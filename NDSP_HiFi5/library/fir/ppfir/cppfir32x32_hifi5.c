/* ------------------------------------------------------------------------ */
/* Copyright (c) 2021 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
/* These coded instructions, statements, and computer programs ('Cadence    */
/* Libraries') are the copyrighted works of Cadence Design Systems Inc.     */
/* Cadence IP is licensed for use with Cadence processor cores only and     */
/* must not be used for any other processors and platforms. Your use of the */
/* Cadence Libraries is subject to the terms of the license agreement you   */
/* have entered into with Cadence Design Systems, or a sublicense granted   */
/* to you by a direct Cadence license.                                     */
/* ------------------------------------------------------------------------ */
/*  IntegrIT, Ltd.   www.integrIT.com, info@integrIT.com                    */
/*                                                                          */
/* NatureDSP_Baseband Library                                               */
/*                                                                          */
/* This library contains copyrighted materials, trade secrets and other     */
/* proprietary information of IntegrIT, Ltd. This software is licensed for  */
/* use with Cadence processor cores only and must not be used for any other */
/* processors and platforms. The license to use these sources was given to  */
/* Cadence, Inc. under Terms and Condition of a Software License Agreement  */
/* between Cadence, Inc. and IntegrIT, Ltd.                                 */
/* ------------------------------------------------------------------------ */
/*          Copyright (C) 2009-2021 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */

/*
  NatureDSP Signal Processing Library. FIR part
    Complex polyphase FIR filter, 32x32-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Internal declarations for polyphase FIR filters. */
#include "ppfir_internal.h"

#define ALIGN_SIZE     (HIFI_SIMD_WIDTH)

#ifndef AE_DSEL32X2_HH_LL
/*
   Equal to:
   a = AE_SEL32_HH(c, d)
   b = AE_SEL32_LL(c, d)
*/
#define AE_DSEL32X2_HH_LL(a,b,c,d) \
{\
    ae_int16x4 aa,bb,cc,dd,sel; \
    sel = AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32X2((7L<<24)|(5L<<16)|(6<<8)|(4), (3L<<24)|(1L<<16)|(2<<8)|(0) )); \
    cc = AE_MOVINT16X4_FROMINT32X2(c); \
    dd = AE_MOVINT16X4_FROMINT32X2(d); \
    AE_DSEL16X4(aa,bb,cc,dd,sel); \
    a = AE_MOVINT32X2_FROMINT16X4(aa); \
    b = AE_MOVINT32X2_FROMINT16X4(bb); \
}
#endif

/*
 * Polyphase FIR filter kernel.
 * This function performs polyphase filtering for a block of M samples
 * specified via the input argument x[M] (type-1 polyphase decomposition).
 * This is accomplished in two steps:
 *   1. Insert M new samples into the current position (input argument p)
 *      of the circular delay line (input/output argument d) of size N*M. Current
 *      position p is updated by M sample positions in a round-robin fashion,
 *      and will be returned to the caller for subsequent calls of this function.
 *   2. Compute output samples of M subfilters, each of N taps (input argument 
 *      h[N*M]), and store results to the output argument y[M].  
 * Domain:
 *   rpp      Real data, real coefficients
 *   cpp      Complex data, real coefficients
 * Data type and precision
 *   pp32x32  32-bit fixed-point data, 32-bit fixed-point coefficients
 *   ppf      Floating point
 * Input:
 *   M        Number of subfilters (or the number of phases)
 *   N        Subfilter length
 *   lsh      (fixed-point variants only) Bi-directional saturating left shift amount
 *            to be applied to input samples prior to storing them into the delay line.
 *            For a negative (i.e. right) shift amount, results are asymmetrically
 *            rounded.
 *   x[M]     A block of input samples, assuming type-1 polyphase decomposition.
 *            That is, x[0] is the oldest sample in a block, and x[M-1] is the 
 *            newest one.
 *   h[N*M]   Subfilters' coefficients organized into an N-by-M matrix, where N is
 *            the number of subfilter taps, and M is the number of subfilters. That
 *            is, N coefficients of the m-th subfilter are stored to the m-th column:
 *            h[(0..N-1)*M+m], where coefficient h[0*M+m] corresponds to the oldest
 *            input sample in subfilter's delay line.
 * Input/Output:
 *   d[N*M]   Subfilters' delay lines organized into an N-by-M matrix, analogously to
 *            subfilter coefficients storage h[N*M]. That is, N-sample delay line of
 *            the m-th subfilter occupies the m-th column: d[(0..N-1)*M+m]. 
 *   p        Current position in the delay line, that is a pointer to one of N rows
 *            of the N-by-M matrix d[N*M]. In a single call of this function, this
 *            row gets filled with M input samples x[M], and the current position 
 *            switches to the next row. After the last row of the delay line is 
 *            filled with new data, the current position wraps around to the 0th row.
 *            The function returns the updated current position to be used for 
 *            subsequent calls.
 * Output:
 *   y[M]     Output samples of M subfilters. Output sample of 0th subfilter is stored
 *            to y[0], output sample of 1st subfilter - to y[1], and so forth, up to
 *            M-1st subfilter's output sample stored to y[M-1].
 * Restrictions:
 *   d,p,h,x  Must not overlap, and must be 16-byte aligned.
 *   N        4..24, multiple of 2
 *   M        32...640, multiple of 32
 */

complex_fract32 * cppfir32x32(complex_fract32 * restrict y,
                              complex_fract32 * restrict d,
                              complex_fract32 * restrict p,
                        const int32_t         * restrict h,
                        const complex_fract32 * restrict x,
                        int M, int N, int lsh)
{
    ae_int32x4 * restrict pY;
    const ae_int32x4 * restrict pX;
          ae_int32x4 * restrict pD0;
          ae_int32x4 * restrict pD1;
    const ae_int32x2 * restrict pH;
    ae_int32x2 h001, h101, h201, h301, h401, h501, h601, h701, scl;
    ae_int32x2 x00, x10, x20, x30, x40, x50, x60, x70;
    ae_int32x2 x01, x11, x21, x31, x41, x51, x61, x71;
    ae_int32x2 y0, y1, y2, y3, y4, y5, y6, y7;
    ae_int64 Y0re, Y1re, Y2re, Y3re, Y4re, Y5re, Y6re, Y7re;
    ae_int64 Y0im, Y1im, Y2im, Y3im, Y4im, Y5im, Y6im, Y7im;
    ae_valignx2 alY;
    int m, n;

    NASSERT_ALIGN(d, ALIGN_SIZE);
    NASSERT_ALIGN(p, ALIGN_SIZE);
    NASSERT_ALIGN(h, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(0==(M%32) && 32<=M && M<=640);
    NASSERT(0==(N%2) && 4<=N && N<=24);

    WUR_AE_SAR(lsh);
    WUR_AE_CBEGIN0((uintptr_t)(d));
    WUR_AE_CEND0((uintptr_t)(d + M*N));
    /* Update the delay line and circularly update the current position. */
    pX = (const ae_int32x4 *)(x);
    pD0 = (     ae_int32x4 *)(p);
    if (lsh >= 0)
    {
        scl = AE_MOVDA32(1L<<lsh);
        __Pragma("loop_count min=4, factor=2");
        for (m = 0; m < (M>>1); m++)
        {
            AE_L32X2X2_IP(x00, x01, pX, sizeof(ae_int32x4));
            AE_MUL2P32X4S(x00, x01, x00, x01, scl, scl);
            AE_S32X2X2_XC(x00, x01, pD0, sizeof(ae_int32x4));
        }
    }
    else
    {
        scl = AE_MOVDA32(1L<<(31+lsh));
        __Pragma("loop_count min=4, factor=2");
        for (m = 0; m < (M>>1); m++)
        {
            AE_L32X2X2_IP(x00, x01, pX, sizeof(ae_int32x4));
            AE_MULF2P32X4RAS(x00, x01, x00, x01, scl, scl);
            AE_S32X2X2_XC(x00, x01, pD0, sizeof(ae_int32x4));
        }
    }
    p = (complex_fract32 *)pD0;

    /* Compute output sample for each of M subfilters. */
    pY = (ae_int32x4 *)y;
    alY = AE_ZALIGN128();
    __Pragma("loop_count min=1");
    for (m=0; m<M; m+=8)
    {
        pD0 = (ae_int32x4 *)(p + m);
        pD1 = pD0; AE_ADDCIRC_XC(castxcc(ae_int64,pD1), M*sizeof(complex_fract32));
        pH = (const ae_int32x2 *)(h + m*2);
        Y0re = Y1re = Y2re = Y3re = Y4re = Y5re = Y6re = Y7re = AE_ZERO64();
        Y0im = Y1im = Y2im = Y3im = Y4im = Y5im = Y6im = Y7im = AE_ZERO64();
        __Pragma("loop_count min=1");
        __Pragma("ymemory(pD0)");
        __Pragma("ymemory(pD1)");
        for (n=0; n<(N>>1); n++)
        {
            AE_L32X2_XP(h001, pH, 2*sizeof(int32_t));
            AE_L32X2_XP(h101, pH, 2*sizeof(int32_t));
            AE_L32X2_XP(h201, pH, 2*sizeof(int32_t));
            AE_L32X2_XP(h301, pH, 2*sizeof(int32_t));
            AE_L32X2_XP(h401, pH, 2*sizeof(int32_t));
            AE_L32X2_XP(h501, pH, 2*sizeof(int32_t));
            AE_L32X2_XP(h601, pH, 2*sizeof(int32_t));
            AE_L32X2_XP(h701, pH, (M*2-2*7)*sizeof(int32_t));

            AE_L32X2X2_IP(x00, x10, pD0, 2*sizeof(complex_fract32));
            AE_L32X2X2_IP(x20, x30, pD0, 2*sizeof(complex_fract32));
            AE_L32X2X2_IP(x40, x50, pD0, 2*sizeof(complex_fract32));
            AE_L32X2X2_XC(x60, x70, pD0, (M*2-6)*sizeof(complex_fract32));
            AE_L32X2X2_IP(x01, x11, pD1, 2*sizeof(complex_fract32));
            AE_L32X2X2_IP(x21, x31, pD1, 2*sizeof(complex_fract32));
            AE_L32X2X2_IP(x41, x51, pD1, 2*sizeof(complex_fract32));
            AE_L32X2X2_XC(x61, x71, pD1, (M*2-6)*sizeof(complex_fract32));
            AE_DSEL32X2_HH_LL(x00, x01, x00, x01);
            AE_DSEL32X2_HH_LL(x10, x11, x10, x11);
            AE_DSEL32X2_HH_LL(x20, x21, x20, x21);
            AE_DSEL32X2_HH_LL(x30, x31, x30, x31);
            AE_DSEL32X2_HH_LL(x40, x41, x40, x41);
            AE_DSEL32X2_HH_LL(x50, x51, x50, x51);
            AE_DSEL32X2_HH_LL(x60, x61, x60, x61);
            AE_DSEL32X2_HH_LL(x70, x71, x70, x71);

            AE_MULAAF2D32RA_HH_LL(Y0re, Y0im, h001, h001, x00, x01);
            AE_MULAAF2D32RA_HH_LL(Y1re, Y1im, h101, h101, x10, x11);
            AE_MULAAF2D32RA_HH_LL(Y2re, Y2im, h201, h201, x20, x21);
            AE_MULAAF2D32RA_HH_LL(Y3re, Y3im, h301, h301, x30, x31);
            AE_MULAAF2D32RA_HH_LL(Y4re, Y4im, h401, h401, x40, x41);
            AE_MULAAF2D32RA_HH_LL(Y5re, Y5im, h501, h501, x50, x51);
            AE_MULAAF2D32RA_HH_LL(Y6re, Y6im, h601, h601, x60, x61);
            AE_MULAAF2D32RA_HH_LL(Y7re, Y7im, h701, h701, x70, x71);
        }
        y0 = AE_ROUND32X2F48SASYM(Y0re, Y0im);
        y1 = AE_ROUND32X2F48SASYM(Y1re, Y1im);
        y2 = AE_ROUND32X2F48SASYM(Y2re, Y2im);
        y3 = AE_ROUND32X2F48SASYM(Y3re, Y3im);
        y4 = AE_ROUND32X2F48SASYM(Y4re, Y4im);
        y5 = AE_ROUND32X2F48SASYM(Y5re, Y5im);
        y6 = AE_ROUND32X2F48SASYM(Y6re, Y6im);
        y7 = AE_ROUND32X2F48SASYM(Y7re, Y7im);
        AE_SA32X2X2_IP(y0, y1, alY, pY);
        AE_SA32X2X2_IP(y2, y3, alY, pY);
        AE_SA32X2X2_IP(y4, y5, alY, pY);
        AE_SA32X2X2_IP(y6, y7, alY, pY);
    }
    AE_SA128POS_FP(alY, pY);
    return p;
} /* cppfir32x32() */
