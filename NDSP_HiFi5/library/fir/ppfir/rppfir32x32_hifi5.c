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
    Real polyphase FIR filter, 32x32-bit
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

int32_t * rppfir32x32(int32_t * restrict y,
                      int32_t * restrict d,
                      int32_t * restrict p,
                const int32_t * restrict h,
                const int32_t * restrict x,
                int M, int N, int lsh)
{
    ae_int32x4 * restrict pY;
    const ae_int32x4 * restrict pX;
          ae_int32x4 * restrict pD0;
          ae_int32x4 * restrict pD1;
    const ae_int32x4 * restrict pH0;
    const ae_int32x4 * restrict pH1;
    ae_int32x2 h0, h1, h2, h3, h4, h5, h6, h7, scl;
    ae_int32x2 h8, h9, hA, hB, hC, hD, hE, hF;
    ae_int32x2 x010, x230, x450, x670, x890, xAB0, xCD0, xEF0;
    ae_int32x2 x011, x231, x451, x671, x891, xAB1, xCD1, xEF1;
    ae_int32x2 x0, x1, x2, x3, x4, x5, x6, x7;
    ae_int32x2 x8, x9, xA, xB, xC, xD, xE, xF;
    ae_int32x2 y01, y23, y45, y67, y89, yAB, yCD, yEF;
    ae_int64 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
    ae_int64 Y8, Y9, YA, YB, YC, YD, YE, YF;
    ae_valignx2 alY;
    int m, n, Md;

    NASSERT_ALIGN(d, ALIGN_SIZE);
    NASSERT_ALIGN(p, ALIGN_SIZE);
    NASSERT_ALIGN(h, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(0==(M%32) && 32<=M && M<=640);
    NASSERT(0==(N%2) && 4<=N && N<=24);

    Md = M + DELAY_PAD_I32;
    WUR_AE_SAR(lsh);
    WUR_AE_CBEGIN0((uintptr_t)(d));
    WUR_AE_CEND0((uintptr_t)(d + Md*N));
    /* Update the delay line and circularly update the current position. */
    pX = (const ae_int32x4 *)(x);
    pD0 = (     ae_int32x4 *)(p);
    if (lsh >= 0)
    {
        scl = AE_MOVDA32(1L<<lsh);
        __Pragma("loop_count min=4, factor=2");
        for (m = 0; m < (M>>2); m++)
        {
            AE_L32X2X2_IP(x010, x230, pX, sizeof(ae_int32x4));
            AE_MUL2P32X4S(x010, x230, x010, x230, scl, scl);
            AE_S32X2X2_XC(x010, x230, pD0, sizeof(ae_int32x4));
        }
        x010 = AE_ZERO32();
        AE_S32X2X2_XC(x010, x010, pD0, sizeof(ae_int32x4));
    }
    else
    {
        scl = AE_MOVDA32(1L<<(31+lsh));
        __Pragma("loop_count min=4, factor=2");
        for (m = 0; m < (M>>2); m++)
        {
            AE_L32X2X2_IP(x010, x230, pX, sizeof(ae_int32x4));
            AE_MULF2P32X4RAS(x010, x230, x010, x230, scl, scl);
            AE_S32X2X2_XC(x010, x230, pD0, sizeof(ae_int32x4));
        }
        x010 = AE_ZERO32();
        AE_S32X2X2_XC(x010, x010, pD0, sizeof(ae_int32x4));
    }
    p = (int32_t *)pD0;

    /* Compute output sample for each of M subfilters. */
    pY = (ae_int32x4 *)y;
    alY = AE_ZALIGN128();
    __Pragma("loop_count min=1");
    for (m=0; m<M; m+=16)
    {
        pD0 = (ae_int32x4 *)(p + m);
        pD1 = pD0; AE_ADDCIRC_XC(castxcc(ae_int64,pD1), Md*sizeof(int32_t));
        pH0 = (const ae_int32x4 *)(h + m*2);
        pH1 = pH0 + 1;
        Y0 = Y1 = Y2 = Y3 = Y4 = Y5 = Y6 = Y7 = AE_ZERO64();
        Y8 = Y9 = YA = YB = YC = YD = YE = YF = AE_ZERO64();
        __Pragma("loop_count min=1");
        for (n=0; n<(N>>1); n++)
        {
            AE_L32X2X2_IP(h0, h1, pH0, 8*sizeof(int32_t));
            AE_L32X2X2_IP(h2, h3, pH1, 8*sizeof(int32_t));
            AE_L32X2X2_IP(h4, h5, pH0, 8*sizeof(int32_t));
            AE_L32X2X2_IP(h6, h7, pH1, 8*sizeof(int32_t));
            AE_L32X2X2_IP(h8, h9, pH0, 8*sizeof(int32_t));
            AE_L32X2X2_IP(hA, hB, pH1, 8*sizeof(int32_t));
            AE_L32X2X2_XP(hC, hD, pH0, (M*2-8*3)*sizeof(int32_t));
            AE_L32X2X2_XP(hE, hF, pH1, (M*2-8*3)*sizeof(int32_t));

            AE_L32X2X2_IP(x010, x230, pD0, 4*sizeof(int32_t));
            AE_L32X2X2_IP(x450, x670, pD0, 4*sizeof(int32_t));
            AE_L32X2X2_IP(x890, xAB0, pD0, 4*sizeof(int32_t));
            AE_L32X2X2_XC(xCD0, xEF0, pD0, (Md*2-4*3)*sizeof(int32_t));
            AE_L32X2X2_IP(x011, x231, pD1, 4*sizeof(int32_t));
            AE_L32X2X2_IP(x451, x671, pD1, 4*sizeof(int32_t));
            AE_L32X2X2_IP(x891, xAB1, pD1, 4*sizeof(int32_t));
            AE_L32X2X2_XC(xCD1, xEF1, pD1, (Md*2-4*3)*sizeof(int32_t));
            AE_DSEL32X2_HH_LL(x0, x1, x010, x011);
            AE_DSEL32X2_HH_LL(x2, x3, x230, x231);
            AE_DSEL32X2_HH_LL(x4, x5, x450, x451);
            AE_DSEL32X2_HH_LL(x6, x7, x670, x671);
            AE_DSEL32X2_HH_LL(x8, x9, x890, x891);
            AE_DSEL32X2_HH_LL(xA, xB, xAB0, xAB1);
            AE_DSEL32X2_HH_LL(xC, xD, xCD0, xCD1);
            AE_DSEL32X2_HH_LL(xE, xF, xEF0, xEF1);

            AE_MULAAF2D32RA_HH_LL(Y0, Y1, h0, h1, x0, x1);
            AE_MULAAF2D32RA_HH_LL(Y2, Y3, h2, h3, x2, x3);
            AE_MULAAF2D32RA_HH_LL(Y4, Y5, h4, h5, x4, x5);
            AE_MULAAF2D32RA_HH_LL(Y6, Y7, h6, h7, x6, x7);
            AE_MULAAF2D32RA_HH_LL(Y8, Y9, h8, h9, x8, x9);
            AE_MULAAF2D32RA_HH_LL(YA, YB, hA, hB, xA, xB);
            AE_MULAAF2D32RA_HH_LL(YC, YD, hC, hD, xC, xD);
            AE_MULAAF2D32RA_HH_LL(YE, YF, hE, hF, xE, xF);
        }
        y01 = AE_ROUND32X2F48SASYM(Y0, Y1);
        y23 = AE_ROUND32X2F48SASYM(Y2, Y3);
        y45 = AE_ROUND32X2F48SASYM(Y4, Y5);
        y67 = AE_ROUND32X2F48SASYM(Y6, Y7);
        y89 = AE_ROUND32X2F48SASYM(Y8, Y9);
        yAB = AE_ROUND32X2F48SASYM(YA, YB);
        yCD = AE_ROUND32X2F48SASYM(YC, YD);
        yEF = AE_ROUND32X2F48SASYM(YE, YF);
        AE_SA32X2X2_IP(y01, y23, alY, pY);
        AE_SA32X2X2_IP(y45, y67, alY, pY);
        AE_SA32X2X2_IP(y89, yAB, alY, pY);
        AE_SA32X2X2_IP(yCD, yEF, alY, pY);
    }
    AE_SA128POS_FP(alY, pY);
    return p;
} /* rppfir32x32() */
