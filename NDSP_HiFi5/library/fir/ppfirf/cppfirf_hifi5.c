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
    Complex polyphase FIR filter, floating point
    C code optimized for HiFi5 with FPU/VFPU
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "common_fpu.h"
/* Internal declarations for polyphase FIR filters. */
#include "ppfir_internal.h"

#define ALIGN_SIZE     (HIFI_SIMD_WIDTH)

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

#if HAVE_VFPU
#ifndef AE_DSEL32_LH_SX2
/*
   Equal to:
   a = AE_SEL32_LH_SX2(c, c)
   b = AE_SEL32_LH_SX2(d, d)
*/
#define AE_DSEL32_LH_SX2(a,b,c,d) \
{\
    ae_int16x4 aa,bb,cc,dd,sel; \
    sel = AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32X2((5L<<24)|(1L<<16)|(4<<8)|(0), (7L<<24)|(3L<<16)|(6<<8)|(2) )); \
    cc = AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(c)); \
    dd = AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(d)); \
    AE_DSEL16X4(aa,bb,cc,dd,sel); \
    a = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(aa)); \
    b = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(bb)); \
}
#endif

complex_float * cppfirf(complex_float * restrict y,
                        complex_float * restrict d,
                        complex_float * restrict p,
                  const float32_t     * restrict h,
                  const complex_float * restrict x,
                  int M, int N)
{
    xtfloatx4 * restrict pY;
    const xtfloatx4 * restrict pX;
          xtfloatx4 * restrict pD0;
          xtfloatx4 * restrict pD1;
    const xtfloatx4 * restrict pH0;
    const xtfloatx4 * restrict pH1;
    const xtfloatx4 * restrict pH2;
    const xtfloatx4 * restrict pH3;
    xtfloatx2 h01, h23, h45, h67, h89, hAB, hCD, hEF;
    xtfloatx2 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, xA, xB, xC, xD, xE, xF;
    xtfloatx2 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, YA, YB, YC, YD, YE, YF;
    ae_valignx2 alY;
    int m, n;

    NASSERT_ALIGN(d, ALIGN_SIZE);
    NASSERT_ALIGN(p, ALIGN_SIZE);
    NASSERT_ALIGN(h, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(0==(M%32) && 32<=M && M<=640);
    NASSERT(0==(N%2) && 4<=N && N<=24);

    WUR_AE_CBEGIN0((uintptr_t)(d));
    WUR_AE_CEND0((uintptr_t)(d + M*N));
    /* Update the delay line and circularly update the current position. */
    pX = (const xtfloatx4 *)(x);
    pD0 = (     xtfloatx4 *)(p);
    __Pragma("loop_count min=1");
    for (m = 0; m < (M>>1); m++)
    {
        AE_LSX2X2_IP(x0, x1, pX , sizeof(xtfloatx4));
        AE_SSX2X2_XC(x0, x1, pD0, sizeof(xtfloatx4));
    }
    p = (complex_float *)pD0;

    /* Compute output sample for each of M subfilters. */
    pY = (xtfloatx4 *)y;
    alY = AE_ZALIGN128();
    __Pragma("loop_count min=1");
    for (m=0; m<M; m+=16)
    {
        pD0 = (xtfloatx4 *)(p + m);
        pD1 = pD0 + 1;
        pH0 = (const xtfloatx4 *)(h + m);
        pH1 = pH0 + 1;
        pH2 = pH0 + 2;
        pH3 = pH0 + 3;
        CONST_SX2X2(Y0, Y1, 0);
        CONST_SX2X2(Y2, Y3, 0);
        CONST_SX2X2(Y4, Y5, 0);
        CONST_SX2X2(Y6, Y7, 0);
        CONST_SX2X2(Y8, Y9, 0);
        CONST_SX2X2(YA, YB, 0);
        CONST_SX2X2(YC, YD, 0);
        CONST_SX2X2(YE, YF, 0);
        __Pragma("loop_count min=1, factor=2");
        for (n=0; n<N; n++)
        {
            AE_LSX2X2_XP(h01, h23, pH0, M*sizeof(float32_t));
            AE_LSX2X2_XP(h45, h67, pH1, M*sizeof(float32_t));
            AE_LSX2X2_XP(h89, hAB, pH2, M*sizeof(float32_t));
            AE_LSX2X2_XP(hCD, hEF, pH3, M*sizeof(float32_t));
            __Pragma("no_reorder");

            AE_LSX2X2_IP(x0, x1, pD0, 4*sizeof(complex_float));
            AE_LSX2X2_IP(x2, x3, pD1, 4*sizeof(complex_float));
            AE_LSX2X2_IP(x4, x5, pD0, 4*sizeof(complex_float));
            AE_LSX2X2_IP(x6, x7, pD1, 4*sizeof(complex_float));
            AE_LSX2X2_IP(x8, x9, pD0, 4*sizeof(complex_float));
            AE_LSX2X2_IP(xA, xB, pD1, 4*sizeof(complex_float));
            AE_LSX2X2_XC(xC, xD, pD0, (M-12)*sizeof(complex_float));
            AE_LSX2X2_XC(xE, xF, pD1, (M-12)*sizeof(complex_float));

            MADDMUX_SX2X2(Y0, Y2, h01, h23, x0, x2, 0);
            MADDMUX_SX2X2(Y1, Y3, h01, h23, x1, x3, 5);
            MADDMUX_SX2X2(Y4, Y6, h45, h67, x4, x6, 0);
            MADDMUX_SX2X2(Y5, Y7, h45, h67, x5, x7, 5);
            MADDMUX_SX2X2(Y8, YA, h89, hAB, x8, xA, 0);
            MADDMUX_SX2X2(Y9, YB, h89, hAB, x9, xB, 5);
            MADDMUX_SX2X2(YC, YE, hCD, hEF, xC, xE, 0);
            MADDMUX_SX2X2(YD, YF, hCD, hEF, xD, xF, 5);
        }
        AE_DSEL32_LH_SX2(Y1, Y3, Y1, Y3);
        AE_DSEL32_LH_SX2(Y5, Y7, Y5, Y7);
        AE_DSEL32_LH_SX2(Y9, YB, Y9, YB);
        AE_DSEL32_LH_SX2(YD, YF, YD, YF);
        AE_SASX2X2_IP(Y0, Y1, alY, pY);
        AE_SASX2X2_IP(Y2, Y3, alY, pY);
        AE_SASX2X2_IP(Y4, Y5, alY, pY);
        AE_SASX2X2_IP(Y6, Y7, alY, pY);
        AE_SASX2X2_IP(Y8, Y9, alY, pY);
        AE_SASX2X2_IP(YA, YB, alY, pY);
        AE_SASX2X2_IP(YC, YD, alY, pY);
        AE_SASX2X2_IP(YE, YF, alY, pY);
    }
    AE_SA128POS_FP(alY, pY);
    return p;
} /* cppfirf() */
#elif HAVE_FPU
complex_float * cppfirf(complex_float * restrict y,
                        complex_float * restrict d,
                        complex_float * restrict p,
                  const float32_t     * restrict h,
                  const complex_float * restrict x,
                  int M, int N)
{
    xtfloat * restrict py;
    const xtfloat * restrict ph;
    const xtfloat * restrict pd;
    const ae_int32x4 * restrict pX;
          ae_int32x4 * restrict pD;
    int m, n;
    ae_int32x2 t0, t1;
    float32_t h0, h1, h2, h3;
    float32_t d0re, d1re, d2re, d3re;
    float32_t d0im, d1im, d2im, d3im;
    float32_t y0re, y1re, y2re, y3re;
    float32_t y0im, y1im, y2im, y3im;

    NASSERT_ALIGN(d, ALIGN_SIZE);
    NASSERT_ALIGN(p, ALIGN_SIZE);
    NASSERT_ALIGN(h, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(0==(M%32) && 32<=M && M<=640);
    NASSERT(0==(N%2) && 4<=N && N<=24);

    WUR_AE_CBEGIN0((uintptr_t)(d));
    WUR_AE_CEND0((uintptr_t)(d + M*N));
    /* Update the delay line and circularly update the current position. */
    pX = (const ae_int32x4 *)(x);
    pD = (      ae_int32x4 *)(p);
    __Pragma("loop_count min=1");
    for (m = 0; m < (M>>1); m++)
    {
        AE_L32X2X2_IP(t0, t1, pX, 2*sizeof(complex_float));
        AE_S32X2X2_XC(t0, t1, pD, 2*sizeof(complex_float));
    }
    p = (complex_float *)pD;

    /* Compute output sample for each of M subfilters. */
    py = (xtfloat *)(y);
    __Pragma("loop_count min=1");
    for (m=0; m<M; m+=4)
    {
        y0re = y1re = y2re = y3re = 0.0f;
        y0im = y1im = y2im = y3im = 0.0f;
        ph = (const xtfloat *)(h + m);
        pd = (const xtfloat *)(p + m);
        __Pragma("loop_count min=1");
        for (n=0; n<N; n++)
        {
            XT_LSIP(h0, ph, sizeof(float32_t));
            XT_LSIP(h1, ph, sizeof(float32_t));
            XT_LSIP(h2, ph, sizeof(float32_t));
            XT_LSXP(h3, ph, (M-3)*sizeof(float32_t));

            XT_LSIP(d0re, pd, sizeof(float32_t));
            XT_LSIP(d0im, pd, sizeof(float32_t));
            XT_LSIP(d1re, pd, sizeof(float32_t));
            XT_LSIP(d1im, pd, sizeof(float32_t));
            XT_LSIP(d2re, pd, sizeof(float32_t));
            XT_LSIP(d2im, pd, sizeof(float32_t));
            XT_LSIP(d3re, pd, sizeof(float32_t));
            XT_LSXC(d3im, pd, (M*2-7)*sizeof(float32_t));

            XT_MADD_S(y0re, h0, d0re);
            XT_MADD_S(y1re, h1, d1re);
            XT_MADD_S(y2re, h2, d2re);
            XT_MADD_S(y3re, h3, d3re);
            XT_MADD_S(y0im, h0, d0im);
            XT_MADD_S(y1im, h1, d1im);
            XT_MADD_S(y2im, h2, d2im);
            XT_MADD_S(y3im, h3, d3im);
        }
        XT_SSIP(y0re, py, sizeof(xtfloat));
        XT_SSIP(y0im, py, sizeof(xtfloat));
        XT_SSIP(y1re, py, sizeof(xtfloat));
        XT_SSIP(y1im, py, sizeof(xtfloat));
        XT_SSIP(y2re, py, sizeof(xtfloat));
        XT_SSIP(y2im, py, sizeof(xtfloat));
        XT_SSIP(y3re, py, sizeof(xtfloat));
        XT_SSIP(y3im, py, sizeof(xtfloat));
    }
    return p;
} /* cppfirf() */
#endif /* HAVE_FPU || HAVE_VFPU */

