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
    Real polyphase FIR filter, floating point
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
float32_t * rppfirf(float32_t * restrict y,
                    float32_t * restrict d,
                    float32_t * restrict p,
              const float32_t * restrict h,
              const float32_t * restrict x,
              int M, int N)
{
    xtfloatx4 * restrict pY;
    const xtfloatx4 * restrict pX;
          xtfloatx4 * restrict pD;
    const xtfloatx4 * restrict pH;
    xtfloatx2 h001, h023, h045, h067, h089, h0AB, h0CD, h0EF;
    xtfloatx2 x001, x023, x045, x067, x089, x0AB, x0CD, x0EF;
    xtfloatx2 Y001, Y023, Y045, Y067, Y089, Y0AB, Y0CD, Y0EF;
    xtfloatx2 h101, h123, h145, h167, h189, h1AB, h1CD, h1EF;
    xtfloatx2 x101, x123, x145, x167, x189, x1AB, x1CD, x1EF;
    xtfloatx2 Y101, Y123, Y145, Y167, Y189, Y1AB, Y1CD, Y1EF;
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
    pD = (      xtfloatx4 *)(p);
    __Pragma("loop_count min=1");
    for (m = 0; m < (M>>2); m++)
    {
        AE_LSX2X2_IP(x001, x023, pX, sizeof(xtfloatx4));
        AE_SSX2X2_XC(x001, x023, pD, sizeof(xtfloatx4));
    }
    p = (float32_t *)pD;

    /* Compute output sample for each of M subfilters. */
    pY = (xtfloatx4 *)y;
    alY = AE_ZALIGN128();
    __Pragma("loop_count min=1");
    for (m=0; m<M; m+=32)
    {
        pD = (xtfloatx4 *)(p + m);
        pH = (const xtfloatx4 *)(h + m);
        CONST_SX2X2(Y001, Y023, 0);
        CONST_SX2X2(Y045, Y067, 0);
        CONST_SX2X2(Y089, Y0AB, 0);
        CONST_SX2X2(Y0CD, Y0EF, 0);
        CONST_SX2X2(Y101, Y123, 0);
        CONST_SX2X2(Y145, Y167, 0);
        CONST_SX2X2(Y189, Y1AB, 0);
        CONST_SX2X2(Y1CD, Y1EF, 0);
        __Pragma("loop_count min=1, factor=2");
        for (n=0; n<N; n++)
        {
            AE_LSX2X2_IP(h001, h023, pH, 4*sizeof(float32_t));
            AE_LSX2X2_IP(h045, h067, pH, 4*sizeof(float32_t));
            AE_LSX2X2_IP(h089, h0AB, pH, 4*sizeof(float32_t));
            AE_LSX2X2_IP(h0CD, h0EF, pH, 4*sizeof(float32_t));
            AE_LSX2X2_IP(h101, h123, pH, 4*sizeof(float32_t));
            AE_LSX2X2_IP(h145, h167, pH, 4*sizeof(float32_t));
            AE_LSX2X2_IP(h189, h1AB, pH, 4*sizeof(float32_t));
            AE_LSX2X2_XP(h1CD, h1EF, pH, (M-7*4)*sizeof(float32_t));

            AE_LSX2X2_IP(x001, x023, pD, 4*sizeof(float32_t));
            AE_LSX2X2_IP(x045, x067, pD, 4*sizeof(float32_t));
            AE_LSX2X2_IP(x089, x0AB, pD, 4*sizeof(float32_t));
            AE_LSX2X2_IP(x0CD, x0EF, pD, 4*sizeof(float32_t));
            AE_LSX2X2_IP(x101, x123, pD, 4*sizeof(float32_t));
            AE_LSX2X2_IP(x145, x167, pD, 4*sizeof(float32_t));
            AE_LSX2X2_IP(x189, x1AB, pD, 4*sizeof(float32_t));
            AE_LSX2X2_XC(x1CD, x1EF, pD, (M-7*4)*sizeof(float32_t));

            MADD_SX2X2(Y001, Y023, h001, h023, x001, x023);
            MADD_SX2X2(Y045, Y067, h045, h067, x045, x067);
            MADD_SX2X2(Y089, Y0AB, h089, h0AB, x089, x0AB);
            MADD_SX2X2(Y0CD, Y0EF, h0CD, h0EF, x0CD, x0EF);
            MADD_SX2X2(Y101, Y123, h101, h123, x101, x123);
            MADD_SX2X2(Y145, Y167, h145, h167, x145, x167);
            MADD_SX2X2(Y189, Y1AB, h189, h1AB, x189, x1AB);
            MADD_SX2X2(Y1CD, Y1EF, h1CD, h1EF, x1CD, x1EF);
        }
        AE_SASX2X2_IP(Y001, Y023, alY, pY);
        AE_SASX2X2_IP(Y045, Y067, alY, pY);
        AE_SASX2X2_IP(Y089, Y0AB, alY, pY);
        AE_SASX2X2_IP(Y0CD, Y0EF, alY, pY);
        AE_SASX2X2_IP(Y101, Y123, alY, pY);
        AE_SASX2X2_IP(Y145, Y167, alY, pY);
        AE_SASX2X2_IP(Y189, Y1AB, alY, pY);
        AE_SASX2X2_IP(Y1CD, Y1EF, alY, pY);
    }
    AE_SA128POS_FP(alY, pY);
    return p;
} /* rppfirf() */
#elif HAVE_FPU
float32_t * rppfirf(float32_t * restrict y,
                    float32_t * restrict d,
                    float32_t * restrict p,
              const float32_t * restrict h,
              const float32_t * restrict x,
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
    float32_t d0, d1, d2, d3;
    float32_t y0, y1, y2, y3;

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
    for (m = 0; m < (M>>2); m++)
    {
        AE_L32X2X2_IP(t0, t1, pX, 4*sizeof(float32_t));
        AE_S32X2X2_XC(t0, t1, pD, 4*sizeof(float32_t));
    }
    p = (float32_t *)pD;

    /* Compute output sample for each of M subfilters. */
    py = (xtfloat *)(y);
    __Pragma("loop_count min=1");
    for (m=0; m<M; m+=4)
    {
        y0 = y1 = y2 = y3 = 0.0f;
        ph = (const xtfloat *)(h + m);
        pd = (const xtfloat *)(p + m);
        __Pragma("loop_count min=1");
        for (n=0; n<N; n++)
        {
            XT_LSIP(h0, ph, sizeof(float32_t));
            XT_LSIP(h1, ph, sizeof(float32_t));
            XT_LSIP(h2, ph, sizeof(float32_t));
            XT_LSXP(h3, ph, (M-3)*sizeof(float32_t));

            XT_LSIP(d0, pd, sizeof(float32_t));
            XT_LSIP(d1, pd, sizeof(float32_t));
            XT_LSIP(d2, pd, sizeof(float32_t));
            XT_LSXC(d3, pd, (M-3)*sizeof(float32_t));

            XT_MADD_S(y0, h0, d0);
            XT_MADD_S(y1, h1, d1);
            XT_MADD_S(y2, h2, d2);
            XT_MADD_S(y3, h3, d3);
        }
        XT_SSIP(y0, py, sizeof(xtfloat));
        XT_SSIP(y1, py, sizeof(xtfloat));
        XT_SSIP(y2, py, sizeof(xtfloat));
        XT_SSIP(y3, py, sizeof(xtfloat));
    }
    return p;
} /* rppfirf() */
#endif /* HAVE_FPU || HAVE_VFPU */
