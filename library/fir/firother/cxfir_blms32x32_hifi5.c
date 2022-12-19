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
    Blockwise Adaptive LMS Algorithm for Complex Data, 32x32-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

/*-------------------------------------------------------------------------
  Blockwise Adaptive LMS Algorithm for Real Data
  Blockwise LMS algorithm performs filtering of reference samples x[N+M-1],
  computation of error e[N] over a block of input samples r[N] and makes
  blockwise update of IR to minimize the error output.
  Algorithm includes FIR filtering, calculation of correlation between the 
  error output e[N] and reference signal x[N+M-1] and IR taps update based
  on that correlation.
NOTES: 
  1. The algorithm must be provided with the normalization factor, which is
     the power of the reference signal times N - the number of samples in a
     data block. This can be calculated using the vec_power32x32() or 
     vec_power16x16() function. In order to avoid the saturation of the 
     normalization factor, it may be biased, i.e. shifted to the right.
     If it's the case, then the adaptation coefficient must be also shifted
     to the right by the same number of bit positions.
  2. This algorithm consumes less CPU cycles per block than single 
     sample algorithm at similar convergence rate.
  3. Right selection of N depends on the change rate of impulse response:
     on static or slow varying channels convergence rate depends on
     selected mu and M, but not on N.
  4. 16x16 routine may converge slower on small errors due to roundoff 
     errors. In that cases, 16x32 routine will give better results although
     convergence rate on bigger errors is the same.
  5. Terms near-end and far-end come from echo cancellation theory where the 
     LMS is used widely. For echo cancellation them term far-end means the 
     output of speakerphone (far end designates that the origin of it is 
     somewhere outside say came from the remote speaker). The near-end is 
     a signal from the local microphone representing a sum of the echo, 
     speech of local speaker and the noise. The LMS is used to estimate the 
     equivalent impulse response of the echopath further compensation and 
     removal the echo from the near-end signal.

  Precision: 
  16x16    16-bit coefficients, 16-bit data, 16-bit output
  16x32    32-bit coefficients, 16-bit data, 16-bit output
  32x32    32-bit coefficients, 32-bit data, 32-bit output, complex and real
  32x32ep  the same as above but using 72-bit accumulator for intermediate 
           computations
  f        floating point, complex and real
  Input:
  h[M]     impulse response, Q15, Q31 or floating point
  r[N]	   input data vector (near end). First in time value is in 
           r[0], Q15, Q31 or floating point
  x[N+M-1] reference data vector (far end). First in time value is in x[0],  
           Q15, Q31 or floating point
  norm     normalization factor: power of signal multiplied by N, Q15, Q31  
           or floating point
           Fixed-point format for the 32x16-bit variant: Q(2*x+1-bias)
  mu       adaptation coefficient (LMS step), Q(31-bias) or Q(15-bias)
  N        length of data block
  M        length of h
  Output:
  e[N]     estimated error, Q15, Q31 or floating point
  h[M]     updated impulse response, Q15, Q31 or floating point

  Restriction:
  x,r,h,e  should not overlap
  x,r,h,e  aligned on a 16-bytes boundary
  N,M      multiples of 8 and >0
-------------------------------------------------------------------------*/

void cxfir_blms32x32( complex_fract32 * restrict e,
                      complex_fract32 * restrict h,
                const complex_fract32 * restrict r,
                const complex_fract32 * restrict x,
                int32_t norm, int32_t mu,
                int N, int M)
{
    const ae_int64   * restrict pX;
    const ae_int32x2 * restrict pR;
          ae_int32x2 * restrict pE;
    const ae_int32x2 * restrict pH;
          ae_int32x2 * restrict pHw;
    ae_int32x2 s_frac;
    int        s_exp;
    int n, m;

    NASSERT(e);
    NASSERT(h);
    NASSERT(r);
    NASSERT(x);
    NASSERT_ALIGN(e, 16);
    NASSERT_ALIGN(h, 16);
    NASSERT_ALIGN(r, 16);
    NASSERT_ALIGN(x, 16);
    NASSERT(N > 0 && M > 0);
    NASSERT(M % 8 == 0 && N % 8 == 0);

    //
    // Pass the reference signal x[] through the adaptive filter to obtain the
    // predicted signal and calculate the error, i.e. the distance to the
    // actual input signal r[].
    //
    pR = (const ae_int32x2*)(r);
    pE = (      ae_int32x2*)(e);
    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 3); n++)
    {
        ae_f64 q0_re, q1_re, q2_re, q3_re;
        ae_f64 q0_im, q1_im, q2_im, q3_im;
        ae_f64 q4_re, q5_re, q6_re, q7_re;
        ae_f64 q4_im, q5_im, q6_im, q7_im;
        ae_int32x2 h0, r0, r1;
        ae_int32x2 x0, x1, x2, x3, x4, x5, x6, x7;
        ae_int64 xx;

        pX = (const ae_int64   *)(x + 8 * n);
        pH = (const ae_int32x2 *)(h + M - 1);
        q0_re = q1_re = q2_re = q3_re = AE_ZERO64();
        q0_im = q1_im = q2_im = q3_im = AE_ZERO64();
        q4_re = q5_re = q6_re = q7_re = AE_ZERO64();
        q4_im = q5_im = q6_im = q7_im = AE_ZERO64();
        /* preload data from x */
        AE_L64_IP(xx, pX, sizeof(complex_fract32)); x0 = AE_MOVINT32X2_FROMINT64(xx);
        AE_L64_IP(xx, pX, sizeof(complex_fract32)); x1 = AE_MOVINT32X2_FROMINT64(xx);
        AE_L64_IP(xx, pX, sizeof(complex_fract32)); x2 = AE_MOVINT32X2_FROMINT64(xx);
        AE_L64_IP(xx, pX, sizeof(complex_fract32)); x3 = AE_MOVINT32X2_FROMINT64(xx);
        AE_L64_IP(xx, pX, sizeof(complex_fract32)); x4 = AE_MOVINT32X2_FROMINT64(xx);
        AE_L64_IP(xx, pX, sizeof(complex_fract32)); x5 = AE_MOVINT32X2_FROMINT64(xx);
        AE_L64_IP(xx, pX, sizeof(complex_fract32)); x6 = AE_MOVINT32X2_FROMINT64(xx);

        __Pragma("loop_count min=8, factor=8");
        for (m = 0; m < M; m++)
        {
            /* load data from x */
            AE_L64_IP(xx, pX, sizeof(complex_fract32)); x7 = AE_MOVINT32X2_FROMINT64(xx);
            /* load data from y */
            AE_L32X2_XP(h0, pH, -(int)sizeof(complex_fract32));
            /* compute correlation of 8 values */
            AE_MULAFC32RA(q0_im, q0_re, x0, h0);
            AE_MULAFC32RA(q1_im, q1_re, x1, h0);
            AE_MULAFC32RA(q2_im, q2_re, x2, h0);
            AE_MULAFC32RA(q3_im, q3_re, x3, h0);
            AE_MULAFC32RA(q4_im, q4_re, x4, h0);
            AE_MULAFC32RA(q5_im, q5_re, x5, h0);
            AE_MULAFC32RA(q6_im, q6_re, x6, h0);
            AE_MULAFC32RA(q7_im, q7_re, x7, h0);
            /* shift input line for the next iteration */
            x0 = x1; x1 = x2; x2 = x3;
            x3 = x4; x4 = x5; x5 = x6; x6 = x7;
        }
        /* save computed samples */
        AE_L32X2X2_IP(r0, r1, castxcc(ae_int32x4, pR), 2 * sizeof(complex_fract32));
        r0 = AE_SUB32S(r0, AE_ROUND32X2F48SASYM(q0_re, q0_im));
        r1 = AE_SUB32S(r1, AE_ROUND32X2F48SASYM(q1_re, q1_im));
        AE_S32X2X2_IP(r0, r1, castxcc(ae_int32x4, pE), 2 * sizeof(complex_fract32));

        AE_L32X2X2_IP(r0, r1, castxcc(ae_int32x4, pR), 2 * sizeof(complex_fract32));
        r0 = AE_SUB32S(r0, AE_ROUND32X2F48SASYM(q2_re, q2_im));
        r1 = AE_SUB32S(r1, AE_ROUND32X2F48SASYM(q3_re, q3_im));
        AE_S32X2X2_IP(r0, r1, castxcc(ae_int32x4, pE), 2 * sizeof(complex_fract32));

        AE_L32X2X2_IP(r0, r1, castxcc(ae_int32x4, pR), 2 * sizeof(complex_fract32));
        r0 = AE_SUB32S(r0, AE_ROUND32X2F48SASYM(q4_re, q4_im));
        r1 = AE_SUB32S(r1, AE_ROUND32X2F48SASYM(q5_re, q5_im));
        AE_S32X2X2_IP(r0, r1, castxcc(ae_int32x4, pE), 2 * sizeof(complex_fract32));

        AE_L32X2X2_IP(r0, r1, castxcc(ae_int32x4, pR), 2 * sizeof(complex_fract32));
        r0 = AE_SUB32S(r0, AE_ROUND32X2F48SASYM(q6_re, q6_im));
        r1 = AE_SUB32S(r1, AE_ROUND32X2F48SASYM(q7_re, q7_im));
        AE_S32X2X2_IP(r0, r1, castxcc(ae_int32x4, pE), 2 * sizeof(complex_fract32));
    }

    //
    // Compute the reciprocal for the normalization factor.
    //
    {
        ae_int32x2  mu_v;
        int         mu_exp;
        xtbool2     sgn, inf;
        ae_int32x2  X;
        ae_f32x2    Y, E;

        X = norm;
        sgn = AE_LT32(X, AE_ZERO32());
        inf = AE_EQ32(X, AE_ZERO32());
        X = AE_ABS32S(X);
        s_exp = AE_NSAZ32_L(X);
        X = AE_SLAA32(X, s_exp);//x in 0.5..1 ,Q31
        // first approximation 3-2x,Q30
        Y = AE_SUB32(AE_MOVDA32X2(0xBAEC0000, 0xBAEC0000), X);
        // 1-st iteration
        E = AE_MOVDA32X2(0x40000000, 0x40000000);
        AE_MULSFP32X2RAS(E, X, Y);
        E = AE_SLAI32(E, 1);
        AE_MULAFP32X2RAS(Y, Y, E);
        // 2-st iteration
        E = AE_MOVDA32X2(0x40000000, 0x40000000);
        AE_MULSFP32X2RAS(E, X, Y);
        E = AE_SLAI32(E, 1);
        AE_MULAFP32X2RAS(Y, Y, E);
        // 3-st iteration
        E = AE_MOVDA32X2(0x40000000, 0x40000000);
        AE_MULSFP32X2RAS(E, X, Y);
        E = AE_SLAI32(E, 1);
        AE_MULAFP32X2RAS(Y, Y, E);
        // apply sign and move right values if input is zero
        AE_MOVT32X2(Y, AE_NEG32(Y), sgn);
        AE_MOVT32X2(Y, AE_MOVDA32X2(0x7FFFFFFF, 0x7FFFFFFF), inf);

        mu_exp = AE_NSAZ32_L(mu);
        mu_v = AE_SLAA32(mu, mu_exp);

        s_frac = AE_MULFP32X2RAS(Y, mu_v);
        s_exp -= mu_exp;
        ASSERT(s_exp >= -31 && s_exp <= 31);
    }

    //
    // Calculate the cross-correlation between the error signal and the
    // reference signal. Scale the result and update the estimation of the
    // impulse response.
    //
    pH = (const ae_int32x2 *)(h + M - 2);
    pHw= (      ae_int32x2 *)pH;
    __Pragma("loop_count min=1");
    for (m = 0; m < (M >> 3); m++)
    {
        ae_f64 q0_re, q1_re, q2_re, q3_re;
        ae_f64 q0_im, q1_im, q2_im, q3_im;
        ae_f64 q4_re, q5_re, q6_re, q7_re;
        ae_f64 q4_im, q5_im, q6_im, q7_im;
        ae_int32x2 e0, h0, h1;
        ae_int32x2 x0, x1, x2, x3, x4, x5, x6, x7;
        ae_int64 xx;

        pX = (const ae_int64   *)(x + 8 * m);
        pE = (      ae_int32x2 *)(e);
        q0_re = q1_re = q2_re = q3_re = AE_ZERO64();
        q0_im = q1_im = q2_im = q3_im = AE_ZERO64();
        q4_re = q5_re = q6_re = q7_re = AE_ZERO64();
        q4_im = q5_im = q6_im = q7_im = AE_ZERO64();
        /* preload data from x */
        AE_L64_IP(xx, pX, sizeof(complex_fract32)); x0 = AE_MOVINT32X2_FROMINT64(xx);
        AE_L64_IP(xx, pX, sizeof(complex_fract32)); x1 = AE_MOVINT32X2_FROMINT64(xx);
        AE_L64_IP(xx, pX, sizeof(complex_fract32)); x2 = AE_MOVINT32X2_FROMINT64(xx);
        AE_L64_IP(xx, pX, sizeof(complex_fract32)); x3 = AE_MOVINT32X2_FROMINT64(xx);
        AE_L64_IP(xx, pX, sizeof(complex_fract32)); x4 = AE_MOVINT32X2_FROMINT64(xx);
        AE_L64_IP(xx, pX, sizeof(complex_fract32)); x5 = AE_MOVINT32X2_FROMINT64(xx);
        AE_L64_IP(xx, pX, sizeof(complex_fract32)); x6 = AE_MOVINT32X2_FROMINT64(xx);

        __Pragma("loop_count min=8, factor=8");
        for (n = 0; n < N; n++)
        {
            /* load data from x */
            AE_L64_IP(xx, pX, sizeof(complex_fract32)); x7 = AE_MOVINT32X2_FROMINT64(xx);
            /* load data from y */
            AE_L32X2_IP(e0, pE, sizeof(complex_fract32));
            /* compute correlation of 8 values */
            AE_MULAFC32RA(q0_im, q0_re, x0, e0);
            AE_MULAFC32RA(q1_im, q1_re, x1, e0);
            AE_MULAFC32RA(q2_im, q2_re, x2, e0);
            AE_MULAFC32RA(q3_im, q3_re, x3, e0);
            AE_MULAFC32RA(q4_im, q4_re, x4, e0);
            AE_MULAFC32RA(q5_im, q5_re, x5, e0);
            AE_MULAFC32RA(q6_im, q6_re, x6, e0);
            AE_MULAFC32RA(q7_im, q7_re, x7, e0);
            /* shift input line for the next iteration */
            x0 = x1; x1 = x2; x2 = x3;
            x3 = x4; x4 = x5; x5 = x6; x6 = x7;
        }
        /* save computed samples */
        x0 = AE_ROUND32X2F48SASYM(q0_re, q0_im);
        x1 = AE_ROUND32X2F48SASYM(q1_re, q1_im);
        AE_MULFP32X2S_HH_LL(q0_re, q0_im, x0, s_frac);
        AE_MULFP32X2S_HH_LL(q1_re, q1_im, x1, s_frac);
        x0 = AE_TRUNCA32X2F64S(q0_re, q0_im, s_exp + 1);
        x1 = AE_TRUNCA32X2F64S(q1_re, q1_im, s_exp + 1);
        AE_L32X2X2_IP(h1, h0, castxcc(ae_int32x4, pH), -2 * (int)sizeof(complex_fract32));
        h0 = AE_ADD32S(h0, x0);
        h1 = AE_ADD32S(h1, x1);
        AE_S32X2X2_IP(h1, h0, castxcc(ae_int32x4, pHw), -2 * (int)sizeof(complex_fract32));

        x0 = AE_ROUND32X2F48SASYM(q2_re, q2_im);
        x1 = AE_ROUND32X2F48SASYM(q3_re, q3_im);
        AE_MULFP32X2S_HH_LL(q0_re, q0_im, x0, s_frac);
        AE_MULFP32X2S_HH_LL(q1_re, q1_im, x1, s_frac);
        x0 = AE_TRUNCA32X2F64S(q0_re, q0_im, s_exp + 1);
        x1 = AE_TRUNCA32X2F64S(q1_re, q1_im, s_exp + 1);
        AE_L32X2X2_IP(h1, h0, castxcc(ae_int32x4, pH), -2 * (int)sizeof(complex_fract32));
        h0 = AE_ADD32S(h0, x0);
        h1 = AE_ADD32S(h1, x1);
        AE_S32X2X2_IP(h1, h0, castxcc(ae_int32x4, pHw), -2 * (int)sizeof(complex_fract32));

        x0 = AE_ROUND32X2F48SASYM(q4_re, q4_im);
        x1 = AE_ROUND32X2F48SASYM(q5_re, q5_im);
        AE_MULFP32X2S_HH_LL(q0_re, q0_im, x0, s_frac);
        AE_MULFP32X2S_HH_LL(q1_re, q1_im, x1, s_frac);
        x0 = AE_TRUNCA32X2F64S(q0_re, q0_im, s_exp + 1);
        x1 = AE_TRUNCA32X2F64S(q1_re, q1_im, s_exp + 1);
        AE_L32X2X2_IP(h1, h0, castxcc(ae_int32x4, pH), -2 * (int)sizeof(complex_fract32));
        h0 = AE_ADD32S(h0, x0);
        h1 = AE_ADD32S(h1, x1);
        AE_S32X2X2_IP(h1, h0, castxcc(ae_int32x4, pHw), -2 * (int)sizeof(complex_fract32));

        x0 = AE_ROUND32X2F48SASYM(q6_re, q6_im);
        x1 = AE_ROUND32X2F48SASYM(q7_re, q7_im);
        AE_MULFP32X2S_HH_LL(q0_re, q0_im, x0, s_frac);
        AE_MULFP32X2S_HH_LL(q1_re, q1_im, x1, s_frac);
        x0 = AE_TRUNCA32X2F64S(q0_re, q0_im, s_exp + 1);
        x1 = AE_TRUNCA32X2F64S(q1_re, q1_im, s_exp + 1);
        AE_L32X2X2_IP(h1, h0, castxcc(ae_int32x4, pH), -2 * (int)sizeof(complex_fract32));
        h0 = AE_ADD32S(h0, x0);
        h1 = AE_ADD32S(h1, x1);
        AE_S32X2X2_IP(h1, h0, castxcc(ae_int32x4, pHw), -2 * (int)sizeof(complex_fract32));
    }
} // cxfir_blms32x32()
