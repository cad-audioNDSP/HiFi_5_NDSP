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
    Blockwise Adaptive LMS Algorithm for Real Data, 16x16-bit
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

void fir_blms16x16( int16_t * restrict e,
                    int16_t * restrict h,
              const int16_t * restrict r,
              const int16_t * restrict x,
              int16_t norm, int16_t mu,
              int N, int M)
{
    const ae_int16x4 * pX;
    const ae_int16x4 * S0;
    const ae_int16x4 * S1;
    const ae_int16x8 * pR;
    const ae_int16x4 * pEr;
          ae_int16x8 * pEw;
    const ae_int64   * pHr;
          ae_int64   * restrict pHw;

    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_int16x4 x0, x1, x2;
    ae_int32x2 t0, t1;
    ae_int16x4 h0, h1;
    ae_int16x4 r0, r1;
    ae_int16x4 e0;

    ae_int32x2 s_frac;
    int s_exp;

    int m, n;

    NASSERT(e);
    NASSERT(h);
    NASSERT(r);
    NASSERT(x);
    NASSERT_ALIGN(e, 16);
    NASSERT_ALIGN(h, 16);
    NASSERT_ALIGN(r, 16);
    NASSERT_ALIGN(x, 16);
    NASSERT((N>0) && (M>0));
    NASSERT(((M%8) == 0) && ((N%8) == 0));

    //
    // Pass the reference signal x[] through the adaptive filter to obtain the
    // predicted signal and calculate the error, i.e. the distance to the
    // actual input signal r[].
    //
    pX = (const ae_int16x4 *)x;
    pR = (const ae_int16x8 *)r;
    pEw = (      ae_int16x8 *)e;
    for (n = 0; n < (N >> 4); n++)
    {
        ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
        ae_int16x4 x3, x4;
        ae_int64 t;

        pHr = (ae_int64 *)(h + M - 4);
        /* load in reversed order */
        AE_L64_IP(t, pHr, -4 * (int)sizeof(int16_t));
        h0 = AE_MOVINT16X4_FROMINT64(t);

        AE_L16X4_IP(x0, pX, 4 * sizeof(int16_t));
        S0 = pX;
        AE_L16X4_IP(x1, pX, 4 * sizeof(int16_t));
        AE_L16X4_IP(x2, pX, 4 * sizeof(int16_t));
        S1 = pX;
        AE_L16X4_IP(x3, pX, 4 * sizeof(int16_t));
        x4 = AE_L16X4_I(pX, 0);

        AE_MULFQ16X2_FIR_3(q0, q1, x0, x1, h0);
        AE_MULFQ16X2_FIR_1(q2, q3, x0, x1, h0);
        AE_MULFQ16X2_FIR_3(q4, q5, x1, x2, h0);
        AE_MULFQ16X2_FIR_1(q6, q7, x1, x2, h0);
        AE_MULFQ16X2_FIR_3(q8, q9, x2, x3, h0);
        AE_MULFQ16X2_FIR_1(qa, qb, x2, x3, h0);
        AE_MULFQ16X2_FIR_3(qc, qd, x3, x4, h0);
        AE_MULFQ16X2_FIR_1(qe, qf, x3, x4, h0);

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2) - 1; m++)
        {
            /* load in reversed order */
            AE_L64_IP(t, pHr, -4 * (int)sizeof(int16_t));
            h0 = AE_MOVINT16X4_FROMINT64(t);

            AE_L16X4_IP(x0, S0, 4 * sizeof(int16_t));
            x1 = AE_L16X4_I(S0, 0);
            AE_L16X4_IP(x2, S1, 4 * sizeof(int16_t));
            x3 = AE_L16X4_I(S1, 0);
            x4 = AE_L16X4_I(S1, 4 * sizeof(int16_t));

            AE_MULAFQ16X2_FIR_3(q0, q1, x0, x1, h0);
            AE_MULAFQ16X2_FIR_1(q2, q3, x0, x1, h0);
            AE_MULAFQ16X2_FIR_3(q4, q5, x1, x2, h0);
            AE_MULAFQ16X2_FIR_1(q6, q7, x1, x2, h0);
            AE_MULAFQ16X2_FIR_3(q8, q9, x2, x3, h0);
            AE_MULAFQ16X2_FIR_1(qa, qb, x2, x3, h0);
            AE_MULAFQ16X2_FIR_3(qc, qd, x3, x4, h0);
            AE_MULAFQ16X2_FIR_1(qe, qf, x3, x4, h0);
        }

        t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
        t1 = AE_SAT32X2(q2, q3);
        x0 = AE_ROUND16X4F32SASYM(t0, t1);
        t0 = AE_TRUNCA32X2F64S(q4, q5, 32);
        t1 = AE_SAT32X2(q6, q7);
        x1 = AE_ROUND16X4F32SASYM(t0, t1);

        AE_L16X4X2_IP(r0, r1, pR, 8 * sizeof(int16_t));
        x0 = AE_SUB16S(r0, x0);
        x1 = AE_SUB16S(r1, x1);
        AE_S16X4X2_IP(x0, x1, pEw, 8 * sizeof(int16_t));

        t0 = AE_TRUNCA32X2F64S(q8, q9, 32);
        t1 = AE_SAT32X2(qa, qb);
        x0 = AE_ROUND16X4F32SASYM(t0, t1);
        t0 = AE_TRUNCA32X2F64S(qc, qd, 32);
        t1 = AE_SAT32X2(qe, qf);
        x1 = AE_ROUND16X4F32SASYM(t0, t1);

        AE_L16X4X2_IP(r0, r1, pR, 8 * sizeof(int16_t));
        x0 = AE_SUB16S(r0, x0);
        x1 = AE_SUB16S(r1, x1);
        AE_S16X4X2_IP(x0, x1, pEw, 8 * sizeof(int16_t));
    }
    if (N & 8)
    {
        pHr = (ae_int64 *)(h + M - 4);
        S0 = pX;
        S1 = (const ae_int16x4 *)XT_ADDX2(4, (uintptr_t)pX);
        pX = (const ae_int16x4 *)XT_ADDX2(8, (uintptr_t)pX);
        q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

        __Pragma("loop_count min=2, factor=2");
        for (m = 0; m < (M >> 2); m++)
        {
            /* load in reversed order */
            ae_int64 t;
            AE_L64_IP(t, pHr, -4 * (int)sizeof(int16_t));
            h0 = AE_MOVINT16X4_FROMINT64(t);

            AE_L16X4_XP(x0, S0, 4 * sizeof(int16_t));
            AE_L16X4_XP(x1, S1, 4 * sizeof(int16_t));
            x2 = AE_L16X4_I(S1, 0);

            AE_MULAFQ16X2_FIR_3(q0, q1, x0, x1, h0);
            AE_MULAFQ16X2_FIR_1(q2, q3, x0, x1, h0);
            AE_MULAFQ16X2_FIR_3(q4, q5, x1, x2, h0);
            AE_MULAFQ16X2_FIR_1(q6, q7, x1, x2, h0);
        }

        t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
        t1 = AE_SAT32X2(q2, q3);
        x0 = AE_ROUND16X4F32SASYM(t0, t1);
        t0 = AE_TRUNCA32X2F64S(q4, q5, 32);
        t1 = AE_SAT32X2(q6, q7);
        x1 = AE_ROUND16X4F32SASYM(t0, t1);

        AE_L16X4X2_IP(r0, r1, pR, 8 * sizeof(int16_t));
        x0 = AE_SUB16S(r0, x0);
        x1 = AE_SUB16S(r1, x1);
        AE_S16X4X2_IP(x0, x1, pEw, 8 * sizeof(int16_t));
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

        X = ((int32_t)norm) << 16;
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

        mu_exp = AE_NSAZ16_0(mu);
        mu_v = AE_SLAA32(((int32_t)mu) << 16, mu_exp);

        s_frac = AE_MULFP32X2RAS(Y, mu_v);
        s_exp -= mu_exp;
        ASSERT(s_exp >= -15 && s_exp <= 31);
    }

    //
    // Calculate the cross-correlation between the error signal and the
    // reference signal. Scale the result and update the estimation of the
    // impulse response.
    //
    pX = (const ae_int16x4 *)x;
    pHw = (ae_int64 *)(h + M - 8);
    pHr = pHw;
    for (m = 0; m < (M >> 4); m++)
    {
        ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
        ae_int16x4 x3, x4;

        pEr = (ae_int16x4 *)(e);
        AE_L16X4_IP(e0, pEr, 4 * sizeof(int16_t));

        AE_L16X4_IP(x0, pX, 4 * sizeof(int16_t));
        S0 = pX;
        AE_L16X4_IP(x1, pX, 4 * sizeof(int16_t));
        AE_L16X4_IP(x2, pX, 4 * sizeof(int16_t));
        S1 = pX;
        AE_L16X4_IP(x3, pX, 4 * sizeof(int16_t));
        x4 = AE_L16X4_I(pX, 0);

        AE_MULFQ16X2_FIR_3(q0, q1, x0, x1, e0);
        AE_MULFQ16X2_FIR_1(q2, q3, x0, x1, e0);
        AE_MULFQ16X2_FIR_3(q4, q5, x1, x2, e0);
        AE_MULFQ16X2_FIR_1(q6, q7, x1, x2, e0);
        AE_MULFQ16X2_FIR_3(q8, q9, x2, x3, e0);
        AE_MULFQ16X2_FIR_1(qa, qb, x2, x3, e0);
        AE_MULFQ16X2_FIR_3(qc, qd, x3, x4, e0);
        AE_MULFQ16X2_FIR_1(qe, qf, x3, x4, e0);

        __Pragma("loop_count min=1");
        for (n = 0; n < (N >> 2) - 1; n++)
        {
            AE_L16X4_IP(e0, pEr, 4 * sizeof(int16_t));

            AE_L16X4_IP(x0, S0, 4 * sizeof(int16_t));
            x1 = AE_L16X4_I(S0, 0);
            AE_L16X4_IP(x2, S1, 4 * sizeof(int16_t));
            x3 = AE_L16X4_I(S1, 0);
            x4 = AE_L16X4_I(S1, 4 * sizeof(int16_t));

            AE_MULAFQ16X2_FIR_3(q0, q1, x0, x1, e0);
            AE_MULAFQ16X2_FIR_1(q2, q3, x0, x1, e0);
            AE_MULAFQ16X2_FIR_3(q4, q5, x1, x2, e0);
            AE_MULAFQ16X2_FIR_1(q6, q7, x1, x2, e0);
            AE_MULAFQ16X2_FIR_3(q8, q9, x2, x3, e0);
            AE_MULAFQ16X2_FIR_1(qa, qb, x2, x3, e0);
            AE_MULAFQ16X2_FIR_3(qc, qd, x3, x4, e0);
            AE_MULAFQ16X2_FIR_1(qe, qf, x3, x4, e0);
        }

        t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
        t1 = AE_SAT32X2(q2, q3);
        AE_MULFP32X2S_HH_LL(q0, q1, s_frac, t0);
        AE_MULFP32X2S_HH_LL(q2, q3, s_frac, t1);
        t0 = AE_TRUNCA32X2F64S(q0, q1, s_exp - 15);
        t1 = AE_TRUNCA32X2F64S(q2, q3, s_exp - 15);
        x0 = AE_SAT16X4(t0, t1);

        t0 = AE_TRUNCA32X2F64S(q4, q5, 32);
        t1 = AE_SAT32X2(q6, q7);
        AE_MULFP32X2S_HH_LL(q4, q5, s_frac, t0);
        AE_MULFP32X2S_HH_LL(q6, q7, s_frac, t1);
        t0 = AE_TRUNCA32X2F64S(q4, q5, s_exp - 15);
        t1 = AE_TRUNCA32X2F64S(q6, q7, s_exp - 15);
        x1 = AE_SAT16X4(t0, t1);

        t0 = AE_TRUNCA32X2F64S(q8, q9, 32);
        t1 = AE_SAT32X2(qa, qb);
        AE_MULFP32X2S_HH_LL(q8, q9, s_frac, t0);
        AE_MULFP32X2S_HH_LL(qa, qb, s_frac, t1);
        t0 = AE_TRUNCA32X2F64S(q8, q9, s_exp - 15);
        t1 = AE_TRUNCA32X2F64S(qa, qb, s_exp - 15);
        x2 = AE_SAT16X4(t0, t1);

        t0 = AE_TRUNCA32X2F64S(qc, qd, 32);
        t1 = AE_SAT32X2(qe, qf);
        AE_MULFP32X2S_HH_LL(qc, qd, s_frac, t0);
        AE_MULFP32X2S_HH_LL(qe, qf, s_frac, t1);
        t0 = AE_TRUNCA32X2F64S(qc, qd, s_exp - 15);
        t1 = AE_TRUNCA32X2F64S(qe, qf, s_exp - 15);
        x3 = AE_SAT16X4(t0, t1);

        {
            ae_int64 th0, th1;
            /* load in reversed order */
            AE_L64X2_IP(th1, th0, castxcc(ae_int64x2, pHr), -8 * (int)sizeof(int16_t));
            h0 = AE_MOVINT16X4_FROMINT64(th0);
            h1 = AE_MOVINT16X4_FROMINT64(th1);
            h0 = AE_ADD16S(x0, h0);
            h1 = AE_ADD16S(x1, h1);
            th0 = AE_MOVINT64_FROMINT16X4(h0);
            th1 = AE_MOVINT64_FROMINT16X4(h1);
            AE_S64X2_IP(th1, th0, castxcc(ae_int64x2, pHw), -8 * (int)sizeof(int16_t));
            /* load in reversed order */
            AE_L64X2_IP(th1, th0, castxcc(ae_int64x2, pHr), -8 * (int)sizeof(int16_t));
            h0 = AE_MOVINT16X4_FROMINT64(th0);
            h1 = AE_MOVINT16X4_FROMINT64(th1);
            h0 = AE_ADD16S(x2, h0);
            h1 = AE_ADD16S(x3, h1);
            th0 = AE_MOVINT64_FROMINT16X4(h0);
            th1 = AE_MOVINT64_FROMINT16X4(h1);
            AE_S64X2_IP(th1, th0, castxcc(ae_int64x2, pHw), -8 * (int)sizeof(int16_t));
        }
    }
    if (M & 8)
    {
        pEr = (ae_int16x4 *)(e);
        S0 = pX;
        S1 = (const ae_int16x4 *)XT_ADDX2(4, (uintptr_t)pX);
        pX = (const ae_int16x4 *)XT_ADDX2(8, (uintptr_t)pX);
        q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

        __Pragma("loop_count min=2, factor=2");
        for (n = 0; n < (N >> 2); n++)
        {
            AE_L16X4_IP(e0, pEr, 4 * sizeof(int16_t));

            AE_L16X4_XP(x0, S0, 4 * sizeof(int16_t));
            AE_L16X4_XP(x1, S1, 4 * sizeof(int16_t));
            x2 = AE_L16X4_I(S1, 0);

            AE_MULAFQ16X2_FIR_3(q0, q1, x0, x1, e0);
            AE_MULAFQ16X2_FIR_1(q2, q3, x0, x1, e0);
            AE_MULAFQ16X2_FIR_3(q4, q5, x1, x2, e0);
            AE_MULAFQ16X2_FIR_1(q6, q7, x1, x2, e0);
        }

        t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
        t1 = AE_SAT32X2(q2, q3);
        AE_MULFP32X2S_HH_LL(q0, q1, s_frac, t0);
        AE_MULFP32X2S_HH_LL(q2, q3, s_frac, t1);
        t0 = AE_TRUNCA32X2F64S(q0, q1, s_exp - 15);
        t1 = AE_TRUNCA32X2F64S(q2, q3, s_exp - 15);
        x0 = AE_SAT16X4(t0, t1);

        t0 = AE_TRUNCA32X2F64S(q4, q5, 32);
        t1 = AE_SAT32X2(q6, q7);
        AE_MULFP32X2S_HH_LL(q4, q5, s_frac, t0);
        AE_MULFP32X2S_HH_LL(q6, q7, s_frac, t1);
        t0 = AE_TRUNCA32X2F64S(q4, q5, s_exp - 15);
        t1 = AE_TRUNCA32X2F64S(q6, q7, s_exp - 15);
        x1 = AE_SAT16X4(t0, t1);

        {
            ae_int64 th0, th1;
            /* load in reversed order */
            AE_L64X2_IP(th1, th0, castxcc(ae_int64x2, pHr), -8 * (int)sizeof(int16_t));
            h0 = AE_MOVINT16X4_FROMINT64(th0);
            h1 = AE_MOVINT16X4_FROMINT64(th1);
            h0 = AE_ADD16S(x0, h0);
            h1 = AE_ADD16S(x1, h1);
            th0 = AE_MOVINT64_FROMINT16X4(h0);
            th1 = AE_MOVINT64_FROMINT16X4(h1);
            AE_S64X2_IP(th1, th0, castxcc(ae_int64x2, pHw), -8 * (int)sizeof(int16_t));
        }
    }
} /* fir_blms16x16() */
