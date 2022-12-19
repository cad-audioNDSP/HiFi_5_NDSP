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
    Blockwise Adaptive LMS Algorithm for Real Data, 16x32-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

#define AE_L32X2X2_RIP(h0, h1, pH, ofs)                    \
{                                                          \
    ae_int64 tmp0, tmp1;                                   \
    AE_L64X2_IP(tmp1, tmp0, castxcc(ae_int64x2, pH), ofs); \
    h0 = AE_MOVINT32X2_FROMINT64(tmp0);                    \
    h1 = AE_MOVINT32X2_FROMINT64(tmp1);                    \
}

#define AE_S32X2X2_RIP(h0, h1, pH, ofs)                    \
{                                                          \
    ae_int64 tmp0, tmp1;                                   \
    tmp0 = AE_MOVINT64_FROMINT32X2(h0);                    \
    tmp1 = AE_MOVINT64_FROMINT32X2(h1);                    \
    AE_S64X2_IP(tmp1, tmp0, castxcc(ae_int64x2, pH), ofs); \
}

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

void fir_blms16x32( int32_t * restrict e,
                    int32_t * restrict h,
              const int16_t * restrict r,
              const int16_t * restrict x,
              int32_t norm, int16_t mu,
              int N, int M )
{
    const ae_int32x4 *          pH;
          ae_int32x4 *          pHw;
    const ae_int32x4 *          S0;
    const ae_int32x4 *          S1;
    const ae_int16x4 *          pX;
    const ae_int16x4 *          pX0;
    const ae_int16x4 *          pX1;
    const ae_int16x8 *          pR;
          ae_int32x4 * restrict pE;

    ae_int32x2 s_frac;  // scale factor (fraction part)
    int s_exp;          // scale factor (exponent)
    int m, n;

#if 1
    static const ALIGN(16) int8_t Sel[8] = { 1, 3, 5, 7, 0, 2, 4, 6 };
    ae_int16x4 sel;
    sel = AE_L16X4_I((ae_int16x4*)&Sel, 0);
#endif

    NASSERT(e && h && r && x);
    NASSERT_ALIGN(e, 16);
    NASSERT_ALIGN(h, 16);
    NASSERT_ALIGN(r, 16);
    NASSERT_ALIGN(x, 16);
    NASSERT(mu >= 0 && norm > 0);
    NASSERT(N >= 0 && N % 8 == 0);
    NASSERT(M >= 0 && M % 8 == 0);

    //
    // Pass the reference signal x[] through the adaptive filter to obtain the
    // predicted signal and calculate the error, i.e. the distance to the
    // actual input signal r[].
    //
    pH = (const ae_int32x4 *)(h + M - 4);
    pX = (const ae_int16x4 *)x;
    pR = (const ae_int16x8 *)r;
    pE = (      ae_int32x4 *)e;

    for (n = 0; n < (N >> 4); n++)
    {
        ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
        ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
        ae_int32x2 h0, h1, h2, h3;
        ae_int16x4 x0, x1, x2, x3;
        ae_int32x2 e0, e1, e2, e3;

        S0 = pH;
        S1 = pH;
        pX0 = pX;
        pX1 = pX + 2;
        pX += 4;

        AE_L16X4X2_X(x0, x1, (ae_int16x8 *)pX0, M * sizeof(int16_t));
        AE_L16X4X2_X(x2, x3, (ae_int16x8 *)pX1, M * sizeof(int16_t));
        AE_L32X2X2_RIP(h0, h1, h, 0);
        h2 = h3 = AE_ZERO32();
        AE_MUL2Q32X16_FIR_L(q3, q2, h0, h1, h2, x0);
        AE_MUL2Q32X16_FIR_L(q1, q0, h1, h2, h3, x0);
        AE_MUL2Q32X16_FIR_L(q7, q6, h0, h1, h2, x1);
        AE_MUL2Q32X16_FIR_L(q5, q4, h1, h2, h3, x1);
        AE_MUL2Q32X16_FIR_L(qb, qa, h0, h1, h2, x2);
        AE_MUL2Q32X16_FIR_L(q9, q8, h1, h2, h3, x2);
        AE_MUL2Q32X16_FIR_L(qf, qe, h0, h1, h2, x3);
        AE_MUL2Q32X16_FIR_L(qd, qc, h1, h2, h3, x3);
        h0 = h1 = AE_ZERO32();

        __Pragma("loop_count min=2, factor=2");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L16X4_IP(x0, pX0, 4 * sizeof(int16_t));
            x1 = AE_L16X4_X(pX0, 0 * sizeof(int16_t));
            AE_L16X4_IP(x2, pX1, 4 * sizeof(int16_t));
            x3 = AE_L16X4_X(pX1, 0 * sizeof(int16_t));
            AE_L32X2X2_RIP(h2, h3, S0, -4 * (int)sizeof(int32_t));
            AE_MULA2Q32X16_FIR_L(q3, q2, h0, h1, h2, x0);
            AE_MULA2Q32X16_FIR_L(q1, q0, h1, h2, h3, x0);
            AE_MULA2Q32X16_FIR_L(q7, q6, h0, h1, h2, x1);
            AE_MULA2Q32X16_FIR_L(q5, q4, h1, h2, h3, x1);
            AE_MULA2Q32X16_FIR_L(qb, qa, h0, h1, h2, x2);
            AE_MULA2Q32X16_FIR_L(q9, q8, h1, h2, h3, x2);
            AE_MULA2Q32X16_FIR_L(qf, qe, h0, h1, h2, x3);
            AE_MULA2Q32X16_FIR_L(qd, qc, h1, h2, h3, x3);
            AE_L32X2X2_RIP(h0, h1, S1, -4 * (int)sizeof(int32_t));
        }

        AE_PKSR32(e0, q0, 1); AE_PKSR32(e0, q1, 1);
        AE_PKSR32(e1, q2, 1); AE_PKSR32(e1, q3, 1);
        AE_PKSR32(e2, q4, 1); AE_PKSR32(e2, q5, 1);
        AE_PKSR32(e3, q6, 1); AE_PKSR32(e3, q7, 1);

        AE_L16X4X2_IP(x0, x1, pR, 8 * sizeof(int16_t));
#if 0
        h0 = AE_CVT32X2F16_32(x0);
        h1 = AE_CVT32X2F16_10(x0);
        h2 = AE_CVT32X2F16_32(x1);
        h3 = AE_CVT32X2F16_10(x1);
#else
        AE_DSEL16X4(x2, x3, AE_ZERO16(), x0, sel);
        h0 = AE_MOVINT32X2_FROMINT16X4(x2);
        h1 = AE_MOVINT32X2_FROMINT16X4(x3);
        AE_DSEL16X4(x2, x3, AE_ZERO16(), x1, sel);
        h2 = AE_MOVINT32X2_FROMINT16X4(x2);
        h3 = AE_MOVINT32X2_FROMINT16X4(x3);
#endif
        e0 = AE_SUB32S(h0, e0);
        e1 = AE_SUB32S(h1, e1);
        e2 = AE_SUB32S(h2, e2);
        e3 = AE_SUB32S(h3, e3);
        AE_S32X2X2_IP(e0, e1, pE, 4 * sizeof(int32_t));
        AE_S32X2X2_IP(e2, e3, pE, 4 * sizeof(int32_t));

        AE_PKSR32(e0, q8, 1); AE_PKSR32(e0, q9, 1);
        AE_PKSR32(e1, qa, 1); AE_PKSR32(e1, qb, 1);
        AE_PKSR32(e2, qc, 1); AE_PKSR32(e2, qd, 1);
        AE_PKSR32(e3, qe, 1); AE_PKSR32(e3, qf, 1);

        AE_L16X4X2_IP(x0, x1, pR, 8 * sizeof(int16_t));
#if 0
        h0 = AE_CVT32X2F16_32(x0);
        h1 = AE_CVT32X2F16_10(x0);
        h2 = AE_CVT32X2F16_32(x1);
        h3 = AE_CVT32X2F16_10(x1);
#else
        AE_DSEL16X4(x2, x3, AE_ZERO16(), x0, sel);
        h0 = AE_MOVINT32X2_FROMINT16X4(x2);
        h1 = AE_MOVINT32X2_FROMINT16X4(x3);
        AE_DSEL16X4(x2, x3, AE_ZERO16(), x1, sel);
        h2 = AE_MOVINT32X2_FROMINT16X4(x2);
        h3 = AE_MOVINT32X2_FROMINT16X4(x3);
#endif
        e0 = AE_SUB32S(h0, e0);
        e1 = AE_SUB32S(h1, e1);
        e2 = AE_SUB32S(h2, e2);
        e3 = AE_SUB32S(h3, e3);
        AE_S32X2X2_IP(e0, e1, pE, 4 * sizeof(int32_t));
        AE_S32X2X2_IP(e2, e3, pE, 4 * sizeof(int32_t));
    }
    if (N & 8)
    {
        ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
        ae_int32x2 h0, h1, h2, h3;
        ae_int16x4 x0, x1;
        ae_int16x4 x2, x3;
        ae_int32x2 e0, e1, e2, e3;

        S0 = pH;
        S1 = pH;
        pX0 = pX;
        pX += 2;

        AE_L16X4X2_X(x0, x1, (ae_int16x8 *)pX0, M * sizeof(int16_t));
        AE_L32X2X2_RIP(h0, h1, h, 0);
        h2 = h3 = AE_ZERO32();
        AE_MUL2Q32X16_FIR_L(q3, q2, h0, h1, h2, x0);
        AE_MUL2Q32X16_FIR_L(q1, q0, h1, h2, h3, x0);
        AE_MUL2Q32X16_FIR_L(q7, q6, h0, h1, h2, x1);
        AE_MUL2Q32X16_FIR_L(q5, q4, h1, h2, h3, x1);
        h0 = h1 = AE_ZERO32();

        __Pragma("loop_count min=2, factor=2");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L16X4_IP(x0, pX0, 4 * sizeof(int16_t));
            x1 = AE_L16X4_X(pX0, 0);
            AE_L32X2X2_RIP(h2, h3, S0, -4 * (int)sizeof(int32_t));
            AE_MULA2Q32X16_FIR_L(q3, q2, h0, h1, h2, x0);
            AE_MULA2Q32X16_FIR_L(q1, q0, h1, h2, h3, x0);
            AE_MULA2Q32X16_FIR_L(q7, q6, h0, h1, h2, x1);
            AE_MULA2Q32X16_FIR_L(q5, q4, h1, h2, h3, x1);
            AE_L32X2X2_RIP(h0, h1, S1, -4 * (int)sizeof(int32_t));
        }

#if 0
        q0 = AE_SLAI64(q0, 1); q1 = AE_SLAI64(q1, 1);
        e0 = AE_ROUND32X2F48SASYM(q0, q1);
        AE_PKSR32(e1, q2, 1); AE_PKSR32(e1, q3, 1);
        q4 = AE_SLAI64(q4, 1); q5 = AE_SLAI64(q5, 1);
        e2 = AE_ROUND32X2F48SASYM(q4, q5);
        AE_PKSR32(e3, q6, 1); AE_PKSR32(e3, q7, 1);
#elif 1
        AE_PKSR32(e0, q0, 1); AE_PKSR32(e0, q1, 1);
        AE_PKSR32(e1, q2, 1); AE_PKSR32(e1, q3, 1);
        AE_PKSR32(e2, q4, 1); AE_PKSR32(e2, q5, 1);
        AE_PKSR32(e3, q6, 1); AE_PKSR32(e3, q7, 1);
#else
        q0 = AE_SLAI64(q0, 1); q1 = AE_SLAI64(q1, 1);
        e0 = AE_ROUND32X2F48SASYM(q0, q1);
        q2 = AE_SLAI64(q2, 1); q3 = AE_SLAI64(q3, 1);
        e1 = AE_ROUND32X2F48SASYM(q2, q3);
        q4 = AE_SLAI64(q4, 1); q5 = AE_SLAI64(q5, 1);
        e2 = AE_ROUND32X2F48SASYM(q4, q5);
        q6 = AE_SLAI64(q6, 1); q7 = AE_SLAI64(q7, 1);
        e3 = AE_ROUND32X2F48SASYM(q6, q7);
#endif

        AE_L16X4X2_IP(x0, x1, pR, 8 * sizeof(int16_t));
#if 0
        h0 = AE_CVT32X2F16_32(x0);
        h1 = AE_CVT32X2F16_10(x0);
        h2 = AE_CVT32X2F16_32(x1);
        h3 = AE_CVT32X2F16_10(x1);
#else
        AE_DSEL16X4(x2, x3, AE_ZERO16(), x0, sel);
        h0 = AE_MOVINT32X2_FROMINT16X4(x2);
        h1 = AE_MOVINT32X2_FROMINT16X4(x3);
        AE_DSEL16X4(x2, x3, AE_ZERO16(), x1, sel);
        h2 = AE_MOVINT32X2_FROMINT16X4(x2);
        h3 = AE_MOVINT32X2_FROMINT16X4(x3);
#endif
        e0 = AE_SUB32S(h0, e0);
        e1 = AE_SUB32S(h1, e1);
        e2 = AE_SUB32S(h2, e2);
        e3 = AE_SUB32S(h3, e3);
        AE_S32X2X2_IP(e0, e1, pE, 4 * sizeof(int32_t));
        AE_S32X2X2_IP(e2, e3, pE, 4 * sizeof(int32_t));
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
    pE = (      ae_int32x4 *)e;
    pH = (const ae_int32x4 *)(h + M - 4);
    pHw= (      ae_int32x4 *)pH;

    for (m = 0; m < (M >> 4); m++)
    {
        ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
        ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
        ae_int32x2 h0, h1;
        ae_int16x4 x0, x1, x2, x3;
        ae_int32x2 e0, e1, e2, e3;

        S0 = pE;
        S1 = pE;
        pX0 = pX;
        pX1 = pX + 2;
        pX += 4;

        AE_L16X4X2_X(x0, x1, (ae_int16x8 *)pX0, N * sizeof(int16_t));
        AE_L16X4X2_X(x2, x3, (ae_int16x8 *)pX1, N * sizeof(int16_t));
        AE_L32X2X2_X(e0, e1, pE, (N - 4) * sizeof(int32_t));
        e2 = e3 = AE_ZERO32();
        AE_MUL2Q32X16_FIR_L(q3, q2, e0, e1, e2, x0);
        AE_MUL2Q32X16_FIR_L(q1, q0, e1, e2, e3, x0);
        AE_MUL2Q32X16_FIR_L(q7, q6, e0, e1, e2, x1);
        AE_MUL2Q32X16_FIR_L(q5, q4, e1, e2, e3, x1);
        AE_MUL2Q32X16_FIR_L(qb, qa, e0, e1, e2, x2);
        AE_MUL2Q32X16_FIR_L(q9, q8, e1, e2, e3, x2);
        AE_MUL2Q32X16_FIR_L(qf, qe, e0, e1, e2, x3);
        AE_MUL2Q32X16_FIR_L(qd, qc, e1, e2, e3, x3);
        e0 = e1 = AE_ZERO32();

        __Pragma("loop_count min=2, factor=2");
        for (n = 0; n < (N >> 2); n++)
        {
            AE_L16X4_IP(x0, pX0, 4 * sizeof(int16_t));
            x1 = AE_L16X4_X(pX0, 0 * sizeof(int16_t));
            AE_L16X4_IP(x2, pX1, 4 * sizeof(int16_t));
            x3 = AE_L16X4_X(pX1, 0 * sizeof(int16_t));
            AE_L32X2X2_IP(e2, e3, S0, 4 * sizeof(int32_t));
            AE_MULA2Q32X16_FIR_L(q3, q2, e0, e1, e2, x0);
            AE_MULA2Q32X16_FIR_L(q1, q0, e1, e2, e3, x0);
            AE_MULA2Q32X16_FIR_L(q7, q6, e0, e1, e2, x1);
            AE_MULA2Q32X16_FIR_L(q5, q4, e1, e2, e3, x1);
            AE_MULA2Q32X16_FIR_L(qb, qa, e0, e1, e2, x2);
            AE_MULA2Q32X16_FIR_L(q9, q8, e1, e2, e3, x2);
            AE_MULA2Q32X16_FIR_L(qf, qe, e0, e1, e2, x3);
            AE_MULA2Q32X16_FIR_L(qd, qc, e1, e2, e3, x3);
            AE_L32X2X2_IP(e0, e1, S1, 4 * sizeof(int32_t));
        }

        AE_PKSR32(e0, q0, 1); AE_PKSR32(e0, q1, 1);
        AE_PKSR32(e1, q2, 1); AE_PKSR32(e1, q3, 1);
        AE_PKSR32(e2, q4, 1); AE_PKSR32(e2, q5, 1);
        AE_PKSR32(e3, q6, 1); AE_PKSR32(e3, q7, 1);

        AE_MULFP32X2S_HH_LL(q0, q1, e0, s_frac);
        AE_MULFP32X2S_HH_LL(q2, q3, e1, s_frac);
        e0 = AE_TRUNCA32X2F64S(q0, q1, s_exp + 1);
        e1 = AE_TRUNCA32X2F64S(q2, q3, s_exp + 1);
        AE_L32X2X2_RIP(h0, h1, pH, -4 * (int)sizeof(int32_t));
        h0 = AE_ADD32S(h0, e0);
        h1 = AE_ADD32S(h1, e1);
        AE_S32X2X2_RIP(h0, h1, pHw, -4 * (int)sizeof(int32_t));

        AE_MULFP32X2S_HH_LL(q0, q1, e2, s_frac);
        AE_MULFP32X2S_HH_LL(q2, q3, e3, s_frac);
        e0 = AE_TRUNCA32X2F64S(q0, q1, s_exp + 1);
        e1 = AE_TRUNCA32X2F64S(q2, q3, s_exp + 1);
        AE_L32X2X2_RIP(h0, h1, pH, -4 * (int)sizeof(int32_t));
        h0 = AE_ADD32S(h0, e0);
        h1 = AE_ADD32S(h1, e1);
        AE_S32X2X2_RIP(h0, h1, pHw, -4 * (int)sizeof(int32_t));

        AE_PKSR32(e0, q8, 1); AE_PKSR32(e0, q9, 1);
        AE_PKSR32(e1, qa, 1); AE_PKSR32(e1, qb, 1);
        AE_PKSR32(e2, qc, 1); AE_PKSR32(e2, qd, 1);
        AE_PKSR32(e3, qe, 1); AE_PKSR32(e3, qf, 1);

        AE_MULFP32X2S_HH_LL(q0, q1, e0, s_frac);
        AE_MULFP32X2S_HH_LL(q2, q3, e1, s_frac);
        e0 = AE_TRUNCA32X2F64S(q0, q1, s_exp + 1);
        e1 = AE_TRUNCA32X2F64S(q2, q3, s_exp + 1);
        AE_L32X2X2_RIP(h0, h1, pH, -4 * (int)sizeof(int32_t));
        h0 = AE_ADD32S(h0, e0);
        h1 = AE_ADD32S(h1, e1);
        AE_S32X2X2_RIP(h0, h1, pHw, -4 * (int)sizeof(int32_t));

        AE_MULFP32X2S_HH_LL(q0, q1, e2, s_frac);
        AE_MULFP32X2S_HH_LL(q2, q3, e3, s_frac);
        e0 = AE_TRUNCA32X2F64S(q0, q1, s_exp + 1);
        e1 = AE_TRUNCA32X2F64S(q2, q3, s_exp + 1);
        AE_L32X2X2_RIP(h0, h1, pH, -4 * (int)sizeof(int32_t));
        h0 = AE_ADD32S(h0, e0);
        h1 = AE_ADD32S(h1, e1);
        AE_S32X2X2_RIP(h0, h1, pHw, -4 * (int)sizeof(int32_t));
    }
    if (M & 8)
    {
        ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
        ae_int32x2 h0, h1;
        ae_int16x4 x0, x1;
        ae_int32x2 e0, e1, e2, e3;

        S0 = pE;
        S1 = pE;
        pX0 = pX;
        pX += 2;

        AE_L16X4X2_X(x0, x1, (ae_int16x8 *)pX0, N * sizeof(int16_t));
        AE_L32X2X2_X(e0, e1, pE, (N - 4) * sizeof(int32_t));
        e2 = e3 = AE_ZERO32();
        AE_MUL2Q32X16_FIR_L(q3, q2, e0, e1, e2, x0);
        AE_MUL2Q32X16_FIR_L(q1, q0, e1, e2, e3, x0);
        AE_MUL2Q32X16_FIR_L(q7, q6, e0, e1, e2, x1);
        AE_MUL2Q32X16_FIR_L(q5, q4, e1, e2, e3, x1);
        e0 = e1 = AE_ZERO32();

        __Pragma("loop_count min=2, factor=2");
        for (n = 0; n < (N >> 2); n++)
        {
            AE_L16X4_IP(x0, pX0, 4 * sizeof(int16_t));
            x1 = AE_L16X4_X(pX0, 0);
            AE_L32X2X2_IP(e2, e3, S0, 4 * sizeof(int32_t));
            AE_MULA2Q32X16_FIR_L(q3, q2, e0, e1, e2, x0);
            AE_MULA2Q32X16_FIR_L(q1, q0, e1, e2, e3, x0);
            AE_MULA2Q32X16_FIR_L(q7, q6, e0, e1, e2, x1);
            AE_MULA2Q32X16_FIR_L(q5, q4, e1, e2, e3, x1);
            AE_L32X2X2_IP(e0, e1, S1, 4 * sizeof(int32_t));
        }

        AE_PKSR32(e0, q0, 1); AE_PKSR32(e0, q1, 1);
        AE_PKSR32(e1, q2, 1); AE_PKSR32(e1, q3, 1);
        AE_PKSR32(e2, q4, 1); AE_PKSR32(e2, q5, 1);
        AE_PKSR32(e3, q6, 1); AE_PKSR32(e3, q7, 1);

        AE_MULFP32X2S_HH_LL(q0, q1, e0, s_frac);
        AE_MULFP32X2S_HH_LL(q2, q3, e1, s_frac);
        e0 = AE_TRUNCA32X2F64S(q0, q1, s_exp + 1);
        e1 = AE_TRUNCA32X2F64S(q2, q3, s_exp + 1);
        AE_L32X2X2_RIP(h0, h1, pH, -4 * (int)sizeof(int32_t));
        h0 = AE_ADD32S(h0, e0);
        h1 = AE_ADD32S(h1, e1);
        AE_S32X2X2_RIP(h0, h1, pHw, -4 * (int)sizeof(int32_t));

        AE_MULFP32X2S_HH_LL(q0, q1, e2, s_frac);
        AE_MULFP32X2S_HH_LL(q2, q3, e3, s_frac);
        e0 = AE_TRUNCA32X2F64S(q0, q1, s_exp + 1);
        e1 = AE_TRUNCA32X2F64S(q2, q3, s_exp + 1);
        AE_L32X2X2_RIP(h0, h1, pH, -4 * (int)sizeof(int32_t));
        h0 = AE_ADD32S(h0, e0);
        h1 = AE_ADD32S(h1, e1);
        AE_S32X2X2_RIP(h0, h1, pHw, -4 * (int)sizeof(int32_t));
    }
}
