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
    Blockwise Adaptive LMS Algorithm for Real Data, floating point
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common helper macros. */
#include "common.h"
#include "common_fpu.h"

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void,fir_blmsf,( float32_t * e, float32_t * h, const float32_t * r,
                const float32_t * x, float32_t norm, float32_t mu, int N, int M ))
#elif (HAVE_VFPU)

#define AE_DSELSX2(out0, out1, in0, in1, sel_mask)                       \
{                                                                        \
    ae_int16x4 tmp0, tmp1;                                               \
    tmp0 = AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(in0));  \
    tmp1 = AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(in1));  \
    AE_DSEL16X4(tmp0, tmp1, tmp0, tmp1, sel_mask);                       \
    out0 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(tmp0)); \
    out1 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(tmp1)); \
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

void fir_blmsf( float32_t * e,
                float32_t * h,
          const float32_t * r,
          const float32_t * x,
          float32_t norm, float32_t mu,
          int N, int M )
{
    const xtfloatx4 * restrict pX0;
    const xtfloatx4 * restrict pX1;
    const xtfloatx4 * restrict pX2;
    const xtfloatx4 * restrict pX3;
    const xtfloatx4 * restrict pX4;
    const xtfloatx4 * restrict pH;
    const xtfloatx4 * restrict pR;
          xtfloatx4 * restrict pE;
          xtfloatx4 * restrict pHw;
    xtfloatx2 r00, r10, r20, r30, r40, r50, r60, r70;
    xtfloatx2 r01, r11, r21, r31, r41, r51, r61, r71;
    xtfloatx2 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9;
    xtfloatx2 h0, h1, e0, e1;
    xtfloatx2 r0, r1, r2, r3;
    ae_int64 t0, t1;
    ae_valignx2 aX3, aX4;
    xtfloatx2 b;
    int m, n;

    static const ALIGN(16) int16_t Sel[4] = { 0x0705, 0x0604, 0x0301, 0x0200 };
    ae_int16x4 sel;
    sel = AE_L16X4_I((ae_int16x4*)&Sel, 0);

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

    /* estimate error */
    pR = (const xtfloatx4*)r;
    pE = (      xtfloatx4*)e;
    __Pragma("loop_count min=1");
    for (n = 0; n < N; n += 8)
    {
        pH  = (const xtfloatx4 *)(h + M - 4);
        pX0 = (const xtfloatx4 *)(x + n);
        pX1 = (const xtfloatx4 *)(x + n + 4);
        pX2 = (const xtfloatx4 *)(x + n + 8);
        pX3 = (const xtfloatx4 *)(x + n + 3);
        pX4 = (const xtfloatx4 *)(x + n + 7);
        aX3 = AE_LA128_PP(pX3);
        aX4 = AE_LA128_PP(pX4);

        CONST_SX2X2(r00, r01, 0); CONST_SX2X2(r10, r11, 0);
        CONST_SX2X2(r20, r21, 0); CONST_SX2X2(r30, r31, 0);
        CONST_SX2X2(r40, r41, 0); CONST_SX2X2(r50, r51, 0);
        CONST_SX2X2(r60, r61, 0); CONST_SX2X2(r70, r71, 0);

        /* preload data from x */
        AE_LSX2X2_IP(x0, x2, pX0, 4 * sizeof(float32_t));
        x1 = XT_SEL32_LH_SX2(x0, x2);

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            /* load data from x */
            AE_LSX2X2_IP(x4, x6, pX1, 4 * sizeof(float32_t));
            AE_LSX2IP(x8, castxcc(xtfloatx2, pX2), 4 * sizeof(float32_t));
            AE_LASX2X2_IP(x3, x5, aX3, pX3);
            AE_LASX2X2_IP(x7, x9, aX4, pX4);
            /* load data from y */
            AE_L64X2_IP(t1, t0, castxcc(ae_int64x2, pH), -4 * (int)sizeof(float32_t));
            h0 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT64(t0));
            h1 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT64(t1));
            /* compute correlation of 8 values */
            MADD_SX2X2(r00, r10, x0, x1, h0, h0);
            MADD_SX2X2(r20, r30, x2, x3, h0, h0);
            MADD_SX2X2(r40, r50, x4, x5, h0, h0);
            MADD_SX2X2(r60, r70, x6, x7, h0, h0);
            MADD_SX2X2(r01, r11, x2, x3, h1, h1);
            MADD_SX2X2(r21, r31, x4, x5, h1, h1);
            MADD_SX2X2(r41, r51, x6, x7, h1, h1);
            MADD_SX2X2(r61, r71, x8, x9, h1, h1);
            /* reload data and shift input line for the next iteration */
            AE_LSX2X2_IP(x0, x2, pX0, 4 * sizeof(float32_t));
            x1 = x5;
        }
        /* reduction add */
        ADD_SX2X2(r00, r10, r00, r10, r01, r11);
        ADD_SX2X2(r20, r30, r20, r30, r21, r31);
        AE_DSELSX2(r0, r1, r00, r10, sel);
        AE_DSELSX2(r2, r3, r20, r30, sel);
        ADD_SX2X2(e0, e1, r0, r2, r1, r3);
        /* compute and save error */
        AE_LSX2X2_IP(r0, r1, pR, 4 * sizeof(float32_t));
        SUB_SX2X2(e0, e1, r0, r1, e0, e1);
        AE_SSX2X2_IP(e0, e1, pE, 4 * sizeof(float32_t));

        /* reduction add */
        ADD_SX2X2(r40, r50, r40, r50, r41, r51);
        ADD_SX2X2(r60, r70, r60, r70, r61, r71);
        AE_DSELSX2(r0, r1, r40, r50, sel);
        AE_DSELSX2(r2, r3, r60, r70, sel);
        ADD_SX2X2(e0, e1, r0, r2, r1, r3);
        /* compute and save error */
        AE_LSX2X2_IP(r0, r1, pR, 4 * sizeof(float32_t));
        SUB_SX2X2(e0, e1, r0, r1, e0, e1);
        AE_SSX2X2_IP(e0, e1, pE, 4 * sizeof(float32_t));
    }

    /* update impluse response */
    b = XT_MUL_SX2(mu, XT_RECIP_SX2(norm));
    pHw = (xtfloatx4*)(h + M - 4);
    pH = pHw;
    __Pragma("loop_count min=1");
    for (m = 0; m < M; m += 8)
    {
        pE = (xtfloatx4 *)e;
        pX0 = (const xtfloatx4 *)(x + m);
        pX1 = (const xtfloatx4 *)(x + m + 4);
        pX2 = (const xtfloatx4 *)(x + m + 8);
        pX3 = (const xtfloatx4 *)(x + m + 3);
        pX4 = (const xtfloatx4 *)(x + m + 7);
        aX3 = AE_LA128_PP(pX3);
        aX4 = AE_LA128_PP(pX4);

        CONST_SX2X2(r00, r01, 0); CONST_SX2X2(r10, r11, 0);
        CONST_SX2X2(r20, r21, 0); CONST_SX2X2(r30, r31, 0);
        CONST_SX2X2(r40, r41, 0); CONST_SX2X2(r50, r51, 0);
        CONST_SX2X2(r60, r61, 0); CONST_SX2X2(r70, r71, 0);

        /* preload data from x */
        AE_LSX2X2_IP(x0, x2, pX0, 4 * sizeof(float32_t));
        x1 = XT_SEL32_LH_SX2(x0, x2);

        __Pragma("loop_count min=1");
        for (n = 0; n < (N >> 2); n++)
        {
            /* load data from x */
            AE_LSX2X2_IP(x4, x6, pX1, 4 * sizeof(float32_t));
            AE_LSX2IP(x8, castxcc(xtfloatx2, pX2), 4 * sizeof(float32_t));
            AE_LASX2X2_IP(x3, x5, aX3, pX3);
            AE_LASX2X2_IP(x7, x9, aX4, pX4);
            /* load data from y */
            AE_LSX2X2_IP(e0, e1, pE, 4 * sizeof(float32_t));
            /* compute correlation of 8 values */
            MADD_SX2X2(r00, r10, x0, x1, e0, e0);
            MADD_SX2X2(r20, r30, x2, x3, e0, e0);
            MADD_SX2X2(r40, r50, x4, x5, e0, e0);
            MADD_SX2X2(r60, r70, x6, x7, e0, e0);
            MADD_SX2X2(r01, r11, x2, x3, e1, e1);
            MADD_SX2X2(r21, r31, x4, x5, e1, e1);
            MADD_SX2X2(r41, r51, x6, x7, e1, e1);
            MADD_SX2X2(r61, r71, x8, x9, e1, e1);
            /* reload data and shift input line for the next iteration */
            AE_LSX2X2_IP(x0, x2, pX0, 4 * sizeof(float32_t));
            x1 = x5;
        }
        /* reduction add */
        ADD_SX2X2(r00, r10, r00, r10, r01, r11);
        ADD_SX2X2(r20, r30, r20, r30, r21, r31);
        AE_DSELSX2(r0, r1, r10, r00, sel);
        AE_DSELSX2(r2, r3, r30, r20, sel);
        ADD_SX2X2(r0, r1, r0, r2, r1, r3);
        /* update and save IR */
        AE_LSX2X2_IP(h0, h1, pH, -4 * (int)sizeof(float32_t));
        MADD_SX2X2(h0, h1, b, b, r1, r0);
        AE_SSX2X2_IP(h0, h1, pHw, -4 * (int)sizeof(float32_t));

        /* reduction add */
        ADD_SX2X2(r40, r50, r40, r50, r41, r51);
        ADD_SX2X2(r60, r70, r60, r70, r61, r71);
        AE_DSELSX2(r0, r1, r50, r40, sel);
        AE_DSELSX2(r2, r3, r70, r60, sel);
        ADD_SX2X2(r0, r1, r0, r2, r1, r3);
        /* update and save IR */
        AE_LSX2X2_IP(h0, h1, pH, -4 * (int)sizeof(float32_t));
        MADD_SX2X2(h0, h1, b, b, r1, r0);
        AE_SSX2X2_IP(h0, h1, pHw, -4 * (int)sizeof(float32_t));
    }
} /* fir_blmsf() */
#else
// for scalar FPU

void fir_blmsf( float32_t * e, float32_t * h, const float32_t * r,
                const float32_t * x, float32_t norm, float32_t mu, int N, int M )
{
  float32_t b;
  int m, n;

  const xtfloat*  restrict pX = (const xtfloat*)x;
  const xtfloat*  restrict pR = (const xtfloat*)r;
  xtfloat*  restrict pE = (      xtfloat*)e;
  const xtfloat*  restrict pe = (      xtfloat*)e;
  const xtfloat*    restrict ph = (const xtfloat*)h;
  xtfloat*    restrict pH = (xtfloat*)h;
  xtfloat*    restrict pH_wr = (xtfloat*)h;
  NASSERT(e);
  NASSERT(h);
  NASSERT(r);
  NASSERT(x);
  NASSERT_ALIGN(e, 8);
  NASSERT_ALIGN(h, 8);
  NASSERT_ALIGN(r, 8);
  NASSERT_ALIGN(x, 8);
  NASSERT(N>0 && M>0);
  NASSERT(M % 8 == 0 && N % 8 == 0);

  if (N <= 0 || M <= 0) return;

  /* estimate error */
  for (n = 0; n<(N); n+=4)
  {
    xtfloat r0, r1, r2, r3;
    xtfloat s0, s1, s2, s3;
    xtfloat x0, x1, x2, x3,x4;
    s0 = XT_CONST_S(0);
    s1 = XT_CONST_S(0);
    s2 = XT_CONST_S(0);
    s3 = XT_CONST_S(0);
    pX = (xtfloat*)((uintptr_t)(x+n));
    XT_LSIP(x0, pX, 4);
    XT_LSIP(x1, pX, 4);
    XT_LSIP(x2, pX, 4);
    XT_LSIP(x3, pX, 4);
    for (m = 0; m<M; m+=2)
    {
      xtfloat hm;
      hm = ph[M - 1 - m];
      XT_LSIP(x4, pX, 4);
      XT_MADD_S(s0, hm, x0);
      XT_MADD_S(s1, hm, x1);
      XT_MADD_S(s2, hm, x2);
      XT_MADD_S(s3, hm, x3);
      hm = ph[M - 1 - m - 1];
      XT_MADD_S(s0, hm, x1);
      XT_MADD_S(s1, hm, x2);
      XT_MADD_S(s2, hm, x3);
      XT_MADD_S(s3, hm, x4);
      x0 = x2; x1 = x3; x2 = x4;
      XT_LSIP(x3, pX, 4);
    }
    XT_LSIP(r0, pR, 4);
    XT_LSIP(r1, pR, 4);
    XT_LSIP(r2, pR, 4);
    XT_LSIP(r3, pR, 4);
    s0=XT_SUB_S(r0, s0);
    s1=XT_SUB_S(r1, s1);
    s2=XT_SUB_S(r2, s2);
    s3=XT_SUB_S(r3, s3);
    XT_SSIP(s0, pE, 4);
    XT_SSIP(s1, pE, 4);
    XT_SSIP(s2, pE, 4);
    XT_SSIP(s3, pE, 4);
  }
  /* update impluse response */
  b = mu / norm;
  for (m = 0; m<M; m+=4)
  {
    xtfloat h0, h1, h2, h3;
    xtfloat s0, s1, s2, s3;
    xtfloat x0, x1, x2, x3, x4;
    s0 = XT_CONST_S(0);
    s1 = XT_CONST_S(0);
    s2 = XT_CONST_S(0);
    s3 = XT_CONST_S(0);
    pX = (xtfloat*)((uintptr_t)(x + m));
    XT_LSIP(x0, pX, 4);
    XT_LSIP(x1, pX, 4);
    XT_LSIP(x2, pX, 4);
    XT_LSIP(x3, pX, 4);

    for (n = 0; n<N; n+=2)
    {
      xtfloat en;
      en = pe[n];
      XT_LSIP(x4, pX, 4);
      XT_MADD_S(s0, en, x0);
      XT_MADD_S(s1, en, x1);
      XT_MADD_S(s2, en, x2);
      XT_MADD_S(s3, en, x3);
      en = pe[n+1];
      XT_MADD_S(s0, en, x1);
      XT_MADD_S(s1, en, x2);
      XT_MADD_S(s2, en, x3);
      XT_MADD_S(s3, en, x4);
      x0 = x2; x1 = x3; x2 = x4;
      XT_LSIP(x3, pX, 4);
    }
    pH = (xtfloat*)((uintptr_t)(h + M - 1 - m));
    XT_LSXP(h0, pH, -4);
    XT_LSXP(h1, pH, -4);
    XT_LSXP(h2, pH, -4);
    XT_LSXP(h3, pH, -4);
    XT_MADD_S(h0, b, s0);
    XT_MADD_S(h1, b, s1);
    XT_MADD_S(h2, b, s2);
    XT_MADD_S(h3, b, s3);
    pH_wr = (xtfloat*)((uintptr_t)(h + M - 1 - m));
    XT_SSXP(h0, pH_wr, -4);
    XT_SSXP(h1, pH_wr, -4);
    XT_SSXP(h2, pH_wr, -4);
    XT_SSXP(h3, pH_wr, -4);
  }
} /* fir_blmsf() */
#endif
