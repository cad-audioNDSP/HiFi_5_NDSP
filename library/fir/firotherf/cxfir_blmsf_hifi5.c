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
    Blockwise Adaptive LMS Algorithm for Complex Data, floating point
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
DISCARD_FUN(void,cxfir_blmsf,( complex_float * e, complex_float * h,
                const complex_float * r,
                const complex_float * x,
                float32_t norm, float32_t mu,
                int          N, int       M ))

#elif (HAVE_VFPU)

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

void cxfir_blmsf ( complex_float * e,
                   complex_float * h,
             const complex_float * r,
             const complex_float * x,
             float32_t norm, float32_t mu,
             int N, int M )
{
    const xtfloatx4 * restrict pX;
    const xtfloatx4 * restrict S0;
    const xtfloatx4 * restrict S1;
    const xtfloatx4 * restrict pH;
    const xtfloatx4 * restrict pR;
          xtfloatx4 * restrict pE;
          xtfloatx4 * restrict pHw;
    xtfloatx2 r00, r10, r20, r30;
    xtfloatx2 r01, r11, r21, r31;
    xtfloatx2 r02, r12, r22, r32;
    xtfloatx2 r03, r13, r23, r33;
    xtfloatx2 x0, x1, x2, x3, x4;
    xtfloatx2 h0, h1, e0, e1;
    xtfloatx2 b;
    int m, n;

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

    WUR_AE_CBEGIN0((uintptr_t)(x));
    WUR_AE_CEND0((uintptr_t)(x + N + M));

    /* estimate error */
    pX = (const xtfloatx4 *)x;
    pR = (const xtfloatx4 *)r;
    pE = (      xtfloatx4 *)e;
    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 2); n++)
    {
        pH = (const xtfloatx4 *)(h + M - 2);
        /* preload data from x */
        AE_LSX2X2_IP(x0, x1, pX, 2 * sizeof(complex_float));
        S0 = pX; S1 = pX;
        pX += 1;

        AE_LSX2X2_IP(r00, r10, pR, 2 * sizeof(complex_float));
        AE_LSX2X2_IP(r20, r30, pR, 2 * sizeof(complex_float));
        CONST_SX2X2(r01, r11, 0); CONST_SX2X2(r21, r31, 0);
        CONST_SX2X2(r02, r12, 0); CONST_SX2X2(r22, r32, 0);
        CONST_SX2X2(r03, r13, 0); CONST_SX2X2(r23, r33, 0);

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 1); m++)
        {
            /* load data from x */
            AE_LSX2X2_IP(x2, x3, S1, 2 * sizeof(complex_float));
            x4 = AE_LSX2I((xtfloatx2*)S1, 0);
            /* load data from y */
            AE_LSX2X2_IP(h1, h0, pH, -2 * (int)sizeof(complex_float));
            /* compute correlation of 4 values */
            MADDMUX_SX2X2(r00, r10, x0, x1, h0, h0, 6);
            MADDMUX_SX2X2(r01, r11, x0, x1, h0, h0, 7);
            MADDMUX_SX2X2(r20, r30, x2, x3, h0, h0, 6);
            MADDMUX_SX2X2(r21, r31, x2, x3, h0, h0, 7);
            MADDMUX_SX2X2(r02, r12, x1, x2, h1, h1, 6);
            MADDMUX_SX2X2(r03, r13, x1, x2, h1, h1, 7);
            MADDMUX_SX2X2(r22, r32, x3, x4, h1, h1, 6);
            MADDMUX_SX2X2(r23, r33, x3, x4, h1, h1, 7);
            /* reload data for the next iteration */
            AE_LSX2X2_XC(x0, x1, S0, 2 * sizeof(complex_float));
        }
        /* reduction add */
        ADD_SX2X2(r00, r10, r00, r10, r01, r11);
        ADD_SX2X2(r20, r30, r20, r30, r21, r31);
        ADD_SX2X2(r02, r12, r02, r12, r03, r13);
        ADD_SX2X2(r22, r32, r22, r32, r23, r33);
        ADD_SX2X2(r00, r10, r00, r10, r02, r12);
        ADD_SX2X2(r20, r30, r20, r30, r22, r32);
        /* save computed samples */
        AE_SSX2X2_IP(r00, r10, pE, 2 * sizeof(complex_float));
        AE_SSX2X2_IP(r20, r30, pE, 2 * sizeof(complex_float));
    }

    /* update impluse response */
    b = XT_MUL_SX2(mu, XT_RECIP_SX2(norm));
    pX = (const xtfloatx4 *)x;
    pHw= (      xtfloatx4 *)(h + M - 2);
    pH = pHw;
    __Pragma("loop_count min=1");
    for (m = 0; m < (M >> 2); m++)
    {
        pE = (xtfloatx4 *)e;
        /* preload data from x */
        AE_LSX2X2_IP(x0, x1, pX, 2 * sizeof(complex_float));
        S0 = pX; S1 = pX;
        pX += 1;

        CONST_SX2X2(r00, r10, 0); CONST_SX2X2(r20, r30, 0);
        CONST_SX2X2(r01, r11, 0); CONST_SX2X2(r21, r31, 0);
        CONST_SX2X2(r02, r12, 0); CONST_SX2X2(r22, r32, 0);
        CONST_SX2X2(r03, r13, 0); CONST_SX2X2(r23, r33, 0);

        __Pragma("loop_count min=1");
        for (n = 0; n < (N >> 1); n++)
        {
            /* load data from x */
            AE_LSX2X2_IP(x2, x3, S1, 2 * sizeof(complex_float));
            x4 = AE_LSX2I((xtfloatx2*)S1, 0);
            /* load data from y */
            AE_LSX2X2_IP(e0, e1, pE, 2 * sizeof(complex_float));
            /* compute correlation of 4 values */
            MADDMUX_SX2X2(r00, r10, x0, x1, e0, e0, 4);
            MADDMUX_SX2X2(r01, r11, x0, x1, e0, e0, 5);
            MADDMUX_SX2X2(r20, r30, x2, x3, e0, e0, 4);
            MADDMUX_SX2X2(r21, r31, x2, x3, e0, e0, 5);
            MADDMUX_SX2X2(r02, r12, x1, x2, e1, e1, 4);
            MADDMUX_SX2X2(r03, r13, x1, x2, e1, e1, 5);
            MADDMUX_SX2X2(r22, r32, x3, x4, e1, e1, 4);
            MADDMUX_SX2X2(r23, r33, x3, x4, e1, e1, 5);
            /* reload data for the next iteration */
            AE_LSX2X2_XC(x0, x1, S0, 2 * sizeof(complex_float));
        }
        /* reduction add */
        ADD_SX2X2(r00, r10, r00, r10, r01, r11);
        ADD_SX2X2(r20, r30, r20, r30, r21, r31);
        ADD_SX2X2(r02, r12, r02, r12, r03, r13);
        ADD_SX2X2(r22, r32, r22, r32, r23, r33);
        ADD_SX2X2(r00, r10, r00, r10, r02, r12);
        ADD_SX2X2(r20, r30, r20, r30, r22, r32);
        /* update and save IR */
        AE_LSX2X2_IP(h0, h1, pH, -2 * (int)sizeof(complex_float));
        MADD_SX2X2(h0, h1, b, b, r10, r00);
        AE_SSX2X2_IP(h0, h1, pHw, -2 * (int)sizeof(complex_float));
        AE_LSX2X2_IP(h0, h1, pH, -2 * (int)sizeof(complex_float));
        MADD_SX2X2(h0, h1, b, b, r30, r20);
        AE_SSX2X2_IP(h0, h1, pHw, -2 * (int)sizeof(complex_float));
    }
} /* cxfir_blmsf() */
#else
/* for scalar FPU */
void cxfir_blmsf    ( complex_float * e, complex_float * h,
                const complex_float * r,
                const complex_float * x,
                float32_t norm, float32_t mu,
                int          N, int       M )
{
          xtfloat* restrict pH;
    const xtfloat* restrict pX;
    const xtfloat* restrict pR;
          xtfloat* restrict pE;
    xtfloat b;
    int m,n;
    NASSERT(e);
    NASSERT(h);
    NASSERT(r);
    NASSERT(x);
    NASSERT_ALIGN(e,8);
    NASSERT_ALIGN(h,8);
    NASSERT_ALIGN(r,8);
    NASSERT_ALIGN(x,8);
    NASSERT(N>0 && M>0);
    NASSERT(M%8==0 && N%8==0);

    /* estimate error */
    pR=(const xtfloat*)r;
    pE=(      xtfloat*)e;
    for (n=0; n<N; n+=2)
    {
        xtfloat s0_re,s0_im,h0_r,h0_i,x0_r,x0_i;
        xtfloat s1_re,s1_im;
        XT_LSIP(s0_re,pR,sizeof(float32_t));
        XT_LSIP(s0_im,pR,sizeof(float32_t));
        XT_LSIP(s1_re,pR,sizeof(float32_t));
        XT_LSIP(s1_im,pR,sizeof(float32_t));
        pX=(const xtfloat*)(x+n);
        pH=(      xtfloat*)(h+M-1);
        XT_LSIP(x0_r,pX,sizeof(float32_t));
        XT_LSIP(x0_i,pX,sizeof(float32_t));
        __Pragma("loop_count min=8")
        for (m=0; m<M; m++)
        {
            h0_i=XT_LSI (pH,sizeof(float32_t));
            XT_LSXP(h0_r,pH,-2*(int)sizeof(float32_t));
            XT_MSUB_S(s0_re,h0_r,x0_r); XT_MSUB_S(s0_re,h0_i,x0_i);
            XT_MADD_S(s0_im,h0_i,x0_r); XT_MSUB_S(s0_im,h0_r,x0_i);
            XT_LSIP(x0_r,pX,sizeof(float32_t));
            XT_LSIP(x0_i,pX,sizeof(float32_t));
            XT_MSUB_S(s1_re,h0_r,x0_r); XT_MSUB_S(s1_re,h0_i,x0_i);
            XT_MADD_S(s1_im,h0_i,x0_r); XT_MSUB_S(s1_im,h0_r,x0_i);
        }
        XT_SSIP(s0_re,pE,sizeof(float32_t));
        XT_SSIP(s0_im,pE,sizeof(float32_t));
        XT_SSIP(s1_re,pE,sizeof(float32_t));
        XT_SSIP(s1_im,pE,sizeof(float32_t));
    }
    /* update impluse response */
    b=XT_MUL_S(mu,XT_RECIP_S(norm));
    pH=(xtfloat*)(h+(M-2));
    for (m=0; m<M; m+=2)
    {
        xtfloat s0_re,s0_im,x0_r,x0_i,e0_r,e0_i;
        xtfloat s1_re,s1_im;
        xtfloat h0_re,h0_im;
        xtfloat h1_re,h1_im;
        s0_re=0.f;        s0_im=0.f;
        s1_re=0.f;        s1_im=0.f;
        pX=(const xtfloat*)(x+m);
        pE=(      xtfloat*)(e);
        XT_LSIP(x0_r,pX,sizeof(xtfloat));
        XT_LSIP(x0_i,pX,sizeof(xtfloat));
        __Pragma("loop_count min=8")
        for (n=0; n<N; n++)
        {
            XT_LSIP(e0_r,pE,sizeof(xtfloat));
            XT_LSIP(e0_i,pE,sizeof(xtfloat));
            XT_MADD_S(s0_re,x0_r,e0_r); XT_MADD_S(s0_re,x0_i,e0_i);
            XT_MSUB_S(s0_im,x0_r,e0_i); XT_MADD_S(s0_im,x0_i,e0_r);
            XT_LSIP(x0_r,pX,sizeof(xtfloat));
            XT_LSIP(x0_i,pX,sizeof(xtfloat));
            XT_MADD_S(s1_re,x0_r,e0_r); XT_MADD_S(s1_re,x0_i,e0_i);
            XT_MSUB_S(s1_im,x0_r,e0_i); XT_MADD_S(s1_im,x0_i,e0_r);
        }
        h1_re=XT_LSI (pH,0*sizeof(float32_t));
        h1_im=XT_LSI (pH,1*sizeof(float32_t));
        h0_re=XT_LSI (pH,2*sizeof(float32_t));
        h0_im=XT_LSI (pH,3*sizeof(float32_t));
        XT_MADD_S(h1_re,s1_re,b);
        XT_MADD_S(h1_im,s1_im,b);
        XT_MADD_S(h0_re,s0_re,b);
        XT_MADD_S(h0_im,s0_im,b);
        XT_SSI (h1_im,pH,1*sizeof(float32_t));
        XT_SSI (h0_re,pH,2*sizeof(float32_t));
        XT_SSI (h0_im,pH,3*sizeof(float32_t));
        XT_SSXP(h1_re,pH,-4*(int)sizeof(float32_t));
    }
} /* cxfir_blmsf() */
#endif
