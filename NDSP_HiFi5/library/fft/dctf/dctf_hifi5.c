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
   NatureDSP Signal Processing Library. FFT part
    Discrete Cosine Transform, Type II 
    C code optimized for HiFi4
   Integrit, 2006-2019
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fft.h"
#include "common.h"
#include "dct2_twd.h"

#include "common_fpu.h"

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(int,dctf,( float32_t *  y,float32_t * x,dct_handle_t h))
#elif HAVE_VFPU

#ifndef AE_DSEL32_HL_LH_SX2
/*
   Equal to:
   a = AE_SEL32_HL_SX2(c, d)
   b = AE_SEL32_LH_SX2(c, d)
*/
#define AE_DSEL32_HL_LH_SX2(a,b,c,d) \
{\
    ae_int16x4 aa,bb,cc,dd,sel; \
    sel = AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32X2((7L<<24)|(5L<<16)|(6<<8)|(4), (1L<<24)|(3L<<16)|(0<<8)|(2) )); \
    cc = AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(c)); \
    dd = AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(d)); \
    AE_DSEL16X4(aa,bb,cc,dd,sel); \
    a = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(aa)); \
    b = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(bb)); \
}
#endif

#ifndef AE_DSEL32_HH_LLSWP_SX2
/*
   Equal to:
   a = AE_SEL32_HH_SX2(c, d)
   b = AE_SEL32_LL_SX2(d, c)
*/
#define AE_DSEL32_HH_LLSWP_SX2(a,b,c,d) \
{\
    ae_int16x4 aa,bb,cc,dd,sel; \
    sel = AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32X2((7L<<24)|(1L<<16)|(6<<8)|(0), (3L<<24)|(5L<<16)|(2<<8)|(4) )); \
    cc = AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(c)); \
    dd = AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(d)); \
    AE_DSEL16X4(aa,bb,cc,dd,sel); \
    a = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(aa)); \
    b = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(bb)); \
}
#endif

/*
Complex FFT for length N==32

  Input:
    x[N] - complex signal
    ptwd[N*3/4] - twiddle factors
  Output:
    y[N] - complex spectrum
*/
static void fft_cplxf_N32(complex_float *y, complex_float *x, complex_float *twd);

/*-------------------------------------------------------------------------
  Discrete Cosine Transform.
  These functions apply DCT (Type II, Type IV) to input.
  Scaling:
      +-----------------------+--------------------------------------+
      |      Function         |           Scaling options            |
      +-----------------------+--------------------------------------+
      |       dct_16x16       |  3 - fixed scaling before each stage |
      |       dct_32x16       |  3 - fixed scaling before each stage |
      |       dct_32x32       |  3 - fixed scaling before each stage |
      |       dct4_32x16      |  3 - fixed scaling before each stage |
      |       dct4_32x32      |  3 - fixed scaling before each stage |
      +-----------------------+--------------------------------------+
  NOTES:
     1. DCT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after 
     the call.
     2. N - DCT size (depends on selected DCT handle)

  Precision: 
  16x16  16-bit input/outputs, 16-bit twiddles
  32x16  32-bit input/outputs, 16-bit twiddles
  32x32  32-bit input/outputs, 32-bit twiddles
  f      floating point

  Input:
  x[N]        input signal
  h           DCT handle
  scalingOpt  scaling option (see table above) 
              not applicable to the floating point function
  Output:
  y[N]        transform output
  
  Returned value:
              total number of right shifts occurred during scaling 
              procedure 
  Restriction:
  x,y         should not overlap
  x,y         aligned on 16-bytes boundary
-------------------------------------------------------------------------*/
#define SZ_CF32 (sizeof(complex_float))
#define SZ_F32  (sizeof(float32_t))
/* 1/sqrt(2.0) */
static const union ufloat32uint32 _invsqrt2f_ = { 0x3f3504f3 };
int dctf     ( float32_t * restrict y,float32_t * restrict x,dct_handle_t h)
/*
    Reference Matlab code:
    function y=dctf(x)
    N=numel(x);
    y(1:N/2)     =x(1:2:N);
    y(N:-1:N/2+1)=x(2:2:N);
    % take fft of N/2
    y=fft(y(1:2:N)+j*y(2:2:N));
    w=exp(i*pi/2*(0:N-1)/N);
    % DCT split algorithm
    Y0=y(1);
    T0=real(Y0)+imag(Y0);
    T1=real(Y0)-imag(Y0);
    z(1      )= real(T0);%*sqrt(2)/2;
    z(N/2+1  )= real(T1)*sqrt(2)/2;
    for k=2:N/4
        Y0=y(k);
        Y1=y(N/2+2-k);
        COSI=(w(4*(k-1)+1));
        W1=w(k);
        W2=w(N/2+2-k);
        S=Y0+Y1;
        D=Y0-Y1;
        T0=i*real(D)+imag(S);
        T1=i*imag(D)+real(S);
        Y0=  ( imag(T0)*imag(COSI)-real(T0)*real(COSI)) + ...
           i*( real(T0)*imag(COSI)+imag(T0)*real(COSI));
        T0=0.5*(T1-Y0);
        T1=0.5*(T1+Y0);
        z(k      )= real(T0)*real(W1)+imag(T0)*imag(W1);
        z(N+2-k  )= real(T0)*imag(W1)-imag(T0)*real(W1);
        z(N/2+2-k)= real(T1)*real(W2)-imag(T1)*imag(W2);
        z(N/2+k  )= real(T1)*imag(W2)+imag(T1)*real(W2);
    end
    W1=w(N/4+1);
    T0=y(N/4+1);
    z(N/4+1  )= real(T0)*real(W1)-imag(T0)*imag(W1);
    z(N+1-N/4)= real(T0)*imag(W1)+imag(T0)*real(W1);
    y=z;
*/
{
    const union ufloat32uint32 *rfft_split_twd; 
    const tdct2_twd *descr=(const tdct2_twd *)h;
    int N;
    const xtfloatx2 *restrict p0_twd;
    const xtfloatx2 *restrict p1_twd;
    const xtfloatx2 *restrict p2_twd;
    const xtfloatx4 *restrict p0_ld;
    const xtfloatx2 *restrict p1_ld;
          xtfloatx4 *restrict p0_stx2;
          xtfloatx4 *restrict p1_stx2;
          xtfloat   *restrict p0_st;
          xtfloat   *restrict p1_st;
          xtfloat   *restrict p2_st;
          xtfloat   *restrict p3_st;
    int k, n;
    int N2, N4;

    NASSERT_ALIGN(x,16);
    NASSERT_ALIGN(y,16);
    NASSERT(x!=y);
    NASSERT(descr->magic==MAGIC_DCT2_F);
    NASSERT(descr->N==32 || descr->N==64);
    N=descr->N;
    rfft_split_twd=(const union ufloat32uint32 *)descr->rfft_split_twd;
    N2 = N>>1;
    N4 = N2>>1;

    /* permute inputs */
    p0_ld  = (const xtfloatx4 *)x;
    p0_stx2 = (xtfloatx4 *)(y);
    p1_stx2 = (xtfloatx4 *)(y+N-4);
    __Pragma("loop_count min=2");
    for (n=0; n<(N4>>1); n++)
    {
      xtfloatx2 t0, t1, t2, t3, y0, y1, y2, y3;
      /* y[n]    =x[2*n+0] */
      /* y[N-1-n]=x[2*n+1] */
      AE_LSX2X2_IP(t0, t1, p0_ld, 2*SZ_CF32);
      AE_LSX2X2_IP(t2, t3, p0_ld, 2*SZ_CF32);
      AE_DSEL32_HH_LLSWP_SX2(y0, y1, t0, t1);
      AE_DSEL32_HH_LLSWP_SX2(y2, y3, t2, t3);
      AE_SSX2X2_IP(y0, y2, p0_stx2, 2*SZ_CF32);
      AE_SSX2X2_IP(y3, y1, p1_stx2, -2*(int)SZ_CF32);
    }

    /* compute fft(N/2) */
    if(N == 32)
    {
      xtfloatx2 a0, a1, a2, a3, a4, a5, a6, a7;
      xtfloatx2 a8, a9, aA, aB, aC, aD, aE, aF;
      xtfloatx2 b0, b1, b2, b3, b4, b5, b6, b7;
      xtfloatx2 b8, b9, bA, bB, bC, bD, bE, bF;
      xtfloatx2 tw1, tw2, tw3, tw5, tw6, tw7;
      xtfloatx2 tw9, twA, twB, twD, twE, twF;
      p0_twd = (const xtfloatx2 *)descr->fft_twd;
      p0_ld = (const xtfloatx4 *)(y);
      p0_stx2 = (xtfloatx4 *)(x);

      /* FFT, radix-4 butterflies, fully unrolled */
      {
        /* load input samples */
        AE_LSX2X2_IP(a0, a4, p0_ld, 2*SZ_CF32);
        AE_LSX2X2_IP(a8, aC, p0_ld, 2*SZ_CF32);
        AE_LSX2X2_IP(a1, a5, p0_ld, 2*SZ_CF32);
        AE_LSX2X2_IP(a9, aD, p0_ld, 2*SZ_CF32);
        AE_LSX2X2_IP(a2, a6, p0_ld, 2*SZ_CF32);
        AE_LSX2X2_IP(aA, aE, p0_ld, 2*SZ_CF32);
        AE_LSX2X2_IP(a3, a7, p0_ld, 2*SZ_CF32);
        AE_LSX2X2_IP(aB, aF, p0_ld, 2*SZ_CF32);
        /* compute butterfly */
        /* 1-st substage */
        ADD_SX2X2(b0, b1, a0, a1, a2, a3);
        SUB_SX2X2(b2, b3, a0, a1, a2, a3);
        ADD_SX2X2(b4, b5, a4, a5, a6, a7);
        SUB_SX2X2(b6, b7, a4, a5, a6, a7);
        ADD_SX2X2(b8, b9, a8, a9, aA, aB);
        SUB_SX2X2(bA, bB, a8, a9, aA, aB);
        ADD_SX2X2(bC, bD, aC, aD, aE, aF);
        SUB_SX2X2(bE, bF, aC, aD, aE, aF);
        /* 2-nd substage */
        ADD_SX2X2(a0, a4, b0, b4, b1, b5);
        ADDANDSUBJC_SX2(a3, a1, b2, b3);
        SUB_SX2X2(a2, a6, b0, b4, b1, b5);
        ADDANDSUBJC_SX2(a7, a5, b6, b7);
        ADD_SX2X2(a8, aC, b8, bC, b9, bD);
        ADDANDSUBJC_SX2(aB, a9, bA, bB);
        SUB_SX2X2(aA, aE, b8, bC, b9, bD);
        ADDANDSUBJC_SX2(aF, aD, bE, bF);
        /* multiply by twiddle factors */
        AE_LSX2X2_IP(tw1, tw2, castxcc(xtfloatx4,p0_twd), 2*SZ_CF32);
        AE_LSX2X2_IP(tw3, tw5, castxcc(xtfloatx4,p0_twd), 2*SZ_CF32);
        AE_LSX2X2_IP(tw6, tw7, castxcc(xtfloatx4,p0_twd), 2*SZ_CF32);
        AE_LSX2X2_IP(tw9, twA, castxcc(xtfloatx4,p0_twd), 2*SZ_CF32);
        AE_LSX2X2_IP(twB, twD, castxcc(xtfloatx4,p0_twd), 2*SZ_CF32);
        AE_LSX2X2_IP(twE, twF, castxcc(xtfloatx4,p0_twd), 2*SZ_CF32);

        b0 = a0;  b1 = a1;
        b2 = a2;  b3 = a3;
        b4 = a4;
        MULC_SX2(b5, b7, a5, a7, tw5, tw7);
        b6 = MULC_S(a6, tw6);
        b8 = a8;  bC = aC;
        MULC_SX2(b9, bD, a9, aD, tw9, twD);
        MULC_SX2(bA, bE, aA, aE, twA, twE);
        MULC_SX2(bB, bF, aB, aF, twB, twF);
        /* Two middle quartiles are swapped on all but the last stage to use the bit reversal
         * permutation instead of the digit reverse. */
        a0 = b0; a1 = b4; 
        a2 = b8; a3 = bC; 
        a4 = b2; a5 = b6; 
        a6 = bA; a7 = bE; 
        a8 = b1; a9 = b5; 
        aA = b9; aB = bD; 
        aC = b3; aD = b7; 
        aE = bB; aF = bF; 
        /* Last stage (radix-4) with bit reversal
         * permutation. */
        /* Last stage is fully unrolled */
        /* 1-st substage */
        ADD_SX2X2(b0, b1, a0, a1, a2, a3);
        SUB_SX2X2(b2, b3, a0, a1, a2, a3);
        ADD_SX2X2(b4, b5, a4, a5, a6, a7);
        SUB_SX2X2(b6, b7, a4, a5, a6, a7);
        ADD_SX2X2(b8, b9, a8, a9, aA, aB);
        SUB_SX2X2(bA, bB, a8, a9, aA, aB);
        ADD_SX2X2(bC, bD, aC, aD, aE, aF);
        SUB_SX2X2(bE, bF, aC, aD, aE, aF);
        /* 2-nd substage */
        ADDANDSUBC_SX2 (a0, a2, b0, b1);
        ADDANDSUBJC_SX2(a3, a1, b2, b3);
        ADDANDSUBC_SX2 (a4, a6, b4, b5);
        ADDANDSUBJC_SX2(a7, a5, b6, b7);
        ADDANDSUBC_SX2 (a8, aA, b8, b9);
        ADDANDSUBJC_SX2(aB, a9, bA, bB);
        ADDANDSUBC_SX2 (aC, aE, bC, bD);
        ADDANDSUBJC_SX2(aF, aD, bE, bF);
        /* Save results with bit-reversal permutation */
        AE_SSX2X2_IP(a0, a8, p0_stx2, 2*SZ_CF32);
        AE_SSX2X2_IP(a4, aC, p0_stx2, 2*SZ_CF32);
        AE_SSX2X2_IP(a1, a9, p0_stx2, 2*SZ_CF32);
        AE_SSX2X2_IP(a5, aD, p0_stx2, 2*SZ_CF32);
        AE_SSX2X2_IP(a2, aA, p0_stx2, 2*SZ_CF32);
        AE_SSX2X2_IP(a6, aE, p0_stx2, 2*SZ_CF32);
        AE_SSX2X2_IP(a3, aB, p0_stx2, 2*SZ_CF32);
        AE_SSX2X2_IP(a7, aF, p0_stx2, 2*SZ_CF32);
      }
    }
    else
    {
        fft_cplxf_N32((complex_float*)x,(complex_float*)y,(complex_float*)descr->fft_twd);
    }

    /* make final DCT transformation of FFT outputs */
    {
      xtfloatx2 t00, t01, t10, t11, y00, y01, y10, y11;
      xtfloatx2 w0, w1, w2, w3, s0, s1, d0, d1, cos1, cos2, c05;
      xtfloat b0, b1, re, im, invsqrt2f;
      ae_int32x2 t32x2;
      ae_valign al_p2, al_p3;
      p0_ld  = (const xtfloatx4 *)x;
      p1_ld  = (const xtfloatx2 *)x+N2-1;
      p0_twd = (const xtfloatx2 *)rfft_split_twd+4;
      p1_twd = (const xtfloatx2 *)rfft_split_twd+1;
      p2_twd = (const xtfloatx2 *)rfft_split_twd+(N2-1);
      p0_st = (xtfloat *)y;
      p2_st = p0_st+N2-1;
      p1_st = p0_st+N2;
      p3_st = p0_st+N-1;

      /* Load constants */
      c05 = (xtfloatx2)0.5f;/* 0.5 */
      invsqrt2f = AE_LSI((xtfloat *)&_invsqrt2f_, 0);/* 1/sqrt(2) */

      AE_LSX2IP(y00, castxcc(xtfloatx2,p0_ld), SZ_CF32);
      /* b0 = y0.re + y0.im */
      /* b1 = y0.re - y0.im */
      re = HIGH_S(y00);
      im = LOW_S (y00);
      ADDANDSUB_S(b0, b1, re, im);
      b1 = MUL_S(b1, invsqrt2f);
      AE_SSIP(b0, p0_st, SZ_F32);
      AE_SSIP(b1, p1_st, SZ_F32);

      {
        AE_LSX2IP(y00, castxcc(xtfloatx2,p0_ld), SZ_CF32);
        AE_LSX2XP(y10, p1_ld, -(int)SZ_CF32);
        AE_LSX2IP(cos1, p0_twd, 4*SZ_CF32);
        AE_LSX2IP(w0  , p1_twd,   SZ_CF32);
        AE_LSX2XP(w1  , p2_twd, -(int)SZ_CF32);

        ADDANDSUB_SX2(s0, d0, y00, y10);
        /* t0.re = s.im; t0.im = d.re */
        /* t1.re = s.re; t1.im = d.im */
        AE_DSEL32_HL_LH_SX2(t01, t00, s0, d0);

        /* y0 = conj(t0*cosi)      */
        /* t0 = 0.5*w0*conj(t1+y0) */
        /* t1 = 0.5*w1*    (t1-y0) */
        y00 = MULMUX_S(t00, cos1, 1);
        MADDMUX_S(y00, t00, cos1, 7);
        ADDANDSUB_SX2(t00, t01, t01, y00);
        MUL_SX2X2(w0, w1, w0, w1, c05, c05);
        t00 = CONJC_SX2(t00);
        MULC_SX2(t00, t01, w0, w1, t00, t01);

        /* y[k    ] = t0.re */
        /* y[N-k  ] = t0.im */
        /* y[N/2-k] = t1.re */
        /* y[N/2+k] = t1.im */
        re = HIGH_S(t00);
        t32x2 = AE_MOVINT32X2_FROMXTFLOATX2(t00);
        AE_SSIP(re, p0_st, SZ_F32);/* save real part */
        AE_S32_L_IP(t32x2, castxcc(ae_int32,p3_st), -(int)SZ_F32);/* save imag part */
        re = HIGH_S(t01);
        t32x2 = AE_MOVINT32X2_FROMXTFLOATX2(t01);
        AE_SSIP(re, p2_st, -(int)SZ_F32);/* save real part */
        AE_S32_L_IP(t32x2, castxcc(ae_int32,p1_st), SZ_F32);/* save imag part*/
      }
      al_p3 = AE_ZALIGN64();
      al_p2 = AE_ZALIGN64();
      __Pragma("loop_count min=2")
      for (k=1; k<(N4>>1); k++)
      {
        AE_LSX2X2_IP(y00, y01, p0_ld, 2*SZ_CF32);
        AE_LSX2XP(y10, p1_ld, -(int)SZ_CF32);
        AE_LSX2XP(y11, p1_ld, -(int)SZ_CF32);
        AE_LSX2IP(cos1, p0_twd, 4*SZ_CF32);
        AE_LSX2IP(cos2, p0_twd, 4*SZ_CF32);
        AE_LSX2X2_IP(w0, w2, castxcc(xtfloatx4,p1_twd), 2*SZ_CF32);
        AE_LSX2XP(w1, p2_twd, -(int)SZ_CF32);
        AE_LSX2XP(w3, p2_twd, -(int)SZ_CF32);
      
        ADDANDSUB_SX2(s0, d0, y00, y10);
        ADDANDSUB_SX2(s1, d1, y01, y11);
        /* t0.re = s.im; t0.im = d.re */
        /* t1.re = s.re; t1.im = d.im */
        AE_DSEL32_HL_LH_SX2(t01, t00, s0, d0);
        AE_DSEL32_HL_LH_SX2(t11, t10, s1, d1);

        /* y0 = conj(t0*cosi)      */
        /* t0 = 0.5*w0*conj(t1+y0) */
        /* t1 = 0.5*w1*    (t1-y0) */
        MULMUX_SX2X2(y00, y01, t00, t10, cos1, cos2, 1);
        MADDMUX_SX2X2(y00, y01, t00, t10, cos1, cos2, 7);
        ADDANDSUB_SX2(t00, t01, t01, y00);
        ADDANDSUB_SX2(t10, t11, t11, y01);
        MUL_SX2X2(w0, w2, w0, w2, c05, c05);
        MUL_SX2X2(w1, w3, w1, w3, c05, c05);
        MULCCONJ_SX2(t00, t10, w0, w2, t00, t10);
        MULC_SX2(t01, t11, w1, w3, t01, t11);

        /* y[k    ]= t0.re */
        /* y[N-k  ]= t0.im */
        /* y[N/2-k]= t1.re */
        /* y[N/2+k]= t1.im */
        y00 = AE_SEL32_HH_SX2(t00, t10);
        y11 = AE_SEL32_LL_SX2(t00, t10);
        y10 = AE_SEL32_HH_SX2(t01, t11);
        y01 = AE_SEL32_LL_SX2(t01, t11);
        AE_SSX2IP(y00, castxcc(xtfloatx2,p0_st), 2*SZ_F32);/* save real part */
        AE_SSX2IP(y01, castxcc(xtfloatx2,p1_st), 2*SZ_F32);/* save imag part*/
        AE_SASX2RIP(y10, al_p2, castxcc(xtfloatx2,p2_st));/* save real part */
        AE_SASX2RIP(y11, al_p3, castxcc(xtfloatx2,p3_st));/* save imag part */
      }
      AE_SA64NEG_FP(al_p3, p3_st);
      AE_SA64NEG_FP(al_p2, p2_st);

      t00 = AE_LSX2I((const xtfloatx2 *)p0_ld, 0);
      w0 = AE_LSX2I(p1_twd, 0);
      t00 = MULC_S(t00, w0);
      re = HIGH_S(t00);
      im = LOW_S (t00);
      AE_SSI(re, p0_st, 0);
      AE_SSI(im, p3_st, 0);
    }

    return 0;
} /* dctf() */

/*
Complex FFT for length N==32

  Input:
    x[N] - complex signal
    ptwd[N*3/4] - twiddle factors
  Output:
    y[N] - complex spectrum
*/
void fft_cplxf_N32(complex_float *y, complex_float *x, complex_float *twd)
{
  const int N = 32;
  const int logN = 5;
  const xtfloatx4 *restrict p_twd;
  const xtfloatx4 *restrict p0_ld;
  const xtfloatx4 *restrict p1_ld;
  const xtfloatx4 *restrict p2_ld;
  const xtfloatx4 *restrict p3_ld;
        xtfloatx4 *restrict p0_st;
        xtfloatx4 *restrict p1_st;
        xtfloatx4 *restrict p2_st;
        xtfloatx4 *restrict p3_st;
        xtfloatx4 *restrict p4_st;
        xtfloatx4 *restrict p5_st;
        xtfloatx4 *restrict p6_st;
        xtfloatx4 *restrict p7_st;
  xtfloatx2 tw01, tw02, tw03, tw11, tw12, tw13;
  int N4;
  int n;
  unsigned int idx, bitrevstride;

  NASSERT( x );
  NASSERT( y );
  NASSERT( twd );
  NASSERT( x != y );
  NASSERT_ALIGN( x, 16 );
  NASSERT_ALIGN( y, 16 );
  NASSERT_ALIGN( twd, 16 );

  N4 = N>>2;
  /* Set the pointer to the twiddle table */
  p_twd = (const xtfloatx4 *)(twd);

  /*
   * Perform the first stage. We use DIF, all permutations are deferred
   * until the last stage.
   */
  {
    p0_st = (xtfloatx4 *)(x);
    p1_st = (xtfloatx4 *)((complex_float *)p0_st + N4);
    p2_st = (xtfloatx4 *)((complex_float *)p1_st + N4);
    p3_st = (xtfloatx4 *)((complex_float *)p2_st + N4);
    p0_ld = p0_st;
    p1_ld = p1_st;
    p2_ld = p2_st;
    p3_ld = p3_st;
    /* Radix-4 butterfly */
    __Pragma("ymemory (p_twd)")
    for ( n=0; n<(N4>>1); n++ )
    {
        xtfloatx2 a01, a11, a21, a31;
        xtfloatx2 b01, b11, b21, b31;
        xtfloatx2 a00, a10, a20, a30;
        xtfloatx2 b00, b10, b20, b30;

        /* load input samples */
        AE_LSX2X2_IP(a00, a01, p0_ld, 2*SZ_CF32);
        AE_LSX2X2_IP(a10, a11, p1_ld, 2*SZ_CF32);
        AE_LSX2X2_IP(a20, a21, p2_ld, 2*SZ_CF32);
        AE_LSX2X2_IP(a30, a31, p3_ld, 2*SZ_CF32);

        /* compute butterfly */
        /* 1-st substage */
        ADDANDSUBC_SX2(b00, b20, a00, a20);
        ADDANDSUBC_SX2(b10, b30, a10, a30);
        ADDANDSUBC_SX2(b01, b21, a01, a21);
        ADDANDSUBC_SX2(b11, b31, a11, a31);

        /* 2-nd substage */
        ADDANDSUBC_SX2 (a00, a20, b00, b10);
        ADDANDSUBJC_SX2(a30, a10, b20, b30);
        ADDANDSUBC_SX2 (a01, a21, b01, b11);
        ADDANDSUBJC_SX2(a31, a11, b21, b31);

        /* multiply by twiddle factors */
        AE_LSX2X2_IP(tw01, tw02, p_twd, 2*SZ_CF32);
        AE_LSX2X2_IP(tw03, tw11, p_twd, 2*SZ_CF32);
        AE_LSX2X2_IP(tw12, tw13, p_twd, 2*SZ_CF32);
        b00 = a00;    b01 = a01;
        MULC_SX2(b10, b11, a10, a11, tw01, tw11);
        MULC_SX2(b20, b21, a20, a21, tw02, tw12);
        MULC_SX2(b30, b31, a30, a31, tw03, tw13);

        /* Two middle quartiles are swapped on all but the last stage to use the bit reversal
        * permutation instead of the digit reverse. */
        AE_SSX2IP(b00, castxcc(xtfloatx2, p0_st), SZ_CF32);
        AE_SSX2IP(b01, castxcc(xtfloatx2, p0_st), SZ_CF32);
        AE_SSX2X2_IP(b20, b21, p1_st, 2 * SZ_CF32);
        AE_SSX2X2_IP(b10, b11, p2_st, 2 * SZ_CF32);
        AE_SSX2X2_IP(b30, b31, p3_st, 2 * SZ_CF32);
    }
  }
  __Pragma("no_reorder");
  /*
    Last stage (radix-8) with bit reversal permutation.
  */
  {
    xtfloatx2 a0, a1, a2, a3, a4, a5, a6, a7;
    xtfloatx2 b0, b1, b2, b3, b4, b5, b6, b7;
    xtfloatx2 a3_, a7_, invsqrt2_;
    int N8 = N4 >> 1;
    idx = 0;
    bitrevstride = 0x80000000U >> (logN-1);
    
    invsqrt2_ = XT_LSI((xtfloat *)&_invsqrt2f_, 0);/* 1.0/sqrt(2.0) */
    p0_ld = (const xtfloatx4 *)(x);
    p1_ld = p0_ld+2;

    p0_st = (xtfloatx4 *)(y);
    p1_st = (xtfloatx4 *)((complex_float *)p0_st+N8);
    p2_st = (xtfloatx4 *)((complex_float *)p1_st+N8);
    p3_st = (xtfloatx4 *)((complex_float *)p2_st+N8);
    p4_st = (xtfloatx4 *)((complex_float *)p3_st+N8);
    p5_st = (xtfloatx4 *)((complex_float *)p4_st+N8);
    p6_st = (xtfloatx4 *)((complex_float *)p5_st+N8);
    p7_st = (xtfloatx4 *)((complex_float *)p6_st+N8);
    
    for ( n=0; n<N8; n++ )
    {
      /* Load input samples */
      AE_LSX2X2_IP(a0, a1, castxcc(xtfloatx4, p0_ld), 2*SZ_CF32);
      AE_LSX2X2_IP(a2, a3, castxcc(xtfloatx4, p0_ld), 6*SZ_CF32);
      AE_LSX2X2_IP(a4, a5, castxcc(xtfloatx4, p1_ld), 2*SZ_CF32);
      AE_LSX2X2_IP(a6, a7, castxcc(xtfloatx4, p1_ld), 6*SZ_CF32);

      /* Compute butterfly */
      /* 1-st substage */
      ADDANDSUBC_SX2(b0, b4, a0, a4);
      ADDANDSUBC_SX2(b1, b5, a1, a5);
      ADDANDSUBC_SX2(b2, b6, a2, a6);
      ADDANDSUBC_SX2(b3, b7, a3, a7);
      /* 2-nd substage */
      ADDANDSUBC_SX2 (a0, a4, b0, b2);
      ADDANDSUBC_SX2 (a1, a5, b1, b3);
      ADDANDSUBJC_SX2(a6, a2, b4, b6);
      ADDANDSUBJC_SX2(a7, a3, b5, b7);
      /* 3-rd substage */
      MULJC_SX2X2(a3_, a7_, a3, a7);
      a3 = a3 - a3_;
      a7 = a7 + a7_;
      MUL_SX2X2(a3, a7, a3, a7, invsqrt2_, invsqrt2_);
      ADDANDSUBC_SX2 (b0, b4, a0, a1);
      ADDANDSUBC_SX2 (b1, b5, a2, a3);
      ADDANDSUBJC_SX2(b6, b2, a4, a5);
      ADDANDSUBC_SX2 (b7, b3, a6, a7);

      /* Store samples */
      XT_SSX2X(b0, (xtfloatx2 *)p0_st, idx);
      XT_SSX2X(b1, (xtfloatx2 *)p1_st, idx);
      XT_SSX2X(b2, (xtfloatx2 *)p2_st, idx);
      XT_SSX2X(b3, (xtfloatx2 *)p3_st, idx);
      XT_SSX2X(b4, (xtfloatx2 *)p4_st, idx);
      XT_SSX2X(b5, (xtfloatx2 *)p5_st, idx);
      XT_SSX2X(b6, (xtfloatx2 *)p6_st, idx);
      XT_SSX2X(b7, (xtfloatx2 *)p7_st, idx);

      idx = AE_ADDBRBA32(idx, bitrevstride);
    }
  }
} /* fft_cplxf_N32() */


#else

int dctf     ( float32_t * restrict y,float32_t * restrict x,dct_handle_t h)
{
    const tdct2_twd *descr=(const tdct2_twd *)h;
    const xtfloat * restrict pX0;
    const xtfloat * restrict pX1;
    const xtfloat * restrict pCosi;
    const xtfloat * restrict pTwd1;
    const xtfloat * restrict pTwd2;
    xtfloat y0_re,y1_re,t0_re,t1_re,cosi_re,w1_re,w2_re,s_re,d_re;
    xtfloat y0_im,y1_im,t0_im,t1_im,cosi_im,w1_im,w2_im,s_im,d_im;
    const ae_int32x2 * restrict p0_ld;
          ae_int32x2 * restrict p0_stx2;
          ae_int32x2 * restrict p1_stx2;

    int N, k, n, twd_stride;
    const complex_float* dct_twd;
    const complex_float *fft_twd;

    NASSERT_ALIGN(x,8);
    NASSERT_ALIGN(y,8);
    NASSERT(x!=y);
    NASSERT(descr->magic==MAGIC_DCT2_F);
    N = descr->N;
    ASSERT(N==32 || N==64);

    /* permute inputs */
    p0_ld  = (const ae_int32x2 *)x;
    p0_stx2 = (ae_int32x2 *)y;
    p1_stx2 = (ae_int32x2 *)(y+N-2);
    __Pragma("loop_count min=1")
    for (n=0; n<(N>>2); n++)
    {
        ae_int32x2 t0,t1,y0,y1;
      /* y[n]    =x[2*n+0] */
      /* y[N-1-n]=x[2*n+1] */
      AE_L32X2_IP(t0, p0_ld, sizeof(ae_int32x2));
      AE_L32X2_IP(t1, p0_ld, sizeof(ae_int32x2));
      y0 = AE_SEL32_HH(t0, t1);
      y1 = AE_SEL32_LL(t1, t0);
      AE_S32X2_IP(y0, p0_stx2,       sizeof(ae_int32x2));
      AE_S32X2_XP(y1, p1_stx2, -(int)sizeof(ae_int32x2));
    }
    /* compute fft(N/2) */
    dct_twd=(const complex_float*)descr->rfft_split_twd;
    fft_twd=(const complex_float *)descr->fft_twd;
    twd_stride=1;
    fft_cplxf_ie((complex_float*)x,(complex_float*)y,fft_twd,twd_stride,N/2);
    /* make final DCT transformation of FFT outputs */
    pX0=(const xtfloat *)(x);
    pX1=(const xtfloat *)(x+N-2);
    pCosi=(const xtfloat *)(dct_twd+4*twd_stride);
    pTwd1=(const xtfloat *)(dct_twd+  twd_stride);
    pTwd2=(const xtfloat *)(dct_twd+(N/2-1)*twd_stride);
    y0_im=XT_LSI(pX0,sizeof(xtfloat)); XT_LSIP(y0_re,pX0, 2*sizeof(xtfloat));
    s_re=XT_ADD_S(y0_re,y0_im);
    s_im=XT_SUB_S(y0_re,y0_im);
    y[0  ]= s_re;
    y[N/2]= XT_MUL_S(s_im,0.707106781f);
    for (k=1; k<(N>>2); k++)
    {
        y0_im=XT_LSI(pX0,sizeof(xtfloat)); XT_LSIP(y0_re,pX0, 2*sizeof(xtfloat));
        y1_im=XT_LSI(pX1,sizeof(xtfloat)); XT_LSXP(y1_re,pX1,-2*(int)sizeof(xtfloat));

        w1_im=XT_LSI(pTwd1,sizeof(xtfloat)); XT_LSXP(w1_re,pTwd1, twd_stride*2*sizeof(xtfloat));
        w2_im=XT_LSI(pTwd2,sizeof(xtfloat)); XT_LSXP(w2_re,pTwd2,-twd_stride*2*sizeof(xtfloat));
        cosi_im=XT_LSI(pCosi,sizeof(xtfloat)); XT_LSXP(cosi_re,pCosi,4*twd_stride*2*sizeof(xtfloat));

        s_re=XT_ADD_S(y0_re,y1_re);
        s_im=XT_ADD_S(y0_im,y1_im);
        d_re=XT_SUB_S(y0_re,y1_re);
        d_im=XT_SUB_S(y0_im,y1_im);

        t0_re=s_im; t0_im=d_re;
        t1_re=s_re; t1_im=d_im;
        y0_re=XT_MUL_S(t0_im,cosi_im);XT_MSUB_S(y0_re,t0_re,cosi_re);
        y0_im=XT_MUL_S(t0_re,cosi_im);XT_MADD_S(y0_im,t0_im,cosi_re);
        t0_re=XT_MUL_S(0.5f,XT_SUB_S(t1_re,y0_re));
        t0_im=XT_MUL_S(0.5f,XT_SUB_S(t1_im,y0_im));
        t1_re=XT_MUL_S(0.5f,XT_ADD_S(t1_re,y0_re));
        t1_im=XT_MUL_S(0.5f,XT_ADD_S(t1_im,y0_im));
        s_re=XT_MUL_S(t0_re,w1_re); XT_MADD_S(s_re,t0_im,w1_im);
        s_im=XT_MUL_S(t0_re,w1_im); XT_MSUB_S(s_im,t0_im,w1_re);
        d_re=XT_MUL_S(t1_re,w2_re); XT_MSUB_S(d_re,t1_im,w2_im);
        d_im=XT_MUL_S(t1_re,w2_im); XT_MADD_S(d_im,t1_im,w2_re);
        y[k         ]= s_re;
        y[(N>>1)+k  ]= d_im;
        y[N-k       ]= s_im;
        y[N-k-(N>>1)]= d_re;
    }
    t0_im=XT_LSI(pX0,sizeof(xtfloat)); XT_LSIP(t0_re,pX0, 2*sizeof(xtfloat));
    w1_im=XT_LSI(pTwd1,sizeof(xtfloat)); XT_LSXP(w1_re,pTwd1, twd_stride*2*sizeof(xtfloat));
    d_re=XT_MUL_S(t0_re,w1_re); XT_MSUB_S(d_re,t0_im,w1_im);
    d_im=XT_MUL_S(t0_re,w1_im); XT_MADD_S(d_im,t0_im,w1_re);
    y[k]  = d_re;
    y[N-k]= d_im;

    return 0;
}
#endif
