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
    DCT Type II: 32-bit data, 16-bit twiddle factors
    C code optimized for HiFi4
    IntegrIT, 2006-2019
*/

/* Signal Processing Library API. */
#include "NatureDSP_Signal_fft.h"
/* Common utility macros. */
#include "common.h"
/* Twiddle factor tables for DCTs. */
#include "dct2_twd.h"

#ifndef AE_DSEL32X2_HH_LLSWP
/*
   Equal to:
   a = AE_SEL32_HH(c, d)
   b = AE_SEL32_LL(d, c)
*/
#define AE_DSEL32X2_HH_LLSWP(a,b,c,d) \
{\
    ae_int16x4 aa,bb,cc,dd,sel; \
    sel = AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32X2((7L<<24)|(1L<<16)|(6<<8)|(0), (3L<<24)|(5L<<16)|(2<<8)|(4) )); \
    cc = AE_MOVINT16X4_FROMINT32X2(c); \
    dd = AE_MOVINT16X4_FROMINT32X2(d); \
    AE_DSEL16X4(aa,bb,cc,dd,sel); \
    a = AE_MOVINT32X2_FROMINT16X4(aa); \
    b = AE_MOVINT32X2_FROMINT16X4(bb); \
}
#endif

/*
   scaled fft with reordering
   NOTE: y is input and output, x - temporary
*/
void fft16_32x16(int32_t *y, int32_t *x, const int16_t *ptwd);/* N=16 */
void fft32_32x16(int32_t *y, int32_t *x, const int16_t *ptwd);/* N=32 */
/* pointer to complex fft function with reordering */
typedef void(*cfftProc_func)(int32_t *y, int32_t *x, const int16_t *ptwd);

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
int dct_32x16( int32_t * restrict y, int32_t * restrict x, dct_handle_t h, int scalingOpt )
{
  //
  // MATLAB reference code:
  //
  //  function y = dct_ref(x)
  //  % Compute DCT Type II by conjugate-pair split-radix algorithm.
  //  % Relates to MATLAB's implementation as dct(x)*(numel(x)/2)^0.5 == dct_ref(x).
  //  N = numel(x);
  //  % Reorder input data as required by the split-radix DCT-II algorithm
  //  x = [x(1:2:N-1)',wrev(x(2:2:N)')];
  //  % Real-valued FFT of size N through the complex-valued FFT of size N/2.
  //  s = fft(x(1:2:N-1)+1j*x(2:2:N));
  //  % Real-to-complex spectrum conversion twiddle factors
  //  spc = -1j*exp(-2*pi*1j*(1:N/4)/N);
  //  % Conjugate-pair split-radix algorithm twiddle factors
  //  spr = reshape([2^-0.5,exp(-2*pi*1j*(1:N/4)/(4*N)); ...
  //                      0,exp(-2*pi*1j*(N/2-1:-1:N/4)/(4*N))],1,N/2+2);
  //  y = [spr(1)*(real(s(1))+imag(s(1)));zeros(N/2-1,1); ...
  //       spr(1)*(real(s(1))-imag(s(1)));zeros(N/2-1,1)];
  //  for n=2:N/4
  //    a0 = s(n);
  //    a1 = s(N/2+2-n);
  //    b0 = 1/2*(a0+conj(a1));
  //    b1 = 1/2*(a0-conj(a1))*spc(n-1);
  //    a0 = (b0+b1)*spr((n-1)*2+1);
  //    a1 = conj(b0-b1)*spr((n-1)*2+2);
  //    y(n) = real(a0);
  //    y(N/2+2-n) = real(a1);
  //    y(N/2+n) = -imag(a1);
  //    y(N+2-n) = -imag(a0);
  //  end
  //  y(N/4+1) = real(conj(s(N/4+1))*spr(N/2+1));
  //  y(3*N/4+1) = -imag(conj(s(N/4+1))*spr(N/2+2));
  //

    const tdct2_twd * descr=(const tdct2_twd *)h;
    const ae_int32x4 * restrict p_y0;
    const ae_int32x4 * restrict p_y1;
    ae_int32x2       * restrict p_z0;
    ae_int32x2       * restrict p_z1;
    ae_int32x2       * restrict p_z2;
    ae_int32x2       * restrict p_z3;
    const ae_int16x4 * restrict p_cos;
    const ae_int16x8 * restrict p_w;
    ae_valign al_z0, al_z2;
    ae_int16x4  vC0, vW0, vW1;
    ae_int32x2  vA0, vA1, vA2, vA3, vB0, vB1, vB2, vB3, vT0;
    int N, i;

    NASSERT_ALIGN16( y );
    NASSERT_ALIGN16( x );
    NASSERT(scalingOpt==3);
    NASSERT(descr->magic==MAGIC_DCT2_16);
    N = descr->N;
    NASSERT( N==32 || N==64 );

    /*
     * FFT of half-size with reordering
     */
    switch (N)
    {
    case 32:
        fft16_32x16(x, y, (const int16_t *)descr->fft_twd);
        break;
    case 64:
        fft32_32x16(x, y, (const int16_t *)descr->fft_twd);
        break;
    default:
        NASSERT(0);
    }

    /*
     * DCT split algorithm
     */
    p_y0  = (const ae_int32x4 *)(x);
    p_y1  = (const ae_int32x4 *)(x+N-4);
    p_cos = (const ae_int16x4 *)(descr->rfft_split_twd);
    p_w   = (const ae_int16x8 *)(descr->dct_twd);
    p_z0  = (      ae_int32x2 *)(y);
    p_z1  = (      ae_int32x2 *)(y+N-2);
    p_z2  = (      ae_int32x2 *)(y+N/2);
    p_z3  = (      ae_int32x2 *)(y+N/2-2);

    WUR_AE_SAR(1);
    /*** process 1st sample ***/
    AE_L32X2X2_IP(vA0, vT0, p_y0, 2*sizeof(ae_int32x2));
    vW0 = AE_MOVINT16X4_FROMINT32X2(AE_L32X2_X((ae_int32x2 *)p_w, N*sizeof(int16_t)-2*sizeof(ae_int16x4)));
    vB0 = AE_SUBADD32S_HL_LH(vA0, vA0);
    vB0 = AE_MULFP32X16X2RAS_L(vB0, vW0);

    AE_S32_H_IP(vB0, castxcc(ae_int32,p_z2), sizeof(ae_int32));
    AE_S32_L_IP(vB0, castxcc(ae_int32,p_z0), sizeof(ae_int32));

    al_z0 = al_z2 = AE_ZALIGN64();
    __Pragma("loop_count min=2");
    __Pragma("no_unroll");
    for (i = 1; i < (N>>3); i++)
    {
      /*
        Y0=Y(k);
        Y1=Y(N/2+2-k);
        COSI=cosi(k);
        W1=w(k);
        W2=w(N/2+2-k);
        S=Y0+Y1;
        D=Y0-Y1;
        T0=i*real(D)+imag(S);
        T1=i*imag(D)+real(S);
        T1=T1/2;
        Y0=  (imag(T0)*imag(COSI)-real(T0)*real(COSI)) + ...
              i*(real(T0)*imag(COSI)+imag(T0)*real(COSI));
        Y0=Y0/2;
        T0=T1-Y0;
        T1=T1+Y0;
      */
      vA0 = vT0;
      AE_L32X2X2_IP(vA2, vT0, p_y0, sizeof(ae_int32x4));
      AE_L32X2X2_IP(vA3, vA1, p_y1, -(int)sizeof(ae_int32x4));

      vB0 = AE_SUBADD32S(vA0, vA1);
      vB1 = AE_ADDSUB32S(vA0, vA1);
      vB2 = AE_SUBADD32S(vA2, vA3);
      vB3 = AE_ADDSUB32S(vA2, vA3);
      AE_L16X4_IP(vC0, p_cos, sizeof(ae_int16x4));
      AE_MULFPC32X16X2RAS(vB0, vB2, vB0, vB2, vC0);

      AE_ADDANDSUBRNG32(vA1, vA0, vB1, vB0);
      AE_ADDANDSUBRNG32(vA3, vA2, vB3, vB2);

      /*
        z(k      )= real(T0)*real(W1)-imag(T0)*imag(W1);
        z(N+2-k  )=-real(T0)*imag(W1)-imag(T0)*real(W1);
        z(N/2+2-k)= real(T1)*real(W2)+imag(T1)*imag(W2);
        z(N/2+k  )=-real(T1)*imag(W2)+imag(T1)*real(W2);
      */
      AE_L16X4X2_IP(vW0, vW1, p_w, 2*sizeof(ae_int16x4));
      AE_MULFPC32X16X2RAS(vB0, vB2, vA0, vA2, vW0);
      AE_MULFPC32X16X2RAS(vB1, vB3, vA1, vA3, vW1);
      AE_DSEL32X2_HH_LLSWP(vA0, vA1, vB0, vB2);
      AE_DSEL32X2_HH_LLSWP(vA3, vA2, vB3, vB1);
      vA1 = AE_NEG32S(vA1);

      AE_SA32X2_IP(vA0, al_z0, p_z0);
      AE_S32X2_XP (vA1, p_z1, -(int)sizeof(ae_int32x2));
      AE_SA32X2_IP(vA2, al_z2, p_z2);
      AE_S32X2_XP (vA3, p_z3, -(int)sizeof(ae_int32x2));
    }

    /*** last iteration ***/
    vA0 = vT0;
    AE_L32X2X2_I(vA3, vA1, p_y1, 0);

    vB0 = AE_SUBADD32S(vA0, vA1);
    vB1 = AE_ADDSUB32S(vA0, vA1);
    AE_L16X4_IP(vC0, p_cos, sizeof(ae_int16x4));
    vB0 = AE_MULFC32X16RAS_H(vB0, vC0);

    AE_ADDANDSUBRNG32(vA1, vA0, vB1, vB0);

    AE_L16X4X2_IP(vW0, vW1, p_w, 2*sizeof(ae_int16x4));
    vB0 = AE_MULFC32X16RAS_H(vA0, vW0);
    vB1 = AE_MULFC32X16RAS_H(vA1, vW1);
    vT0 = AE_NEG32S(vB0);

    /*** middle sample ***/
    /*
      W1=w(N/4+1);
      T0=Y(N/4+1);
      z(N/4+1  )= real(T0)*real(W1)+imag(T0)*imag(W1);
      z(N+1-N/4)=-real(T0)*imag(W1)+imag(T0)*real(W1);
    */
    vB3 = AE_MULFC32X16RAS_L(vA3, vW1);

    vA0 = AE_SEL32_HH(vB0, vB3);
    vA1 = AE_SEL32_LL(vB3, vT0);
    AE_SA32X2_IP(vA0, al_z0, p_z0);
    AE_S32X2_I  (vA1, p_z1, 0);
    AE_SA64POS_FP(al_z2, p_z2);
    AE_S32_L_I(vB1, (ae_int32 *)p_z2, 0);
    AE_S32_H_I(vB1, (ae_int32 *)p_z3, sizeof(ae_int32));
    AE_SA64POS_FP(al_z0, p_z0);

    return 30-NSA(N);
} // dct_32x16()
