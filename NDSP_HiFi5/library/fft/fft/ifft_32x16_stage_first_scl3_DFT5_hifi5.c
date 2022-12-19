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
    Complex-valued FFT stages with butterflies radix-2, radix-3, radix-5
    with static data scaling: 32-bit data, 16-bit twiddle factors
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/
#if 0 /* not used in library for now */

#include "NatureDSP_types.h"
#include "common.h"
#include "fft_32x16_stages.h"

/*
DFT5 algorithm:
x - input complex vector
y - output complex vector
y = fft(x)
w1 =  exp(-1j*2*pi/5);
w2 =  exp(-1j*2*pi*2/5);

y = zeros(5,1);
s1 = (x1+x4);
s2 = (x2+x3);
d1 = (x1-x4);
d2 = (x2-x3);

y(1) = x0 + s1 + s2;
y(2) = x0 + (s1*real(w1) + s2*real(w2)) + 1j*(d1*imag(w1) + d2*imag(w2));
y(5) = x0 + (s1*real(w1) + s2*real(w2)) - 1j*(d1*imag(w1) + d2*imag(w2));
y(3) = x0 + (s1*real(w2) + s2*real(w1)) + 1j*(d1*imag(w2) - d2*imag(w1));
y(4) = x0 + (s1*real(w2) + s2*real(w1)) - 1j*(d1*imag(w2) - d2*imag(w1));

*/
#define DFT5X1(x0, x1, x2, x3, x4)\
{ \
  ae_int32x2 s1, s2, d1, d2;       \
  ae_int32x2 y0, y1;               \
  ae_int32x2 t0, t1, t2, t3, wtmp; \
  ae_int16x4 w1, w2;               \
  wtmp = AE_MOVDA32X2(0x278E278E,0x00008644); \
  w1 = AE_MOVINT16X4_FROMINT32X2(wtmp);       \
  wtmp = AE_MOVDA32X2(0x98729872,0x0000B4C3); \
  w2 = AE_MOVINT16X4_FROMINT32X2(wtmp);       \
  AE_ADDANDSUB32S(s1, d1, x1, x4);            \
  AE_ADDANDSUB32S(s2, d2, x2, x3);            \
  t0 = AE_MULFP32X16X2RAS_H(s1, w1);          \
  t1 = AE_MULFP32X16X2RAS_H(s2, w2);          \
  t2 = AE_MULFP32X16X2RAS_H(s1, w2);          \
  t3 = AE_MULFP32X16X2RAS_H(s2, w1);          \
  y0 = AE_ADD32S(x0, AE_ADD32S(t0, t1));      \
  y1 = AE_ADD32S(x0, AE_ADD32S(t2, t3));      \
  t0 = AE_MULFC32X16RAS_L(d1, w1);            \
  t1 = AE_MULFC32X16RAS_L(d2, w2);            \
  t2 = AE_MULFC32X16RAS_L(d2, w1);            \
  t3 = AE_MULFC32X16RAS_L(d1, w2);            \
  t0 = AE_ADD32S(t0, t1);                     \
  t1 = AE_SUB32S(t3, t2);                     \
  x0 = AE_ADD32S(x0, AE_ADD32S(s1, s2));      \
  AE_ADDANDSUB32S(x1, x4, y0, t0);            \
  AE_ADDANDSUB32S(x2, x3, y1, t1);            \
}

/*
 *  32x16 IFFT first stage Radix 5, scalingOption=3
 */
int ifft_32x16_stage_first_scl3_DFT5(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  int shift;
  int i;
  const int R = 5; // stage radix
  const int stride = N / R;
  ae_int32x2 * restrict px0;
  ae_int32x2 * restrict py0;
  const ae_int16x4 * restrict ptwd;
  ae_int16x4 tw12, tw34;
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  shift = 3;
  ptwd = (const ae_int16x4 *)tw;
  py0 = (ae_int32x2 *)y;
  px0 = (ae_int32x2 *)x + N - 4*stride;
  /*
  ifft(x) = fft( [ x(1), x(end:-1:2)] )
  */
  {
    /* First butterfly radix 4 */
    ae_int32x2 x0, x1, x2, x3, x4;

    AE_L16X4_IP(tw12, ptwd, sizeof(ae_int16x4));
    AE_L16X4_IP(tw34, ptwd, sizeof(ae_int16x4));

    AE_L32X2_XP(x4, px0, stride*sizeof(ae_int32x2));
    AE_L32X2_XP(x3, px0, stride*sizeof(ae_int32x2));
    AE_L32X2_XP(x2, px0, stride*sizeof(ae_int32x2));
    AE_L32X2_XP(x1, px0, stride*sizeof(ae_int32x2));
    x0 = AE_L32X2_X((ae_int32x2*)x, 0);
    x0 = AE_SRAI32R(x0, 3);
    x1 = AE_SRAI32R(x1, 3);
    x2 = AE_SRAI32R(x2, 3);
    x3 = AE_SRAI32R(x3, 3);
    x4 = AE_SRAI32R(x4, 3);

    DFT5X1(x0, x1, x2, x3, x4);
    x1 = AE_MULFC32X16RAS_H(x1, tw12);
    x2 = AE_MULFC32X16RAS_L(x2, tw12);
    x3 = AE_MULFC32X16RAS_H(x3, tw34);
    x4 = AE_MULFC32X16RAS_L(x4, tw34);

    AE_S32X2_IP(x0, py0, sizeof(ae_int32x2));
    AE_S32X2_IP(x1, py0, sizeof(ae_int32x2));
    AE_S32X2_IP(x2, py0, sizeof(ae_int32x2));
    AE_S32X2_IP(x3, py0, sizeof(ae_int32x2));
    AE_S32X2_IP(x4, py0, sizeof(ae_int32x2));
  }
  px0 = (ae_int32x2 *)x + N - 4*stride - 1;
  __Pragma("loop_count min=2");
  for (i = 1; i <stride; i++)
  {
    ae_int32x2 x0, x1, x2, x3, x4;

    AE_L16X4_IP(tw12, ptwd, sizeof(ae_int16x4));
    AE_L16X4_IP(tw34, ptwd, sizeof(ae_int16x4));

    AE_L32X2_XP(x4, px0, stride*sizeof(ae_int32x2));
    AE_L32X2_XP(x3, px0, stride*sizeof(ae_int32x2));
    AE_L32X2_XP(x2, px0, stride*sizeof(ae_int32x2));
    AE_L32X2_XP(x1, px0, stride*sizeof(ae_int32x2));
    AE_L32X2_XP(x0, px0, (-4 * stride - 1)*sizeof(ae_int32x2));
    x0 = AE_SRAI32R(x0, 3);
    x1 = AE_SRAI32R(x1, 3);
    x2 = AE_SRAI32R(x2, 3);
    x3 = AE_SRAI32R(x3, 3);
    x4 = AE_SRAI32R(x4, 3);

    DFT5X1(x0, x1, x2, x3, x4);
    x1 = AE_MULFC32X16RAS_H(x1, tw12);
    x2 = AE_MULFC32X16RAS_L(x2, tw12);
    x3 = AE_MULFC32X16RAS_H(x3, tw34);
    x4 = AE_MULFC32X16RAS_L(x4, tw34);

    AE_S32X2_IP(x0, py0, sizeof(ae_int32x2));
    AE_S32X2_IP(x1, py0, sizeof(ae_int32x2));
    AE_S32X2_IP(x2, py0, sizeof(ae_int32x2));
    AE_S32X2_IP(x3, py0, sizeof(ae_int32x2));
    AE_S32X2_IP(x4, py0, sizeof(ae_int32x2));
  }
  *v *= R;
  return shift; 
} /* ifft_32x16_stage_first_scl3_DFT5() */
#endif
