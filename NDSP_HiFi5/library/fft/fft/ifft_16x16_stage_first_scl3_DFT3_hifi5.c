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
    with static data scaling: 16-bit data, 16-bit twiddle factors
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

#if 0 /* not used in library for now */

#include "NatureDSP_types.h"
#include "common.h"
#include "fft_16x16_stages.h"

/*
DFT3 algorithm:
x - input complex vector
y - output complex vector
y = fft(x)
y = [ x(1) + x(2)  + x(3);
x(1) + (x(2) + x(3))*cos(2*pi/3) - 1j*(x(2) - x(3))*sin(2*pi/3);
x(1) + (x(2) + x(3))*cos(2*pi/3) + 1j*(x(2) - x(3))*sin(2*pi/3) ]

*/
#define DFT3X2(x0, x1, x2)\
{\
    ae_int16x4 s0, s1, d0, c;\
    ae_int32x2 c32;          \
    c32 = AE_MOVDA32(0x6EDA9126);      \
    c = AE_MOVINT16X4_FROMINT32X2(c32);\
    s0 = AE_ADD16S(x1, x2);            \
    d0 = AE_SUB16S(x1, x2);            \
    s1 = AE_ADD16S(x0, s0);            \
    s0 = AE_SRAI16(s0, 1);             \
    s0 = AE_SUB16S(x0, s0);            \
    d0 = AE_MULFP16X4RAS(d0, c);       \
    d0 = AE_SEL16_2301(d0, d0);        \
    x0 = s1;                           \
    x1 = AE_SUB16S(s0, d0);            \
    x2 = AE_ADD16S(s0, d0);            \
}

int ifft_16x16_stage_first_scl3_DFT3(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
  #define shift 3
  int i;
  const ae_int16x4 * restrict px0;
  const ae_int32   * restrict px0s;
        ae_int16x4 * restrict py0;
  const ae_int16x4 * restrict ptwd;
    const int N3=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),715827883));//N/3
  const int stride = N3;
  ae_int16x4 x0, x1, x2;
  ae_int32x2 x00, x01, x02;
  ae_int32x2 x10, x11, x12;
  ae_int32x2 x20, x21, x22;
  ae_int16x4 tw012, tw112, c1;
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  ptwd = (const ae_int16x4 *)tw;
  py0 = (ae_int16x4 *)y;
  c1 = 1;
  /*
  ifft(x) = fft( [ x(1), x(end:-1:2)] )
  */
  px0s = (ae_int32 *)x + N - 2*stride;
  {
    /* First butterfly radix 3 */
    AE_L32_XP(x02, px0s, stride*sizeof(ae_int32));
    x01 = AE_L32_I(px0s, 0);
    x00 = AE_L32_I((ae_int32*)x, 0);
    x0 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x00));
    x1 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x01));
    x2 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x02));

    x0 = AE_SRAI16R(x0, shift);
    x1 = AE_SRAI16R(x1, shift);
    x2 = AE_SRAI16R(x2, shift);
    DFT3X2(x0, x1, x2);
    x00 = AE_SEXT32X2D16_32(x0);
    x01 = AE_SEXT32X2D16_32(x1);
    x02 = AE_SEXT32X2D16_32(x2);

    AE_L16X4_IP(tw012, ptwd, sizeof(ae_int16x4));
    x01 = AE_MULFC32X16RAS_H(x01, tw012);
    x02 = AE_MULFC32X16RAS_L(x02, tw012);
  }
  px0 = (ae_int16x4 *)((complex_fract16*)x + N-2 - 2*stride);
  __Pragma("loop_count min=1");
  for (i = 1; i < (stride>>1); i++)
  {
    AE_L16X4_XP(x2, px0, stride * sizeof(complex_fract16));
    AE_L16X4_XP(x1, px0, stride * sizeof(complex_fract16));
    AE_L16X4_XP(x0, px0, (-2*stride - 2) * sizeof(complex_fract16));

    x0 = AE_SRAI16R(x0, shift);
    x1 = AE_SRAI16R(x1, shift);
    x2 = AE_SRAI16R(x2, shift);
    DFT3X2(x0, x1, x2);
    x20 = AE_SEXT32X2D16_32(x0);
    x10 = AE_SEXT32X2D16_10(x0);
    AE_MUL16X4(x21, x11, x1, c1);
    AE_MUL16X4(x22, x12, x2, c1);

    AE_L16X4_IP(tw112, ptwd, sizeof(ae_int16x4));
    AE_L16X4_IP(tw012, ptwd, sizeof(ae_int16x4));
    x11 = AE_MULFC32X16RAS_H(x11, tw112);
    x12 = AE_MULFC32X16RAS_L(x12, tw112);
    x21 = AE_MULFC32X16RAS_H(x21, tw012);
    x22 = AE_MULFC32X16RAS_L(x22, tw012);

    x0 = AE_SAT16X4(x00,x01);
    x1 = AE_SAT16X4(x02,x10);
    x2 = AE_SAT16X4(x11,x12);
    AE_S16X4_IP(x0, py0, sizeof(ae_int16x4));
    AE_S16X4_IP(x1, py0, sizeof(ae_int16x4));
    AE_S16X4_IP(x2, py0, sizeof(ae_int16x4));

    x00 = x20;
    x01 = x21;
    x02 = x22;
  }
  /* Last butterfly radix 3 */
  px0s = (ae_int32 *)x + 1;
  {
    AE_L32_XP(x12, px0s, stride * sizeof(complex_fract16));
    AE_L32_XP(x11, px0s, stride * sizeof(complex_fract16));
    x10 = AE_L32_I(px0s, 0);
    x0 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x10));
    x1 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x11));
    x2 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x12));

    x0 = AE_SRAI16R(x0, shift);
    x1 = AE_SRAI16R(x1, shift);
    x2 = AE_SRAI16R(x2, shift);
    DFT3X2(x0, x1, x2);
    x10 = AE_SEXT32X2D16_10(x0);
    x11 = AE_SEXT32X2D16_10(x1);
    x12 = AE_SEXT32X2D16_10(x2);

    AE_L16X4_IP(tw112, ptwd, sizeof(ae_int16x4));
    x11 = AE_MULFC32X16RAS_H(x11, tw112);
    x12 = AE_MULFC32X16RAS_L(x12, tw112);

    x0 = AE_SAT16X4(x00,x01);
    x1 = AE_SAT16X4(x02,x10);
    x2 = AE_SAT16X4(x11,x12);
    AE_S16X4_IP(x0, py0, sizeof(ae_int16x4));
    AE_S16X4_IP(x1, py0, sizeof(ae_int16x4));
    AE_S16X4_IP(x2, py0, sizeof(ae_int16x4));
  }

  *v *= 3;
  return shift;
  #undef shift
} /* ifft_16x16_stage_first_scl3_DFT3() */

#endif
