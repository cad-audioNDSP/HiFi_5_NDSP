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

ALIGN(32) static const int16_t __dft5_tw[] =
{
    (int16_t)0x278E, (int16_t)0x278E, (int16_t)0x278E, (int16_t)0x278E,
    (int16_t)0x8644, (int16_t)0x79BC, (int16_t)0x8644, (int16_t)0x79BC,
    (int16_t)0x9872, (int16_t)0x9872, (int16_t)0x9872, (int16_t)0x9872,
    (int16_t)0xB4C3, (int16_t)0x4B3D, (int16_t)0xB4C3, (int16_t)0x4B3D
};
/* twiddles should be loaded from the table above */
#define DFT5X2(x0, x1, x2, x3, x4, w1, w2, w3, w4)\
{ \
    ae_int16x4 s1, s2, d1, d2;             \
    ae_int16x4 t0, t1, t2, t3;             \
    ae_int16x4 y0, y1, y2, y3;             \
    s1 = AE_ADD16S(x1, x4);                \
    s2 = AE_ADD16S(x2, x3);                \
    d1 = AE_SUB16S(x1, x4);                \
    d2 = AE_SUB16S(x2, x3);                \
    \
    t0 = AE_MULFP16X4RAS(s1, w1);         \
    t1 = AE_MULFP16X4RAS(s2, w3);         \
    t2 = AE_MULFP16X4RAS(s1, w3);         \
    t3 = AE_MULFP16X4RAS(s2, w1);         \
    y0 = AE_ADD16S(x0, AE_ADD16S(t0, t1)); \
    y1 = AE_ADD16S(x0, AE_ADD16S(t2, t3)); \
    \
    t0 = AE_MULFP16X4RAS(d1, w2); \
    t1 = AE_MULFP16X4RAS(d2, w4); \
    t2 = AE_MULFP16X4RAS(d2, w2); \
    t3 = AE_MULFP16X4RAS(d1, w4); \
    y2 = AE_ADD16S(t0, t1);      \
    y3 = AE_SUB16S(t3, t2);      \
    y2 = AE_SEL16_2301(y2, y2);  \
    y3 = AE_SEL16_2301(y3, y3);  \
    \
    x0 = AE_ADD16S(x0, AE_ADD16S(s1, s2)); \
    x1 = AE_ADD16S(y0, y2);               \
    x2 = AE_ADD16S(y1, y3);               \
    x3 = AE_SUB16S(y1, y3);               \
    x4 = AE_SUB16S(y0, y2);               \
}

/*
 *  First stage of IFFT 16x16, radix-5, static scaling
 */
int ifft_16x16_stage_first_scl3_DFT5(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
  #define shift 3
  int i;
  const int N5=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),429496730));//N/5
  const int stride = N5;
  const ae_int16x4 * restrict px0;
  const ae_int32   * restrict px0s;
        ae_int16x4 * restrict py0;
  const ae_int16x4 * restrict ptwd;
  const ae_int16x4 * restrict ptwd_dft;
  ae_int16x4 x0, x1, x2, x3, x4;
  ae_int32x2 x00, x01, x02, x03, x04;
  ae_int32x2 x10, x11, x12, x13, x14;
  ae_int32x2 x20, x21, x22, x23, x24;
  ae_int16x4 tw012, tw034, tw112, tw134, c1;
  ae_int16x4 w1, w2, w3, w4;
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  ptwd     = (const ae_int16x4 *)tw;
  ptwd_dft = (const ae_int16x4 *)__dft5_tw;
  py0 = (ae_int16x4 *)y;
  AE_L16X4_IP(w1, ptwd_dft, sizeof(ae_int16x4));
  AE_L16X4_IP(w2, ptwd_dft, sizeof(ae_int16x4));
  AE_L16X4_IP(w3, ptwd_dft, sizeof(ae_int16x4));
  AE_L16X4_XP(w4, ptwd_dft, -3*(int)sizeof(ae_int16x4));
  c1 = 1;
  /*
  ifft(x) = fft( [ x(1), x(end:-1:2)] )
  */
  px0s = (ae_int32 *)x + N - 4*stride;
  {
    /* First butterfly radix 4 */
    AE_L32_XP(x04, px0s, stride*sizeof(ae_int32));
    AE_L32_XP(x03, px0s, stride*sizeof(ae_int32));
    AE_L32_XP(x02, px0s, stride*sizeof(ae_int32));
    x01 = AE_L32_I(px0s, 0);
    x00 = AE_L32_I((ae_int32*)x, 0);
    x0 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x00));
    x1 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x01));
    x2 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x02));
    x3 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x03));
    x4 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x04));
    x0 = AE_SRAI16R(x0, shift);
    x1 = AE_SRAI16R(x1, shift);
    x2 = AE_SRAI16R(x2, shift);
    x3 = AE_SRAI16R(x3, shift);
    x4 = AE_SRAI16R(x4, shift);
    DFT5X2(x0, x1, x2, x3, x4,  w1, w2, w3, w4);
    x00 = AE_SEXT32X2D16_32(x0);
    x01 = AE_SEXT32X2D16_32(x1);
    x02 = AE_SEXT32X2D16_32(x2);
    x03 = AE_SEXT32X2D16_32(x3);
    x04 = AE_SEXT32X2D16_32(x4);

    AE_L16X4_IP(tw012, ptwd, sizeof(ae_int16x4));
    AE_L16X4_IP(tw034, ptwd, sizeof(ae_int16x4));
    x01 = AE_MULFC32X16RAS_H(x01, tw012);
    x02 = AE_MULFC32X16RAS_L(x02, tw012);
    x03 = AE_MULFC32X16RAS_H(x03, tw034);
    x04 = AE_MULFC32X16RAS_L(x04, tw034);
  }
  px0 = (ae_int16x4 *)((complex_fract16 *)x + N - 2 - 4*stride);
  __Pragma("loop_count min=2");
  for (i = 1; i < (stride>>1); i++)
  {
    x0 = AE_L16X4_X(px0, 4*stride*sizeof(complex_fract16));
    x1 = AE_L16X4_X(px0, 3*stride*sizeof(complex_fract16));
    x2 = AE_L16X4_X(px0, 2*stride*sizeof(complex_fract16));
    x3 = AE_L16X4_X(px0, 1*stride*sizeof(complex_fract16));
    AE_L16X4_XP(x4, px0, -2*(int)sizeof(complex_fract16));

    x0 = AE_SRAI16R(x0, shift);
    x1 = AE_SRAI16R(x1, shift);
    x2 = AE_SRAI16R(x2, shift);
    x3 = AE_SRAI16R(x3, shift);
    x4 = AE_SRAI16R(x4, shift);

    DFT5X2(x0, x1, x2, x3, x4,  w1, w2, w3, w4);
    x20 = AE_SEXT32X2D16_32(x0);
    x10 = AE_SEXT32X2D16_10(x0);
    x21 = AE_SEXT32X2D16_32(x1);
    x11 = AE_SEXT32X2D16_10(x1);
    x22 = AE_SEXT32X2D16_32(x2);
    x12 = AE_SEXT32X2D16_10(x2);
    x23 = AE_SEXT32X2D16_32(x3);
    x13 = AE_SEXT32X2D16_10(x3);
    x24 = AE_SEXT32X2D16_32(x4);
    x14 = AE_SEXT32X2D16_10(x4);

    AE_L16X4_IP(tw112, ptwd, sizeof(ae_int16x4));
    AE_L16X4_IP(tw134, ptwd, sizeof(ae_int16x4));
    AE_L16X4_IP(tw012, ptwd, sizeof(ae_int16x4));
    AE_L16X4_IP(tw034, ptwd, sizeof(ae_int16x4));
    x11 = AE_MULFC32X16RAS_H(x11, tw112);
    x12 = AE_MULFC32X16RAS_L(x12, tw112);
    x13 = AE_MULFC32X16RAS_H(x13, tw134);
    x14 = AE_MULFC32X16RAS_L(x14, tw134);
    x21 = AE_MULFC32X16RAS_H(x21, tw012);
    x22 = AE_MULFC32X16RAS_L(x22, tw012);
    x23 = AE_MULFC32X16RAS_H(x23, tw034);
    x24 = AE_MULFC32X16RAS_L(x24, tw034);

    x0 = AE_SAT16X4(x00,x01);
    x1 = AE_SAT16X4(x02,x03);
    x2 = AE_SAT16X4(x04,x10);
    x3 = AE_SAT16X4(x11,x12);
    x4 = AE_SAT16X4(x13,x14);
    AE_S16X4_IP(x0, py0, sizeof(ae_int16x4));
    AE_S16X4_IP(x1, py0, sizeof(ae_int16x4));
    AE_S16X4_IP(x2, py0, sizeof(ae_int16x4));
    AE_S16X4_IP(x3, py0, sizeof(ae_int16x4));
    AE_S16X4_IP(x4, py0, sizeof(ae_int16x4));

    x00 = x20;
    x01 = x21;
    x02 = x22;
    x03 = x23;
    x04 = x24;
  }
  px0s = (ae_int32 *)x + 1;
  {
    AE_L32_XP(x14, px0s, stride*sizeof(ae_int32));
    AE_L32_XP(x13, px0s, stride*sizeof(ae_int32));
    AE_L32_XP(x12, px0s, stride*sizeof(ae_int32));
    AE_L32_XP(x11, px0s, stride*sizeof(ae_int32));
    x10 = AE_L32_I(px0s, 0);
    x0 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x10));
    x1 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x11));
    x2 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x12));
    x3 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x13));
    x4 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(x14));
    x0 = AE_SRAI16R(x0, shift);
    x1 = AE_SRAI16R(x1, shift);
    x2 = AE_SRAI16R(x2, shift);
    x3 = AE_SRAI16R(x3, shift);
    x4 = AE_SRAI16R(x4, shift);
    DFT5X2(x0, x1, x2, x3, x4,  w1, w2, w3, w4);
    x10 = AE_SEXT32X2D16_10(x0);
    x11 = AE_SEXT32X2D16_10(x1);
    x12 = AE_SEXT32X2D16_10(x2);
    x13 = AE_SEXT32X2D16_10(x3);
    x14 = AE_SEXT32X2D16_10(x4);

    AE_L16X4_IP(tw112, ptwd, sizeof(ae_int16x4));
    AE_L16X4_IP(tw134, ptwd, sizeof(ae_int16x4));
    x11 = AE_MULFC32X16RAS_H(x11, tw112);
    x12 = AE_MULFC32X16RAS_L(x12, tw112);
    x13 = AE_MULFC32X16RAS_H(x13, tw134);
    x14 = AE_MULFC32X16RAS_L(x14, tw134);

    x0 = AE_SAT16X4(x00,x01);
    x1 = AE_SAT16X4(x02,x03);
    x2 = AE_SAT16X4(x04,x10);
    x3 = AE_SAT16X4(x11,x12);
    x4 = AE_SAT16X4(x13,x14);
    AE_S16X4_IP(x0, py0, sizeof(ae_int16x4));
    AE_S16X4_IP(x1, py0, sizeof(ae_int16x4));
    AE_S16X4_IP(x2, py0, sizeof(ae_int16x4));
    AE_S16X4_IP(x3, py0, sizeof(ae_int16x4));
    AE_S16X4_IP(x4, py0, sizeof(ae_int16x4));
  }
  *v *= 5;
  return shift;
  #undef shift
} /* ifft_16x16_stage_first_scl3_DFT5() */

#endif
