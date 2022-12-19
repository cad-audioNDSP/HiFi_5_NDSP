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
#include "NatureDSP_types.h"
#include "common.h"
#include "fft_32x16_stages.h"


/* radix-3 butterfly with data scaling using AE.SAR register & 'shift' value */
#define DFT3X1_RNG_shift(x0, x1, x2, shift)\
{\
    ae_int32x2 s0, s1, d0, c32; \
    ae_int16x4 c;               \
    c32 = AE_MOVDA32X2(0x40004000,0x00006EDA); \
    c = AE_MOVINT16X4_FROMINT32X2(c32);        \
    x0 = AE_SRAI32(x0, shift);                 \
    AE_ADDANDSUBRNG32(s0, d0, x1, x2);         \
    s1 = AE_ADD32S(x0, s0);                    \
    AE_MULSFP32X16X2RAS_H(x0, s0, c);          \
    d0 = AE_MULFC32X16RAS_L(d0, c);            \
    s0 = x0;                                   \
    x0 = s1;                                   \
    AE_ADDANDSUB32S(x2, x1, s0, d0);           \
}

/*
 *  32x16 FFT/IFFT last stage Radix 3, scalingOption=3
 */
int fft_32x16_stage_last_scl3_DFT3(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
#define shift 2
  ae_int32x2 * restrict px0;
  ae_int32x2 * restrict py0;
  const int R = 3; // stage radix
  const int stride = N / R;
  int j, _v;

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  _v = *v;
  px0 = (ae_int32x2 *)x;
  py0 = (ae_int32x2 *)y;

  WUR_AE_SAR(shift);

  __Pragma("loop_count min=1");
  for (j = 0; j < (_v >> 1); j++)
  {
      /*  11 cycles per pipeline stage in steady state with unroll=2 */
      ae_int32x2 x00, x10, x20;
      ae_int32x2 x01, x11, x21;

      AE_L32X2X2_X(x10, x11, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride);
      AE_L32X2X2_X(x20, x21, (ae_int32x4*)px0, 2 * sizeof(ae_int32x2)* stride);
      AE_L32X2X2_IP(x00, x01, castxcc(ae_int32x4, px0), 2 * sizeof(ae_int32x2));

      DFT3X1_RNG_shift(x00, x10, x20, shift);
      DFT3X1_RNG_shift(x01, x11, x21, shift);

      AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), (-2 * stride + 2) * sizeof(ae_int32x2));
  }

  return shift;
#undef shift
} /* fft_32x16_stage_last_scl3_DFT3() */
