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
 *  32x16 FFT first stage Radix 3, scalingOption=3
 */
int fft_32x16_stage_first_scl3_DFT3(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  int i, shift;
  const ae_int32x2 * restrict px0;
  const ae_int32x2 * restrict px1;
  const ae_int32x2 * restrict px2;
        ae_int32x2 * restrict py0;
  const ae_int16x4 * restrict ptwd;
  const int R = 3; // stage radix
  const int stride = N / R;
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  px0 = (ae_int32x2 *)x;
  px1 = px0 + stride;
  px2 = px1 + stride;
  py0 = (ae_int32x2 *)y;
  ptwd = (const ae_int16x4 *)tw;

  shift = 3;
  WUR_AE_SAR(shift);

  __Pragma("loop_count min=2");
  for (i = 0; i <stride; i++)
  {
    ae_int32x2 x0, x1, x2;
    ae_int16x4 tw12;

    AE_L32X2_IP(x0, px0, sizeof(ae_int32x2));
    AE_L32X2_IP(x1, px1, sizeof(ae_int32x2));
    AE_L32X2_IP(x2, px2, sizeof(ae_int32x2));

    AE_L16X4_IP(tw12, ptwd, sizeof(ae_int16x4));

    DFT3X1_RNG_shift(x0, x1, x2, 3);
    x1 = AE_MULFC32X16RAS_H(x1, tw12);
    x2 = AE_MULFC32X16RAS_L(x2, tw12);

    AE_S32X2_IP(x0, py0, sizeof(ae_int32x2));
    AE_S32X2_IP(x1, py0, sizeof(ae_int32x2));
    AE_S32X2_IP(x2, py0, sizeof(ae_int32x2));
  }

  *v = *v * R;
  return shift;
} /* fft_32x16_stage_first_scl3_DFT3() */
#endif
