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
    Complex-valued FFT stages with butterflies radix-4, radix-8
    with static data scaling: 32-bit data, 16-bit twiddle factors
    C code optimized for HiFi4
  IntegrIT, 2006-2019
*/
#include "NatureDSP_types.h"
#include "common.h"
#include "fft_32x16_stages.h"
#include "fft_x16_common.h"

/* radix-4 butterfly with normalization applied twice */
#define DFT4X1RNG2(x0, x1, x2, x3)\
{   \
    ae_int32x2 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG32(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32(s1, d1, x1, x3); \
    d1 = AE_MUL32JS(d1);               \
    AE_ADDANDSUBRNG32(x0, x2, s0, s1); \
    AE_ADDANDSUBRNG32(x3, x1, d0, d1); \
}

/*
 *  32x16 FFT/IFFT intermediate stage Radix 4, scalingOption=3
 */
int fft_32x16_stage_inner_scl3_DFT4(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{ 
    int i, j, shift, _v = *v;
        ae_int32x2 * restrict px0; 
        ae_int32x2 * restrict py0;
  const ae_int16x4 * restrict ptwd;
  const int R = 4; // stage radix
  const int stride = (N >> 2);

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT(tw_step==1);
  NASSERT((stride / _v)%2==0); 

  shift = 2;
  WUR_AE_SAR(shift>>1);
  
  for (j = 0; j < _v; j++)
  {
    ae_int32x2 x0, x1, x2, x3;
    ptwd = (const ae_int16x4 *)tw;
    px0 = (ae_int32x2 *)x + j;
    py0 = (ae_int32x2 *)y + j;
    __Pragma("loop_count min=1");
    for (i = 0; i < ((stride/_v)>>1); i++)
    {
        /* 10 cycles per pipeline stage in steady state with unroll=1 */
      ae_int16x4 tw0102, tw0311, tw1213;

      AE_L16X4_IP(tw0102, ptwd, sizeof(ae_int16x4));
      AE_L16X4_IP(tw0311, ptwd, sizeof(ae_int16x4));
      AE_L16X4_IP(tw1213, ptwd, sizeof(ae_int16x4));

      x1 = AE_L32X2_X(px0, sizeof(ae_int32x2) * stride);
      x3 = AE_L32X2_X(px0, sizeof(ae_int32x2) * stride * 3);
      x2 = AE_L32X2_X(px0, sizeof(ae_int32x2) * stride * 2);
      AE_L32X2_XP(x0, px0, _v * sizeof(ae_int32x2));

      DFT4X1RNG2(x0, x1, x2, x3);
      x1 = AE_MULFC32X16RAS_H(x1, tw0102);
      x2 = AE_MULFC32X16RAS_L(x2, tw0102);
      x3 = AE_MULFC32X16RAS_H(x3, tw0311);

      AE_S32X2_XP(x0, py0, _v * sizeof(ae_int32x2));
      AE_S32X2_XP(x1, py0, _v * sizeof(ae_int32x2));
      AE_S32X2_XP(x2, py0, _v * sizeof(ae_int32x2));
      AE_S32X2_XP(x3, py0, _v * sizeof(ae_int32x2));

      x1 = AE_L32X2_X(px0, sizeof(ae_int32x2) * stride);
      x3 = AE_L32X2_X(px0, sizeof(ae_int32x2) * stride * 3);
      x2 = AE_L32X2_X(px0, sizeof(ae_int32x2) * stride * 2);
      AE_L32X2_XP(x0, px0, _v * sizeof(ae_int32x2));

      DFT4X1RNG2(x0, x1, x2, x3);
      x1 = AE_MULFC32X16RAS_L(x1, tw0311);
      x2 = AE_MULFC32X16RAS_H(x2, tw1213);
      x3 = AE_MULFC32X16RAS_L(x3, tw1213);

      AE_S32X2_XP(x0, py0, _v * sizeof(ae_int32x2));
      AE_S32X2_XP(x1, py0, _v * sizeof(ae_int32x2));
      AE_S32X2_XP(x2, py0, _v * sizeof(ae_int32x2));
      AE_S32X2_XP(x3, py0, _v * sizeof(ae_int32x2));
    }
  }

  *v *= R;
  return shift;
} /* fft_32x16_stage_inner_scl3_DFT4() */
