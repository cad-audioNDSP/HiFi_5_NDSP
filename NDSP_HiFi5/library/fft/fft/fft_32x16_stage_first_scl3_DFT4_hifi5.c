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

/* radix-4 butterfly with normalization */
#define DFT4X1RNG(x0, x1, x2, x3)\
{   \
    ae_int32x2 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG32_L(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32_L(s1, d1, x1, x3); \
    AE_ADDANDSUB32S(x0, x2, s0, s1);   \
    AE_ADDANDSUB32JS(x3, x1, d0, d1);   \
}
/*
 *  32x16 FFT first stage Radix 4, scalingOption=3
 */
int fft_32x16_stage_first_scl3_DFT4(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  int i;
  const ae_int32x2 * restrict px0;
  const ae_int32x2 * restrict px1;
  const ae_int32x2 * restrict px2;
  const ae_int32x2 * restrict px3;
        ae_int32x2 * restrict py0;
  const ae_int16x4 * restrict ptwd;
  const int R = 4; // stage radix
  const int stride = (N >> 2);
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT(tw_step == 1);
  NASSERT(0==(N&7)); /* N must be a multiple of 8 */

  px0 = (const ae_int32x2 *)x;
  px1 = (const ae_int32x2 *)px0 + stride;
  px2 = (const ae_int32x2 *)px1 + stride;
  px3 = (const ae_int32x2 *)px2 + stride;
  py0 = (      ae_int32x2 *)y;

  ptwd = (const ae_int16x4 *)tw;
  WUR_AE_SAR(3);
  __Pragma("loop_count min=2");
  for (i = 0; i < (stride>>1); i++)
  {
      /* 13 cycles per pipeline stage in steady state with unroll=2  */
      ae_int32x2 x00, x10, x20, x30;
      ae_int32x2 x01, x11, x21, x31;
    ae_int16x4 tw0102, tw0311, tw1213;

    AE_L16X4_IP(tw0102, ptwd, sizeof(ae_int16x4));
    AE_L16X4_IP(tw0311, ptwd, sizeof(ae_int16x4));
    AE_L16X4_IP(tw1213, ptwd, sizeof(ae_int16x4));

    AE_L32X2X2_IP(x00, x01, castxcc(ae_int32x4, px0), 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(x10, x11, castxcc(ae_int32x4, px1), 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(x20, x21, castxcc(ae_int32x4, px2), 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(x30, x31, castxcc(ae_int32x4, px3), 2*sizeof(ae_int32x2));


    DFT4X1RNG(x00, x10, x20, x30);
    DFT4X1RNG(x01, x11, x21, x31);

    x10 = AE_MULFC32X16RAS_H(x10, tw0102);
    x20 = AE_MULFC32X16RAS_H(x20, tw0311);
    x30 = AE_MULFC32X16RAS_H(x30, tw1213);
    x11 = AE_MULFC32X16RAS_L(x11, tw0102);
    x21 = AE_MULFC32X16RAS_L(x21, tw0311);
    x31 = AE_MULFC32X16RAS_L(x31, tw1213);

    AE_S32X2X2_IP(x00, x10, castxcc(ae_int32x4, py0), 2 * sizeof(ae_int32x2));
    AE_S32X2X2_IP(x20, x30, castxcc(ae_int32x4, py0), 2 * sizeof(ae_int32x2));
    AE_S32X2X2_IP(x01, x11, castxcc(ae_int32x4, py0), 2 * sizeof(ae_int32x2));
    AE_S32X2X2_IP(x21, x31, castxcc(ae_int32x4, py0), 2 * sizeof(ae_int32x2));
  }

  *v = v[0] * R;
  return 3;
} /* fft_32x16_stage_first_scl3_DFT4() */
