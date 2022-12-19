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
    with dynamic data scaling: 32-bit data, 16-bit twiddle factors
    C code optimized for HiFi5
*/
#include "NatureDSP_types.h"
#include "common.h"
#include "fft_32x16_stages.h"
#include "fft_x16_common.h"

/*
 *  32x16 FFT/IFFT penultimate stage Radix 4 unrolled 4 times, scalingOption=2
 *  Restriction: last stage must be radix 4
 */
int fft_32x16_stage_penultimate_scl2_DFT4x4(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
#if 1
    return fft_32x16_stage_inner_scl2_DFT4x2(tw, x, y, N, v, tw_step, bexp); 
#else
  int shift;
  const ae_int32x2 * restrict px0;
  const ae_int32x2 * restrict px1;
  const ae_int32x2 * restrict px2;
  const ae_int32x2 * restrict px3;
        ae_int32x2 * restrict py0;
        ae_int32x2 * restrict py1;
        ae_int32x2 * restrict py2;
        ae_int32x2 * restrict py3;
  const ae_int16x4 * restrict ptwd;
  int j, _v;
  ae_int32x2 x0, x1, x2, x3;
  ae_int16x4 tw0311, tw1213;
  ae_int16x4 tw2122, tw2331;
  ae_int16x4 tw3233;
  const int R = 4; // stage radix
  const int stride = (N >> 2);

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT(tw_step==1);

  shift = 3 - *bexp;
  WUR_AE_SAR( shift );
  _v = *v;
  NASSERT(stride/_v == 4);

  px0 = (ae_int32x2 *)x;
  px1 = px0 + stride;
  px2 = px1 + stride;
  px3 = px2 + stride;
  py0 = (ae_int32x2 *)y;
  py1 = py0 + _v;
  py2 = py1 + _v;
  py3 = py2 + _v;

  ptwd = (const ae_int16x4 *)tw;
  tw0311 = AE_L16X4_I(ptwd, 1*sizeof(ae_int16x4));
  tw1213 = AE_L16X4_I(ptwd, 2*sizeof(ae_int16x4));
  tw2122 = AE_L16X4_I(ptwd, 3*sizeof(ae_int16x4));
  tw2331 = AE_L16X4_I(ptwd, 4*sizeof(ae_int16x4));
  tw3233 = AE_L16X4_I(ptwd, 5*sizeof(ae_int16x4));

  __Pragma("loop_count min=4, factor=4");
  for (j = 0; j < _v; j++)
  {
    AE_L32X2_XP(x0, px0, _v*sizeof(ae_int32x2));
    AE_L32X2_XP(x1, px1, _v*sizeof(ae_int32x2));
    AE_L32X2_XP(x2, px2, _v*sizeof(ae_int32x2));
    AE_L32X2_XP(x3, px3, _v*sizeof(ae_int32x2));
    DFT4X1RNG(x0, x1, x2, x3);
    AE_S32X2RNG_XP(x0, py0, 4*_v*sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x1, py1, 4*_v*sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x2, py2, 4*_v*sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x3, py3, 4*_v*sizeof(ae_int32x2));

    AE_L32X2_XP(x0, px0, _v*sizeof(ae_int32x2));
    AE_L32X2_XP(x1, px1, _v*sizeof(ae_int32x2));
    AE_L32X2_XP(x2, px2, _v*sizeof(ae_int32x2));
    AE_L32X2_XP(x3, px3, _v*sizeof(ae_int32x2));
    DFT4X1RNG(x0, x1, x2, x3);
    x1 = AE_MULFC32X16RAS_L(x1, tw0311);
    x2 = AE_MULFC32X16RAS_H(x2, tw1213);
    x3 = AE_MULFC32X16RAS_L(x3, tw1213);
    AE_S32X2RNG_XP(x0, py0, 4*_v*sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x1, py1, 4*_v*sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x2, py2, 4*_v*sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x3, py3, 4*_v*sizeof(ae_int32x2));

    AE_L32X2_XP(x0, px0, _v*sizeof(ae_int32x2));
    AE_L32X2_XP(x1, px1, _v*sizeof(ae_int32x2));
    AE_L32X2_XP(x2, px2, _v*sizeof(ae_int32x2));
    AE_L32X2_XP(x3, px3, _v*sizeof(ae_int32x2));
    DFT4X1RNG(x0, x1, x2, x3);
    x1 = AE_MULFC32X16RAS_H(x1, tw2122);
    x2 = AE_MULFC32X16RAS_L(x2, tw2122);
    x3 = AE_MULFC32X16RAS_H(x3, tw2331);
    AE_S32X2RNG_XP(x0, py0, 4*_v*sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x1, py1, 4*_v*sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x2, py2, 4*_v*sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x3, py3, 4*_v*sizeof(ae_int32x2));

    AE_L32X2_XP(x0, px0, (1-3*_v)*sizeof(ae_int32x2));
    AE_L32X2_XP(x1, px1, (1-3*_v)*sizeof(ae_int32x2));
    AE_L32X2_XP(x2, px2, (1-3*_v)*sizeof(ae_int32x2));
    AE_L32X2_XP(x3, px3, (1-3*_v)*sizeof(ae_int32x2));
    DFT4X1RNG(x0, x1, x2, x3);
    x1 = AE_MULFC32X16RAS_L(x1, tw2331);
    x2 = AE_MULFC32X16RAS_H(x2, tw3233);
    x3 = AE_MULFC32X16RAS_L(x3, tw3233);
    AE_S32X2RNG_XP(x0, py0, (1-3*4*_v)*sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x1, py1, (1-3*4*_v)*sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x2, py2, (1-3*4*_v)*sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x3, py3, (1-3*4*_v)*sizeof(ae_int32x2));

  }

  AE_CALCRNG3();
  *bexp = 3 - RUR_AE_SAR();
  *v *= R;
  return shift;
  #endif
} /* fft_32x16_stage_penultimate_scl2_DFT4x4() */
