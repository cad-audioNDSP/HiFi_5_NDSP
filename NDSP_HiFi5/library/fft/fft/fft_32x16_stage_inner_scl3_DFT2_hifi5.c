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

/*
 *  32x16 FFT/IFFT intermediate stage Radix 2, scalingOption=3
 */
int fft_32x16_stage_inner_scl3_DFT2(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{ 
  int i, j;
  int shift, _v;
  ae_int32x2 * restrict px0;
  ae_int32x2 * restrict py0;
  const ae_int32 * restrict ptwd;
  const int R = 2; // stage radix
  const int stride = (N >> 1);

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  shift = 1;
  WUR_AE_SAR(shift);
  _v = *v;
  for (j = 0; j < _v; j++)
  {
    ae_int32x2 x0, x1, y0, y1;
    ptwd = (const ae_int32 *)tw;
    px0 = (ae_int32x2 *)x + j;
    py0 = (ae_int32x2 *)y + j;

    __Pragma("loop_count min=2");
    for (i = 0; i < (stride / _v); i++)
    {
      ae_int16x4 tw1;
      ae_int32x2 tw1_tmp;

      AE_L32_XP(tw1_tmp, ptwd, tw_step*2*sizeof(int16_t));
      tw1 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(tw1_tmp));

      x1 = AE_L32X2_X(px0, sizeof(ae_int32x2) * stride);
      AE_L32X2_XP(x0, px0, _v * sizeof(ae_int32x2));

      AE_ADDANDSUBRNG32(y0, y1, x0, x1);
      y1 = AE_MULFC32X16RAS_L(y1, tw1);

      AE_S32X2_XP(y0, py0, _v * sizeof(ae_int32x2));
      AE_S32X2_XP(y1, py0, _v * sizeof(ae_int32x2));

    }
  }

  *v *= R;
  return shift;
} /* fft_32x16_stage_inner_scl3_DFT2() */
