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
 *  32x16 FFT/IFFT last stage Radix 2, scalingOption=3
 */
int fft_32x16_stage_last_scl3_DFT2(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  int j;
  const ae_int32x2 * restrict px0;
  const ae_int32x2 * restrict px1;
        ae_int32x2 * restrict py0;
        ae_int32x2 * restrict py1;
  const int stride = (N >> 1);
  int shift;

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  shift = 1;
  WUR_AE_SAR(shift);
  px0 = (const ae_int32x2 *)x;
  px1 = (const ae_int32x2 *)x+stride;
  py0 = (      ae_int32x2 *)y;
  py1 = (      ae_int32x2 *)y+stride;

  __Pragma("loop_count min=3");
  for (j = 0; j < stride; j++)
  {
    ae_int32x2 x0, x1, y0, y1;

    AE_L32X2_IP(x0, px0, sizeof(ae_int32x2));
    AE_L32X2_IP(x1, px1, sizeof(ae_int32x2));

    AE_ADDANDSUBRNG32(y0, y1, x0, x1);

    AE_S32X2_IP(y0, py0, sizeof(ae_int32x2));
    AE_S32X2_IP(y1, py1, sizeof(ae_int32x2));
  }

  return shift;
} /* fft_32x16_stage_last_scl3_DFT2() */

#endif
