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
    with dynamic data scaling: 32-bit data, 16-bit twiddle factors
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/
#if 0 /* not used in library for now */
#include "NatureDSP_types.h"
#include "common.h"
#include "fft_32x16_stages.h"

/*
 *  32x16 IFFT first stage Radix 2, scalingOption=2
 */
int ifft_32x16_stage_first_scl2_DFT2(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{ 
  int i;
  int shift, shiftl, shiftr;
        ae_int32x2 * restrict px0;
        ae_int32x2 * restrict py0;
  const ae_int32   * restrict ptwd;
  ae_int32x2 scl;
  ae_int32x2 x0, x1, y0, y1;
  const int R = 2; // stage radix
  const int stride = (N >> 1);
  const int min_shift = 2;

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  shift = min_shift - *bexp;
  shiftl = XT_MAX(0, -shift);
  shiftr = XT_MAX(0,  shift);
  NASSERT(shift>-32 && shift<4);

  scl = 1 << shiftl;
  WUR_AE_SAR(shiftr);

  ptwd = (const ae_int32 *)tw;
  px0 = (ae_int32x2*)x;
  py0 = (ae_int32x2 *)(y);
  {
    /* First butterfly radix 2 */
    ae_int16x4 tw1;
    ae_int32x2 tw1_tmp;

    AE_L32_IP(tw1_tmp, ptwd, 2*sizeof(int16_t));
    tw1 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(tw1_tmp));

    x0 = AE_L32X2_I(px0, 0);
    x1 = AE_L32X2_X(px0, (N-stride)*sizeof(ae_int32x2));

    x0 = AE_MULP32X2(x0, scl);
    x1 = AE_MULP32X2(x1, scl);
    AE_ADDANDSUBRNG32(y0, y1, x0, x1);
    y1 = AE_MULFC32X16RAS_L(y1, tw1);

    AE_S32X2RNG_IP(y0, py0, sizeof(ae_int32x2));
    AE_S32X2RNG_IP(y1, py0, sizeof(ae_int32x2));
  }
  px0 = (ae_int32x2*)x + N-stride-1;
  for (i = 1; i < stride; i++)
  { 
    ae_int16x4 tw1;
    ae_int32x2 tw1_tmp;

    AE_L32_IP(tw1_tmp, ptwd, 2*sizeof(int16_t));
    tw1 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(tw1_tmp));

    AE_L32X2_XP(x1, px0, stride * sizeof(ae_int32x2));
    AE_L32X2_XP(x0, px0, (-stride - 1) * sizeof(ae_int32x2));

    x0 = AE_MULP32X2(x0, scl);
    x1 = AE_MULP32X2(x1, scl);
    AE_ADDANDSUBRNG32(y0, y1, x0, x1);
    y1 = AE_MULFC32X16RAS_L(y1, tw1);

    AE_S32X2RNG_IP(y0, py0, sizeof(ae_int32x2));
    AE_S32X2RNG_IP(y1, py0, sizeof(ae_int32x2));
  }

  AE_CALCRNG3();
  *bexp = 3 - RUR_AE_SAR();
  *v *= R;
  return shift;
} /* ifft_32x16_stage_first_scl2_DFT2() */
#endif
