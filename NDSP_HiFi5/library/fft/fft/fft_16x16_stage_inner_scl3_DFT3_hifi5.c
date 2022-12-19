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
#include "NatureDSP_types.h"
#include "common.h"
#include "fft_16x16_stages.h"

/* Assumed that SAR is 0x204*/
inline_  void __DFT3xI2_s2(ae_int16x4 *x0, ae_int16x4 *x1, ae_int16x4 *x2, const ae_int16x4 *c, int shift)
{
    ae_int16x4 s0, s1, d0;

    s0 = *x1;
    d0 = *x2;
    *x0 = AE_SRAI16R(*x0, 2);
    AE_ADDANDSUBRNG16RAS_S2(s0, d0);

    s1 = AE_ADD16S(*x0, s0);
    s0 = AE_SRAI16(s0, 1);
    s0 = AE_SUB16S(*x0, s0);
    d0 = AE_MULFC16RAS(*c, d0);
    *x0 = s1;
    
    AE_ADDANDSUBRNG16RAS_S1(s0, d0);
    *x1 = d0;
    *x2 = s0;
}
#define DFT3XI2_S2(__x0, __x1, __x2, __c, __s) __DFT3xI2_s2(&__x0, &__x1, &__x2, (const ae_int16x4*)&__c, __s)

/*
 *  Intermediate stage of FFT/IFFT 16x16, radix-3, static scaling
 */
int fft_16x16_stage_inner_scl3_DFT3(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{ 
  int i, j, _v;
  const int N3=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),715827883));//N/3
  const int stride = N3;
  const ae_int16x4 * restrict px0;
        ae_int16x4 * restrict py0;
  const ae_int16x4 * restrict ptwd;
  ae_int16x4 x00, x10, x20;
  ae_int16x4 x01, x11, x21;
  const int s2 = 2;
  const int shift = s2;

  ae_int16x4 sel = AE_MOVINT16X4_FROMINT64(0x0705060403010200);

  WUR_AE_SAR(s2 * 0x102);
  /* Internal twiddle for DFT3XI2_S2 */
  ae_int16x4  c = AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32(0x00006EDA));

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  _v = *v;
  NASSERT((_v&1)==0);
  __Pragma("loop_count min=1");
  for (j = 0; j < (_v>>1); j+=2)
  {
    ptwd = (const ae_int16x4 *)tw;
    px0 = (const ae_int16x4 *)x + j;
    py0 = (      ae_int16x4 *)y + j;
    /* 8 cycles per pipeline stage in steady state with unroll=1 */
    __Pragma("loop_count min=1");
    for (i = 0; i < (stride/ _v); i++)
    {
      ae_int16x4 tw1, tw2, tmp;

      AE_L16X4_IP(tmp, ptwd, sizeof(tmp)); 
      AE_DSEL16X4(tw1, tw2, tmp, tmp, sel);

      AE_L16X4X2_XP(x00, x01, castxcc(ae_int16x8, px0), stride*sizeof(complex_fract16));
      AE_L16X4X2_XP(x10, x11, castxcc(ae_int16x8, px0), stride*sizeof(complex_fract16));
      AE_L16X4X2_XP(x20, x21, castxcc(ae_int16x8, px0), (_v - 2 * stride)*sizeof(complex_fract16));

      DFT3XI2_S2(x00, x10, x20, c, s2);
      DFT3XI2_S2(x01, x11, x21, c, s2);

      x10 = AE_MULFC16RAS(x10, tw1);
      x20 = AE_MULFC16RAS(x20, tw2);
      x11 = AE_MULFC16RAS(x11, tw1);
      x21 = AE_MULFC16RAS(x21, tw2);

      AE_S16X4X2_XP(x00, x01, castxcc(ae_int16x8, py0), _v * sizeof(complex_fract16));
      AE_S16X4X2_XP(x10, x11, castxcc(ae_int16x8, py0), _v * sizeof(complex_fract16));
      AE_S16X4X2_XP(x20, x21, castxcc(ae_int16x8, py0), _v * sizeof(complex_fract16));
    }
  }
  *v *= 3;
  return shift;
} /* fft_16x16_stage_inner_scl3_DFT3() */
