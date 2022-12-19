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
 *  32x16 FFT/IFFT intermediate stage Radix 4, scalingOption=3, v=4
 */
int fft_32x16_stage_inner_scl3_DFT4x2_v4(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  const int shift = 2;
  //const int R = 4; /* stage radix */
  const int stride = (N >> 2);
  int i;
  int ninner = stride >> 2; 
  const ae_int32x2 * restrict px0 = (ae_int32x2 *)x + 0 * 2; 
        ae_int32x2 * restrict py0 = (ae_int32x2 *)y + 0 * 2;
  const ae_int32x2 * restrict px1 = (ae_int32x2 *)x + 1 * 2; 
        ae_int32x2 * restrict py1 = (ae_int32x2 *)y + 1 * 2;
  const ae_int16x4 * restrict ptwd = (const ae_int16x4 *)tw;

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT(v[0] == 4);
  WUR_AE_SAR(shift>>1);
  
  ae_int32x2 x00, x10, x20, x30;
  ae_int32x2 x01, x11, x21, x31;

  ae_int16x4 tw0102, tw0311, tw1213;

  __Pragma("loop_count min=1");
  for (i = 0; i < (ninner >> 1); i++)
  {
     
      AE_L16X4_IP(tw0102, ptwd, sizeof(ae_int16x4));
      AE_L16X4_IP(tw0311, ptwd, sizeof(ae_int16x4));
      AE_L16X4_IP(tw1213, ptwd, sizeof(ae_int16x4));

      AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (4 - 3 * stride)* sizeof(ae_int32x2));

      DFT4X1RNG2(x00, x10, x20, x30);
      DFT4X1RNG2(x01, x11, x21, x31);
      x10 = AE_MULFC32X16RAS_H(x10, tw0102);
      x20 = AE_MULFC32X16RAS_L(x20, tw0102);
      x30 = AE_MULFC32X16RAS_H(x30, tw0311);
      x11 = AE_MULFC32X16RAS_H(x11, tw0102);
      x21 = AE_MULFC32X16RAS_L(x21, tw0102);
      x31 = AE_MULFC32X16RAS_H(x31, tw0311);

      AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));

      AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (4 - 3 * stride)* sizeof(ae_int32x2));

      DFT4X1RNG2(x00, x10, x20, x30);
      DFT4X1RNG2(x01, x11, x21, x31);
      x10 = AE_MULFC32X16RAS_L(x10, tw0311);
      x20 = AE_MULFC32X16RAS_H(x20, tw1213);
      x30 = AE_MULFC32X16RAS_L(x30, tw1213);
      x11 = AE_MULFC32X16RAS_L(x11, tw0311);
      x21 = AE_MULFC32X16RAS_H(x21, tw1213);
      x31 = AE_MULFC32X16RAS_L(x31, tw1213);

      AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));

      /////////////////////

      AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px1), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px1), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px1), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px1), (4 - 3 * stride)* sizeof(ae_int32x2));

      DFT4X1RNG2(x00, x10, x20, x30);
      DFT4X1RNG2(x01, x11, x21, x31);
      x10 = AE_MULFC32X16RAS_H(x10, tw0102);
      x20 = AE_MULFC32X16RAS_L(x20, tw0102);
      x30 = AE_MULFC32X16RAS_H(x30, tw0311);
      x11 = AE_MULFC32X16RAS_H(x11, tw0102);
      x21 = AE_MULFC32X16RAS_L(x21, tw0102);
      x31 = AE_MULFC32X16RAS_H(x31, tw0311);

      AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));

      AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px1), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px1), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px1), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px1), (4 - 3 * stride)* sizeof(ae_int32x2));

      DFT4X1RNG2(x00, x10, x20, x30);
      DFT4X1RNG2(x01, x11, x21, x31);
      x10 = AE_MULFC32X16RAS_L(x10, tw0311);
      x20 = AE_MULFC32X16RAS_H(x20, tw1213);
      x30 = AE_MULFC32X16RAS_L(x30, tw1213);
      x11 = AE_MULFC32X16RAS_L(x11, tw0311);
      x21 = AE_MULFC32X16RAS_H(x21, tw1213);
      x31 = AE_MULFC32X16RAS_L(x31, tw1213);

      AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
  } /* for (i = 0; i < (ninner >> 1); i++) */
  if (ninner & 1)
  {
      AE_L16X4_IP(tw0102, ptwd, sizeof(ae_int16x4));
      tw0311 = AE_SHORTSWAP(AE_MOVINT16X4_FROMF32X2(AE_L32_I((ae_int32*)ptwd, 0)));

      AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (4 - 3 * stride)* sizeof(ae_int32x2));

      DFT4X1RNG2(x00, x10, x20, x30);
      DFT4X1RNG2(x01, x11, x21, x31);
      x10 = AE_MULFC32X16RAS_H(x10, tw0102);
      x20 = AE_MULFC32X16RAS_L(x20, tw0102);
      x30 = AE_MULFC32X16RAS_H(x30, tw0311);
      x11 = AE_MULFC32X16RAS_H(x11, tw0102);
      x21 = AE_MULFC32X16RAS_L(x21, tw0102);
      x31 = AE_MULFC32X16RAS_H(x31, tw0311);

      AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));

      ///////

      AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px1), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px1), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px1), stride * sizeof(ae_int32x2));
      AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px1), (4 - 3 * stride)* sizeof(ae_int32x2));

      DFT4X1RNG2(x00, x10, x20, x30);
      DFT4X1RNG2(x01, x11, x21, x31);
      x10 = AE_MULFC32X16RAS_H(x10, tw0102);
      x20 = AE_MULFC32X16RAS_L(x20, tw0102);
      x30 = AE_MULFC32X16RAS_H(x30, tw0311);
      x11 = AE_MULFC32X16RAS_H(x11, tw0102);
      x21 = AE_MULFC32X16RAS_L(x21, tw0102);
      x31 = AE_MULFC32X16RAS_H(x31, tw0311);

      AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
      AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
  }

    *v = 16;
    return shift;
} /* fft_32x16_stage_inner_scl3_DFT4x2() */
