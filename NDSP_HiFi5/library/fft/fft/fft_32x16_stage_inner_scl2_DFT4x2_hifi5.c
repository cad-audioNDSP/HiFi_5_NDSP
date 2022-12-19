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



#define DFT4X1(x0, x1, x2, x3)\
{   \
\
    ae_int32x2 s0, s1, d0, d1;       \
    AE_ADDANDSUB32S(s0, d0, x0, x2); \
    AE_ADDANDSUB32S(s1, d1, x1, x3); \
    d1 = AE_MUL32JS(d1);             \
    AE_ADDANDSUB32S(x0, x2, s0, s1); \
    AE_ADDANDSUB32S(x3, x1, d0, d1); \
}
/* radix-4 butterfly with normalization */
#define DFT4X1RNG(x0, x1, x2, x3)\
{   \
    ae_int32x2 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG32(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32(s1, d1, x1, x3); \
    d1 = AE_MUL32JS(d1);               \
    AE_ADDANDSUB32S(x0, x2, s0, s1);   \
    AE_ADDANDSUB32S(x3, x1, d0, d1);   \
}

/*
 *  32x16 FFT/IFFT intermediate stage Radix 4 unrolled 2 times, scalingOption=2
 */
int fft_32x16_stage_inner_scl2_DFT4x2(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  int i, j, _v = *v;
  const ae_int32x2 * restrict px0; 
        ae_int32x2 * restrict py0;
  const ae_int16x4 * restrict ptwd;
  int shift;
  const int min_shift = 3;
  const int R = 4; /* stage radix */
  const int stride = (N >> 2);
  int ninner = stride / _v;
  ptwd = (const ae_int16x4 *)tw;

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  shift = XT_MAX(0, (min_shift - *bexp));


  WUR_AE_SAR(shift);
  if (ninner > _v)
  {
      __Pragma("loop_count min=1");
      for (j = 0; j < (_v >> 1); j++)
      {
          ae_int32x2 x00, x10, x20, x30;
          ae_int32x2 x01, x11, x21, x31;

          ae_int16x4 tw0102, tw0311, tw1213;
          ptwd = (const ae_int16x4 *)tw;
          px0 = (ae_int32x2 *)x + j * 2;
          py0 = (ae_int32x2 *)y + j * 2;

          __Pragma("loop_count min=1");
          for (i = 0; i < (ninner >> 1); i++)
          {
              /* 12 cycles per pipeline stage in steady state with unroll=1 */
              AE_L16X4_IP(tw0102, ptwd, sizeof(ae_int16x4));
              AE_L16X4_IP(tw0311, ptwd, sizeof(ae_int16x4));
              AE_L16X4_IP(tw1213, ptwd, sizeof(ae_int16x4));

              AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
              AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
              AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
              AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (_v - 3 * stride)* sizeof(ae_int32x2));

              DFT4X1RNG(x00, x10, x20, x30);
              DFT4X1RNG(x01, x11, x21, x31);
              x10 = AE_MULFC32X16RAS_H(x10, tw0102);
              x20 = AE_MULFC32X16RAS_L(x20, tw0102);
              x30 = AE_MULFC32X16RAS_H(x30, tw0311);
              x11 = AE_MULFC32X16RAS_H(x11, tw0102);
              x21 = AE_MULFC32X16RAS_L(x21, tw0102);
              x31 = AE_MULFC32X16RAS_H(x31, tw0311);

              AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
              AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
              AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
              AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));

              AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
              AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
              AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
              AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (_v - 3 * stride)* sizeof(ae_int32x2));

              DFT4X1RNG(x00, x10, x20, x30);
              DFT4X1RNG(x01, x11, x21, x31);
              x10 = AE_MULFC32X16RAS_L(x10, tw0311);
              x20 = AE_MULFC32X16RAS_H(x20, tw1213);
              x30 = AE_MULFC32X16RAS_L(x30, tw1213);
              x11 = AE_MULFC32X16RAS_L(x11, tw0311);
              x21 = AE_MULFC32X16RAS_H(x21, tw1213);
              x31 = AE_MULFC32X16RAS_L(x31, tw1213);

              AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
              AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
              AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
              AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
          }
          if (ninner & 1)
          {
              AE_L16X4_IP(tw0102, ptwd, sizeof(ae_int16x4));
              tw0311 = AE_SHORTSWAP(AE_MOVINT16X4_FROMF32X2(AE_L32_I((ae_int32*)ptwd, 0)));

              AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
              AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
              AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
              AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (_v - 3 * stride)* sizeof(ae_int32x2));

              DFT4X1RNG(x00, x10, x20, x30);
              DFT4X1RNG(x01, x11, x21, x31);
              x10 = AE_MULFC32X16RAS_H(x10, tw0102);
              x20 = AE_MULFC32X16RAS_L(x20, tw0102);
              x30 = AE_MULFC32X16RAS_H(x30, tw0311);
              x11 = AE_MULFC32X16RAS_H(x11, tw0102);
              x21 = AE_MULFC32X16RAS_L(x21, tw0102);
              x31 = AE_MULFC32X16RAS_H(x31, tw0311);

              AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
              AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
              AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
              AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
          }
      }
  }
  else
  {
      __Pragma("loop_count min=1");
      for (i = 0; i < ninner; i++)
      {
          ae_int32x2 x00, x10, x20, x30;
          ae_int32x2 x01, x11, x21, x31;
          ae_int32x2 t1, t2, t3;
          ae_int16x4 tw1, tw2, tw3;

          AE_L32_IP(t1, castxcc(ae_int32, ptwd), 4);
          AE_L32_IP(t2, castxcc(ae_int32, ptwd), 4);
          AE_L32_IP(t3, castxcc(ae_int32, ptwd), 4);

          tw1 = AE_SHORTSWAP(AE_MOVINT16X4_FROMF32X2(t1));
          tw2 = AE_SHORTSWAP(AE_MOVINT16X4_FROMF32X2(t2));
          tw3 = AE_SHORTSWAP(AE_MOVINT16X4_FROMF32X2(t3));

          px0 = (ae_int32x2 *)x + i*_v;
          py0 = (ae_int32x2 *)y + i * 4 * _v;

          __Pragma("loop_count min=1");
          for (j = 0; j < (_v >> 1); j++)
          {
              /* 6(?) cycles per pipeline stage in steady state with unroll=1 */
              //      px0 = (ae_int32x2 *)x + j * 2 + i*_v;
              //      py0 = (ae_int32x2 *)y + j * 2 + i*4*_v;

              AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
              AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
              AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
              AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (2 - 3 * stride)* sizeof(ae_int32x2));

              DFT4X1RNG(x00, x10, x20, x30);
              DFT4X1RNG(x01, x11, x21, x31);
              x10 = AE_MULFC32X16RAS_H(x10, tw1);
              x20 = AE_MULFC32X16RAS_H(x20, tw2);
              x30 = AE_MULFC32X16RAS_H(x30, tw3);
              x11 = AE_MULFC32X16RAS_H(x11, tw1);
              x21 = AE_MULFC32X16RAS_H(x21, tw2);
              x31 = AE_MULFC32X16RAS_H(x31, tw3);

              AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
              AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
              AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
              AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), (2 - 3 * _v) * sizeof(ae_int32x2));
          }  /* for (j = 0; j < (_v >> 1); j++) */
      } /* for (i = 0; i < ninner; i++) */
  }

    AE_CALCRNG3();
    *bexp = 3 - RUR_AE_SAR();
    *v *= R;
    return shift;
  
} /* fft_32x16_stage_inner_scl2_DFT4x2() */
