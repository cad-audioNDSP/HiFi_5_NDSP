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
    with dynamic data scaling: 32-bit data, 32-bit twiddle factors
    C code optimized for HiFi4
  IntegrIT, 2006-2019
*/
#include "NatureDSP_Signal_fft.h"
#include "fft_twiddles32x32.h"
#include "common.h"

inline_ void _cmult32x32(ae_int32x2 *result, ae_int32x2 *x, ae_int32x2 *y)
{
  ae_f32x2 z;

  z = AE_MULFC32RAS(AE_MOVF32X2_FROMINT32X2(*x), AE_MOVF32X2_FROMINT32X2(*y));
  *result = AE_MOVINT32X2_FROMF32X2(z);
}

#define DFT4X1(x0, x1, x2, x3)\
{   \
\
    ae_int32x2 s0, s1, d0, d1;       \
    AE_ADDANDSUB32S(s0, d0, x0, x2); \
    AE_ADDANDSUB32S(s1, d1, x1, x3); \
    AE_ADDANDSUB32S(x0, x2, s0, s1); \
    AE_ADDANDSUB32JS(x3, x1, d0, d1); \
}
/* radix-4 butterfly with normalization */
#define DFT4X1RNG(x0, x1, x2, x3)\
{   \
    ae_int32x2 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG32(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32(s1, d1, x1, x3); \
    AE_ADDANDSUB32S(x0, x2, s0, s1);   \
    AE_ADDANDSUB32JS(x3, x1, d0, d1);   \
}

/*
 *  32x32 FFT first stage Radix 4, scalingOption=2
 */
int fft_stageS2_DFT4_first_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  int i;
  int shift;
  ae_int32x2 * restrict px0;
  ae_int32x2 * restrict px1;
  ae_int32x2 * restrict px2;
  ae_int32x2 * restrict px3;
  ae_int32x2 * restrict py0;
  const ae_int32x2 * restrict ptwd;
  const int R = 4; // stage radix
  const int stride = (N >> 2);
 // const int min_shift = 3;
  shift = 3 - *bexp;
 // shiftl = XT_MAX(0, -shift);
 // shiftr = XT_MAX(0,  shift);

  NASSERT(shift>-32 && shift<4);
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

 // scl = 1 << shiftl;
  WUR_AE_SAR(3);
  px0 = (ae_int32x2 *)x;
  px1 = px0 + stride;
  px2 = px1 + stride;
  px3 = px2 + stride;
  py0 = (ae_int32x2 *)y;
  ptwd = (const ae_int32x2 *)tw;
  if ((stride & 1) != 0 || tw_step != 1)
  {
      __Pragma("loop_count min=3");
      for (i = 0; i < stride; i++)
      {
        ae_int32x2 x0, x1, x2, x3;
        ae_int32x2 tw1, tw2, tw3;

        AE_L32X2_IP(tw1, ptwd, sizeof(ae_int32x2));
        AE_L32X2_IP(tw2, ptwd, sizeof(ae_int32x2));
        AE_L32X2_XP(tw3, ptwd, (3 * tw_step - 2) * sizeof(ae_int32x2));

        AE_L32X2_IP(x0, px0, sizeof(ae_int32x2));
        AE_L32X2_IP(x1, px1, sizeof(ae_int32x2));
        AE_L32X2_IP(x2, px2, sizeof(ae_int32x2));
        AE_L32X2_IP(x3, px3, sizeof(ae_int32x2));

        x0 = AE_SLAA32(x0, *bexp);
        x1 = AE_SLAA32(x1, *bexp);
        x2 = AE_SLAA32(x2, *bexp);
        x3 = AE_SLAA32(x3, *bexp);

        DFT4X1RNG(x0, x1, x2, x3);
        x1 = AE_MULFC32RAS(x1, tw1);
        x2 = AE_MULFC32RAS(x2, tw2);
        x3 = AE_MULFC32RAS(x3, tw3);

        AE_S32X2RNG_IP(x0, py0, sizeof(ae_int32x2));
        AE_S32X2RNG_IP(x1, py0, sizeof(ae_int32x2));
        AE_S32X2RNG_IP(x2, py0, sizeof(ae_int32x2));
        AE_S32X2RNG_IP(x3, py0, sizeof(ae_int32x2));
      }
  }
  else
  {
      __Pragma("loop_count min=1");
      for (i = 0; i < (stride >> 1); i++)
      {
          /* 6 cycles per pipeline stage in steady state with unroll=1
          2 pipeline stages */
          ae_int32x2 x00, x10, x20, x30;
          ae_int32x2 x01, x11, x21, x31;
          ae_int32x2 tw10, tw20, tw30;
          ae_int32x2 tw11, tw21, tw31;

          AE_L32X2X2_IP(tw10, tw20, castxcc(ae_int32x4, ptwd), sizeof(ae_int32x4));
          AE_L32X2X2_IP(tw30, tw11, castxcc(ae_int32x4, ptwd), sizeof(ae_int32x4));
          AE_L32X2X2_IP(tw21, tw31, castxcc(ae_int32x4, ptwd), sizeof(ae_int32x4));

          AE_L32X2X2_IP(x00, x01, castxcc(ae_int32x4, px0), sizeof(ae_int32x4));
          AE_L32X2X2_IP(x10, x11, castxcc(ae_int32x4, px1), sizeof(ae_int32x4));
          AE_L32X2X2_IP(x20, x21, castxcc(ae_int32x4, px2), sizeof(ae_int32x4));
          AE_L32X2X2_IP(x30, x31, castxcc(ae_int32x4, px3), sizeof(ae_int32x4));

          x00 = AE_SLAA32(x00, *bexp);
          x10 = AE_SLAA32(x10, *bexp);
          x20 = AE_SLAA32(x20, *bexp);
          x30 = AE_SLAA32(x30, *bexp);
          x01 = AE_SLAA32(x01, *bexp);
          x11 = AE_SLAA32(x11, *bexp);
          x21 = AE_SLAA32(x21, *bexp);
          x31 = AE_SLAA32(x31, *bexp);

          DFT4X1RNG(x00, x10, x20, x30);
          DFT4X1RNG(x01, x11, x21, x31);

          x10 = AE_MULFC32RAS(x10, tw10);
          x20 = AE_MULFC32RAS(x20, tw20);
          x30 = AE_MULFC32RAS(x30, tw30);
          x11 = AE_MULFC32RAS(x11, tw11);
          x21 = AE_MULFC32RAS(x21, tw21);
          x31 = AE_MULFC32RAS(x31, tw31);

          AE_S32X2X2RNG_IP(x00, x10, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
          AE_S32X2X2RNG_IP(x20, x30, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
          AE_S32X2X2RNG_IP(x01, x11, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
          AE_S32X2X2RNG_IP(x21, x31, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
      } /* for (i = 0; i < (stride>>1); i++) */
  }
  AE_CALCRNG3();
  *bexp = 3 - RUR_AE_SAR();
  *v *= R;
  return shift;

} /* fft_stageS2_DFT4_first_32x32() */

/*
 *  32x32 FFT last stage Radix 4, scalingOption=2
 */
int fft_stageS2_DFT4_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  const ae_int32x2 * restrict px0;
  const ae_int32x2 * restrict px1;
  const ae_int32x2 * restrict px2;
  const ae_int32x2 * restrict px3;
        ae_int32x2 * restrict py0;
        ae_int32x2 * restrict py1;
        ae_int32x2 * restrict py2;
        ae_int32x2 * restrict py3;
  const int min_shift = 3;
  const int stride = (N >> 2);
  int shift;
  int j;

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT(*bexp >= 0);

  shift = (min_shift-*bexp);
  WUR_AE_SAR(shift);
  px0 = (const ae_int32x2 *)x;
  px1 = px0 + stride;
  px2 = px1 + stride;
  px3 = px2 + stride;
  py0 = (ae_int32x2 *)y;
  py1 = py0 + stride;
  py2 = py1 + stride;
  py3 = py2 + stride;

  if (stride & 1)
  {
    __Pragma("loop_count min=3"); 
    for (j = 0; j < stride; j++)
    {
      ae_int32x2 x0, x1, x2, x3;
  
      AE_L32X2_IP(x0, px0, sizeof(ae_int32x2));
      AE_L32X2_IP(x1, px1, sizeof(ae_int32x2));
      AE_L32X2_IP(x2, px2, sizeof(ae_int32x2));
      AE_L32X2_IP(x3, px3, sizeof(ae_int32x2));

      DFT4X1RNG(x0, x1, x2, x3);

      AE_S32X2RNG_IP(x0, py0, sizeof(ae_int32x2));
      AE_S32X2RNG_IP(x1, py1, sizeof(ae_int32x2));
      AE_S32X2RNG_IP(x2, py2, sizeof(ae_int32x2));
      AE_S32X2RNG_IP(x3, py3, sizeof(ae_int32x2));
    }
  }
  else
  {
      ae_int32x2 x0, x1, x2, x3;
      ae_int32x2 x4, x5, x6, x7;
    __Pragma("loop_count min=2"); 
    for (j = 0; j < (stride>>1); j++)
    {
        AE_L32X2X2_XP(x0, x4, castxcc(ae_int32x4, px0), stride*sizeof(complex_fract32));
        AE_L32X2X2_XP(x1, x5, castxcc(ae_int32x4, px0), stride*sizeof(complex_fract32));
        AE_L32X2X2_XP(x2, x6, castxcc(ae_int32x4, px0), stride*sizeof(complex_fract32));
        AE_L32X2X2_XP(x3, x7, castxcc(ae_int32x4, px0), (2 - 3 * stride)*sizeof(complex_fract32));

        DFT4X1RNG(x0, x1, x2, x3);
        DFT4X1RNG(x4, x5, x6, x7);

        AE_S32X2X2RNG_XP(x0, x4, castxcc(ae_int32x4, py0), stride*sizeof(complex_fract32));
        AE_S32X2X2RNG_XP(x1, x5, castxcc(ae_int32x4, py0), stride*sizeof(complex_fract32));
        AE_S32X2X2RNG_XP(x2, x6, castxcc(ae_int32x4, py0), stride*sizeof(complex_fract32));
        AE_S32X2X2RNG_XP(x3, x7, castxcc(ae_int32x4, py0), (2 - 3 * stride)*sizeof(complex_fract32));
    }
  }

  return shift;
} /* fft_stageS2_DFT4_last_32x32() */


/*
*  32x32 IFFT last stage Radix 4, scalingOption=2
*/
int ifft_stageS2_DFT4_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    const ae_int32x2 * restrict px0;
    const ae_int32x2 * restrict px1;
    const ae_int32x2 * restrict px2;
    const ae_int32x2 * restrict px3;
    ae_int32x2 * restrict py0;
    ae_int32x2 * restrict py1;
    ae_int32x2 * restrict py2;
    ae_int32x2 * restrict py3;
    const int min_shift = 3;
    const int stride = (N >> 2);
    int shift;
    int j;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(*bexp >= 0);

    shift = XT_MAX(0, (min_shift - *bexp));
    WUR_AE_SAR(shift);
    px0 = (const ae_int32x2 *)x;
    px1 = px0 + stride;
    px2 = px1 + stride;
    px3 = px2 + stride;
    py0 = (ae_int32x2 *)y;
    py1 = py0 + stride;
    py2 = py1 + stride;
    py3 = py2 + stride;

    if (stride & 1)
    {
        __Pragma("loop_count min=3");
        for (j = 0; j < stride; j++)
        {
            ae_int32x2 x0, x1, x2, x3;

            AE_L32X2_IP(x0, px0, sizeof(ae_int32x2));
            AE_L32X2_IP(x1, px1, sizeof(ae_int32x2));
            AE_L32X2_IP(x2, px2, sizeof(ae_int32x2));
            AE_L32X2_IP(x3, px3, sizeof(ae_int32x2));

            DFT4X1RNG(x0, x1, x2, x3);

            AE_S64_IP(AE_MOVINT64_FROMINT32X2(x0), castxcc(ae_int64, py0), sizeof(ae_int32x2));
            AE_S64_IP(AE_MOVINT64_FROMINT32X2(x1), castxcc(ae_int64, py1), sizeof(ae_int32x2));
            AE_S64_IP(AE_MOVINT64_FROMINT32X2(x2), castxcc(ae_int64, py2), sizeof(ae_int32x2));
            AE_S64_IP(AE_MOVINT64_FROMINT32X2(x3), castxcc(ae_int64, py3), sizeof(ae_int32x2));
        }
    }
    else
    {
        ae_int32x2 x0, x1, x2, x3;
        ae_int32x2 x4, x5, x6, x7;
        __Pragma("loop_count min=2");
        for (j = 0; j < (stride >> 1); j++)
        {
            AE_L32X2X2_XP(x0, x4, castxcc(ae_int32x4, px0), stride*sizeof(complex_fract32));
            AE_L32X2X2_XP(x1, x5, castxcc(ae_int32x4, px0), stride*sizeof(complex_fract32));
            AE_L32X2X2_XP(x2, x6, castxcc(ae_int32x4, px0), stride*sizeof(complex_fract32));
            AE_L32X2X2_XP(x3, x7, castxcc(ae_int32x4, px0), (2 - 3 * stride)*sizeof(complex_fract32));

            DFT4X1RNG(x0, x1, x2, x3);
            DFT4X1RNG(x4, x5, x6, x7);

            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x0), AE_MOVINT64_FROMINT32X2(x4), castxcc(ae_int64x2, py0), stride*sizeof(complex_fract32));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x1), AE_MOVINT64_FROMINT32X2(x5), castxcc(ae_int64x2, py0), stride*sizeof(complex_fract32));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x2), AE_MOVINT64_FROMINT32X2(x6), castxcc(ae_int64x2, py0), stride*sizeof(complex_fract32));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x3), AE_MOVINT64_FROMINT32X2(x7), castxcc(ae_int64x2, py0), (2 - 3 * stride)*sizeof(complex_fract32));
        }
    }

    return shift;
} /* ifft_stageS2_DFT4_last_32x32() */

ALIGN(32) static const int32_t __fft8_tw1[] =
{
    (int32_t)0x2D413CCD, (int32_t)0xD2BEC333,
    (int32_t)0x00000000, (int32_t)0xC0000000,
    (int32_t)0xD2BEC333, (int32_t)0xD2BEC333,
};

/*
 *  32x32 FFT last stage Radix 8, scalingOption=2
 */
int fft_stageS2_DFT8_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    int shift;
    int j;
          ae_int32x2 * restrict px0;
          ae_int32x2 * restrict py0;
    const ae_int32x2 * restrict ptwd;

    const int stride = (N >> 3);
    const int min_shift = 3;
    ae_int32x2 tw1, tw2, tw3;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(*bexp >= 0);

    shift = XT_MAX(0, min_shift - *bexp);
    WUR_AE_SAR(shift);
    px0 = (ae_int32x2 *)x + 7*stride;
    py0 = (ae_int32x2 *)y;
    ptwd = (const ae_int32x2 *)__fft8_tw1;
    tw1 = AE_L32X2_I(ptwd, 0 * sizeof(ae_int32x2));
    tw2 = AE_L32X2_I(ptwd, 1 * sizeof(ae_int32x2));
    tw3 = AE_L32X2_I(ptwd, 2 * sizeof(ae_int32x2));    

    if ((stride & 1) > 0)
    {
        /* Odd loop count, used in the some mixed radix FFT */
        __Pragma("loop_count min=1");
        for (j = 0; j < stride; j++)
        {
            ae_int32x2 x0, x1, x2, x3;
            ae_int32x2 x4, x5, x6, x7;

            AE_L32X2_XP(x7, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x6, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x5, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x4, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x3, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x2, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x1, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x0, px0, sizeof(ae_int32x2)*(7 * stride + 1));

            DFT4X1RNG(x0, x2, x4, x6);
            DFT4X1RNG(x1, x3, x5, x7);

            _cmult32x32(&x3, &x3, &tw1);
            _cmult32x32(&x5, &x5, &tw2);
            _cmult32x32(&x7, &x7, &tw3);

            {
                ae_int32x2 s0, s1, s2, s3;
                ae_int32x2 d0, d1, d2, d3;

                x0 = AE_SRAI32(x0, 1);
                x1 = AE_SRAI32(x1, 1);
                x2 = AE_SRAI32(x2, 1);
                x4 = AE_SRAI32(x4, 1);
                x6 = AE_SRAI32(x6, 1);

                AE_ADDANDSUB32S(s0, d0, x0, x1);
                AE_ADDANDSUB32S(s1, d1, x2, x3);
                AE_ADDANDSUB32S(s2, d2, x4, x5);
                AE_ADDANDSUB32S(s3, d3, x6, x7);

                x0 = s0;        x4 = d0;
                x1 = s1;        x5 = d1;
                x2 = s2;        x6 = d2;
                x3 = s3;        x7 = d3;
            }

            AE_S32X2RNG_XP(x0, py0, stride * sizeof(ae_int32x2));
            AE_S32X2RNG_XP(x1, py0, stride * sizeof(ae_int32x2));
            AE_S32X2RNG_XP(x2, py0, stride * sizeof(ae_int32x2));
            AE_S32X2RNG_XP(x3, py0, stride * sizeof(ae_int32x2));
            AE_S32X2RNG_XP(x4, py0, stride * sizeof(ae_int32x2));
            AE_S32X2RNG_XP(x5, py0, stride * sizeof(ae_int32x2));
            AE_S32X2RNG_XP(x6, py0, stride * sizeof(ae_int32x2));
            AE_S32X2RNG_XP(x7, py0, (-7 * stride + 1)* sizeof(ae_int32x2));
        }
    }
    else
    {
        ae_int32x2 x00, x10, x20, x30;
        ae_int32x2 x40, x50, x60, x70;
        ae_int32x2 x01, x11, x21, x31;
        ae_int32x2 x41, x51, x61, x71;

        /* Even loop count,  unrolled twice  */
        __Pragma("loop_count min=1");
        for (j = 0; j < (stride >> 1); j++)
        {
            AE_L32X2X2_XP(x70, x71, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x60, x61, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x50, x51, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x40, x41, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), sizeof(ae_int32x2)*(7 * stride + 2));
            /**************** DFT8 ***************/
            DFT4X1RNG(x00, x20, x40, x60);
            DFT4X1RNG(x10, x30, x50, x70);

            x30 = AE_MULFC32RAS(x30, tw1);
            x50 = AE_MULFC32RAS(x50, tw2);
            x70 = AE_MULFC32RAS(x70, tw3);

            {
                ae_int32x2 s0, s1, s2, s3;
                ae_int32x2 d0, d1, d2, d3;

                x00 = AE_SRAI32(x00, 1);
                x10 = AE_SRAI32(x10, 1);
                x20 = AE_SRAI32(x20, 1);
                x40 = AE_SRAI32(x40, 1);
                x60 = AE_SRAI32(x60, 1);

                AE_ADDANDSUB32S(s0, d0, x00, x10);
                AE_ADDANDSUB32S(s1, d1, x20, x30);
                AE_ADDANDSUB32S(s2, d2, x40, x50);
                AE_ADDANDSUB32S(s3, d3, x60, x70);

                x00 = s0;        x40 = d0;
                x10 = s1;        x50 = d1;
                x20 = s2;        x60 = d2;
                x30 = s3;        x70 = d3;
            }
            /**************** DFT8 ***************/
            DFT4X1RNG(x01, x21, x41, x61);
            DFT4X1RNG(x11, x31, x51, x71);

            x31 = AE_MULFC32RAS(x31, tw1);
            x51 = AE_MULFC32RAS(x51, tw2);
            x71 = AE_MULFC32RAS(x71, tw3);

            {
                ae_int32x2 s0, s1, s2, s3;
                ae_int32x2 d0, d1, d2, d3;

                x01 = AE_SRAI32(x01, 1);
                x11 = AE_SRAI32(x11, 1);
                x21 = AE_SRAI32(x21, 1);
                x41 = AE_SRAI32(x41, 1);
                x61 = AE_SRAI32(x61, 1);

                AE_ADDANDSUB32S(s0, d0, x01, x11);
                AE_ADDANDSUB32S(s1, d1, x21, x31);
                AE_ADDANDSUB32S(s2, d2, x41, x51);
                AE_ADDANDSUB32S(s3, d3, x61, x71);

                x01 = s0;        x41 = d0;
                x11 = s1;        x51 = d1;
                x21 = s2;        x61 = d2;
                x31 = s3;        x71 = d3;
            }

            AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x40, x41, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x50, x51, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x60, x61, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x70, x71, castxcc(ae_int32x4, py0), (-7 * stride + 2)* sizeof(ae_int32x2));
        }  /* for (j = 0; j < (stride >> 1); j++) */
    }
    return shift+1;
} /* fft_stageS2_DFT8_last_32x32() */

/*
*  32x32 FFT last stage Radix 8, scalingOption=2
*/
int ifft_stageS2_DFT8_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    int shift;
    int j;
    ae_int32x2 * restrict px0;
    ae_int32x2 * restrict py0;
    const ae_int32x2 * restrict ptwd;

    const int stride = (N >> 3);
    const int min_shift = 3;
    ae_int32x2 tw1, tw2, tw3;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(*bexp >= 0);

    shift = XT_MAX(0, min_shift - *bexp);
    WUR_AE_SAR(shift);
    px0 = (ae_int32x2 *)x + 7 * stride;
    py0 = (ae_int32x2 *)y;
    ptwd = (const ae_int32x2 *)__fft8_tw1;
    tw1 = AE_L32X2_I(ptwd, 0 * sizeof(ae_int32x2));
    tw2 = AE_L32X2_I(ptwd, 1 * sizeof(ae_int32x2));
    tw3 = AE_L32X2_I(ptwd, 2 * sizeof(ae_int32x2));

    if ((stride & 1) > 0)
    {
        /* Odd loop count, used in the some mixed radix FFT */
        __Pragma("loop_count min=1");
        for (j = 0; j < stride; j++)
        {
            ae_int32x2 x0, x1, x2, x3;
            ae_int32x2 x4, x5, x6, x7;

            AE_L32X2_XP(x7, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x6, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x5, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x4, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x3, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x2, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x1, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x0, px0, sizeof(ae_int32x2)*(7 * stride + 1));

            DFT4X1RNG(x0, x2, x4, x6);
            DFT4X1RNG(x1, x3, x5, x7);

            _cmult32x32(&x3, &x3, &tw1);
            _cmult32x32(&x5, &x5, &tw2);
            _cmult32x32(&x7, &x7, &tw3);

            {
                ae_int32x2 s0, s1, s2, s3;
                ae_int32x2 d0, d1, d2, d3;

                x0 = AE_SRAI32(x0, 1);
                x1 = AE_SRAI32(x1, 1);
                x2 = AE_SRAI32(x2, 1);
                x4 = AE_SRAI32(x4, 1);
                x6 = AE_SRAI32(x6, 1);

                AE_ADDANDSUB32S(s0, d0, x0, x1);
                AE_ADDANDSUB32S(s1, d1, x2, x3);
                AE_ADDANDSUB32S(s2, d2, x4, x5);
                AE_ADDANDSUB32S(s3, d3, x6, x7);

                x0 = s0;        x4 = d0;
                x1 = s1;        x5 = d1;
                x2 = s2;        x6 = d2;
                x3 = s3;        x7 = d3;
            }

            AE_S64_XP(AE_MOVINT64_FROMINT32X2(x0), castxcc(ae_int64, py0), stride * sizeof(ae_int32x2));
            AE_S64_XP(AE_MOVINT64_FROMINT32X2(x1), castxcc(ae_int64, py0), stride * sizeof(ae_int32x2));
            AE_S64_XP(AE_MOVINT64_FROMINT32X2(x2), castxcc(ae_int64, py0), stride * sizeof(ae_int32x2));
            AE_S64_XP(AE_MOVINT64_FROMINT32X2(x3), castxcc(ae_int64, py0), stride * sizeof(ae_int32x2));
            AE_S64_XP(AE_MOVINT64_FROMINT32X2(x4), castxcc(ae_int64, py0), stride * sizeof(ae_int32x2));
            AE_S64_XP(AE_MOVINT64_FROMINT32X2(x5), castxcc(ae_int64, py0), stride * sizeof(ae_int32x2));
            AE_S64_XP(AE_MOVINT64_FROMINT32X2(x6), castxcc(ae_int64, py0), stride * sizeof(ae_int32x2));
            AE_S64_XP(AE_MOVINT64_FROMINT32X2(x7), castxcc(ae_int64, py0), (-7 * stride + 1)* sizeof(ae_int32x2));
        }
    }
    else
    {
        ae_int32x2 x00, x10, x20, x30;
        ae_int32x2 x40, x50, x60, x70;
        ae_int32x2 x01, x11, x21, x31;
        ae_int32x2 x41, x51, x61, x71;

        /* Even loop count,  unrolled twice  */
        __Pragma("loop_count min=1");
        for (j = 0; j < (stride >> 1); j++)
        {
            AE_L32X2X2_XP(x70, x71, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x60, x61, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x50, x51, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x40, x41, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), sizeof(ae_int32x2)*(7 * stride + 2));
            /**************** DFT8 ***************/
            DFT4X1RNG(x00, x20, x40, x60);
            DFT4X1RNG(x10, x30, x50, x70);

            x30 = AE_MULFC32RAS(x30, tw1);
            x50 = AE_MULFC32RAS(x50, tw2);
            x70 = AE_MULFC32RAS(x70, tw3);

            {
                ae_int32x2 s0, s1, s2, s3;
                ae_int32x2 d0, d1, d2, d3;

                x00 = AE_SRAI32(x00, 1);
                x10 = AE_SRAI32(x10, 1);
                x20 = AE_SRAI32(x20, 1);
                x40 = AE_SRAI32(x40, 1);
                x60 = AE_SRAI32(x60, 1);

                AE_ADDANDSUB32S(s0, d0, x00, x10);
                AE_ADDANDSUB32S(s1, d1, x20, x30);
                AE_ADDANDSUB32S(s2, d2, x40, x50);
                AE_ADDANDSUB32S(s3, d3, x60, x70);

                x00 = s0;        x40 = d0;
                x10 = s1;        x50 = d1;
                x20 = s2;        x60 = d2;
                x30 = s3;        x70 = d3;
            }
            /**************** DFT8 ***************/
            DFT4X1RNG(x01, x21, x41, x61);
            DFT4X1RNG(x11, x31, x51, x71);

            x31 = AE_MULFC32RAS(x31, tw1);
            x51 = AE_MULFC32RAS(x51, tw2);
            x71 = AE_MULFC32RAS(x71, tw3);

            {
                ae_int32x2 s0, s1, s2, s3;
                ae_int32x2 d0, d1, d2, d3;

                x01 = AE_SRAI32(x01, 1);
                x11 = AE_SRAI32(x11, 1);
                x21 = AE_SRAI32(x21, 1);
                x41 = AE_SRAI32(x41, 1);
                x61 = AE_SRAI32(x61, 1);

                AE_ADDANDSUB32S(s0, d0, x01, x11);
                AE_ADDANDSUB32S(s1, d1, x21, x31);
                AE_ADDANDSUB32S(s2, d2, x41, x51);
                AE_ADDANDSUB32S(s3, d3, x61, x71);

                x01 = s0;        x41 = d0;
                x11 = s1;        x51 = d1;
                x21 = s2;        x61 = d2;
                x31 = s3;        x71 = d3;
            }

            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x00), AE_MOVINT64_FROMINT32X2(x01), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x10), AE_MOVINT64_FROMINT32X2(x11), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x20), AE_MOVINT64_FROMINT32X2(x21), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x30), AE_MOVINT64_FROMINT32X2(x31), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x40), AE_MOVINT64_FROMINT32X2(x41), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x50), AE_MOVINT64_FROMINT32X2(x51), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x60), AE_MOVINT64_FROMINT32X2(x61), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x70), AE_MOVINT64_FROMINT32X2(x71), castxcc(ae_int64x2, py0), (-7 * stride + 2)* sizeof(ae_int32x2));
        }  /* for (j = 0; j < (stride >> 1); j++) */
    }
    return shift + 1;
} /* ifft_stageS2_DFT8_last_32x32() */

/*
 *  32x32 FFT/IFFT intermediate stage Radix 4 unrolled 2 times, scalingOption=2
 */
int fft_stageS2_DFT4x2_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  const ae_int32x2 * restrict px0; 
        ae_int32x2 * restrict py0;
  const ae_int32x2 * restrict ptwd;
  int shift;
  const int min_shift = 3;
  const int R = 4; /* stage radix */
  const int stride = (N >> 2);

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  shift = (min_shift-*bexp);
  WUR_AE_SAR(shift);
  {
    int i, j, _v;
    _v = *v;
    NASSERT(_v>=2 && (_v&1)==0);
    
    if (stride / _v > (_v >> 1))
    {
        __Pragma("loop_count min=1");
        for (j = 0; j < (_v>>1); j++)
        {
            ae_int32x2 x00, x10, x20, x30;
            ae_int32x2 x01, x11, x21, x31;
            ae_int32x2 tw1, tw2, tw3;
            ptwd = (const ae_int32x2 *)tw;
            px0 = (ae_int32x2 *)x + j*2;
            py0 = (ae_int32x2 *)y + j*2;

          __Pragma("loop_count min=2");
          for (i = 0; i < (stride/_v); i++)
          {
            AE_L32X2_IP(tw1, ptwd, sizeof(ae_int32x2));
            AE_L32X2_IP(tw2, ptwd, sizeof(ae_int32x2));
            AE_L32X2_XP(tw3, ptwd, (3 * (tw_step-1)+1)*sizeof(ae_int32x2));
           
            AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (_v - 3 * stride) * sizeof(ae_int32x2));

            DFT4X1RNG(x00, x10, x20, x30);
            x10 = AE_MULFC32RAS(x10, tw1);
            x20 = AE_MULFC32RAS(x20, tw2);
            x30 = AE_MULFC32RAS(x30, tw3);

            DFT4X1RNG(x01, x11, x21, x31);
            x11 = AE_MULFC32RAS(x11, tw1);
            x21 = AE_MULFC32RAS(x21, tw2);
            x31 = AE_MULFC32RAS(x31, tw3);

            AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
          }
        }
    }  
    else /*  if (stride / _v > (_v >> 1)) */
    {
        ptwd = (const ae_int32x2 *)tw;
        __Pragma("loop_count min=1");
        for (i = 0; i < (stride / _v); i++)
        {
            ae_int32x2 tw1, tw2, tw3;

            px0 = (ae_int32x2 *)x + _v * i;
            py0 = (ae_int32x2 *)y + 4 * _v * i;

            AE_L32X2_IP(tw1, ptwd, sizeof(ae_int32x2));
            AE_L32X2_IP(tw2, ptwd, sizeof(ae_int32x2));
            AE_L32X2_XP(tw3, ptwd, (3 * (tw_step - 1) + 1)*sizeof(ae_int32x2));
            __Pragma("loop_count min=1");
            for (j = 0; j < (_v >> 1); j++)
            {
                /* 11 cycles per pipeline stage in steady state with unroll=2 */
                ae_int32x2 x00, x10, x20, x30;
                ae_int32x2 x01, x11, x21, x31;

                AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
                AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
                AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
                AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (2 - 3 * stride) * sizeof(ae_int32x2));

                DFT4X1RNG(x00, x10, x20, x30);
                x10 = AE_MULFC32RAS(x10, tw1);
                x20 = AE_MULFC32RAS(x20, tw2);
                x30 = AE_MULFC32RAS(x30, tw3);

                DFT4X1RNG(x01, x11, x21, x31);
                x11 = AE_MULFC32RAS(x11, tw1);
                x21 = AE_MULFC32RAS(x21, tw2);
                x31 = AE_MULFC32RAS(x31, tw3);

                AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), (2 - 3 * _v) * sizeof(ae_int32x2));
            }
        }

    }   /* if (stride / _v > (_v >> 1)) else */
    AE_CALCRNG3();
    *bexp = 3 - RUR_AE_SAR();
    *v *= R;
    return shift;
  }
} /* fft_stageS2_DFT4x2_32x32 */

/*
 *  32x32 FFT/IFFT intermediate stage Radix 4 unrolled 4 times, scalingOption=2
 *  Restriction: N must be a power of 2
 */
int fft_stageS2_DFT4x4_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  const ae_int32x2 * restrict px0;
  const ae_int32x2 * restrict px1;
        ae_int32x2 * restrict py0;
        ae_int32x2 * restrict py1;
  const ae_int32x2 * restrict ptwd;
  int shift;
  const int R = 4; /* stage radix */
  const int stride = (N >> 2);
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  NASSERT(*bexp >= 0);
  shift = 3 - *bexp;
  WUR_AE_SAR( shift );
  {
    int i, j, _v;
    _v = *v;
    NASSERT((_v&3)==0);
    __Pragma( "loop_count min=1" );
    for (j = 0; j < (_v>>2); j++)
    {
      ae_int32x2 tw1, tw2, tw3;
      ptwd = (const ae_int32x2 *)tw;
      px0 = (const ae_int32x2 *)x+j*4;
      py0 = (ae_int32x2 *)y+j*4;
      px1 = px0 + 2;
      py1 = py0 + 2;

      __Pragma( "loop_count min=1" );
      for (i = 0; i < (stride/_v); i++)
      {
          ae_int32x2 x00, x10, x20, x30;
          ae_int32x2 x01, x11, x21, x31;
        AE_L32X2_IP(tw1, ptwd, sizeof(ae_int32x2));
        AE_L32X2_IP(tw2, ptwd, sizeof(ae_int32x2));
        AE_L32X2_XP(tw3, ptwd, (3*tw_step-2)*sizeof(ae_int32x2));
           
        AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (_v - 3 * stride) * sizeof(ae_int32x2));

        DFT4X1RNG(x00, x10, x20, x30);
        _cmult32x32(&x10, &x10, &tw1);
        _cmult32x32(&x20, &x20, &tw2);
        _cmult32x32(&x30, &x30, &tw3);

        DFT4X1RNG(x01, x11, x21, x31);
        _cmult32x32(&x11, &x11, &tw1);
        _cmult32x32(&x21, &x21, &tw2);
        _cmult32x32(&x31, &x31, &tw3);

        AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));

        AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px1), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px1), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px1), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px1), (_v - 3 * stride) * sizeof(ae_int32x2));

        DFT4X1RNG(x00, x10, x20, x30);
        _cmult32x32(&x10, &x10, &tw1);
        _cmult32x32(&x20, &x20, &tw2);
        _cmult32x32(&x30, &x30, &tw3);

        DFT4X1RNG(x01, x11, x21, x31);
        _cmult32x32(&x11, &x11, &tw1);
        _cmult32x32(&x21, &x21, &tw2);
        _cmult32x32(&x31, &x31, &tw3);

        AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py1), _v * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py1), _v * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py1), _v * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py1), _v * sizeof(ae_int32x2));
      }
    }
    AE_CALCRNG3();
    *bexp = 3 - RUR_AE_SAR();
    *v *= R;
    return shift;
  }
} /* fft_stageS2_DFT4x4_32x32 */

/*
 *  32x32 FFT/IFFT penultimate stage Radix 4 unrolled 2(4) times, scalingOption=2
 *  Restriction: N must be a power of 2
 */
int fft_stageS2_DFT4x4_penultimate_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  int shift;
  const ae_int32x2 * restrict px0;
        ae_int32x2 * restrict py0;
  const ae_int32x2 * restrict ptwd;
  int j, _v;
  ae_int32x2 tw11, tw12, tw13;
  ae_int32x2 tw21, tw22, tw23;
  ae_int32x2 tw31, tw32, tw33;
  const int R = 4; // stage radix
  const int stride = (N >> 2);

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT(0==(N&(N-1)));

  shift = 3 - *bexp;
  WUR_AE_SAR( shift );
  _v = *v;
  NASSERT((stride/_v == 4) || (stride/_v == 2));

  px0 = (ae_int32x2 *)x;
  py0 = (ae_int32x2 *)y;

  if (stride==(_v<<2))
  {
    ptwd = (const ae_int32x2 *)tw + 3*tw_step;
    tw11 = AE_L32X2_I(ptwd, 0*sizeof(ae_int32x2));
    tw12 = AE_L32X2_I(ptwd, 1*sizeof(ae_int32x2));
    tw13 = AE_L32X2_I(ptwd, 2*sizeof(ae_int32x2));
    ptwd = ptwd + 3*tw_step;
    tw21 = AE_L32X2_I(ptwd, 0*sizeof(ae_int32x2));
    tw22 = AE_L32X2_I(ptwd, 1*sizeof(ae_int32x2));
    tw23 = AE_L32X2_I(ptwd, 2*sizeof(ae_int32x2));
    ptwd = ptwd + 3*tw_step;
    tw31 = AE_L32X2_I(ptwd, 0*sizeof(ae_int32x2));
    tw32 = AE_L32X2_I(ptwd, 1*sizeof(ae_int32x2));
    tw33 = AE_L32X2_I(ptwd, 2*sizeof(ae_int32x2));

    __Pragma("loop_count min=1");
    for (j = 0; j < (_v >> 1); j++)
    {
        /* 20 cycles per pipeline stage in steady state with unroll=1
        2 pipeline stages */
        px0 = (ae_int32x2 *)x + j * 2;
        py0 = (ae_int32x2 *)y + j * 2;

        {
            ae_int32x2 x00, x10, x20, x30;
            ae_int32x2 x01, x11, x21, x31;


            AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (_v - 3 * stride) * sizeof(ae_int32x2));

            DFT4X1RNG(x00, x10, x20, x30);
            DFT4X1RNG(x01, x11, x21, x31);

            AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
        }
        {
            ae_int32x2 x00, x10, x20, x30;
            ae_int32x2 x01, x11, x21, x31;

            AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (_v - 3 * stride) * sizeof(ae_int32x2));

            DFT4X1RNG(x00, x10, x20, x30);
            DFT4X1RNG(x01, x11, x21, x31);

            x10 = AE_MULFC32RAS(x10, tw11);
            x20 = AE_MULFC32RAS(x20, tw21);
            x30 = AE_MULFC32RAS(x30, tw31);
            x11 = AE_MULFC32RAS(x11, tw11);
            x21 = AE_MULFC32RAS(x21, tw21);
            x31 = AE_MULFC32RAS(x31, tw31);

            AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
        }
        {
            ae_int32x2 x00, x10, x20, x30;
            ae_int32x2 x01, x11, x21, x31;

            AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (_v - 3 * stride) * sizeof(ae_int32x2));

            DFT4X1RNG(x00, x10, x20, x30);
            DFT4X1RNG(x01, x11, x21, x31);

            x10 = AE_MULFC32RAS(x10, tw12);
            x20 = AE_MULFC32RAS(x20, tw22);
            x30 = AE_MULFC32RAS(x30, tw32);
            x11 = AE_MULFC32RAS(x11, tw12);
            x21 = AE_MULFC32RAS(x21, tw22);
            x31 = AE_MULFC32RAS(x31, tw32);

            AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
        }
        {
            ae_int32x2 x00, x10, x20, x30;
            ae_int32x2 x01, x11, x21, x31;

            AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (_v - 3 * stride) * sizeof(ae_int32x2));

            DFT4X1RNG(x00, x10, x20, x30);
            DFT4X1RNG(x01, x11, x21, x31);

            x11 = AE_MULFC32RAS(x11, tw13);
            x21 = AE_MULFC32RAS(x21, tw23);
            x31 = AE_MULFC32RAS(x31, tw33);
            x10 = AE_MULFC32RAS(x10, tw13);
            x20 = AE_MULFC32RAS(x20, tw23);
            x30 = AE_MULFC32RAS(x30, tw33);

            AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
        }
    }
  }
  else /* if (stride==(_v<<2)) */
  {
    ae_int32x2 tw1, tw2, tw3;
    ptwd = (const ae_int32x2 *)tw + 3*tw_step;
    tw1 = AE_L32X2_I(ptwd, 0*sizeof(ae_int32x2));
    tw2 = AE_L32X2_I(ptwd, 1*sizeof(ae_int32x2));
    tw3 = AE_L32X2_I(ptwd, 2*sizeof(ae_int32x2));

    __Pragma("loop_count min=1");
    for (j = 0; j < (_v>>1); j++)
    {
        px0 = (ae_int32x2 *)x + j * 2;
        py0 = (ae_int32x2 *)y + j * 2;

        {
            ae_int32x2 x00, x10, x20, x30;
            ae_int32x2 x01, x11, x21, x31;

            AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (_v - 3 * stride) * sizeof(ae_int32x2));

            DFT4X1RNG(x00, x10, x20, x30);
            DFT4X1RNG(x01, x11, x21, x31);

            AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
        }
        {
            ae_int32x2 x00, x10, x20, x30;
            ae_int32x2 x01, x11, x21, x31;

            AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (_v - 3 * stride) * sizeof(ae_int32x2));

            DFT4X1RNG(x00, x10, x20, x30);
            x10 = AE_MULFC32RAS(x10, tw1);
            x20 = AE_MULFC32RAS(x20, tw2);
            x30 = AE_MULFC32RAS(x30, tw3);

            DFT4X1RNG(x01, x11, x21, x31);
            x11 = AE_MULFC32RAS(x11, tw1);
            x21 = AE_MULFC32RAS(x21, tw2);
            x31 = AE_MULFC32RAS(x31, tw3);

            AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
        }
    }
  } /*if (stride==(_v<<2)) else ... */

  AE_CALCRNG3();
  *bexp = 3 - RUR_AE_SAR();
  *v *= R;
  return shift;
} /* fft_stageS2_DFT4x4_penultimate_32x32 */

/*
 *  32x32 IFFT first stage Radix 4, scalingOption=2
 */
int ifft_stageS2_DFT4_first_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    int i;
    int shift;
    ae_int32x2 * restrict px0;
    ae_int32x2 * restrict px1;
    ae_int32x2 * restrict px2;
    ae_int32x2 * restrict px3;
    ae_int32x2 * restrict py0;
    const ae_int32x2 * restrict ptwd;
    const int R = 4; // stage radix
    const int stride = (N >> 2);
    shift = 3 - *bexp;
    ae_int64 tmp00, tmp10, tmp20, tmp30;
    ae_int64 tmp01, tmp11, tmp21, tmp31;

    NASSERT(shift>-32 && shift<4);
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    AE_CALCRNG3();
    // scl = 1 << shiftl;
    WUR_AE_SAR(3);
    px0 = (ae_int32x2 *)x;
    px1 = px0 + stride;
    px2 = px1 + stride;
    px3 = px2 + stride;
    py0 = (ae_int32x2 *)y;
    ptwd = (const ae_int32x2 *)tw;
    if ((stride & 1) != 0 || tw_step != 1)
    {
        __Pragma("loop_count min=3");
        for (i = 0; i < stride; i++)
        {
            ae_int32x2 x0, x1, x2, x3;
            ae_int32x2 tw1, tw2, tw3;

            AE_L32X2_IP(tw1, ptwd, sizeof(ae_int32x2));
            AE_L32X2_IP(tw2, ptwd, sizeof(ae_int32x2));
            AE_L32X2_XP(tw3, ptwd, (3 * tw_step - 2) * sizeof(ae_int32x2));

            AE_L64_IP(tmp00, castxcc(ae_int64, px0), sizeof(ae_int32x2));
            AE_L64_IP(tmp10, castxcc(ae_int64, px1), sizeof(ae_int32x2));
            AE_L64_IP(tmp20, castxcc(ae_int64, px2), sizeof(ae_int32x2));
            AE_L64_IP(tmp30, castxcc(ae_int64, px3), sizeof(ae_int32x2));

            x0 = AE_MOVINT32X2_FROMINT64(tmp00);
            x1 = AE_MOVINT32X2_FROMINT64(tmp10);
            x2 = AE_MOVINT32X2_FROMINT64(tmp20);
            x3 = AE_MOVINT32X2_FROMINT64(tmp30);

            x0 = AE_SLAA32S(x0, *bexp);
            x1 = AE_SLAA32S(x1, *bexp);
            x2 = AE_SLAA32S(x2, *bexp);
            x3 = AE_SLAA32S(x3, *bexp);

            DFT4X1RNG(x0, x1, x2, x3);
            x1 = AE_MULFC32RAS(x1, tw1);
            x2 = AE_MULFC32RAS(x2, tw2);
            x3 = AE_MULFC32RAS(x3, tw3);

            AE_S32X2RNG_IP(x0, py0, sizeof(ae_int32x2));
            AE_S32X2RNG_IP(x1, py0, sizeof(ae_int32x2));
            AE_S32X2RNG_IP(x2, py0, sizeof(ae_int32x2));
            AE_S32X2RNG_IP(x3, py0, sizeof(ae_int32x2));
        }
    }
    else
    {
        __Pragma("loop_count min=1");
        for (i = 0; i < (stride >> 1); i++)
        {
            /* 6 cycles per pipeline stage in steady state with unroll=1
            2 pipeline stages */
            ae_int32x2 x00, x10, x20, x30;
            ae_int32x2 x01, x11, x21, x31;
            ae_int32x2 tw10, tw20, tw30;
            ae_int32x2 tw11, tw21, tw31;


            AE_L32X2X2_IP(tw10, tw20, castxcc(ae_int32x4, ptwd), sizeof(ae_int32x4));
            AE_L32X2X2_IP(tw30, tw11, castxcc(ae_int32x4, ptwd), sizeof(ae_int32x4));
            AE_L32X2X2_IP(tw21, tw31, castxcc(ae_int32x4, ptwd), sizeof(ae_int32x4));

            AE_L64X2_IP(tmp00, tmp01, castxcc(ae_int64x2, px0), 2 * sizeof(ae_int32x2));
            AE_L64X2_IP(tmp10, tmp11, castxcc(ae_int64x2, px1), 2 * sizeof(ae_int32x2));
            AE_L64X2_IP(tmp20, tmp21, castxcc(ae_int64x2, px2), 2 * sizeof(ae_int32x2));
            AE_L64X2_IP(tmp30, tmp31, castxcc(ae_int64x2, px3), 2 * sizeof(ae_int32x2));

            x00 = AE_MOVINT32X2_FROMINT64(tmp00);
            x10 = AE_MOVINT32X2_FROMINT64(tmp10);
            x20 = AE_MOVINT32X2_FROMINT64(tmp20);
            x30 = AE_MOVINT32X2_FROMINT64(tmp30);
            x01 = AE_MOVINT32X2_FROMINT64(tmp01);
            x11 = AE_MOVINT32X2_FROMINT64(tmp11);
            x21 = AE_MOVINT32X2_FROMINT64(tmp21);
            x31 = AE_MOVINT32X2_FROMINT64(tmp31);

            x00 = AE_SLAA32S(x00, *bexp);
            x10 = AE_SLAA32S(x10, *bexp);
            x20 = AE_SLAA32S(x20, *bexp);
            x30 = AE_SLAA32S(x30, *bexp);
            x01 = AE_SLAA32S(x01, *bexp);
            x11 = AE_SLAA32S(x11, *bexp);
            x21 = AE_SLAA32S(x21, *bexp);
            x31 = AE_SLAA32S(x31, *bexp);

            DFT4X1RNG(x00, x10, x20, x30);
            DFT4X1RNG(x01, x11, x21, x31);

            x10 = AE_MULFC32RAS(x10, tw10);
            x20 = AE_MULFC32RAS(x20, tw20);
            x30 = AE_MULFC32RAS(x30, tw30);
            x11 = AE_MULFC32RAS(x11, tw11);
            x21 = AE_MULFC32RAS(x21, tw21);
            x31 = AE_MULFC32RAS(x31, tw31);

            AE_S32X2X2RNG_IP(x00, x10, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
            AE_S32X2X2RNG_IP(x20, x30, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
            AE_S32X2X2RNG_IP(x01, x11, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
            AE_S32X2X2RNG_IP(x21, x31, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
        } /* for (i = 0; i < (stride>>1); i++) */
    }
    AE_CALCRNG3();
    *bexp = 3 - RUR_AE_SAR();
    *v *= R;
    return shift;

} /* ifft_stageS2_DFT4_first_32x32() */

int ifft_stageS2_DFT4_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp){ return 0; }
