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
    with static data scaling: 32-bit data, 32-bit twiddle factors
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
/* radix-4 butterfly with normalization */
#define DFT4X1RNG(x0, x1, x2, x3)\
{   \
    ae_int32x2 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG32(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32(s1, d1, x1, x3); \
    AE_ADDANDSUB32S(x0, x2, s0, s1);   \
    AE_ADDANDSUB32JS(x3, x1, d0, d1);   \
}

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

int fft_stageS3_DFT4_first_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  int i;
  const ae_int32x2 * restrict px0;
  const ae_int32x2 * restrict px1;
  const ae_int32x2 * restrict px2;
  const ae_int32x2 * restrict px3;
        ae_int32x2 * restrict py0;
  const ae_int32x2 * restrict ptwd;
  const int R = 4; // stage radix
  const int stride = (N >> 2);

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT_ALIGN16(tw);
  //ASSERT(tw_step == 1);

  px0 = (const ae_int32x2 *)x;
  px1 = (const ae_int32x2 *)px0 + stride;
  px2 = (const ae_int32x2 *)px1 + stride;
  px3 = (const ae_int32x2 *)px2 + stride;
  py0 = (      ae_int32x2 *)y;

  ptwd = (const ae_int32x2 *)tw;
  WUR_AE_SAR(3);
 
  if ( (stride & 1) != 0 || tw_step != 1)
  {
      __Pragma("loop_count min=3");
      for (i = 0; i < stride; i++)
      {
          /* 11 cycles per pipeline stage in steady state with unroll=2 */
          ae_int32x2 x00, x10, x20, x30;
          ae_int32x2 tw10, tw20, tw30;

          AE_L32X2_IP(tw10, ptwd, sizeof(ae_int32x2));
          AE_L32X2_IP(tw20, ptwd, sizeof(ae_int32x2));
          AE_L32X2_XP(tw30, ptwd, (3 * tw_step - 2) * sizeof(ae_int32x2));

          AE_L32X2_IP(x00, px0, sizeof(ae_int32x2));
          AE_L32X2_IP(x10, px1, sizeof(ae_int32x2));
          AE_L32X2_IP(x20, px2, sizeof(ae_int32x2));
          AE_L32X2_IP(x30, px3, sizeof(ae_int32x2));

          DFT4X1RNG(x00, x10, x20, x30);
          x10 = AE_MULFC32RAS(x10, tw10);
          x20 = AE_MULFC32RAS(x20, tw20);
          x30 = AE_MULFC32RAS(x30, tw30);

          AE_S32X2_IP(x00, py0, sizeof(ae_int32x2));
          AE_S32X2_IP(x10, py0, sizeof(ae_int32x2));
          AE_S32X2_IP(x20, py0, sizeof(ae_int32x2));
          AE_S32X2_IP(x30, py0, sizeof(ae_int32x2));
        }
  }
  else
  {
      __Pragma("loop_count min=1");
       for (i = 0; i < (stride>>1); i++)
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

            DFT4X1RNG(x00, x10, x20, x30);
            DFT4X1RNG(x01, x11, x21, x31);

            x10 = AE_MULFC32RAS(x10, tw10);
            x20 = AE_MULFC32RAS(x20, tw20);
            x30 = AE_MULFC32RAS(x30, tw30);
            x11 = AE_MULFC32RAS(x11, tw11);
            x21 = AE_MULFC32RAS(x21, tw21);
            x31 = AE_MULFC32RAS(x31, tw31);

            AE_S32X2X2_IP(x00, x10, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
            AE_S32X2X2_IP(x20, x30, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
            AE_S32X2X2_IP(x01, x11, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
            AE_S32X2X2_IP(x21, x31, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
        } /* for (i = 0; i < (stride>>1); i++) */
    } /*if(stride&1) else... */

  *v = v[0] * R;
  return 3;
} /* fft_stageS3_DFT4_first_32x32() */

int fft_stageS3_DFT4_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  const ae_int32x2 * restrict px0;
  const ae_int32x2 * restrict px1;
  const ae_int32x2 * restrict px2;
  const ae_int32x2 * restrict px3;
        ae_int32x2 * restrict py0;
        ae_int32x2 * restrict py1;
        ae_int32x2 * restrict py2;
        ae_int32x2 * restrict py3;
  const int stride = (N >> 2);
  int j, shift;
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  
  shift = 2;
  WUR_AE_SAR(shift);

  px0 = (const ae_int32x2 *)x;
  px1 = (const ae_int32x2 *)px0 + stride;
  px2 = (const ae_int32x2 *)px1 + stride;
  px3 = (const ae_int32x2 *)px2 + stride;
  py0 = (      ae_int32x2 *)y;
  py1 = (      ae_int32x2 *)py0 + stride;
  py2 = (      ae_int32x2 *)py1 + stride;
  py3 = (      ae_int32x2 *)py2 + stride;

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

      AE_S32X2_IP(x0, py0, sizeof(ae_int32x2));
      AE_S32X2_IP(x1, py1, sizeof(ae_int32x2));
      AE_S32X2_IP(x2, py2, sizeof(ae_int32x2));
      AE_S32X2_IP(x3, py3, sizeof(ae_int32x2));
    }
  }
  else
  {
    ae_int32x2 x0, x1, x2, x3;
    ae_int32x2 x4, x5, x6, x7;

      __Pragma("loop_count min=1");
    for (j = 0; j < (stride>>1); j++)
    {
        AE_L32X2X2_XP(x0, x4, castxcc(ae_int32x4, px0), stride*sizeof(complex_fract32));
        AE_L32X2X2_XP(x1, x5, castxcc(ae_int32x4, px0), stride*sizeof(complex_fract32));
        AE_L32X2X2_XP(x2, x6, castxcc(ae_int32x4, px0), stride*sizeof(complex_fract32));
        AE_L32X2X2_XP(x3, x7, castxcc(ae_int32x4, px0), (2 - 3 * stride)*sizeof(complex_fract32));

        DFT4X1RNG(x0, x1, x2, x3);
        DFT4X1RNG(x4, x5, x6, x7);

        AE_S32X2X2_XP(x0, x4, castxcc(ae_int32x4, py0), stride*sizeof(complex_fract32));
        AE_S32X2X2_XP(x1, x5, castxcc(ae_int32x4, py0), stride*sizeof(complex_fract32));
        AE_S32X2X2_XP(x2, x6, castxcc(ae_int32x4, py0), stride*sizeof(complex_fract32));
        AE_S32X2X2_XP(x3, x7, castxcc(ae_int32x4, py0), (2 - 3 * stride)*sizeof(complex_fract32));
    }
  }

  return shift;

} /* fft_stageS3_DFT4_last_32x32() */


int ifft_stageS3_DFT4_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    const ae_int32x2 * restrict px0;
    const ae_int32x2 * restrict px1;
    const ae_int32x2 * restrict px2;
    const ae_int32x2 * restrict px3;
    ae_int32x2 * restrict py0;
    ae_int32x2 * restrict py1;
    ae_int32x2 * restrict py2;
    ae_int32x2 * restrict py3;
    const int stride = (N >> 2);
    int j, shift;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);

    shift = 2;
    WUR_AE_SAR(shift);

    px0 = (const ae_int32x2 *)x;
    px1 = (const ae_int32x2 *)px0 + stride;
    px2 = (const ae_int32x2 *)px1 + stride;
    px3 = (const ae_int32x2 *)px2 + stride;
    py0 = (ae_int32x2 *)y;
    py1 = (ae_int32x2 *)py0 + stride;
    py2 = (ae_int32x2 *)py1 + stride;
    py3 = (ae_int32x2 *)py2 + stride;

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

        __Pragma("loop_count min=1");
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

} /* ifft_stageS3_DFT4_last_32x32() */

int fft_stageS3_DFT4x2_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  const ae_int32x2 * restrict px0; 
        ae_int32x2 * restrict py0;
  const ae_int32x2 * restrict ptwd;
  const int shift = 2;
  const int R = 4; /* stage radix */
  const int stride = (N >> 2);
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  WUR_AE_SAR(shift);
  {
    int i, j, _v;
    _v = v[0];
    NASSERT(_v>=2 && (_v&1)==0);
    if (stride / _v > (_v >> 1))
    {
        __Pragma("loop_count min=1");
        for (j = 0; j < (_v>>1); j++)
        {
          px0 = (ae_int32x2 *)x + j * 2;
          py0 = (ae_int32x2 *)y + j * 2;
          ptwd = (const ae_int32x2 *)tw;
          /* 6 cycles per pipeline stage in steady state with unroll=1 */
          __Pragma( "loop_count min=1" );
          for (i = 0; i <   (stride/ _v); i++)
          {
            ae_int32x2 x00, x10, x20, x30;
            ae_int32x2 x01, x11, x21, x31;
            ae_int32x2 tw1, tw2, tw3;

          //  px0 = (ae_int32x2 *)x + j * 2 + _v * i;
          //  py0 = (ae_int32x2 *)y + j * 2 + 4 * _v * i;

            AE_L32X2_IP(tw1, ptwd, sizeof(ae_int32x2));
            AE_L32X2_IP(tw2, ptwd, sizeof(ae_int32x2));
            AE_L32X2_XP(tw3, ptwd, (3 * (tw_step-1)+1)*sizeof(ae_int32x2));

            AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
            AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (_v-3*stride) * sizeof(ae_int32x2));

            DFT4X1RNG(x00, x10, x20, x30);
            x10 = AE_MULFC32RAS(x10, tw1);
            x20 = AE_MULFC32RAS(x20, tw2); 
            x30 = AE_MULFC32RAS(x30, tw3);

            DFT4X1RNG(x01, x11, x21, x31);
            x11 = AE_MULFC32RAS(x11, tw1);
            x21 = AE_MULFC32RAS(x21, tw2);
            x31 = AE_MULFC32RAS(x31, tw3);

            AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
          }
        }
    }
    else
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

                AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), (2-3*_v) * sizeof(ae_int32x2));
            }
        }
    }
    *v *= R;
    return shift;
  }
} /* fft_stageS3_DFT4x2_32x32 */

int fft_stageS3_DFT4x4_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{ 
  const ae_int32x2 * restrict px0;
        ae_int32x2 * restrict py0;
  const ae_int32x2 * restrict px1;
        ae_int32x2 * restrict py1;
  const ae_int32x2 * restrict ptwd;
  const int shift = 2;
  const int R = 4; /* stage radix */
  const int stride = (N >> 2);
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  WUR_AE_SAR(shift);
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

      __Pragma( "loop_count min=2" );
      for (i = 0; i < (stride/_v); i++)
      {
          /* 11 cycles per pipeline stage in steady state with unroll=1 */
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

        AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
        AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
        AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
        AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));

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

        AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py1), _v * sizeof(ae_int32x2));
        AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py1), _v * sizeof(ae_int32x2));
        AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py1), _v * sizeof(ae_int32x2));
        AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py1), _v * sizeof(ae_int32x2));
      } /* for (i = 0; i < (stride>>log_v); i++) */
    }
    *v *= R;
    return shift;
  }
} /* fft_stageS3_DFT4x4_32x32 */

/*
 *  32x32 FFT/IFFT penultimate stage Radix 4 unrolled 2(4) times, scalingOption=3
 *  Restriction: N must be a power of 2
 */
int fft_stageS3_DFT4x4_penultimate_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
#if 1
    const ae_int32x2 * restrict px0;
    ae_int32x2 * restrict py0;
    const ae_int32x2 * restrict ptwd;
    const int shift = 2;
    const int R = 4; /* stage radix */
    const int stride = (N >> 2);
    int j, _v;
    _v = v[0];

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT((stride / _v == 4) || (stride / _v == 2));

    WUR_AE_SAR(shift);
    if (stride == (_v<<2))
    {
        ae_int32x2 tw11, tw21, tw31;
        ae_int32x2 tw12, tw22, tw32;
        ae_int32x2 tw13, tw23, tw33;

        ptwd = 3*tw_step + (const ae_int32x2 *)tw;

        AE_L32X2_IP(tw11, ptwd, sizeof(ae_int32x2));
        AE_L32X2_IP(tw21, ptwd, sizeof(ae_int32x2));
        AE_L32X2_XP(tw31, ptwd, (3 * (tw_step - 1) + 1)*sizeof(ae_int32x2));
        AE_L32X2_IP(tw12, ptwd, sizeof(ae_int32x2));
        AE_L32X2_IP(tw22, ptwd, sizeof(ae_int32x2));
        AE_L32X2_XP(tw32, ptwd, (3 * (tw_step - 1) + 1)*sizeof(ae_int32x2));
        AE_L32X2_IP(tw13, ptwd, sizeof(ae_int32x2));
        AE_L32X2_IP(tw23, ptwd, sizeof(ae_int32x2));
        AE_L32X2_XP(tw33, ptwd, (3 * (tw_step - 1) + 1)*sizeof(ae_int32x2));

        NASSERT(_v >= 2 && (_v & 1) == 0);
        for (j = 0; j < (_v>>1); j++)
        {
            /* 20 cycles per pipeline stage in steady state with unroll=1 
               2 pipeline stages */
            px0 = (ae_int32x2 *)x + j*2;
            py0 = (ae_int32x2 *)y + j*2;

            {
                ae_int32x2 x00, x10, x20, x30;
                ae_int32x2 x01, x11, x21, x31;


                AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
                AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
                AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
                AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (_v-3*stride) * sizeof(ae_int32x2));

                DFT4X1RNG(x00, x10, x20, x30);
                DFT4X1RNG(x01, x11, x21, x31);

                AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
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

                AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
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

                AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
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

                AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            }
        }
        *v *= R;
        return shift;
    }
    else
    {
        ptwd = 3 * tw_step + (const ae_int32x2 *)tw;
        ae_int32x2 tw1, tw2, tw3;

        AE_L32X2_IP(tw1, ptwd, sizeof(ae_int32x2));
        AE_L32X2_IP(tw2, ptwd, sizeof(ae_int32x2));
        AE_L32X2_XP(tw3, ptwd, (3 * (tw_step - 1) + 1)*sizeof(ae_int32x2));

        for (j = 0; j < (_v>>1); j++)
        {
            
            px0 = (ae_int32x2 *)x + j*2;
            py0 = (ae_int32x2 *)y + j*2;

            {
                ae_int32x2 x00, x10, x20, x30;
                ae_int32x2 x01, x11, x21, x31;

                AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
                AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
                AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
                AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (_v-3*stride) * sizeof(ae_int32x2));

                DFT4X1RNG(x00, x10, x20, x30);
                DFT4X1RNG(x01, x11, x21, x31);

                AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
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

                AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
                AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            }
        }
        *v *= R;
        return shift;

    }
#else
  const ae_int32x2 * restrict px0;
  const ae_int32x2 * restrict px1;
  const ae_int32x2 * restrict px2;
  const ae_int32x2 * restrict px3;
        ae_int32x2 * restrict py0;
        ae_int32x2 * restrict py1;
        ae_int32x2 * restrict py2;
        ae_int32x2 * restrict py3;
  const ae_int32x2 * restrict ptwd;
  int j, _v;
  ae_int32x2 x0, x1, x2, x3;
  ae_int32x2 tw11, tw12, tw13;
  ae_int32x2 tw21, tw22, tw23;
  ae_int32x2 tw31, tw32, tw33;
  const int shift = 2;
  const int R = 4; // stage radix
  const int stride = (N >> 2);

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT(0==(N&(N-1)));

  WUR_AE_SAR( shift>>1 );
  _v = *v;
  NASSERT((stride/_v == 4) || (stride/_v == 2));

  px0 = (ae_int32x2 *)x;
  px1 = px0 + stride;
  px2 = px1 + stride;
  px3 = px2 + stride;
  py0 = (ae_int32x2 *)y;
  py1 = py0 + _v;
  py2 = py1 + _v;
  py3 = py2 + _v;

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

    __Pragma("loop_count min=4, factor=4");
    for (j = 0; j < _v; j++)
    {
      AE_L32X2_XP(x0, px0, _v*sizeof(ae_int32x2));
      AE_L32X2_XP(x1, px1, _v*sizeof(ae_int32x2));
      AE_L32X2_XP(x2, px2, _v*sizeof(ae_int32x2));
      AE_L32X2_XP(x3, px3, _v*sizeof(ae_int32x2));
      DFT4X1RNG2(x0, x1, x2, x3);
      AE_S32X2_XP(x0, py0, 4*_v*sizeof(ae_int32x2));
      AE_S32X2_XP(x1, py1, 4*_v*sizeof(ae_int32x2));
      AE_S32X2_XP(x2, py2, 4*_v*sizeof(ae_int32x2));
      AE_S32X2_XP(x3, py3, 4*_v*sizeof(ae_int32x2));

      AE_L32X2_XP(x0, px0, _v*sizeof(ae_int32x2));
      AE_L32X2_XP(x1, px1, _v*sizeof(ae_int32x2));
      AE_L32X2_XP(x2, px2, _v*sizeof(ae_int32x2));
      AE_L32X2_XP(x3, px3, _v*sizeof(ae_int32x2));
      DFT4X1RNG2(x0, x1, x2, x3);
      _cmult32x32(&x1, &x1, &tw11);
      _cmult32x32(&x2, &x2, &tw12);
      _cmult32x32(&x3, &x3, &tw13);
      AE_S32X2_XP(x0, py0, 4*_v*sizeof(ae_int32x2));
      AE_S32X2_XP(x1, py1, 4*_v*sizeof(ae_int32x2));
      AE_S32X2_XP(x2, py2, 4*_v*sizeof(ae_int32x2));
      AE_S32X2_XP(x3, py3, 4*_v*sizeof(ae_int32x2));

      AE_L32X2_XP(x0, px0, _v*sizeof(ae_int32x2));
      AE_L32X2_XP(x1, px1, _v*sizeof(ae_int32x2));
      AE_L32X2_XP(x2, px2, _v*sizeof(ae_int32x2));
      AE_L32X2_XP(x3, px3, _v*sizeof(ae_int32x2));
      DFT4X1RNG2(x0, x1, x2, x3);
      _cmult32x32(&x1, &x1, &tw21);
      _cmult32x32(&x2, &x2, &tw22);
      _cmult32x32(&x3, &x3, &tw23);
      AE_S32X2_XP(x0, py0, 4*_v*sizeof(ae_int32x2));
      AE_S32X2_XP(x1, py1, 4*_v*sizeof(ae_int32x2));
      AE_S32X2_XP(x2, py2, 4*_v*sizeof(ae_int32x2));
      AE_S32X2_XP(x3, py3, 4*_v*sizeof(ae_int32x2));

      AE_L32X2_XP(x0, px0, (1-3*_v)*sizeof(ae_int32x2));
      AE_L32X2_XP(x1, px1, (1-3*_v)*sizeof(ae_int32x2));
      AE_L32X2_XP(x2, px2, (1-3*_v)*sizeof(ae_int32x2));
      AE_L32X2_XP(x3, px3, (1-3*_v)*sizeof(ae_int32x2));
      DFT4X1RNG2(x0, x1, x2, x3);
      _cmult32x32(&x1, &x1, &tw31);
      _cmult32x32(&x2, &x2, &tw32);
      _cmult32x32(&x3, &x3, &tw33);
      AE_S32X2_XP(x0, py0, (1-3*4*_v)*sizeof(ae_int32x2));
      AE_S32X2_XP(x1, py1, (1-3*4*_v)*sizeof(ae_int32x2));
      AE_S32X2_XP(x2, py2, (1-3*4*_v)*sizeof(ae_int32x2));
      AE_S32X2_XP(x3, py3, (1-3*4*_v)*sizeof(ae_int32x2));

    }
  }
  else
  {
    ptwd = (const ae_int32x2 *)tw + 3*tw_step;
    tw11 = AE_L32X2_I(ptwd, 0*sizeof(ae_int32x2));
    tw12 = AE_L32X2_I(ptwd, 1*sizeof(ae_int32x2));
    tw13 = AE_L32X2_I(ptwd, 2*sizeof(ae_int32x2));

    __Pragma("loop_count min=4, factor=4");
    for (j = 0; j < _v; j++)
    {
      AE_L32X2_XP(x0, px0, _v*sizeof(ae_int32x2));
      AE_L32X2_XP(x1, px1, _v*sizeof(ae_int32x2));
      AE_L32X2_XP(x2, px2, _v*sizeof(ae_int32x2));
      AE_L32X2_XP(x3, px3, _v*sizeof(ae_int32x2));
      DFT4X1RNG2(x0, x1, x2, x3);
      AE_S32X2_XP(x0, py0, 4*_v*sizeof(ae_int32x2));
      AE_S32X2_XP(x1, py1, 4*_v*sizeof(ae_int32x2));
      AE_S32X2_XP(x2, py2, 4*_v*sizeof(ae_int32x2));
      AE_S32X2_XP(x3, py3, 4*_v*sizeof(ae_int32x2));

      AE_L32X2_XP(x0, px0, (1-_v)*sizeof(ae_int32x2));
      AE_L32X2_XP(x1, px1, (1-_v)*sizeof(ae_int32x2));
      AE_L32X2_XP(x2, px2, (1-_v)*sizeof(ae_int32x2));
      AE_L32X2_XP(x3, px3, (1-_v)*sizeof(ae_int32x2));
      DFT4X1RNG2(x0, x1, x2, x3);
      _cmult32x32(&x1, &x1, &tw11);
      _cmult32x32(&x2, &x2, &tw12);
      _cmult32x32(&x3, &x3, &tw13);
      AE_S32X2_XP(x0, py0, (1-4*_v)*sizeof(ae_int32x2));
      AE_S32X2_XP(x1, py1, (1-4*_v)*sizeof(ae_int32x2));
      AE_S32X2_XP(x2, py2, (1-4*_v)*sizeof(ae_int32x2));
      AE_S32X2_XP(x3, py3, (1-4*_v)*sizeof(ae_int32x2));

    }
  }

  *v *= R;
  return shift;
#endif
} /* fft_stageS3_DFT4x4_penultimate_32x32 */

int ifft_stageS3_DFT4_first_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    int i;
    const ae_int32x2 * restrict px0;
    const ae_int32x2 * restrict px1;
    const ae_int32x2 * restrict px2;
    const ae_int32x2 * restrict px3;
    ae_int32x2 * restrict py0;
    const ae_int32x2 * restrict ptwd;
    const int R = 4; // stage radix
    const int stride = (N >> 2);
    ae_int64 tmp00, tmp10, tmp20, tmp30;
    ae_int64 tmp01, tmp11, tmp21, tmp31;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT_ALIGN16(tw);
    //ASSERT(tw_step == 1);

    px0 = (const ae_int32x2 *)x;
    px1 = (const ae_int32x2 *)px0 + stride;
    px2 = (const ae_int32x2 *)px1 + stride;
    px3 = (const ae_int32x2 *)px2 + stride;
    py0 = (ae_int32x2 *)y;

    ptwd = (const ae_int32x2 *)tw;
    WUR_AE_SAR(3);

    if ((stride & 1) != 0 || tw_step != 1)
    {
        __Pragma("loop_count min=3");
        for (i = 0; i < stride; i++)
        {
            /* 11 cycles per pipeline stage in steady state with unroll=2 */
            ae_int32x2 x00, x10, x20, x30;
            ae_int32x2 tw10, tw20, tw30;

            AE_L32X2_IP(tw10, ptwd, sizeof(ae_int32x2));
            AE_L32X2_IP(tw20, ptwd, sizeof(ae_int32x2));
            AE_L32X2_XP(tw30, ptwd, (3 * tw_step - 2) * sizeof(ae_int32x2));

            AE_L64_IP(tmp00, castxcc(ae_int64, px0), sizeof(ae_int32x2));
            AE_L64_IP(tmp10, castxcc(ae_int64, px1), sizeof(ae_int32x2));
            AE_L64_IP(tmp20, castxcc(ae_int64, px2), sizeof(ae_int32x2));
            AE_L64_IP(tmp30, castxcc(ae_int64, px3), sizeof(ae_int32x2));

            x00 = AE_MOVINT32X2_FROMINT64(tmp00);
            x10 = AE_MOVINT32X2_FROMINT64(tmp10);
            x20 = AE_MOVINT32X2_FROMINT64(tmp20);
            x30 = AE_MOVINT32X2_FROMINT64(tmp30);

            DFT4X1RNG(x00, x10, x20, x30);
            x10 = AE_MULFC32RAS(x10, tw10);
            x20 = AE_MULFC32RAS(x20, tw20);
            x30 = AE_MULFC32RAS(x30, tw30);

            AE_S32X2_IP(x00, py0, sizeof(ae_int32x2));
            AE_S32X2_IP(x10, py0, sizeof(ae_int32x2));
            AE_S32X2_IP(x20, py0, sizeof(ae_int32x2));
            AE_S32X2_IP(x30, py0, sizeof(ae_int32x2));
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

            AE_L64X2_IP(tmp00, tmp01, castxcc(ae_int64x2, px0), sizeof(ae_int32x4));
            AE_L64X2_IP(tmp10, tmp11, castxcc(ae_int64x2, px1), sizeof(ae_int32x4));
            AE_L64X2_IP(tmp20, tmp21, castxcc(ae_int64x2, px2), sizeof(ae_int32x4));
            AE_L64X2_IP(tmp30, tmp31, castxcc(ae_int64x2, px3), sizeof(ae_int32x4));

            x00 = AE_MOVINT32X2_FROMINT64(tmp00);
            x10 = AE_MOVINT32X2_FROMINT64(tmp10);
            x20 = AE_MOVINT32X2_FROMINT64(tmp20);
            x30 = AE_MOVINT32X2_FROMINT64(tmp30);
            x01 = AE_MOVINT32X2_FROMINT64(tmp01);
            x11 = AE_MOVINT32X2_FROMINT64(tmp11);
            x21 = AE_MOVINT32X2_FROMINT64(tmp21);
            x31 = AE_MOVINT32X2_FROMINT64(tmp31);

            DFT4X1RNG(x00, x10, x20, x30);
            DFT4X1RNG(x01, x11, x21, x31);

            x10 = AE_MULFC32RAS(x10, tw10);
            x20 = AE_MULFC32RAS(x20, tw20);
            x30 = AE_MULFC32RAS(x30, tw30);
            x11 = AE_MULFC32RAS(x11, tw11);
            x21 = AE_MULFC32RAS(x21, tw21);
            x31 = AE_MULFC32RAS(x31, tw31);

            AE_S32X2X2_IP(x00, x10, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
            AE_S32X2X2_IP(x20, x30, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
            AE_S32X2X2_IP(x01, x11, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
            AE_S32X2X2_IP(x21, x31, castxcc(ae_int32x4, py0), sizeof(ae_int32x4));
        } /* for (i = 0; i < (stride>>1); i++) */
    } /*if(stride&1) else... */

    *v = v[0] * R;
    return 3;
} /* ifft_stageS3_DFT4_first_32x32() */

#if 0
ALIGN(32) static const int32_t __fft8_tw1[] =
{
  (int32_t)0x7FFFFFFF, (int32_t)0x00000000,
  (int32_t)0x2D413CCD, (int32_t)0xD2BEC333,
  (int32_t)0x00000000, (int32_t)0xC0000000,
  (int32_t)0xD2BEC333, (int32_t)0xD2BEC333,
};
#else
/* exp(-1j*2*pi/8*(1:3)) */
ALIGN(32) static const int32_t __fft8_tw1[] =
{
    (int32_t)0x5A82799A, (int32_t)0xA57D8666, 
    (int32_t)0x00000000, (int32_t)0x80000000, 
    (int32_t)0xA57D8666, (int32_t)0xA57D8666
}; 
#endif

int fft_stageS3_DFT8_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  const int R = 8; // stage radix
  int j;
  const int stride = N / R;
  ae_int32x2 * restrict px0;
  ae_int32x2 * restrict py0;
  ae_int32x2 tw1, tw2, tw3;
  const ae_int32x2 * restrict ptwd;
  int shift;
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT(v[0] == stride);

  shift = 2;
  WUR_AE_SAR(shift>>1);
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
        ae_int32x2 s0, s1, s2, s3;
        ae_int32x2 d0, d1, d2, d3;

        AE_L32X2_XP(x7, px0, -1*(sizeof(ae_int32x2) * stride) );
        AE_L32X2_XP(x6, px0, -1*(sizeof(ae_int32x2) * stride) );
        AE_L32X2_XP(x5, px0, -1*(sizeof(ae_int32x2) * stride) );
        AE_L32X2_XP(x4, px0, -1*(sizeof(ae_int32x2) * stride) );
        AE_L32X2_XP(x3, px0, -1*(sizeof(ae_int32x2) * stride) );
        AE_L32X2_XP(x2, px0, -1*(sizeof(ae_int32x2) * stride) );
        AE_L32X2_XP(x1, px0, -1*(sizeof(ae_int32x2) * stride) );
        AE_L32X2_XP(x0, px0, sizeof(ae_int32x2) *( 7*stride +1 ));

        DFT4X1RNG2(x0, x2, x4, x6); 
        DFT4X1RNG2(x1, x3, x5, x7);

        x3 = AE_MULFC32RAS(x3, tw1); 
        x5 = AE_MULFC32RAS(x5, tw2);
        x7 = AE_MULFC32RAS(x7, tw3);

        AE_ADDANDSUBRNG32(s0, d0, x0, x1);
        AE_ADDANDSUBRNG32(s1, d1, x2, x3);
        AE_ADDANDSUBRNG32(s2, d2, x4, x5);
        AE_ADDANDSUBRNG32(s3, d3, x6, x7);

        x0 = s0; x4 = d0; x2 = s2; x6 = d2;
        x1 = s1; x5 = d1; x3 = s3; x7 = d3;

        AE_S32X2_XP(x0, py0, stride * sizeof(ae_int32x2));
        AE_S32X2_XP(x1, py0, stride * sizeof(ae_int32x2));
        AE_S32X2_XP(x2, py0, stride * sizeof(ae_int32x2));
        AE_S32X2_XP(x3, py0, stride * sizeof(ae_int32x2));
        AE_S32X2_XP(x4, py0, stride * sizeof(ae_int32x2));
        AE_S32X2_XP(x5, py0, stride * sizeof(ae_int32x2));
        AE_S32X2_XP(x6, py0, stride * sizeof(ae_int32x2));
        AE_S32X2_XP(x7, py0, (-7 * stride + 1)* sizeof(ae_int32x2));
      }
  }
  else /* if ((stride & 1) > 0) */
  {
      /* Even loop count,  unrolled twice  */
      __Pragma("loop_count min=1");
      for (j = 0; j < (stride>>1); j++)
      {
          /*
          12 cycles per pipeline stage in steady state with unroll=1
          2 pipeline stages
          */
          ae_int32x2 x00, x10, x20, x30;
          ae_int32x2 x40, x50, x60, x70;
          ae_int32x2 x01, x11, x21, x31;
          ae_int32x2 x41, x51, x61, x71;

          ae_int32x2 s0, s1, s2, s3;
          ae_int32x2 d0, d1, d2, d3;
                        
          AE_L32X2X2_XP(x70, x71, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
          AE_L32X2X2_XP(x60, x61, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
          AE_L32X2X2_XP(x50, x51, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
          AE_L32X2X2_XP(x40, x41, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
          AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
          AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
          AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
          AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), sizeof(ae_int32x2)*(7 * stride + 2));
          /**************** DFT8 ***************/
          DFT4X1RNG2(x00, x20, x40, x60);
          DFT4X1RNG2(x10, x30, x50, x70);

          x30 = AE_MULFC32RAS(x30, tw1);
          x50 = AE_MULFC32RAS(x50, tw2);
          x70 = AE_MULFC32RAS(x70, tw3);

          AE_ADDANDSUBRNG32(s0, d0, x00, x10);
          AE_ADDANDSUBRNG32(s1, d1, x20, x30);
          AE_ADDANDSUBRNG32(s2, d2, x40, x50);
          AE_ADDANDSUBRNG32(s3, d3, x60, x70);

          x00 = s0; x40 = d0; x20 = s2; x60 = d2;
          x10 = s1; x50 = d1; x30 = s3; x70 = d3;
          /**************** DFT8 ***************/
          DFT4X1RNG2(x01, x21, x41, x61);
          DFT4X1RNG2(x11, x31, x51, x71);

          x31 = AE_MULFC32RAS(x31, tw1);
          x51 = AE_MULFC32RAS(x51, tw2);
          x71 = AE_MULFC32RAS(x71, tw3);

          AE_ADDANDSUBRNG32(s0, d0, x01, x11);
          AE_ADDANDSUBRNG32(s1, d1, x21, x31);
          AE_ADDANDSUBRNG32(s2, d2, x41, x51);
          AE_ADDANDSUBRNG32(s3, d3, x61, x71);

          x01 = s0; x41 = d0; x21 = s2; x61 = d2;
          x11 = s1; x51 = d1; x31 = s3; x71 = d3;

          AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
          AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
          AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
          AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
          AE_S32X2X2_XP(x40, x41, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
          AE_S32X2X2_XP(x50, x51, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
          AE_S32X2X2_XP(x60, x61, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
          AE_S32X2X2_XP(x70, x71, castxcc(ae_int32x4, py0), (-7 * stride + 2)* sizeof(ae_int32x2));
      }
  }    /* if ((stride & 1) > 0) else..  */

  *v *= R;
  return (shift + 1);
} /* fft_stageS3_DFT8_last_32x32() */

int ifft_stageS3_DFT8_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    const int R = 8; // stage radix
    int j;
    const int stride = N / R;
    ae_int32x2 * restrict px0;
    ae_int32x2 * restrict py0;
    ae_int32x2 tw1, tw2, tw3;
    const ae_int32x2 * restrict ptwd;
    int shift;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(v[0] == stride);

    shift = 2;
    WUR_AE_SAR(shift >> 1);
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
            ae_int32x2 s0, s1, s2, s3;
            ae_int32x2 d0, d1, d2, d3;

            AE_L32X2_XP(x7, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x6, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x5, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x4, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x3, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x2, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x1, px0, -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2_XP(x0, px0, sizeof(ae_int32x2)*(7 * stride + 1));

            DFT4X1RNG2(x0, x2, x4, x6);
            DFT4X1RNG2(x1, x3, x5, x7);

            x3 = AE_MULFC32RAS(x3, tw1);
            x5 = AE_MULFC32RAS(x5, tw2);
            x7 = AE_MULFC32RAS(x7, tw3);

            AE_ADDANDSUBRNG32(s0, d0, x0, x1);
            AE_ADDANDSUBRNG32(s1, d1, x2, x3);
            AE_ADDANDSUBRNG32(s2, d2, x4, x5);
            AE_ADDANDSUBRNG32(s3, d3, x6, x7);

            x0 = s0; x4 = d0; x2 = s2; x6 = d2;
            x1 = s1; x5 = d1; x3 = s3; x7 = d3;

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
    else /* if ((stride & 1) > 0) */
    {
        /* Even loop count,  unrolled twice  */
        __Pragma("loop_count min=1");
        for (j = 0; j < (stride >> 1); j++)
        {
            /*
            12 cycles per pipeline stage in steady state with unroll=1
            2 pipeline stages
            */
            ae_int32x2 x00, x10, x20, x30;
            ae_int32x2 x40, x50, x60, x70;
            ae_int32x2 x01, x11, x21, x31;
            ae_int32x2 x41, x51, x61, x71;

            ae_int32x2 s0, s1, s2, s3;
            ae_int32x2 d0, d1, d2, d3;

            AE_L32X2X2_XP(x70, x71, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x60, x61, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x50, x51, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x40, x41, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), -1 * (sizeof(ae_int32x2)* stride));
            AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), sizeof(ae_int32x2)*(7 * stride + 2));
            /**************** DFT8 ***************/
            DFT4X1RNG2(x00, x20, x40, x60);
            DFT4X1RNG2(x10, x30, x50, x70);

            x30 = AE_MULFC32RAS(x30, tw1);
            x50 = AE_MULFC32RAS(x50, tw2);
            x70 = AE_MULFC32RAS(x70, tw3);

            AE_ADDANDSUBRNG32(s0, d0, x00, x10);
            AE_ADDANDSUBRNG32(s1, d1, x20, x30);
            AE_ADDANDSUBRNG32(s2, d2, x40, x50);
            AE_ADDANDSUBRNG32(s3, d3, x60, x70);

            x00 = s0; x40 = d0; x20 = s2; x60 = d2;
            x10 = s1; x50 = d1; x30 = s3; x70 = d3;
            /**************** DFT8 ***************/
            DFT4X1RNG2(x01, x21, x41, x61);
            DFT4X1RNG2(x11, x31, x51, x71);

            x31 = AE_MULFC32RAS(x31, tw1);
            x51 = AE_MULFC32RAS(x51, tw2);
            x71 = AE_MULFC32RAS(x71, tw3);

            AE_ADDANDSUBRNG32(s0, d0, x01, x11);
            AE_ADDANDSUBRNG32(s1, d1, x21, x31);
            AE_ADDANDSUBRNG32(s2, d2, x41, x51);
            AE_ADDANDSUBRNG32(s3, d3, x61, x71);

            x01 = s0; x41 = d0; x21 = s2; x61 = d2;
            x11 = s1; x51 = d1; x31 = s3; x71 = d3;

            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x00), AE_MOVINT64_FROMINT32X2(x01), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x10), AE_MOVINT64_FROMINT32X2(x11), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x20), AE_MOVINT64_FROMINT32X2(x21), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x30), AE_MOVINT64_FROMINT32X2(x31), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x40), AE_MOVINT64_FROMINT32X2(x41), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x50), AE_MOVINT64_FROMINT32X2(x51), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x60), AE_MOVINT64_FROMINT32X2(x61), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x70), AE_MOVINT64_FROMINT32X2(x71), castxcc(ae_int64x2, py0), (-7 * stride + 2)* sizeof(ae_int32x2));
        }
    }    /* if ((stride & 1) > 0) else..  */

    *v *= R;
    return (shift + 1);
} /* ifft_stageS3_DFT8_last_32x32() */

int ifft_stageS3_DFT4_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp){ return 0; }

ALIGN(32) static const int32_t __fft16_tw1[] =
{
    (int32_t)0x7FFFFFFF, (int32_t)0x00000000, (int32_t)0x7FFFFFFF, (int32_t)0x00000000, (int32_t)0x7FFFFFFF, (int32_t)0x00000000, 
    (int32_t)0x7641AF3D, (int32_t)0xCF043AB3, //t5,
    (int32_t)0x5A82799A, (int32_t)0xA57D8666, //t9,
    (int32_t)0x30FBC54D, (int32_t)0x89BE50C3, //t13,
    (int32_t)0x5A82799A, (int32_t)0xA57D8666, //t6,
    (int32_t)0x00000000, (int32_t)0x80000000, //t10,
    (int32_t)0xA57D8666, (int32_t)0xA57D8666, //t14,
    (int32_t)0x30FBC54D, (int32_t)0x89BE50C3, //t7,
    (int32_t)0xA57D8666, (int32_t)0xA57D8666, //t11,
    (int32_t)0x89BE50C3, (int32_t)0x30FBC54D, //t15,
};

/* radix-4 butterfly with normalization */
#define DFT4X1RNG_L(x0, x1, x2, x3)\
{   \
    ae_int32x2 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG32_L(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32_L(s1, d1, x1, x3); \
    AE_ADDANDSUB32S(x0, x2, s0, s1);   \
    AE_ADDANDSUB32JS(x3, x1, d0, d1);   \
}

/* radix-4 butterfly with normalization */
#define DFT4X1RNG_H(x0, x1, x2, x3)\
{   \
    ae_int32x2 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG32_H(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32_H(s1, d1, x1, x3); \
    AE_ADDANDSUB32S(x0, x2, s0, s1);   \
    AE_ADDANDSUB32JS(x3, x1, d0, d1);   \
}

int fft_stageS3_DFT16_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    const ae_int32x2 * restrict px0;
    ae_int32x2 * restrict py0;
    ae_int32x2 * restrict ptw = (ae_int32x2*)(__fft16_tw1 + 6);

    const int stride = (N >> 4);
    int j;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);

    px0 = (const ae_int32x2 *)x;
    py0 = (ae_int32x2 *)y;

    ae_int32x2 x0, x1, x2, x3;
    ae_int32x2 x4, x5, x6, x7;
    ae_int32x2 x8, x9, x10, x11;
    ae_int32x2 x12, x13, x14, x15;
    ae_int32x2 t5, t9,  t13;
    ae_int32x2 t6, t10, t14;
    ae_int32x2 t7, t11, t15;

    AE_L32X2_IP(t5, ptw, sizeof(complex_fract32)); 
    AE_L32X2_IP(t9, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t13, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t6, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t10, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t14, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t7, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t11, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t15, ptw, sizeof(complex_fract32));
    
    WUR_AE_SAR(2); 

    __Pragma("loop_count min=1");
    for (j = 0; j < stride; j++)
    {
        AE_L32X2_XP(x0,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x1,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x2,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x3,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x4,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x5,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x6,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x7,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x8,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x9,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x10,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x11,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x12,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x13,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x14,  px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x15,  px0, (1 - 15 * stride)*sizeof(complex_fract32));

        DFT4X1RNG(x0, x4, x8 , x12);
        DFT4X1RNG(x1, x5, x9 , x13);
        DFT4X1RNG(x2, x6, x10, x14);
        DFT4X1RNG(x3, x7, x11, x15);

        x5  = AE_MULFC32RAS(x5, t5);
        x6  = AE_MULFC32RAS(x6, t6);
        x7  = AE_MULFC32RAS(x7, t7);
        x9  = AE_MULFC32RAS(x9, t6);
        x10 = AE_MULFC32RAS(x10, t10);
        x11 = AE_MULFC32RAS(x11, t11);
        x13 = AE_MULFC32RAS(x13, t7);
        x14 = AE_MULFC32RAS(x14, t11);
        x15 = AE_MULFC32RAS(x15, t15);

        DFT4X1RNG(x0, x1, x2, x3);
        DFT4X1RNG(x4, x5, x6, x7);
        DFT4X1RNG(x8, x9, x10, x11);
        DFT4X1RNG(x12, x13, x14, x15);

        AE_S32X2_XP(x0 , py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x4 , py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x8 , py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x12, py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x1 , py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x5 , py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x9 , py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x13, py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x2 , py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x6 , py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x10,  py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x14,  py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x3 ,  py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x7 ,  py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x11,  py0, stride*sizeof(complex_fract32));
        AE_S32X2_XP(x15,  py0, (1 - 15 * stride)*sizeof(complex_fract32));
    }

    return 4;

} /* fft_stageS3_DFT16_last_32x32() */

int ifft_stageS3_DFT16_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    const ae_int32x2 * restrict px0;
    ae_int64 * restrict py0;
    ae_int32x2 * restrict ptw = (ae_int32x2*)(__fft16_tw1 + 6);

    const int stride = (N >> 4);
    int j;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);

    px0 = (const ae_int32x2 *)x;
    py0 = (ae_int64 *)y;

    ae_int32x2 x0, x1, x2, x3;
    ae_int32x2 x4, x5, x6, x7;
    ae_int32x2 x8, x9, x10, x11;
    ae_int32x2 x12, x13, x14, x15;
    ae_int32x2 t5, t9, t13;
    ae_int32x2 t6, t10, t14;
    ae_int32x2 t7, t11, t15;

    AE_L32X2_IP(t5, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t9, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t13, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t6, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t10, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t14, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t7, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t11, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t15, ptw, sizeof(complex_fract32));

    WUR_AE_SAR(2);

    __Pragma("loop_count min=1");
    for (j = 0; j < stride; j++)
    {
        AE_L32X2_XP(x0, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x1, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x2, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x3, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x4, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x5, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x6, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x7, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x8, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x9, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x10, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x11, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x12, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x13, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x14, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x15, px0, (1 - 15 * stride)*sizeof(complex_fract32));

        DFT4X1RNG(x0, x4, x8, x12);
        DFT4X1RNG(x1, x5, x9, x13);
        DFT4X1RNG(x2, x6, x10, x14);
        DFT4X1RNG(x3, x7, x11, x15);

        x5 = AE_MULFC32RAS(x5, t5);
        x6 = AE_MULFC32RAS(x6, t6);
        x7 = AE_MULFC32RAS(x7, t7);
        x9 = AE_MULFC32RAS(x9, t6);
        x10 = AE_MULFC32RAS(x10, t10);
        x11 = AE_MULFC32RAS(x11, t11);
        x13 = AE_MULFC32RAS(x13, t7);
        x14 = AE_MULFC32RAS(x14, t11);
        x15 = AE_MULFC32RAS(x15, t15);

        DFT4X1RNG(x0, x1, x2, x3);
        DFT4X1RNG(x4, x5, x6, x7);
        DFT4X1RNG(x8, x9, x10, x11);
        DFT4X1RNG(x12, x13, x14, x15);

        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x0), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x4), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x8), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x12), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x1), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x5), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x9), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x13), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x2), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x6), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x10), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x14), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x3), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x7), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x11), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x15), py0, (1 - 15 * stride)*sizeof(complex_fract32));
    }

    return 4;

} /* ifft_stageS3_DFT16_last_32x32() */


