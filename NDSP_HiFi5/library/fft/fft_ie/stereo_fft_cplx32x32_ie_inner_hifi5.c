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
    C code optimized for HiFi4
    Integrit, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"

#define SWAP_PTR(_x, _y) {complex_fract32 *tmp = _x; _x = _y ; _y = tmp; } 

/* Radix-4 butterfly with normalization   *
 * x0, x1, x2, x3 - input/output samples  *
 * AE_SAR - contains normalization factor */
#define DFT4X1RNG(x0, x1, x2, x3)\
{   \
    ae_int32x2 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG32(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32(s1, d1, x1, x3); \
    d1 = AE_MUL32JS(d1);               \
    AE_ADDANDSUB32S(x0, x2, s0, s1);   \
    AE_ADDANDSUB32S(x3, x1, d0, d1);   \
}

/* Radix-4 butterfly with normalization applied twice  *
 * x0, x1, x2, x3 - input/output samples               *
 * AE_SAR - contains normalization factor (used twice) */
#define DFT4X1RNG2(x0, x1, x2, x3)\
{   \
    ae_int32x2 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG32(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32(s1, d1, x1, x3); \
    d1 = AE_MUL32JS(d1);               \
    AE_ADDANDSUBRNG32(x0, x2, s0, s1); \
    AE_ADDANDSUBRNG32(x3, x1, d0, d1); \
}

/*-------------------------------------------------------------------------
  Internal stages of functions stereo_fft_cplx32x32_ie(),
  stereo_ifft_cplx32x32_ie()

  Performs all stages of FFT except the first one. It is assumed that
  the first stage performs radix-4 DIF butterflies with normal input/output
  order.
  Static scaling is applied.
 
  Input/Output:
  x[N]                intermediate signal after computing the first
                      stage of FFT/IFFT.

  Input:
  twd[N*tw_step*3/4]  twiddle factor table of a complex-valued FFT of 
                      size N*twdstep
  tw_step             twiddle step 
  N                   FFT size

  Temporary:
  y[N]                temporary buffer

  Returned value:     total number of right shifts occurred during 
                      scaling procedure

  Restrictions:
  x, y - should not overlap and must be aligned on 16-bytes boundary
-------------------------------------------------------------------------*/
int stereo_fft_cplx32x32_innerstage_scl3(complex_fract32 *x, complex_fract32 *y, const complex_fract32 *tw, int tw_step, int N)
{
  const ae_int32x2 * restrict px0;
  const ae_int32x2 * restrict px1;
  const ae_int32x2 * restrict px2;
  const ae_int32x2 * restrict px3;
        ae_int32x2 * restrict py0;
        ae_int32x2 * restrict py1;
        ae_int32x2 * restrict py2;
        ae_int32x2 * restrict py3;
  const ae_int32x2 * restrict ptwd;
  ae_int32x2 x0, x1, y0, y1;
  ae_int32x2 tw1, tw2, tw3;
  int shift = 0;
  int stridey;
  const int stridex = N>>2;
  const int log2N = 30 - NSA(N);
  int i, j, M;
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  stridey = 4;/* first stage should be Radix-4 */
  M = N>>2;
  WUR_AE_SAR(1);

  /*
   * Perform all stages except the first and last
   */
  while ( stridey < (N >> 2) )
  {
    M = M>>2;
    for (j = 0; j < (stridey>>1); j++)
    {
      ptwd = (const ae_int32x2 *)tw;
      px0 = (ae_int32x2 *)x + j*2*2;
      px1 = px0 + stridex*2;
      px2 = px1 + stridex*2;
      px3 = px2 + stridex*2;
      py0 = (ae_int32x2 *)y + j*2*2;
      py2 = py0 + 2;

      __Pragma( "loop_count min=2" );
      for (i = 0; i < M; i++)
      {
          ae_int32x2 x00, x10, x20, x30;
          ae_int32x2 x01, x11, x21, x31;
        /* load twiddles */
        AE_L32X2_IP(tw1, ptwd, sizeof(ae_int32x2));
        AE_L32X2_IP(tw2, ptwd, sizeof(ae_int32x2));
        AE_L32X2_XP(tw3, ptwd, (3*tw_step-2)*sizeof(ae_int32x2));

        /* butterfly 0, left and right channels */
        AE_L32X2X2_IP(x00, x01, castxcc(ae_int32x4,px0), 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x10, x11, castxcc(ae_int32x4,px1), 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x20, x21, castxcc(ae_int32x4,px2), 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x30, x31, castxcc(ae_int32x4,px3), 2*sizeof(ae_int32x2));

        DFT4X1RNG2(x00, x10, x20, x30);
        x10 = AE_MULFC32RAS(x10, tw1);
        x20 = AE_MULFC32RAS(x20, tw2);
        x30 = AE_MULFC32RAS(x30, tw3);

        DFT4X1RNG2(x01, x11, x21, x31);
        x11 = AE_MULFC32RAS(x11, tw1);
        x21 = AE_MULFC32RAS(x21, tw2);
        x31 = AE_MULFC32RAS(x31, tw3);

        AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), stridey * 2 * sizeof(ae_int32x2));
        AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), stridey * 2 * sizeof(ae_int32x2));
        AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), stridey * 2 * sizeof(ae_int32x2));
        AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), stridey * 2 * sizeof(ae_int32x2));

        /* butterfly 1, left and right channels */
        AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), (stridey * 2 - 2) * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px1), (stridey * 2 - 2) * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px2), (stridey * 2 - 2) * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px3), (stridey * 2 - 2) * sizeof(ae_int32x2));

        DFT4X1RNG2(x00, x10, x20, x30);
        x10 = AE_MULFC32RAS(x10, tw1);
        x20 = AE_MULFC32RAS(x20, tw2);
        x30 = AE_MULFC32RAS(x30, tw3);
        DFT4X1RNG2(x01, x11, x21, x31);
        x11 = AE_MULFC32RAS(x11, tw1);
        x21 = AE_MULFC32RAS(x21, tw2);
        x31 = AE_MULFC32RAS(x31, tw3);

        AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py2), stridey * 2 * sizeof(ae_int32x2));
        AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py2), stridey * 2 * sizeof(ae_int32x2));
        AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py2), stridey * 2 * sizeof(ae_int32x2));
        AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py2), stridey * 2 * sizeof(ae_int32x2));
      }
    }
    stridey = stridey<<2;
    tw_step = tw_step<<2;
    SWAP_PTR(x, y);
    shift += 2;
  }

  /*
   * Last stage
   */
  if ((log2N-1) & 2)
  {
    /* Execute the last stage inplace */
    y = x;
  }
  px0 = (const ae_int32x2 *)x;
  py0 = (      ae_int32x2 *)y;

  /* Radix-2 */
  if (log2N & 1)
  {
    const int stride = N;
    shift += 1;
    WUR_AE_SAR(1);
    __Pragma("loop_count min=6");
    for (j = 0; j < stride; j++)
    {
      x1 = AE_L32X2_X(px0, sizeof(ae_int32x2) * stride);
      AE_L32X2_IP(x0, px0, sizeof(ae_int32x2));

      AE_ADDANDSUBRNG32(y0, y1, x0, x1);

      AE_S32X2_XP(y0, py0, stride * sizeof(ae_int32x2));
      AE_S32X2_XP(y1, py0, (1-stride)* sizeof(ae_int32x2));
    }
  }
  /* Radix-4 */
  else
  {
    const int stride = N>>1;
    shift += 2;
    WUR_AE_SAR(1);
    px1 = px0 + stride;
    px2 = px1 + stride;
    px3 = px2 + stride;
    py1 = py0 + stride;
    py2 = py1 + stride;
    py3 = py2 + stride;
    __Pragma("loop_count min=1");
    for (j = 0; j < (stride>>1); j++)
    {
        ae_int32x2 x00, x10, x20, x30;
        ae_int32x2 x01, x11, x21, x31;
        AE_L32X2X2_IP(x00, x01, castxcc(ae_int32x4, px0), 2 * sizeof(ae_int32x2));
        AE_L32X2X2_IP(x10, x11, castxcc(ae_int32x4, px1), 2 * sizeof(ae_int32x2));
        AE_L32X2X2_IP(x20, x21, castxcc(ae_int32x4, px2), 2 * sizeof(ae_int32x2));
        AE_L32X2X2_IP(x30, x31, castxcc(ae_int32x4, px3), 2 * sizeof(ae_int32x2));

        DFT4X1RNG2(x00, x10, x20, x30);
        DFT4X1RNG2(x01, x11, x21, x31);

        AE_S32X2X2_IP(x00, x01, castxcc(ae_int32x4, py0), 2 * sizeof(ae_int32x2));
        AE_S32X2X2_IP(x10, x11, castxcc(ae_int32x4, py1), 2 * sizeof(ae_int32x2));
        AE_S32X2X2_IP(x20, x21, castxcc(ae_int32x4, py2), 2 * sizeof(ae_int32x2));
        AE_S32X2X2_IP(x30, x31, castxcc(ae_int32x4, py3), 2 * sizeof(ae_int32x2));
    }
  }

  return shift;
} /* stereo_fft_cplx32x32_innerstage_scl3() */


/*-------------------------------------------------------------------------
  Internal stages of functions stereo_fft_cplx32x32_ie(),
  stereo_ifft_cplx32x32_ie()

  Performs all stages of FFT except the first one. It is assumed that
  the first stage performs radix-4 DIF butterflies with normal input/output
  order.
  Autoscaling is applied.
 
  Input/Output:
  x[N]                intermediate signal after computing the first
                      stage of FFT/IFFT.

  Input:
  twd[N*tw_step*3/4]  twiddle factor table of a complex-valued FFT of 
                      size N*twdstep
  tw_step             twiddle step 
  N                   FFT size

  Temporary:
  y[N]                temporary buffer

  Returned value:     total number of right shifts occurred during 
                      scaling procedure

  Restrictions:
  x, y - should not overlap and must be aligned on 16-bytes boundary
-------------------------------------------------------------------------*/
int stereo_fft_cplx32x32_innerstage_scl2(complex_fract32 *x, complex_fract32 *y, const complex_fract32 *tw, int tw_step, int N, int bexp)
{
  const ae_int32x2 * restrict px0;
  const ae_int32x2 * restrict px1;
  const ae_int32x2 * restrict px2;
  const ae_int32x2 * restrict px3;
        ae_int32x2 * restrict py0;
        ae_int32x2 * restrict py1;
        ae_int32x2 * restrict py2;
        ae_int32x2 * restrict py3;
  const ae_int32x2 * restrict ptwd;
  int i, j, M;
  int stridey;
  const int stride = N>>2;
  int shift_curr, shift = 0;
  const int min_shift = 3;
  const int log2N = 30 - NSA(N); 

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  stridey = 4;/* first stage should be Radix-4 */
  M = N>>2;

  /*
   * Perform all stages except the first and last
   */
  while ( stridey < (N>>2) )
  {
    M = M>>2;
    shift_curr = min_shift-bexp;
    WUR_AE_SAR(shift_curr);
    shift += shift_curr;
    __Pragma("loop_count min=1");
    for (j = 0; j < (stridey>>1); j++)
    {
      ae_int32x2 tw1, tw2, tw3;
      ptwd = (const ae_int32x2 *)tw;
      px0 = (ae_int32x2 *)x + j*2*2;
      px1 = px0 + stride*2;
      px2 = px1 + stride*2;
      px3 = px2 + stride*2;
      py0 = (ae_int32x2 *)y + j*2*2;
      py1 = py0 + 1;
      py2 = py1 + 1;
      py3 = py2 + 1;

      __Pragma("loop_count min=2");
      for (i = 0; i < M; i++)
      {
          ae_int32x2 x00, x10, x20, x30;
          ae_int32x2 x01, x11, x21, x31;
          /* load twiddles */
          AE_L32X2_IP(tw1, ptwd, sizeof(ae_int32x2));
          AE_L32X2_IP(tw2, ptwd, sizeof(ae_int32x2));
          AE_L32X2_XP(tw3, ptwd, (3 * tw_step - 2)*sizeof(ae_int32x2));

          /* butterfly 0, left and right channels */
          AE_L32X2X2_IP(x00, x01, castxcc(ae_int32x4, px0), 2 * sizeof(ae_int32x2));
          AE_L32X2X2_IP(x10, x11, castxcc(ae_int32x4, px1), 2 * sizeof(ae_int32x2));
          AE_L32X2X2_IP(x20, x21, castxcc(ae_int32x4, px2), 2 * sizeof(ae_int32x2));
          AE_L32X2X2_IP(x30, x31, castxcc(ae_int32x4, px3), 2 * sizeof(ae_int32x2));

          DFT4X1RNG(x00, x10, x20, x30);
          x10 = AE_MULFC32RAS(x10, tw1);
          x20 = AE_MULFC32RAS(x20, tw2);
          x30 = AE_MULFC32RAS(x30, tw3);

          DFT4X1RNG(x01, x11, x21, x31);
          x11 = AE_MULFC32RAS(x11, tw1);
          x21 = AE_MULFC32RAS(x21, tw2);
          x31 = AE_MULFC32RAS(x31, tw3);

          AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), stridey * 2 * sizeof(ae_int32x2));
          AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), stridey * 2 * sizeof(ae_int32x2));
          AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), stridey * 2 * sizeof(ae_int32x2));
          AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), stridey * 2 * sizeof(ae_int32x2));

          /* butterfly 1, left and right channels */
          AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), (stridey * 2 - 2) * sizeof(ae_int32x2));
          AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px1), (stridey * 2 - 2) * sizeof(ae_int32x2));
          AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px2), (stridey * 2 - 2) * sizeof(ae_int32x2));
          AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px3), (stridey * 2 - 2) * sizeof(ae_int32x2));

          DFT4X1RNG(x00, x10, x20, x30);
          x10 = AE_MULFC32RAS(x10, tw1);
          x20 = AE_MULFC32RAS(x20, tw2);
          x30 = AE_MULFC32RAS(x30, tw3);
          DFT4X1RNG(x01, x11, x21, x31);
          x11 = AE_MULFC32RAS(x11, tw1);
          x21 = AE_MULFC32RAS(x21, tw2);
          x31 = AE_MULFC32RAS(x31, tw3);

          AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py2), stridey * 2 * sizeof(ae_int32x2));
          AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py2), stridey * 2 * sizeof(ae_int32x2));
          AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py2), stridey * 2 * sizeof(ae_int32x2));
          AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py2), stridey * 2 * sizeof(ae_int32x2));
      }
    }

    AE_CALCRNG3();
    bexp = 3 - RUR_AE_SAR();
    stridey *= 4;
    tw_step *= 4;
    SWAP_PTR(x, y);
  }

  /*
   * Last stage
   */
  if ((log2N-1) & 2)
  {
    /* Execute the last stage inplace */
    y = x;
  }
  px0 = (const ae_int32x2 *)x;
  py0 = (      ae_int32x2 *)y;

  /* radix-2 */
  if (log2N & 1)
  {
    const int stride = N;
    const int min_shift = 1;
    shift_curr = XT_MAX( 0, min_shift-bexp );
    WUR_AE_SAR(shift_curr);
    shift += shift_curr;

    __Pragma("loop_count min=3");
    for (j = 0; j < stride; j++)
    {
      ae_int32x2 x0, x1, y0, y1;

      AE_L32X2_XP(x0, px0, stride*sizeof(ae_int32x2));
      AE_L32X2_XP(x1, px0, (1-stride)*sizeof(ae_int32x2));

      AE_ADDANDSUBRNG32(y0, y1, x0, x1);

      AE_S32X2RNG_XP(y0, py0, stride*sizeof(ae_int32x2));
      AE_S32X2RNG_XP(y1, py0, (1-stride)*sizeof(ae_int32x2));
    }
  }
  /* radix-4 */
  else
  {
    const int stride = N>>1;
    const int min_shift = 3;
    shift_curr = XT_MAX( 0, min_shift-bexp );
    WUR_AE_SAR(shift_curr);
    shift += shift_curr;
    px1 = px0 + stride;
    px2 = px1 + stride;
    px3 = px2 + stride;
    py1 = py0 + stride;
    py2 = py1 + stride;
    py3 = py2 + stride;
    __Pragma("loop_count min=1");
    for (j = 0; j < (stride >> 1); j++)
    {
        ae_int32x2 x00, x10, x20, x30;
        ae_int32x2 x01, x11, x21, x31;
        AE_L32X2X2_IP(x00, x01, castxcc(ae_int32x4, px0), 2 * sizeof(ae_int32x2));
        AE_L32X2X2_IP(x10, x11, castxcc(ae_int32x4, px1), 2 * sizeof(ae_int32x2));
        AE_L32X2X2_IP(x20, x21, castxcc(ae_int32x4, px2), 2 * sizeof(ae_int32x2));
        AE_L32X2X2_IP(x30, x31, castxcc(ae_int32x4, px3), 2 * sizeof(ae_int32x2));

        DFT4X1RNG(x00, x10, x20, x30);
        DFT4X1RNG(x01, x11, x21, x31);

        AE_S32X2X2RNG_IP(x00, x01, castxcc(ae_int32x4, py0), 2 * sizeof(ae_int32x2));
        AE_S32X2X2RNG_IP(x10, x11, castxcc(ae_int32x4, py1), 2 * sizeof(ae_int32x2));
        AE_S32X2X2RNG_IP(x20, x21, castxcc(ae_int32x4, py2), 2 * sizeof(ae_int32x2));
        AE_S32X2X2RNG_IP(x30, x31, castxcc(ae_int32x4, py3), 2 * sizeof(ae_int32x2));
    }
  }

  return shift;
} /* stereo_fft_cplx32x32_innerstage_scl2() */
