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
DFT3 algorithm:
x - input complex vector
y - output complex vector
y = fft(x)
y = [ x(1) + x(2)  + x(3);
x(1) + (x(2) + x(3))*cos(2*pi/3) - 1j*(x(2) - x(3))*sin(2*pi/3);
x(1) + (x(2) + x(3))*cos(2*pi/3) + 1j*(x(2) - x(3))*sin(2*pi/3) ]

*/
#define DFT3X1(x0, x1, x2)\
{\
    ae_int32x2 s0, s1, d0, c32; \
    ae_int16x4 c;               \
    c32 = AE_MOVDA32X2(0x40004000,0x00006EDA); \
    c = AE_MOVINT16X4_FROMINT32X2(c32);        \
    AE_ADDANDSUB32S(s0, d0, x1, x2);           \
    s1 = AE_ADD32S(x0, s0);                    \
    AE_MULSFP32X16X2RAS_H(x0, s0, c);          \
    d0 = AE_MULFC32X16RAS_L(d0, c);            \
    s0 = x0;                                   \
    x0 = s1;                                   \
    AE_ADDANDSUB32S(x2, x1, s0, d0);           \
}

/*
 *  32x16 FFT first stage Radix 3, scalingOption=2
 */
int fft_32x16_stage_first_scl2_DFT3(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  const int R = 3; // stage radix
  const int stride = N / R;
  int shift;
  ae_int32x2 * restrict px0;
  ae_int32x2 * restrict px1;
  ae_int32x2 * restrict px2;
  ae_int32x2 * restrict py0;
  const ae_int16x4 * restrict ptwd;
  int min_shift = 3;
  int i;
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  px0 = (ae_int32x2 *)x;
  px1 = px0 + stride;
  px2 = px1 + stride;
  py0 = (ae_int32x2 *)y;
  ptwd = (const ae_int16x4 *)tw;

  shift = min_shift - *bexp;
  NASSERT(shift>-32 && shift<32);
  WUR_AE_SAR(0);

  __Pragma("loop_count min=2");
  for (i = 0; i <stride; i++)
  {
    ae_int32x2 x0, x1, x2;
    ae_int16x4 tw12;

    AE_L16X4_IP(tw12, ptwd, sizeof(ae_int16x4));

    AE_L32X2_IP(x0, px0, sizeof(ae_int32x2));
    AE_L32X2_IP(x1, px1, sizeof(ae_int32x2));
    AE_L32X2_IP(x2, px2, sizeof(ae_int32x2));
    x0 = AE_SRAA32RS(x0, shift);
    x1 = AE_SRAA32RS(x1, shift);
    x2 = AE_SRAA32RS(x2, shift);

    DFT3X1(x0, x1, x2);
    x1 = AE_MULFC32X16RAS_H(x1, tw12);
    x2 = AE_MULFC32X16RAS_L(x2, tw12);

    AE_S32X2RNG_IP(x0, py0, sizeof(ae_int32x2));
    AE_S32X2RNG_IP(x1, py0, sizeof(ae_int32x2));
    AE_S32X2RNG_IP(x2, py0, sizeof(ae_int32x2));
  }

  AE_CALCRNG3();
  *bexp = 3 - RUR_AE_SAR();
  *v *= R;
  return shift;
} /* fft_32x16_stage_first_scl2_DFT3() */
#endif
