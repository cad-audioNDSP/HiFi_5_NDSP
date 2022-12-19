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
    with dynamic data scaling: 32-bit data, 32-bit twiddle factors
    C code optimized for HiFi4
  IntegrIT, 2006-2019
*/
#include "NatureDSP_Signal_fft.h"
#include "fft_twiddles32x32.h"
#include "common.h"
#include <xtensa/config/core-isa.h>


inline_ void _cmult32x32(ae_int32x2 *result, ae_int32x2 *x, ae_int32x2 *y)
{
  ae_f32x2 z;

  z = AE_MULFC32RAS(AE_MOVF32X2_FROMINT32X2(*x), AE_MOVF32X2_FROMINT32X2(*y));
  *result = AE_MOVINT32X2_FROMF32X2(z);
}

/*
DFT3 algorithm:
x - input complex vector
y - output complex vector
y = fft(x)
y = [ x(1) + x(2)  + x(3);
x(1) + (x(2) + x(3))*cos(2*pi/3) - 1j*(x(2) - x(3))*sin(2*pi/3);
x(1) + (x(2) + x(3))*cos(2*pi/3) + 1j*(x(2) - x(3))*sin(2*pi/3) ]

*/
/* DFT3 without scaling */
#if 1
inline_ void __dft3x1(ae_int32x2 *x0, ae_int32x2 *x1, ae_int32x2 *x2)
{
    ae_int32x2 s0, s1, d0;
    ae_int32x2 c;
    c = AE_MOVDA32X2(0x0, 0x6ED9EBA1);
    AE_ADDANDSUB32S(s0, d0, *x1, *x2);
    s1 = AE_ADD32S(*x0, s0);
    s0 = AE_SRAI32(s0, 1);
    _cmult32x32(&d0, &d0, &c);
    s0 = AE_SUB32S(*x0, s0);
    *x0 = s1;
    AE_ADDANDSUB32S(*x2, *x1, s0, d0);
}
#define DFT3X1(_x0, _x1, _x2) __dft3x1(&_x0, &_x1, &_x2 )
#else
#define DFT3X1(x0, x1, x2)\
{\
    ae_int32x2 s0, s1, d0;            \
    ae_int32x2 c;                     \
    c = AE_MOVDA32X2(0x0,0x6ED9EBA1); \
    AE_ADDANDSUB32S(s0, d0, x1, x2);  \
    s1 = AE_ADD32S(x0, s0);           \
    s0 = AE_SRAI32(s0, 1);            \
    _cmult32x32(&d0, &d0, &c);        \
    s0 = AE_SUB32S(x0, s0);           \
    x0 = s1;                          \
    AE_ADDANDSUB32S(x2, x1, s0, d0);  \
}
#endif
/* radix-3 butterfly with data scaling using AE.SAR register & 'shift' value */
#define DFT3X1_RNG_shift(x0, x1, x2, shift)\
{\
    ae_int32x2 s0, s1, d0;            \
    ae_int32x2 c;                     \
    c = AE_MOVDA32X2(0x0,0x6ED9EBA1); \
    x0 = AE_SRAA32(x0, shift);        \
    AE_ADDANDSUBRNG32(s0, d0, x1, x2);\
    s1 = AE_ADD32S(x0, s0);           \
    s0 = AE_SRAI32(s0, 1);            \
    _cmult32x32(&d0, &d0, &c);        \
    s0 = AE_SUB32S(x0, s0);           \
    x0 = s1;                          \
    AE_ADDANDSUB32S(x2, x1, s0, d0);  \
}

/*
DFT5 algorithm, No scaling:
x - input complex vector
y - output complex vector
y = fft(x)
w1 =  exp(-1j*2*pi/5);
w2 =  exp(-1j*2*pi*2/5);

y = zeros(5,1);
s1 = (x1+x4);
s2 = (x2 + x3);
d1 = (x1-x4);
d2 = (x2-x3);

y(1) = x0 + s1 + s2;
y(2) = x0 + (s1*real(w1) + s2*real(w2)) + 1j*(d1*imag(w1) + d2*imag(w2));
y(5) = x0 + (s1*real(w1) + s2*real(w2)) - 1j*(d1*imag(w1) + d2*imag(w2));
y(3) = x0 + (s1*real(w2) + s2*real(w1)) + 1j*(d1*imag(w2)  - d2*imag(w1));
y(4) = x0 + (s1*real(w2) + s2*real(w1)) - 1j*(d1*imag(w2)  - d2*imag(w1));

*/
#define DFT5X1(x0, x1, x2, x3, x4)\
{ \
  ae_int32x2 s1, s2, d1, d2;                                 \
  ae_int32x2 y0, y1, y2        ;                             \
  ae_int32x2 t0, t1, t2, t3;                                 \
  ae_int32x2 real_w1, jimag_w1, real_w2, jimag_w2;           \
  real_w1 = AE_MOVDA32X2(0x278DDE6E, 0x0);                   \
  jimag_w1 = AE_MOVDA32X2(0x0,0x8643C7B3);                   \
  real_w2 = AE_MOVDA32X2(0x98722192,0x0);                    \
  jimag_w2 = AE_MOVDA32X2(0x0,0xB4C373EE);                   \
  AE_ADDANDSUB32S(s1, d1, x1, x4);                           \
  AE_ADDANDSUB32S(s2, d2, x2, x3);                           \
  _cmult32x32(&t0, &s1, &real_w1);                           \
  _cmult32x32(&t1, &s2, &real_w2);                           \
  _cmult32x32(&t2, &s1, &real_w2);                           \
  _cmult32x32(&t3, &s2, &real_w1);                           \
  y0 = AE_ADD32S(x0, AE_ADD32S(s1, s2));                     \
  y1 = AE_ADD32S(x0, AE_ADD32S(t0, t1));                     \
  y2 = AE_ADD32S(x0, AE_ADD32S(t2, t3));                     \
  _cmult32x32(&t0, &d1, &jimag_w1);                          \
  _cmult32x32(&t1, &d2, &jimag_w2);                          \
  _cmult32x32(&t2, &d2, &jimag_w1);                          \
  _cmult32x32(&t3, &d1, &jimag_w2);                          \
  t0 = AE_ADD32S(t0, t1);                                    \
  t1 = AE_SUB32S(t3, t2);                                    \
  x0 = y0;                                                   \
  AE_ADDANDSUB32S(x1, x4, y1, t0);                           \
  AE_ADDANDSUB32S(x2, x3, y2, t1);                           \
}                                                       



/*
DFT5 algorithm with scaling:
x - input complex vector
y - output complex vector
y = fft(x)
w1 =  exp(-1j*2*pi/5);
w2 =  exp(-1j*2*pi*2/5);

y = zeros(5,1);
s1 = (x1+x4);
s2 = (x2+x3);
d1 = (x1-x4);
d2 = (x2-x3);

y(1) = x0 + s1 + s2;
y(2) = x0 + (s1*real(w1) + s2*real(w2)) + 1j*(d1*imag(w1) + d2*imag(w2));
y(5) = x0 + (s1*real(w1) + s2*real(w2)) - 1j*(d1*imag(w1) + d2*imag(w2));
y(3) = x0 + (s1*real(w2) + s2*real(w1)) + 1j*(d1*imag(w2)  - d2*imag(w1));
y(4) = x0 + (s1*real(w2) + s2*real(w1)) - 1j*(d1*imag(w2)  - d2*imag(w1));

*/
#define DFT5X1S(x0, x1, x2, x3, x4, shift)                       \
{                                                                \
    ae_int32x2 s1, s2, d1, d2;                                   \
    ae_int32x2 y0, y1, y2;                                       \
    ae_int32x2 t0, t1, t2, t3;                                   \
    ae_int32x2 real_w1, jimag_w1, real_w2, jimag_w2;             \
    real_w1 = AE_MOVDA32X2(0x278DDE6E, 0x0);                     \
    jimag_w1 = AE_MOVDA32X2(0x0, 0x8643C7B3);                    \
    real_w2 = AE_MOVDA32X2(0x98722192, 0x0);                     \
    jimag_w2 = AE_MOVDA32X2(0x0, 0xB4C373EE);                    \
    x0 = AE_SRAA32RS(x0, shift);                                  \
    AE_ADDANDSUBRNG32(s1, d1, x1, x4);                           \
    AE_ADDANDSUBRNG32(s2, d2, x2, x3);                           \
    _cmult32x32(&t0, &s1, &real_w1);                             \
    _cmult32x32(&t1, &s2, &real_w2);                             \
    _cmult32x32(&t2, &s1, &real_w2);                             \
    _cmult32x32(&t3, &s2, &real_w1);                             \
    y0 = AE_ADD32S(x0, AE_ADD32S(s1, s2));                       \
    y1 = AE_ADD32S(x0, AE_ADD32S(t0, t1));                       \
    y2 = AE_ADD32S(x0, AE_ADD32S(t2, t3));                       \
    _cmult32x32(&t0, &d1, &jimag_w1);                            \
    _cmult32x32(&t1, &d2, &jimag_w2);                            \
    _cmult32x32(&t2, &d2, &jimag_w1);                            \
    _cmult32x32(&t3, &d1, &jimag_w2);                            \
    t0 = AE_ADD32S(t0, t1);                                      \
    t1 = AE_SUB32S(t3, t2);                                      \
    x0 = y0;                                                     \
    AE_ADDANDSUB32S(x1, x4, y1, t0);                             \
    AE_ADDANDSUB32S(x2, x3, y2, t1);                             \
}

/*
 *  32x32 FFT first stage Radix 2, scalingOption=2
 */
int fft_stageS2_DFT2_first_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  const ae_int32x2 * restrict px0;
  const ae_int32x2 * restrict px1;
        ae_int32x2 * restrict py0;
  const ae_int32x2 * restrict ptwd;
  const int stride = (N >> 1);
  int shift;
  const int R = 2; // stage radix
  const int min_shift = 2;
  int i;

  shift = min_shift - *bexp;

  NASSERT(shift>-32 && shift<4);
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT_ALIGN16(tw);
  NASSERT((stride&1)==1); 
  NASSERT((stride >> 1)>=3); 

  //scl = 1 << shiftl;
  WUR_AE_SAR(min_shift);

  px0 = (ae_int32x2 *)x;
  px1 = px0 + stride;
  py0 = (ae_int32x2 *)y;
  ptwd = (const ae_int32x2 *)tw;

  __Pragma("loop_count min=3"); 
  for (i = 0; i < (stride>>1); i++)
  {
    ae_int32x2 x00, x10, y00, y10, tw10;
    ae_int32x2 x01, x11, y01, y11, tw11;

    AE_L32X2X2_IP(tw10, tw11, castxcc(ae_int32x4, ptwd), 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(x00, x01, castxcc(ae_int32x4, px0), 2 * sizeof(ae_int32x2));
    /* px1 always unaligned by 16 bytes */
    AE_L32X2_IP(x10, px1, sizeof(ae_int32x2));
    AE_L32X2_IP(x11, px1, sizeof(ae_int32x2));

    x00 = AE_SLAA32(x00, *bexp); 
    x10 = AE_SLAA32(x10, *bexp); 
    x01 = AE_SLAA32(x01, *bexp);
    x11 = AE_SLAA32(x11, *bexp);

    AE_ADDANDSUBRNG32(y00, y10, x00, x10);
    AE_ADDANDSUBRNG32(y01, y11, x01, x11);

    y11 = AE_MULFC32RAS(y11, tw11);
    y10 = AE_MULFC32RAS(y10, tw10);

    AE_S32X2X2RNG_IP(y00, y10, castxcc(ae_int32x4, py0), 2*sizeof(ae_int32x2));
    AE_S32X2X2RNG_IP(y01, y11, castxcc(ae_int32x4, py0), 2*sizeof(ae_int32x2));
  }

  {
      ae_int32x2 x0, x1, y0, y1;
      ae_int32x2 tw1;

      AE_L32X2_IP(tw1, ptwd, sizeof(ae_int32x2));

      AE_L32X2_IP(x0, px0, sizeof(ae_int32x2));
      AE_L32X2_IP(x1, px1, sizeof(ae_int32x2));

      x0 = AE_SLAA32(x0, *bexp);
      x1 = AE_SLAA32(x1, *bexp);

      AE_ADDANDSUBRNG32(y0, y1, x0, x1);
      y1 = AE_MULFC32RAS(y1, tw1);

      AE_S32X2X2RNG_IP(y0, y1, castxcc(ae_int32x4, py0), 2*sizeof(ae_int32x2));

  }

  AE_CALCRNG3();
  *bexp = 3 - RUR_AE_SAR();
  *v *= R;
  return shift;
} /* fft_stageS2_DFT2_first_32x32() */

/*
 *  32x32 FFT last stage Radix 2, scalingOption=2
 */
int fft_stageS2_DFT2_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  const int stride = (N >> 1);
  int shift;
  ae_int32x2 * restrict px0;
  ae_int32x2 * restrict py0;
  int j;

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT(*bexp >= 0);

  shift = XT_MAX( 0, 1 - *bexp );
  WUR_AE_SAR(shift);
  px0 = (ae_int32x2 *)(x);
  py0 = (ae_int32x2 *)(y);

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

  return shift;
} /* fft_stageS2_DFT2_last_32x32() */


/*
*  32x32 FFT last stage Radix 2, scalingOption=2
*/
int ifft_stageS2_DFT2_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    const int stride = (N >> 1);
    int shift;
    ae_int32x2 * restrict px0;
    ae_int32x2 * restrict py0;
    int j;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(*bexp >= 0);

    shift = XT_MAX(0, 1 - *bexp);
    WUR_AE_SAR(shift);
    px0 = (ae_int32x2 *)(x);
    py0 = (ae_int32x2 *)(y);

    __Pragma("loop_count min=3");
    for (j = 0; j < stride; j++)
    {
        ae_int32x2 x0, x1, y0, y1;
        ae_int64 t0, t1;

        AE_L64_XP(t0, castxcc(ae_int64, px0), stride*sizeof(ae_int32x2));
        AE_L64_XP(t1, castxcc(ae_int64, px0), (1 - stride)*sizeof(ae_int32x2));

        x0 = AE_MOVINT32X2_FROMINT64(t0);
        x1 = AE_MOVINT32X2_FROMINT64(t1);

        AE_ADDANDSUBRNG32(y0, y1, x0, x1);

        AE_S32X2RNG_XP(y0, py0, stride*sizeof(ae_int32x2));
        AE_S32X2RNG_XP(y1, py0, (1 - stride)*sizeof(ae_int32x2));
    }

    return shift;
} /* fft_stageS2_DFT2_last_32x32() */

/*
 *  32x32 FFT/IFFT intermediate stage Radix 2, scalingOption=2
 */
int fft_stageS2_DFT2_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{ 
  int i, j, _v;
  int shift;
  const ae_int32x2 * restrict px0;
        ae_int32x2 * restrict py0;
  const ae_int32x2 * restrict ptwd;
  const int stride = (N >> 1);
  const int R = 2; // stage radix

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  shift = XT_MAX( 0, 2 - *bexp );
  WUR_AE_SAR(shift);

  _v = *v;
  for (j = 0; j < _v; j++)
  {
    ae_int32x2 x0, x1, y0, y1;
    ae_int32x2 tw1;
    ptwd = (const ae_int32x2 *)tw;
    px0 = (ae_int32x2 *)x + j;
    py0 = (ae_int32x2 *)y + j;
    __Pragma("loop_count min=1");
    for (i = 0; i < (stride / _v); i++)
    {
      AE_L32X2_XP(tw1, ptwd, tw_step*sizeof(ae_int32x2));

      x1 = AE_L32X2_X(px0, sizeof(ae_int32x2) * stride);
      AE_L32X2_XP(x0, px0, _v * sizeof(ae_int32x2));
      AE_ADDANDSUBRNG32(y0, y1, x0, x1);
      _cmult32x32(&y1, &y1, &tw1);
      AE_S32X2RNG_XP(y0, py0, _v * sizeof(ae_int32x2));
      AE_S32X2RNG_XP(y1, py0, _v * sizeof(ae_int32x2));
    }
  }
  AE_CALCRNG3();
  *bexp = 3 - RUR_AE_SAR();
  *v *= R;
  return shift;
} /* fft_stageS2_DFT2_32x32() */

/*
 *  32x32 FFT first stage Radix 3, scalingOption=2
 */
int fft_stageS2_DFT3_first_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  const int R = 3; // stage radix
  const int stride = N / R;
  int shift;
  ae_int32x2 * restrict px0;
  ae_int32x2 * restrict px1;
  ae_int32x2 * restrict px2;
  ae_int32x2 * restrict py0;
  const ae_int32x2 * restrict ptwd;
  int min_shift = 3;
  int i;
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  px0 = (ae_int32x2 *)x;
  px1 = px0 + stride;
  px2 = px1 + stride;
  py0 = (ae_int32x2 *)y;
  ptwd = (const ae_int32x2 *)tw;

  shift = min_shift - *bexp;
  NASSERT(shift>-32 && shift<32);
  WUR_AE_SAR(0);

  __Pragma("loop_count min=2");
  for (i = 0; i <stride; i++)
  {
    ae_int32x2 x0, x1, x2;
    ae_int32x2 tw1, tw2;

    AE_L32X2_IP(tw1, ptwd, sizeof(ae_int32x2));
    AE_L32X2_IP(tw2, ptwd, sizeof(ae_int32x2));

    AE_L32X2_IP(x0, px0, sizeof(ae_int32x2));
    AE_L32X2_IP(x1, px1, sizeof(ae_int32x2));
    AE_L32X2_IP(x2, px2, sizeof(ae_int32x2));
    x0 = AE_SRAA32RS(x0, shift);
    x1 = AE_SRAA32RS(x1, shift);
    x2 = AE_SRAA32RS(x2, shift);

    DFT3X1(x0, x1, x2);
    _cmult32x32(&x1, &x1, &tw1);
    _cmult32x32(&x2, &x2, &tw2);

    AE_S32X2RNG_IP(x0, py0, sizeof(ae_int32x2));
    AE_S32X2RNG_IP(x1, py0, sizeof(ae_int32x2));
    AE_S32X2RNG_IP(x2, py0, sizeof(ae_int32x2));
  }

  AE_CALCRNG3();
  *bexp = 3 - RUR_AE_SAR();
  *v *= R;
  return shift;
} /* fft_stageS2_DFT3_first_32x32() */

/*
 *  32x32 FFT last stage Radix 3, scalingOption=2
 */
int fft_stageS2_DFT3_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  const int R = 3; // stage radix
  const int stride = N / R;
  int shift;
  ae_int32x2 * restrict px0;
  ae_int32x2 * restrict py0;
  int j, _v;

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  _v = *v;
  px0 = (ae_int32x2 *)x;
  py0 = (ae_int32x2 *)y;

  shift = XT_MAX(0, 3 - *bexp);
  NASSERT(shift>=0 && shift<4);
  NASSERT((_v&1)==0); 

  WUR_AE_SAR(shift);

  __Pragma("loop_count min=1");
  for (j = 0; j < (_v >> 1); j++)
  {
      /*  11 cycles per pipeline stage in steady state with unroll=2 */
      ae_int32x2 x00, x10, x20;
      ae_int32x2 x01, x11, x21;

      AE_L32X2X2_X(x10, x11, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride);
      AE_L32X2X2_X(x20, x21, (ae_int32x4*)px0, 2 * sizeof(ae_int32x2)* stride);
      AE_L32X2X2_IP(x00, x01, castxcc(ae_int32x4, px0), 2 * sizeof(ae_int32x2));

      DFT3X1_RNG_shift(x00, x10, x20, shift);
      DFT3X1_RNG_shift(x01, x11, x21, shift);

      AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
      AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
      AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), (-2 * stride + 2) * sizeof(ae_int32x2));
  }

  return shift;
} /* fft_stageS2_DFT3_last_32x32() */


/*
*  32x32 FFT last stage Radix 3, scalingOption=2
*/
int ifft_stageS2_DFT3_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    const int R = 3; // stage radix
    const int stride = N / R;
    int shift;
    ae_int32x2 * restrict px0;
    ae_int32x2 * restrict py0;
    int j, _v;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);

    _v = *v;
    px0 = (ae_int32x2 *)x;
    py0 = (ae_int32x2 *)y;

    shift = XT_MAX(0, 3 - *bexp);
    NASSERT(shift >= 0 && shift<4);
    NASSERT((_v & 1) == 0);

    WUR_AE_SAR(shift);

    __Pragma("loop_count min=1");
    for (j = 0; j < (_v >> 1); j++)
    {
        /*  11 cycles per pipeline stage in steady state with unroll=2 */
        ae_int32x2 x00, x10, x20;
        ae_int32x2 x01, x11, x21;

        AE_L32X2X2_X(x10, x11, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride);
        AE_L32X2X2_X(x20, x21, (ae_int32x4*)px0, 2 * sizeof(ae_int32x2)* stride);
        AE_L32X2X2_IP(x00, x01, castxcc(ae_int32x4, px0), 2 * sizeof(ae_int32x2));

        DFT3X1_RNG_shift(x00, x10, x20, shift);
        DFT3X1_RNG_shift(x01, x11, x21, shift);

        AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x00), AE_MOVINT64_FROMINT32X2(x01), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
        AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x10), AE_MOVINT64_FROMINT32X2(x11), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
        AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x20), AE_MOVINT64_FROMINT32X2(x21), castxcc(ae_int64x2, py0), (-2 * stride + 2) * sizeof(ae_int32x2));
    }

    return shift;
} /* ifft_stageS2_DFT3_last_32x32() */

/*
 *  32x32 FFT/IFFT intermediate stage Radix 3, scalingOption=2
 */
int fft_stageS2_DFT3_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{ 
  const int R = 3; // stage radix
  const int stride = N / R;
  int shift;
        ae_int32x2 * restrict px0;
        ae_int32x2 * restrict py0;
  const ae_int32x2 * restrict ptwd;

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  shift = XT_MAX(0, 3 - *bexp);
  NASSERT(shift>=0 && shift<4);
  WUR_AE_SAR(shift);

  {
    int i, j;
    int _v;
    _v = *v;
    if ((_v & 1) == 0)
    {
        ae_int32x2 tw1, tw2;
        ae_int32x2 x00, x10, x20;
        ae_int32x2 x01, x11, x21;
        int flag = 0;

        py0 = (ae_int32x2 *)y;
        px0 = (ae_int32x2 *)x;
        ptwd = (const ae_int32x2 *)tw;

        __Pragma("loop_count min=2");
        for (i = 0; i< (stride >> 1); i++)
        {
            /*  8 cycles per pipeline stage in steady state with unroll=1 */
            int py_inc;
            int tw_inc;
            py_inc = (2 - _v * 2) * sizeof(ae_int32x2); //(-2 * _v + 1)* sizeof(ae_int32x2);
            tw_inc = 0;
        #if (XCHAL_HW_VERSION >= 281090)
            AE_ADDICIRC((int*)flag, _v<<2, 1<<3);
        #else
            AE_ADDICIRC(flag,  _v<<2, 1<<3);
        #endif 
            XT_MOVEQZ(py_inc, 2 * sizeof(ae_int32x2), flag);
            XT_MOVEQZ(tw_inc, 2 * sizeof(ae_int32x2), flag);

            AE_L32X2X2_XP(tw1, tw2, castxcc(ae_int32x4, ptwd), tw_inc);
            AE_L32X2X2_X(x10, x11, (ae_int32x4*)px0, 1 * stride*sizeof(ae_int32x2));
            AE_L32X2X2_X(x20, x21, (ae_int32x4*)px0, 2 * stride*sizeof(ae_int32x2));
            AE_L32X2X2_IP(x00, x01, castxcc(ae_int32x4, px0), 2 * sizeof(ae_int32x2));
            ////////////////////
            DFT3X1_RNG_shift(x00, x10, x20, shift);
            x10 = AE_MULFC32RAS(x10, tw1);
            x20 = AE_MULFC32RAS(x20, tw2);
            ///////////////////
            DFT3X1_RNG_shift(x01, x11, x21, shift);
            x11 = AE_MULFC32RAS(x11, tw1);
            x21 = AE_MULFC32RAS(x21, tw2);

            AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
            AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), py_inc);
        }
    }
    else /* if ((_v & 3) == 0) */
    {
        __Pragma("loop_count min=1");
        for (j = 0; j < _v; j++)
        {
            ae_int32x2 x0, x1, x2;
            ae_int32x2 tw1, tw2;
            ptwd = (const ae_int32x2 *)tw;
            px0 = (ae_int32x2 *)x + j;
            py0 = (ae_int32x2 *)y + j;
            __Pragma("loop_count min=1");
            for (i = 0; i < (stride / _v); i++)
            {
                AE_L32X2_IP(tw1, ptwd, sizeof(ae_int32x2));
                AE_L32X2_IP(tw2, ptwd, sizeof(ae_int32x2));

                AE_L32X2_XP(x0, px0, stride * sizeof(ae_int32x2));
                AE_L32X2_XP(x1, px0, stride * sizeof(ae_int32x2));
                AE_L32X2_XP(x2, px0, (_v - 2 * stride) * sizeof(ae_int32x2));

                DFT3X1_RNG_shift(x0, x1, x2, shift);
                _cmult32x32(&x1, &x1, &tw1);
                _cmult32x32(&x2, &x2, &tw2);

                AE_S32X2RNG_XP(x0, py0, _v * sizeof(ae_int32x2));
                AE_S32X2RNG_XP(x1, py0, _v * sizeof(ae_int32x2));
                AE_S32X2RNG_XP(x2, py0, _v * sizeof(ae_int32x2));
            }
        }
    }  /* if ((_v & 3) == 0) else .. */
  }

  AE_CALCRNG3();
  *bexp = 3 - RUR_AE_SAR();
  *v *= R;
  return shift;
} /* fft_stageS2_DFT3_32x32() */


/*
 *  32x32 FFT last stage Radix 5, scalingOption=2
 */
int fft_stageS2_DFT5_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  const int R = 5; // stage radix
  const int stride = N / R;
  const int shift = XT_MAX(0, 3 - *bexp);
  ae_int32x2 * restrict px0;
  ae_int32x2 * restrict py0;
  int j, _v;

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);

  _v = *v;
  NASSERT(shift>=0 && shift<4);
  
  px0 = (ae_int32x2 *)x;
  py0 = (ae_int32x2 *)y;
  /* Set scaling for DFT5X1S */
  WUR_AE_SAR(shift);
  if ((_v & 1) == 0)
  {
      __Pragma("loop_count min=1");
      for (j = 0; j < (_v >> 1); j++)
      {
          ae_int32x2 x00, x10, x20, x30, x40;
          ae_int32x2 x01, x11, x21, x31, x41;
          /* 15 cycles per pipeline stage in steady state with unroll=1 */
          AE_L32X2X2_X(x10, x11, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride);
          AE_L32X2X2_X(x20, x21, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride * 2);
          AE_L32X2X2_X(x30, x31, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride * 3);
          AE_L32X2X2_X(x40, x41, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride * 4);
          AE_L32X2X2_IP(x00, x01, castxcc(ae_int32x4, px0), 2 * sizeof(ae_int32x2));

          DFT5X1S(x00, x10, x20, x30, x40, shift);
          DFT5X1S(x01, x11, x21, x31, x41, shift);

          AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
          AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
          AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
          AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
          AE_S32X2X2RNG_XP(x40, x41, castxcc(ae_int32x4, py0), (-4 * stride + 2) * sizeof(ae_int32x2));
      }
      *v *= R;
      return shift;
  }

  /* Used in rfft N=30 */
  __Pragma("loop_count min=1");
  for (j = 0; j < _v; j++)
  {
    ae_int32x2 x0, x1, x2, x3, x4;

    x1 = AE_L32X2_X(px0, sizeof(ae_int32x2) * stride);
    x2 = AE_L32X2_X(px0, sizeof(ae_int32x2) * stride * 2);
    x3 = AE_L32X2_X(px0, sizeof(ae_int32x2) * stride * 3);
    x4 = AE_L32X2_X(px0, sizeof(ae_int32x2) * stride * 4);
    AE_L32X2_IP(x0, px0, sizeof(ae_int32x2));

    DFT5X1S(x0, x1, x2, x3, x4, shift);

    AE_S32X2RNG_XP(x0, py0, stride * sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x1, py0, stride * sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x2, py0, stride * sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x3, py0, stride * sizeof(ae_int32x2));
    AE_S32X2RNG_XP(x4, py0, (-4*stride+1) * sizeof(ae_int32x2));
  }

  return shift;
} /* fft_stageS2_DFT5_last_32x32() */


/*
*  32x32 IFFT last stage Radix 5, scalingOption=2
*/
int ifft_stageS2_DFT5_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    const int R = 5; // stage radix
    const int stride = N / R;
    const int shift = XT_MAX(0, 3 - *bexp);
    ae_int32x2 * restrict px0;
    ae_int32x2 * restrict py0;
    int j, _v;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);

    _v = *v;
    NASSERT(shift >= 0 && shift<4);

    px0 = (ae_int32x2 *)x;
    py0 = (ae_int32x2 *)y;
    /* Set scaling for DFT5X1S */
    WUR_AE_SAR(shift);
    if ((_v & 1) == 0)
    {
        __Pragma("loop_count min=1");
        for (j = 0; j < (_v >> 1); j++)
        {
            ae_int32x2 x00, x10, x20, x30, x40;
            ae_int32x2 x01, x11, x21, x31, x41;
            /* 15 cycles per pipeline stage in steady state with unroll=1 */
            AE_L32X2X2_X(x10, x11, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride);
            AE_L32X2X2_X(x20, x21, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride * 2);
            AE_L32X2X2_X(x30, x31, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride * 3);
            AE_L32X2X2_X(x40, x41, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride * 4);
            AE_L32X2X2_IP(x00, x01, castxcc(ae_int32x4, px0), 2 * sizeof(ae_int32x2));

            DFT5X1S(x00, x10, x20, x30, x40, shift);
            DFT5X1S(x01, x11, x21, x31, x41, shift);

            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x00), AE_MOVINT64_FROMINT32X2(x01), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x10), AE_MOVINT64_FROMINT32X2(x11), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x20), AE_MOVINT64_FROMINT32X2(x21), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x30), AE_MOVINT64_FROMINT32X2(x31), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
            AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x40), AE_MOVINT64_FROMINT32X2(x41), castxcc(ae_int64x2, py0), (-4 * stride + 2) * sizeof(ae_int32x2));
        }
        *v *= R;
        return shift;
    }

    /* Used in rfft N=30 */
    __Pragma("loop_count min=1");
    for (j = 0; j < _v; j++)
    {
        ae_int32x2 x0, x1, x2, x3, x4;

        x1 = AE_L32X2_X(px0, sizeof(ae_int32x2)* stride);
        x2 = AE_L32X2_X(px0, sizeof(ae_int32x2)* stride * 2);
        x3 = AE_L32X2_X(px0, sizeof(ae_int32x2)* stride * 3);
        x4 = AE_L32X2_X(px0, sizeof(ae_int32x2)* stride * 4);
        AE_L32X2_IP(x0, px0, sizeof(ae_int32x2));

        DFT5X1S(x0, x1, x2, x3, x4, shift);

        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x0), castxcc(ae_int64, py0), stride * sizeof(ae_int32x2));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x1), castxcc(ae_int64, py0), stride * sizeof(ae_int32x2));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x2), castxcc(ae_int64, py0), stride * sizeof(ae_int32x2));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x3), castxcc(ae_int64, py0), stride * sizeof(ae_int32x2));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x4), castxcc(ae_int64, py0), (-4 * stride + 1) * sizeof(ae_int32x2));
    }

    return shift;
} /* ifft_stageS2_DFT5_last_32x32() */


/*
 *  32x32 FFT/IFFT intermediate stage Radix 5, scalingOption=2
 */
int fft_stageS2_DFT5_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
  const int R = 5; // stage radix
  const int stride = N / R;
  int shift;
  int i, j, _v = *v;
        ae_int32x2 * restrict px0;
        ae_int32x2 * restrict py0;
  const ae_int32x2 * restrict ptwd;
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT((_v & 1) == 0);

  shift = XT_MAX(0, 3 - *bexp);
  NASSERT(shift>=0 && shift<5);
  /* Set scaling for AE_ADDANDSUBRNG32 */
  WUR_AE_SAR(shift);

  __Pragma("loop_count min=1");
  for (j = 0; j < _v; j += 2)
  {
      ae_int32x2 x00, x10, x20, x30, x40;
      ae_int32x2 x01, x11, x21, x31, x41;
      ae_int32x2 tw1, tw2, tw3, tw4;
      ptwd = (const ae_int32x2 *)tw;
      px0 = (ae_int32x2 *)x + j;
      py0 = (ae_int32x2 *)y + j;

      __Pragma("loop_count min=1");
      for (i = 0; i < (stride / _v); i++)
      {
          /* 18 cycles per pipeline stage in steady state with unroll=1 */
          AE_L32X2X2_IP(tw1, tw2, castxcc(ae_int32x4, ptwd), 2 * sizeof(ae_int32x2));
          AE_L32X2X2_IP(tw3, tw4, castxcc(ae_int32x4, ptwd), 2 * sizeof(ae_int32x2));

          AE_L32X2X2_X(x10, x11, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride);
          AE_L32X2X2_X(x20, x21, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride * 2);
          AE_L32X2X2_X(x30, x31, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride * 3);
          AE_L32X2X2_X(x40, x41, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride * 4);
          AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), _v*sizeof(ae_int32x2));

          DFT5X1S(x00, x10, x20, x30, x40, shift);
          DFT5X1S(x01, x11, x21, x31, x41, shift);

          x10 = AE_MULFC32RAS(x10, tw1);
          x20 = AE_MULFC32RAS(x20, tw2);
          x30 = AE_MULFC32RAS(x30, tw3);
          x40 = AE_MULFC32RAS(x40, tw4);
          x11 = AE_MULFC32RAS(x11, tw1);
          x21 = AE_MULFC32RAS(x21, tw2);
          x31 = AE_MULFC32RAS(x31, tw3);
          x41 = AE_MULFC32RAS(x41, tw4);

          AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
          AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
          AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
          AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
          AE_S32X2X2RNG_XP(x40, x41, castxcc(ae_int32x4, py0), _v * sizeof(ae_int32x2));
      }
  }

  AE_CALCRNG3();
  *bexp = 3 - RUR_AE_SAR();
  *v *= R;
  return shift;
} /* fft_stageS2_DFT5_32x32() */

/*
 *  32x32 IFFT first stage Radix 2, scalingOption=2
 */
int ifft_stageS2_DFT2_first_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{ 
    const ae_int32x2 * restrict px0;
    const ae_int32x2 * restrict px1;
    ae_int32x2 * restrict py0;
    const ae_int32x2 * restrict ptwd;
    const int stride = (N >> 1);
    int shift;
    const int R = 2; // stage radix
    const int min_shift = 2;
    int i;
    ae_int64 t0, t1, t2, t3;

    shift = min_shift - *bexp;

    NASSERT(shift>-32 && shift<4);
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT_ALIGN16(tw);
    NASSERT((stride & 1) == 1);
    NASSERT((stride >> 1) >= 3);

    //scl = 1 << shiftl;
    WUR_AE_SAR(min_shift);

    px0 = (ae_int32x2 *)x;
    px1 = px0 + stride;
    py0 = (ae_int32x2 *)y;
    ptwd = (const ae_int32x2 *)tw;

    __Pragma("loop_count min=3");
    for (i = 0; i < (stride >> 1); i++)
    {
        ae_int32x2 x00, x10, y00, y10, tw10;
        ae_int32x2 x01, x11, y01, y11, tw11;

        AE_L32X2X2_IP(tw10, tw11, castxcc(ae_int32x4, ptwd), 2 * sizeof(ae_int32x2));
        AE_L64X2_IP(t0, t1, castxcc(ae_int64x2, px0), 2 * sizeof(ae_int32x2));
        /* px1 always unaligned by 16 bytes */
        AE_L64_IP(t2, castxcc(ae_int64, px1), sizeof(ae_int32x2));
        AE_L64_IP(t3, castxcc(ae_int64, px1), sizeof(ae_int32x2));

        x00 = AE_MOVINT32X2_FROMINT64(t0);
        x01 = AE_MOVINT32X2_FROMINT64(t1);
        x10 = AE_MOVINT32X2_FROMINT64(t2);
        x11 = AE_MOVINT32X2_FROMINT64(t3);

        x00 = AE_SLAA32(x00, *bexp);
        x10 = AE_SLAA32(x10, *bexp);
        x01 = AE_SLAA32(x01, *bexp);
        x11 = AE_SLAA32(x11, *bexp);

        AE_ADDANDSUBRNG32(y00, y10, x00, x10);
        AE_ADDANDSUBRNG32(y01, y11, x01, x11);

        y11 = AE_MULFC32RAS(y11, tw11);
        y10 = AE_MULFC32RAS(y10, tw10);

        AE_S32X2X2RNG_IP(y00, y10, castxcc(ae_int32x4, py0), 2 * sizeof(ae_int32x2));
        AE_S32X2X2RNG_IP(y01, y11, castxcc(ae_int32x4, py0), 2 * sizeof(ae_int32x2));
    }

    {
        ae_int32x2 x0, x1, y0, y1;
        ae_int32x2 tw1;

        AE_L32X2_IP(tw1, ptwd, sizeof(ae_int32x2));

        AE_L64_IP(t0, castxcc(ae_int64, px0), sizeof(ae_int32x2));
        AE_L64_IP(t1, castxcc(ae_int64, px1), sizeof(ae_int32x2));
        x0 = AE_MOVINT32X2_FROMINT64(t0);
        x1 = AE_MOVINT32X2_FROMINT64(t1);

        x0 = AE_SLAA32(x0, *bexp);
        x1 = AE_SLAA32(x1, *bexp);

        AE_ADDANDSUBRNG32(y0, y1, x0, x1);
        y1 = AE_MULFC32RAS(y1, tw1);

        AE_S32X2X2RNG_IP(y0, y1, castxcc(ae_int32x4, py0), 2 * sizeof(ae_int32x2));

    }

    AE_CALCRNG3();
    *bexp = 3 - RUR_AE_SAR();
    *v *= R;
    return shift;
} /* ifft_stageS2_DFT2_first_32x32() */

/*
 *  32x32 IFFT first stage Radix 3, scalingOption=2
 */
int ifft_stageS2_DFT3_first_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{ 
    const int R = 3; // stage radix
    const int stride = N / R;
    int shift;
    ae_int32x2 * restrict px0;
    ae_int32x2 * restrict px1;
    ae_int32x2 * restrict px2;
    ae_int32x2 * restrict py0;
    const ae_int32x2 * restrict ptwd;
    int min_shift = 3;
    int i;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);

    px0 = (ae_int32x2 *)x;
    px1 = px0 + stride;
    px2 = px1 + stride;
    py0 = (ae_int32x2 *)y;
    ptwd = (const ae_int32x2 *)tw;

    shift = min_shift - *bexp;
    NASSERT(shift>-32 && shift<32);
    WUR_AE_SAR(0);

    __Pragma("loop_count min=2");
    for (i = 0; i <stride; i++)
    {
        ae_int32x2 x0, x1, x2;
        ae_int64 t0, t1, t2;
        ae_int32x2 tw1, tw2;

        AE_L32X2_IP(tw1, ptwd, sizeof(ae_int32x2));
        AE_L32X2_IP(tw2, ptwd, sizeof(ae_int32x2));

        AE_L64_IP(t0, castxcc(ae_int64, px0), sizeof(ae_int32x2));
        AE_L64_IP(t1, castxcc(ae_int64, px1), sizeof(ae_int32x2));
        AE_L64_IP(t2, castxcc(ae_int64, px2), sizeof(ae_int32x2));

        x0 = AE_MOVINT32X2_FROMINT64(t0);
        x1 = AE_MOVINT32X2_FROMINT64(t1);
        x2 = AE_MOVINT32X2_FROMINT64(t2);


        x0 = AE_SRAA32RS(x0, shift);
        x1 = AE_SRAA32RS(x1, shift);
        x2 = AE_SRAA32RS(x2, shift);

        DFT3X1(x0, x1, x2);
        _cmult32x32(&x1, &x1, &tw1);
        _cmult32x32(&x2, &x2, &tw2);

        AE_S32X2RNG_IP(x0, py0, sizeof(ae_int32x2));
        AE_S32X2RNG_IP(x1, py0, sizeof(ae_int32x2));
        AE_S32X2RNG_IP(x2, py0, sizeof(ae_int32x2));
    }

    AE_CALCRNG3();
    *bexp = 3 - RUR_AE_SAR();
    *v *= R;
    return shift;
} /* ifft_stageS2_DFT3_first_32x32() */

int fft_stageS2_DFT6_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    const int shift =  XT_MAX(0, 3 - *bexp);
    ae_int32x2 * restrict px0;
    ae_int32x2 * restrict py0;
    const int R = 6; // stage radix
    const int stride = N / R;
    int j, _v;
    _v = *v;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT((_v & 1) == 0);

    px0 = (ae_int32x2 *)x;
    py0 = (ae_int32x2 *)y;

    WUR_AE_SAR(shift);

    __Pragma("loop_count min=1");
    for (j = 0; j < (_v >> 1); j++)
    {
        /*  ?? cycles per pipeline stage in steady state with unroll=2 */
        ae_int32x2 x00, x10, x20, x30, x40, x50;
        ae_int32x2 x01, x11, x21, x31, x41, x51;

        AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x40, x41, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x50, x51, castxcc(ae_int32x4, px0), (2 - 5 * stride) * sizeof(ae_int32x2));
        {
            /* DFT6. Prime-factor FFT algorithm */
            ae_int32x2 _s0, _s1, _s2, _d0, _d1, _d2;

            _s0 = x00;    _d0 = x30;    _s1 = x40;
            _s2 = x20;    _d2 = x50;    _d1 = x10;
            /* Stage1: DFT2 with scaling */
            AE_ADDANDSUBRNG32(_s0, _d0, _s0, _d0);
            AE_ADDANDSUBRNG32(_s1, _d1, _s1, _d1);
            AE_ADDANDSUBRNG32(_s2, _d2, _s2, _d2);
            /* Stage2: DFT3 without scaling */
            DFT3X1(_s0, _s1, _s2);
            DFT3X1(_d0, _d1, _d2);

            x00 = _s0;    x10 = _d2;    x20 = _s1;
            x40 = _s2;    x50 = _d1;    x30 = _d0;
        }
        {
            /* DFT6. Prime-factor FFT algorithm */
            ae_int32x2 _s0, _s1, _s2, _d0, _d1, _d2;

            _s0 = x01;    _d0 = x31;    _s1 = x41;
            _s2 = x21;    _d2 = x51;    _d1 = x11;
            /* Stage1: DFT2 with scaling */
            AE_ADDANDSUBRNG32(_s0, _d0, _s0, _d0);
            AE_ADDANDSUBRNG32(_s1, _d1, _s1, _d1);
            AE_ADDANDSUBRNG32(_s2, _d2, _s2, _d2);
            /* Stage2: DFT3 without scaling */
            DFT3X1(_s0, _s1, _s2);
            DFT3X1(_d0, _d1, _d2);

            x01 = _s0;    x11 = _d2;    x21 = _s1;
            x41 = _s2;    x51 = _d1;    x31 = _d0;
        }

        AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x40, x41, castxcc(ae_int32x4, py0), stride * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x50, x51, castxcc(ae_int32x4, py0), (2 - 5 * stride) * sizeof(ae_int32x2));
    }

    *v = *v * R;
    return shift;
} /* fft_stageS3_DFT6_last_32x32() */

int ifft_stageS2_DFT6_last_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    const int shift = XT_MAX(0, 3 - *bexp);
    ae_int32x2 * restrict px0;
    ae_int32x2 * restrict py0;
    const int R = 6; // stage radix
    const int stride = N / R;
    int j, _v;
    _v = *v;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT((_v & 1) == 0);

    px0 = (ae_int32x2 *)x;
    py0 = (ae_int32x2 *)y;

    WUR_AE_SAR(shift);

    __Pragma("loop_count min=1");
    for (j = 0; j < (_v >> 1); j++)
    {
        /*  ?? cycles per pipeline stage in steady state with unroll=2 */
        ae_int32x2 x00, x10, x20, x30, x40, x50;
        ae_int32x2 x01, x11, x21, x31, x41, x51;

        AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x40, x41, castxcc(ae_int32x4, px0), stride * sizeof(ae_int32x2));
        AE_L32X2X2_XP(x50, x51, castxcc(ae_int32x4, px0), (2 - 5 * stride) * sizeof(ae_int32x2));
        {
            /* DFT6. Prime-factor FFT algorithm */
            ae_int32x2 _s0, _s1, _s2, _d0, _d1, _d2;

            _s0 = x00;    _d0 = x30;    _s1 = x40;
            _s2 = x20;    _d2 = x50;    _d1 = x10;
            /* Stage1: DFT2 with scaling */
            AE_ADDANDSUBRNG32(_s0, _d0, _s0, _d0);
            AE_ADDANDSUBRNG32(_s1, _d1, _s1, _d1);
            AE_ADDANDSUBRNG32(_s2, _d2, _s2, _d2);
            /* Stage2: DFT3 without scaling */
            DFT3X1(_s0, _s1, _s2);
            DFT3X1(_d0, _d1, _d2);

            x00 = _s0;    x10 = _d2;    x20 = _s1;
            x40 = _s2;    x50 = _d1;    x30 = _d0;
        }
        {
            /* DFT6. Prime-factor FFT algorithm */
            ae_int32x2 _s0, _s1, _s2, _d0, _d1, _d2;

            _s0 = x01;    _d0 = x31;    _s1 = x41;
            _s2 = x21;    _d2 = x51;    _d1 = x11;
            /* Stage1: DFT2 with scaling */
            AE_ADDANDSUBRNG32(_s0, _d0, _s0, _d0);
            AE_ADDANDSUBRNG32(_s1, _d1, _s1, _d1);
            AE_ADDANDSUBRNG32(_s2, _d2, _s2, _d2);
            /* Stage2: DFT3 without scaling */
            DFT3X1(_s0, _s1, _s2);
            DFT3X1(_d0, _d1, _d2);

            x01 = _s0;    x11 = _d2;    x21 = _s1;
            x41 = _s2;    x51 = _d1;    x31 = _d0;
        }

        AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x00), AE_MOVINT64_FROMINT32X2(x01), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
        AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x10), AE_MOVINT64_FROMINT32X2(x11), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
        AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x20), AE_MOVINT64_FROMINT32X2(x21), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
        AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x30), AE_MOVINT64_FROMINT32X2(x31), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
        AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x40), AE_MOVINT64_FROMINT32X2(x41), castxcc(ae_int64x2, py0), stride * sizeof(ae_int32x2));
        AE_S64X2_XP(AE_MOVINT64_FROMINT32X2(x50), AE_MOVINT64_FROMINT32X2(x51), castxcc(ae_int64x2, py0), (2 - 5 * stride) * sizeof(ae_int32x2));
    }

    *v = *v * R;
    return shift;
} /* ifft_stageS3_DFT6_last_32x32() */

int ifft_stageS2_DFT2_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{ return 0; }
int ifft_stageS2_DFT3_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp){ return 0; }
int ifft_stageS2_DFT5_32x32(const int32_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp){ return 0; }
