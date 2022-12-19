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
    Optimized for HiFi5
	Integrit, 2006-2019
*/
#include "NatureDSP_Signal_fft.h"
#include "NatureDSP_Signal_vector.h"
#include "common.h"
/* Twiddle factor tables and FFT descriptor structure. */
#include "fft_x16_common.h"
#include "fft_16x16_stages.h"

#if 0
#include "baseop.h"
static uint16_t _HI(uint32_t x) { return (uint16_t)(x >> 16); }
static uint16_t _LO(uint32_t x) { return (uint16_t)(x); }
static uint32_t _PACK(uint16_t hi, uint16_t lo) { return (((uint32_t)hi) << 16) | (lo); }
static uint32_t _SHR2(uint32_t x, int rsh)
{
  int16_t hi = _HI(x), lo = _LO(x);
  hi >>= rsh; lo >>= rsh;
  return _PACK(hi, lo);
}
static uint32_t _SHR2R(uint32_t x, int rsh)
{
  int32_t hi = (int32_t)(int16_t)_HI(x) + (1L<<(rsh-1));
  int32_t lo = (int32_t)(int16_t)_LO(x) + (1L<<(rsh-1));
  hi >>= rsh; lo >>= rsh;
  return _PACK(S_sature_l(hi), S_sature_l(lo));
}
static uint32_t _ADD2S(uint32_t x, uint32_t y)
{
  return _PACK(S_add_ss(_HI(x), _HI(y)), S_add_ss(_LO(x), _LO(y)));
}
static uint32_t _SUB2S(uint32_t x, uint32_t y)
{
  return _PACK(S_sub_ss(_HI(x), _HI(y)), S_sub_ss(_LO(x), _LO(y)));
}
/* Two-way bidirectional arithmetic (sign-extending) right shift */
static uint32_t _SHR2R_BIDIR(uint32_t x, int rsh)
{
    int16_t hi = _HI(x), lo = _LO(x);
    if (rsh>0)
    {
        int32_t _hi = (int32_t)hi + (1L<<(rsh-1));
        int32_t _lo = (int32_t)lo + (1L<<(rsh-1));
        _hi >>= rsh; _lo >>= rsh;
        hi = S_sature_l(_hi); lo = S_sature_l(_lo);
    }
    else
    {
        hi <<= -rsh; lo <<= -rsh;
    }
    return _PACK(hi, lo);
}

static uint32_t cmult16x16(uint32_t x, uint32_t y)
{
    int32_t re = S_round_l(L_sub_ll(L_mpy_ss(_LO(x), _LO(y)), L_mpy_ss(_HI(x), _HI(y))));
    int32_t im = S_round_l(L_add_ll(L_mpy_ss(_LO(x), _HI(y)), L_mpy_ss(_HI(x), _LO(y))));
    return _PACK(im, re);
}

static void _cmult16x16(uint32_t *result, uint32_t *x, uint32_t *y)
{
    int16_t re = S_extract_l(L_sub_ll(L_mpy_ss(_LO(*x), _LO(*y)), L_mpy_ss(_HI(*x), _HI(*y))));
    int16_t im = S_extract_l(L_add_ll(L_mpy_ss(_LO(*x), _HI(*y)), L_mpy_ss(_HI(*x), _LO(*y))));
    *result = _PACK(im, re);
}
#endif
/*
 *  Intermediate stage of FFT/IFFT 16x16, radix-2, dynamic scaling
 */
int fft_16x16_stage_inner_scl2_DFT2_v4(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{ 
  int i;//, j;
  int Ninner;
  const int shift = XT_MAX(0, 2 - *bexp);
  const int stride = N>>1;
  const ae_int16x8 * restrict px0;
  const ae_int16x8 * restrict px1;
        ae_int16x8 * restrict py0;
  const ae_int16x4 * restrict ptwd;
  
  ae_int16x4 y0, y1;
  ae_int16x4 y0_, y1_;
  ae_int16x4 tw01;
  ae_valignx2 alx1;
  
  NASSERT(v[0]==4);
  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT(tw_step == 1);
  NASSERT(stride%4==0);
  Ninner=stride>>2;

  {
      /* Reset RANGE register */
      int a, b;
      AE_CALCRNG16(a, b, 0, 3);
      (void)b;
      (void)a;

  }

  WUR_AE_SAR( shift*0x102 ); 
  ptwd = (const ae_int16x4 *)tw;
  px0 = (const ae_int16x8 *)((complex_fract16 *)x);
  px1 = (const ae_int16x8 *)((complex_fract16 *)px0 + stride);
  py0 = (      ae_int16x8 *)((complex_fract16 *)y);
  alx1 = AE_LA128_PP(px1);
  __Pragma("loop_count min=1");
  for (i = 0; i < (Ninner>>1); i++)
  { /* 8 cycles per pipeline stage in steady state with unroll=1 */
      ae_int32x2 t1; 
      /* butterfly 0 */
      AE_L16X4X2_IP(y0, y0_, px0, 4*sizeof(complex_fract16));
      AE_L16X4X2_IP(y1, y1_, px1, 4*sizeof(complex_fract16));

      AE_ADDANDSUBRNG16RAS_S2(y0, y1); 
      AE_ADDANDSUBRNG16RAS_S2(y0_, y1_); 

      AE_L32_XP(t1, castxcc(ae_int32, ptwd), sizeof(complex_fract16));
      tw01 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t1));
      y1 = AE_MULFC16RAS(y1, tw01);
      y1_ = AE_MULFC16RAS(y1_, tw01);

      AE_S16X4X2RNG_IP(y0, y0_, py0, 4*sizeof(complex_fract16));
      AE_S16X4X2RNG_IP(y1, y1_, py0, 4*sizeof(complex_fract16));

      /* butterfly 1 */
      AE_L16X4X2_IP(y0, y0_, px0, 4*sizeof(complex_fract16));
      AE_L16X4X2_IP(y1, y1_, px1, 4*sizeof(complex_fract16));

      AE_ADDANDSUBRNG16RAS_S2(y0, y1); 
      AE_ADDANDSUBRNG16RAS_S2(y0_, y1_);

      AE_L32_XP(t1, castxcc(ae_int32, ptwd), sizeof(complex_fract16));
      tw01 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t1));
      y1 = AE_MULFC16RAS(y1, tw01);
      y1_ = AE_MULFC16RAS(y1_, tw01);

      AE_S16X4X2RNG_IP(y0, y0_, py0, 4*sizeof(complex_fract16));
      AE_S16X4X2RNG_IP(y1, y1_, py0, 4*sizeof(complex_fract16));
  }
  if (Ninner&1)
  {
      ae_int32x2 t1; 
      /* butterfly 0 */
      AE_L16X4X2_IP(y0, y0_, px0, 4*sizeof(complex_fract16));
      AE_L16X4X2_IP(y1, y1_, px1, 4*sizeof(complex_fract16));

      AE_ADDANDSUBRNG16RAS_S2(y0, y1); 
      AE_ADDANDSUBRNG16RAS_S2(y0_, y1_); 

      AE_L32_XP(t1, castxcc(ae_int32, ptwd), sizeof(complex_fract16));
      tw01 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t1));
      y1 = AE_MULFC16RAS(y1, tw01);
      y1_ = AE_MULFC16RAS(y1_, tw01);

      AE_S16X4X2RNG_IP(y0, y0_, py0, 4*sizeof(complex_fract16));
      AE_S16X4X2RNG_IP(y1, y1_, py0, 4*sizeof(complex_fract16));
  }

  {
      int a, b;
      AE_CALCRNG16(a, b, 0, 3);
      *bexp = 3 - a;
      (void)b;

  }

  *v *= 2;
  return shift;
} /* fft_16x16_stage_inner_scl2_DFT2() */
