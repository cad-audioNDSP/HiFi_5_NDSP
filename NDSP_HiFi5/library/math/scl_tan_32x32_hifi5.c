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
#include "common.h"
#include "NatureDSP_Signal_math.h"
#include "vec_tan32x32_table.h"
/*
  NatureDSP Signal Processing Library. Scalar Mathematics
   Tangent
    Optimized code for HiFi5
  IntegrIT, 2006-2019
*/
/*-------------------------------------------------------------------------
  Tangent 
  Fixed point functions calculate tan(pi*x) for number written in Q31. 
  Floating point functions compute tan(x)
  
  Precision: 
  32x32  32-bit inputs, 32-bit outputs. Accuracy: (1.3e-4*y+1LSB)
                                        if abs(y)<=464873(14.19 in Q15) 
                                        or abs(x)<pi*0.4776
  f      floating point.                Accuracy: 2 ULP

  NOTE:
  1.  Scalar floating point function is compatible with standard ANSI C routines 
      and set errno and exception flags accordingly
  2.  Floating point functions limit the range of allowable input values: [-9099, 9099]
      Whenever the input value does not belong to this range, the result is set to NaN.

  Input:
  x[N]   input data,Q31 or floating point
  N      length of vectors
  Output:
  y[N]   result, Q16.15 or floating point

  Restriction:
  x,y - should not overlap

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  x,z - aligned on 16-byte boundary
  N   - multiple of 2

  Scalar versions:
  ----------------
  Return result, Q16.15 or floating point
-------------------------------------------------------------------------*/

int32_t scl_tan32x32 (int32_t x)
{
  ae_int32x2 x0, x2_0, r0;
  ae_int32x2 n0, j0, t0;
  ae_int32x2 yt0, yc0, sg0;
  xtbool2 sign0, st0;
  xtbool2 sx0;
  ae_int16x4 nsa;
  x0 = x;
  sg0 = x0;
  x0 = AE_ABS32S(x0);
  /*
  * argument reduction
  */
  n0 = AE_SRAI32(x0, 29);
  n0 = AE_ADD32(n0, 1);
  j0 = AE_AND32(n0, -2);
  AE_MULSP32X2(x0, j0, 0x20000000);
  x2_0 = AE_MULFP32X2RAS(x0, x0);
  /*Tangent*/
  yt0 = AE_L32_I((const ae_int32 *)Pt, 4 * 0);
  t0 = AE_L32_I((const ae_int32 *)Pt, 4 * 1);
  AE_MULAFP32X2RAS(t0, x2_0, yt0); yt0 = t0; 
  t0 = AE_L32_I((const ae_int32 *)Pt, 4 * 2);
  AE_MULAFP32X2RAS(t0, x2_0, yt0); yt0 = t0; 
  t0 = AE_L32_I((const ae_int32 *)Pt, 4 * 3);
  AE_MULAFP32X2RAS(t0, x2_0, yt0); yt0 = t0; 
  t0 = AE_L32_I((const ae_int32 *)Pt, 4 * 4);
  AE_MULAFP32X2RAS(t0, x2_0, yt0); yt0 = t0; 
  t0 = AE_SRAI32(x0, 5);
  yt0 = AE_MULFP32X2RAS(t0, yt0); //Q16.15
  /*Cotangent*/
  yc0 = AE_L32_I((const ae_int32 *)Qt, 4 * 0);
  t0 = AE_L32_I((const ae_int32 *)Qt, 4 * 1);
  AE_MULAFP32X2RAS(t0, x2_0, yc0); yc0 = t0;
  t0 = AE_L32_I((const ae_int32 *)Qt, 4 * 2);
  AE_MULAFP32X2RAS(t0, x2_0, yc0); yc0 = t0;
  t0 = AE_L32_I((const ae_int32 *)Qt, 4 * 3);
  AE_MULAFP32X2RAS(t0, x2_0, yc0); yc0 = t0;
  t0 = AE_L32_I((const ae_int32 *)Qt, 4 * 4);
  AE_MULAFP32X2RAS(t0, x2_0, yc0); yc0 = t0;
  t0 = AE_L32_I((const ae_int32 *)Qt, 4 * 5);
  AE_MULAFP32X2RAS(t0, x2_0, yc0); yc0 = t0;

  /* 1/x */
  {
    ae_int32x2 y0, _0x40000000;
    ae_int32x2 e0, xx0;

    _0x40000000 = AE_MOVDA32(0x40000000);
    nsa = AE_NSA32X4(x0, x0);
    xx0 = AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(nsa)));
    sx0 = AE_LT32(xx0, AE_ZERO32());
    xx0 = AE_INT32X2_ABS32S(xx0);
   
    y0 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), xx0);
   
    t0 = _0x40000000;
    AE_MULSFP32X2RAS(t0, xx0, y0); e0 = t0;
    e0 = AE_ADD32(e0, e0);
    t0 = y0;
    AE_MULAFP32X2RAS(t0, e0, y0); y0 = t0;
    t0 = _0x40000000;
    AE_MULSFP32X2RAS(t0, xx0, y0); e0 = t0;
    e0 = AE_ADD32(e0, e0);
    t0 = y0;
    AE_MULAFP32X2RAS(t0, e0, y0); r0 = t0;
  }

  
  t0 = AE_INT32X2_NEG32S(yc0);
  AE_MOVT32X2(yc0, t0, sx0);
  yc0 = AE_MULFP32X2RAS(r0, yc0); //Q(24-nsa-1) 

  nsa = AE_SUB16(8, nsa);
  yc0 = AE_SRAV32RS(yc0, AE_SEXT32X2D16_32(nsa));

  /* adjust sign */
  n0 = AE_SRAI32(n0, 1);
  n0 = AE_SLLI32(n0, 31);
  
  sg0 = AE_AND32(sg0, 0x80000000);
  sg0 = AE_XOR32(sg0, n0);
  st0 = AE_EQ32(n0, 0x80000000);
  sign0 = AE_EQ32(sg0, 0x80000000);

  AE_MOVT32X2(yt0, yc0, st0);

  t0 = AE_NEG32S(yt0);
  AE_MOVT32X2(yt0, t0, sign0);
  return AE_MOVAD32_H(yt0);
} /* scl_tan32x32() */
