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
#include "sine_table32.h"
#include "common.h"
#include "NatureDSP_Signal_math.h"

/*
  NatureDSP Signal Processing Library. Scalar Mathematics
    Sine
    Optimized code for HiFi5
  IntegrIT, 2006-2019
*/

/*-------------------------------------------------------------------------
  Sine/Cosine 
  Fixed-point functions calculate sin(pi*x) or cos(pi*x) for numbers written 
  in Q31 or Q15 format. Return results in the same format. 
  Floating point functions compute sin(x) or cos(x)
  Two versions of functions available: regular version (vec_sine32x32, 
  vec_cosine32x32, , vec_sinef, vec_cosinef) 
  with arbitrary arguments and faster version (vec_sine32x32_fast, 
  vec_cosine32x32_fast) that apply some restrictions.
  NOTE:
  1.  Scalar floating point functions are compatible with standard ANSI C
      routines and set errno and exception flags accordingly
  2.  Floating point functions limit the range of allowable input values:
      [-102940.0, 102940.0] Whenever the input value does not belong to this
      range, the result is set to NaN.

  Precision: 
  32x32  32-bit inputs, 32-bit output. Accuracy: 1700 (7.9e-7)
  f      floating point. Accuracy 2 ULP

  Input:
  x[N]  input data,Q31 or floating point
  N     length of vectors
  Output:
  y[N]  output data,Q31 or floating point

  Restriction:
  Regular versions (vec_sine32x32, vec_cosine32x32, vec_sinef, 
  vec_cosinef):
  x,y - should not overlap

  Faster versions (vec_sine32x32_fast, vec_cosine32x32_fast):
  x,y - should not overlap
  x,y - aligned on 16-byte boundary
  N   - multiple of 2

  Scalar versions:
  ----------------
  return result in Q31 or floating point
-------------------------------------------------------------------------*/
int32_t scl_sine32x32   (int32_t x)
{
  ae_int32x2 x0, x0_2, x0_3, t0, y0, z0;
  ae_int32x2 p0, p1, p2, p3, p4, p5;
  x0 = AE_MOVDA32X2(x, x);
  p0 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 0);
  p1 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 1);
  p2 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 2);
  p3 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 3);
  p4 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 4);
  p5 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 5);

  x0_2 = AE_MULFP32X2RAS(x0, x0);//x^2
  x0_3 = AE_MULFP32X2RAS(x0_2, x0);//x^3

  t0 = p0; y0 = p1;
  AE_MULAFP32X2RAS(y0, x0_2, t0); t0 = y0;
  y0 = p2;        
  AE_MULAFP32X2RAS(y0, x0_2, t0); t0 = y0;
  y0 = p3;        
  AE_MULAFP32X2RAS(y0, x0_2, t0); t0 = y0;
  y0 = p4;        
  AE_MULAFP32X2RAS(y0, x0_2, t0); t0 = y0;
  y0 = AE_MULFP32X2RAS(x0, p5);
  AE_MULAFP32X2RAS(y0, x0_3, t0); //Q28
  z0 = AE_SLAI32S(y0, 3);
  return  AE_MOVAD32_H(z0);
} /* scl_sine32x32() */
