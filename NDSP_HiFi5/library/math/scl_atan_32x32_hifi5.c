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
#include "atan_table.h"
#include "common.h"
#include "NatureDSP_Signal_math.h"

/*
  NatureDSP Signal Processing Library. Scalar Mathematics
   Arctangent
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/*-------------------------------------------------------------------------
  Arctangent 
  Functions calculate arctangent of number. Fixed point functions 
  scale output to pi so it is always in range -0x20000000 : 0x20000000 
  which corresponds to the real phases +pi/4 and represent input and output 
  in Q31
  NOTE:
  1.  Scalar floating point function is compatible with standard ANSI C
      routines and sets errno and exception flags accordingly

  Accuracy:
  24 bit version: 74000 (3.4e-5) 
  32 bit version: 42    (2.0e-8)
  floating point: 2 ULP

  Precision: 
  32x32  32-bit inputs, 32-bit output.
  f      floating point

  Input:
  x[N]   input data, Q31 or floating point
  N      length of vectors
  Output:
  z[N]   result, Q31 or floating point

  Restriction:
  x,z should not overlap

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  x,z - aligned on 16-byte boundary
  N   - multiple of 2

  Scalar versions:
  ----------------
  return result, Q31 or floating point
-------------------------------------------------------------------------*/
int32_t scl_atan32x32(int32_t x) 
{
  ae_int32x2 x0, x0_2, x0_3, z0;
  ae_int32x2 y0, y2, y4, t0, t2, t4;
  ae_int32x2 p0, p1, p2, p3, p4, p5;
  ae_int32x2 p6, p7, p8, p9, p10;
  xtbool2 sx0;
  x0 = AE_MOVDA32X2(x, x);
  sx0 = AE_LT32(x0, AE_ZERO32());

  p0  = AE_L32_I((const ae_int32 *)atan_table, 4 * 0);
  p1  = AE_L32_I((const ae_int32 *)atan_table, 4 * 1);
  p2  = AE_L32_I((const ae_int32 *)atan_table, 4 * 2);
  p3  = AE_L32_I((const ae_int32 *)atan_table, 4 * 3);
  p4  = AE_L32_I((const ae_int32 *)atan_table, 4 * 4);
  p5  = AE_L32_I((const ae_int32 *)atan_table, 4 * 5);
  p6  = AE_L32_I((const ae_int32 *)atan_table, 4 * 6);
  p7  = AE_L32_I((const ae_int32 *)atan_table, 4 * 7);
  p8  = AE_L32_X((const ae_int32 *)atan_table, 4 * 8);
  p9  = AE_L32_X((const ae_int32 *)atan_table, 4 * 9);
  p10 = AE_L32_X((const ae_int32 *)atan_table, 4 * 10);

  x0 = AE_ABS32S(x0);
  x0 = AE_SUB32(x0, 0x40000000);

  x0_2 = AE_MULFP32X2RAS(x0, x0);//x^2
  x0_3 = AE_MULFP32X2RAS(x0_2, x0);//x^3

  y0 = p3;
  AE_MULAFP32X2RAS(y0, x0_3, p0); t0 = y0;
  y2 = p4;
  AE_MULAFP32X2RAS(y2, x0_3, p1); t2 = y2;
  y4 = p5;
  AE_MULAFP32X2RAS(y4, x0_3, p2); t4 = y4;

  y0 = p6;
  AE_MULAFP32X2RAS(y0, x0_3, t0); t0 = y0;
  y2 = p7;
  AE_MULAFP32X2RAS(y2, x0_3, t2); t2 = y2;
  y4 = p8;
  AE_MULAFP32X2RAS(y4, x0_3, t4); t4 = y4;

  y0 = p9;
  AE_MULAFP32X2RAS(y0, x0_3, t0); t0 = y0;
  y2 = p10;
  AE_MULAFP32X2RAS(y2, x0_3, t2); t2 = y2;

  AE_MULAFP32X2RAS(y2, x0_2, y4);

  AE_MULAFP32X2RAS(y2, x0, y0); y0 = y2;

  z0 = AE_NEG32S(y0);
  AE_MOVT32X2(y0, z0, sx0);
  return  AE_MOVAD32_H(y0);

} /* scl_atan32x32() */
