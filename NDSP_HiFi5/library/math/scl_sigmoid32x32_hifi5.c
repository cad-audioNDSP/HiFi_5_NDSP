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
  NatureDSP Signal Processing Library. Scalar Mathematics
    Sigmoid
    Optimized code for HiFi5
  IntegrIT, 2006-2019
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Signal_math.h"
#include "common.h"

/*-------------------------------------------------------------------------
  Sigmoid
  The functions compute the sigmoid of input argument. 32-bit fixed-point 
  functions accept inputs in Q6.25 and form outputs in Q16.15 format.

  Precision:
  32x32  32-bit inputs, 32-bit output. Accuracy: 2 LSB.
  f      single precision floating-point. Accuracy 2 ULP
  fp16   half precision floating-point. Accuracy 2 ULP
  Input:
  x[N]   input data, Q6.25 or floating point
  N      length of vectors
  Output:
  y[N]   result, Q16.15 or floating point

  Restriction:
  x,y should not overlap

  Scalar versions:
  ----------------
  return result, Q16.15 or floating point
-------------------------------------------------------------------------*/
static const int32_t ALIGN(32) polypow2[] = { 14685057, -114217091, 514075394, -1488269031, 2147475316 };

int32_t scl_sigmoid32x32(int32_t x)
{
  ae_int32x2 x0, e0, z0, d0, y0, t0;
  xtbool2 sign0;
  x0 = AE_MOVDA32X2(x, x);
  sign0 = AE_LT32(x0, 0);

  z0 = AE_MULFP32X2RAS(x0, 774541002);
  x0 = AE_ABS32S(z0);
  e0 = AE_SRAI32(x0, 23);

  x0 = AE_MOVDEXT(x0, 23,8);//Q31
  y0 = AE_L32_I((const ae_int32 *)polypow2, 4 * 0);
  t0 = AE_L32_I((const ae_int32 *)polypow2, 4 * 1);
  AE_MULAFP32X2RAS(t0, x0, y0); y0 = t0;
  t0 = AE_L32_I((const ae_int32 *)polypow2, 4 * 2);
  AE_MULAFP32X2RAS(t0, x0, y0); y0 = t0;
  t0 = AE_L32_I((const ae_int32 *)polypow2, 4 * 3);
  AE_MULAFP32X2RAS(t0, x0, y0); y0 = t0;
  t0 = AE_L32_I((const ae_int32 *)polypow2, 4 * 4);
  AE_MULAFP32X2RAS(t0, x0, y0); y0 = t0;
  x0 = AE_SRAV32RS(y0, e0);

  /* 0.96-x/2 */
  z0 = 1030792151;
  AE_MULSFP32X2RAS(z0, x0, 0x20000000);//Q30
  t0 = z0;
  AE_MULAFP32X2RAS(t0, z0, x0);
  d0 = AE_SUB32(1073741824, t0);//Q30
  t0 = AE_SRAI32(z0, 1);
  AE_MULAFP32X2RAS(t0, z0, d0); z0 = t0;
  AE_MULAFP32X2RAS(t0, z0, x0);
  d0 = AE_SUB32(536870912, t0);//Q29
  t0 = AE_MULFP32X2RAS(z0,0x20000000);
  AE_MULAFP32X2RAS(t0, z0, d0); z0 = t0;

  z0 = AE_SRAA32RS(z0, 12);
  x0 = AE_SUB32(32768, z0);

  AE_MOVT32X2(z0, x0, sign0);

  return AE_MOVAD32_H(z0);
} /* scl_sigmoid32x32() */
