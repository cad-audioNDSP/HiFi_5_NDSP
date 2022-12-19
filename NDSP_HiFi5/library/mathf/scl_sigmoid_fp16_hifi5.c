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
  NatureDSP Signal Processing Library. Scalar matematics
    Sigmoid, 
    Code optimized for HiFi5 core
  IntegrIT, 2006-2019
*/
#include "NatureDSP_Signal_math.h"
#include "NatureDSP_types.h"
#include "common.h"
#include "common_fpu.h"


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
#if !HAVE_HPFPU
DISCARD_FUN(float16_t,scl_sigmoid_fp16,(float16_t x))
#else
float16_t scl_sigmoid_fp16(float16_t x)
{
  static const union ufloat16uint16 ALIGN(32) cnt[] = { { 0xcc29 }, /* -16.6355==-log(1/double(eps(it_half(0)))-1) */
  { 0x3dc5 },   /* 1.4424 */
  { 0x0d1e }};   /* 0.00031233 */
  
  /* polynonial coefficients for 2^x, x=-0.5...0.5 */
  static const union ufloat16uint16 ALIGN(32) p[] = { { 0x2b27 }, { 0x33c1 }, { 0x398c }, { 0x3c00 } };
  const xthalf * restrict pP = (const xthalf  *)p;
  const xthalf  * restrict pC = (const xthalf   *)cnt;

  int16_t n0, n1, nn;
  ae_int16x4 y0, y1;
 
  xthalf y, z, d, t;
  int sx = x;
  xthalf minsigmoid_fp16 = AE_LHI(pC, 0);
  xthalf log2e0 = AE_LHI(pC, 2);
  xthalf log2e1 = AE_LHI(pC, 4);
  
  xthalf fx = AE_LHI((xthalf *)(&x), 0);
  {
    xtbool4 b0;
    xthalf x = fx;
    xthalf t0, t1;
    x = NEG_H(ABS_H(x));
    x = MAX_H(minsigmoid_fp16, x);
    /* compute d+n=log2(e)*x */
    y = FIROUND_H(MUL_H(x, log2e0));
    d = y; MSUBN_H(d, x, log2e0);
    MSUBN_H(d, x, log2e1);
    nn = TRUNC16_H(y, 0);
    /* approx 2^d */
    z = AE_LHI(pP, 0);
    t = AE_LHI(pP, 2); MSUBN_H(t, d, z); z = t;
    t = AE_LHI(pP, 4); MSUBN_H(t, d, z); z = t;
    t = CONST_H(1); MSUBN_H(t, d, z); z = t;
    /* simplified ldexpf */
    n0 = AE_SRAI16(nn, 1);
    n1 = AE_SUB16(nn, n0);
    n0 = AE_SLLI16S(AE_ADD16(n0, 15), 10);
    n1 = AE_SLLI16S(AE_ADD16(n1, 15), 10);
    t0 = AE_LHI((xthalf *)(&n1), 0);
    t1 = AE_LHI((xthalf *)(&n0), 0);
    x = MUL_H(t0, MUL_H(t1,z));
    /* approx y=1/(1+x); */
    t = ADD_H(CONST_H(1), x);
    y = RECIP0_H(t);
    t = SUB_H(CONST_H(1), y);
    MSUBN_H(t, x, y);
    MADDN_H(y, t, y);

    t = MUL_H(y, x);
    b0 = AE_LT16(sx, AE_ZERO16());
    y0 = AE_L16_I((ae_int16 *)(&y), 0);
    y1 = AE_L16_I((ae_int16 *)(&t), 0);
    AE_MOVT16X4(y0, y1, b0);
  }
  return AE_MOVAD16_0(y0);

} /* scl_sigmoid_fp16() */
#endif
