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
  NatureDSP Signal Processing Library. Mathematics
    Vector operations
    Optimized code for HiFi5
  IntegrIT, 2006-2019
*/
/* DSP Library API */
#include "NatureDSP_Signal_math.h"
#include "common.h"

/*===========================================================================
  Scalar matematics:
  scl_recip            Reciprocal on Q31/Q15 Numbers
===========================================================================*/

/*-------------------------------------------------------------------------
  Reciprocal on Q63/Q31/Q15 Numbers
  These routines return the fractional and exponential portion of the 
  reciprocal of a vector x of Q31 or Q15 numbers. Since the reciprocal is 
  always greater than 1, it returns fractional portion frac in Q(31-exp) 
  or Q(15-exp) format and exponent exp so true reciprocal value in the 
  Q0.31/Q0.15 may be found by shifting fractional part left by exponent 
  value.

  Mantissa accuracy is 1 LSB, so relative accuracy is:
  vec_recip16x16, scl_recip16x16                   6.2e-5 
  scl_recip32x32                                   2.4e-7 
  vec_recip32x32                                   9.2e-10
  vec_recip64x64                                   2.2e-19

  Precision: 
  64x64  64-bit input, 64-bit output. 
  32x32  32-bit input, 32-bit output. 
  16x16  16-bit input, 16-bit output. 

  Input:
  x[N]    input data, Q63, Q31 or Q15
  N       length of vectors

  Output:
  frac[N] fractional part of result, Q(63-exp), Q(31-exp) or Q(15-exp)
  exp[N]  exponent of result 

  Restriction:
  x,frac,exp should not overlap

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  frac,x - aligned on 16-byte boundary
  N      - multiple of 4 and >4

  Scalar versions:
  ----------------
  Return packed value: 
  scl_recip64x64():
  bits 55:0 fractional part
  bits 63:56 exponent
  scl_recip32x32():
  bits 23:0 fractional part
  bits 31:24 exponent
  scl_recip16x16():
  bits 15:0 fractional part
  bits 31:16 exponent
-------------------------------------------------------------------------*/
uint32_t scl_recip16x16 (int16_t x)
{
    ae_int16x4 x0,mx0,e0,y0;
    ae_int16x4 yy;
    xtbool4 sx0;
    x0=AE_MOVDA16(x);
    e0=AE_NSA16X4(x0);
    yy=AE_ADD16(e0,1);
    x0=AE_SRAV16RS(x0,AE_NEG16S(e0));
    AE_MOVT16X4(x0,0x4000,AE_EQ16(x0,0));
    sx0=AE_LT16(x0,0);
    x0=AE_ABS16S(x0);
    y0=AE_SUB16((int16_t)47852,x0);
    mx0=AE_NEG16S(x0);
    e0=AE_MULFD16X16X4RAS(AE_MOVDA16(0x7fff),mx0,AE_MOVDA16(0x4000),y0);
    e0=AE_ADD16(e0,e0);
    y0=AE_MULFD16X16X4RAS(y0,e0,AE_MOVDA16(0x7fff),y0);

    e0=AE_MULFD16X16X4RAS(AE_MOVDA16(0x7fff),mx0,AE_MOVDA16(0x4000),y0);
    e0=AE_ADD16(e0,e0);
    y0=AE_MULFD16X16X4RAS(y0,e0,AE_MOVDA16(0x7fff),y0);

    e0=AE_SUB16(0x4000,AE_MULFP16X4S(x0,y0));
    e0=AE_ADD16(e0,e0);
    y0=AE_MULFD16X16X4RAS(y0,e0,AE_MOVDA16(0x7fff),y0);
    AE_MOVT16X4(y0,AE_NEG16S(y0),sx0);
    yy=AE_SEL16_4321(yy,y0);
    return AE_MOVAD32_H(AE_MOVINT32X2_FROMINT16X4(yy));
} /* scl_recip16x16() */

