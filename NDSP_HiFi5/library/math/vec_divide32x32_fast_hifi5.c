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
#include "NatureDSP_Signal_math.h"
#include "common.h"

/*
  NatureDSP Signal Processing Library. Vector Mathematics
   Division
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/


/*-------------------------------------------------------------------------
  Division
  These routines perform pair-wise division of vectors written in Q63, Q31 or 
  Q15 format. They return the fractional and exponential portion of the division 
  result. Since the division may generate result greater than 1, it returns 
  fractional portion frac in Q(63-exp), Q(31-exp) or Q(15-exp) format and 
  exponent exp so true division result in the Q0.31 may be found by shifting 
  fractional part left by exponent value.
  Additional routine makes integer division of 64-bit number to 32-bit 
  denominator forming 32-bit result. If result is overflown, 0x7fffffff 
  or 0x80000000 is returned depending on the signs of inputs.
  For division to 0, the result is not defined.

  Two versions of routines are available: regular versions (vec_divide64x32i,
  vec_divide64x64, vec_divide32x32, vec_divide16x16) work 
  with arbitrary arguments, faster versions (vec_divide32x32_fast, 
  vec_divide16x16_fast) apply some restrictions.

  Accuracy is measured as accuracy of fractional part (mantissa):
  vec_divide64x32i, scl_divide64x32                      :  1 LSB   
  vec_divide64x64                                        :  2 LSB 
  vec_divide32x32, vec_divide32x32_fast                  :  2 LSB (1.8e-9) 
  scl_divide32x32                                        :  2 LSB (4.8e-7) 
  vec_divide16x16, scl_divide16x16, vec_divide16x16_fast :  2 LSB (1.2e-4)

  Precision: 
  64x32i integer division, 64-bit nominator, 32-bit denominator, 32-bit output. 
  64x64  fractional division, 64-bit inputs, 64-bit output. 
  32x32  fractional division, 32-bit inputs, 32-bit output. 
  16x16  fractional division, 16-bit inputs, 16-bit output. 

  Input:
  x[N]    nominator, 64-bit integer, Q63, Q31 or Q15
  y[N]    denominator, 32-bit integer, Q63, Q31 or Q15
  N       length of vectors
  Output:
  frac[N] fractional parts of result, Q(63-exp), Q(31-exp) or Q(15-exp)
  exp[N]  exponents of result 

  Restriction:
  For regular versions (vec_divide64x32i, vec_divide64x64, vec_divide32x32,
  vec_divide16x16) :
  x,y,frac,exp should not overlap

  For faster versions (vec_divide32x32_fast, vec_divide16x16_fast) :
  x,y,frac,exp  should not overlap
  x,y,frac      to be aligned by 16-byte boundary, N - multiple of 4.

  Scalar versions:
  ----------------
  scl_divide64x32(): integer remainder
  Return packed value: 
  scl_divide64x64():
  bits 55:0 fractional part
  bits 63:56 exponent
  scl_divide32x32():
  bits 23:0 fractional part
  bits 31:24 exponent
  scl_divide16x16():
  bits 15:0 fractional part
  bits 31:16 exponent
-------------------------------------------------------------------------*/

void vec_divide32x32_fast 
(
  int32_t * restrict        frac,
  int16_t * restrict        exp,
  const int32_t * restrict  x,
  const int32_t * restrict  y,
  int                       N
)
{
   int n;
        ae_int32x4 * restrict pF = (      ae_int32x4 *)frac;
        ae_int16x4 * restrict pE = (      ae_int16x4 *)exp;
  const ae_int32x4 * restrict pX = (const ae_int32x4 *)x;
  const ae_int32x4 * restrict pY = (const ae_int32x4 *)y;
  ae_valignx2 aX, aY, aF;
  ae_valign aE;
  aX = AE_LA128_PP(pX);
  aY = AE_LA128_PP(pY);
  aE = AE_ZALIGN64();
  aF = AE_ZALIGN128(); 

  if (N <= 0) return;
  NASSERT_ALIGN8(x);
  NASSERT_ALIGN8(y);
  NASSERT_ALIGN8(frac);
  NASSERT(N % 4 == 0);
  pX = (const ae_int32x4 *)x;
  pF = (ae_int32x4*)frac;
  __Pragma("loop_count min=1")
  for (n = 0; n<(N >> 2); n++)
  {
    ae_int32x2 sx0, sx1, x0, x1, t0, t1, y0, y1, e0, e1, z0, z1;
    ae_int32x2 _0x40000000;
    ae_int16x4 nsax, nsay, expxy;
    _0x40000000 = AE_MOVDA32(0x40000000);
    AE_LA32X2X2_IP(x0, x1, aX, pX);
    AE_LA32X2X2_IP(y0, y1, aY, pY);

    nsax = AE_NSA32X4(x0, x1);
    nsay = AE_NSA32X4(y0, y1);
    x0 = AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(nsax)));
    x1 = AE_SRAV32RS(x1, AE_SEXT32X2D16_10(AE_NEG16S(nsax)));
    y0 = AE_SRAV32RS(y0, AE_SEXT32X2D16_32(AE_NEG16S(nsay)));
    y1 = AE_SRAV32RS(y1, AE_SEXT32X2D16_10(AE_NEG16S(nsay)));

    expxy = AE_SUB16S(nsay, nsax);
    expxy = AE_INT16X4_ADD_INT16X4(expxy, 1);
    AE_SA16X4_IP(expxy, aE, pE);
    z0 = x0; z1 = x1;

    sx0 = y0;
    sx1 = y1;
    x0 = AE_INT32X2_ABS32S(y0);
    x1 = AE_INT32X2_ABS32S(y1);
    y0 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), x0);
    y1 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), x1);
    /* 4 iterations to achieve 1 LSB accuracy in mantissa */
    AE_MOVD32X4(t0,t1,_0x40000000,_0x40000000);
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    t0 = y0; t1 = y1;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
    AE_MOVD32X4(t0,t1,_0x40000000,_0x40000000);
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    t0 = y0; t1 = y1;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
    AE_MOVD32X4(t0,t1,_0x40000000,_0x40000000);
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    t0 = y0; t1 = y1;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
    AE_MOVD32X4(t0,t1,_0x40000000,_0x40000000);
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    t0 = y0; t1 = y1;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;

    /* restore original sign */
    y0=AE_MOVNEG32S_T(y0, sx0);
    y1=AE_MOVNEG32S_T(y1, sx1);

    AE_MULF2P32X4RAS(y0, y1, z0, z1, y0, y1);

    AE_SA32X2X2_IP(y0, y1, aF, pF);
  }
  AE_SA128POS_FP(aF, pF);
  AE_SA64POS_FP(aE, pE);
} /* vec_divide32x32_fast() */
