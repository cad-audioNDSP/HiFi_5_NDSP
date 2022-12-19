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
#include "vec_log_table.h"
#include "common.h"
#include "NatureDSP_Signal_math.h"

/*
  NatureDSP Signal Processing Library. Scalar Mathematics
   Logarithm, base 10
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/


/*-------------------------------------------------------------------------
  Logarithm:
  Different kinds of logarithm (base 2, natural, base 10). Fixed point 
  functions represent results in Q25 format or return 0x80000000 on negative 
  of zero input.

  Precision:
  32x32  32-bit inputs, 32-bit outputs
  f      floating point

  Accuracy :
  vec_log2_32x32,scl_log2_32x32              730 (2.2e-5)
  vec_logn_32x32,scl_logn_32x32              510 (1.5e-5)
  vec_log10_32x32,scl_log10_32x32            230 (6.9e-6)
  floating point                             2 ULP

  NOTES:
  1.  Although 32 and 24 bit functions provide the same accuracy, 32-bit 
      functions have better input/output resolution (dynamic range)
  2.  Scalar Floating point functions are compatible with standard ANSI C routines 
      and set errno and exception flags accordingly.
  3.  Floating point functions limit the range of allowable input values:
      A) If x<0, the result is set to NaN. In addition, scalar floating point
         functions assign the value EDOM to errno and raise the "invalid" 
         floating-point exception.
      B) If x==0, the result is set to minus infinity. Scalar floating  point
         functions assign the value ERANGE to errno and raise the "divide-by-zero"
         floating-point exception.

  Input:
  x[N]  input data, Q16.15 or floating point 
  N     length of vectors
  Output:
  y[N]  result, Q25 or floating point 

  Restriction:
  x,y should not overlap

  Scalar versions:
  ----------------
  return result in Q25 or floating point
-------------------------------------------------------------------------*/
int32_t scl_log10_32x32(int32_t x)
{

  ae_int32x2 x0, x2, e0, d0;
  ae_int32x2 t0, y0, t1, t2, t3, y2;
  int32_t     y;
  xtbool2     inf0, i0;
  ae_int32x2 vsw, vhw, vzw;
  ae_int16x4 nsa01;
  ae_int32x2 c0, c1, c2, c3, c4, c5;
  vsw = AE_MOVDA32X2(646456993, 646456993); /*log10(2)*/
  vhw = AE_MOVDA32X2(0x5A82799A, 0x5A82799A);
  vzw = AE_MOVI(0);
  x0 = AE_MOVDA32X2(x, x);/*Q.16.15*/
  inf0 = AE_LE32(x0, vzw);
  /* Normalize x*/
  nsa01 = AE_NSA32X4(x0, x0);
  x0 = AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(nsa01)));
  nsa01 = AE_SUB16(16, nsa01);

  i0 = AE_LE32(x0, vhw);
  e0 = AE_SEXT32X2D16_32(nsa01);
  AE_MOVT32X2(e0, AE_SUB32(e0, 1), i0);
  e0 = AE_SLAI32(e0, (31 - 6));
  e0 = AE_MULFP32X2RAS(e0, vsw);
  c0 = AE_L32_I((const ae_int32 *)log10_tbl_Q31, 4 * 0);
  c1 = AE_L32_I((const ae_int32 *)log10_tbl_Q31, 4 * 1);
  c2 = AE_L32_I((const ae_int32 *)log10_tbl_Q31, 4 * 2);
  c3 = AE_L32_I((const ae_int32 *)log10_tbl_Q31, 4 * 3);
  c4 = AE_L32_I((const ae_int32 *)log10_tbl_Q31, 4 * 4);
  c5 = AE_L32_I((const ae_int32 *)log10_tbl_Q31, 4 * 5);

  t0 = 0x80000000; d0 = 1;
  AE_MOVT32X2(d0, 2, i0);

  AE_MULSP32X2(t0, x0, d0); x0 = t0;

  x2 = AE_MULFP32X2RAS(x0, x0);
  y0 = c0; t0 = c2;
  y2 = c1; t2 = c3;
  AE_MULAFP32X2RAS(t0, x2, y0); y0 = t0;
  AE_MULAFP32X2RAS(t2, x2, y2); y2 = t2;
  t0 = t1 = c4;
  t2 = t3 = c5;
  AE_MULAFP32X2RAS(t0, x2, y0); y0 = t0;
  AE_MULAFP32X2RAS(t2, x2, y2); y2 = t2;
  AE_MULAFP32X2RAS(y2, x0, y0); y0 = y2;


  y0 = AE_MULFP32X2RAS(x0, y0);
  x0 = AE_SUB32(y0, x0);
  AE_MULAFP32X2RAS(e0, x0, 0x02000000);

  AE_MOVT32X2(e0, 0x80000000, inf0);
  y = AE_MOVAD32_H(e0);
  return y;
} /* scl_log10_32x32() */
