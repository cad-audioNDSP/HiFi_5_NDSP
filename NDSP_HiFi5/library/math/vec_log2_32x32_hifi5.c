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
  NatureDSP Signal Processing Library. Vector Mathematics
   Logarithm, base 2
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

void vec_log2_32x32 (int32_t * restrict y, const int32_t * restrict x, int N)
{
  int n;
  const ae_int32x4 * restrict pX = (const ae_int32x4 *)x;
  const ae_int32x4 * restrict pX1 = (const ae_int32x4 *)x;
        ae_int32x4 * restrict pY = (      ae_int32x4 *)y;
  ae_int32x2 y0, y1, y2, y3, t0, t1, t2, t3;
  ae_int32x2 vmw, vsw, vhw, vzw;
  xtbool2     inf0, inf1, i0, i1;
  ae_int16x4 nsa01;
  ae_valignx2 aX, aX1, aY;
  aX = AE_LA128_PP(pX);
  aX1 = AE_LA128_PP(pX1);
  aY = AE_ZALIGN128();

  if (N<=0) return;
  vmw = AE_MOVDA32X2(0x7fffff, 0x7fffff);
  vsw = AE_MOVDA32X2(0x400000, 0x400000);
  vhw = AE_MOVDA32X2(0x5A82799A, 0x5A82799A);
  vzw = AE_MOVI(0);
  for (n=0; n<(N>>2); n++)
  {
    ae_int32x2 x0, x1, x2, x3, e0, e1, d0, d1;
    ae_int32x2 c0, c1, c2, c3, c4, c5;
    AE_LA32X2X2_IP(x0, x1, aX, pX);
    inf0 = AE_LE32(x0, vzw);
    inf1 = AE_LE32(x1, vzw);
    
    /* Normalize x*/
    nsa01 = AE_NSA32X4(x0, x1);
    x0=AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(nsa01)));
    x1=AE_SRAV32RS(x1, AE_SEXT32X2D16_10(AE_NEG16S(nsa01)));
    nsa01 = AE_SUB16(16, nsa01);

    i0 = AE_LE32(x0, vhw);
    i1 = AE_LE32(x1, vhw);
    
    e0 = AE_SEXT32X2D16_32(nsa01);
    e1 = AE_SEXT32X2D16_10(nsa01);

    AE_MOVT32X2(e0, AE_SUB32(e0, 1), i0);
    AE_MOVT32X2(e1, AE_SUB32(e1, 1), i1);

    e0 = AE_SLAI32(e0, (31 - 6));
    e1 = AE_SLAI32(e1, (31 - 6));

    c0 = AE_L32_I((const ae_int32 *)log2_tbl_Q31, 4 * 0);
    c1 = AE_L32_I((const ae_int32 *)log2_tbl_Q31, 4 * 1);
    c2 = AE_L32_I((const ae_int32 *)log2_tbl_Q31, 4 * 2);
    c3 = AE_L32_I((const ae_int32 *)log2_tbl_Q31, 4 * 3);
    c4 = AE_L32_I((const ae_int32 *)log2_tbl_Q31, 4 * 4);
    c5 = AE_L32_I((const ae_int32 *)log2_tbl_Q31, 4 * 5);

    t0 = t1 = 0x80000000; d0 = d1 = 1;
    AE_MOVT32X2(d0, 2, i0);
    AE_MOVT32X2(d1, 2, i1);

    AE_MULS2P32X4(t0, t1, x0, x1, d0, d1); x0 = t0; x1 = t1;

    AE_MULF2P32X4RAS(x2, x3, x0, x1, x0, x1);
    y0 = y1 = c0; t0 = t1 = c2;
    y2 = y3 = c1; t2 = t3 = c3;
    AE_MULAF2P32X4RAS(t0, t1, x2, x3, y0, y1); y0 = t0; y1 = t1;
    AE_MULAF2P32X4RAS(t2, t3, x2, x3, y2, y3); y2 = t2; y3 = t3;
    t0 = t1 = c4;
    t2 = t3 = c5;
    AE_MULAF2P32X4RAS(t0, t1, x2, x3, y0, y1); y0 = t0; y1 = t1;
    AE_MULAF2P32X4RAS(t2, t3, x2, x3, y2, y3); y2 = t2; y3 = t3;
    AE_MULAF2P32X4RAS(y2, y3, x0, x1, y0, y1); y0 = y2; y1 = y3;


    AE_MULF2P32X4RAS(y0, y1, x0, x1, y0, y1);
    x0 = AE_SUB32(y0, x0);
    x1 = AE_SUB32(y1, x1);
    AE_MULAF2P32X4RAS(e0, e1, x0, x1, 0x02000000, 0x02000000);
    y0 = e0; y1 = e1;

    AE_MOVT32X2(y0, 0x80000000, inf0);
    AE_MOVT32X2(y1, 0x80000000, inf1);

    AE_SA32X2X2_IP(y0, y1, aY, pY);
  }
  AE_SA128POS_FP(aY, pY);
  x += (N&~3);
  y += (N&~3);
  N &= 3;
  if (N>0)
  {
    ae_int32x2 x0, x1, x2, x3, e0, e1, d0, d1;
    ae_int32x2 c0, c1, c2, c3, c4, c5;
    int32_t ALIGN(32) scratch[4];
    ae_int32x4 *pScr;
    pScr = (ae_int32x4*)scratch;
    pX = (const ae_int32x4*)x;
    pY = (      ae_int32x4*)y;
    AE_S32X2X2_I(0, 0, pScr, 0);
    __Pragma("no_unroll")
    for (n = 0; n<N; n++)
    {
      ae_int32x2 t;
      AE_L32_IP(t, castxcc(ae_int32, pX), sizeof(int32_t));
      AE_S32_L_IP(t, castxcc(ae_int32, pScr), sizeof(int32_t));
    }
    pScr = (ae_int32x4*)scratch;
    AE_L32X2X2_I(x0, x1, pScr, 0 * sizeof(ae_int32x4));

    inf0 = AE_LE32(x0, vzw);
    inf1 = AE_LE32(x1, vzw);

    /* Normalize x*/
    nsa01 = AE_NSA32X4(x0, x1);
    x0 = AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(nsa01)));
    x1 = AE_SRAV32RS(x1, AE_SEXT32X2D16_10(AE_NEG16S(nsa01)));
    nsa01 = AE_SUB16(16, nsa01);

    i0 = AE_LE32(x0, vhw);
    i1 = AE_LE32(x1, vhw);

    e0 = AE_SEXT32X2D16_32(nsa01);
    e1 = AE_SEXT32X2D16_10(nsa01);

    AE_MOVT32X2(e0, AE_SUB32(e0, 1), i0);
    AE_MOVT32X2(e1, AE_SUB32(e1, 1), i1);

    e0 = AE_SLAI32(e0, (31 - 6));
    e1 = AE_SLAI32(e1, (31 - 6));

    c0 = AE_L32_I((const ae_int32 *)log2_tbl_Q31, 4 * 0);
    c1 = AE_L32_I((const ae_int32 *)log2_tbl_Q31, 4 * 1);
    c2 = AE_L32_I((const ae_int32 *)log2_tbl_Q31, 4 * 2);
    c3 = AE_L32_I((const ae_int32 *)log2_tbl_Q31, 4 * 3);
    c4 = AE_L32_I((const ae_int32 *)log2_tbl_Q31, 4 * 4);
    c5 = AE_L32_I((const ae_int32 *)log2_tbl_Q31, 4 * 5);

    t0 = t1 = 0x80000000; d0 = d1 = 1;
    AE_MOVT32X2(d0, 2, i0);
    AE_MOVT32X2(d1, 2, i1);

    AE_MULS2P32X4(t0, t1, x0, x1, d0, d1); x0 = t0; x1 = t1;

    AE_MULF2P32X4RAS(x2, x3, x0, x1, x0, x1);
    y0 = y1 = c0; t0 = t1 = c2;
    y2 = y3 = c1; t2 = t3 = c3;
    AE_MULAF2P32X4RAS(t0, t1, x2, x3, y0, y1); y0 = t0; y1 = t1;
    AE_MULAF2P32X4RAS(t2, t3, x2, x3, y2, y3); y2 = t2; y3 = t3;
    t0 = t1 = c4;
    t2 = t3 = c5;
    AE_MULAF2P32X4RAS(t0, t1, x2, x3, y0, y1); y0 = t0; y1 = t1;
    AE_MULAF2P32X4RAS(t2, t3, x2, x3, y2, y3); y2 = t2; y3 = t3;
    AE_MULAF2P32X4RAS(y2, y3, x0, x1, y0, y1); y0 = y2; y1 = y3;


    AE_MULF2P32X4RAS(y0, y1, x0, x1, y0, y1);
    x0 = AE_SUB32(y0, x0);
    x1 = AE_SUB32(y1, x1);
    AE_MULAF2P32X4RAS(e0, e1, x0, x1, 0x02000000, 0x02000000);
    y0 = e0; y1 = e1;

    AE_MOVT32X2(y0, 0x80000000, inf0);
    AE_MOVT32X2(y1, 0x80000000, inf1);
    AE_S32X2X2_I(y0, y1, pScr, 0 * sizeof(ae_int32x4));
    __Pragma("no_unroll")
    for (n = 0; n<N; n++)
    {
      ae_int32x2 t;
      AE_L32_IP(t, castxcc(ae_int32, pScr), sizeof(int32_t));
      AE_S32_L_IP(t, castxcc(ae_int32, pY), sizeof(int32_t));
    }
  }
} /* vec_log2_32x32() */
