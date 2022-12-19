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
#include "vec_alog_table.h"
#include "common.h"
#include "NatureDSP_Signal_math.h"

/*
  NatureDSP Signal Processing Library. Vector Mathematics
   Antilogarithm, base 10
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/
/*-------------------------------------------------------------------------
  Antilogarithm
  These routines calculate antilogarithm (base2, natural and base10). 
  Fixed-point functions accept inputs in Q25 and form outputs in Q16.15 
  format and return 0x7FFFFFFF in case of overflow and 0 in case of 
  underflow.

  Precision:
  32x32  32-bit inputs, 32-bit outputs. Accuracy: 8*e-6*y+1LSB
  f      floating point: Accuracy: 2 ULP
  NOTE:
  1.  Although 32 and 24 bit functions provide the similar accuracy, 32-bit
      functions have better input/output resolution (dynamic range).
  2.  Floating point functions are compatible with standard ANSI C routines 
      and set errno and exception flags accordingly.

  Input:
  x[N]  input data,Q25 or floating point 
  N     length of vectors
  Output:
  y[N]  output data,Q16.15 or floating point  

  Restriction:
  x,y should not overlap

  Scalar versions:
  ----------------
  fixed point functions return result in Q16.15

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  x,y - aligned on 16-byte boundary
  N   - multiple of 2
-------------------------------------------------------------------------*/
void vec_antilog10_32x32(int32_t * restrict y, const int32_t* restrict x, int N)
{
  int n;
        ae_int32x4 * restrict pY = (      ae_int32x4 *)y;
  const ae_int32x4 * restrict pX = (const ae_int32x4 *)x;
  ae_int32x2 y0, y1, t0, t1;
  ae_valignx2 aX, aY;
  aX = AE_LA128_PP(pX);
  aY = AE_ZALIGN128();

  if (N<=0) return;

  for (n=0; n<(N>>2); n++)
  {
    ae_int32x2 x0, x1, e0, e1;
    ae_int64 ah0, al0;
    ae_int64 ah1, al1;

    AE_LA32X2X2_IP(x0, x1, aX, pX);
    AE_MUL32X2S_HH_LL(ah0, al0, 1783446566, x0);
    AE_MUL32X2S_HH_LL(ah1, al1, 1783446566, x1);
    x0 = AE_TRUNCI32X2F64S(ah0, al0, 3); x1 = AE_TRUNCI32X2F64S(ah1, al1, 3);

    e0 = AE_SRAI32(x0, 25);
    e1 = AE_SRAI32(x1, 25);
    e0 = AE_SUB32(15, e0);
    e1 = AE_SUB32(15, e1);

    x0 = AE_SLLI32(x0, 6);
    x1 = AE_SLLI32(x1, 6);
    x0 = AE_OR32(x0, 0x80000000);
    x1 = AE_OR32(x1, 0x80000000);

    y0 = y1 = AE_L32_I((const ae_int32 *)pow2poly, 4 * 0);
    t0 = t1 = AE_L32_I((const ae_int32 *)pow2poly, 4 * 1);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = AE_L32_I((const ae_int32 *)pow2poly, 4 * 2);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = AE_L32_I((const ae_int32 *)pow2poly, 4 * 3);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = AE_L32_I((const ae_int32 *)pow2poly, 4 * 4);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    y0 = AE_SRAV32RS(y0, e0);
    y1 = AE_SRAV32RS(y1, e1);

    AE_SA32X2X2_IP(y0, y1, aY, pY);
  }
  AE_SA128POS_FP(aY, pY);
  x += (N&~3);
  y += (N&~3);
  N &= 3;
  if (N>0)
  {
    ae_int32x2 x0, x1, e0, e1;
    ae_int64 ah0, al0;
    ae_int64 ah1, al1;
    int32_t ALIGN(32) scratch[4];
    ae_int32x4 *pScr;
    pScr=(      ae_int32x4*)scratch;
    pX  =(const ae_int32x4*)x;
    pY  =(      ae_int32x4*)y;
    AE_S32X2X2_I(0,0,pScr,0);
    __Pragma("no_unroll")
    for(n=0; n<N; n++) 
    {
      ae_int32x2 t;
      AE_L32_IP(t,castxcc(ae_int32,pX),sizeof(int32_t));
      AE_S32_L_IP(t,castxcc(ae_int32,pScr),sizeof(int32_t));
    }
    pScr=(ae_int32x4*)scratch;
    AE_L32X2X2_I (x0,x1,pScr,0*sizeof(ae_int32x4));
    AE_MUL32X2S_HH_LL(ah0, al0, 1783446566, x0);
    AE_MUL32X2S_HH_LL(ah1, al1, 1783446566, x1);
    x0 = AE_TRUNCI32X2F64S(ah0, al0, 3); x1 = AE_TRUNCI32X2F64S(ah1, al1, 3);

    e0 = AE_SRAI32(x0, 25);
    e1 = AE_SRAI32(x1, 25);
    e0 = AE_SUB32(15, e0);
    e1 = AE_SUB32(15, e1);

    x0 = AE_SLLI32(x0, 6);
    x1 = AE_SLLI32(x1, 6);
    x0 = AE_OR32(x0, 0x80000000);
    x1 = AE_OR32(x1, 0x80000000);

    y0 = y1 = AE_L32_I((const ae_int32 *)pow2poly, 4 * 0);
    t0 = t1 = AE_L32_I((const ae_int32 *)pow2poly, 4 * 1);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = AE_L32_I((const ae_int32 *)pow2poly, 4 * 2);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = AE_L32_I((const ae_int32 *)pow2poly, 4 * 3);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = AE_L32_I((const ae_int32 *)pow2poly, 4 * 4);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    y0 = AE_SRAV32RS(y0, e0);
    y1 = AE_SRAV32RS(y1, e1);
    AE_S32X2X2_I(y0, y1, pScr, 0 * sizeof(ae_int32x4));
    __Pragma("no_unroll")
    for (n = 0; n<N; n++)
    {
      ae_int32x2 t;
      AE_L32_IP(t, castxcc(ae_int32, pScr), sizeof(int32_t));
      AE_S32_L_IP(t, castxcc(ae_int32, pY), sizeof(int32_t));
    }
  }
} /* vec_antilog10_32x32() */
