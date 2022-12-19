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
  NatureDSP Signal Processing Library. Vector Mathematics
   Sine
    Code optimized for HiFi5
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

void vec_sine32x32 (  int32_t * restrict y, const int32_t * restrict x, int N)
{
  int n;
        ae_int32x4 * restrict pY = (      ae_int32x4 *)y;
  const ae_int32x4 * restrict pX = (const ae_int32x4 *)x;
  ae_int32x2 y0, y1, z0, z1, t0, t1;
  ae_int32x2 p0, p1, p2, p3, p4, p5;
  ae_int32x2 x0, x1, x0_2, x1_2, x0_3, x1_3;
  ae_valignx2 aX, aY;
  aX = AE_LA128_PP(pX);
  aY = AE_ZALIGN128();

  if (N<=0) return;

  for (n=0; n<(N>>2); n++)
  {
    AE_LA32X2X2_IP(x0, x1, aX, pX);

    p0 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 0);
    p1 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 1);
    p2 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 2);
    p3 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 3);
    p4 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 4);
    p5 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 5);

    AE_MULF2P32X4RAS(x0_2, x1_2, x0, x1, x0, x1);//x^2
    AE_MULF2P32X4RAS(x0_3, x1_3, x0_2, x1_2, x0, x1);//x^3

    t0 = t1 = p0; y0 = y1 = p1;
    AE_MULAF2P32X4RAS(y0, y1, x0_2, x1_2, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p2;
    AE_MULAF2P32X4RAS(y0, y1, x0_2, x1_2, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p3;
    AE_MULAF2P32X4RAS(y0, y1, x0_2, x1_2, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p4;
    AE_MULAF2P32X4RAS(y0, y1, x0_2, x1_2, t0, t1); t0 = y0; t1 = y1;
    AE_MULF2P32X4RAS(y0, y1, x0, x1, p5, p5);
    AE_MULAF2P32X4RAS(y0, y1, x0_3, x1_3, t0, t1); //Q28
    AE_MUL2P32X4S(z0, z1, y0, y1, 0x8, 0x8); //Q31
    AE_SA32X2X2_IP(z0, z1, aY, pY);
  }
  AE_SA128POS_FP(aY, pY);
  x += (N&~3);
  y += (N&~3);
  N &= 3;
  if (N>0)
  {
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
    p0 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 0);
    p1 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 1);
    p2 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 2);
    p3 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 3);
    p4 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 4);
    p5 = AE_L32_I((const ae_int32 *)sine_table32, 4 * 5);

    AE_MULF2P32X4RAS(x0_2, x1_2, x0, x1, x0, x1);//x^2
    AE_MULF2P32X4RAS(x0_3, x1_3, x0_2, x1_2, x0, x1);//x^3

    t0 = t1 = p0; y0 = y1 = p1;
    AE_MULAF2P32X4RAS(y0, y1, x0_2, x1_2, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p2;
    AE_MULAF2P32X4RAS(y0, y1, x0_2, x1_2, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p3;
    AE_MULAF2P32X4RAS(y0, y1, x0_2, x1_2, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p4;
    AE_MULAF2P32X4RAS(y0, y1, x0_2, x1_2, t0, t1); t0 = y0; t1 = y1;
    AE_MULF2P32X4RAS(y0, y1, x0, x1, p5, p5);
    AE_MULAF2P32X4RAS(y0, y1, x0_3, x1_3, t0, t1); //Q28
    AE_MUL2P32X4S(z0, z1, y0, y1, 0x8, 0x8); //Q31
    AE_S32X2X2_I(z0, z1, pScr, 0 * sizeof(ae_int32x4));
    __Pragma("no_unroll")
    for (n = 0; n<N; n++)
    {
      ae_int32x2 t;
      AE_L32_IP(t, castxcc(ae_int32, pScr), sizeof(int32_t));
      AE_S32_L_IP(t, castxcc(ae_int32, pY), sizeof(int32_t));
    }
  }
} /* vec_sine32x32() */
