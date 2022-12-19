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
  NatureDSP Signal Processing Library. Vector Mathematics
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

#include "common.h"

void vec_atan32x32 (int32_t * restrict z, 
              const int32_t * restrict x, 
              int N )
{
  int n;
        ae_int32x4 * restrict pZ = (      ae_int32x4 *)z;
  const ae_int32x4 * restrict pX = (const ae_int32x4 *)x;
  const ae_int32x4 * restrict pX1 = (const ae_int32x4 *)x;

  ae_int32x2 y0, y1, y2, y3, y4, y5;
  ae_int32x2 t0, t1, t2, t3, t4, t5;
  ae_int32x2 p0, p1, p2, p3, p4, p5;
  ae_int32x2 p6, p7, p8, p9, p10;
  ae_int32x2 x0, x1, x0_2, x1_2, x0_3, x1_3;
  ae_valignx2 aX, aZ, aX1;
  aX = AE_LA128_PP(pX);
  aX1 = AE_LA128_PP(pX1);
  aZ = AE_ZALIGN128();

  if (N <= 0) return;

  for (n = 0; n<(N >> 2); n++)
  {
    AE_LA32X2X2_IP(x0, x1, aX, pX);

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
    x1 = AE_ABS32S(x1);
    x0 = AE_SUB32(x0, 0x40000000);
    x1 = AE_SUB32(x1, 0x40000000);
    #if 0
    t0 = t1 = p0; y0 = y1 = p1;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, p0, p0); t0 = y0; t1 = y1;   
    y0 = y1 = p2;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p3;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p4;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p5;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p6;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p7;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p8;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p9;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p10;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
  
    #elif 0
    AE_MULF2P32X4RAS(x0_2, x1_2, x0, x1, x0, x1);//x^2
    AE_MULF2P32X4RAS(x0_3, x1_3, x0_2, x1_2, x0_2, x1_2);//x^4

    t0 = t1 = p0; y0 = y1 = p4;
    AE_MULAF2P32X4RAS(y0, y1, x0_3, x1_3, p0, p0); t0 = y0; t1 = y1;
    t2 = t3 = p1; y2 = y3 = p5;
    AE_MULAF2P32X4RAS(y2, y3, x0_3, x1_3, p1, p1); t2 = y2; t3 = y3;
    t4 = t5 = p2; y4 = y5 = p6;
    AE_MULAF2P32X4RAS(y4, y5, x0_3, x1_3, p2, p2); t4 = y4; t5 = y5;

    y0 = y1 = p8;
    AE_MULAF2P32X4RAS(y0, y1, x0_3, x1_3, t0, t1); t0 = y0; t1 = y1;
    y2 = y3 = p9;
    AE_MULAF2P32X4RAS(y2, y3, x0_3, x1_3, t2, t3); t2 = y2; t3 = y3;
    y4 = y5 = p10;
    AE_MULAF2P32X4RAS(y4, y5, x0_3, x1_3, t4, t5); t4 = y4; t5 = y5;


    t6 = t7 = p3; y6 = y7 = p7;
    AE_MULAF2P32X4RAS(y6, y7, x0_3, x1_3, p3, p3); t6 = y6; t7 = y7;
    AE_MULF2P32X4RAS(y6, y7, x0_2, x1_2, y6, y7);

    AE_MULAF2P32X4RAS(y4, y5, x0, x1, y6, y7);

    AE_MULAF2P32X4RAS(y4, y5, x0_2, x1_2, y0, y1);
    AE_MULAF2P32X4RAS(y4, y5, x0, x1, y2, y3);y0 = y4; y1 = y5;
    #else
    AE_MULF2P32X4RAS(x0_2, x1_2, x0, x1, x0, x1);//x^2
    AE_MULF2P32X4RAS(x0_3, x1_3, x0_2, x1_2, x0, x1);//x^3

    t0 = t1 = p0; y0 = y1 = p3;
    AE_MULAF2P32X4RAS(y0, y1, x0_3, x1_3, p0, p0); t0 = y0; t1 = y1;
    t2 = t3 = p1; y2 = y3 = p4;
    AE_MULAF2P32X4RAS(y2, y3, x0_3, x1_3, p1, p1); t2 = y2; t3 = y3;
    t4 = t5 = p2; y4 = y5 = p5;
    AE_MULAF2P32X4RAS(y4, y5, x0_3, x1_3, p2, p2); t4 = y4; t5 = y5;

    y0 = y1 = p6;
    AE_MULAF2P32X4RAS(y0, y1, x0_3, x1_3, t0, t1); t0 = y0; t1 = y1;
    y2 = y3 = p7;
    AE_MULAF2P32X4RAS(y2, y3, x0_3, x1_3, t2, t3); t2 = y2; t3 = y3;
    y4 = y5 = p8;
    AE_MULAF2P32X4RAS(y4, y5, x0_3, x1_3, t4, t5); t4 = y4; t5 = y5;

    y0 = y1 = p9;
    AE_MULAF2P32X4RAS(y0, y1, x0_3, x1_3, t0, t1); t0 = y0; t1 = y1;
    y2 = y3 = p10;
    AE_MULAF2P32X4RAS(y2, y3, x0_3, x1_3, t2, t3); t2 = y2; t3 = y3;

    //y = L_add_ll(satQ31(((int64_t)(x)*y1 + (int64_t)(x2)*y3 + (1LL << 30)) >> 31), y2);
    AE_MULAF2P32X4RAS(y2, y3, x0_2, x1_2, y4, y5); 

    AE_MULAF2P32X4RAS(y2, y3, x0, x1, y0, y1); y0 = y2; y1 = y3;
    #endif
    
    AE_LA32X2X2_IP(x0, x1, aX1, pX1);
//    sx0 = AE_LT32(x0, AE_ZERO32());
//    sx1 = AE_LT32(x1, AE_ZERO32());
//
//    z0 = AE_NEG32S(y0);
//    z1 = AE_NEG32S(y1);
//    AE_MOVT32X2(y0, z0, sx0);
//    AE_MOVT32X2(y1, z1, sx1);
    y0=AE_MOVNEG32S_T(y0,x0);
    y1=AE_MOVNEG32S_T(y1,x1);

    AE_SA32X2X2_IP(y0, y1, aZ, pZ);
  }
  AE_SA128POS_FP(aZ, pZ);
  x += (N&~3);
  z += (N&~3);
  N &= 3;
  if (N>0)
  {
    int32_t ALIGN(32) scratch[4];
    ae_int32x4 *pScr;
    pScr = (ae_int32x4*)scratch;
    pX = (const ae_int32x4*)x;
    pZ = (      ae_int32x4*)z;
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
    p0 = AE_L32_I((const ae_int32 *)atan_table, 4 * 0);
    p1 = AE_L32_I((const ae_int32 *)atan_table, 4 * 1);
    p2 = AE_L32_I((const ae_int32 *)atan_table, 4 * 2);
    p3 = AE_L32_I((const ae_int32 *)atan_table, 4 * 3);
    p4 = AE_L32_I((const ae_int32 *)atan_table, 4 * 4);
    p5 = AE_L32_I((const ae_int32 *)atan_table, 4 * 5);
    p6 = AE_L32_I((const ae_int32 *)atan_table, 4 * 6);
    p7 = AE_L32_I((const ae_int32 *)atan_table, 4 * 7);
    p8 = AE_L32_X((const ae_int32 *)atan_table, 4 * 8);
    p9 = AE_L32_X((const ae_int32 *)atan_table, 4 * 9);
    p10 = AE_L32_X((const ae_int32 *)atan_table, 4 * 10);

    x0 = AE_ABS32S(x0);
    x1 = AE_ABS32S(x1);
    x0 = AE_SUB32(x0, 0x40000000);
    x1 = AE_SUB32(x1, 0x40000000);

    AE_MULF2P32X4RAS(x0_2, x1_2, x0, x1, x0, x1);//x^2
    AE_MULF2P32X4RAS(x0_3, x1_3, x0_2, x1_2, x0, x1);//x^3

    t0 = t1 = p0; y0 = y1 = p3;
    AE_MULAF2P32X4RAS(y0, y1, x0_3, x1_3, p0, p0); t0 = y0; t1 = y1;
    t2 = t3 = p1; y2 = y3 = p4;
    AE_MULAF2P32X4RAS(y2, y3, x0_3, x1_3, p1, p1); t2 = y2; t3 = y3;
    t4 = t5 = p2; y4 = y5 = p5;
    AE_MULAF2P32X4RAS(y4, y5, x0_3, x1_3, p2, p2); t4 = y4; t5 = y5;

    y0 = y1 = p6;
    AE_MULAF2P32X4RAS(y0, y1, x0_3, x1_3, t0, t1); t0 = y0; t1 = y1;
    y2 = y3 = p7;
    AE_MULAF2P32X4RAS(y2, y3, x0_3, x1_3, t2, t3); t2 = y2; t3 = y3;
    y4 = y5 = p8;
    AE_MULAF2P32X4RAS(y4, y5, x0_3, x1_3, t4, t5); t4 = y4; t5 = y5;

    y0 = y1 = p9;
    AE_MULAF2P32X4RAS(y0, y1, x0_3, x1_3, t0, t1); t0 = y0; t1 = y1;
    y2 = y3 = p10;
    AE_MULAF2P32X4RAS(y2, y3, x0_3, x1_3, t2, t3); t2 = y2; t3 = y3;

    //y = L_add_ll(satQ31(((int64_t)(x)*y1 + (int64_t)(x2)*y3 + (1LL << 30)) >> 31), y2);
    AE_MULAF2P32X4RAS(y2, y3, x0_2, x1_2, y4, y5);

    AE_MULAF2P32X4RAS(y2, y3, x0, x1, y0, y1); y0 = y2; y1 = y3;
    AE_L32X2X2_I(x0, x1, pScr, 0 * sizeof(ae_int32x4));
    y0=AE_MOVNEG32S_T(y0,x0);
    y1=AE_MOVNEG32S_T(y1,x1);

    AE_S32X2X2_I(y0, y1, pScr, 0 * sizeof(ae_int32x4));
    __Pragma("no_unroll")
    for (n = 0; n<N; n++)
    {
      ae_int32x2 t;
      AE_L32_IP(t, castxcc(ae_int32, pScr), sizeof(int32_t));
      AE_S32_L_IP(t, castxcc(ae_int32, pZ), sizeof(int32_t));
    }
  }
} /* vec_atan32x32() */
