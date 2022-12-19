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
  NatureDSP Signal Processing Library. Vector Operations
    Common Exponent 
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_vector.h"
#include "common.h"
#include "common_fpu.h"




/*-------------------------------------------------------------------------
  Batch Computation of Vector Dot products
  These routines take a set of input vectors and compute their dot product 
  with specific reference data.
  Two versions of routines are available: 
  - regular versions (vec_dot_batch8x8, vec_dot_batch8x16, vec_dot_batch16x16, 
    vec_dot_batchf). They work with arbitratry arguments
  - fast versions (vec_dot_batch8x8_fast, vec_dot_batch8x16_fast, 
    vec_dot_batch16x16_fast, vec_dot_batchf_fast) apply some restrictions.  

  Precision: 
  8x8    8x8-bit data, 16-bit output (fractional multiply Q7xQ7->Q15)
  8x16   8x16-bit data, 16-bit output (fractional multiply Q7xQ15->Q15)
  16x16  16x16-bit data, 32-bit output (fractional multiply Q15xQ15->Q31)
  f      single precision floating point
  fp16   half precision floating point

  Input:
  x[N]     input (reference) data, Q7, Q15 or floating point
  y[M][N]  pointers to M input vectors, Q7, Q15 or floating point
  N        length of vectors
  M        number of vectors
  rsh      right shift for output (for fixed point API only!)
  Output:
  z[M]     dot products between references and M inputs, Q15, Q31 or 
           floating point

  Restrictions:
  Regular versions:
    none
  Faster versions:
    x,y[m] - aligned on 16-byte boundary
    N      - multiple of 8
    M        multiple of 4
-------------------------------------------------------------------------*/
#if !HAVE_HPFPU
DISCARD_FUN(void,vec_dot_batch_fp16_fast,(float16_t *z, const float16_t   * x,const cfloat16ptr_t * y        , int N, int M))
#else
void vec_dot_batch_fp16_fast    (float16_t *z, const float16_t   * x,const cfloat16ptr_t * y        , int N, int M)
{
  xthalfx4 x0, y0, y1, y2, y3, vacc0, vacc1, vacc2, vacc3;
  xthalfx4 x1, y4, y5, y6, y7, vacc4, vacc5, vacc6, vacc7;
  xthalfx4 res0, res1, res2, res3, res4, res5, res6, res7;
  xthalfx4 one;
  const xthalfx8 * restrict px = (const xthalfx8 *)x;
  const xthalfx8 * restrict py0;
  const xthalfx8 * restrict py1;
  const xthalfx8 * restrict py2;
  const xthalfx8 * restrict py3;
  ae_int16 * restrict pz = (ae_int16 *)z;
  int m, n;
  NASSERT(z);  NASSERT_ALIGN(z, HIFI_SIMD_WIDTH);
  NASSERT(x);  NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
  NASSERT(y);
  NASSERT(M % 4 == 0 && N % 8 == 0);
  if (M <= 0) return;
  if (N == 0)
  {
    xthalf zero = CONST_H(0);
    for (m = 0; m < M; ++m)
      AE_SHIP(zero, castxcc(xthalf, pz), sizeof(xthalf));
    return;
  }
  one = CONST_HX4(1);
  for (m = 0; m < (M>>2); m++)
  {
    px  = (const xthalfx8 *)x;
    py0 = (const xthalfx8 *)(y[4*m+0]);
    py1 = (const xthalfx8 *)(y[4*m+1]);
    py2 = (const xthalfx8 *)(y[4*m+2]);
    py3 = (const xthalfx8 *)(y[4*m+3]);

    MOV_HX4X2(vacc0, vacc1, ZERO_HX4(), ZERO_HX4());
    MOV_HX4X2(vacc2, vacc3, ZERO_HX4(), ZERO_HX4());
    MOV_HX4X2(vacc4, vacc5, ZERO_HX4(), ZERO_HX4());
    MOV_HX4X2(vacc6, vacc7, ZERO_HX4(), ZERO_HX4());
    for (n = 0; n < (N >> 3); n++)
    {
      AE_LHX4X2_IP(x0, x1,  px, 8 * sizeof(float16_t));           
      AE_LHX4X2_IP(y0, y4, py0, 8 * sizeof(float16_t));
      AE_LHX4X2_IP(y1, y5, py1, 8 * sizeof(float16_t));
      AE_LHX4X2_IP(y2, y6, py2, 8 * sizeof(float16_t));
      AE_LHX4X2_IP(y3, y7, py3, 8 * sizeof(float16_t));
      MADDQ_H(vacc0, vacc1, y0, y1, x0);
      MADDQ_H(vacc2, vacc3, y2, y3, x0);
      MADDQ_H(vacc4, vacc5, y4, y5, x1);
      MADDQ_H(vacc6, vacc7, y6, y7, x1);
    }
    ADD_HX4X2(vacc0, vacc1, vacc0, vacc1, vacc4, vacc5);
    ADD_HX4X2(vacc2, vacc3, vacc2, vacc3, vacc6, vacc7);
    
    MULCNVH_HX4X2(res0, res4, one, vacc0, vacc0);
    MULACNVL_HX4X2(res0, res4, one, vacc0, vacc0);
    MULCNVH_HX4X2(res1, res5, one, vacc1, vacc1);
    MULACNVL_HX4X2(res1, res5, one, vacc1, vacc1);
    ADD_HX4X2(res0, res1, res0, res1, res4, res5);

    AE_S16_0_IP(AE_MOVINT16X4_FROMF16X4(AE_MOVF16X4_FROMHALFX4(res0)), pz, sizeof(xthalf));
    AE_S16_0_IP(AE_MOVINT16X4_FROMF16X4(AE_MOVF16X4_FROMHALFX4(res1)), pz, sizeof(xthalf));

    
    MULCNVH_HX4X2(res2, res6, one, vacc2, vacc2);
    MULACNVL_HX4X2(res2, res6, one, vacc2, vacc2);
    MULCNVH_HX4X2(res3, res7, one, vacc3, vacc3);
    MULACNVL_HX4X2(res3, res7, one, vacc3, vacc3);
    ADD_HX4X2(res2, res3, res2, res3, res6, res7);
    AE_S16_0_IP(AE_MOVINT16X4_FROMF16X4(AE_MOVF16X4_FROMHALFX4(res2)), pz, sizeof(xthalf));
    AE_S16_0_IP(AE_MOVINT16X4_FROMF16X4(AE_MOVF16X4_FROMHALFX4(res3)), pz, sizeof(xthalf));

  }
} /* vec_dot_batch_fp16_fast() */
#endif
