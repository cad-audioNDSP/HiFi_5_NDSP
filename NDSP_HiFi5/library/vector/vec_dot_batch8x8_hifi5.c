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

/* Library API */
#include "NatureDSP_Signal_vector.h"
#include "NatureDSP_types.h"
#include "common.h"

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

void vec_dot_batch8x8  (int16_t   *restrict z, const int8_t *restrict x,const cint8ptr_t * restrict y,int rsh, int N, int M)
{
  int m, n, sh;

  ae_f64      q0, q1, q2, q3;
  ae_f64      q4, q5, q6, q7;
  ae_int16x4  vai0, vai1;
  ae_int8x8  x0, x1;
  ae_int8x8  y0_0, y0_1, y1_0, y1_1;
  ae_int8x8  y2_0, y2_1, y3_0, y3_1;
  ae_valignx2 ax, ay0, ay1, ay2, ay3;
  ae_valign   az;

  const ae_int8x16 * restrict px;
  const ae_int8x16 * restrict py0;
  const ae_int8x16 * restrict py1;
  const ae_int8x16 * restrict py2;
  const ae_int8x16 * restrict py3;
        ae_int16x4 * restrict pz;

  NASSERT(x);
  NASSERT(y);
  NASSERT(z);

  pz = (ae_int16x4 *)z;
  az = AE_ZALIGN64();

  if (M <= 0) return;

  if (N <= 0)
  {
    vai0 = 0;
    for (m = 0; m < M - 3; m += 4)
    {
      AE_SA16X4_IP(vai0, az, castxcc(ae_int16x4, pz));
    }
    AE_SA64POS_FP(az, castxcc(ae_int16x4, pz));
    for (; m < M; m++)
    {
      z[m] = 0;
    }
    return;
  }

  rsh = 48 - rsh;

  for (m = 0; m < (M >> 2); m++)
  {
    NASSERT(y[4 * m + 0]);
    NASSERT(y[4 * m + 1]);
    NASSERT(y[4 * m + 2]);
    NASSERT(y[4 * m + 3]);

    // Process a batch of 4 vectors;
    px  = (const ae_int8x16 *)x;
    py0 = (const ae_int8x16 *)y[4 * m + 0];
    py1 = (const ae_int8x16 *)y[4 * m + 1];
    py2 = (const ae_int8x16 *)y[4 * m + 2];
    py3 = (const ae_int8x16 *)y[4 * m + 3];
    q0 = AE_ZERO64();
    q1 = AE_ZERO64();
    q2 = AE_ZERO64();
    q3 = AE_ZERO64();
    q4 = AE_ZERO64();
    q5 = AE_ZERO64();
    q6 = AE_ZERO64();
    q7 = AE_ZERO64();
    sh = 0;
    sh = (((uintptr_t)(x)) & 15);
    sh = XT_MAX(0, (XT_MIN(16, N) - sh));
    if (sh > 0)
    {
      sh = sh & 15;
#if defined(AE_LAV8X8X2_XP)   
      {
        ax = AE_LA128_PP(px);
        AE_LAV8X8X2_XP(x0, x1, ax, px, sh);

        ay0 = AE_LA128_PP(py0);
        ay1 = AE_LA128_PP(py1);
        ay2 = AE_LA128_PP(py2);
        ay3 = AE_LA128_PP(py3);
        AE_LAV8X8X2_XP(y0_0, y0_1, ay0, py0, sh);
        AE_LAV8X8X2_XP(y1_0, y1_1, ay1, py1, sh);
        AE_LAV8X8X2_XP(y2_0, y2_1, ay2, py2, sh);
        AE_LAV8X8X2_XP(y3_0, y3_1, ay3, py3, sh);
        AE_MULAAAA2Q8(q0, q1, y0_0, x0);
        AE_MULAAAA2Q8(q2, q3, y1_0, x0);
        AE_MULAAAA2Q8(q4, q5, y2_0, x0);
        AE_MULAAAA2Q8(q6, q7, y3_0, x0);
        AE_MULAAAA2Q8(q0, q1, y0_1, x1);
        AE_MULAAAA2Q8(q2, q3, y1_1, x1);
        AE_MULAAAA2Q8(q4, q5, y2_1, x1);
        AE_MULAAAA2Q8(q6, q7, y3_1, x1);
      }
#else
      for (n = 0; n < (sh); n++)
      {        
        AE_L8_IP(x0, castxcc(ae_int8, px), sizeof(ae_int8));
        AE_L8_IP(y0_0, castxcc(ae_int8, py0), sizeof(ae_int8));
        AE_L8_IP(y1_0, castxcc(ae_int8, py1), sizeof(ae_int8));
        AE_L8_IP(y2_0, castxcc(ae_int8, py2), sizeof(ae_int8));
        AE_L8_IP(y3_0, castxcc(ae_int8, py3), sizeof(ae_int8));
        AE_MOVT8X16_L(y0_0, y0_1, y0_0, AE_MOVDA8(0), (0x0000FFFF >> 1));
        AE_MOVT8X16_L(y1_0, y1_1, y1_0, AE_MOVDA8(0), (0x0000FFFF >> 1));
        AE_MOVT8X16_L(y2_0, y2_1, y2_0, AE_MOVDA8(0), (0x0000FFFF >> 1));
        AE_MOVT8X16_L(y3_0, y3_1, y3_0, AE_MOVDA8(0), (0x0000FFFF >> 1));
        
        AE_MULAAAA2Q8(q0, q1, y0_0, x0);
        AE_MULAAAA2Q8(q2, q3, y1_0, x0);
        AE_MULAAAA2Q8(q4, q5, y2_0, x0);
        AE_MULAAAA2Q8(q6, q7, y3_0, x0);
      }
#endif
    }
    ax  = AE_LA128_PP(px);
    ay0 = AE_LA128_PP(py0);
    ay1 = AE_LA128_PP(py1);
    ay2 = AE_LA128_PP(py2);
    ay3 = AE_LA128_PP(py3);
    for (n = 0; n < ((N - sh) >> 4); n++)
    {
      AE_L8X8X2_IP(x0, x1, px, 2 * sizeof(ae_int8x8));
      AE_LA8X8X2_IP(y0_0, y0_1, ay0, py0);
      AE_LA8X8X2_IP(y1_0, y1_1, ay1, py1);
      AE_LA8X8X2_IP(y2_0, y2_1, ay2, py2);
      AE_LA8X8X2_IP(y3_0, y3_1, ay3, py3);
      AE_MULAAAA2Q8(q0, q1, y0_0, x0);
      AE_MULAAAA2Q8(q2, q3, y1_0, x0);
      AE_MULAAAA2Q8(q4, q5, y2_0, x0);
      AE_MULAAAA2Q8(q6, q7, y3_0, x0);
      AE_MULAAAA2Q8(q0, q1, y0_1, x1);
      AE_MULAAAA2Q8(q2, q3, y1_1, x1);
      AE_MULAAAA2Q8(q4, q5, y2_1, x1);
      AE_MULAAAA2Q8(q6, q7, y3_1, x1);
    }
#if defined(AE_LAV16X4X2_XP)   
    if (((N - sh) & 15))
    {
      sh = ((N - sh) & 15);
      ax = AE_LA128_PP(px);
      AE_LAV8X8X2_XP(x0, x1, ax, px, sh);

      ay0 = AE_LA128_PP(py0);
      ay1 = AE_LA128_PP(py1);
      ay2 = AE_LA128_PP(py2);
      ay3 = AE_LA128_PP(py3);
      AE_LAV8X8X2_XP(y0_0, y0_1, ay0, py0, sh);
      AE_LAV8X8X2_XP(y1_0, y1_1, ay1, py1, sh);
      AE_LAV8X8X2_XP(y2_0, y2_1, ay2, py2, sh);
      AE_LAV8X8X2_XP(y3_0, y3_1, ay3, py3, sh);
      AE_MULAAAA2Q8(q0, q1, y0_0, x0);
      AE_MULAAAA2Q8(q2, q3, y1_0, x0);
      AE_MULAAAA2Q8(q4, q5, y2_0, x0);
      AE_MULAAAA2Q8(q6, q7, y3_0, x0);
      AE_MULAAAA2Q8(q0, q1, y0_1, x1);
      AE_MULAAAA2Q8(q2, q3, y1_1, x1);
      AE_MULAAAA2Q8(q4, q5, y2_1, x1);
      AE_MULAAAA2Q8(q6, q7, y3_1, x1);
    }
#else
    {
      ae_int8x8 tmp;
      sh = ((N - sh) & 15);
      ax = AE_LA128_PP(px);
      AE_LA8X8X2_IP(x0, x1, ax, px);

      ay0 = AE_LA128_PP(py0);
      ay1 = AE_LA128_PP(py1);
      ay2 = AE_LA128_PP(py2);
      ay3 = AE_LA128_PP(py3);
      AE_LA8X8X2_IP(y0_0, y0_1, ay0, py0);
      AE_LA8X8X2_IP(y1_0, y1_1, ay1, py1);
      AE_LA8X8X2_IP(y2_0, y2_1, ay2, py2);
      AE_LA8X8X2_IP(y3_0, y3_1, ay3, py3);

      AE_MOVT8X16_H(y0_0, tmp, y0_0, AE_MOVDA8(0), (0xFFFF0000 >> sh));
      AE_MOVT8X16_H(y1_0, tmp, y1_0, AE_MOVDA8(0), (0xFFFF0000 >> sh));
      AE_MOVT8X16_H(y2_0, tmp, y2_0, AE_MOVDA8(0), (0xFFFF0000 >> sh));
      AE_MOVT8X16_H(y3_0, tmp, y3_0, AE_MOVDA8(0), (0xFFFF0000 >> sh));
      AE_MOVT8X16_H(tmp, y0_1, y0_1, AE_MOVDA8(0), (0xFFFF0000 >> sh));
      AE_MOVT8X16_H(tmp, y1_1, y1_1, AE_MOVDA8(0), (0xFFFF0000 >> sh));
      AE_MOVT8X16_H(tmp, y2_1, y2_1, AE_MOVDA8(0), (0xFFFF0000 >> sh));
      AE_MOVT8X16_H(tmp, y3_1, y3_1, AE_MOVDA8(0), (0xFFFF0000 >> sh));

      AE_MULAAAA2Q8(q0, q1, y0_0, x0);
      AE_MULAAAA2Q8(q2, q3, y1_0, x0);
      AE_MULAAAA2Q8(q4, q5, y2_0, x0);
      AE_MULAAAA2Q8(q6, q7, y3_0, x0);
      AE_MULAAAA2Q8(q0, q1, y0_1, x1);
      AE_MULAAAA2Q8(q2, q3, y1_1, x1);
      AE_MULAAAA2Q8(q4, q5, y2_1, x1);
      AE_MULAAAA2Q8(q6, q7, y3_1, x1);
    }
#endif
    q0 = AE_ADD64(q0, q1);
    q1 = AE_ADD64(q2, q3);
    q2 = AE_ADD64(q4, q5);
    q3 = AE_ADD64(q6, q7);

    vai0 = AE_TRUNCA16X4F64S(q0, q0, rsh);
    vai1 = AE_TRUNCA16X4F64S(q1, q1, rsh);
    AE_S16_0_IP(vai0, castxcc(ae_int16, pz), sizeof(ae_int16));
    AE_S16_0_IP(vai1, castxcc(ae_int16, pz), sizeof(ae_int16));
    vai0 = AE_TRUNCA16X4F64S(q2, q2, rsh);
    vai1 = AE_TRUNCA16X4F64S(q3, q3, rsh);
    AE_S16_0_IP(vai0, castxcc(ae_int16, pz), sizeof(ae_int16));
    AE_S16_0_IP(vai1, castxcc(ae_int16, pz), sizeof(ae_int16));
  }

  if (M & 3)
  {
    ae_valign   x_align, y0_align;
    for (m = (M&(~3)); m < M; m++)
    {
      px = (const ae_int8x16 *)x;
      py0 = (const ae_int8x16 *)y[m];

      x_align = AE_LA64_PP(px);
      y0_align = AE_LA64_PP(py0);

      q0 = q1 = AE_ZERO64();
      for (n = 0; n < N - 7; n += 8)
      {
        AE_LA8X8_IP(x0, x_align, castxcc(ae_int8x8, px));
        AE_LA8X8_IP(y0_0, y0_align, castxcc(ae_int8x8, py0));

        // Q16.47 <- [ Q16*Q16 + 1] + [ Q16*Q16 + 1 ] + [ Q16*Q16 + 1] + [ Q16*Q16 + 1 ]
        AE_MULAAAA2Q8(q0, q1, y0_0, x0);
      }
      // tail
      for (n = 0; n < (N&7); n ++)
      {
        AE_L8_IP(x0, castxcc(ae_int8, px), sizeof(ae_int8));
        AE_L8_IP(y0_0, castxcc(ae_int8, py0), sizeof(ae_int8));
        AE_MOVT8X16_L(y0_0, y0_1, y0_0, AE_MOVDA8(0), (0x0000FFFF >> 1));
        AE_MULAAAA2Q8(q0, q2, y0_0, x0);
      }
      
      q0 = AE_ADD64(q0, q1);
      vai0 = AE_TRUNCA16X4F64S(q0, q0, rsh);
      z[m] = AE_MOVAD16_0(vai0);
    }
  }
} /* vec_dot_batch8x8() */
