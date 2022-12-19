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

void vec_dot_batch16x16(int32_t   * restrict z, const int16_t * restrict x,const cint16ptr_t * restrict y,int rsh, int N, int M)
{
  int m, n, sh;

  ae_f64      q0, q1, q2, q3;
  ae_f64      q4, q5, q6, q7;
  ae_int32x2  vai0, vai1;
  ae_int16x4  vxh;
  ae_int16x4  vyh0;
  ae_int16x4  x0, x1;
  ae_int16x4  y0_0, y0_1, y1_0, y1_1;
  ae_int16x4  y2_0, y2_1, y3_0, y3_1;
  ae_valignx2 ax, ay0, ay1, ay2, ay3;
  ae_valign   az;
  xtbool4     bmask;

  const ae_int16x8 * restrict px;
  const ae_int16x8 * restrict py0;
  const ae_int16x8 * restrict py1;
  const ae_int16x8 * restrict py2;
  const ae_int16x8 * restrict py3;
        ae_int32x2 * restrict pz;

  NASSERT(x);
  NASSERT(y);
  NASSERT(z);

  pz = (ae_int32x2 *)z;
  if (M <= 0) return;
  az = AE_ZALIGN64();
  if (N <= 0)
  {
    ae_valign ay1;
    ay0 = AE_ZALIGN128();
    ay1 = AE_ZALIGN64();
    vai0 = 0;
    for (m = 0; m < M - 3; m+=4)
    {
      AE_SA32X2X2_IP(vai0, vai0, ay0, castxcc(ae_int32x4, pz));
    }
    AE_SA128POS_FP(ay0, pz);
    if (M & 2) 
    {
      AE_SA32X2_IP(vai0, ay1, pz);
    }
    AE_SA64POS_FP(ay1, pz);
    if (M & 1)
    {
      z[M - 1] = 0;
    }
    return;
  }

  rsh = 32 - rsh;

  for (m = 0; m < (M>>2); m ++)
  {
    NASSERT(y[4 * m + 0]);
    NASSERT(y[4 * m + 1]);
    NASSERT(y[4 * m + 2]);
    NASSERT(y[4 * m + 3]);

    // Process a batch of 4 vectors;
    px = (const ae_int16x8 *)x;
    py0 = (const ae_int16x8 *)y[4 * m + 0];
    py1 = (const ae_int16x8 *)y[4 * m + 1];
    py2 = (const ae_int16x8 *)y[4 * m + 2];
    py3 = (const ae_int16x8 *)y[4 * m + 3];
    q0 = AE_ZERO64();
    q1 = AE_ZERO64();
    q2 = AE_ZERO64();
    q3 = AE_ZERO64();
    q4 = AE_ZERO64();
    q5 = AE_ZERO64();
    q6 = AE_ZERO64();
    q7 = AE_ZERO64();
    sh = 0;
    sh = (((uintptr_t)(x)) & 15) >> 1;
    sh = XT_MAX(0, (XT_MIN(8, N) - sh));

    if (sh > 0)
    {
      sh = sh & 7;
#if defined(AE_LAV16X4X2_XP)   
      {
        ax = AE_LA128_PP(px);
        AE_LAV16X4X2_XP(x0, x1, ax, px, 2 * sh);

        ay0 = AE_LA128_PP(py0);
        ay1 = AE_LA128_PP(py1);
        ay2 = AE_LA128_PP(py2);
        ay3 = AE_LA128_PP(py3);
        AE_LAV16X4X2_XP(y0_0, y0_1, ay0, py0, 2 * sh);
        AE_LAV16X4X2_XP(y1_0, y1_1, ay1, py1, 2 * sh);
        AE_LAV16X4X2_XP(y2_0, y2_1, ay2, py2, 2 * sh);
        AE_LAV16X4X2_XP(y3_0, y3_1, ay3, py3, 2 * sh);
        AE_MULAAAA2Q16(q0, q1, y0_0, y0_1, x0, x1);
        AE_MULAAAA2Q16(q2, q3, y1_0, y1_1, x0, x1);
        AE_MULAAAA2Q16(q4, q5, y2_0, y2_1, x0, x1);
        AE_MULAAAA2Q16(q6, q7, y3_0, y3_1, x0, x1);        
      }
#else
      // tail
      for (n = 0; n < (sh); n++)
      {
        AE_L16_IP(x0, castxcc(ae_int16, px), sizeof(int16_t));
        AE_L16_IP(y0_0, castxcc(ae_int16,py0), sizeof(int16_t));
        AE_L16_IP(y1_0, castxcc(ae_int16,py1), sizeof(int16_t));
        AE_L16_IP(y2_0, castxcc(ae_int16,py2), sizeof(int16_t));
        AE_L16_IP(y3_0, castxcc(ae_int16,py3), sizeof(int16_t));

        AE_MULA16_00(q0, x0, y0_0);
        AE_MULA16_00(q2, x0, y1_0);
        AE_MULA16_00(q4, x0, y2_0);
        AE_MULA16_00(q6, x0, y3_0);
      }
#endif
    }

    ax = AE_LA128_PP(px);
    ay0 = AE_LA128_PP(py0);
    ay1 = AE_LA128_PP(py1);
    ay2 = AE_LA128_PP(py2);
    ay3 = AE_LA128_PP(py3);
    for (n = 0; n < ((N - sh) >> 3); n++)
    {
      AE_L16X4X2_IP(x0, x1, px, 2 * sizeof(ae_int16x4));
      AE_LA16X4X2_IP(y0_0, y0_1, ay0, py0);
      AE_LA16X4X2_IP(y1_0, y1_1, ay1, py1);
      AE_LA16X4X2_IP(y2_0, y2_1, ay2, py2);
      AE_LA16X4X2_IP(y3_0, y3_1, ay3, py3);
      // Q16.47 <- [ Q16*Q16 + 1] + [ Q16*Q16 + 1 ] + [ Q16*Q16 + 1] + [ Q16*Q16 + 1 ]
      AE_MULAAAA2Q16(q0, q1, y0_0, y0_1, x0, x1);
      AE_MULAAAA2Q16(q2, q3, y1_0, y1_1, x0, x1);
      AE_MULAAAA2Q16(q4, q5, y2_0, y2_1, x0, x1);
      AE_MULAAAA2Q16(q6, q7, y3_0, y3_1, x0, x1);

    }
#if defined(AE_LAV16X4X2_XP)   
    if (((N - sh) & 7))
    {
      sh = ((N - sh) & 7);
      ax = AE_LA128_PP(px);
      AE_LAV16X4X2_XP(x0, x1, ax, px, 2 * sh);

      ay0 = AE_LA128_PP(py0);
      ay1 = AE_LA128_PP(py1);
      ay2 = AE_LA128_PP(py2);
      ay3 = AE_LA128_PP(py3);
      AE_LAV16X4X2_XP(y0_0, y0_1, ay0, py0, 2 * sh);
      AE_LAV16X4X2_XP(y1_0, y1_1, ay1, py1, 2 * sh);
      AE_LAV16X4X2_XP(y2_0, y2_1, ay2, py2, 2 * sh);
      AE_LAV16X4X2_XP(y3_0, y3_1, ay3, py3, 2 * sh);
      AE_MULAAAA2Q16(q0, q1, y0_0, y0_1, x0, x1);
      AE_MULAAAA2Q16(q2, q3, y1_0, y1_1, x0, x1);
      AE_MULAAAA2Q16(q4, q5, y2_0, y2_1, x0, x1);
      AE_MULAAAA2Q16(q6, q7, y3_0, y3_1, x0, x1);
    }
#else
    {
      ae_valign ax, ay0, ay1, ay2, ay3;
      // tail
      ax = AE_LA64_PP(px);
      ay0 = AE_LA64_PP(py0);
      ay1 = AE_LA64_PP(py1);
      ay2 = AE_LA64_PP(py2);
      ay3 = AE_LA64_PP(py3);
      if (((N - sh)& 7) >= 4)
      {
        AE_LA16X4_IP(x0, ax, castxcc(ae_int16x4, px));
        AE_LA16X4_IP(y0_0, ay0, castxcc(ae_int16x4, py0));
        AE_LA16X4_IP(y1_0, ay1, castxcc(ae_int16x4, py1));
        AE_LA16X4_IP(y2_0, ay2, castxcc(ae_int16x4, py2));
        AE_LA16X4_IP(y3_0, ay3, castxcc(ae_int16x4, py3));

        AE_MULAAAAQ16(q0, y0_0, x0);
        AE_MULAAAAQ16(q2, y1_0, x0);
        AE_MULAAAAQ16(q4, y2_0, x0);
        AE_MULAAAAQ16(q6, y3_0, x0);
      }
      if (((N - sh) & 3))
      {
        bmask = AE_MOVBA4(0xF0 >> ((N - sh)& 3));
        AE_LA16X4_IP(x0, ax, castxcc(ae_int16x4, px));
        AE_LA16X4_IP(y0_0, ay0, castxcc(ae_int16x4, py0));
        AE_LA16X4_IP(y1_0, ay1, castxcc(ae_int16x4, py1));
        AE_LA16X4_IP(y2_0, ay2, castxcc(ae_int16x4, py2));
        AE_LA16X4_IP(y3_0, ay3, castxcc(ae_int16x4, py3));

        AE_MOVF16X4(y0_0, AE_ZERO16(), bmask);
        AE_MOVF16X4(y1_0, AE_ZERO16(), bmask);
        AE_MOVF16X4(y2_0, AE_ZERO16(), bmask);
        AE_MOVF16X4(y3_0, AE_ZERO16(), bmask);

        AE_MULAAAAQ16(q0, y0_0, x0);
        AE_MULAAAAQ16(q2, y1_0, x0);
        AE_MULAAAAQ16(q4, y2_0, x0);
        AE_MULAAAAQ16(q6, y3_0, x0);
      }
    }
    
#endif
    q0 = AE_ADD64(q0, q1);
    q1 = AE_ADD64(q2, q3);
    q2 = AE_ADD64(q4, q5);
    q3 = AE_ADD64(q6, q7);
    // Q(15 - rsh) <- [ Q16.47 - w/ rounding and saturation ]
    vai0 = AE_TRUNCA32X2F64S(q0, q1, rsh);
    vai1 = AE_TRUNCA32X2F64S(q2, q3, rsh);
    AE_SA32X2_IP(vai0, az, pz);
    AE_SA32X2_IP(vai1, az, pz);
  }
  AE_SA64POS_FP(az, pz);

  if (M & 3)
  {
    ae_valign   x_align, y0_align;
    for (m = (M&(~3)); m < M; m++)
    {
      px = (const ae_int16x8 *)x;
      py0 = (const ae_int16x8 *)y[m];

      x_align = AE_LA64_PP(px);
      y0_align = AE_LA64_PP(py0);

      q0 = AE_ZERO64();

      for (n = 0; n < N - 3; n += 4)
      {
        AE_LA16X4_IP(vxh, x_align, castxcc(ae_int16x4, px));
        AE_LA16X4_IP(vyh0, y0_align, castxcc(ae_int16x4, py0));

        // Q16.47 <- [ Q16*Q16 + 1] + [ Q16*Q16 + 1 ] + [ Q16*Q16 + 1] + [ Q16*Q16 + 1 ]
        AE_MULAAAAQ16(q0, vyh0, vxh);
      }
      // tail
      bmask = AE_MOVBA4(0xF0 >> (N & 3));

      AE_LA16X4_IP(vxh, x_align, castxcc(ae_int16x4, px));
      AE_LA16X4_IP(vyh0, y0_align, castxcc(ae_int16x4, py0));

      AE_MOVF16X4(vyh0, AE_ZERO16(), bmask);

      // Q16.47 <- [ Q16*Q16 + 1] + [ Q16*Q16 + 1 ] + [ Q16*Q16 + 1] + [ Q16*Q16 + 1 ]
      AE_MULAAAAQ16(q0, vyh0, vxh);

      // Q(15 - rsh) <- [ Q16.47 - w/ rounding and saturation ]
      q0 = AE_SLAA64S(q0, rsh);
      z[m] = AE_TRUNCA32Q64(q0);
    }
  }
} /* vec_dot_batch16x16() */
