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
   Matrix transpose
   Code optimized for HiFi5 core 
   */
/* Code optimized for HiFi5 core */
/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Library API */
#include "NatureDSP_Signal_matop.h"
/* Common helper macros. */
#include "common.h"
#include <stdio.h>

/*-------------------------------------------------------------------------
  Matrix Transpose
  These functions transpose matrices.

  Precision: 
  32x32 32-bit input, 32-bit output
  16x16 16-bit input, 16-bit output
  8x8   8-bit inputs, 8-bit output
  f     floating point

  Input:
  x[M][N] input matrix,Q31,Q15,Q7 or floating point
  M       number of rows in matrix x
  N       number of columns in matrix x
  Output:
  y[N][M] output vector,Q31,Q15,Q7 or floating point

  Restriction:
  For regular routines (mtx_transpose_32x32, mtx_transpose_16x16, 
  mtx_transpose_8x8, mtx_transposef):
  x,y should not overlap

  For faster routines (mtx_transpose 32x32_fast, mtx_transpose 16x16_fast, 
  mtx_transpose_8x8_fast, mtx_transposef_fast)
  x,y   should not overlap
  x,y   aligned on 16-byte boundary
  N and M are multiples of 4
-------------------------------------------------------------------------*/
void mtx_transpose32x32_fast(int32_t  *  y, const int32_t*     x, int M, int N)
{
    const ae_int32x4 * restrict pX;
          ae_int32x4 * restrict pY;
    int m, n;
    NASSERT(x);
    NASSERT(y);
    NASSERT(x != y);
    NASSERT(N % 4 == 0 && M % 4 == 0);
    if (M <= 0 || N <= 0) return;
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    pX = (const ae_int32x4*)(x);

    for (m = 0; m < (M>>2); m++)
    {
        pY = (ae_int32x4*)(y);
        __Pragma("loop_count min=1")
        for (n = 0; n < (N>>2); n++)
        {
            ae_int32x2 x0,x1,x2,x3,x4,x5,x6,x7;
            ae_int32x2 y0,y1,y2,y3,y4,y5,y6,y7;
            AE_L32X2X2_X (x2, x3, pX, 1 * N*sizeof(int32_t));
            AE_L32X2X2_X (x4, x5, pX, 2 * N*sizeof(int32_t));
            AE_L32X2X2_X (x6, x7, pX, 3 * N*sizeof(int32_t));
            AE_L32X2X2_IP(x0, x1, pX, sizeof(ae_int32x4));
            y0 = AE_SEL32_HH(x0,x2);
            y1 = AE_SEL32_LL(x0,x2);
            y2 = AE_SEL32_HH(x4,x6);
            y3 = AE_SEL32_LL(x4,x6);
            y4 = AE_SEL32_HH(x1,x3);
            y5 = AE_SEL32_LL(x1,x3);
            y6 = AE_SEL32_HH(x5,x7);
            y7 = AE_SEL32_LL(x5,x7);
            AE_S32X2X2_XP(y0,y2, pY, M*sizeof(int32_t));
            AE_S32X2X2_XP(y1,y3, pY, M*sizeof(int32_t));
            AE_S32X2X2_XP(y4,y6, pY, M*sizeof(int32_t));
            AE_S32X2X2_XP(y5,y7, pY, M*sizeof(int32_t));
        }
        pX += 3 * (N >> 2);
        y  +=4;
    }
}
