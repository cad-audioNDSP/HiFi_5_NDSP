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
  NatureDSP Signal Processing Library. FIR part
    2D Convolution, 16x16-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

/*-------------------------------------------------------------------------
  2D convolution
  Functions compute the two-dimensional convolution of input matrix x[M][N]
  and y[P][Q] and store the result in matrix z[M+P-1][N+Q-1]
  Additional parameter rsh allows to contro l fixed point representation of 
  output data.
  Two versions of functions available: 
  - generic version with _gen_ suffix. 
    These functions work with arbitrary arguments.
  - fast version with no _gen_ suffix. 
    These functions expose some additional restrictions on argument


  Precision: 
  8x8      8-bit coefficients, 8-bit data, 8-bit output, Q7
  8x16     8-bit coefficients Q7, 16-bit data, 16-bit output, Q15
  16x16    16-bit coefficients, 16-bit data, 16-bit output, Q15
  f        single precision floating point data
  fp16     half precision floating point data

  Input:
  x[M][N]   input data Q15, Q7, floating point
  y[P][Q]   input data Q15, Q7, floating point
  M         number of rows in the matrix x
  N         number of columns in the matrix x
  P         number of rows in the matrix y
  Q         number of columns in the matrix y
  rsh       additional right shift (for fixed point API only)

  Output:
  z	[M+P-1][N+Q-1] output data, Q(7-rsh), Q(15-rsh)

  Temporary:
  pScr     scratch data. Should have size at least as requested by 
           corresponding scratch allocation function

  Restrictions:
  For regular routines:
  x,y,z        should not overlap
  pScr         aligned on a 16-bytes boundary
  P, Q	       >0

  For fast routines:
  x,y,z        should not overlap
  x,y,z,pScr   aligned on a 16-bytes boundary
  P, Q	       >0 and multiplies of 8
-------------------------------------------------------------------------*/

void conv2d_5x5_16x16(void * pScr, int16_t * z, const int16_t * x, const int16_t * y, int rsh, int P, int Q)
{
#define M 5
#define N 5
    int16_t* w;// [9][2*8]: weights - reordered x
    const ae_int16x4 * restrict pY0;
    const ae_int16x4 * restrict pY1;
    const ae_int16x4 * restrict pY2;
    const ae_int16x4 * restrict pY3;
    const ae_int16x4 * restrict pY4;
          ae_int16x4 * restrict pZ = (ae_int16x4 *)z;
          ae_int16x4 * restrict pW;
    int i, j;

    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT_ALIGN(z, 16);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0 && P % 8 == 0);
    NASSERT(Q >= 0 && Q % 8 == 0);
    if (P <= 0 || Q <= 0) return;

    w = (int16_t*)pScr;

    /* fill coefficients in right order
    order:
    x[4] x[3] x[2] x[1] x[0] 0    0    0
    0    0    0    0    x[3] x[2] x[1] x[0]
    .... next 5 coefficients ....
    */
    pW = (ae_int16x4*)w;
    {
        int m;
        ae_int16x4 * pX;
        ae_valign aX;
        ae_int16x4 x0, x1;
        /* pad coefficients with zeroes for processing N-1 first and last rows */
        for (m = 0; m < (M - 1); m++)
        {
            AE_S16X4X2_IP(AE_ZERO16(), AE_ZERO16(), castxcc(ae_int16x8, pW), 8 * sizeof(int16_t));
            AE_S16X4X2_IP(AE_ZERO16(), AE_ZERO16(), castxcc(ae_int16x8, pW), 8 * sizeof(int16_t));
        }
        pW = (ae_int16x4 *)XT_ADDX2(9 * 8, (uintptr_t)pW);
        pX = (ae_int16x4 *)x;
        for (m = 0; m < M; m++)
        {
            aX = AE_LA64_PP(pX);
            AE_LA16X4_IP(x1, aX, pX);
            AE_L16_IP(x0, castxcc(ae_int16, pX), sizeof(int16_t));
            AE_S64X2_IP(AE_ZERO64(), AE_MOVINT64_FROMINT16X4(x1), castxcc(ae_int64x2, pW), -8 * (int)sizeof(int16_t));
            x0 = AE_SEL16_6543(x1, x0);
            x1 = AE_SEL16_6543(AE_ZERO16(), x1);
            AE_S64X2_IP(AE_MOVINT64_FROMINT16X4(x0), AE_MOVINT64_FROMINT16X4(x1), castxcc(ae_int64x2, pW), -8 * (int)sizeof(int16_t));
        }
    }
    WAE_CBEGIN0((uintptr_t)(w));
    WAE_CEND0  ((uintptr_t)(w + (M + M - 1) * 2 * 8));

    for (i = 0; i < (P + M - 1); i++)
    {
        int wrow, yrow;
        ae_int16x4 y00, y01, y10, y11, y20, y21, y30, y31, y40, y41;
        ae_int16x4 w00, w01, w10, w11, w20, w21, w30, w31, w40, w41;
        ae_int64 S0, S1, S2, S3;
        ae_int16x4 r;

        yrow = XT_MAX(XT_MIN(i, P - 1) - (M - 1), 0);
        wrow = 2 * (M - 1) - i + yrow;
        NASSERT(yrow >= 0 && yrow <= (P - M));
        NASSERT(wrow >= 0 && wrow <= 8);
        pY0 = (ae_int16x4 *)(y + yrow*Q);
        pY1 = (ae_int16x4 *)((int16_t *)pY0 + Q);
        pY2 = (ae_int16x4 *)((int16_t *)pY1 + Q);
        pY3 = (ae_int16x4 *)((int16_t *)pY2 + Q);
        pY4 = (ae_int16x4 *)((int16_t *)pY3 + Q);

        /* First N-1 samples */
        {
            pW = (ae_int16x4*)(w + 2 * 8 * wrow) + 3;
            AE_L16X4_XC(w01, pW, 4 * 4 * sizeof(int16_t));
            AE_L16X4_XC(w11, pW, 4 * 4 * sizeof(int16_t));
            AE_L16X4_XC(w21, pW, 4 * 4 * sizeof(int16_t));
            AE_L16X4_XC(w31, pW, 4 * 4 * sizeof(int16_t));
            AE_L16X4_IP(w41, pW, 4 * 4 * sizeof(int16_t));

            AE_L16X4_IP(y01, pY0, 4 * sizeof(int16_t));
            AE_L16X4_IP(y11, pY1, 4 * sizeof(int16_t));
            AE_L16X4_IP(y21, pY2, 4 * sizeof(int16_t));
            AE_L16X4_IP(y31, pY3, 4 * sizeof(int16_t));
            AE_L16X4_IP(y41, pY4, 4 * sizeof(int16_t));

            AE_MULFQ16X2_FIR_2(S0, S1, AE_ZERO16(), y01, w01);
            AE_MULFQ16X2_FIR_0(S2, S3, AE_ZERO16(), y01, w01);
            AE_MULAFQ16X2_FIR_2(S0, S1, AE_ZERO16(), y11, w11);
            AE_MULAFQ16X2_FIR_0(S2, S3, AE_ZERO16(), y11, w11);
            AE_MULAFQ16X2_FIR_2(S0, S1, AE_ZERO16(), y21, w21);
            AE_MULAFQ16X2_FIR_0(S2, S3, AE_ZERO16(), y21, w21);
            AE_MULAFQ16X2_FIR_2(S0, S1, AE_ZERO16(), y31, w31);
            AE_MULAFQ16X2_FIR_0(S2, S3, AE_ZERO16(), y31, w31);
            AE_MULAFQ16X2_FIR_2(S0, S1, AE_ZERO16(), y41, w41);
            AE_MULAFQ16X2_FIR_0(S2, S3, AE_ZERO16(), y41, w41);

            r = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
            AE_S16X4_IP(r, pZ, 4 * sizeof(int16_t));

            y00 = y01;
            y10 = y11;
            y20 = y21;
            y30 = y31;
            y40 = y41;
        }
        /* Next samples */
        pW = (ae_int16x4*)(w + 2 * 8 * wrow);
        AE_L16X4X2_XC(w00, w01, castxcc(ae_int16x8, pW), 4 * 4 * sizeof(int16_t));
        AE_L16X4X2_XC(w10, w11, castxcc(ae_int16x8, pW), 4 * 4 * sizeof(int16_t));
        AE_L16X4X2_XC(w20, w21, castxcc(ae_int16x8, pW), 4 * 4 * sizeof(int16_t));
        AE_L16X4X2_XC(w30, w31, castxcc(ae_int16x8, pW), 4 * 4 * sizeof(int16_t));
        AE_L16X4X2_IP(w40, w41, castxcc(ae_int16x8, pW), 4 * 4 * sizeof(int16_t));
        __Pragma("loop_count min=1");
        for (j = 0; j < ((Q + N - 1) >> 2) - 2; j++)
        {
            AE_L16X4_IP(y01, pY0, 4 * sizeof(int16_t));
            AE_L16X4_IP(y11, pY1, 4 * sizeof(int16_t));
            AE_L16X4_IP(y21, pY2, 4 * sizeof(int16_t));
            AE_L16X4_IP(y31, pY3, 4 * sizeof(int16_t));
            AE_L16X4_IP(y41, pY4, 4 * sizeof(int16_t));

            AE_MULFQ16X2_FIR_3(S0, S1, y00, y01, w00);
            AE_MULFQ16X2_FIR_1(S2, S3, y00, y01, w00);
            AE_MULAFQ16X2_FIR_3(S0, S1, y01, AE_ZERO16(), w01);
            AE_MULAFQ16X2_FIR_1(S2, S3, y01, AE_ZERO16(), w01);
            AE_MULAFQ16X2_FIR_3(S0, S1, y10, y11, w10);
            AE_MULAFQ16X2_FIR_1(S2, S3, y10, y11, w10);
            AE_MULAFQ16X2_FIR_3(S0, S1, y11, AE_ZERO16(), w11);
            AE_MULAFQ16X2_FIR_1(S2, S3, y11, AE_ZERO16(), w11);
            AE_MULAFQ16X2_FIR_3(S0, S1, y20, y21, w20);
            AE_MULAFQ16X2_FIR_1(S2, S3, y20, y21, w20);
            AE_MULAFQ16X2_FIR_3(S0, S1, y21, AE_ZERO16(), w21);
            AE_MULAFQ16X2_FIR_1(S2, S3, y21, AE_ZERO16(), w21);
            AE_MULAFQ16X2_FIR_3(S0, S1, y30, y31, w30);
            AE_MULAFQ16X2_FIR_1(S2, S3, y30, y31, w30);
            AE_MULAFQ16X2_FIR_3(S0, S1, y31, AE_ZERO16(), w31);
            AE_MULAFQ16X2_FIR_1(S2, S3, y31, AE_ZERO16(), w31);
            AE_MULAFQ16X2_FIR_3(S0, S1, y40, y41, w40);
            AE_MULAFQ16X2_FIR_1(S2, S3, y40, y41, w40);
            AE_MULAFQ16X2_FIR_3(S0, S1, y41, AE_ZERO16(), w41);
            AE_MULAFQ16X2_FIR_1(S2, S3, y41, AE_ZERO16(), w41);

            r = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
            AE_S16X4_IP(r, pZ, 4 * sizeof(int16_t));

            y00 = AE_L16X4_I(pY0, -4 * (int)sizeof(int16_t));
            y10 = AE_L16X4_I(pY1, -4 * (int)sizeof(int16_t));
            y20 = AE_L16X4_I(pY2, -4 * (int)sizeof(int16_t));
            y30 = AE_L16X4_I(pY3, -4 * (int)sizeof(int16_t));
            y40 = AE_L16X4_I(pY4, -4 * (int)sizeof(int16_t));
        }
        /* Last N-1 samples */
        {
            AE_MULFQ16X2_FIR_3(S0, S1, y00, AE_ZERO16(), w00);
            AE_MULFQ16X2_FIR_1(S2, S3, y00, AE_ZERO16(), w00);
            AE_MULAFQ16X2_FIR_3(S0, S1, y10, AE_ZERO16(), w10);
            AE_MULAFQ16X2_FIR_1(S2, S3, y10, AE_ZERO16(), w10);
            AE_MULAFQ16X2_FIR_3(S0, S1, y20, AE_ZERO16(), w20);
            AE_MULAFQ16X2_FIR_1(S2, S3, y20, AE_ZERO16(), w20);
            AE_MULAFQ16X2_FIR_3(S0, S1, y30, AE_ZERO16(), w30);
            AE_MULAFQ16X2_FIR_1(S2, S3, y30, AE_ZERO16(), w30);
            AE_MULAFQ16X2_FIR_3(S0, S1, y40, AE_ZERO16(), w40);
            AE_MULAFQ16X2_FIR_1(S2, S3, y40, AE_ZERO16(), w40);

            r = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
            AE_S16X4_IP(r, pZ, 4 * sizeof(int16_t));
        }
    }
#undef M
#undef N
}

// scratch allocatation functions. return required scratch size in bytes
size_t conv2d_5x5_16x16_getScratchSize(int P, int Q)
{
    return 2*9*8*sizeof(int16_t);/* weights (reordered x[5][5]) */
} // MxN=5x5
