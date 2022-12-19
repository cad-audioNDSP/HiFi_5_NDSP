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
    2D Convolution, 8x8-bit
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

#if defined(AE_MULA8Q8X8CNV_L) // NN extension
void conv2d_5x5_8x8(void* pScr, int8_t* z, const int8_t* x, const int8_t* y, int rsh, int P, int Q)
{
#define M 5
#define N 5
    const ae_int8x8* restrict pY;
    const ae_int8x8* restrict pY0;
    const ae_int32x4* restrict pT_read;
    ae_int8x8* restrict pW;
    ae_int32x4* restrict pT;
    ae_int8x8* restrict pZ0;
    ae_int32x2 S0, S1, S2, S3;
    ae_int16x4 T0, T1;
    ae_int8x8 R0;
    const ae_int8* pX = (const ae_int8*)x;

    ae_valign alZ;
    ae_int8x8 Y00, Y01;

    int i, j, mmin, mmax, xstart, ystart, m;

    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT_ALIGN(z, 16);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0 && P % 8 == 0);
    NASSERT(Q >= 0 && Q % 8 == 0);
    if (P <= 0 || Q <= 0) return;

    ae_int8x8 W0, W1, W2, W3, W4;

    alZ = AE_ZALIGN64();

    pW = (ae_int8x8*)pScr;
    // preload X and save in scratch
    {
        ae_int8x8 w0, w1, w2, w3, w4, w5;
        w0 = AE_MOVINT8X8_FROMINT16(0);
        AE_L8_IP(w5, pX, sizeof(ae_int8));
        AE_L8_IP(w4, pX, sizeof(ae_int8));
        AE_L8_IP(w3, pX, sizeof(ae_int8));
        AE_L8_IP(w2, pX, sizeof(ae_int8));
        AE_L8_IP(w1, pX, sizeof(ae_int8));

        W4 = w0;
        W4 = AE_SEL8X8I(W4, w1, 16);
        W4 = AE_SEL8X8I(W4, w2, 16);
        W4 = AE_SEL8X8I(W4, w3, 16);
        W4 = AE_SEL8X8I(W4, w4, 16);
        W4 = AE_SEL8X8I(W4, w5, 16);
        W4 = AE_SEL8X8I(W4, w0, 17);

        AE_L8_IP(w5, pX, sizeof(ae_int8));
        AE_L8_IP(w4, pX, sizeof(ae_int8));
        AE_L8_IP(w3, pX, sizeof(ae_int8));
        AE_L8_IP(w2, pX, sizeof(ae_int8));
        AE_L8_IP(w1, pX, sizeof(ae_int8));

        W3 = w0;
        W3 = AE_SEL8X8I(W3, w1, 16);
        W3 = AE_SEL8X8I(W3, w2, 16);
        W3 = AE_SEL8X8I(W3, w3, 16);
        W3 = AE_SEL8X8I(W3, w4, 16);
        W3 = AE_SEL8X8I(W3, w5, 16);
        W3 = AE_SEL8X8I(W3, w0, 17);

        AE_L8_IP(w5, pX, sizeof(ae_int8));
        AE_L8_IP(w4, pX, sizeof(ae_int8));
        AE_L8_IP(w3, pX, sizeof(ae_int8));
        AE_L8_IP(w2, pX, sizeof(ae_int8));
        AE_L8_IP(w1, pX, sizeof(ae_int8));

        W2 = w0;
        W2 = AE_SEL8X8I(W2, w1, 16);
        W2 = AE_SEL8X8I(W2, w2, 16);
        W2 = AE_SEL8X8I(W2, w3, 16);
        W2 = AE_SEL8X8I(W2, w4, 16);
        W2 = AE_SEL8X8I(W2, w5, 16);
        W2 = AE_SEL8X8I(W2, w0, 17);

        AE_L8_IP(w5, pX, sizeof(ae_int8));
        AE_L8_IP(w4, pX, sizeof(ae_int8));
        AE_L8_IP(w3, pX, sizeof(ae_int8));
        AE_L8_IP(w2, pX, sizeof(ae_int8));
        AE_L8_IP(w1, pX, sizeof(ae_int8));

        W1 = w0;
        W1 = AE_SEL8X8I(W1, w1, 16);
        W1 = AE_SEL8X8I(W1, w2, 16);
        W1 = AE_SEL8X8I(W1, w3, 16);
        W1 = AE_SEL8X8I(W1, w4, 16);
        W1 = AE_SEL8X8I(W1, w5, 16);
        W1 = AE_SEL8X8I(W1, w0, 17);

        AE_L8_IP(w5, pX, sizeof(ae_int8));
        AE_L8_IP(w4, pX, sizeof(ae_int8));
        AE_L8_IP(w3, pX, sizeof(ae_int8));
        AE_L8_IP(w2, pX, sizeof(ae_int8));
        AE_L8_IP(w1, pX, sizeof(ae_int8));

        W0 = w0;
        W0 = AE_SEL8X8I(W0, w1, 16);
        W0 = AE_SEL8X8I(W0, w2, 16);
        W0 = AE_SEL8X8I(W0, w3, 16);
        W0 = AE_SEL8X8I(W0, w4, 16);
        W0 = AE_SEL8X8I(W0, w5, 16);
        W0 = AE_SEL8X8I(W0, w0, 17);

        AE_S8X8_IP(W0, pW, sizeof(ae_int8x8));
        AE_S8X8_IP(W1, pW, sizeof(ae_int8x8));
        AE_S8X8_IP(W2, pW, sizeof(ae_int8x8));
        AE_S8X8_IP(W3, pW, sizeof(ae_int8x8));
        AE_S8X8_IP(W4, pW, sizeof(ae_int8x8));
    }

    for (i = 0; i < P + M - 1; i += 1)
    {
        mmin = i < M - 1 ? M - 1 - i : 0;
        mmax = i >= P ? M - 1 + P - i : M;
        xstart = i < M - 1 ? (mmin) : 0;
        ystart = i < M - 1 ? 0 : i - (M - 1);

        pW = (ae_int8x8*)pScr + xstart;
        pY0 = (ae_int8x8*)(y + ystart * Q);
        pZ0 = (ae_int8x8*)(z + i * (Q + N - 1));

        if (mmax - mmin != 1)
        {
            pY = (ae_int8x8*)((ae_int8*)pY0 + Q);
            // first iteration
            {
                pT_read = pT = (ae_int32x4*)pScr + 3;

                Y00 = AE_MOVINT8X8_FROMINT16(0);
                AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));

                for (j = 0; j < Q; j += 8)
                {
                    AE_L8X8_IP(Y01, pY0, sizeof(ae_int8x8));

                    AE_MUL8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                    AE_MUL8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                    Y00 = Y01;

                    AE_S32X2X2_IP(S0, S1, pT, sizeof(ae_int32x4));
                    AE_S32X2X2_IP(S2, S3, pT, sizeof(ae_int32x4));
                }
                //last N-1 elements
                {
                    Y01 = AE_MOVINT8X8_FROMINT16(0);

                    AE_MUL8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                    AE_MUL8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                    AE_S32X2X2_IP(S0, S1, pT, sizeof(ae_int32x4));
                    AE_S32X2X2_IP(S2, S3, pT, sizeof(ae_int32x4));
                }
            }
            pY0 = pY;
            pY = (ae_int8x8*)((ae_int8*)pY + Q);
            __Pragma("loop_count min=0, max=2");
            for (m = mmin + 1; m < mmax - 1; ++m)
            {
                pT_read = pT = (ae_int32x4*)pScr + 3;

                Y00 = AE_MOVINT8X8_FROMINT16(0);
                AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));

                for (j = 0; j < Q; j += 8)
                {
                    AE_L32X2X2_IP(S0, S1, pT_read, sizeof(ae_int32x4));
                    AE_L32X2X2_IP(S2, S3, pT_read, sizeof(ae_int32x4));

                    AE_L8X8_IP(Y01, pY0, sizeof(ae_int8x8));

                    AE_MULA8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                    AE_MULA8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                    Y00 = Y01;

                    AE_S32X2X2_IP(S0, S1, pT, sizeof(ae_int32x4));
                    AE_S32X2X2_IP(S2, S3, pT, sizeof(ae_int32x4));

                }
                //last N-1 elements
                {
                    AE_L32X2X2_IP(S0, S1, pT_read, sizeof(ae_int32x4));
                    AE_L32X2X2_IP(S2, S3, pT_read, sizeof(ae_int32x4));

                    Y01 = AE_MOVINT8X8_FROMINT16(0);

                    AE_MULA8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                    AE_MULA8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                    AE_S32X2X2_IP(S0, S1, pT, sizeof(ae_int32x4));
                    AE_S32X2X2_IP(S2, S3, pT, sizeof(ae_int32x4));
                }
                pY0 = pY;
                pY = (ae_int8x8*)((ae_int8*)pY + Q);
            }
            //last iteration
            {
                pT_read = pT = (ae_int32x4*)pScr + 3;

                Y00 = AE_MOVINT8X8_FROMINT16(0);

                AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));

                for (j = 0; j < Q; j += 8)
                {
                    AE_L32X2X2_IP(S0, S1, pT_read, sizeof(ae_int32x4));
                    AE_L32X2X2_IP(S2, S3, pT_read, sizeof(ae_int32x4));

                    AE_L8X8_IP(Y01, pY0, sizeof(ae_int8x8));
                    AE_MULA8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                    AE_MULA8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                    Y00 = Y01;

                    T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                    AE_SA8X8_IP(AE_ROUND8X8F16SSYM(T0, T1), alZ, pZ0);
                }
                AE_SA64POS_FP(alZ, pZ0);
                //last N-1 elements
                {
                    AE_L32X2X2_IP(S0, S1, pT_read, sizeof(ae_int32x4));
                    AE_L32X2X2_IP(S2, S3, pT_read, sizeof(ae_int32x4));

                    Y01 = AE_MOVINT8X8_FROMINT16(0);

                    AE_MULA8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                    AE_MULA8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                    T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                    R0 = AE_ROUND8X8F16SSYM(T0, T1);

                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
        }
        else
        {
            Y00 = AE_MOVINT8X8_FROMINT16(0);

            AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));

            for (j = 0; j < Q; j += 8)
            {
                AE_L8X8_IP(Y01, pY0, sizeof(ae_int8x8));
                AE_MUL8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                AE_MUL8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                Y00 = Y01;

                T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                AE_SA8X8_IP(AE_ROUND8X8F16SSYM(T0, T1), alZ, pZ0);
            }
            AE_SA64POS_FP(alZ, pZ0);
            //last N-1 elements
            {
                Y01 = AE_MOVINT8X8_FROMINT16(0);

                AE_MUL8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                AE_MUL8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                R0 = AE_ROUND8X8F16SSYM(T0, T1);

                R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
            }
        }
    }
#undef M
#undef N
}

size_t conv2d_5x5_8x8_getScratchSize(int P, int Q)
{
    return 6 * 8 * sizeof(int8_t) + (Q + 16) * sizeof(int32_t);
} // MxN=5x5

#else // no NN extension

void conv2d_5x5_8x8(void * pScr, int8_t * z, const int8_t * x, const int8_t * y, int rsh, int P, int Q)
{
#define M 5
#define N 5
    int16_t* w;// [9][2*8]: weights - reordered x
    const int8_t     * restrict pY0;
    const int8_t     * restrict pY1;
    const int8_t     * restrict pY2;
    const int8_t     * restrict pY3;
    const int8_t     * restrict pY4;
          ae_int32   * restrict pZ = (ae_int32 *)z;
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
        int8_t * pX;
        ae_valign aX;
        ae_int16x4 x0, x1;
        ae_int8x8 x0_;
        /* pad coefficients with zeroes for processing N-1 first and last rows */
        for (m = 0; m < (M - 1); m++)
        {
            AE_S16X4X2_IP(AE_ZERO16(), AE_ZERO16(), castxcc(ae_int16x8, pW), 8 * sizeof(int16_t));
            AE_S16X4X2_IP(AE_ZERO16(), AE_ZERO16(), castxcc(ae_int16x8, pW), 8 * sizeof(int16_t));
        }
        pW = (ae_int16x4 *)XT_ADDX2(9 * 8, (uintptr_t)pW);
        pX = (int8_t *)x;
        for (m = 0; m < M; m++)
        {
            aX = AE_LA64_PP(pX);
            AE_LA8X4U_IP(x1, aX, pX);
            x1 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(AE_MOVINT8X8_FROMINT16X4(x1), AE_MOVINT8X8_FROMINT16X4(AE_ZERO16()), 30));
            AE_L8_IP(x0_, castxcc(ae_int8, pX), sizeof(int8_t));
            x0 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(x0_, AE_MOVINT8X8_FROMINT16X4(AE_ZERO16()), 21));
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
        ae_int32x2 r;

        yrow = XT_MAX(XT_MIN(i, P - 1) - (M - 1), 0);
        wrow = 2 * (M - 1) - i + yrow;
        NASSERT(yrow >= 0 && yrow <= (P - M));
        NASSERT(wrow >= 0 && wrow <= 8);
        pY0 = y + yrow*Q;
        pY1 = pY0 + Q;
        pY2 = pY1 + Q;
        pY3 = pY2 + Q;
        pY4 = pY3 + Q;

        /* First N-1 samples */
        {
            pW = (ae_int16x4*)(w + 2 * 8 * wrow) + 3;
            AE_L16X4_XC(w01, pW, 4 * 4 * sizeof(int16_t));
            AE_L16X4_XC(w11, pW, 4 * 4 * sizeof(int16_t));
            AE_L16X4_XC(w21, pW, 4 * 4 * sizeof(int16_t));
            AE_L16X4_XC(w31, pW, 4 * 4 * sizeof(int16_t));
            AE_L16X4_IP(w41, pW, 4 * 4 * sizeof(int16_t));

            AE_L8X4F_IP(y01, pY0, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y11, pY1, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y21, pY2, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y31, pY3, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y41, pY4, 4 * sizeof(int8_t));

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

            r = AE_MOVINT32X2_FROMINT8X8(AE_ROUND8X4F32SSYM_L(AE_TRUNCA32X2F64S(S3, S2, 32 - rsh), AE_TRUNCA32X2F64S(S1, S0, 32 - rsh)));
            AE_S32_L_IP(r, pZ, 4 * sizeof(int8_t));

            y00 = y01;
            y10 = y11;
            y20 = y21;
            y30 = y31;
            y40 = y41;
        }
        /* Next samples */
        pW = (ae_int16x4*)(w + 2 * 8 * wrow);
        AE_L16X4_IP(w00, pW, 4 * sizeof(int16_t)); AE_L16X4_XC(w01, pW, 3 * 4 * sizeof(int16_t));
        AE_L16X4_IP(w10, pW, 4 * sizeof(int16_t)); AE_L16X4_XC(w11, pW, 3 * 4 * sizeof(int16_t));
        AE_L16X4_IP(w20, pW, 4 * sizeof(int16_t)); AE_L16X4_XC(w21, pW, 3 * 4 * sizeof(int16_t));
        AE_L16X4_IP(w30, pW, 4 * sizeof(int16_t)); AE_L16X4_XC(w31, pW, 3 * 4 * sizeof(int16_t));
        AE_L16X4_IP(w40, pW, 4 * sizeof(int16_t)); AE_L16X4_XC(w41, pW, 3 * 4 * sizeof(int16_t));
        __Pragma("loop_count min=1");
        for (j = 0; j < ((Q + N - 1) >> 2) - 2; j++)
        {
            AE_L8X4F_IP(y01, pY0, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y11, pY1, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y21, pY2, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y31, pY3, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y41, pY4, 4 * sizeof(int8_t));

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

            r = AE_MOVINT32X2_FROMINT8X8(AE_ROUND8X4F32SSYM_L(AE_TRUNCA32X2F64S(S3, S2, 32 - rsh), AE_TRUNCA32X2F64S(S1, S0, 32 - rsh)));
            AE_S32_L_IP(r, pZ, 4 * sizeof(int8_t));

            y00 = AE_L8X4F_I(pY0, -4 * (int)sizeof(int8_t));
            y10 = AE_L8X4F_I(pY1, -4 * (int)sizeof(int8_t));
            y20 = AE_L8X4F_I(pY2, -4 * (int)sizeof(int8_t));
            y30 = AE_L8X4F_I(pY3, -4 * (int)sizeof(int8_t));
            y40 = AE_L8X4F_I(pY4, -4 * (int)sizeof(int8_t));
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

            r = AE_MOVINT32X2_FROMINT8X8(AE_ROUND8X4F32SSYM_L(AE_TRUNCA32X2F64S(S3, S2, 32 - rsh), AE_TRUNCA32X2F64S(S1, S0, 32 - rsh)));
            AE_S32_L_IP(r, pZ, 4 * sizeof(int8_t));
        }
    }
#undef M
#undef N
}

// scratch allocatation functions. return required scratch size in bytes
size_t conv2d_5x5_8x8_getScratchSize(int P, int Q)
{
    return 2*9*8*sizeof(int16_t);/* weights (reordered x[5][5]) */
} // MxN=5x5

#endif
