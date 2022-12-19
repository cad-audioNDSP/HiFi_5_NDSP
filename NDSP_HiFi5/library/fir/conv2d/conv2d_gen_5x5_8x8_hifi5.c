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
void conv2d_gen_5x5_8x8(void * pScr, int8_t * z, const int8_t * x, const int8_t * y, int rsh, int P, int Q)
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

    ae_valign al0, alZ;
    ae_int8x8 Y00, Y01;

    int i, j, mmin, mmax, xstart, ystart, m;

    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0);
    NASSERT(Q >= 0);
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

    __Pragma("loop_count min=1");
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
                al0 = AE_LA64_PP(pY0);
                AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));

                for (j = 0; j < (Q & ~7); j += 8)
                {
                    AE_LA8X8_IP(Y01, al0, pY0);

                    AE_MUL8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                    AE_MUL8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                    Y00 = Y01;

                    AE_S32X2X2_IP(S0, S1, pT, sizeof(ae_int32x4));
                    AE_S32X2X2_IP(S2, S3, pT, sizeof(ae_int32x4));

                }
                AE_SA64POS_FP(alZ, pZ0);
                // tail
                if (j < Q)
                {

                    AE_LA8X8_IP(Y01, al0, pY0);

                    AE_MUL8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                    AE_MUL8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                    AE_S32X2X2_IP(S0, S1, pT, sizeof(ae_int32x4));
                    AE_S32X2X2_IP(S2, S3, pT, sizeof(ae_int32x4));

                    for (; j < Q; ++j)
                    {
                        AE_DSEL8X8(Y00, Y01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xE6D5C4B3A2918070));
                    }
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
                al0 = AE_LA64_PP(pY0);
                AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));

                for (j = 0; j < (Q & ~7); j += 8)
                {
                    AE_L32X2X2_IP(S0, S1, pT_read, sizeof(ae_int32x4));
                    AE_L32X2X2_IP(S2, S3, pT_read, sizeof(ae_int32x4));

                    AE_LA8X8_IP(Y01, al0, pY0);

                    AE_MULA8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                    AE_MULA8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                    Y00 = Y01;

                    AE_S32X2X2_IP(S0, S1, pT, sizeof(ae_int32x4));
                    AE_S32X2X2_IP(S2, S3, pT, sizeof(ae_int32x4));

                }
                AE_SA64POS_FP(alZ, pZ0);
                // tail
                if (j < Q)
                {
                    AE_L32X2X2_IP(S0, S1, pT_read, sizeof(ae_int32x4));
                    AE_L32X2X2_IP(S2, S3, pT_read, sizeof(ae_int32x4));

                    AE_LA8X8_IP(Y01, al0, pY0);

                    AE_MULA8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                    AE_MULA8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                    AE_S32X2X2_IP(S0, S1, pT, sizeof(ae_int32x4));
                    AE_S32X2X2_IP(S2, S3, pT, sizeof(ae_int32x4));

                    for (; j < Q; ++j)
                    {
                        AE_DSEL8X8(Y00, Y01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xE6D5C4B3A2918070));
                    }
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

                al0 = AE_LA64_PP(pY0);

                AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));

                for (j = 0; j < (Q & ~7); j += 8)
                {
                    AE_L32X2X2_IP(S0, S1, pT_read, sizeof(ae_int32x4));
                    AE_L32X2X2_IP(S2, S3, pT_read, sizeof(ae_int32x4));

                    AE_LA8X8_IP(Y01, al0, pY0);
                    AE_MULA8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                    AE_MULA8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                    Y00 = Y01;

                    T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                    AE_SA8X8_IP(AE_ROUND8X8F16SSYM(T0, T1), alZ, pZ0);
                }
                AE_SA64POS_FP(alZ, pZ0);
                // tail
                if (j < Q)
                {
                    AE_L32X2X2_IP(S0, S1, pT_read, sizeof(ae_int32x4));
                    AE_L32X2X2_IP(S2, S3, pT_read, sizeof(ae_int32x4));

                    AE_LA8X8_IP(Y01, al0, pY0);
                    AE_MULA8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                    AE_MULA8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                    T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                    R0 = AE_ROUND8X8F16SSYM(T0, T1);

                    for (; j < Q; ++j)
                    {
                        R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                        AE_DSEL8X8(Y00, Y01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xE6D5C4B3A2918070));
                        AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                    }
                }
                //last N-1 elements
                {
                    AE_L32X2X2_IP(S0, S1, pT_read, sizeof(ae_int32x4));
                    AE_L32X2X2_IP(S2, S3, pT_read, sizeof(ae_int32x4));

                    Y01 = AE_MOVINT8X8_FROMINT16(0);

                    AE_MULA8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                    AE_MULA8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                    T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                    R0 = AE_ROUND8X8F16SSYM(T0, T1);

                    for (j = Q; j < Q + N - 1; ++j)
                    {
                        R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                        AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                    }
                }
            }
        }
        else
        {
            Y00 = AE_MOVINT8X8_FROMINT16(0);

            al0 = AE_LA64_PP(pY0);

            AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));

            for (j = 0; j < (Q & ~7); j += 8)
            {
                AE_LA8X8_IP(Y01, al0, pY0);
                AE_MUL8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                AE_MUL8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                Y00 = Y01;

                T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                AE_SA8X8_IP(AE_ROUND8X8F16SSYM(T0, T1), alZ, pZ0);
            }
            AE_SA64POS_FP(alZ, pZ0);
            // tail
            if (j < Q)
            {
                AE_LA8X8_IP(Y01, al0, pY0);
                AE_MUL8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                AE_MUL8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                R0 = AE_ROUND8X8F16SSYM(T0, T1);

                for (; j < Q; ++j)
                {
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_DSEL8X8(Y00, Y01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xE6D5C4B3A2918070));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
            //last N-1 elements
            {
                Y01 = AE_MOVINT8X8_FROMINT16(0);

                AE_MUL8Q8X8CNV_L(S0, S1, W0, Y00, Y01);
                AE_MUL8Q8X8CNV_H(S2, S3, W0, Y01, Y01);

                T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                R0 = AE_ROUND8X8F16SSYM(T0, T1);

                for (j = Q; j < Q + N - 1; ++j)
                {
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
        }
    }
#undef M
#undef N
}

size_t conv2d_gen_5x5_8x8_getScratchSize(int P, int Q)
{
    return 6*8*sizeof(int8_t)+(Q+16)*sizeof(int32_t);
} // MxN=5x5

#else // no NN extension


void conv2d_gen_5x5_8x8(void* pScr, int8_t* z, const int8_t* x, const int8_t* y, int rsh, int P, int Q)
{
#define M 5
#define N 5
    const int8_t* restrict pY;
    const int8_t* restrict pY0;
    const ae_int64x2* restrict pT_read;
    int8_t* restrict pW;
    ae_int64x2* restrict pT;
    ae_int8x8* restrict pZ0;
    ae_int64 S0, S1, S2, S3, S4, S5, S6, S7;
    ae_int16x4 T0, T1;
    ae_int8x8 R0;
    const ae_int8* pX = (const ae_int8*)x;

    ae_valign al0, alZ;
    ae_int16x4 Y00, Y01, Y02;

    int i, j, mmin, mmax, xstart, ystart, m;

    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0);
    NASSERT(Q >= 0);
    if (P <= 0 || Q <= 0) return;

    
    ae_int16x4 W00, W01;

    alZ = AE_ZALIGN64();

    pW = (int8_t*)pScr;
    // preload X and save in scratch
    {
        ae_int8x8 W0, W1, W2, W3, W4;
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

        AE_S8X8_IP(W0, castxcc(ae_int8x8,pW), sizeof(ae_int8x8));
        AE_S8X8_IP(W1, castxcc(ae_int8x8,pW), sizeof(ae_int8x8));
        AE_S8X8_IP(W2, castxcc(ae_int8x8,pW), sizeof(ae_int8x8));
        AE_S8X8_IP(W3, castxcc(ae_int8x8,pW), sizeof(ae_int8x8));
        AE_S8X8_IP(W4, castxcc(ae_int8x8,pW), sizeof(ae_int8x8));
    }

    __Pragma("loop_count min=1");
    for (i = 0; i < P + M - 1; i += 1)
    {
        mmin = i < M - 1 ? M - 1 - i : 0;
        mmax = i >= P ? M - 1 + P - i : M;
        xstart = i < M - 1 ? (mmin) : 0;
        ystart = i < M - 1 ? 0 : i - (M - 1);

        pW = (int8_t*)pScr + 8*xstart;
        pY0 = (int8_t*)(y + ystart * Q);
        pZ0 = (ae_int8x8*)(z + i * (Q + N - 1));

        if (mmax - mmin != 1)
        {
            pY = (int8_t*)(pY0 + Q);
            // first iteration
            {
                pT = (ae_int64x2*)pScr + 3;

                Y00 = 0;
                al0 = AE_LA64_PP(pY0);

                AE_L8X4S_IP(W00, pW, 4 * sizeof(ae_int8)); AE_L8X4S_IP(W01, pW, 4 * sizeof(ae_int8));
                for (j = 0; j < (Q & ~7); j += 8)
                {
                    AE_LA8X4S_IP(Y01, al0, pY0); AE_LA8X4S_IP(Y02, al0, pY0);

                    AE_MULFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                    Y00 = Y02;

                    AE_S64X2_IP(S0, S1, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S2, S3, pT, sizeof(ae_int64x2));
                    AE_S64X2_IP(S4, S5, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S6, S7, pT, sizeof(ae_int64x2));

                }
                // tail
                if (j < Q)
                {
                    AE_LA8X4S_IP(Y01, al0, pY0); AE_LA8X4S_IP(Y02, al0, pY0);

                    AE_MULFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                    AE_S64X2_IP(S0, S1, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S2, S3, pT, sizeof(ae_int64x2));
                    AE_S64X2_IP(S4, S5, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S6, S7, pT, sizeof(ae_int64x2));
                    for (; j < Q; ++j)
                    {
                        AE_DSEL16X4(Y00, Y01, Y00, Y01, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                        AE_DSEL16X4(Y01, Y02, Y01, Y02, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    }
                }
                //last N-1 elements
                {
                    Y01 = 0; Y02 = 0;

                    AE_MULFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                    AE_S64X2_IP(S0, S1, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S2, S3, pT, sizeof(ae_int64x2));
                    AE_S64X2_IP(S4, S5, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S6, S7, pT, sizeof(ae_int64x2));
                }
            }
            pY0 = pY;
            pY = (int8_t*)(pY + Q);
            __Pragma("loop_count min=0, max=2");
            for (m = mmin + 1; m < mmax - 1; ++m)
            {
                pT_read = pT = (ae_int64x2*)pScr + 3;

                Y00 = 0;
                al0 = AE_LA64_PP(pY0);

                AE_L8X4S_IP(W00, pW, 4*sizeof(ae_int8)); AE_L8X4S_IP(W01, pW, 4*sizeof(ae_int8));
                for (j = 0; j < (Q & ~7); j += 8)
                {
                    AE_L64X2_IP(S0, S1, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S2, S3, pT_read, sizeof(ae_int64x2));
                    AE_L64X2_IP(S4, S5, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S6, S7, pT_read, sizeof(ae_int64x2));

                    AE_LA8X4S_IP(Y01, al0, pY0); AE_LA8X4S_IP(Y02, al0, pY0);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                    Y00 = Y02;

                    AE_S64X2_IP(S0, S1, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S2, S3, pT, sizeof(ae_int64x2));
                    AE_S64X2_IP(S4, S5, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S6, S7, pT, sizeof(ae_int64x2));

                }
                // tail
                if (j < Q)
                {
                    AE_L64X2_IP(S0, S1, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S2, S3, pT_read, sizeof(ae_int64x2));
                    AE_L64X2_IP(S4, S5, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S6, S7, pT_read, sizeof(ae_int64x2));

                    AE_LA8X4S_IP(Y01, al0, pY0); AE_LA8X4S_IP(Y02, al0, pY0);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                    AE_S64X2_IP(S0, S1, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S2, S3, pT, sizeof(ae_int64x2));
                    AE_S64X2_IP(S4, S5, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S6, S7, pT, sizeof(ae_int64x2));
                    for (; j < Q; ++j)
                    {
                        AE_DSEL16X4(Y00, Y01, Y00, Y01, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                        AE_DSEL16X4(Y01, Y02, Y01, Y02, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    }
                }
                //last N-1 elements
                {
                    AE_L64X2_IP(S0, S1, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S2, S3, pT_read, sizeof(ae_int64x2));
                    AE_L64X2_IP(S4, S5, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S6, S7, pT_read, sizeof(ae_int64x2));

                    Y01 = 0; Y02 = 0;

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                    AE_S64X2_IP(S0, S1, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S2, S3, pT, sizeof(ae_int64x2));
                    AE_S64X2_IP(S4, S5, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S6, S7, pT, sizeof(ae_int64x2));
                }

                pY0 = pY;
                pY = (int8_t*)(pY + Q);
            }
            //last iteration
            {
                pT_read = (ae_int64x2*)pScr + 3;

                Y00 = 0;

                al0 = AE_LA64_PP(pY0);

                AE_L8X4S_IP(W00, pW, 4 * sizeof(ae_int8)); AE_L8X4S_IP(W01, pW, 4 * sizeof(ae_int8));
                for (j = 0; j < (Q & ~7); j += 8)
                {
                    AE_L64X2_IP(S0, S1, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S2, S3, pT_read, sizeof(ae_int64x2));
                    AE_L64X2_IP(S4, S5, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S6, S7, pT_read, sizeof(ae_int64x2));

                    AE_LA8X4S_IP(Y01, al0, pY0); AE_LA8X4S_IP(Y02, al0, pY0);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                    Y00 = Y02;

                    T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                    T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 48 - rsh), AE_TRUNCA32X2F64S(S6, S7, 48 - rsh));
                    R0 = AE_ROUND8X8F16SSYM(T0, T1);

                    AE_SA8X8_IP(R0, alZ, pZ0);
                    
                }
                AE_SA64POS_FP(alZ, pZ0);
                // tail
                if (j < Q)
                {
                    AE_L64X2_IP(S0, S1, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S2, S3, pT_read, sizeof(ae_int64x2));
                    AE_L64X2_IP(S4, S5, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S6, S7, pT_read, sizeof(ae_int64x2));

                    AE_LA8X4S_IP(Y01, al0, pY0); AE_LA8X4S_IP(Y02, al0, pY0);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                    T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                    T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 48 - rsh), AE_TRUNCA32X2F64S(S6, S7, 48 - rsh));
                    R0 = AE_ROUND8X8F16SSYM(T0, T1);
                    for (; j < Q; ++j)
                    {
                        R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                        AE_DSEL16X4(Y00, Y01, Y00, Y01, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                        AE_DSEL16X4(Y01, Y02, Y01, Y02, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                        AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                    }
                }
                //last N-1 elements
                {
                    AE_L64X2_IP(S0, S1, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S2, S3, pT_read, sizeof(ae_int64x2));
                    AE_L64X2_IP(S4, S5, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S6, S7, pT_read, sizeof(ae_int64x2));

                    Y01 = Y02 = 0;

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                    T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                    T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 48 - rsh), AE_TRUNCA32X2F64S(S6, S7, 48 - rsh));
                    R0 = AE_ROUND8X8F16SSYM(T0, T1);


                    for (; j < Q+N-1; ++j)
                    {
                        R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                        AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                    }
                }
            }
        }
        else
        {
            Y00 = 0;
            al0 = AE_LA64_PP(pY0);

            AE_L8X4S_IP(W00, pW, 4 * sizeof(ae_int8)); AE_L8X4S_IP(W01, pW, 4 * sizeof(ae_int8));
            for (j = 0; j < (Q & ~7); j += 8)
            {

                AE_LA8X4S_IP(Y01, al0, pY0); AE_LA8X4S_IP(Y02, al0, pY0);

                AE_MULFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                AE_MULFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                AE_MULFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                AE_MULFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                Y00 = Y02;

                T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 48 - rsh), AE_TRUNCA32X2F64S(S6, S7, 48 - rsh));
                R0 = AE_ROUND8X8F16SSYM(T0, T1);

                AE_SA8X8_IP(R0, alZ, pZ0);

            }
            AE_SA64POS_FP(alZ, pZ0);
            // tail
            if (j < Q)
            {
                AE_LA8X4S_IP(Y01, al0, pY0); AE_LA8X4S_IP(Y02, al0, pY0);

                AE_MULFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                AE_MULFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                AE_MULFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                AE_MULFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 48 - rsh), AE_TRUNCA32X2F64S(S6, S7, 48 - rsh));
                R0 = AE_ROUND8X8F16SSYM(T0, T1);

                for (; j < Q; ++j)
                {
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_DSEL16X4(Y00, Y01, Y00, Y01, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    AE_DSEL16X4(Y01, Y02, Y01, Y02, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
            //last N-1 elements
            {

                Y01 = Y02 = 0;

                AE_MULFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                AE_MULFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                AE_MULFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                AE_MULFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 48 - rsh), AE_TRUNCA32X2F64S(S6, S7, 48 - rsh));
                R0 = AE_ROUND8X8F16SSYM(T0, T1);

                for (; j < Q + N - 1; ++j)
                {
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
        }
    }
#undef M
#undef N
}

size_t conv2d_gen_5x5_8x8_getScratchSize(int P, int Q)
{
    return 6 * 8 * sizeof(int8_t) + (Q + 16) * sizeof(int64_t);
} // MxN=5x5


#endif
