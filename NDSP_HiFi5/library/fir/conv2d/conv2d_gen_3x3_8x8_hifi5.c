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

#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )
#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )

#if defined(AE_MUL4O8X8CNV_H) // NN extension

void conv2d_gen_3x3_8x8(void * pScr, int8_t * z, const int8_t * x, const int8_t * y, int rsh, int P, int Q)
{
#define M 3
#define N 3
    const ae_int8x8* restrict pY0;
    const ae_int8x8* restrict pY1;
    const ae_int8x8* restrict pY2;
          ae_int8x8* restrict pW;
          ae_int8x8* restrict pZ0;
    ae_int8x8 v00, v01, v10, v11, v20, v21;
    ae_int32x2 S0, S1, S2, S3;
    ae_int16x4 T0, T1;
    ae_int8x8 R0;
    const ae_int8* pX = (const ae_int8*)x;

    ae_valign al0, al1, al2, alZ;

    int i, j, mmin, mmax, xstart, ystart;


    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0);
    NASSERT(Q >= 0);
    if (P <= 0 || Q <= 0) return;

    ae_int8x8 Y00,Y01;
    ae_int8x8 Y10,Y11;
    ae_int8x8 Y20,Y21;

    ae_int8x8 W0,W1,W2;
    
    alZ = AE_ZALIGN64();
    
    pW = (ae_int8x8*)pScr;
    // preload X and save in scratch
    {
        ae_int8x8 w0, w1, w2, w3;
        w0 = AE_MOVINT8X8_FROMINT16(0);
        AE_L8_IP(w3, pX, sizeof(ae_int8));
        AE_L8_IP(w2, pX, sizeof(ae_int8));
        AE_L8_IP(w1, pX, sizeof(ae_int8));
            
        W2 = AE_SEL8X8I(w0, w0, 16);
        W2 = AE_SEL8X8I(W2, w1, 16);
        W2 = AE_SEL8X8I(W2,w2,16);
        W2 = AE_SEL8X8I(W2, w3, 16);
        W2 = AE_SEL8X8I(W2,w0,3);

        AE_L8_IP(w3, pX, sizeof(ae_int8));
        AE_L8_IP(w2, pX, sizeof(ae_int8));
        AE_L8_IP(w1, pX, sizeof(ae_int8));

        W1 = AE_SEL8X8I(w0, w0, 16);
        W1 = AE_SEL8X8I(W1, w1, 16);
        W1 = AE_SEL8X8I(W1, w2, 16);
        W1 = AE_SEL8X8I(W1, w3, 16);
        W1 = AE_SEL8X8I(W1, w0, 3);

        AE_L8_IP(w3, pX, sizeof(ae_int8));
        AE_L8_IP(w2, pX, sizeof(ae_int8));
        AE_L8_IP(w1, pX, sizeof(ae_int8));

        W0 = AE_SEL8X8I(w0, w0, 16);
        W0 = AE_SEL8X8I(W0, w1, 16);
        W0 = AE_SEL8X8I(W0, w2, 16);
        W0 = AE_SEL8X8I(W0, w3, 16);
        W0 = AE_SEL8X8I(W0, w0, 3);

        AE_S8X8_IP(W0, pW, sizeof(ae_int8x8));
        AE_S8X8_IP(W1, pW, sizeof(ae_int8x8));
        AE_S8X8_IP(W2, pW, sizeof(ae_int8x8));
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
        pY1 = (ae_int8x8*)((uintptr_t)pY0 + Q);
        pY2 = (ae_int8x8*)((uintptr_t)pY1 + Q);
        pZ0 = (ae_int8x8*)(z + i * (Q + N - 1));

        // all 3 row needed
        if (mmax - mmin == 3)
        {
            Y00 = AE_MOVINT8X8_FROMINT16(0);
            Y10 = AE_MOVINT8X8_FROMINT16(0);
            Y20 = AE_MOVINT8X8_FROMINT16(0);

            al0 = AE_LA64_PP(pY0);
            al1 = AE_LA64_PP(pY1);
            al2 = AE_LA64_PP(pY2);

            AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));
            AE_L8X8_IP(W1, pW, sizeof(ae_int8x8));
            AE_L8X8_IP(W2, pW, sizeof(ae_int8x8));

            for (j = 0; j < (Q & ~7); j += 8)
            {                   
                AE_LA8X8_IP(Y01, al0, pY0);
                AE_LA8X8_IP(Y11, al1, pY1);
                AE_LA8X8_IP(Y21, al2, pY2);

                AE_DSEL8X8(v00, v01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                AE_DSEL8X8(v10, v11, Y10, Y11, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                AE_DSEL8X8(v20, v21, Y20, Y21, AE_MOVINT8X8_FROMINT64(0xA291807060504030));

                AE_MOVD8X16(Y00, Y10, Y01, Y11); Y20 = Y21;
                AE_MUL4O8X8CNV_H(S0, S1, S2, S3, W0, v00, v01);
                AE_MULA4O8X8CNV_H(S0, S1, S2, S3, W1, v10, v11);
                AE_MULA4O8X8CNV_H(S0, S1, S2, S3, W2, v20, v21);

                T0 = AE_TRUNCA16X4F32S(S0,S1,17-rsh); T1 = AE_TRUNCA16X4F32S(S2,S3,17-rsh);
                AE_SA8X8_IP(AE_ROUND8X8F16SSYM(T0, T1), alZ, pZ0);
            }
            AE_SA64POS_FP(alZ, pZ0);
            // tail
            if (j<Q)
            {
                AE_LA8X8_IP(Y01, al0, pY0);
                AE_LA8X8_IP(Y11, al1, pY1);
                AE_LA8X8_IP(Y21, al2, pY2);

                AE_DSEL8X8(v00, v01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                AE_DSEL8X8(v10, v11, Y10, Y11, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                AE_DSEL8X8(v20, v21, Y20, Y21, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
               
                AE_MUL4O8X8CNV_H(S0, S1, S2, S3, W0, v00, v01);
                AE_MULA4O8X8CNV_H(S0, S1, S2, S3, W1, v10, v11);
                AE_MULA4O8X8CNV_H(S0, S1, S2, S3, W2, v20, v21);

                T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                R0 = AE_ROUND8X8F16SSYM(T0, T1);

                for (; j < Q; ++j) 
                {
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_DSEL8X8(Y00, Y01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xE6D5C4B3A2918070));
                    AE_DSEL8X8(Y10, Y11, Y10, Y11, AE_MOVINT8X8_FROMINT64(0xE6D5C4B3A2918070));
                    AE_DSEL8X8(Y20, Y21, Y20, Y21, AE_MOVINT8X8_FROMINT64(0xE6D5C4B3A2918070));
                    AE_S8_0_IP(R0,castxcc(ae_int8,pZ0),sizeof(ae_int8));
                }
            }
            //last N-1 elements
            {
                Y01 = AE_MOVINT8X8_FROMINT16(0);
                Y11 = AE_MOVINT8X8_FROMINT16(0);
                Y21 = AE_MOVINT8X8_FROMINT16(0);

                AE_DSEL8X8(v00, v01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                AE_DSEL8X8(v10, v11, Y10, Y11, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                AE_DSEL8X8(v20, v21, Y20, Y21, AE_MOVINT8X8_FROMINT64(0xA291807060504030));

                AE_MUL4O8X8CNV_H(S0, S1, S2, S3, W0, v00, v01);
                AE_MULA4O8X8CNV_H(S0, S1, S2, S3, W1, v10, v11);
                AE_MULA4O8X8CNV_H(S0, S1, S2, S3, W2, v20, v21);

                T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                R0 = AE_ROUND8X8F16SSYM(T0, T1);

                for (j = Q; j < Q + N - 1; ++j)
                {
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
        }
        // if only one Y row need to 
        else if (mmax - mmin == 1)
        {
            Y00 = AE_MOVINT8X8_FROMINT16(0);

            al0 = AE_LA64_PP(pY0);

            AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));

            for (j = 0; j < (Q & ~7); j += 8)
            {
                AE_LA8X8_IP(Y01, al0, pY0);

                AE_DSEL8X8(v00, v01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                Y00 = Y01;
                AE_MUL4O8X8CNV_H(S0, S1, S2, S3, W0, v00, v01);

                T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                AE_SA8X8_IP(AE_ROUND8X8F16SSYM(T0, T1), alZ, pZ0);
            }
            AE_SA64POS_FP(alZ, pZ0);
            // tail
            if (j < Q)
            {
                AE_LA8X8_IP(Y01, al0, pY0);
                AE_DSEL8X8(v00, v01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                AE_MUL4O8X8CNV_H(S0, S1, S2, S3, W0, v00, v01);
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
                AE_DSEL8X8(v00, v01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                AE_MUL4O8X8CNV_H(S0, S1, S2, S3, W0, v00, v01);
                T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                R0 = AE_ROUND8X8F16SSYM(T0, T1);
                for (j = Q; j < Q + N - 1; ++j)
                {
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
        }
        else if (mmax - mmin == 2)
        {
            Y00 = AE_MOVINT8X8_FROMINT16(0);
            Y10 = AE_MOVINT8X8_FROMINT16(0);

            al0 = AE_LA64_PP(pY0);
            al1 = AE_LA64_PP(pY1);

            AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));
            AE_L8X8_IP(W1, pW, sizeof(ae_int8x8));

            for (j = 0; j < (Q & ~7); j += 8)
            {
                AE_LA8X8_IP(Y01, al0, pY0);
                AE_LA8X8_IP(Y11, al1, pY1);

                AE_DSEL8X8(v00, v01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                AE_DSEL8X8(v10, v11, Y10, Y11, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                AE_MOVD8X16(Y00, Y10, Y01, Y11);
                AE_MUL4O8X8CNV_H(S0, S1, S2, S3, W0, v00, v01);
                AE_MULA4O8X8CNV_H(S0, S1, S2, S3, W1, v10, v11);

                T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                AE_SA8X8_IP(AE_ROUND8X8F16SSYM(T0, T1), alZ, pZ0);
            }
            AE_SA64POS_FP(alZ, pZ0);
            // tail
            if (j < Q)
            {
                AE_LA8X8_IP(Y01, al0, pY0);
                AE_LA8X8_IP(Y11, al1, pY1);

                AE_DSEL8X8(v00, v01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                AE_DSEL8X8(v10, v11, Y10, Y11, AE_MOVINT8X8_FROMINT64(0xA291807060504030));

                AE_MUL4O8X8CNV_H(S0, S1, S2, S3, W0, v00, v01);
                AE_MULA4O8X8CNV_H(S0, S1, S2, S3, W1, v10, v11);

                T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                R0 = AE_ROUND8X8F16SSYM(T0, T1);

                for (; j < Q; ++j)
                {
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_DSEL8X8(Y00, Y01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xE6D5C4B3A2918070));
                    AE_DSEL8X8(Y10, Y11, Y10, Y11, AE_MOVINT8X8_FROMINT64(0xE6D5C4B3A2918070));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
            //last N-1 elements
            {
                Y01 = AE_MOVINT8X8_FROMINT16(0);
                Y11 = AE_MOVINT8X8_FROMINT16(0);

                AE_DSEL8X8(v00, v01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                AE_DSEL8X8(v10, v11, Y10, Y11, AE_MOVINT8X8_FROMINT64(0xA291807060504030));

                AE_MUL4O8X8CNV_H(S0, S1, S2, S3, W0, v00, v01);
                AE_MULA4O8X8CNV_H(S0, S1, S2, S3, W1, v10, v11);

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

// scratch allocatation functions. return required scratch size in bytes
size_t conv2d_gen_3x3_8x8_getScratchSize(int P, int Q)
{
    return 4 * sizeof(ae_int8x8);
}
#elif 1 // no NN extension
void conv2d_gen_3x3_8x8(void* pScr, int8_t* z, const int8_t* x, const int8_t* y, int rsh, int P, int Q)
{
#define M 3
#define N 3
    const int8_t* restrict pY0;
    const int8_t* restrict pY1;
    const int8_t* restrict pY2;
    int8_t* restrict pW;
    ae_int8x8* restrict pZ0;
    const ae_int8* pX = (const ae_int8*)x;

    ae_valign al0, al1, al2, alZ;

    int i, j, mmin, mmax, xstart, ystart;


    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0);
    NASSERT(Q >= 0);
    if (P <= 0 || Q <= 0) return;


    ae_int16x4 Y00, Y01, Y02;
    ae_int16x4 Y10, Y11, Y12;
    ae_int16x4 Y20, Y21, Y22;
    ae_int64 S0, S1, S2, S3;
    ae_int16x4 r0, r1;
    ae_int16x4 W0, W1, W2;
    ae_int8x8 R0;

    alZ = AE_ZALIGN64();

    pW = (int8_t*)pScr;
    // preload X and save in scratch
    {
        ae_int8x8 W0, W1, W2;

        ae_int8x8 w0, w1, w2, w3;
        w0 = AE_MOVINT8X8_FROMINT16(0);
        AE_L8_IP(w3, pX, sizeof(ae_int8));
        AE_L8_IP(w2, pX, sizeof(ae_int8));
        AE_L8_IP(w1, pX, sizeof(ae_int8));

        W2 = AE_SEL8X8I(w0, w0, 16);
        W2 = AE_SEL8X8I(W2, w1, 16);
        W2 = AE_SEL8X8I(W2, w2, 16);
        W2 = AE_SEL8X8I(W2, w3, 16);
        W2 = AE_SEL8X8I(W2, w0, 3);

        AE_L8_IP(w3, pX, sizeof(ae_int8));
        AE_L8_IP(w2, pX, sizeof(ae_int8));
        AE_L8_IP(w1, pX, sizeof(ae_int8));

        W1 = AE_SEL8X8I(w0, w0, 16);
        W1 = AE_SEL8X8I(W1, w1, 16);
        W1 = AE_SEL8X8I(W1, w2, 16);
        W1 = AE_SEL8X8I(W1, w3, 16);
        W1 = AE_SEL8X8I(W1, w0, 3);

        AE_L8_IP(w3, pX, sizeof(ae_int8));
        AE_L8_IP(w2, pX, sizeof(ae_int8));
        AE_L8_IP(w1, pX, sizeof(ae_int8));

        W0 = AE_SEL8X8I(w0, w0, 16);
        W0 = AE_SEL8X8I(W0, w1, 16);
        W0 = AE_SEL8X8I(W0, w2, 16);
        W0 = AE_SEL8X8I(W0, w3, 16);
        W0 = AE_SEL8X8I(W0, w0, 3);

        AE_S8X8_IP(W0, castxcc(ae_int8x8,pW), sizeof(ae_int8x8));
        AE_S8X8_IP(W1, castxcc(ae_int8x8,pW), sizeof(ae_int8x8));
        AE_S8X8_IP(W2, castxcc(ae_int8x8,pW), sizeof(ae_int8x8));
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
        pY1 = (int8_t*)((uintptr_t)pY0 + Q);
        pY2 = (int8_t*)((uintptr_t)pY1 + Q);
        pZ0 = (ae_int8x8*)(z + i * (Q + N - 1));

        // all 3 row needed
        if (mmax - mmin == 3)
        {
            Y00 = 0;
            Y10 = 0;
            Y20 = 0;

            al0 = AE_LA64_PP(pY0);
            al1 = AE_LA64_PP(pY1);
            al2 = AE_LA64_PP(pY2);

            AE_L8X4S_IP(W0, pW, sizeof(ae_int16x4));
            AE_L8X4S_IP(W1, pW, sizeof(ae_int16x4));
            AE_L8X4S_IP(W2, pW, sizeof(ae_int16x4));

            for (j = 0; j < (Q & ~7); j += 8)
            {


                AE_LA8X4S_IP(Y01, al0, pY0); AE_LA8X4S_IP(Y02, al0, pY0);
                AE_LA8X4S_IP(Y11, al1, pY1); AE_LA8X4S_IP(Y12, al1, pY1);
                AE_LA8X4S_IP(Y21, al2, pY2); AE_LA8X4S_IP(Y22, al2, pY2);
                
                AE_MULFQ16X2_FIR_2(S0, S1, Y00, Y01, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y00, Y01, W0);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y10, Y11, W1);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y10, Y11, W1);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y20, Y21, W2);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y20, Y21, W2);
                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, Y01, Y02, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y01, Y02, W0);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y11, Y12, W1);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y11, Y12, W1);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y21, Y22, W2);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y21, Y22, W2);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), alZ, pZ0);

                AE_MOVD16X8(Y00, Y10, Y02, Y12); Y20 = Y22;
            }
            AE_SA64POS_FP(alZ, pZ0);
            // tail
            if (j < Q)
            {
                AE_LA8X4S_IP(Y01, al0, pY0); AE_LA8X4S_IP(Y02, al0, pY0);
                AE_LA8X4S_IP(Y11, al1, pY1); AE_LA8X4S_IP(Y12, al1, pY1);
                AE_LA8X4S_IP(Y21, al2, pY2); AE_LA8X4S_IP(Y22, al2, pY2);

                AE_MULFQ16X2_FIR_2(S0, S1, Y00, Y01, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y00, Y01, W0);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y10, Y11, W1);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y10, Y11, W1);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y20, Y21, W2);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y20, Y21, W2);
                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, Y01, Y02, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y01, Y02, W0);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y11, Y12, W1);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y11, Y12, W1);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y21, Y22, W2);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y21, Y22, W2);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                R0 = AE_ROUND8X8F16SSYM(r0, r1);
                for (; j < Q; ++j)
                {
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_DSEL16X4(Y00, Y01, Y00, Y01, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    AE_DSEL16X4(Y01, Y02, Y01, Y02, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    AE_DSEL16X4(Y10, Y11, Y10, Y11, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    AE_DSEL16X4(Y11, Y12, Y11, Y12, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    AE_DSEL16X4(Y20, Y21, Y20, Y21, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    AE_DSEL16X4(Y21, Y22, Y21, Y22, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
            //last N-1 elements
            {
                Y01 = 0; Y02 = 0;
                Y11 = 0; Y12 = 0;
                Y21 = 0; Y22 = 0;

                AE_MULFQ16X2_FIR_2(S0, S1, Y00, Y01, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y00, Y01, W0);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y10, Y11, W1);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y10, Y11, W1);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y20, Y21, W2);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y20, Y21, W2);
                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, Y01, Y02, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y01, Y02, W0);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y11, Y12, W1);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y11, Y12, W1);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y21, Y22, W2);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y21, Y22, W2);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                R0 = AE_ROUND8X8F16SSYM(r0, r1);
                for (; j < Q+N-1; ++j)
                {
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
        }
        // if only one Y row need to 
        else if (mmax - mmin == 1)
        {
            Y00 = 0;

            al0 = AE_LA64_PP(pY0);

            AE_L8X4S_IP(W0, pW, sizeof(ae_int16x4));

            for (j = 0; j < (Q & ~7); j += 8)
            {

                AE_LA8X4S_IP(Y01, al0, pY0); AE_LA8X4S_IP(Y02, al0, pY0);

                AE_MULFQ16X2_FIR_2(S0, S1, Y00, Y01, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y00, Y01, W0);

                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, Y01, Y02, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y01, Y02, W0);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), alZ, pZ0);

                Y00 = Y02;
            }
            AE_SA64POS_FP(alZ, pZ0);
            // tail
            if (j < Q)
            {
                AE_LA8X4S_IP(Y01, al0, pY0); AE_LA8X4S_IP(Y02, al0, pY0);

                AE_MULFQ16X2_FIR_2(S0, S1, Y00, Y01, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y00, Y01, W0);

                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, Y01, Y02, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y01, Y02, W0);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                R0 = AE_ROUND8X8F16SSYM(r0, r1);
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
                Y01 = 0; Y02 = 0;
                Y11 = 0; Y12 = 0;
                Y21 = 0; Y22 = 0;

                AE_MULFQ16X2_FIR_2(S0, S1, Y00, Y01, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y00, Y01, W0);
                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, Y01, Y02, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y01, Y02, W0);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                R0 = AE_ROUND8X8F16SSYM(r0, r1);
                for (; j < Q + N - 1; ++j)
                {
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
        }
        else if (mmax - mmin == 2)
        {
            Y00 = 0;
            Y10 = 0;

            al0 = AE_LA64_PP(pY0);
            al1 = AE_LA64_PP(pY1);

            AE_L8X4S_IP(W0, pW, sizeof(ae_int16x4));
            AE_L8X4S_IP(W1, pW, sizeof(ae_int16x4));

            for (j = 0; j < (Q & ~7); j += 8)
            {


                AE_LA8X4S_IP(Y01, al0, pY0); AE_LA8X4S_IP(Y02, al0, pY0);
                AE_LA8X4S_IP(Y11, al1, pY1); AE_LA8X4S_IP(Y12, al1, pY1);

                AE_MULFQ16X2_FIR_2(S0, S1, Y00, Y01, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y00, Y01, W0);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y10, Y11, W1);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y10, Y11, W1);
                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, Y01, Y02, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y01, Y02, W0);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y11, Y12, W1);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y11, Y12, W1);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), alZ, pZ0);

                AE_MOVD16X8(Y00, Y10, Y02, Y12);
            }
            AE_SA64POS_FP(alZ, pZ0);
            // tail
            if (j < Q)
            {
                AE_LA8X4S_IP(Y01, al0, pY0); AE_LA8X4S_IP(Y02, al0, pY0);
                AE_LA8X4S_IP(Y11, al1, pY1); AE_LA8X4S_IP(Y12, al1, pY1);

                AE_MULFQ16X2_FIR_2(S0, S1, Y00, Y01, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y00, Y01, W0);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y10, Y11, W1);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y10, Y11, W1);
                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, Y01, Y02, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y01, Y02, W0);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y11, Y12, W1);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y11, Y12, W1);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                R0 = AE_ROUND8X8F16SSYM(r0, r1);
                for (; j < Q; ++j)
                {
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_DSEL16X4(Y00, Y01, Y00, Y01, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    AE_DSEL16X4(Y01, Y02, Y01, Y02, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    AE_DSEL16X4(Y10, Y11, Y10, Y11, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    AE_DSEL16X4(Y11, Y12, Y11, Y12, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
            //last N-1 elements
            {
                Y01 = 0; Y02 = 0;
                Y11 = 0; Y12 = 0;
                Y21 = 0; Y22 = 0;

                AE_MULFQ16X2_FIR_2(S0, S1, Y00, Y01, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y00, Y01, W0);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y10, Y11, W1);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y10, Y11, W1);
                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, Y01, Y02, W0);
                AE_MULFQ16X2_FIR_0(S2, S3, Y01, Y02, W0);
                AE_MULAFQ16X2_FIR_2(S0, S1, Y11, Y12, W1);
                AE_MULAFQ16X2_FIR_0(S2, S3, Y11, Y12, W1);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 48 - rsh), AE_TRUNCA32X2F64S(S2, S3, 48 - rsh));
                R0 = AE_ROUND8X8F16SSYM(r0, r1);
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

// scratch allocatation functions. return required scratch size in bytes
size_t conv2d_gen_3x3_8x8_getScratchSize(int P, int Q)
{
    return 4 * sizeof(ae_int16x4);
}

#else // reference code
void conv2d_gen_3x3_8x8(void* pScr, int8_t* z, const int8_t* x, const int8_t* y, int rsh, int P, int Q)
{
#define M 3
#define N 3
    const int8_t *restrict pY0;
    const int8_t *restrict pY1;
    const int8_t *restrict pY2;
          int8_t *restrict pZ0;
    int8_t *restrict pW;


    int i,j, mmin, mmax, xstart, ystart, rnd;

    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0);
    NASSERT(Q >= 0);
    if (P <= 0 || Q <= 0) return;
    
    rnd = (7 + rsh) > 0 ? (1 << (7 + rsh)) >> 1 : 0;

    int8_t y0_0, y0_1, y0_2,  y0_3,  y0_4,  y0_5,  y0_6,  y0_7;
    int8_t y0_8, y0_9, y0_10, y0_11, y0_12, y0_13, y0_14, y0_15;
    int8_t y1_0, y1_1, y1_2,  y1_3,  y1_4,  y1_5,  y1_6,  y1_7;
    int8_t y1_8, y1_9, y1_10, y1_11, y1_12, y1_13, y1_14, y1_15;
    int8_t y2_0, y2_1, y2_2,  y2_3,  y2_4,  y2_5,  y2_6,  y2_7;
    int8_t y2_8, y2_9, y2_10, y2_11, y2_12, y2_13, y2_14, y2_15;

    int32_t s0, s1, s2, s3, s4, s5, s6, s7;

    int8_t w0, w1, w2, w3; // low bits are eq 0
    pW = (int8_t*) pScr;
    // preload X and save in scratch
    {
        for (i = 0; i < 3; ++i)
        {
            w3 = *x++;
            w2 = *x++;
            w1 = *x++;
            w0 = 0;
            pW[(2 - i) * 4 + 0] = w0;
            pW[(2 - i) * 4 + 1] = w1;
            pW[(2 - i) * 4 + 2] = w2;
            pW[(2 - i) * 4 + 3] = w3;
        }
        // do some other staf
    }

    __Pragma("loop_count min=1");
    for (i = 0; i < P + M - 1; i += 1)
    {
        mmin = i < M - 1 ? M - 1 - i : 0;
        mmax = i >= P ? M - 1 + P - i : M;
        xstart = i < M - 1 ? (mmin) : 0;
        ystart = i < M - 1 ? 0 : i - (M - 1);

        pW = (int8_t*)pScr + 4*xstart;
        pY0 = (int8_t*)(y + ystart * Q);
        pY1 = (int8_t*)(pY0 + Q);
        pY2 = (int8_t*)(pY1 + Q);
        pZ0 = (int8_t*)(z + i * (Q + N - 1));
        
        // all 3 row needed
        if (mmax - mmin == 3)
        {
            {
                y0_0 = y0_1 = y0_2 = y0_3 = y0_4 = y0_5 = y0_6 = y0_7 = 0;
                y1_0 = y1_1 = y1_2 = y1_3 = y1_4 = y1_5 = y1_6 = y1_7 = 0;
                y2_0 = y2_1 = y2_2 = y2_3 = y2_4 = y2_5 = y2_6 = y2_7 = 0;
                

                for (j = 0; j < (Q & ~7); j += 8)
                {
                    
                    y0_8 = *pY0++; y0_9 = *pY0++; y0_10 = *pY0++; y0_11 = *pY0++; y0_12 = *pY0++; y0_13 = *pY0++; y0_14 = *pY0++; y0_15 = *pY0++;
                    w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;

                    s0 = (int32_t)w0 * y0_5; s0 +=  (int32_t)w1 * y0_6; s0 +=  (int32_t)w2 * y0_7; s0 += (int32_t)w3 *  y0_8;
                    s1 = (int32_t)w0 * y0_6; s1 +=  (int32_t)w1 * y0_7; s1 +=  (int32_t)w2 * y0_8; s1 += (int32_t)w3 *  y0_9;
                    s2 = (int32_t)w0 * y0_7; s2 +=  (int32_t)w1 * y0_8; s2 +=  (int32_t)w2 * y0_9; s2 += (int32_t)w3 *  y0_10;
                    s3 = (int32_t)w0 * y0_8; s3 +=  (int32_t)w1 * y0_9; s3 +=  (int32_t)w2 * y0_10; s3 += (int32_t)w3 * y0_11;
                    s4 = (int32_t)w0 * y0_9; s4 +=  (int32_t)w1 * y0_10; s4 += (int32_t)w2 * y0_11; s4 += (int32_t)w3 * y0_12;
                    s5 = (int32_t)w0 * y0_10; s5 += (int32_t)w1 * y0_11; s5 += (int32_t)w2 * y0_12; s5 += (int32_t)w3 * y0_13;
                    s6 = (int32_t)w0 * y0_11; s6 += (int32_t)w1 * y0_12; s6 += (int32_t)w2 * y0_13; s6 += (int32_t)w3 * y0_14;
                    s7 = (int32_t)w0 * y0_12; s7 += (int32_t)w1 * y0_13; s7 += (int32_t)w2 * y0_14; s7 += (int32_t)w3 * y0_15;
                    y0_0 = y0_8; y0_1 = y0_9; y0_2 = y0_10; y0_3 = y0_11; y0_4 = y0_12; y0_5 = y0_13; y0_6 = y0_14; y0_7 = y0_15;

                    w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;
                    y1_8 = *pY1++; y1_9 = *pY1++; y1_10 = *pY1++; y1_11 = *pY1++; y1_12 = *pY1++; y1_13 = *pY1++; y1_14 = *pY1++; y1_15 = *pY1++;
                    s0 += (int32_t)w0 * y1_5; s0 += (int32_t)w1 * y1_6; s0 += (int32_t)w2 * y1_7; s0 += (int32_t)w3 * y1_8;
                    s1 += (int32_t)w0 * y1_6; s1 += (int32_t)w1 * y1_7; s1 += (int32_t)w2 * y1_8; s1 += (int32_t)w3 * y1_9;
                    s2 += (int32_t)w0 * y1_7; s2 += (int32_t)w1 * y1_8; s2 += (int32_t)w2 * y1_9; s2 += (int32_t)w3 * y1_10;
                    s3 += (int32_t)w0 * y1_8; s3 += (int32_t)w1 * y1_9; s3 += (int32_t)w2 * y1_10; s3 += (int32_t)w3 * y1_11;
                    s4 += (int32_t)w0 * y1_9; s4 += (int32_t)w1 * y1_10; s4 += (int32_t)w2 * y1_11; s4 += (int32_t)w3 * y1_12;
                    s5 += (int32_t)w0 * y1_10; s5 += (int32_t)w1 * y1_11; s5 += (int32_t)w2 * y1_12; s5 += (int32_t)w3 * y1_13;
                    s6 += (int32_t)w0 * y1_11; s6 += (int32_t)w1 * y1_12; s6 += (int32_t)w2 * y1_13; s6 += (int32_t)w3 * y1_14;
                    s7 += (int32_t)w0 * y1_12; s7 += (int32_t)w1 * y1_13; s7 += (int32_t)w2 * y1_14; s7 += (int32_t)w3 * y1_15;
                    y1_0 = y1_8; y1_1 = y1_9; y1_2 = y1_10; y1_3 = y1_11; y1_4 = y1_12; y1_5 = y1_13; y1_6 = y1_14; y1_7 = y1_15;

                    w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;
                    y2_8 = *pY2++; y2_9 = *pY2++; y2_10 = *pY2++; y2_11 = *pY2++; y2_12 = *pY2++; y2_13 = *pY2++; y2_14 = *pY2++; y2_15 = *pY2++;
                    s0 += (int32_t)w0 * y2_5; s0 += (int32_t)w1 *  y2_6; s0 += (int32_t)w2 *  y2_7; s0 += (int32_t)w3 *  y2_8;
                    s1 += (int32_t)w0 * y2_6; s1 += (int32_t)w1 *  y2_7; s1 += (int32_t)w2 *  y2_8; s1 += (int32_t)w3 *  y2_9;
                    s2 += (int32_t)w0 * y2_7; s2 += (int32_t)w1 *  y2_8; s2 += (int32_t)w2 *  y2_9; s2 += (int32_t)w3 *  y2_10;
                    s3 += (int32_t)w0 * y2_8; s3 += (int32_t)w1 *  y2_9; s3 += (int32_t)w2 *  y2_10; s3 += (int32_t)w3 * y2_11;
                    s4 += (int32_t)w0 * y2_9; s4 += (int32_t)w1 *  y2_10; s4 += (int32_t)w2 * y2_11; s4 += (int32_t)w3 * y2_12;
                    s5 += (int32_t)w0 * y2_10; s5 += (int32_t)w1 * y2_11; s5 += (int32_t)w2 * y2_12; s5 += (int32_t)w3 * y2_13;
                    s6 += (int32_t)w0 * y2_11; s6 += (int32_t)w1 * y2_12; s6 += (int32_t)w2 * y2_13; s6 += (int32_t)w3 * y2_14;
                    s7 += (int32_t)w0 * y2_12; s7 += (int32_t)w1 * y2_13; s7 += (int32_t)w2 * y2_14; s7 += (int32_t)w3 * y2_15;
                    y2_0 = y2_8; y2_1 = y2_9; y2_2 = y2_10; y2_3 = y2_11; y2_4 = y2_12; y2_5 = y2_13; y2_6 = y2_14; y2_7 = y2_15;

                    pW = (int8_t*)pScr;

                    s0 += rnd; s0 = (7 + rsh) > 0 ? s0 >> (7 + rsh) : s0 << -(7 + rsh); s0 = MIN(127, MAX(s0, -128));
                    s1 += rnd; s1 = (7 + rsh) > 0 ? s1 >> (7 + rsh) : s1 << -(7 + rsh); s1 = MIN(127, MAX(s1, -128));
                    s2 += rnd; s2 = (7 + rsh) > 0 ? s2 >> (7 + rsh) : s2 << -(7 + rsh); s2 = MIN(127, MAX(s2, -128));
                    s3 += rnd; s3 = (7 + rsh) > 0 ? s3 >> (7 + rsh) : s3 << -(7 + rsh); s3 = MIN(127, MAX(s3, -128));
                    s4 += rnd; s4 = (7 + rsh) > 0 ? s4 >> (7 + rsh) : s4 << -(7 + rsh); s4 = MIN(127, MAX(s4, -128));
                    s5 += rnd; s5 = (7 + rsh) > 0 ? s5 >> (7 + rsh) : s5 << -(7 + rsh); s5 = MIN(127, MAX(s5, -128));
                    s6 += rnd; s6 = (7 + rsh) > 0 ? s6 >> (7 + rsh) : s6 << -(7 + rsh); s6 = MIN(127, MAX(s6, -128));
                    s7 += rnd; s7 = (7 + rsh) > 0 ? s7 >> (7 + rsh) : s7 << -(7 + rsh); s7 = MIN(127, MAX(s7, -128));

                    *pZ0++ = (int8_t)s0;
                    *pZ0++ = (int8_t)s1;
                    *pZ0++ = (int8_t)s2;
                    *pZ0++ = (int8_t)s3;
                    *pZ0++ = (int8_t)s4;
                    *pZ0++ = (int8_t)s5;
                    *pZ0++ = (int8_t)s6;
                    *pZ0++ = (int8_t)s7;

                }
                for (; j < Q; ++j)
                {
                    y0_8 = *pY0++;
                    y1_8 = *pY1++;
                    y2_8 = *pY2++;
                    w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;
                    s0 = (int32_t)w0 * y0_5; s0 += (int32_t)w1 * y0_6; s0 += (int32_t)w2 * y0_7; s0 += (int32_t)w3 * y0_8;
                    w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;
                    s0 += (int32_t)w0 * y1_5; s0 += (int32_t)w1 * y1_6; s0 += (int32_t)w2 * y1_7; s0 += (int32_t)w3 * y1_8;
                    w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;
                    s0 += (int32_t)w0 * y2_5; s0 += (int32_t)w1 * y2_6; s0 += (int32_t)w2 * y2_7; s0 += (int32_t)w3 * y2_8;
                    pW = (int8_t*)pScr;

                    s0 += rnd; s0 = (7 + rsh) > 0 ? s0 >> (7 + rsh) : s0 << -(7 + rsh); s0 = MIN(127, MAX(s0, -128));
                    *pZ0++ = (int8_t)s0;
                    y0_5 = y0_6; y0_6 = y0_7; y0_7 = y0_8;
                    y1_5 = y1_6; y1_6 = y1_7; y1_7 = y1_8;
                    y2_5 = y2_6; y2_6 = y2_7; y2_7 = y2_8;
                }
                for (; j < Q + N - 1; ++j)
                {
                    y0_8 = 0;
                    y1_8 = 0;
                    y2_8 = 0;
                    w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;
                    s0 = (int32_t)w0 * y0_5; s0 += (int32_t)w1 * y0_6; s0 += (int32_t)w2 * y0_7; s0 += (int32_t)w3 * y0_8;
                    w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;
                    s0 += (int32_t)w0 * y1_5; s0 += (int32_t)w1 * y1_6; s0 += (int32_t)w2 * y1_7; s0 += (int32_t)w3 * y1_8;
                    w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;
                    s0 += (int32_t)w0 * y2_5; s0 += (int32_t)w1 * y2_6; s0 += (int32_t)w2 * y2_7; s0 += (int32_t)w3 * y2_8;
                    pW = (int8_t*)pScr;

                    s0 += rnd; s0 = (7 + rsh) > 0 ? s0 >> (7 + rsh) : s0 << -(7 + rsh); s0 = MIN(127, MAX(s0, -128));
                    *pZ0++ = (int8_t)s0;
                    y0_5 = y0_6; y0_6 = y0_7; y0_7 = y0_8;
                    y1_5 = y1_6; y1_6 = y1_7; y1_7 = y1_8;
                    y2_5 = y2_6; y2_6 = y2_7; y2_7 = y2_8;
                }
            }
        }
        // if only one Y row need to 
        else if (mmax - mmin == 1)
        {

            
            y0_0 = y0_1 = y0_2 = y0_3 = y0_4 = y0_5 = y0_6 = y0_7 = 0;
            w0 = *pW++;
            w1 = *pW++;
            w2 = *pW++;
            w3 = *pW++;

            for (j = 0; j < (Q&~7); j +=8)
            {
                y0_8 = *pY0++; y0_9 = *pY0++; y0_10 = *pY0++; y0_11 = *pY0++; y0_12 = *pY0++; y0_13 = *pY0++; y0_14 = *pY0++; y0_15 = *pY0++;

                s0 = (int32_t)w0 * y0_5; s0 += (int32_t)w1 * y0_6; s0 += (int32_t)w2 * y0_7; s0 += (int32_t)w3 * y0_8;
                s1 = (int32_t)w0 * y0_6; s1 += (int32_t)w1 * y0_7; s1 += (int32_t)w2 * y0_8; s1 += (int32_t)w3 * y0_9;
                s2 = (int32_t)w0 * y0_7; s2 += (int32_t)w1 * y0_8; s2 += (int32_t)w2 * y0_9; s2 += (int32_t)w3 * y0_10;
                s3 = (int32_t)w0 * y0_8; s3 += (int32_t)w1 * y0_9; s3 += (int32_t)w2 * y0_10; s3 += (int32_t)w3 * y0_11;
                s4 = (int32_t)w0 * y0_9; s4 += (int32_t)w1 * y0_10; s4 += (int32_t)w2 * y0_11; s4 += (int32_t)w3 * y0_12;
                s5 = (int32_t)w0 * y0_10; s5 += (int32_t)w1 * y0_11; s5 += (int32_t)w2 * y0_12; s5 += (int32_t)w3 * y0_13;
                s6 = (int32_t)w0 * y0_11; s6 += (int32_t)w1 * y0_12; s6 += (int32_t)w2 * y0_13; s6 += (int32_t)w3 * y0_14;
                s7 = (int32_t)w0 * y0_12; s7 += (int32_t)w1 * y0_13; s7 += (int32_t)w2 * y0_14; s7 += (int32_t)w3 * y0_15;
                y0_0 = y0_8; y0_1 = y0_9; y0_2 = y0_10; y0_3 = y0_11; y0_4 = y0_12; y0_5 = y0_13; y0_6 = y0_14; y0_7 = y0_15;

                s0 += rnd; s0 = (7 + rsh) > 0 ? s0 >> (7 + rsh) : s0 << -(7 + rsh); s0 = MIN(127, MAX(s0, -128));
                s1 += rnd; s1 = (7 + rsh) > 0 ? s1 >> (7 + rsh) : s1 << -(7 + rsh); s1 = MIN(127, MAX(s1, -128));
                s2 += rnd; s2 = (7 + rsh) > 0 ? s2 >> (7 + rsh) : s2 << -(7 + rsh); s2 = MIN(127, MAX(s2, -128));
                s3 += rnd; s3 = (7 + rsh) > 0 ? s3 >> (7 + rsh) : s3 << -(7 + rsh); s3 = MIN(127, MAX(s3, -128));
                s4 += rnd; s4 = (7 + rsh) > 0 ? s4 >> (7 + rsh) : s4 << -(7 + rsh); s4 = MIN(127, MAX(s4, -128));
                s5 += rnd; s5 = (7 + rsh) > 0 ? s5 >> (7 + rsh) : s5 << -(7 + rsh); s5 = MIN(127, MAX(s5, -128));
                s6 += rnd; s6 = (7 + rsh) > 0 ? s6 >> (7 + rsh) : s6 << -(7 + rsh); s6 = MIN(127, MAX(s6, -128));
                s7 += rnd; s7 = (7 + rsh) > 0 ? s7 >> (7 + rsh) : s7 << -(7 + rsh); s7 = MIN(127, MAX(s7, -128));
                
                *pZ0++ = (int8_t)s0;
                *pZ0++ = (int8_t)s1;
                *pZ0++ = (int8_t)s2;
                *pZ0++ = (int8_t)s3;
                *pZ0++ = (int8_t)s4;
                *pZ0++ = (int8_t)s5;
                *pZ0++ = (int8_t)s6;
                *pZ0++ = (int8_t)s7;
            }
            for (; j < Q; ++j)
            {
                y0_8 = *pY0++;
                s0 = (int32_t)w0 * y0_5; s0 += (int32_t)w1 * y0_6; s0 += (int32_t)w2 * y0_7; s0 += (int32_t)w3 * y0_8;
                s0 += rnd; s0 = (7 + rsh) > 0 ? s0 >> (7 + rsh) : s0 << -(7 + rsh); s0 = MIN(127, MAX(s0, -128));
                *pZ0++ = (int8_t)s0;
                y0_5 = y0_6; y0_6 = y0_7; y0_7 = y0_8;
            }
            for (; j < Q+N-1; ++j)
            {
                y0_8 = 0;
                s0 = (int32_t)w0 * y0_5; s0 += (int32_t)w1 * y0_6; s0 += (int32_t)w2 * y0_7; s0 += (int32_t)w3 * y0_8;
                s0 += rnd; s0 = (7 + rsh) > 0 ? s0 >> (7 + rsh) : s0 << -(7 + rsh); s0 = MIN(127, MAX(s0, -128));
                *pZ0++ = (int8_t)s0;
                y0_5 = y0_6; y0_6 = y0_7; y0_7 = y0_8;
            }
        }
        else if (mmax - mmin == 2)
        {
            y0_0 = y0_1 = y0_2 = y0_3 = y0_4 = y0_5 = y0_6 = y0_7 = 0;
            y1_0 = y1_1 = y1_2 = y1_3 = y1_4 = y1_5 = y1_6 = y1_7 = 0;

            for (j = 0; j < (Q & ~7); j += 8)
            {
               

                y0_8 = *pY0++; y0_9 = *pY0++; y0_10 = *pY0++; y0_11 = *pY0++; y0_12 = *pY0++; y0_13 = *pY0++; y0_14 = *pY0++; y0_15 = *pY0++;
                w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;

                s0 = (int32_t)w0 * y0_5; s0 += (int32_t)w1 * y0_6; s0 += (int32_t)w2 * y0_7; s0 += (int32_t)w3 * y0_8;
                s1 = (int32_t)w0 * y0_6; s1 += (int32_t)w1 * y0_7; s1 += (int32_t)w2 * y0_8; s1 += (int32_t)w3 * y0_9;
                s2 = (int32_t)w0 * y0_7; s2 += (int32_t)w1 * y0_8; s2 += (int32_t)w2 * y0_9; s2 += (int32_t)w3 * y0_10;
                s3 = (int32_t)w0 * y0_8; s3 += (int32_t)w1 * y0_9; s3 += (int32_t)w2 * y0_10; s3 += (int32_t)w3 * y0_11;
                s4 = (int32_t)w0 * y0_9; s4 += (int32_t)w1 * y0_10; s4 += (int32_t)w2 * y0_11; s4 += (int32_t)w3 * y0_12;
                s5 = (int32_t)w0 * y0_10; s5 += (int32_t)w1 * y0_11; s5 += (int32_t)w2 * y0_12; s5 += (int32_t)w3 * y0_13;
                s6 = (int32_t)w0 * y0_11; s6 += (int32_t)w1 * y0_12; s6 += (int32_t)w2 * y0_13; s6 += (int32_t)w3 * y0_14;
                s7 = (int32_t)w0 * y0_12; s7 += (int32_t)w1 * y0_13; s7 += (int32_t)w2 * y0_14; s7 += (int32_t)w3 * y0_15;
                y0_0 = y0_8; y0_1 = y0_9; y0_2 = y0_10; y0_3 = y0_11; y0_4 = y0_12; y0_5 = y0_13; y0_6 = y0_14; y0_7 = y0_15;

                w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;
                y1_8 = *pY1++; y1_9 = *pY1++; y1_10 = *pY1++; y1_11 = *pY1++; y1_12 = *pY1++; y1_13 = *pY1++; y1_14 = *pY1++; y1_15 = *pY1++;
                s0 += (int32_t)w0 * y1_5; s0 += (int32_t)w1 * y1_6; s0 += (int32_t)w2 * y1_7; s0 += (int32_t)w3 * y1_8;
                s1 += (int32_t)w0 * y1_6; s1 += (int32_t)w1 * y1_7; s1 += (int32_t)w2 * y1_8; s1 += (int32_t)w3 * y1_9;
                s2 += (int32_t)w0 * y1_7; s2 += (int32_t)w1 * y1_8; s2 += (int32_t)w2 * y1_9; s2 += (int32_t)w3 * y1_10;
                s3 += (int32_t)w0 * y1_8; s3 += (int32_t)w1 * y1_9; s3 += (int32_t)w2 * y1_10; s3 += (int32_t)w3 * y1_11;
                s4 += (int32_t)w0 * y1_9; s4 += (int32_t)w1 * y1_10; s4 += (int32_t)w2 * y1_11; s4 += (int32_t)w3 * y1_12;
                s5 += (int32_t)w0 * y1_10; s5 += (int32_t)w1 * y1_11; s5 += (int32_t)w2 * y1_12; s5 += (int32_t)w3 * y1_13;
                s6 += (int32_t)w0 * y1_11; s6 += (int32_t)w1 * y1_12; s6 += (int32_t)w2 * y1_13; s6 += (int32_t)w3 * y1_14;
                s7 += (int32_t)w0 * y1_12; s7 += (int32_t)w1 * y1_13; s7 += (int32_t)w2 * y1_14; s7 += (int32_t)w3 * y1_15;
                y1_0 = y1_8; y1_1 = y1_9; y1_2 = y1_10; y1_3 = y1_11; y1_4 = y1_12; y1_5 = y1_13; y1_6 = y1_14; y1_7 = y1_15;
                pW -= 8;

                s0 += rnd; s0 = (7 + rsh) > 0 ? s0 >> (7 + rsh) : s0 << -(7 + rsh); s0 = MIN(127, MAX(s0, -128));
                s1 += rnd; s1 = (7 + rsh) > 0 ? s1 >> (7 + rsh) : s1 << -(7 + rsh); s1 = MIN(127, MAX(s1, -128));
                s2 += rnd; s2 = (7 + rsh) > 0 ? s2 >> (7 + rsh) : s2 << -(7 + rsh); s2 = MIN(127, MAX(s2, -128));
                s3 += rnd; s3 = (7 + rsh) > 0 ? s3 >> (7 + rsh) : s3 << -(7 + rsh); s3 = MIN(127, MAX(s3, -128));
                s4 += rnd; s4 = (7 + rsh) > 0 ? s4 >> (7 + rsh) : s4 << -(7 + rsh); s4 = MIN(127, MAX(s4, -128));
                s5 += rnd; s5 = (7 + rsh) > 0 ? s5 >> (7 + rsh) : s5 << -(7 + rsh); s5 = MIN(127, MAX(s5, -128));
                s6 += rnd; s6 = (7 + rsh) > 0 ? s6 >> (7 + rsh) : s6 << -(7 + rsh); s6 = MIN(127, MAX(s6, -128));
                s7 += rnd; s7 = (7 + rsh) > 0 ? s7 >> (7 + rsh) : s7 << -(7 + rsh); s7 = MIN(127, MAX(s7, -128));

                *pZ0++ = (int8_t)s0;
                *pZ0++ = (int8_t)s1;
                *pZ0++ = (int8_t)s2;
                *pZ0++ = (int8_t)s3;
                *pZ0++ = (int8_t)s4;
                *pZ0++ = (int8_t)s5;
                *pZ0++ = (int8_t)s6;
                *pZ0++ = (int8_t)s7;

            }
            for (; j < Q; ++j)
            {
                y0_8 = *pY0++;
                y1_8 = *pY1++;
                w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;
                s0 = (int32_t)w0 * y0_5; s0 += (int32_t)w1 * y0_6; s0 += (int32_t)w2 * y0_7; s0 += (int32_t)w3 * y0_8;
                w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;
                s0 += (int32_t)w0 * y1_5; s0 += (int32_t)w1 * y1_6; s0 += (int32_t)w2 * y1_7; s0 += (int32_t)w3 * y1_8;
                pW -= 8;
                s0 += rnd; s0 = (7 + rsh) > 0 ? s0 >> (7 + rsh) : s0 << -(7 + rsh); s0 = MIN(127, MAX(s0, -128));
                *pZ0++ = (int8_t)s0;
                y0_5 = y0_6; y0_6 = y0_7; y0_7 = y0_8;
                y1_5 = y1_6; y1_6 = y1_7; y1_7 = y1_8;
            }
            for (; j < Q + N - 1; ++j)
            {
                y0_8 = 0;
                y1_8 = 0;
                w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;
                s0 = (int32_t)w0 * y0_5; s0 += (int32_t)w1 * y0_6; s0 += (int32_t)w2 * y0_7; s0 += (int32_t)w3 * y0_8;
                w0 = *pW++; w1 = *pW++; w2 = *pW++; w3 = *pW++;
                s0 += (int32_t)w0 * y1_5; s0 += (int32_t)w1 * y1_6; s0 += (int32_t)w2 * y1_7; s0 += (int32_t)w3 * y1_8;
                pW -= 8;
                s0 += rnd; s0 = (7 + rsh) > 0 ? s0 >> (7 + rsh) : s0 << -(7 + rsh); s0 = MIN(127, MAX(s0, -128));
                *pZ0++ = (int8_t)s0;
                y0_5 = y0_6; y0_6 = y0_7; y0_7 = y0_8;
                y1_5 = y1_6; y1_6 = y1_7; y1_7 = y1_8;
            }
        }

    }
#undef M
#undef N
}

// scratch allocatation functions. return required scratch size in bytes
size_t conv2d_gen_3x3_8x8_getScratchSize(int P, int Q)
{
    const int N = 3;
    return 4 * 4 * sizeof(int8_t);
}

#endif
