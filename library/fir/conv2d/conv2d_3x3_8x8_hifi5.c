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


#if defined(AE_MUL4O8X8CNV_H) // NN extension

void conv2d_3x3_8x8(void* pScr, int8_t* z, const int8_t* x, const int8_t* y, int rsh, int P, int Q)
{
#define M 3
#define N 3
    const ae_int8x8* restrict pY0;
    const ae_int8x8* restrict pY1;
    const ae_int8x8* restrict pY2;
    const ae_int8x8* restrict pY3;
    ae_int8x8* restrict pW;
    ae_int8x8* restrict pZ0;
    ae_int8x8* restrict pZ1;

    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT_ALIGN(z, 16);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0 && P % 8 == 0);
    NASSERT(Q >= 0 && Q % 8 == 0);
    if (P <= 0 || Q <= 0) return;

    if (Q >= 128)
    {

        ae_int8x8 v00, v01, v10, v11, v20, v21;
        ae_int32x2 S0, S1, S2, S3;
        ae_int16x4 T0, T1;
        ae_int8x8 R0;
        const ae_int8* pX = (const ae_int8*)x;

        ae_valign alZ;

        int i, j, mmin, mmax, xstart, ystart;

        ae_int8x8 Y00, Y01;
        ae_int8x8 Y10, Y11;
        ae_int8x8 Y20, Y21;

        ae_int8x8 W0, W1, W2;

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

            AE_S8X8_IP(W0, pW, sizeof(ae_int8x8));
            AE_S8X8_IP(W1, pW, sizeof(ae_int8x8));
            AE_S8X8_IP(W2, pW, sizeof(ae_int8x8));
        }

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

                AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));
                AE_L8X8_IP(W1, pW, sizeof(ae_int8x8));
                AE_L8X8_IP(W2, pW, sizeof(ae_int8x8));

                for (j = 0; j < Q; j += 8)
                {
                    AE_L8X8_IP(Y01, pY0, sizeof(ae_int8x8));
                    AE_L8X8_IP(Y11, pY1, sizeof(ae_int8x8));
                    AE_L8X8_IP(Y21, pY2, sizeof(ae_int8x8));

                    AE_DSEL8X8(v00, v01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                    AE_DSEL8X8(v10, v11, Y10, Y11, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                    AE_DSEL8X8(v20, v21, Y20, Y21, AE_MOVINT8X8_FROMINT64(0xA291807060504030));

                    AE_MOVD8X16(Y00, Y10, Y01, Y11); Y20 = Y21;
                    AE_MUL4O8X8CNV_H(S0, S1, S2, S3, W0, v00, v01);
                    AE_MULA4O8X8CNV_H(S0, S1, S2, S3, W1, v10, v11);
                    AE_MULA4O8X8CNV_H(S0, S1, S2, S3, W2, v20, v21);

                    T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                    AE_SA8X8_IP(AE_ROUND8X8F16SSYM(T0, T1), alZ, pZ0);
                }
                AE_SA64POS_FP(alZ, pZ0);
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

                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));

                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
            // if only one Y row need to 
            else if (mmax - mmin == 1)
            {
                Y00 = AE_MOVINT8X8_FROMINT16(0);
                AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));
                for (j = 0; j < Q; j += 8)
                {
                    AE_L8X8_IP(Y01, pY0, sizeof(ae_int8x8));

                    AE_DSEL8X8(v00, v01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                    Y00 = Y01;
                    AE_MUL4O8X8CNV_H(S0, S1, S2, S3, W0, v00, v01);

                    T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                    AE_SA8X8_IP(AE_ROUND8X8F16SSYM(T0, T1), alZ, pZ0);
                }
                AE_SA64POS_FP(alZ, pZ0);
                //last N-1 elements
                {
                    Y01 = AE_MOVINT8X8_FROMINT16(0);
                    AE_DSEL8X8(v00, v01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                    AE_MUL4O8X8CNV_H(S0, S1, S2, S3, W0, v00, v01);
                    T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                    R0 = AE_ROUND8X8F16SSYM(T0, T1);
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
            else if (mmax - mmin == 2)
            {
                Y00 = AE_MOVINT8X8_FROMINT16(0);
                Y10 = AE_MOVINT8X8_FROMINT16(0);

                AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));
                AE_L8X8_IP(W1, pW, sizeof(ae_int8x8));

                for (j = 0; j < Q; j += 8)
                {
                    AE_L8X8_IP(Y01, pY0, sizeof(ae_int8x8));
                    AE_L8X8_IP(Y11, pY1, sizeof(ae_int8x8));

                    AE_DSEL8X8(v00, v01, Y00, Y01, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                    AE_DSEL8X8(v10, v11, Y10, Y11, AE_MOVINT8X8_FROMINT64(0xA291807060504030));
                    AE_MOVD8X16(Y00, Y10, Y01, Y11);
                    AE_MUL4O8X8CNV_H(S0, S1, S2, S3, W0, v00, v01);
                    AE_MULA4O8X8CNV_H(S0, S1, S2, S3, W1, v10, v11);

                    T0 = AE_TRUNCA16X4F32S(S0, S1, 17 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 17 - rsh);
                    AE_SA8X8_IP(AE_ROUND8X8F16SSYM(T0, T1), alZ, pZ0);
                }
                AE_SA64POS_FP(alZ, pZ0);
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

                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                    R0 = AE_SEL8X8(R0, R0, AE_MOVINT8X8_FROMINT64(0x0E0D0C0B0A090807));
                    AE_S8_0_IP(R0, castxcc(ae_int8, pZ0), sizeof(ae_int8));
                }
            }
        }
    }
    else
    {
        int i, j;
        ae_valign aZ0, aZ1;
        ae_int16x4 w0, w1, w2;
        ae_int16x4 y00, y01, y02, y10, y11, y12;
        ae_int16x4 y20, y21, y22, y30, y31, y32;
        ae_int64 S0, S1, S2, S3;
        ae_int16x4 r0, r1;



        /* Preload coefficients */
        {
            ae_int16x4 w0_, w1_, w2_;
            ae_int64 tmp0;

            pW = (ae_int8x8*)x;

            AE_L64_IP(tmp0, castxcc(ae_int64, pW), 8 * sizeof(int8_t));
            w2_ = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(AE_MOVINT8X8_FROMINT64(tmp0), AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 21));
            w1_ = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(AE_MOVINT8X8_FROMINT64(tmp0), AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 20));
            w0_ = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(AE_L8_I((ae_int8*)pW, 0), AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 21));

            w2 = AE_SEL16I(w2_, AE_ZERO16(), 5);
            w2 = AE_SEL16I(AE_ZERO16(), w2, 4);

            w1 = AE_SEL16I(w1_, w2_, 0);
            w1 = AE_SEL16I(AE_ZERO16(), w1, 4);

            w0 = AE_SEL16I(w0_, w1_, 4);
            w0 = AE_SEL16I(AE_ZERO16(), w0, 4);
        }

        pZ0 = (ae_int8x8*)z;
        pZ1 = (ae_int8x8*)(z + Q + N - 1);
        aZ0 = AE_ZALIGN64();
        aZ1 = AE_ZALIGN64();

        pY0 = (ae_int8x8*)y;
        pY1 = (ae_int8x8*)((ae_int8*)pY0 + Q);
        pY2 = (ae_int8x8*)((ae_int8*)pY1 + Q);
        pY3 = (ae_int8x8*)((ae_int8*)pY2 + Q);

        {
            y00 = AE_ZERO16();
            y10 = AE_ZERO16();

            /* Process all but last N-1 samples of the 0 and 1 rows */
            __Pragma("loop_count min=1");
            for (j = 0; j < (Q >> 3); j++)
            {
                AE_L8X4F_IP(y01, castxcc(int8_t,pY0), 4 * sizeof(int8_t));
                AE_L8X4F_IP(y11, castxcc(int8_t,pY1), 4 * sizeof(int8_t));
                AE_L8X4F_IP(y02, castxcc(int8_t,pY0), 4 * sizeof(int8_t));
                AE_L8X4F_IP(y12, castxcc(int8_t,pY1), 4 * sizeof(int8_t));

                AE_MULFQ16X2_FIR_2(S0, S1, y00, y01, w2);
                AE_MULFQ16X2_FIR_0(S2, S3, y00, y01, w2);
                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, y01, y02, w2);
                AE_MULFQ16X2_FIR_0(S2, S3, y01, y02, w2);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), aZ0, pZ0);


                AE_MULFQ16X2_FIR_2(S0, S1, y00, y01, w1);
                AE_MULFQ16X2_FIR_0(S2, S3, y00, y01, w1);
                AE_MULAFQ16X2_FIR_2(S0, S1, y10, y11, w2);
                AE_MULAFQ16X2_FIR_0(S2, S3, y10, y11, w2);
                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, y01, y02, w1);
                AE_MULFQ16X2_FIR_0(S2, S3, y01, y02, w1);
                AE_MULAFQ16X2_FIR_2(S0, S1, y11, y12, w2);
                AE_MULAFQ16X2_FIR_0(S2, S3, y11, y12, w2);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), aZ1, pZ1);

                y00 = y02;
                y10 = y12;
            }
            AE_SA64POS_FP(aZ0, pZ0);
            AE_SA64POS_FP(aZ1, pZ1);

            /* Last N-1 samples */
            {
                AE_MULFQ16X2_FIR_2(S0, S1, y00, AE_ZERO16(), w2);

                AE_MULFQ16X2_FIR_2(S2, S3, y00, AE_ZERO16(), w1);
                AE_MULAFQ16X2_FIR_2(S2, S3, y10, AE_ZERO16(), w2);

                r0 = AE_MOVINT16X4_FROMINT8X8(AE_ROUND8X4F32SSYM_L(AE_TRUNCA32X2F64S(S3, S2, 32 - rsh), AE_TRUNCA32X2F64S(S1, S0, 32 - rsh)));

                AE_S16_0_XP(r0, castxcc(ae_int16, pZ0), (Q + (N - 1) + 2) * sizeof(int8_t));
                AE_S16_0_XP(AE_SEL16_4321(r0, r0), castxcc(ae_int16, pZ1), (Q + (N - 1) + 2) * sizeof(int8_t));
            }

            pY0 = (ae_int8x8*)((ae_int8*)pY0 - Q);
            pY1 = (ae_int8x8*)((ae_int8*)pY1 - Q);
        }

        __Pragma("loop_count min=1");
        for (i = (M - 1); i < (P); i += 2)
        {
            y00 = AE_ZERO16();
            y10 = AE_ZERO16();
            y20 = AE_ZERO16();
            y30 = AE_ZERO16();

            /* Process all but last N-1 samples of the i and i+1 rows */
            __Pragma("loop_count min=1");
            for (j = 0; j < (Q >> 3); j++)
            {
                AE_L8X4F_IP(y01, castxcc(int8_t,pY0), 4 * sizeof(int8_t));
                AE_L8X4F_IP(y11, castxcc(int8_t,pY1), 4 * sizeof(int8_t));
                AE_L8X4F_IP(y21, castxcc(int8_t,pY2), 4 * sizeof(int8_t));
                AE_L8X4F_IP(y31, castxcc(int8_t,pY3), 4 * sizeof(int8_t));
                AE_L8X4F_IP(y02, castxcc(int8_t,pY0), 4 * sizeof(int8_t));
                AE_L8X4F_IP(y12, castxcc(int8_t,pY1), 4 * sizeof(int8_t));
                AE_L8X4F_IP(y22, castxcc(int8_t,pY2), 4 * sizeof(int8_t));
                AE_L8X4F_IP(y32, castxcc(int8_t,pY3), 4 * sizeof(int8_t));

                AE_MULFQ16X2_FIR_2(S0, S1, y00, y01, w0);
                AE_MULFQ16X2_FIR_0(S2, S3, y00, y01, w0);
                AE_MULAFQ16X2_FIR_2(S0, S1, y10, y11, w1);
                AE_MULAFQ16X2_FIR_0(S2, S3, y10, y11, w1);
                AE_MULAFQ16X2_FIR_2(S0, S1, y20, y21, w2);
                AE_MULAFQ16X2_FIR_0(S2, S3, y20, y21, w2);
                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, y01, y02, w0);
                AE_MULFQ16X2_FIR_0(S2, S3, y01, y02, w0);
                AE_MULAFQ16X2_FIR_2(S0, S1, y11, y12, w1);
                AE_MULAFQ16X2_FIR_0(S2, S3, y11, y12, w1);
                AE_MULAFQ16X2_FIR_2(S0, S1, y21, y22, w2);
                AE_MULAFQ16X2_FIR_0(S2, S3, y21, y22, w2);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), aZ0, pZ0);


                AE_MULFQ16X2_FIR_2(S0, S1, y10, y11, w0);
                AE_MULFQ16X2_FIR_0(S2, S3, y10, y11, w0);
                AE_MULAFQ16X2_FIR_2(S0, S1, y20, y21, w1);
                AE_MULAFQ16X2_FIR_0(S2, S3, y20, y21, w1);
                AE_MULAFQ16X2_FIR_2(S0, S1, y30, y31, w2);
                AE_MULAFQ16X2_FIR_0(S2, S3, y30, y31, w2);
                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, y11, y12, w0);
                AE_MULFQ16X2_FIR_0(S2, S3, y11, y12, w0);
                AE_MULAFQ16X2_FIR_2(S0, S1, y21, y22, w1);
                AE_MULAFQ16X2_FIR_0(S2, S3, y21, y22, w1);
                AE_MULAFQ16X2_FIR_2(S0, S1, y31, y32, w2);
                AE_MULAFQ16X2_FIR_0(S2, S3, y31, y32, w2);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), aZ1, pZ1);

                y00 = y02;
                y10 = y12;
                y20 = y22;
                y30 = y32;
            }
            AE_SA64POS_FP(aZ0, pZ0);
            AE_SA64POS_FP(aZ1, pZ1);

            /* Last N-1 samples */
            {
                AE_MULFQ16X2_FIR_2(S0, S1, y00, AE_ZERO16(), w0);
                AE_MULAFQ16X2_FIR_2(S0, S1, y10, AE_ZERO16(), w1);
                AE_MULAFQ16X2_FIR_2(S0, S1, y20, AE_ZERO16(), w2);

                AE_MULFQ16X2_FIR_2(S2, S3, y10, AE_ZERO16(), w0);
                AE_MULAFQ16X2_FIR_2(S2, S3, y20, AE_ZERO16(), w1);
                AE_MULAFQ16X2_FIR_2(S2, S3, y30, AE_ZERO16(), w2);

                r0 = AE_MOVINT16X4_FROMINT8X8(AE_ROUND8X4F32SSYM_L(AE_TRUNCA32X2F64S(S3, S2, 32 - rsh), AE_TRUNCA32X2F64S(S1, S0, 32 - rsh)));

                AE_S16_0_XP(r0, castxcc(ae_int16, pZ0), (Q + (N - 1) + 2) * sizeof(int8_t));
                AE_S16_0_XP(AE_SEL16_4321(r0, r0), castxcc(ae_int16, pZ1), (Q + (N - 1) + 2) * sizeof(int8_t));
            }

            pY0 = (ae_int8x8*)((ae_int8*)pY0 + Q);
            pY1 = (ae_int8x8*)((ae_int8*)pY1 + Q);
            pY2 = (ae_int8x8*)((ae_int8*)pY2 + Q);
            pY3 = (ae_int8x8*)((ae_int8*)pY3 + Q);
        }

        {
            y10 = AE_ZERO16();
            y20 = AE_ZERO16();

            /* Process all but last N-1 samples of the (P+M-3) and (P+M-2) rows */
            __Pragma("loop_count min=1");
            for (j = 0; j < (Q >> 3); j++)
            {
                AE_L8X4F_IP(y11, castxcc(int8_t,pY0), 4 * sizeof(int8_t));
                AE_L8X4F_IP(y21, castxcc(int8_t,pY1), 4 * sizeof(int8_t));
                AE_L8X4F_IP(y12, castxcc(int8_t,pY0), 4 * sizeof(int8_t));
                AE_L8X4F_IP(y22, castxcc(int8_t,pY1), 4 * sizeof(int8_t));

                AE_MULFQ16X2_FIR_2(S0, S1, y10, y11, w0);
                AE_MULFQ16X2_FIR_0(S2, S3, y10, y11, w0);
                AE_MULAFQ16X2_FIR_2(S0, S1, y20, y21, w1);
                AE_MULAFQ16X2_FIR_0(S2, S3, y20, y21, w1);
                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, y11, y12, w0);
                AE_MULFQ16X2_FIR_0(S2, S3, y11, y12, w0);
                AE_MULAFQ16X2_FIR_2(S0, S1, y21, y22, w1);
                AE_MULAFQ16X2_FIR_0(S2, S3, y21, y22, w1);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), aZ0, pZ0);


                AE_MULFQ16X2_FIR_2(S0, S1, y20, y21, w0);
                AE_MULFQ16X2_FIR_0(S2, S3, y20, y21, w0);
                r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));

                AE_MULFQ16X2_FIR_2(S0, S1, y21, y22, w0);
                AE_MULFQ16X2_FIR_0(S2, S3, y21, y22, w0);
                r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), aZ1, pZ1);

                y10 = y12;
                y20 = y22;
            }
            AE_SA64POS_FP(aZ0, pZ0);
            AE_SA64POS_FP(aZ1, pZ1);

            /* Last N-1 samples */
            {
                AE_MULFQ16X2_FIR_2(S0, S1, y10, AE_ZERO16(), w0);
                AE_MULAFQ16X2_FIR_2(S0, S1, y20, AE_ZERO16(), w1);

                AE_MULFQ16X2_FIR_2(S2, S3, y20, AE_ZERO16(), w0);

                r0 = AE_MOVINT16X4_FROMINT8X8(AE_ROUND8X4F32SSYM_L(AE_TRUNCA32X2F64S(S3, S2, 32 - rsh), AE_TRUNCA32X2F64S(S1, S0, 32 - rsh)));

                AE_S16_0_XP(r0, castxcc(ae_int16, pZ0), (Q + (N - 1) + 2) * sizeof(int8_t));
                AE_S16_0_XP(AE_SEL16_4321(r0, r0), castxcc(ae_int16, pZ1), (Q + (N - 1) + 2) * sizeof(int8_t));
            }
        }

    }
#undef M
#undef N
}

// scratch allocatation functions. return required scratch size in bytes
size_t conv2d_3x3_8x8_getScratchSize(int P, int Q)
{
    if (Q >= 128)
        return 4 * sizeof(ae_int8x8);
    else
        return 0;
}

#else // no NN extension

void conv2d_3x3_8x8(void * pScr, int8_t * z, const int8_t * x, const int8_t * y, int rsh, int P, int Q)
{
#define M 3
#define N 3
    const int8_t    * restrict pY0;
    const int8_t    * restrict pY1;
    const int8_t    * restrict pY2;
    const int8_t    * restrict pY3;
          ae_int8x8 * restrict pZ0;
          ae_int8x8 * restrict pZ1;
    const ae_int64  * restrict pW;
    int i,j;
    ae_valign aZ0, aZ1;
    ae_int16x4 w0, w1, w2;
    ae_int16x4 y00, y01, y02, y10, y11, y12;
    ae_int16x4 y20, y21, y22, y30, y31, y32;
    ae_int64 S0, S1, S2, S3;
    ae_int16x4 r0, r1;

    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT_ALIGN(z, 16);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0 && P % 8 == 0);
    NASSERT(Q >= 0 && Q % 8 == 0);
    if (P <= 0 || Q <= 0) return;

    /* Preload coefficients */
    {
        ae_int16x4 w0_, w1_, w2_;
        ae_int64 tmp0;

        pW = (const ae_int64 *)x;

        AE_L64_IP(tmp0, pW, 8 * sizeof(int8_t));
        w2_ = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(AE_MOVINT8X8_FROMINT64(tmp0), AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 21));
        w1_ = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(AE_MOVINT8X8_FROMINT64(tmp0), AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 20));
        w0_ = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(AE_L8_I((ae_int8 *)pW, 0), AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 21));

        w2 = AE_SEL16I(w2_, AE_ZERO16(), 5);
        w2 = AE_SEL16I(AE_ZERO16(), w2, 4);

        w1 = AE_SEL16I(w1_, w2_, 0);
        w1 = AE_SEL16I(AE_ZERO16(), w1, 4);

        w0 = AE_SEL16I(w0_, w1_, 4);
        w0 = AE_SEL16I(AE_ZERO16(), w0, 4);
    }

    pZ0 = (ae_int8x8 *)z;
    pZ1 = (ae_int8x8 *)(z + Q + N - 1);
    aZ0 = AE_ZALIGN64();
    aZ1 = AE_ZALIGN64();

    pY0 = y;
    pY1 = pY0 + Q;
    pY2 = pY1 + Q;
    pY3 = pY2 + Q;

    {
        y00 = AE_ZERO16();
        y10 = AE_ZERO16();

        /* Process all but last N-1 samples of the 0 and 1 rows */
        __Pragma("loop_count min=1");
        for (j = 0; j < (Q >> 3); j++)
        {
            AE_L8X4F_IP(y01, pY0, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y11, pY1, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y02, pY0, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y12, pY1, 4 * sizeof(int8_t));

            AE_MULFQ16X2_FIR_2(S0, S1, y00, y01, w2);
            AE_MULFQ16X2_FIR_0(S2, S3, y00, y01, w2);
            r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));

            AE_MULFQ16X2_FIR_2(S0, S1, y01, y02, w2);
            AE_MULFQ16X2_FIR_0(S2, S3, y01, y02, w2);
            r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
            AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), aZ0, pZ0);


            AE_MULFQ16X2_FIR_2(S0, S1, y00, y01, w1);
            AE_MULFQ16X2_FIR_0(S2, S3, y00, y01, w1);
            AE_MULAFQ16X2_FIR_2(S0, S1, y10, y11, w2);
            AE_MULAFQ16X2_FIR_0(S2, S3, y10, y11, w2);
            r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));

            AE_MULFQ16X2_FIR_2(S0, S1, y01, y02, w1);
            AE_MULFQ16X2_FIR_0(S2, S3, y01, y02, w1);
            AE_MULAFQ16X2_FIR_2(S0, S1, y11, y12, w2);
            AE_MULAFQ16X2_FIR_0(S2, S3, y11, y12, w2);
            r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
            AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), aZ1, pZ1);

            y00 = y02;
            y10 = y12;
        }
        AE_SA64POS_FP(aZ0, pZ0);
        AE_SA64POS_FP(aZ1, pZ1);

        /* Last N-1 samples */
        {
            AE_MULFQ16X2_FIR_2(S0, S1, y00, AE_ZERO16(), w2);

            AE_MULFQ16X2_FIR_2(S2, S3, y00, AE_ZERO16(), w1);
            AE_MULAFQ16X2_FIR_2(S2, S3, y10, AE_ZERO16(), w2);

            r0 = AE_MOVINT16X4_FROMINT8X8(AE_ROUND8X4F32SSYM_L(AE_TRUNCA32X2F64S(S3, S2, 32 - rsh), AE_TRUNCA32X2F64S(S1, S0, 32 - rsh)));

            AE_S16_0_XP(r0, castxcc(ae_int16, pZ0), (Q + (N - 1) + 2) * sizeof(int8_t));
            AE_S16_0_XP(AE_SEL16_4321(r0, r0), castxcc(ae_int16, pZ1), (Q + (N - 1) + 2) * sizeof(int8_t));
        }

        pY0 -= Q;
        pY1 -= Q;
    }

    __Pragma("loop_count min=1");
    for (i = (M - 1); i < (P); i += 2)
    {
        y00 = AE_ZERO16();
        y10 = AE_ZERO16();
        y20 = AE_ZERO16();
        y30 = AE_ZERO16();

        /* Process all but last N-1 samples of the i and i+1 rows */
        __Pragma("loop_count min=1");
        for (j = 0; j < (Q >> 3); j++)
        {
            AE_L8X4F_IP(y01, pY0, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y11, pY1, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y21, pY2, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y31, pY3, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y02, pY0, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y12, pY1, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y22, pY2, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y32, pY3, 4 * sizeof(int8_t));

            AE_MULFQ16X2_FIR_2(S0, S1, y00, y01, w0);
            AE_MULFQ16X2_FIR_0(S2, S3, y00, y01, w0);
            AE_MULAFQ16X2_FIR_2(S0, S1, y10, y11, w1);
            AE_MULAFQ16X2_FIR_0(S2, S3, y10, y11, w1);
            AE_MULAFQ16X2_FIR_2(S0, S1, y20, y21, w2);
            AE_MULAFQ16X2_FIR_0(S2, S3, y20, y21, w2);
            r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));

            AE_MULFQ16X2_FIR_2(S0, S1, y01, y02, w0);
            AE_MULFQ16X2_FIR_0(S2, S3, y01, y02, w0);
            AE_MULAFQ16X2_FIR_2(S0, S1, y11, y12, w1);
            AE_MULAFQ16X2_FIR_0(S2, S3, y11, y12, w1);
            AE_MULAFQ16X2_FIR_2(S0, S1, y21, y22, w2);
            AE_MULAFQ16X2_FIR_0(S2, S3, y21, y22, w2);
            r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
            AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), aZ0, pZ0);


            AE_MULFQ16X2_FIR_2(S0, S1, y10, y11, w0);
            AE_MULFQ16X2_FIR_0(S2, S3, y10, y11, w0);
            AE_MULAFQ16X2_FIR_2(S0, S1, y20, y21, w1);
            AE_MULAFQ16X2_FIR_0(S2, S3, y20, y21, w1);
            AE_MULAFQ16X2_FIR_2(S0, S1, y30, y31, w2);
            AE_MULAFQ16X2_FIR_0(S2, S3, y30, y31, w2);
            r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));

            AE_MULFQ16X2_FIR_2(S0, S1, y11, y12, w0);
            AE_MULFQ16X2_FIR_0(S2, S3, y11, y12, w0);
            AE_MULAFQ16X2_FIR_2(S0, S1, y21, y22, w1);
            AE_MULAFQ16X2_FIR_0(S2, S3, y21, y22, w1);
            AE_MULAFQ16X2_FIR_2(S0, S1, y31, y32, w2);
            AE_MULAFQ16X2_FIR_0(S2, S3, y31, y32, w2);
            r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
            AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), aZ1, pZ1);

            y00 = y02;
            y10 = y12;
            y20 = y22;
            y30 = y32;
        }
        AE_SA64POS_FP(aZ0, pZ0);
        AE_SA64POS_FP(aZ1, pZ1);

        /* Last N-1 samples */
        {
            AE_MULFQ16X2_FIR_2(S0, S1, y00, AE_ZERO16(), w0);
            AE_MULAFQ16X2_FIR_2(S0, S1, y10, AE_ZERO16(), w1);
            AE_MULAFQ16X2_FIR_2(S0, S1, y20, AE_ZERO16(), w2);

            AE_MULFQ16X2_FIR_2(S2, S3, y10, AE_ZERO16(), w0);
            AE_MULAFQ16X2_FIR_2(S2, S3, y20, AE_ZERO16(), w1);
            AE_MULAFQ16X2_FIR_2(S2, S3, y30, AE_ZERO16(), w2);

            r0 = AE_MOVINT16X4_FROMINT8X8(AE_ROUND8X4F32SSYM_L(AE_TRUNCA32X2F64S(S3, S2, 32 - rsh), AE_TRUNCA32X2F64S(S1, S0, 32 - rsh)));

            AE_S16_0_XP(r0, castxcc(ae_int16, pZ0), (Q + (N - 1) + 2) * sizeof(int8_t));
            AE_S16_0_XP(AE_SEL16_4321(r0, r0), castxcc(ae_int16, pZ1), (Q + (N - 1) + 2) * sizeof(int8_t));
        }

        pY0 += Q;
        pY1 += Q;
        pY2 += Q;
        pY3 += Q;
    }

    {
        y10 = AE_ZERO16();
        y20 = AE_ZERO16();

        /* Process all but last N-1 samples of the (P+M-3) and (P+M-2) rows */
        __Pragma("loop_count min=1");
        for (j = 0; j < (Q >> 3); j++)
        {
            AE_L8X4F_IP(y11, pY0, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y21, pY1, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y12, pY0, 4 * sizeof(int8_t));
            AE_L8X4F_IP(y22, pY1, 4 * sizeof(int8_t));

            AE_MULFQ16X2_FIR_2(S0, S1, y10, y11, w0);
            AE_MULFQ16X2_FIR_0(S2, S3, y10, y11, w0);
            AE_MULAFQ16X2_FIR_2(S0, S1, y20, y21, w1);
            AE_MULAFQ16X2_FIR_0(S2, S3, y20, y21, w1);
            r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));

            AE_MULFQ16X2_FIR_2(S0, S1, y11, y12, w0);
            AE_MULFQ16X2_FIR_0(S2, S3, y11, y12, w0);
            AE_MULAFQ16X2_FIR_2(S0, S1, y21, y22, w1);
            AE_MULAFQ16X2_FIR_0(S2, S3, y21, y22, w1);
            r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
            AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), aZ0, pZ0);


            AE_MULFQ16X2_FIR_2(S0, S1, y20, y21, w0);
            AE_MULFQ16X2_FIR_0(S2, S3, y20, y21, w0);
            r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));

            AE_MULFQ16X2_FIR_2(S0, S1, y21, y22, w0);
            AE_MULFQ16X2_FIR_0(S2, S3, y21, y22, w0);
            r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
            AE_SA8X8_IP(AE_ROUND8X8F16SSYM(r0, r1), aZ1, pZ1);

            y10 = y12;
            y20 = y22;
        }
        AE_SA64POS_FP(aZ0, pZ0);
        AE_SA64POS_FP(aZ1, pZ1);

        /* Last N-1 samples */
        {
            AE_MULFQ16X2_FIR_2(S0, S1, y10, AE_ZERO16(), w0);
            AE_MULAFQ16X2_FIR_2(S0, S1, y20, AE_ZERO16(), w1);

            AE_MULFQ16X2_FIR_2(S2, S3, y20, AE_ZERO16(), w0);

            r0 = AE_MOVINT16X4_FROMINT8X8(AE_ROUND8X4F32SSYM_L(AE_TRUNCA32X2F64S(S3, S2, 32 - rsh), AE_TRUNCA32X2F64S(S1, S0, 32 - rsh)));

            AE_S16_0_XP(r0, castxcc(ae_int16, pZ0), (Q + (N - 1) + 2) * sizeof(int8_t));
            AE_S16_0_XP(AE_SEL16_4321(r0, r0), castxcc(ae_int16, pZ1), (Q + (N - 1) + 2) * sizeof(int8_t));
        }
    }
#undef M
#undef N
}

// scratch allocatation functions. return required scratch size in bytes
size_t conv2d_3x3_8x8_getScratchSize(int P, int Q)
{
    return 0;
}

#endif

