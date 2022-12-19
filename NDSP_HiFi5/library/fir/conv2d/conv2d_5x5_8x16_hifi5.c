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
    2D Convolution, 8x16-bit
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

#if defined(AE_MUL8Q8X16CNV) & defined(AE_SAV16X4X2_XP)
void conv2d_5x5_8x16(void* pScr, int16_t* z, const int8_t* x, const int16_t* y, int rsh, int P, int Q)
{
#define M 5
#define N 5
    const ae_int16x8* restrict pY;
    const ae_int16x8* restrict pY0;
    const ae_int32x4* restrict pT_read;
    ae_int8x8* restrict pW;
    ae_int32x4* restrict pT;
    ae_int16x8* restrict pZ0;


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

        ae_int32x2 S0, S1, S2, S3;
        ae_int32x2 rnd0, rnd1, rnd2, rnd3, rnd;
        ae_int16x4 T0, T1;
        const ae_int8* pX = (const ae_int8*)x;

        ae_valignx2 alZ;
        ae_int16x4 Y00, Y01, Y02;

        int i, j, mmin, mmax, xstart, ystart, m;

        ae_int8x8 W0, W1, W2, W3, W4;

        alZ = AE_ZALIGN128();
        rnd = (7 + rsh) > 0 ? (1 << (7 + rsh)) >> 1 : 0;
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
            pY0 = (ae_int16x8*)(y + ystart * Q);
            pZ0 = (ae_int16x8*)(z + i * (Q + N - 1));

            if (mmax - mmin != 1)
            {
                pY = (ae_int16x8*)((ae_int16*)pY0 + Q);
                // first iteration
                {
                    pT_read = pT = (ae_int32x4*)pScr + 3;
                    Y00 = 0;

                    AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));

                    for (j = 0; j < Q; j += 8)
                    {
                        AE_L16X4X2_IP(Y01, Y02, pY0, sizeof(ae_int16x8));

                        AE_MUL8Q8X16CNV(S0, S1, W0, Y00, Y01, Y01);
                        AE_MUL8Q8X16CNV(S2, S3, W0, Y01, Y02, Y02);

                        Y00 = Y02;

                        AE_S32X2X2_IP(S0, S1, pT, sizeof(ae_int32x4));
                        AE_S32X2X2_IP(S2, S3, pT, sizeof(ae_int32x4));

                    }
                    //last N-1 elements
                    {
                        Y01 = Y02 = 0;

                        AE_MUL8Q8X16CNV(S0, S1, W0, Y00, Y01, Y01);
                        AE_MUL8Q8X16CNV(S2, S3, W0, Y01, Y02, Y02);

                        AE_S32X2X2_IP(S0, S1, pT, sizeof(ae_int32x4));
                        AE_S32X2X2_IP(S2, S3, pT, sizeof(ae_int32x4));
                    }
                }
                pY0 = pY;
                pY = (ae_int16x8*)((ae_int16*)pY + Q);
                __Pragma("loop_count min=0, max=2");
                for (m = mmin + 1; m < mmax - 1; ++m)
                {
                    pT_read = pT = (ae_int32x4*)pScr + 3;

                    Y00 = 0;
                    AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));

                    for (j = 0; j < Q; j += 8)
                    {
                        AE_L32X2X2_IP(S0, S1, pT_read, sizeof(ae_int32x4));
                        AE_L32X2X2_IP(S2, S3, pT_read, sizeof(ae_int32x4));

                        AE_L16X4X2_IP(Y01, Y02, pY0, sizeof(ae_int16x8));

                        AE_MULA8Q8X16CNV(S0, S1, W0, Y00, Y01, Y01);
                        AE_MULA8Q8X16CNV(S2, S3, W0, Y01, Y02, Y02);

                        Y00 = Y02;

                        AE_S32X2X2_IP(S0, S1, pT, sizeof(ae_int32x4));
                        AE_S32X2X2_IP(S2, S3, pT, sizeof(ae_int32x4));

                    }
                    //last N-1 elements
                    {
                        AE_L32X2X2_IP(S0, S1, pT_read, sizeof(ae_int32x4));
                        AE_L32X2X2_IP(S2, S3, pT_read, sizeof(ae_int32x4));

                        Y01 = Y02 = 0;

                        AE_MULA8Q8X16CNV(S0, S1, W0, Y00, Y01, Y01);
                        AE_MULA8Q8X16CNV(S2, S3, W0, Y01, Y02, Y02);

                        AE_S32X2X2_IP(S0, S1, pT, sizeof(ae_int32x4));
                        AE_S32X2X2_IP(S2, S3, pT, sizeof(ae_int32x4));
                    }

                    pY0 = pY;
                    pY = (ae_int16x8*)((ae_int16*)pY + Q);
                }
                //last iteration
                {
                    pT_read = pT = (ae_int32x4*)pScr + 3;

                    Y00 = 0;

                    AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));

                    for (j = 0; j < Q; j += 8)
                    {
                        AE_L32X2X2_IP(S0, S1, pT_read, sizeof(ae_int32x4));
                        AE_L32X2X2_IP(S2, S3, pT_read, sizeof(ae_int32x4));

                        AE_L16X4X2_IP(Y01, Y02, pY0, sizeof(ae_int16x8));

                        AE_MULA8Q8X16CNV(S0, S1, W0, Y00, Y01, Y01);
                        AE_MULA8Q8X16CNV(S2, S3, W0, Y01, Y02, Y02);

                        Y00 = Y02;

                        rnd0 = AE_MOVNEG32S_T(rnd, S0);
                        rnd1 = AE_MOVNEG32S_T(rnd, S1);
                        rnd2 = AE_MOVNEG32S_T(rnd, S2);
                        rnd3 = AE_MOVNEG32S_T(rnd, S3);

                        S0 = AE_ADD32S(S0, rnd0);
                        S1 = AE_ADD32S(S1, rnd1);
                        S2 = AE_ADD32S(S2, rnd2);
                        S3 = AE_ADD32S(S3, rnd3);

                        T0 = AE_TRUNCA16X4F32S(S0, S1, 9 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 9 - rsh);

                        AE_SA16X4X2_IP(T0, T1, alZ, pZ0);
                    }
                    //last N-1 elements
                    {
                        AE_L32X2X2_IP(S0, S1, pT_read, sizeof(ae_int32x4));
                        AE_L32X2X2_IP(S2, S3, pT_read, sizeof(ae_int32x4));

                        Y01 = Y02 = 0;

                        AE_MULA8Q8X16CNV(S0, S1, W0, Y00, Y01, Y01);
                        AE_MULA8Q8X16CNV(S2, S3, W0, Y01, Y02, Y02);

                        rnd0 = AE_MOVNEG32S_T(rnd, S0);
                        rnd1 = AE_MOVNEG32S_T(rnd, S1);
                        rnd2 = AE_MOVNEG32S_T(rnd, S2);
                        rnd3 = AE_MOVNEG32S_T(rnd, S3);

                        S0 = AE_ADD32S(S0, rnd0);
                        S1 = AE_ADD32S(S1, rnd1);
                        S2 = AE_ADD32S(S2, rnd2);
                        S3 = AE_ADD32S(S3, rnd3);

                        T0 = AE_TRUNCA16X4F32S(S0, S1, 9 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 9 - rsh);

                        AE_SAV16X4X2_XP(T0, T1, alZ, pZ0, 4 * sizeof(ae_int16));
                    }
                    AE_SA128POS_FP(alZ, pZ0);
                }
            }
            else
            {
                Y00 = 0;

                AE_L8X8_IP(W0, pW, sizeof(ae_int8x8));

                for (j = 0; j < Q; j += 8)
                {
                    AE_L16X4X2_IP(Y01, Y02, pY0, sizeof(ae_int16x8));

                    AE_MUL8Q8X16CNV(S0, S1, W0, Y00, Y01, Y01);
                    AE_MUL8Q8X16CNV(S2, S3, W0, Y01, Y02, Y02);

                    Y00 = Y02;

                    rnd0 = AE_MOVNEG32S_T(rnd, S0);
                    rnd1 = AE_MOVNEG32S_T(rnd, S1);
                    rnd2 = AE_MOVNEG32S_T(rnd, S2);
                    rnd3 = AE_MOVNEG32S_T(rnd, S3);

                    S0 = AE_ADD32S(S0, rnd0);
                    S1 = AE_ADD32S(S1, rnd1);
                    S2 = AE_ADD32S(S2, rnd2);
                    S3 = AE_ADD32S(S3, rnd3);

                    T0 = AE_TRUNCA16X4F32S(S0, S1, 9 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 9 - rsh);

                    AE_SA16X4X2_IP(T0, T1, alZ, pZ0);
                }
                //last N-1 elements
                {

                    Y01 = Y02 = 0;

                    AE_MUL8Q8X16CNV(S0, S1, W0, Y00, Y01, Y01);
                    AE_MUL8Q8X16CNV(S2, S3, W0, Y01, Y02, Y02);

                    rnd0 = AE_MOVNEG32S_T(rnd, S0);
                    rnd1 = AE_MOVNEG32S_T(rnd, S1);
                    rnd2 = AE_MOVNEG32S_T(rnd, S2);
                    rnd3 = AE_MOVNEG32S_T(rnd, S3);

                    S0 = AE_ADD32S(S0, rnd0);
                    S1 = AE_ADD32S(S1, rnd1);
                    S2 = AE_ADD32S(S2, rnd2);
                    S3 = AE_ADD32S(S3, rnd3);

                    T0 = AE_TRUNCA16X4F32S(S0, S1, 9 - rsh); T1 = AE_TRUNCA16X4F32S(S2, S3, 9 - rsh);

                    AE_SAV16X4X2_XP(T0, T1, alZ, pZ0, 4 * sizeof(ae_int16));
                }
                AE_SA128POS_FP(alZ, pZ0);
            }
        }
    }
    else
    {
        pW = (ae_int8x8*)pScr;
        ae_int16x4 t0, t1;
        ae_int8x8 t8x8;

        /* just copy original 8-bit coefficients to 16-bit shifting left by 8 bits */
        AE_L8X8_IP(t8x8, castxcc(ae_int8x8, x), 8 * sizeof(int8_t));
        t0 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(t8x8, AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 20));
        t1 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(t8x8, AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 21));
        AE_S16X4X2_IP(t0, t1, castxcc(ae_int16x8,pW), 8 * sizeof(int16_t));
        AE_L8X8_IP(t8x8, castxcc(ae_int8x8, x), 8 * sizeof(int8_t));
        t0 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(t8x8, AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 20));
        t1 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(t8x8, AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 21));
        AE_S16X4X2_IP(t0, t1, castxcc(ae_int16x8, pW), 8 * sizeof(int16_t));
        AE_L8X8_IP(t8x8, castxcc(ae_int8x8, x), 8 * sizeof(int8_t));
        t0 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(t8x8, AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 20));
        t1 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(t8x8, AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 21));
        AE_S16X4X2_IP(t0, t1, castxcc(ae_int16x8, pW), 8 * sizeof(int16_t));
        t0 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(AE_L8_I((ae_int8*)x, 0), AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 21));
        AE_S16_0_I(t0, (ae_int16*)pW, 0);

        conv2d_5x5_16x16((int16_t*)pScr + 32, z, (const int16_t*)pScr, y, rsh, P, Q);
    }
#undef M
#undef N
}
size_t conv2d_5x5_8x16_getScratchSize(int P, int Q)
{
    if (Q>= 128)
        return 6 * 8 * sizeof(int8_t) + (Q + 16) * sizeof(int32_t);
    else
        return 32 * sizeof(int16_t) + conv2d_5x5_16x16_getScratchSize(P, Q);
}

#else

void conv2d_5x5_8x16(void * pScr, int16_t * z, const int8_t * x, const int16_t * y, int rsh, int P, int Q)
{
    ae_int16x8 * restrict pW = (ae_int16x8 *)pScr;
    ae_int16x4 t0, t1;
    ae_int8x8 t8x8;

    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT_ALIGN(z, 16);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0 && P % 8 == 0);
    NASSERT(Q >= 0 && Q % 8 == 0);
    if (P <= 0 || Q <= 0) return;

    /* just copy original 8-bit coefficients to 16-bit shifting left by 8 bits */
    AE_L8X8_IP(t8x8, castxcc(ae_int8x8, x), 8 * sizeof(int8_t));
    t0 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(t8x8, AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 20));
    t1 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(t8x8, AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 21));
    AE_S16X4X2_IP(t0, t1, pW, 8 * sizeof(int16_t));
    AE_L8X8_IP(t8x8, castxcc(ae_int8x8, x), 8 * sizeof(int8_t));
    t0 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(t8x8, AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 20));
    t1 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(t8x8, AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 21));
    AE_S16X4X2_IP(t0, t1, pW, 8 * sizeof(int16_t));
    AE_L8X8_IP(t8x8, castxcc(ae_int8x8, x), 8 * sizeof(int8_t));
    t0 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(t8x8, AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 20));
    t1 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(t8x8, AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 21));
    AE_S16X4X2_IP(t0, t1, pW, 8 * sizeof(int16_t));
    t0 = AE_MOVINT16X4_FROMINT8X8(AE_SEL8X8I(AE_L8_I((ae_int8 *)x, 0), AE_MOVINT8X8_FROMINT64(AE_ZERO64()), 21));
    AE_S16_0_I(t0, (ae_int16 *)pW, 0);

    conv2d_5x5_16x16((int16_t*)pScr + 32, z, (const int16_t*)pScr, y, rsh, P, Q);
}

// scratch allocatation functions. return required scratch size in bytes
size_t conv2d_5x5_8x16_getScratchSize(int P, int Q)
{
    return 32 * sizeof(int16_t) + conv2d_5x5_16x16_getScratchSize(P, Q);
}

#endif
