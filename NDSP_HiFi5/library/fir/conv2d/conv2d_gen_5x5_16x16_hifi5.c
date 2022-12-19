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

void conv2d_gen_5x5_16x16(void * pScr, int16_t * z, const int16_t * x, const int16_t * y, int rsh, int P, int Q)
{
#define M 5
#define N 5
    const ae_int16x8* restrict pY;
    const ae_int16x8* restrict pY0;
    const ae_int16x8* restrict pY1;
    const ae_int16x8* restrict pY2;
    const ae_int64x2* restrict pT_read;
    ae_int16x8* restrict pW;
    ae_int64x2* restrict pT;
    ae_int16x8* restrict pZ0;
    ae_int64 S0, S1, S2, S3, S4, S5, S6, S7;
    ae_int16x4 T0, T1;
    const ae_int16* pX = (const ae_int16*)x;

    ae_valignx2 al0, al1, al2, alZ;
    ae_int16x4 Y00, Y01, Y02;
    ae_int16x4 Y10, Y11, Y12;
    ae_int16x4 Y20, Y21, Y22;
    ae_int16x4 W00, W01, W10, W11, W20, W21;
    int i, j, mmin, mmax, xstart, ystart, m;

    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0);
    NASSERT(Q >= 0);
    if (P <= 0 || Q <= 0) return;



    alZ = AE_ZALIGN128();

    pW = (ae_int16x8*)pScr+4;
    // preload X and save in scratch
    {
        ae_int16x4 w0, w1, w2, w3, w4, w5;
        w0 = 0;
        for (m = 0; m < M; ++m)
        {
            AE_L16_IP(w5, pX, sizeof(ae_int16));
            AE_L16_IP(w4, pX, sizeof(ae_int16));
            AE_L16_IP(w3, pX, sizeof(ae_int16));
            AE_L16_IP(w2, pX, sizeof(ae_int16));
            AE_L16_IP(w1, pX, sizeof(ae_int16));

            W00 = w0;
            W00 = AE_SEL16_6543(W00, w1);
            W00 = AE_SEL16_6543(W00, w2);
            W00 = AE_SEL16_6543(W00, w3);
            W00 = AE_SEL16_6543(W00, w4);
            W01 = AE_SEL16_4321(w5, w0);

            AE_S16X4X2_IP(W00,W01, pW, -(int)sizeof(ae_int16x8));
        }
    }

    __Pragma("loop_count min=1");
    for (i = 0; i < P + M - 1; i += 1)
    {
        mmin = i < M - 1 ? M - 1 - i : 0;
        mmax = i >= P ? M - 1 + P - i : M;
        xstart = i < M - 1 ? (mmin) : 0;
        ystart = i < M - 1 ? 0 : i - (M - 1);

        pW = (ae_int16x8*)pScr + xstart;
        pY0 = (ae_int16x8*)(y + ystart * Q);
        pZ0 = (ae_int16x8*)(z + i * (Q + N - 1));

        if (mmax - mmin == 5)
        {
            pY = (ae_int16x8*)((ae_int16*)pY0 + 2*Q);
            // first 2 iterations
            {
                pT = (ae_int64x2*)pScr + 6;

                pY1 = (ae_int16x8*)((ae_int16*)pY0 + Q);

                Y00 = 0;
                Y10 = 0;
                al0 = AE_LA128_PP(pY0);
                al1 = AE_LA128_PP(pY1);

                AE_L16X4X2_IP(W00, W01, pW, sizeof(ae_int16x8));
                AE_L16X4X2_IP(W10, W11, pW, sizeof(ae_int16x8));
                for (j = 0; j < (Q & ~7); j += 8)
                {
                    AE_LA16X4X2_IP(Y01, Y02, al0, pY0);
                    AE_LA16X4X2_IP(Y11, Y12, al1, pY1);

                    AE_MULFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y10, Y11, W10);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y10, Y11, W10);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y11, Y12, W10);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y11, Y12, W10);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y11, Y11, W11);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y11, Y11, W11);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y12, Y12, W11);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y12, Y12, W11);

                    AE_MOVD16X8(Y00, Y10, Y02, Y12);

                    AE_S64X2_IP(S0, S1, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S2, S3, pT, sizeof(ae_int64x2));
                    AE_S64X2_IP(S4, S5, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S6, S7, pT, sizeof(ae_int64x2));

                }
                // tail
                if (j < Q)
                {
                    AE_LA16X4X2_IP(Y01, Y02, al0, pY0);
                    AE_LA16X4X2_IP(Y11, Y12, al1, pY1);

                    AE_MULFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y10, Y11, W10);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y10, Y11, W10);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y11, Y12, W10);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y11, Y12, W10);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y11, Y11, W11);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y11, Y11, W11);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y12, Y12, W11);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y12, Y12, W11);

                    AE_S64X2_IP(S0, S1, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S2, S3, pT, sizeof(ae_int64x2));
                    AE_S64X2_IP(S4, S5, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S6, S7, pT, sizeof(ae_int64x2));
                    for (; j < Q; ++j)
                    {
                        AE_DSEL16X4(Y00, Y01, Y00, Y01, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                        AE_DSEL16X4(Y01, Y02, Y01, Y02, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                        AE_DSEL16X4(Y10, Y11, Y10, Y11, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                        AE_DSEL16X4(Y11, Y12, Y11, Y12, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    }
                }
                //last N-1 elements
                {
                    Y01 = 0; Y02 = 0;
                    Y11 = 0; Y12 = 0;


                    AE_MULFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y10, Y11, W10);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y10, Y11, W10);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y11, Y12, W10);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y11, Y12, W10);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y11, Y11, W11);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y11, Y11, W11);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y12, Y12, W11);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y12, Y12, W11);

                    AE_S64X2_IP(S0, S1, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S2, S3, pT, sizeof(ae_int64x2));
                    AE_S64X2_IP(S4, S5, pT, sizeof(ae_int64x2)); AE_S64X2_IP(S6, S7, pT, sizeof(ae_int64x2));
                }
            }      
            //last 3 iteration
            {
                pT_read = (ae_int64x2*)pScr + 6;
                
                pY0 = pY;
                pY1 = (ae_int16x8*)((ae_int16*)pY + Q);
                pY2 = (ae_int16x8*)((ae_int16*)pY + 2 * Q);

                Y00 = 0;
                Y10 = 0;
                Y20 = 0;

                al0 = AE_LA128_PP(pY0);
                al1 = AE_LA128_PP(pY1);
                al2 = AE_LA128_PP(pY2);

                AE_L16X4X2_IP(W00, W01, pW, sizeof(ae_int16x8));
                AE_L16X4X2_IP(W10, W11, pW, sizeof(ae_int16x8));
                AE_L16X4X2_IP(W20, W21, pW, sizeof(ae_int16x8));
                for (j = 0; j < (Q & ~7); j += 8)
                {
                    AE_L64X2_IP(S0, S1, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S2, S3, pT_read, sizeof(ae_int64x2));
                    AE_L64X2_IP(S4, S5, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S6, S7, pT_read, sizeof(ae_int64x2));

                    AE_LA16X4X2_IP(Y01, Y02, al0, pY0);
                    AE_LA16X4X2_IP(Y11, Y12, al1, pY1);
                    AE_LA16X4X2_IP(Y21, Y22, al2, pY2);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);


                    AE_MULAFQ16X2_FIR_3(S0, S1, Y10, Y11, W10);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y10, Y11, W10);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y11, Y12, W10);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y11, Y12, W10);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y11, Y11, W11);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y11, Y11, W11);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y12, Y12, W11);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y12, Y12, W11);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y20, Y21, W20);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y20, Y21, W20);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y21, Y22, W20);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y21, Y22, W20);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y21, Y21, W21);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y21, Y21, W21);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y22, Y22, W21);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y22, Y22, W21);

                    AE_MOVD16X8(Y00, Y10, Y02, Y12); Y20 = Y22;

                    T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                    T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 32 - rsh), AE_TRUNCA32X2F64S(S6, S7, 32 - rsh));

                    AE_SA16X4X2_IP(T0, T1, alZ, pZ0);

                }
                AE_SA128POS_FP(alZ, pZ0);
                // tail
                if (j < Q)
                {
                    AE_L64X2_IP(S0, S1, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S2, S3, pT_read, sizeof(ae_int64x2));
                    AE_L64X2_IP(S4, S5, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S6, S7, pT_read, sizeof(ae_int64x2));

                    AE_LA16X4X2_IP(Y01, Y02, al0, pY0);
                    AE_LA16X4X2_IP(Y11, Y12, al1, pY1);
                    AE_LA16X4X2_IP(Y21, Y22, al2, pY2);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);


                    AE_MULAFQ16X2_FIR_3(S0, S1, Y10, Y11, W10);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y10, Y11, W10);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y11, Y12, W10);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y11, Y12, W10);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y11, Y11, W11);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y11, Y11, W11);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y12, Y12, W11);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y12, Y12, W11);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y20, Y21, W20);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y20, Y21, W20);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y21, Y22, W20);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y21, Y22, W20);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y21, Y21, W21);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y21, Y21, W21);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y22, Y22, W21);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y22, Y22, W21);

                    T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                    T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 32 - rsh), AE_TRUNCA32X2F64S(S6, S7, 32 - rsh));

                    for (; j < Q; ++j)
                    {
                        AE_DSEL16X4(T0, T1, T0, T1, AE_MOVINT16X4_FROMINT64(0x0602050104000307));
                        AE_DSEL16X4(Y00, Y01, Y00, Y01, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                        AE_DSEL16X4(Y01, Y02, Y01, Y02, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                        AE_DSEL16X4(Y10, Y11, Y10, Y11, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                        AE_DSEL16X4(Y11, Y12, Y11, Y12, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                        AE_DSEL16X4(Y20, Y21, Y20, Y21, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                        AE_DSEL16X4(Y21, Y22, Y21, Y22, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                        AE_S16_0_IP(T1, castxcc(ae_int16, pZ0), sizeof(ae_int16));
                    }
                }
                //last N-1 elements
                {
                    AE_L64X2_IP(S0, S1, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S2, S3, pT_read, sizeof(ae_int64x2));
                    AE_L64X2_IP(S4, S5, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S6, S7, pT_read, sizeof(ae_int64x2));

                    Y01 = Y02 = 0;
                    Y11 = Y12 = 0;
                    Y21 = Y22 = 0;

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);


                    AE_MULAFQ16X2_FIR_3(S0, S1, Y10, Y11, W10);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y10, Y11, W10);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y11, Y12, W10);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y11, Y12, W10);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y11, Y11, W11);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y11, Y11, W11);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y12, Y12, W11);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y12, Y12, W11);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y20, Y21, W20);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y20, Y21, W20);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y21, Y22, W20);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y21, Y22, W20);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y21, Y21, W21);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y21, Y21, W21);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y22, Y22, W21);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y22, Y22, W21);

                    T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                    T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 32 - rsh), AE_TRUNCA32X2F64S(S6, S7, 32 - rsh));

                    for (; j < Q + N - 1; ++j)
                    {
                        AE_DSEL16X4(T0, T1, T0, T1, AE_MOVINT16X4_FROMINT64(0x0602050104000307));
                        AE_S16_0_IP(T1, castxcc(ae_int16, pZ0), sizeof(ae_int16));
                    }
                }
            }
        }
        else if (mmax - mmin != 1)
        {
            pY = (ae_int16x8*)((ae_int16*)pY0 + Q);
            // first iteration
            {
                pT = (ae_int64x2*)pScr + 6;

                Y00 = 0;
                al0 = AE_LA128_PP(pY0);

                AE_L16X4X2_IP(W00, W01, pW, sizeof(ae_int16x8));
                for (j = 0; j < (Q & ~7); j += 8)
                {
                    AE_LA16X4X2_IP(Y01, Y02, al0, pY0);

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
                    AE_LA16X4X2_IP(Y01, Y02, al0, pY0);

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
            pY = (ae_int16x8*)((ae_int16*)pY0 + Q);
            __Pragma("loop_count min=0, max=2");
            for (m = mmin + 1; m < mmax - 1; ++m)
            {
                pT_read = pT = (ae_int64x2*)pScr + 6;

                Y00 = 0;
                al0 = AE_LA128_PP(pY0);

                AE_L16X4X2_IP(W00, W01, pW, sizeof(ae_int16x8));
                for (j = 0; j < (Q & ~7); j += 8)
                {
                    AE_L64X2_IP(S0, S1, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S2, S3, pT_read, sizeof(ae_int64x2));
                    AE_L64X2_IP(S4, S5, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S6, S7, pT_read, sizeof(ae_int64x2));

                    AE_LA16X4X2_IP(Y01, Y02, al0, pY0);

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

                    AE_LA16X4X2_IP(Y01, Y02, al0, pY0);

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
                pY = (ae_int16x8*)((ae_int16*)pY0 + Q);
            }
            //last iteration
            {
                pT_read = (ae_int64x2*)pScr + 6;

                Y00 = 0;

                al0 = AE_LA128_PP(pY0);

                AE_L16X4X2_IP(W00, W01, pW, sizeof(ae_int16x8));
                for (j = 0; j < (Q & ~7); j += 8)
                {
                    AE_L64X2_IP(S0, S1, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S2, S3, pT_read, sizeof(ae_int64x2));
                    AE_L64X2_IP(S4, S5, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S6, S7, pT_read, sizeof(ae_int64x2));

                    AE_LA16X4X2_IP(Y01, Y02, al0, pY0);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                    Y00 = Y02;

                    T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                    T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 32 - rsh), AE_TRUNCA32X2F64S(S6, S7, 32 - rsh));

                    AE_SA16X4X2_IP(T0, T1, alZ, pZ0);

                }
                AE_SA128POS_FP(alZ, pZ0);
                // tail
                if (j < Q)
                {
                    AE_L64X2_IP(S0, S1, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S2, S3, pT_read, sizeof(ae_int64x2));
                    AE_L64X2_IP(S4, S5, pT_read, sizeof(ae_int64x2)); AE_L64X2_IP(S6, S7, pT_read, sizeof(ae_int64x2));

                    AE_LA16X4X2_IP(Y01, Y02, al0, pY0);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                    AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                    AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                    AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                    T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                    T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 32 - rsh), AE_TRUNCA32X2F64S(S6, S7, 32 - rsh));

                    for (; j < Q; ++j)
                    {
                        AE_DSEL16X4(T0, T1, T0, T1, AE_MOVINT16X4_FROMINT64(0x0602050104000307));
                        AE_DSEL16X4(Y00, Y01, Y00, Y01, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                        AE_DSEL16X4(Y01, Y02, Y01, Y02, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                        AE_S16_0_IP(T1, castxcc(ae_int16, pZ0), sizeof(ae_int16));
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

                    T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                    T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 32 - rsh), AE_TRUNCA32X2F64S(S6, S7, 32 - rsh));

                    for (; j < Q + N - 1; ++j)
                    {
                        AE_DSEL16X4(T0, T1, T0, T1, AE_MOVINT16X4_FROMINT64(0x0602050104000307));
                        AE_S16_0_IP(T1, castxcc(ae_int16, pZ0), sizeof(ae_int16));
                    }
                }
            }
        }
        else
        {
            Y00 = 0;
            al0 = AE_LA128_PP(pY0);

            AE_L16X4X2_IP(W00, W01, pW, sizeof(ae_int16x8));
            for (j = 0; j < (Q & ~7); j += 8)
            {

                AE_LA16X4X2_IP(Y01, Y02, al0, pY0);

                AE_MULFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                AE_MULFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                AE_MULFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                AE_MULFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                Y00 = Y02;

                T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 32 - rsh), AE_TRUNCA32X2F64S(S6, S7, 32 - rsh));

                AE_SA16X4X2_IP(T0, T1, alZ, pZ0);

            }
            AE_SA128POS_FP(alZ, pZ0);
            // tail
            if (j < Q)
            {
                AE_LA16X4X2_IP(Y01, Y02, al0, pY0);

                AE_MULFQ16X2_FIR_3(S0, S1, Y00, Y01, W00);
                AE_MULFQ16X2_FIR_1(S2, S3, Y00, Y01, W00);
                AE_MULFQ16X2_FIR_3(S4, S5, Y01, Y02, W00);
                AE_MULFQ16X2_FIR_1(S6, S7, Y01, Y02, W00);

                AE_MULAFQ16X2_FIR_3(S0, S1, Y01, Y01, W01);
                AE_MULAFQ16X2_FIR_1(S2, S3, Y01, Y01, W01);
                AE_MULAFQ16X2_FIR_3(S4, S5, Y02, Y02, W01);
                AE_MULAFQ16X2_FIR_1(S6, S7, Y02, Y02, W01);

                T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 32 - rsh), AE_TRUNCA32X2F64S(S6, S7, 32 - rsh));

                for (; j < Q; ++j)
                {
                    AE_DSEL16X4(T0, T1, T0, T1, AE_MOVINT16X4_FROMINT64(0x0602050104000307));
                    AE_DSEL16X4(Y00, Y01, Y00, Y01, AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    AE_DSEL16X4(Y01, Y02, Y01, Y02, AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    AE_S16_0_IP(T1, castxcc(ae_int16, pZ0), sizeof(ae_int16));
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

                T0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
                T1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 32 - rsh), AE_TRUNCA32X2F64S(S6, S7, 32 - rsh));

                for (; j < Q + N - 1; ++j)
                {
                    AE_DSEL16X4(T0, T1, T0, T1, AE_MOVINT16X4_FROMINT64(0x0602050104000307));
                    AE_S16_0_IP(T1, castxcc(ae_int16, pZ0), sizeof(ae_int16));
                }
            }
        }
    }
#undef M
#undef N
}

// scratch allocatation functions. return required scratch size in bytes
size_t conv2d_gen_5x5_16x16_getScratchSize(int P, int Q)
{
    return 6 * 8 * sizeof(int16_t) + (Q + 16) * sizeof(int64_t);
} // MxN=5x5
