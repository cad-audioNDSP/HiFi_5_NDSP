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

void conv2d_11x7_16x16(void * pScr, int16_t * z, const int16_t * x, const int16_t * y, int rsh, int P, int Q)
{
#define M 11
#define N 7
    int16_t * w;
    const ae_int16x8 * restrict pX;
    const ae_int16x8 * restrict pY0;
    const ae_int16x8 * restrict pY1;
          ae_int16x4 * restrict pZ;
          ae_int64x2 * restrict pWwr;
    const ae_int16x4 * restrict pW;
    int i, j, n, n0, n1, m1;
    ae_f64 S0, S1, S2, S3, S4, S5, S6, S7;
    ae_int16x4 y0, y1, y2, y3;
    ae_int16x4 w0, w1;
    ae_int16x4 z0;
    ae_valign aZ;

    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT_ALIGN(z, 16);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0 && P % 8 == 0);
    NASSERT(Q >= 0 && Q % 8 == 0);
    if (P <= 0 || Q <= 0) return;

    w = (int16_t *)pScr;
    /* Store coefficients in the next order:
     *   x[6] x[5] x[4] x[3]  x[2] x[1] x[0]  0
     *    0   x[6] x[5] x[4]  x[3] x[2] x[1] x[0]
     * Start from the last row
     */
    pX = (ae_int16x8 *)(x + 10 * 7);
    pWwr = (ae_int64x2 *)w;
#if 0 //switch on if ferret warning will occur
    {
        ae_valign aX;
        ae_int16x4 x0, x1, x2, x3;

        aX = AE_LA64_PP(pX);
        AE_LA16X4_IP(x0, aX, castxcc(ae_int16x4, pX));
        pX = (ae_int16x8 *)XT_ADDX2(-1, (uintptr_t)pX);
        aX = AE_LA64_PP(pX);
        AE_LA16X4_IP(x3, aX, castxcc(ae_int16x4, pX));
        pX = (ae_int16x8 *)XT_ADDX2(-14, (uintptr_t)pX);
        x1 = AE_SEL16_6543(x3, AE_ZERO16());
        x2 = AE_SEL16_4321(AE_ZERO16(), x0);

        AE_S64X2_X(AE_MOVINT64_FROMINT16X4(x1), AE_MOVINT64_FROMINT16X4(x0), pWwr, 8 * M*sizeof(int16_t));
        AE_S64X2_IP(AE_MOVINT64_FROMINT16X4(x3), AE_MOVINT64_FROMINT16X4(x2), pWwr, 8 * sizeof(int16_t));
    }
    for (i = 0; i < M - 1; i++)
#else
    for (i = 0; i < M; i++)
#endif
    {
        ae_valignx2 aX;
        ae_int16x4 x0, x1, x2, x3;

        aX = AE_LA128_PP(pX);
        AE_LA16X4X2_IP(x0, x1, aX, pX);
        pX = (ae_int16x8 *)XT_ADDX2(-15, (uintptr_t)pX);
        AE_MOVF16X4(x1, AE_ZERO16(), AE_MOVBA4(14));
        x2 = AE_SEL16_4321(AE_ZERO16(), x0);
        x3 = AE_SEL16_4321(x0, x1);

        AE_S64X2_X(AE_MOVINT64_FROMINT16X4(x1), AE_MOVINT64_FROMINT16X4(x0), pWwr, 8 * M*sizeof(int16_t));
        AE_S64X2_IP(AE_MOVINT64_FROMINT16X4(x3), AE_MOVINT64_FROMINT16X4(x2), pWwr, 8 * sizeof(int16_t));
    }

    /*
     * Processing of convolution
     */
    pZ = (ae_int16x4 *)z;
    aZ = AE_ZALIGN64();
    __Pragma("loop_count min=1");
    for (i = 0; i < M + P - 1; i++)
    {
        n0 = XT_MAX(i + 1 - M, 0);
        n1 = XT_MIN(i + 1, P);
        m1 = XT_MIN(i + 1, M);
        WAE_CBEGIN0((uintptr_t)(y + n0*Q + 8));
        WAE_CEND0  ((uintptr_t)(y + n1*Q));

        /* First 8 samples of the i-th row */
        {
            pY0 = (ae_int16x8 *)(y + n0*Q);
            pW = (ae_int16x4 *)XT_ADDX2(8 * (2 * M - m1), (uintptr_t)w);
            S0 = S1 = S2 = S3 = S4 = S5 = S6 = S7 = AE_ZERO64();
            __Pragma("loop_count min=1, max=11");
            for (n = 0; n < n1 - n0; n++)
            {
                AE_L16X4X2_XP(y0, y1, pY0, Q * sizeof(int16_t));

                AE_L16X4_IP(w0, pW, 4 * sizeof(int16_t));
                AE_L16X4_IP(w1, pW, 4 * sizeof(int16_t));

                AE_MULAFQ16X2_FIR_2(S0, S1, AE_ZERO16(), y0, w1);
                AE_MULAFQ16X2_FIR_0(S2, S3, AE_ZERO16(), y0, w1);
                AE_MULAFQ16X2_FIR_2(S4, S5, AE_ZERO16(), y0, w0);
                AE_MULAFQ16X2_FIR_2(S4, S5, y0, y1, w1);
                AE_MULAFQ16X2_FIR_0(S6, S7, AE_ZERO16(), y0, w0);
                AE_MULAFQ16X2_FIR_0(S6, S7, y0, y1, w1);
            }
            z0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
            AE_SA16X4_IP(z0, aZ, pZ);
            z0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 32 - rsh), AE_TRUNCA32X2F64S(S6, S7, 32 - rsh));
            AE_SA16X4_IP(z0, aZ, pZ);
        }
        /* Next samples */
        pY0 = (ae_int16x8 *)(y + n0*Q);
        pY1 = pY0 + 1;
        for (j = 0; j < ((Q >> 3) - 1); j++)
        {
            pW = (ae_int16x4 *)XT_ADDX2(8 * (M - m1), (uintptr_t)w);
            S0 = S1 = S2 = S3 = S4 = S5 = S6 = S7 = AE_ZERO64();
            __Pragma("loop_count min=1, max=11");
            for (n = 0; n < n1 - n0; n++)
            {
                AE_L16X4X2_XC(y0, y1, pY0, Q * sizeof(int16_t));
                AE_L16X4X2_XC(y2, y3, pY1, Q * sizeof(int16_t));

                AE_L16X4_IP(w0, pW, 4 * sizeof(int16_t));
                AE_L16X4_IP(w1, pW, 4 * sizeof(int16_t));

                AE_MULAFQ16X2_FIR_1(S0, S1, y0, y1, w0);
                AE_MULAFQ16X2_FIR_1(S0, S1, y1, y2, w1);
                AE_MULAFQ16X2_FIR_3(S2, S3, y1, y2, w0);
                AE_MULAFQ16X2_FIR_3(S2, S3, y2, AE_ZERO16(), w1);
                AE_MULAFQ16X2_FIR_1(S4, S5, y1, y2, w0);
                AE_MULAFQ16X2_FIR_1(S4, S5, y2, y3, w1);
                AE_MULAFQ16X2_FIR_3(S6, S7, y2, y3, w0);
                AE_MULAFQ16X2_FIR_3(S6, S7, y3, AE_ZERO16(), w1);
            }
            z0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
            AE_SA16X4_IP(z0, aZ, pZ);
            z0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S4, S5, 32 - rsh), AE_TRUNCA32X2F64S(S6, S7, 32 - rsh));
            AE_SA16X4_IP(z0, aZ, pZ);
        }
        /* Last N-1 samples of the i-th row */
        {
            pW = (ae_int16x4 *)XT_ADDX2(8 * (M - m1), (uintptr_t)w);
            S0 = S1 = S2 = S3 = S4 = S5 = AE_ZERO64();
            __Pragma("loop_count min=1, max=11");
            for (n = 0; n < n1 - n0; n++)
            {
                AE_L16X4X2_XP(y0, y1, pY0, Q * sizeof(int16_t));

                AE_L16X4_XP(w0, pW, 4 * sizeof(int16_t));
                AE_L16X4_XP(w1, pW, 4 * sizeof(int16_t));

                AE_MULAFQ16X2_FIR_1(S0, S1, y0, y1, w0);
                AE_MULAFQ16X2_FIR_1(S0, S1, y1, AE_ZERO16(), w1);
                AE_MULAFQ16X2_FIR_3(S2, S3, y1, AE_ZERO16(), w0);
                AE_MULAFQ16X2_FIR_1(S4, S5, y1, AE_ZERO16(), w0);
            }
            z0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
            AE_SA16X4_IP(z0, aZ, pZ);
            AE_SA64POS_FP(aZ, pZ);
            z0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S5, S4, 32 - rsh), AE_TRUNCA32X2F64S(S5, S4, 32 - rsh));
            AE_S32_L_IP(AE_MOVINT32X2_FROMINT16X4(z0), castxcc(ae_int32, pZ), 2 * sizeof(int16_t));
        }
    }
#undef M
#undef N
}

// scratch allocatation functions. return required scratch size in bytes
size_t conv2d_11x7_16x16_getScratchSize(int P, int Q)
{
    return 2*11*8*sizeof(int16_t);
} // MxN=11x7
