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

void conv2d_3x3_16x16(void * pScr, int16_t * z, const int16_t * x, const int16_t * y, int rsh, int P, int Q)
{
#define M 3
#define N 3
    const ae_int16x8 * restrict pY0;
    const ae_int16x8 * restrict pY1;
    const ae_int16x8 * restrict pY2;
    const ae_int16x8 * restrict pY3;
          ae_int16x8 * restrict pZ0;
          ae_int16x8 * restrict pZ1;
    const ae_int64x2 * restrict pW;
    int i,j;
    ae_valignx2 aZ0, aZ1;
    ae_int16x4 w0, w1, w2;
    ae_int16x4 y00, y01, y02, y10, y11, y12;
    ae_int16x4 y20, y21, y22, y30, y31, y32;
    ae_int64 S0, S1, S2, S3;
    ae_int16x4 r0, r1;
    ae_int32x2 R;

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
        ae_int64 tmp0, tmp1;

        pW = (const ae_int64x2 *)x;

        AE_L64X2_IP(tmp0, tmp1, pW, 8 * sizeof(int16_t));
        w2_ = AE_MOVINT16X4_FROMINT64(tmp0);
        w1_ = AE_MOVINT16X4_FROMINT64(tmp1);
        w0_ = AE_L16_I((ae_int16 *)pW, 0);

        w2 = AE_SEL16I(w2_, AE_ZERO16(), 5);
        w2 = AE_SEL16I(AE_ZERO16(), w2, 4);

        w1 = AE_SEL16I(w1_, w2_, 0);
        w1 = AE_SEL16I(AE_ZERO16(), w1, 4);

        w0 = AE_SEL16I(w0_, w1_, 4);
        w0 = AE_SEL16I(AE_ZERO16(), w0, 4);
    }

    pZ0 = (ae_int16x8 *)z;
    pZ1 = (ae_int16x8 *)(z + Q + N - 1);
    aZ0 = AE_ZALIGN128();
    aZ1 = AE_ZALIGN128();

    pY0 = (ae_int16x8 *)y;
    pY1 = (ae_int16x8 *)((int16_t *)pY0 + Q);
    pY2 = (ae_int16x8 *)((int16_t *)pY1 + Q);
    pY3 = (ae_int16x8 *)((int16_t *)pY2 + Q);

    {
        y00 = AE_ZERO16();
        y10 = AE_ZERO16();

        /* Process all but last N-1 samples of the 0 and 1 rows */
        __Pragma("loop_count min=1");
        for (j = 0; j < (Q >> 3); j++)
        {
            AE_L16X4X2_IP(y01, y02, pY0, 8 * sizeof(int16_t));
            AE_L16X4X2_IP(y11, y12, pY1, 8 * sizeof(int16_t));

            AE_MULFQ16X2_FIR_2(S0, S1, y00, y01, w2);
            AE_MULFQ16X2_FIR_0(S2, S3, y00, y01, w2);
            r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));

            AE_MULFQ16X2_FIR_2(S0, S1, y01, y02, w2);
            AE_MULFQ16X2_FIR_0(S2, S3, y01, y02, w2);
            r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
            AE_SA16X4X2_IP(r0, r1, aZ0, pZ0);


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
            AE_SA16X4X2_IP(r0, r1, aZ1, pZ1);

            //y00 = y02;
            //y10 = y12;
            y00 = AE_L16X4_I((ae_int16x4 *)pY0, -4 * (int)sizeof(int16_t));
            y10 = AE_L16X4_I((ae_int16x4 *)pY1, -4 * (int)sizeof(int16_t));
        }
        AE_SA128POS_FP(aZ0, pZ0);
        AE_SA128POS_FP(aZ1, pZ1);

        /* Last N-1 samples */
        {
            AE_MULFQ16X2_FIR_2(S0, S1, y00, AE_ZERO16(), w2);

            AE_MULFQ16X2_FIR_2(S2, S3, y00, AE_ZERO16(), w1);
            AE_MULAFQ16X2_FIR_2(S2, S3, y10, AE_ZERO16(), w2);

            R = AE_MOVINT32X2_FROMINT16X4(AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S3, S2, 32 - rsh), AE_TRUNCA32X2F64S(S1, S0, 32 - rsh)));

            AE_S32_L_XP(R, castxcc(ae_int32, pZ0), (Q + (N - 1) + 2) * sizeof(int16_t));
            AE_S32_H_XP(R, castxcc(ae_int32, pZ1), (Q + (N - 1) + 2) * sizeof(int16_t));
        }

        pY0 = (ae_int16x8 *)XT_ADDX2(-Q, (uintptr_t)pY0);
        pY1 = (ae_int16x8 *)XT_ADDX2(-Q, (uintptr_t)pY1);
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
            AE_L16X4X2_IP(y01, y02, pY0, 8 * sizeof(int16_t));
            AE_L16X4X2_IP(y11, y12, pY1, 8 * sizeof(int16_t));
            AE_L16X4X2_IP(y21, y22, pY2, 8 * sizeof(int16_t));
            AE_L16X4X2_IP(y31, y32, pY3, 8 * sizeof(int16_t));

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
            AE_SA16X4X2_IP(r0, r1, aZ0, pZ0);


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
            AE_SA16X4X2_IP(r0, r1, aZ1, pZ1);

            y00 = y02;
            y10 = y12;
            y20 = y22;
            y30 = y32;
        }
        AE_SA128POS_FP(aZ0, pZ0);
        AE_SA128POS_FP(aZ1, pZ1);

        /* Last N-1 samples */
        {
            AE_MULFQ16X2_FIR_2(S0, S1, y00, AE_ZERO16(), w0);
            AE_MULAFQ16X2_FIR_2(S0, S1, y10, AE_ZERO16(), w1);
            AE_MULAFQ16X2_FIR_2(S0, S1, y20, AE_ZERO16(), w2);

            AE_MULFQ16X2_FIR_2(S2, S3, y10, AE_ZERO16(), w0);
            AE_MULAFQ16X2_FIR_2(S2, S3, y20, AE_ZERO16(), w1);
            AE_MULAFQ16X2_FIR_2(S2, S3, y30, AE_ZERO16(), w2);

            R = AE_MOVINT32X2_FROMINT16X4(AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S3, S2, 32 - rsh), AE_TRUNCA32X2F64S(S1, S0, 32 - rsh)));

            AE_S32_L_XP(R, castxcc(ae_int32, pZ0), (Q + (N - 1) + 2) * sizeof(int16_t));
            AE_S32_H_XP(R, castxcc(ae_int32, pZ1), (Q + (N - 1) + 2) * sizeof(int16_t));
        }

        pY0 = (ae_int16x8 *)XT_ADDX2(Q, (uintptr_t)pY0);
        pY1 = (ae_int16x8 *)XT_ADDX2(Q, (uintptr_t)pY1);
        pY2 = (ae_int16x8 *)XT_ADDX2(Q, (uintptr_t)pY2);
        pY3 = (ae_int16x8 *)XT_ADDX2(Q, (uintptr_t)pY3);
    }

    {
        y10 = AE_ZERO16();
        y20 = AE_ZERO16();

        /* Process all but last N-1 samples of the (P+M-3) and (P+M-2) rows */
        __Pragma("loop_count min=1");
        for (j = 0; j < (Q >> 3); j++)
        {
            AE_L16X4X2_IP(y11, y12, pY0, 8 * sizeof(int16_t));
            AE_L16X4X2_IP(y21, y22, pY1, 8 * sizeof(int16_t));

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
            AE_SA16X4X2_IP(r0, r1, aZ0, pZ0);


            AE_MULFQ16X2_FIR_2(S0, S1, y20, y21, w0);
            AE_MULFQ16X2_FIR_0(S2, S3, y20, y21, w0);
            r0 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));

            AE_MULFQ16X2_FIR_2(S0, S1, y21, y22, w0);
            AE_MULFQ16X2_FIR_0(S2, S3, y21, y22, w0);
            r1 = AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S0, S1, 32 - rsh), AE_TRUNCA32X2F64S(S2, S3, 32 - rsh));
            AE_SA16X4X2_IP(r0, r1, aZ1, pZ1);

            y10 = y12;
            y20 = y22;
        }
        AE_SA128POS_FP(aZ0, pZ0);
        AE_SA128POS_FP(aZ1, pZ1);

        /* Last N-1 samples */
        {
            AE_MULFQ16X2_FIR_2(S0, S1, y10, AE_ZERO16(), w0);
            AE_MULAFQ16X2_FIR_2(S0, S1, y20, AE_ZERO16(), w1);

            AE_MULFQ16X2_FIR_2(S2, S3, y20, AE_ZERO16(), w0);

            R = AE_MOVINT32X2_FROMINT16X4(AE_ROUND16X4F32SSYM(AE_TRUNCA32X2F64S(S3, S2, 32 - rsh), AE_TRUNCA32X2F64S(S1, S0, 32 - rsh)));

            AE_S32_L_XP(R, castxcc(ae_int32, pZ0), (Q + (N - 1) + 2) * sizeof(int16_t));
            AE_S32_H_XP(R, castxcc(ae_int32, pZ1), (Q + (N - 1) + 2) * sizeof(int16_t));
        }
    }
#undef M
#undef N
}

// scratch allocatation functions. return required scratch size in bytes
size_t conv2d_3x3_16x16_getScratchSize(int P, int Q)
{
    return 0;
}
