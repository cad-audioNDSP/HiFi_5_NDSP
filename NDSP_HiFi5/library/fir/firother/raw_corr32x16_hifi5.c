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
    Helper for circular correlation/convolution with originally non-aligned inputs
    This function takes aligned inputs (allocated from the scratch) with circularly
    duplicated inputs

    C code optimized for HiFi5
    IntegrIT, 2006-2019

    raw linear correlation
    input:
    x[N+M-1 ]  - input . should be aligned on 16-byte boundary
    y[1+M+7]   - input padded with 1 zero from the left side and 7 zeroes from
                 the right side. Should be aligned on 16-byte boundary
    output:
    r[N]       - correlations
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "raw_corr32x16.h"

void raw_corr32x16( int32_t * restrict r,
              const int32_t * restrict x,
              const int16_t * restrict y,
              int N, int M )
{
    const ae_int32x4 *          pX;
    const ae_int32x4 *          S0;
    const ae_int32x4 *          S1;
    const ae_int32x4 *          S2;
    const ae_int32x4 *          S3;
    const ae_int32x4 *          S4;
    const ae_int16x4 *          pY;
          ae_int32x2 * restrict pR;

    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
    ae_int32x2 x0, x1, x2, x3, x4, x5;
    ae_int32x2 x6, x7, x8, x9;
    ae_int16x4 y0;
    ae_valign ar;

    int n, m;

    NASSERT(r);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(r, 4);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 8);
    NASSERT(N > 0 && M > 0 && N >= M - 1);
    WUR_AE_CBEGIN0((uintptr_t)x);
    WUR_AE_CEND0  ((((uintptr_t)(x+N+M-1))+15)&~15);

    M = (M + 3)&~3;

    pX = (const ae_int32x4 *)x;
    pR = (      ae_int32x2 *)r;
    ar = AE_ZALIGN64();

    for (n = 0; n < (N >> 4); n++)
    {
        pY = (const ae_int16x4 *)y;
        AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
        AE_L32X2X2_IP(x0, x1, pX, 4 * sizeof(int32_t)); S0 = pX;
        AE_L32X2X2_IP(x2, x3, pX, 4 * sizeof(int32_t)); S1 = pX;
        AE_L32X2X2_IP(x4, x5, pX, 4 * sizeof(int32_t)); S2 = pX;
        AE_L32X2X2_IP(x6, x7, pX, 4 * sizeof(int32_t)); S3 = pX; S4 = pX;
        AE_L32X2X2_IP(x8, x9, S4, 4 * sizeof(int32_t));
        AE_MUL2Q32X16_FIR_H(q0, q1, x0, x1, x2, y0);
        AE_MUL2Q32X16_FIR_H(q2, q3, x1, x2, x3, y0);
        AE_MUL2Q32X16_FIR_H(q4, q5, x2, x3, x4, y0);
        AE_MUL2Q32X16_FIR_H(q6, q7, x3, x4, x5, y0);
        AE_MUL2Q32X16_FIR_H(q8, q9, x4, x5, x6, y0);
        AE_MUL2Q32X16_FIR_H(qa, qb, x5, x6, x7, y0);
        AE_MUL2Q32X16_FIR_H(qc, qd, x6, x7, x8, y0);
        AE_MUL2Q32X16_FIR_H(qe, qf, x7, x8, x9, y0);
        for (m = 0; m < (M >> 2) - 1; m++)
        {
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
            AE_L32X2X2_IP(x0, x1, S0, 4 * sizeof(int32_t));
            AE_L32X2X2_IP(x2, x3, S1, 4 * sizeof(int32_t));
            AE_L32X2X2_IP(x4, x5, S2, 4 * sizeof(int32_t));
            AE_L32X2X2_IP(x6, x7, S3, 4 * sizeof(int32_t));
            AE_L32X2X2_IP(x8, x9, S4, 4 * sizeof(int32_t));
            AE_MULA2Q32X16_FIR_H(q0, q1, x0, x1, x2, y0);
            AE_MULA2Q32X16_FIR_H(q2, q3, x1, x2, x3, y0);
            AE_MULA2Q32X16_FIR_H(q4, q5, x2, x3, x4, y0);
            AE_MULA2Q32X16_FIR_H(q6, q7, x3, x4, x5, y0);
            AE_MULA2Q32X16_FIR_H(q8, q9, x4, x5, x6, y0);
            AE_MULA2Q32X16_FIR_H(qa, qb, x5, x6, x7, y0);
            AE_MULA2Q32X16_FIR_H(qc, qd, x6, x7, x8, y0);
            AE_MULA2Q32X16_FIR_H(qe, qf, x7, x8, x9, y0);
        }

        q0 = AE_SLAI64(q0, 1); q1 = AE_SLAI64(q1, 1);
        AE_SA32X2_IP(AE_ROUND32X2F48SASYM(q0, q1), ar, pR);
        AE_PKSR32(x0, q2, 1); AE_PKSR32(x0, q3, 1);
        AE_SA32X2_IP(x0, ar, pR);
        q4 = AE_SLAI64(q4, 1); q5 = AE_SLAI64(q5, 1);
        AE_SA32X2_IP(AE_ROUND32X2F48SASYM(q4, q5), ar, pR);
        AE_PKSR32(x0, q6, 1); AE_PKSR32(x0, q7, 1);
        AE_SA32X2_IP(x0, ar, pR);
        q8 = AE_SLAI64(q8, 1); q9 = AE_SLAI64(q9, 1);
        AE_SA32X2_IP(AE_ROUND32X2F48SASYM(q8, q9), ar, pR);
        AE_PKSR32(x0, qa, 1); AE_PKSR32(x0, qb, 1);
        AE_SA32X2_IP(x0, ar, pR);
        qc = AE_SLAI64(qc, 1); qd = AE_SLAI64(qd, 1);
        AE_SA32X2_IP(AE_ROUND32X2F48SASYM(qc, qd), ar, pR);
        AE_PKSR32(x0, qe, 1); AE_PKSR32(x0, qf, 1);
        AE_SA32X2_IP(x0, ar, pR);
    }
    if (N & 8)
    {
        pY = (const ae_int16x4 *)y;
        S0 = pX;
        pX = (const ae_int32x4 *)XT_ADDX4(4, (uintptr_t)pX);
        S1 = pX;
        pX = (const ae_int32x4 *)XT_ADDX4(4, (uintptr_t)pX);
        S2 = pX;

        q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
            AE_L32X2X2_IP(x0, x1, S0, 4 * sizeof(int32_t));
            AE_L32X2X2_IP(x2, x3, S1, 4 * sizeof(int32_t));
            AE_L32X2X2_IP(x4, x5, S2, 4 * sizeof(int32_t));
            AE_MULA2Q32X16_FIR_H(q0, q1, x0, x1, x2, y0);
            AE_MULA2Q32X16_FIR_H(q2, q3, x1, x2, x3, y0);
            AE_MULA2Q32X16_FIR_H(q4, q5, x2, x3, x4, y0);
            AE_MULA2Q32X16_FIR_H(q6, q7, x3, x4, x5, y0);
        }

        q0 = AE_SLAI64(q0, 1); q1 = AE_SLAI64(q1, 1);
        AE_SA32X2_IP(AE_ROUND32X2F48SASYM(q0, q1), ar, pR);
        AE_PKSR32(x0, q2, 1); AE_PKSR32(x0, q3, 1);
        AE_SA32X2_IP(x0, ar, pR);
        q4 = AE_SLAI64(q4, 1); q5 = AE_SLAI64(q5, 1);
        AE_SA32X2_IP(AE_ROUND32X2F48SASYM(q4, q5), ar, pR);
        AE_PKSR32(x0, q6, 1); AE_PKSR32(x0, q7, 1);
        AE_SA32X2_IP(x0, ar, pR);
    }
    if (N & 4)
    {
        pY = (const ae_int16x4 *)y;
        AE_L32X2X2_IP(x0, x1, pX, 4 * sizeof(int32_t));
        S0 = pX;
        q0 = q1 = q2 = q3 = AE_ZERO64();

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
            AE_L32X2X2_IP(x2, x3, S0, 4 * sizeof(int32_t));
            AE_MULA2Q32X16_FIR_H(q0, q1, x0, x1, x2, y0);
            AE_MULA2Q32X16_FIR_H(q2, q3, x1, x2, x3, y0);
            x0 = x2; x1 = x3;
        }

        q0 = AE_SLAI64(q0, 1); q1 = AE_SLAI64(q1, 1);
        AE_SA32X2_IP(AE_ROUND32X2F48SASYM(q0, q1), ar, pR);
        q2 = AE_SLAI64(q2, 1); q3 = AE_SLAI64(q3, 1);
        AE_SA32X2_IP(AE_ROUND32X2F48SASYM(q2, q3), ar, pR);
    }
    AE_SA64POS_FP(ar, pR);

    N &= 3;
    if (N)
    {
        pY = (const ae_int16x4 *)y;
        AE_L32X2X2_XC(x0, x1, pX, 4 * sizeof(int32_t));
        q0 = q1 = q2 = q3 = AE_ZERO64();

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
            AE_L32X2X2_XC(x2, x3, pX, 4 * sizeof(int32_t));
            AE_MULA2Q32X16_FIR_H(q0, q1, x0, x1, x2, y0);
            AE_MULA2Q32X16_FIR_H(q2, q3, x1, x2, x3, y0);
            x0 = x2; x1 = x3;
        }

        q0 = AE_SLAI64(q0, 1);
        q1 = AE_SLAI64(q1, 1);
        q2 = AE_SLAI64(q2, 1);
        AE_S32RA64S_IP(q0, castxcc(ae_int32, pR), sizeof(int32_t));
        if (N > 1) AE_S32RA64S_IP(q1, castxcc(ae_int32, pR), sizeof(int32_t));
        if (N > 2) AE_S32RA64S_IP(q2, castxcc(ae_int32, pR), sizeof(int32_t));
    }
}
