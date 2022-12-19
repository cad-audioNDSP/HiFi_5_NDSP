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
    y[M+3]     - input padded with at least 3 zeroes. should be aligned on 16-byte boundary
    output:
    r[N]       - correlations
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "raw_corr32x32.h"

void raw_corr32x32( int32_t * restrict r,
              const int32_t * restrict x,
              const int32_t * restrict y,
              int N, int M )
{
    const ae_int32x4 *          pX;
    const ae_int32x4 *          S0;
    const ae_int32x4 *          S1;
    const ae_int32x4 *          S2;
    const ae_int32x2 *          pY;
          ae_int32x2 * restrict pR;

    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_int32x2 x0, x1, x2, x3, x4, x5;
    ae_int32x2 y0, y1;
    ae_valign ar;

    int n, m;

    NASSERT(r);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(r, 4);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 8);
    NASSERT(N > 0 && M > 0 && N >= M - 1);

    M = (M + 3)&~3;

    pX = (const ae_int32x4 *)x;
    pR = (      ae_int32x2 *)r;
    ar = AE_ZALIGN64();

    for (n = 0; n < (N >> 3); n++)
    {
        pY = (const ae_int32x2 *)y;
        S0 = pX;
        pX = (const ae_int32x4 *)XT_ADDX4(4, (uintptr_t)pX);
        S1 = pX;
        pX = (const ae_int32x4 *)XT_ADDX4(4, (uintptr_t)pX);
        S2 = pX;

        q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            AE_L32X2_IP(y1, pY, 2 * sizeof(int32_t));
            AE_L32X2X2_IP(x0, x1, S0, 4 * sizeof(int32_t));
            AE_L32X2X2_IP(x2, x3, S1, 4 * sizeof(int32_t));
            AE_L32X2X2_IP(x4, x5, S2, 4 * sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
            AE_MULAFD32X2RA_FIR_H(q2, q3, x1, x2, y0);
            AE_MULAFD32X2RA_FIR_H(q4, q5, x2, x3, y0);
            AE_MULAFD32X2RA_FIR_H(q6, q7, x3, x4, y0);
            AE_MULAFD32X2RA_FIR_H(q0, q1, x1, x2, y1);
            AE_MULAFD32X2RA_FIR_H(q2, q3, x2, x3, y1);
            AE_MULAFD32X2RA_FIR_H(q4, q5, x3, x4, y1);
            AE_MULAFD32X2RA_FIR_H(q6, q7, x4, x5, y1);
        }

        x0 = AE_ROUND32X2F48SASYM(q0, q1);
        x1 = AE_ROUND32X2F48SASYM(q2, q3);
        x2 = AE_ROUND32X2F48SASYM(q4, q5);
        x3 = AE_ROUND32X2F48SASYM(q6, q7);
        AE_SA32X2_IP(x0, ar, pR);
        AE_SA32X2_IP(x1, ar, pR);
        AE_SA32X2_IP(x2, ar, pR);
        AE_SA32X2_IP(x3, ar, pR);
    }
    if (N & 4)
    {
        pY = (const ae_int32x2 *)y;
        S0 = pX;
        pX = (const ae_int32x4 *)XT_ADDX4(4, (uintptr_t)pX);
        S1 = pX;
        q0 = q1 = q2 = q3 = AE_ZERO64();

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            AE_L32X2_IP(y1, pY, 2 * sizeof(int32_t));
            AE_L32X2X2_IP(x0, x1, S0, 4 * sizeof(int32_t));
            AE_L32X2X2_IP(x2, x3, S1, 4 * sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
            AE_MULAFD32X2RA_FIR_H(q2, q3, x1, x2, y0);
            AE_MULAFD32X2RA_FIR_H(q0, q1, x1, x2, y1);
            AE_MULAFD32X2RA_FIR_H(q2, q3, x2, x3, y1);
        }

        x0 = AE_ROUND32X2F48SASYM(q0, q1);
        x1 = AE_ROUND32X2F48SASYM(q2, q3);
        AE_SA32X2_IP(x0, ar, pR);
        AE_SA32X2_IP(x1, ar, pR);
    }
    AE_SA64POS_FP(ar, pR);

    N &= 3;
    if (N)
    {
        pY = (const ae_int32x2 *)y;
        S0 = pX;
        pX = (const ae_int32x4 *)XT_ADDX4(4, (uintptr_t)pX);
        S1 = pX;
        q0 = q1 = q2 = q3 = AE_ZERO64();

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            AE_L32X2_IP(y1, pY, 2 * sizeof(int32_t));
            AE_L32X2X2_IP(x0, x1, S0, 4 * sizeof(int32_t));
            AE_L32X2X2_IP(x2, x3, S1, 4 * sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
            AE_MULAFD32X2RA_FIR_H(q2, q3, x1, x2, y0);
            AE_MULAFD32X2RA_FIR_H(q0, q1, x1, x2, y1);
            AE_MULAFD32X2RA_FIR_H(q2, q3, x2, x3, y1);
        }

        AE_S32RA64S_IP(q0, castxcc(ae_int32, pR), sizeof(int32_t));
        if (N > 1) AE_S32RA64S_IP(q1, castxcc(ae_int32, pR), sizeof(int32_t));
        if (N > 2) AE_S32RA64S_IP(q2, castxcc(ae_int32, pR), sizeof(int32_t));
    }
}
