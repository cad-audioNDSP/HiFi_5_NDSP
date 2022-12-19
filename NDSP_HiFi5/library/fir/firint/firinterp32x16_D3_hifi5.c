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
    Interpolating block real FIR filter, 32x16-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "firinterp32x16_common.h"

/*-----------------------------------------------------------------------------
 * Data processing function of a particular interpolating filter. Stores a
 * block of N input samples to the circular delay line buffer and computes
 * N*D samples of interpolating FIR filter's response.
 * Input:
 *   delayLine - circular delay line buffer start address
 *   delayLen  - Delay line buffer length
 *   wrIx    - next position in the buffer to be filled with an input sample
 *   x[N]    - input samples
 *   h[]     - decimating FIR filter coefficients, array layout varies
 * Output:
 *   y[N*D]  - output samples
 *   retval  - updated index of the oldest sample
 * Notes and restrictions:
 *   1. Most of data processing functions feature a single, hard-coded
 *      interpolation factor, so they expect a determined value for parameter D.
 *   2. All pointers with the exception of y[N] must be aligned on an 16-bytes
 *      boundary.
 *   3. N - must be a multiple of 8.
 *   4. M - must be a multiple of 4.
 -----------------------------------------------------------------------------*/

/* Data processing function for a factor 3 interpolating FIR filter. */
int firinterp32x16_D3_proc( int32_t * restrict y,
                            int32_t * delayLine, int delayLen,
                      const int32_t * restrict x,
                      const int16_t * restrict h,
                      int wrIx, int D, int N, int M )
{
    const ae_int32x4 * restrict pX;
          ae_int32x4 * restrict pDw;
    const ae_int32x4 * restrict pDr;
    const ae_int32x4 * restrict S0;
    const ae_int32x4 * restrict S1;
    const ae_int32x4 * restrict S2;
    const ae_int16x4 * restrict pH;
          ae_int32x4 * restrict pY;

    ae_valignx2 aY;

    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
    ae_f64 q10, q11, q12, q13, q14, q15, q16, q17;
    ae_int32x2 d0, d1, d2, d3, d4, d5;
    ae_int16x4 h0, h1, h2;

    int m, n;

    NASSERT(y && delayLine && x && h);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(h, 16);
    NASSERT_ALIGN(delayLine, 16);
    NASSERT(D == 3);
    NASSERT(M > 0 && M % 4 == 0);
    NASSERT(N > 0 && N % 8 == 0);

    //
    // Setup pointers and circular delay line buffer.
    //
    pX  = (const ae_int32x4 *)x;
    pY  = (      ae_int32x4 *)y;
    pDw = (      ae_int32x4 *)(delayLine + wrIx);
    WUR_AE_CBEGIN0((uintptr_t)(delayLine));
    WUR_AE_CEND0  ((uintptr_t)(delayLine + delayLen));
    aY = AE_ZALIGN128();

    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 3); n++)
    {
        pDr = pDw;
        AE_L32X2X2_IP(d0, d1, pX, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(d2, d3, pX, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d0, d1, pDw, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d2, d3, pDw, 4 * sizeof(int32_t));

        pH = (const ae_int16x4 *)h;

        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 8 * sizeof(int32_t)); S0 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(int32_t)); S1 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(int32_t)); S2 = pDr;

        q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();
        q8 = q9 = qa = qb = qc = qd = qe = qf = AE_ZERO64();
        q10 = q11 = q12 = q13 = q14 = q15 = q16 = q17 = AE_ZERO64();

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            h2 = AE_L16X4_X(pH, 2 * M * sizeof(int16_t));
            h1 = AE_L16X4_X(pH, 1 * M * sizeof(int16_t));
            AE_L16X4_IP(h0, pH, 4 * sizeof(int16_t));

            AE_L32X2X2_XC(d0, d1, S0, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d2, d3, S1, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d4, d5, S2, 4 * sizeof(int32_t));

            AE_MULA2Q32X16_FIR_L(q0, q1, d0, d1, d2, h0);
            AE_MULA2Q32X16_FIR_L(q2, q3, d1, d2, d3, h0);
            AE_MULA2Q32X16_FIR_L(q4, q5, d2, d3, d4, h0);
            AE_MULA2Q32X16_FIR_L(q6, q7, d3, d4, d5, h0);

            AE_MULA2Q32X16_FIR_L(q8, q9, d0, d1, d2, h1);
            AE_MULA2Q32X16_FIR_L(qa, qb, d1, d2, d3, h1);
            AE_MULA2Q32X16_FIR_L(qc, qd, d2, d3, d4, h1);
            AE_MULA2Q32X16_FIR_L(qe, qf, d3, d4, d5, h1);

            AE_MULA2Q32X16_FIR_L(q10, q11, d0, d1, d2, h2);
            AE_MULA2Q32X16_FIR_L(q12, q13, d1, d2, d3, h2);
            AE_MULA2Q32X16_FIR_L(q14, q15, d2, d3, d4, h2);
            AE_MULA2Q32X16_FIR_L(q16, q17, d3, d4, d5, h2);
        }

        d0 = AE_ROUND32X2F48SASYM(AE_SLAI64(q0, 1), AE_SLAI64(q8, 1));
        AE_PKSR32(d1, q10, 1); AE_PKSR32(d1, q1, 1);
        AE_MUL2P32X4S(d0, d1, d0, d1, 3, 3);
        AE_SA32X2X2_IP(d0, d1, aY, pY);

        d0 = AE_ROUND32X2F48SASYM(AE_SLAI64(q9, 1), AE_SLAI64(q11, 1));
        AE_PKSR32(d1, q2, 1); AE_PKSR32(d1, qa, 1);
        AE_MUL2P32X4S(d0, d1, d0, d1, 3, 3);
        AE_SA32X2X2_IP(d0, d1, aY, pY);

        d0 = AE_ROUND32X2F48SASYM(AE_SLAI64(q12, 1), AE_SLAI64(q3, 1));
        AE_PKSR32(d1, qb, 1); AE_PKSR32(d1, q13, 1);
        AE_MUL2P32X4S(d0, d1, d0, d1, 3, 3);
        AE_SA32X2X2_IP(d0, d1, aY, pY);

        d0 = AE_ROUND32X2F48SASYM(AE_SLAI64(q4, 1), AE_SLAI64(qc, 1));
        AE_PKSR32(d1, q14, 1); AE_PKSR32(d1, q5, 1);
        AE_MUL2P32X4S(d0, d1, d0, d1, 3, 3);
        AE_SA32X2X2_IP(d0, d1, aY, pY);

        d0 = AE_ROUND32X2F48SASYM(AE_SLAI64(qd, 1), AE_SLAI64(q15, 1));
        AE_PKSR32(d1, q6, 1); AE_PKSR32(d1, qe, 1);
        AE_MUL2P32X4S(d0, d1, d0, d1, 3, 3);
        AE_SA32X2X2_IP(d0, d1, aY, pY);

        d0 = AE_ROUND32X2F48SASYM(AE_SLAI64(q16, 1), AE_SLAI64(q7, 1));
        AE_PKSR32(d1, qf, 1); AE_PKSR32(d1, q17, 1);
        AE_MUL2P32X4S(d0, d1, d0, d1, 3, 3);
        AE_SA32X2X2_IP(d0, d1, aY, pY);
    }
    AE_SA128POS_FP(aY, pY);
    return (int)((int32_t *)pDw - delayLine);
}
