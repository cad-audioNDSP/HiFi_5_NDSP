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
    Interpolating block real FIR filter, 32x32-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "firinterp32x32_common.h"

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

/* Generic data processing function for a interpolating FIR filter. */
int firinterp32x32_DX_proc( int32_t * restrict y,
                            int32_t * delayLine, int delayLen,
                      const int32_t * restrict x,
                      const int32_t * restrict h,
                      int wrIx, int D, int N, int M )
{
    const ae_int32x4 * restrict pX;
          ae_int32x4 * restrict pDw;
    const ae_int32x4 * restrict pDr;
    const ae_int32x4 * restrict S0;
    const ae_int32x4 * restrict S1;
    const ae_int32x4 * restrict S2;
    const ae_int32x2 * restrict pH;
          ae_int32   * restrict pY;

    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_int32x2 d0, d1, d2, d3, d4, d5;
    ae_int32x2 h0, h1;

    int m, n, d;

    NASSERT(y && delayLine && x && h);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(h, 16);
    NASSERT_ALIGN(delayLine, 16);
    NASSERT(D > 1);
    NASSERT(M > 0 && M % 4 == 0);
    NASSERT(N > 0 && N % 8 == 0);

    //
    // Setup pointers and circular delay line buffer.
    //
    pX  = (const ae_int32x4 *)x;
    pY  = (      ae_int32   *)y;
    pDw = (      ae_int32x4 *)(delayLine + wrIx);
    WUR_AE_CBEGIN0((uintptr_t)(delayLine));
    WUR_AE_CEND0  ((uintptr_t)(delayLine + delayLen));

    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 3); n++)
    {
        AE_L32X2X2_IP(d0, d1, pX, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(d2, d3, pX, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d0, d1, pDw, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d2, d3, pDw, 4 * sizeof(int32_t));

        pH = (const ae_int32x2 *)h;

        __Pragma("loop_count min=1");
        for (d = 0; d < D; d++)
        {
            pDr = pDw;
            S0 = pDr;
            AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(int32_t)); S1 = pDr;
            AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(int32_t)); S2 = pDr;

            q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

            __Pragma("loop_count min=1");
            for (m = 0; m < (M >> 2); m++)
            {
                AE_L32X2_IP(h0, pH, 2 * sizeof(int32_t));
                AE_L32X2_IP(h1, pH, 2 * sizeof(int32_t));
                AE_L32X2X2_XC(d0, d1, S0, 4 * sizeof(int32_t));
                AE_L32X2X2_XC(d2, d3, S1, 4 * sizeof(int32_t));
                AE_L32X2X2_XC(d4, d5, S2, 4 * sizeof(int32_t));
                AE_MULAFD32X2RA_FIR_L(q0, q1, d0, d1, h0);
                AE_MULAFD32X2RA_FIR_L(q2, q3, d1, d2, h0);
                AE_MULAFD32X2RA_FIR_L(q4, q5, d2, d3, h0);
                AE_MULAFD32X2RA_FIR_L(q6, q7, d3, d4, h0);
                AE_MULAFD32X2RA_FIR_L(q0, q1, d1, d2, h1);
                AE_MULAFD32X2RA_FIR_L(q2, q3, d2, d3, h1);
                AE_MULAFD32X2RA_FIR_L(q4, q5, d3, d4, h1);
                AE_MULAFD32X2RA_FIR_L(q6, q7, d4, d5, h1);
            }

            d0 = AE_ROUND32X2F48SASYM(q0, q1);
            d1 = AE_ROUND32X2F48SASYM(q2, q3);
            AE_MUL2P32X4S(d0, d1, d0, d1, D, D);
            AE_S32_H_XP(d0, pY, D*sizeof(int32_t));
            AE_S32_L_XP(d0, pY, D*sizeof(int32_t));
            AE_S32_H_XP(d1, pY, D*sizeof(int32_t));
            AE_S32_L_XP(d1, pY, D*sizeof(int32_t));

            d0 = AE_ROUND32X2F48SASYM(q4, q5);
            d1 = AE_ROUND32X2F48SASYM(q6, q7);
            AE_MUL2P32X4S(d0, d1, d0, d1, D, D);
            AE_S32_H_XP(d0, pY, D*sizeof(int32_t));
            AE_S32_L_XP(d0, pY, D*sizeof(int32_t));
            AE_S32_H_XP(d1, pY, D*sizeof(int32_t));
            AE_S32_L_XP(d1, pY, (1 - 7 * D)*sizeof(int32_t));
        }
        pY = (ae_int32 *)XT_ADDX4(7 * D, (uintptr_t)pY);
    }
    return (int)((int32_t *)pDw - delayLine);
}
