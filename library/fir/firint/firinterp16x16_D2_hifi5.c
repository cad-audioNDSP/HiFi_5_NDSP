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
    Interpolating block real FIR filter, 16x16-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "firinterp16x16_common.h"

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

/* Data processing function for a factor 2 interpolating FIR filter. */
int firinterp16x16_D2_proc( int16_t * restrict y,
                            int16_t * delayLine, int delayLen,
                      const int16_t * restrict x,
                      const int16_t * restrict h,
                      int wrIx, int D, int N, int M )
{
    const ae_int16x4 *          pX;
          ae_int16x4 * restrict pDw;
    const ae_int16x4 *          pDr;
    const ae_int16x4 *          S0;
    const ae_int16x4 *          S1;
    const ae_int16x4 *          pH;
          ae_int16x4 * restrict pY;

    ae_valign aY;

    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
    ae_int16x4 d0, d1, d2;
    ae_int16x4 h0, h1;
    ae_int32x2 t0, t1;

    int m, n;

    NASSERT(y && delayLine && x && h);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(h, 16);
    NASSERT_ALIGN(delayLine, 16);
    NASSERT(D == 2);
    NASSERT(M > 0 && M % 4 == 0);
    NASSERT(N > 0 && N % 8 == 0);

    //
    // Setup pointers and circular delay line buffer.
    //
    pX  = (const ae_int16x4 *)x;
    pY  = (      ae_int16x4 *)y;
    pDw = (      ae_int16x4 *)(delayLine + wrIx);
    WUR_AE_CBEGIN0((uintptr_t)(delayLine));
    WUR_AE_CEND0  ((uintptr_t)(delayLine + delayLen));
    aY = AE_ZALIGN64();

    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 3); n++)
    {
        pDr = pDw;
        // Load 8 input samples.
        AE_L16X4_IP(d0, pX, 4 * sizeof(int16_t));
        AE_L16X4_IP(d1, pX, 4 * sizeof(int16_t));
        // Store 8 samples to the delay line buffer with circular address update.
        AE_S16X4_XC(d0, pDw, 4 * sizeof(int16_t));
        AE_S16X4_XC(d1, pDw, 4 * sizeof(int16_t));

        // Reset the coefficients pointer. Now it looks at the tap corresponding
        // to the oldest sample in the delay line.
        pH = (const ae_int16x4 *)h;

        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 8 * sizeof(int16_t));
        S0 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(int16_t));
        S1 = pDr;

        q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();
        q8 = q9 = qa = qb = qc = qd = qe = qf = AE_ZERO64();

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L16X4_XC(d0, S0, 4 * sizeof(int16_t));
            AE_L16X4_XC(d1, S1, 4 * sizeof(int16_t));
            d2 = AE_L16X4_I(S1, 0);

            h1 = AE_L16X4_X(pH, M * sizeof(int16_t));
            AE_L16X4_IP(h0, pH, 4 * sizeof(int16_t));

            AE_MULAFQ16X2_FIR_2(q0, q1, d0, d1, h0);
            AE_MULAFQ16X2_FIR_0(q2, q3, d0, d1, h0);
            AE_MULAFQ16X2_FIR_2(q4, q5, d1, d2, h0);
            AE_MULAFQ16X2_FIR_0(q6, q7, d1, d2, h0);

            AE_MULAFQ16X2_FIR_2(q8, q9, d0, d1, h1);
            AE_MULAFQ16X2_FIR_0(qa, qb, d0, d1, h1);
            AE_MULAFQ16X2_FIR_2(qc, qd, d1, d2, h1);
            AE_MULAFQ16X2_FIR_0(qe, qf, d1, d2, h1);
        }

        t0 = AE_TRUNCA32X2F64S(q0, q8, 33);
        t1 = AE_TRUNCA32X2F64S(q1, q9, 33);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
        t0 = AE_TRUNCA32X2F64S(q2, qa, 33);
        t1 = AE_TRUNCA32X2F64S(q3, qb, 33);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
        t0 = AE_TRUNCA32X2F64S(q4, qc, 33);
        t1 = AE_TRUNCA32X2F64S(q5, qd, 33);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
        t0 = AE_TRUNCA32X2F64S(q6, qe, 33);
        t1 = AE_TRUNCA32X2F64S(q7, qf, 33);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
    }
    AE_SA64POS_FP(aY, pY);
    return (int)((int16_t *)pDw - delayLine);
}
