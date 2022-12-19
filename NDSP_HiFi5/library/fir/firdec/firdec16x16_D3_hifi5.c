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
    Decimating block real FIR filter, 16x16-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "firdec16x16_common.h"

/*-----------------------------------------------------------------------------
 * Data processing function of a particular decimating filter. Stores a
 * block of input samples to the circular delay line buffer and computes
 * decimating FIR filter's response.
 * Input:
 *   delayLine - circular delay line buffer start address
 *   delayLen  - Delay line buffer length
 *   wrIx    - next position in the buffer to be filled with an input sample
 *   x[N*D]  - input samples
 *   h[]     - decimating FIR filter coefficients, array layout varies
 * Output:
 *   y[N]    - output samples
 *   retval  - updated index of the oldest sample
 * Notes and restrictions:
 *   1. Most of data processing functions feature a single, hard-coded
 *      decimation factor, so they expect a determined value for parameter D.
 *   2. All pointers with the exception of y[N] must be aligned on an 16-bytes
 *      boundary.
 *   3. N - must be a multiple of 8.
 *   4. M - must be a multiple of 8.
 -----------------------------------------------------------------------------*/

int fir16x16_D3_proc(int16_t * restrict y,
                     int16_t * delayLine, int delayLen,
               const int16_t * restrict x,
               const int16_t * restrict h,
               int wrIx, int D, int N, int M)
{
#if 0
    const ae_int16x4  *          pX;
          ae_int16x4  * restrict pDw;
    const ae_int16x8  *          pDr;
    const ae_int16x8  *          pD0;
    const ae_int16x8  *          pD1;
    const ae_int16x8  *          pD2;
    const ae_int16x8  *          pD3;
    const ae_int16x8  *          pH;
          ae_int16x4  * restrict pY;

    ae_valignx2 aD0, aD1, aD2, aD3;
    ae_valign aY;

    ae_f64 q0, q1, q2, q3;
    ae_int16x4 d0, d1, d2, d3, d4, d5, d6, d7;
    ae_int16x4 h0, h1;
    ae_int32x2 t0, t1;

    int m, n;

    NASSERT(y && delayLine && x && h);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(h, 16);
    NASSERT_ALIGN(delayLine, 16);
    NASSERT(D == 3);
    NASSERT(M > 0 && M % 8 == 0);
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

    //
    // Break the input signal into 4*D-samples blocks. For each block, store
    // 4*D samples to the delay line buffer, and compute 4 samples of decimated
    // response signal.
    //
    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 2); n++)
    {
        pDr = (const ae_int16x8 *)pDw;
        AE_L16X4_IP(d0, pX, 4 * sizeof(int16_t));
        AE_S16X4_XC(d0, pDw, 4 * sizeof(int16_t));
        AE_L16X4_IP(d0, pX, 4 * sizeof(int16_t));
        AE_S16X4_XC(d0, pDw, 4 * sizeof(int16_t));
        AE_L16X4_IP(d0, pX, 4 * sizeof(int16_t));
        AE_S16X4_XC(d0, pDw, 4 * sizeof(int16_t));

        //
        // Setup 4-way delay line reading pointers, one per an accumulator.
        //
        pH = (const ae_int16x8 *)h;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 24 * sizeof(int16_t)); pD0 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr),  3 * sizeof(int16_t)); pD1 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr),  3 * sizeof(int16_t)); pD2 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr),  3 * sizeof(int16_t)); pD3 = pDr;

        AE_LA16X4X2POS_PC(aD0, pD0);
        AE_LA16X4X2POS_PC(aD1, pD1);
        AE_LA16X4X2POS_PC(aD2, pD2);
        AE_LA16X4X2POS_PC(aD3, pD3);

        q0 = q1 = q2 = q3 = AE_ZERO64();

        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 4) + 1; m++)
        {
            AE_L16X4X2_IP(h0, h1, pH, 8 * sizeof(int16_t));
            AE_LA16X4X2_IC(d0, d4, aD0, pD0);
            AE_LA16X4X2_IC(d1, d5, aD1, pD1);
            AE_LA16X4X2_IC(d2, d6, aD2, pD2);
            AE_LA16X4X2_IC(d3, d7, aD3, pD3);
            AE_MULAAAA2Q16(q0, q1, d0, d1, h0, h0);
            AE_MULAAAA2Q16(q2, q3, d2, d3, h0, h0);
            AE_MULAAAA2Q16(q0, q1, d4, d5, h1, h1);
            AE_MULAAAA2Q16(q2, q3, d6, d7, h1, h1);

            AE_L16X4X2_IP(h0, h1, pH, 8 * sizeof(int16_t));
            AE_LA16X4X2_IC(d0, d4, aD0, pD0);
            AE_LA16X4X2_IC(d1, d5, aD1, pD1);
            AE_LA16X4X2_IC(d2, d6, aD2, pD2);
            AE_LA16X4X2_IC(d3, d7, aD3, pD3);
            AE_MULAAAA2Q16(q0, q1, d0, d1, h0, h0);
            AE_MULAAAA2Q16(q2, q3, d2, d3, h0, h0);
            AE_MULAAAA2Q16(q0, q1, d4, d5, h1, h1);
            AE_MULAAAA2Q16(q2, q3, d6, d7, h1, h1);
        }
        t0 = AE_TRUNCA32X2F64S(q0, q1, 33);
        t1 = AE_TRUNCA32X2F64S(q2, q3, 33);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
    }
    AE_SA64POS_FP(aY, pY);

    return (int)((int16_t *)pDw - delayLine);
#else
    const ae_int16x8  *          pX;
          ae_int16x8  * restrict pDw;
    const ae_int16x8  *          pDr;
    const ae_int16x8  *          pD0;
    const ae_int16x8  *          pH;
    const ae_int16x8  *          pH1;
    const ae_int16x8  *          pH2;
    const ae_int16x8  *          pH3;
          ae_int16x4  * restrict pY;

    ae_valignx2 aH1, aH2, aH3;
    ae_valign aY;

    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_int16x4 d0, d1, d2, d3, d4, d5, d6, d7;
    ae_int16x4 h0_0, h1_0, h0_1, h1_1, h0_2, h1_2, h0_3, h1_3;
    ae_int32x2 t0, t1;

    int m, n;

    NASSERT(y && delayLine && x && h);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(h, 16);
    NASSERT_ALIGN(delayLine, 16);
    NASSERT(D == 3);
    NASSERT(M > 0 && M % 8 == 0);
    NASSERT(N > 0 && N % 8 == 0);

    //
    // Setup pointers and circular delay line buffer.
    //
    pX  = (const ae_int16x8 *)x;
    pY  = (      ae_int16x4 *)y;
    pDw = (      ae_int16x8 *)(delayLine + wrIx);
    WUR_AE_CBEGIN0((uintptr_t)(delayLine));
    WUR_AE_CEND0  ((uintptr_t)(delayLine + delayLen));
    aY = AE_ZALIGN64();

    //
    // Break the input signal into 8*D-samples blocks. For each block, store
    // 8*D samples to the delay line buffer, and compute 8 samples of decimated
    // response signal.
    //
    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 3); n++)
    {
        pDr = (const ae_int16x8 *)pDw;
        AE_L16X4X2_IP(d0, d1, pX, 8 * sizeof(int16_t));
        AE_L16X4X2_IP(d2, d3, pX, 8 * sizeof(int16_t));
        AE_L16X4X2_IP(d4, d5, pX, 8 * sizeof(int16_t));
        AE_S16X4X2_XC(d0, d1, pDw, 8 * sizeof(int16_t));
        AE_S16X4X2_XC(d2, d3, pDw, 8 * sizeof(int16_t));
        AE_S16X4X2_XC(d4, d5, pDw, 8 * sizeof(int16_t));

        pH = (const ae_int16x8 *)h;
        pH1 = (const ae_int16x8 *)(h - 1);
        pH2 = (const ae_int16x8 *)(h - 2);
        pH3 = (const ae_int16x8 *)(h - 3);
        aH1 = AE_LA128_PP(pH1);
        aH2 = AE_LA128_PP(pH2);
        aH3 = AE_LA128_PP(pH3);
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 24 * sizeof(int16_t));

        q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

        AE_L16X4X2_XC(d0, d1, pDr, 8 * sizeof(int16_t));
        AE_L16X4X2_XC(d2, d3, pDr, 8 * sizeof(int16_t)); pD0 = pDr;
        AE_L16X4X2_XC(d4, d5, pDr, 8 * sizeof(int16_t));
        __Pragma("loop_count min=2, factor=2");
        for (m = 0; m < (M >> 3) + 1; m++)
        {
            AE_L16X4X2_IP(h0_0, h1_0, pH, 8 * sizeof(int16_t));
            AE_LA16X4X2_IP(h0_1, h1_1, aH1, pH1);
            AE_LA16X4X2_IP(h0_2, h1_2, aH2, pH2);
            AE_LA16X4X2_IP(h0_3, h1_3, aH3, pH3);

            AE_L16X4X2_XC(d6, d7, pDr, 8 * sizeof(int16_t));

            AE_MULAAAA2Q16(q0, q1, d0, d0, h0_0, h0_3);
            AE_MULAAAA2Q16(q2, q3, d1, d2, h0_2, h0_1);
            AE_MULAAAA2Q16(q4, q5, d3, d3, h0_0, h0_3);
            AE_MULAAAA2Q16(q6, q7, d4, d5, h0_2, h0_1);

            AE_MULAAAA2Q16(q0, q1, d1, d1, h1_0, h1_3);
            AE_MULAAAA2Q16(q2, q3, d2, d3, h1_2, h1_1);
            AE_MULAAAA2Q16(q4, q5, d4, d4, h1_0, h1_3);
            AE_MULAAAA2Q16(q6, q7, d5, d6, h1_2, h1_1);

            d0 = d2; d1 = d3;
            AE_L16X4X2_XC(d2, d3, pD0, 8 * sizeof(int16_t));
            d4 = d6; d5 = d7;
        }
        t0 = AE_TRUNCA32X2F64S(q0, q1, 33);
        t1 = AE_TRUNCA32X2F64S(q2, q3, 33);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
        t0 = AE_TRUNCA32X2F64S(q4, q5, 33);
        t1 = AE_TRUNCA32X2F64S(q6, q7, 33);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
    }
    AE_SA64POS_FP(aY, pY);

    return (int)((int16_t *)pDw - delayLine);
#endif
} /* fir16x16_D3_proc() */
