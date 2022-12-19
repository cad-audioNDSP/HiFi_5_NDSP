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
    Decimating block real FIR filter, 32x16-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "firdec32x16_common.h"

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
 *   4. M - must be a multiple of 4.
 -----------------------------------------------------------------------------*/

/* Generic data processing function for a decimating FIR filter. */
int fir32x16_DX_proc(int32_t * restrict y,
                     int32_t * delayLine, int delayLen,
               const int32_t * restrict x,
               const int16_t * restrict h,
               int wrIx, int D, int N, int M)
{
    const ae_int32x4  *          pX;
          ae_int32x4  * restrict pDw;
    const ae_int64    *          pDr;
    const ae_int32x4  *          pD0;
    const ae_int32x4  *          pD1;
    const ae_int32x4  *          pD2;
    const ae_int32x4  *          pD3;
    const ae_int16x8  *          pH;
          ae_int32x4  * restrict pY;

    ae_valignx2 aD1, aD2, aD3;
    ae_valignx2 aY;

    ae_f64 q0, q1, q2, q3;
    ae_int32x2 d0, d1, d2, d3, d4, d5, d6, d7;
    ae_int16x4 h0, h1;

    int m, n;

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
    pY  = (      ae_int32x4 *)y;
    pDw = (      ae_int32x4 *)(delayLine + wrIx);
    WUR_AE_CBEGIN0((uintptr_t)(delayLine));
    WUR_AE_CEND0  ((uintptr_t)(delayLine + delayLen));
    aY = AE_ZALIGN128();

    //
    // Break the input signal into 4*D-samples blocks. For each block, store
    // 4*D samples to the delay line buffer, and compute 4 samples of decimated
    // response signal.
    //
    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 2); n++)
    {
        // Reset the coefficients pointer. Now it looks at the tap corresponding
        // to the oldest sample in the delay line.
        pH = (const ae_int16x8 *)h;

        __Pragma("loop_count min=5");
        for (m = 0; m < D; m++)
        {
            AE_L32X2X2_IP(d0, d1, pX, 4 * sizeof(int32_t));
            AE_S32X2X2_XC(d0, d1, pDw, 4 * sizeof(int32_t));
        }

        //
        // Setup 4-way delay line reading pointers, one per an accumulator.
        //
        pDr = (const ae_int64 *)pDw;
        AE_ADDCIRC_XC(pDr, 4 * D * sizeof(int32_t)); pD0 = (const ae_int32x4 *)pDr;
        AE_ADDCIRC_XC(pDr,     D * sizeof(int32_t)); pD1 = (const ae_int32x4 *)pDr;
        AE_ADDCIRC_XC(pDr,     D * sizeof(int32_t)); pD2 = (const ae_int32x4 *)pDr;
        AE_ADDCIRC_XC(pDr,     D * sizeof(int32_t)); pD3 = (const ae_int32x4 *)pDr;

        AE_LA32X2X2POS_PC(aD1, pD1);
        AE_LA32X2X2POS_PC(aD2, pD2);
        AE_LA32X2X2POS_PC(aD3, pD3);

        q0 = q1 = q2 = q3 = AE_ZERO64();

        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 4) + 1; m++)
        {
            AE_L16X4X2_IP(h0, h1, pH, 8 * sizeof(int16_t));
            AE_L32X2X2_XC(d0, d4, pD0, 4 * sizeof(int32_t));
            AE_LA32X2X2_IC(d1, d5, aD1, pD1);
            AE_LA32X2X2_IC(d2, d6, aD2, pD2);
            AE_LA32X2X2_IC(d3, d7, aD3, pD3);
            AE_MULAAAAFQ32X16(q0, d0, d4, h0);
            AE_MULAAAAFQ32X16(q1, d1, d5, h0);
            AE_MULAAAAFQ32X16(q2, d2, d6, h0);
            AE_MULAAAAFQ32X16(q3, d3, d7, h0);
            AE_L32X2X2_XC(d0, d4, pD0, 4 * sizeof(int32_t));
            AE_LA32X2X2_IC(d1, d5, aD1, pD1);
            AE_LA32X2X2_IC(d2, d6, aD2, pD2);
            AE_LA32X2X2_IC(d3, d7, aD3, pD3);
            AE_MULAAAAFQ32X16(q0, d0, d4, h1);
            AE_MULAAAAFQ32X16(q1, d1, d5, h1);
            AE_MULAAAAFQ32X16(q2, d2, d6, h1);
            AE_MULAAAAFQ32X16(q3, d3, d7, h1);

            AE_L16X4X2_IP(h0, h1, pH, 8 * sizeof(int16_t));
            AE_L32X2X2_XC(d0, d4, pD0, 4 * sizeof(int32_t));
            AE_LA32X2X2_IC(d1, d5, aD1, pD1);
            AE_LA32X2X2_IC(d2, d6, aD2, pD2);
            AE_LA32X2X2_IC(d3, d7, aD3, pD3);
            AE_MULAAAAFQ32X16(q0, d0, d4, h0);
            AE_MULAAAAFQ32X16(q1, d1, d5, h0);
            AE_MULAAAAFQ32X16(q2, d2, d6, h0);
            AE_MULAAAAFQ32X16(q3, d3, d7, h0);
            AE_L32X2X2_XC(d0, d4, pD0, 4 * sizeof(int32_t));
            AE_LA32X2X2_IC(d1, d5, aD1, pD1);
            AE_LA32X2X2_IC(d2, d6, aD2, pD2);
            AE_LA32X2X2_IC(d3, d7, aD3, pD3);
            AE_MULAAAAFQ32X16(q0, d0, d4, h1);
            AE_MULAAAAFQ32X16(q1, d1, d5, h1);
            AE_MULAAAAFQ32X16(q2, d2, d6, h1);
            AE_MULAAAAFQ32X16(q3, d3, d7, h1);
        }
        d0 = AE_ROUND32X2F48SASYM(q0, q1);
        d1 = AE_ROUND32X2F48SASYM(q2, q3);
        AE_SA32X2X2_IP(d0, d1, aY, pY);
    }
    AE_SA128POS_FP(aY, pY);

    return (int)((int32_t *)pDw - delayLine);
}
