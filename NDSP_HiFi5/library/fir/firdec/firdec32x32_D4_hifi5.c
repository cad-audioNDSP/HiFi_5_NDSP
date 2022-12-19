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
    Decimating block real FIR filter, 32x32-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "firdec32x32_common.h"

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

int fir32x32_D4_proc(int32_t * restrict y,
                     int32_t * delayLine, int delayLen,
               const int32_t * restrict x,
               const int32_t * restrict h,
               int wrIx, int D, int N, int M)
{
    const ae_int32x4  *          pX;
          ae_int32x4  * restrict pDw;
    const ae_int64    *          pDr;
    const ae_int32x4  *          pD0;
    const ae_int32x4  *          pD1;
    const ae_int32x4  *          pD2;
    const ae_int32x4  *          pD3;
    const ae_int32x4  *          pH;
          ae_int32x4  * restrict pY;

    ae_valignx2 aY;

    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_int32x2 d0, d1, d2, d3, d4, d5, d6, d7;
    ae_int32x2 d8, d9, da, db, dc, dd, de, df, d10, d11;
    ae_int32x2 h0, h1;

    int m, n;

    NASSERT(y && delayLine && x && h);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(h, 16);
    NASSERT_ALIGN(delayLine, 16);
    NASSERT(D == 4);
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
    // Break the input signal into 8*D-samples blocks. For each block, store
    // 8*D samples to the delay line buffer, and compute 8 samples of decimated
    // response signal.
    //
    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 3); n++)
    {
        pDr = (const ae_int64 *)pDw;
        AE_L32X2X2_IP(d0, d1, pX, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(d2, d3, pX, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(d4, d5, pX, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(d6, d7, pX, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d0, d1, pDw, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d2, d3, pDw, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d4, d5, pDw, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d6, d7, pDw, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(d0, d1, pX, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(d2, d3, pX, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(d4, d5, pX, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(d6, d7, pX, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d0, d1, pDw, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d2, d3, pDw, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d4, d5, pDw, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d6, d7, pDw, 4 * sizeof(int32_t));

        pH = (const ae_int32x4 *)h;
        AE_ADDCIRC_XC(pDr, 32 * sizeof(int32_t)); pD0 = (const ae_int32x4 *)pDr;
        AE_ADDCIRC_XC(pDr,  8 * sizeof(int32_t)); pD1 = (const ae_int32x4 *)pDr;
        AE_ADDCIRC_XC(pDr,  8 * sizeof(int32_t)); pD2 = (const ae_int32x4 *)pDr;
        AE_ADDCIRC_XC(pDr,  8 * sizeof(int32_t)); pD3 = (const ae_int32x4 *)pDr;

        q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

        AE_L32X2X2_XC(d0, d1, pD0, 4 * sizeof(int32_t));
        AE_L32X2X2_XC(d2, d3, pD0, 4 * sizeof(int32_t));
        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 3) + 1; m++)
        {
            AE_L32X2X2_XC(d4, d5, pD1, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d6, d7, pD1, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d8, d9, pD2, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(da, db, pD2, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(dc, dd, pD3, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(de, df, pD3, 4 * sizeof(int32_t));
            AE_L32X2X2_I(d10, d11, pD3, 0);

            AE_L32X2X2_IP(h0, h1, pH, 4 * sizeof(int32_t));
            AE_MULAAF2D32RA_HH_LL(q0, q1, d0, d2, h0, h0);
            AE_MULAAF2D32RA_HH_LL(q2, q3, d4, d6, h0, h0);
            AE_MULAAF2D32RA_HH_LL(q4, q5, d8, da, h0, h0);
            AE_MULAAF2D32RA_HH_LL(q6, q7, dc, de, h0, h0);
            AE_MULAAF2D32RA_HH_LL(q0, q1, d1, d3, h1, h1);
            AE_MULAAF2D32RA_HH_LL(q2, q3, d5, d7, h1, h1);
            AE_MULAAF2D32RA_HH_LL(q4, q5, d9, db, h1, h1);
            AE_MULAAF2D32RA_HH_LL(q6, q7, dd, df, h1, h1);

            AE_L32X2X2_IP(h0, h1, pH, 4 * sizeof(int32_t));
            AE_MULAAF2D32RA_HH_LL(q0, q1, d2, d4, h0, h0);
            AE_MULAAF2D32RA_HH_LL(q2, q3, d6, d8, h0, h0);
            AE_MULAAF2D32RA_HH_LL(q4, q5, da, dc, h0, h0);
            AE_MULAAF2D32RA_HH_LL(q6, q7, de, d10, h0, h0);
            AE_MULAAF2D32RA_HH_LL(q0, q1, d3, d5, h1, h1);
            AE_MULAAF2D32RA_HH_LL(q2, q3, d7, d9, h1, h1);
            AE_MULAAF2D32RA_HH_LL(q4, q5, db, dd, h1, h1);
            AE_MULAAF2D32RA_HH_LL(q6, q7, df, d11, h1, h1);

            AE_L32X2X2_XC(d0, d1, pD0, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d2, d3, pD0, 4 * sizeof(int32_t));
        }

        d0 = AE_ROUND32X2F48SASYM(q0, q1);
        d1 = AE_ROUND32X2F48SASYM(q2, q3);
        d2 = AE_ROUND32X2F48SASYM(q4, q5);
        d3 = AE_ROUND32X2F48SASYM(q6, q7);
        AE_SA32X2X2_IP(d0, d1, aY, pY);
        AE_SA32X2X2_IP(d2, d3, aY, pY);
    }
    AE_SA128POS_FP(aY, pY);

    return (int)((int32_t *)pDw - delayLine);
}
