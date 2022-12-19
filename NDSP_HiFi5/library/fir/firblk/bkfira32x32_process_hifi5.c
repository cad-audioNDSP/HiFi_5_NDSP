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
    Real block FIR filter, 32x32-bit, unaligned data and arbitrary M/N allowed
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "bkfira32x32_common.h"

#define SMALLER_CODESIZE 1

/*-------------------------------------------------------------------------
  Real FIR filter.
  Computes a real FIR filter (direct-form) using IR stored in vector h. The 
  real data input is stored in vector x. The filter output result is stored 
  in vector y. The filter calculates N output samples using M coefficients 
  and requires last M-1 samples in the delay line.
  These functions implement FIR filter with no limitation on size of data
  block, alignment and length of impulse response at the cost of increased
  processing complexity.
  NOTE: 
  User application is not responsible for management of delay lines.

  Precision: 
  16x16    16-bit data, 16-bit coefficients, 16-bit outputs
  32x16    32-bit data, 16-bit coefficients, 32-bit outputs
  32x32    32-bit data, 32-bit coefficients, 32-bit outputs
  32x32ep  the same as above but using 72-bit accumulator for intermediate 
           computations
  f        floating point
  Input:
  x[N]     input samples, Q15, Q31, floating point
  h[M]     filter coefficients in normal order, Q15, Q31, floating point
  N        length of sample block
  M        length of filter
  Output:
  y[N]     input samples, Q15, Q31, floating point 

  Restrictions:
  x,y      should not be overlapping
-------------------------------------------------------------------------*/

void bkfira32x32_process( bkfira32x32_handle_t handle,
                          int32_t * restrict       y,
                    const int32_t * restrict       x, int N )
{
    bkfira32x32_t * bkfir = (bkfira32x32_ptr_t)handle;

    const ae_int32x4 * restrict pX;
    const ae_int32x4 * restrict pDr;
          ae_int32x4 * restrict pDw;
    const ae_int32x4 * restrict S0;
    const ae_int32x4 * restrict S1;
    const ae_int32x4 * restrict S2;
    const ae_int32x2 * restrict pH;
          ae_int32x2 * restrict pY;
    ae_valignx2 aX, aD0, aD1;
    ae_valign aY;
    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_int32x2 d0, d1, d2, d3, d4, d5;
    ae_int32x2 h0, h1;

    int M;
    int n, m;
    int wrIx, ua_cnt;

    M = bkfir->M;
    wrIx = bkfir->wrIx;

    NASSERT(bkfir && bkfir->magic == BKFIRA32X32_MAGIC && y && x);
    NASSERT_ALIGN((bkfir->coef), 16);
    NASSERT_ALIGN((bkfir->delayLine), 16);
    NASSERT_ALIGN((bkfir->delayLine + bkfir->delayLen), 16);
    if (N <= 0) return;

    // Setup pointers and circular delay line buffer.
    pX = (const ae_int32x4 *)x;
    pY = (      ae_int32x2 *)y;
    pDw= (      ae_int32x4 *)(bkfir->delayLine + wrIx);
    WUR_AE_CBEGIN0((uintptr_t)(bkfir->delayLine));
    WUR_AE_CEND0  ((uintptr_t)(bkfir->delayLine + bkfir->delayLen));

    ua_cnt = (4 - wrIx) & 3;
    ua_cnt = XT_MIN(ua_cnt, N);
    //process first 0..3 samples (until pDw is aligned to 8 bytes)
    if (ua_cnt)
    {
        __Pragma("loop_count min=1, max=3");
        for (n = 0; n < ua_cnt; n++)
        {
            AE_L32_IP(d0, castxcc(ae_int32, pX), sizeof(int32_t));
            AE_S32_L_XC(d0, castxcc(ae_int32, pDw), sizeof(int32_t));
        }

        pH = (const ae_int32x2 *)bkfir->coef;
        pDr = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), (8 - ua_cnt) * sizeof(int32_t));
        S0 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(int32_t));
        S1 = pDr;
        AE_LA32X2X2POS_PC(aD0, S0);
        AE_LA32X2X2POS_PC(aD1, S1);

        q0 = q1 = q2 = q3 = AE_ZERO64();

        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            AE_L32X2_IP(h0, pH, 2 * sizeof(int32_t));
            AE_L32X2_IP(h1, pH, 2 * sizeof(int32_t));
            AE_LA32X2X2_IC(d0, d1, aD0, S0);
            AE_LA32X2X2_IC(d2, d3, aD1, S1);
            AE_MULAFD32X2RA_FIR_H(q0, q1, d0, d1, h0);
            AE_MULAFD32X2RA_FIR_H(q2, q3, d1, d2, h0);
            AE_MULAFD32X2RA_FIR_H(q0, q1, d1, d2, h1);
            AE_MULAFD32X2RA_FIR_H(q2, q3, d2, d3, h1);
        }

        AE_S32RA64S_IP(q0, castxcc(ae_int32, pY), sizeof(int32_t));
        if (ua_cnt > 1) AE_S32RA64S_IP(q1, castxcc(ae_int32, pY), sizeof(int32_t));
        if (ua_cnt > 2) AE_S32RA64S_IP(q2, castxcc(ae_int32, pY), sizeof(int32_t));
        N -= ua_cnt;
    }

    aX = AE_LA128_PP(pX);
    aY = AE_ZALIGN64();
    for (n = 0; n < (N >> 3); n++)
    {
        AE_LA32X2X2_IP(d0, d1, aX, pX);
        AE_LA32X2X2_IP(d2, d3, aX, pX);
        AE_S32X2X2_XC(d0, d1, pDw, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d2, d3, pDw, 4 * sizeof(int32_t));

        pH = (const ae_int32x2 *)bkfir->coef;
        pDr = pDw;
        S0 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(int32_t));
        S1 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(int32_t));
        S2 = pDr;

        q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            AE_L32X2_IP(h0, pH, 2 * sizeof(int32_t));
            AE_L32X2_IP(h1, pH, 2 * sizeof(int32_t));
            AE_L32X2X2_XC(d0, d1, S0, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d2, d3, S1, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d4, d5, S2, 4 * sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_H(q0, q1, d0, d1, h0);
            AE_MULAFD32X2RA_FIR_H(q2, q3, d1, d2, h0);
            AE_MULAFD32X2RA_FIR_H(q4, q5, d2, d3, h0);
            AE_MULAFD32X2RA_FIR_H(q6, q7, d3, d4, h0);
            AE_MULAFD32X2RA_FIR_H(q0, q1, d1, d2, h1);
            AE_MULAFD32X2RA_FIR_H(q2, q3, d2, d3, h1);
            AE_MULAFD32X2RA_FIR_H(q4, q5, d3, d4, h1);
            AE_MULAFD32X2RA_FIR_H(q6, q7, d4, d5, h1);
        }

        d0 = AE_ROUND32X2F48SASYM(q0, q1);
        d1 = AE_ROUND32X2F48SASYM(q2, q3);
        d2 = AE_ROUND32X2F48SASYM(q4, q5);
        d3 = AE_ROUND32X2F48SASYM(q6, q7);
        AE_SA32X2_IP(d0, aY, pY);
        AE_SA32X2_IP(d1, aY, pY);
        AE_SA32X2_IP(d2, aY, pY);
        AE_SA32X2_IP(d3, aY, pY);
    }
#if SMALLER_CODESIZE
    AE_SA64POS_FP(aY, pY);

    for (n = (8 * n); n < N; n++)
    {
        ae_valign aS, aH;

        AE_L32_IP(d0, castxcc(ae_int32, pX), sizeof(int32_t));
        AE_S32_L_XC(d0, castxcc(ae_int32, pDw), sizeof(int32_t));

        pH = (const ae_int32x2 *)(bkfir->coef + 1);
        aH = AE_LA64_PP(pH);
        S0 = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, S0), 8 * sizeof(int32_t));
        AE_LA32X2POS_PC(aS, castxcc(ae_int32x2, S0));

        q0 = q1 = AE_ZERO64();

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_LA32X2_IP(h0, aH, pH);
            AE_LA32X2_IP(h1, aH, pH);
            AE_LA32X2_IC(d0, aS, castxcc(ae_int32x2, S0));
            AE_LA32X2_IC(d1, aS, castxcc(ae_int32x2, S0));
            AE_MULAFD32X2RA_FIR_H(q0, q1, d0, d0, h0);
            AE_MULAFD32X2RA_FIR_H(q0, q1, d1, d1, h1);
        }

        AE_S32RA64S_IP(q0, castxcc(ae_int32, pY), sizeof(int32_t));
    }
#else
    if (N & 4)
    {
        AE_LA32X2X2_IP(d0, d1, aX, pX);
        AE_S32X2X2_XC(d0, d1, pDw, 4 * sizeof(int32_t));

        pH = (const ae_int32x2 *)bkfir->coef;
        pDr = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(int32_t));
        S0 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(int32_t));
        S1 = pDr;

        q0 = q1 = q2 = q3 = AE_ZERO64();

        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            AE_L32X2_IP(h0, pH, 2 * sizeof(int32_t));
            AE_L32X2_IP(h1, pH, 2 * sizeof(int32_t));
            AE_L32X2X2_XC(d0, d1, S0, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d2, d3, S1, 4 * sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_H(q0, q1, d0, d1, h0);
            AE_MULAFD32X2RA_FIR_H(q2, q3, d1, d2, h0);
            AE_MULAFD32X2RA_FIR_H(q0, q1, d1, d2, h1);
            AE_MULAFD32X2RA_FIR_H(q2, q3, d2, d3, h1);
        }

        d0 = AE_ROUND32X2F48SASYM(q0, q1);
        d1 = AE_ROUND32X2F48SASYM(q2, q3);
        AE_SA32X2_IP(d0, aY, pY);
        AE_SA32X2_IP(d1, aY, pY);
    }
    AE_SA64POS_FP(aY, pY);

    N &= 3;
    if (N)
    {
        __Pragma("loop_count min=1, max=3");
        for (n = 0; n < N; n++)
        {
            AE_L32_IP(d0, castxcc(ae_int32, pX), sizeof(int32_t));
            AE_S32_L_XC(d0, castxcc(ae_int32, pDw), sizeof(int32_t));
        }

        pH = (const ae_int32x2 *)bkfir->coef;
        pDr = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), (8 - N) * sizeof(int32_t));
        S0 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(int32_t));
        S1 = pDr;

        q0 = q1 = q2 = q3 = AE_ZERO64();

        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            AE_L32X2_IP(h0, pH, 2 * sizeof(int32_t));
            AE_L32X2_IP(h1, pH, 2 * sizeof(int32_t));
            AE_L32X2X2_XC(d0, d1, S0, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d2, d3, S1, 4 * sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_H(q0, q1, d0, d1, h0);
            AE_MULAFD32X2RA_FIR_H(q2, q3, d1, d2, h0);
            AE_MULAFD32X2RA_FIR_H(q0, q1, d1, d2, h1);
            AE_MULAFD32X2RA_FIR_H(q2, q3, d2, d3, h1);
        }

        AE_S32RA64S_IP(q0, castxcc(ae_int32, pY), sizeof(int32_t));
        if (N > 1) AE_S32RA64S_IP(q1, castxcc(ae_int32, pY), sizeof(int32_t));
        if (N > 2) AE_S32RA64S_IP(q2, castxcc(ae_int32, pY), sizeof(int32_t));
    }
#endif
    bkfir->wrIx = (int)((int32_t *)pDw - bkfir->delayLine);
} // bkfira32x32_process()
