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
    Real block FIR filter, 32x32-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "bkfir32x32_common.h"

/*-------------------------------------------------------------------------
  Real FIR filter.
  Computes a real FIR filter (direct-form) using IR stored in vector h. The 
  real data input is stored in vector x. The filter output result is stored 
  in vector y. The filter calculates N output samples using M coefficients 
  and requires last M-1 samples in the delay line.
  NOTE: 
  1. User application is not responsible for management of delay lines
  2. User has an option to set IR externally or copy from original location 
     (i.e. from the slower constant memory). In the first case, user is 
     responsible for right alignment, ordering and zero padding of filter 
     coefficients - usually array is composed from zeroes (left padding), 
     reverted IR and right zero padding.


  Precision: 
  16x16    16-bit data, 16-bit coefficients, 16-bit outputs. Ordinary variant 
           and stereo
  32x16    32-bit data, 16-bit coefficients, 32-bit outputs
  32x32    32-bit data, 32-bit coefficients, 32-bit outputs. Ordinary variant 
           and stereo
  32x32ep  the same as above but using 72-bit accumulator for intermediate 
           computations
  f        floating point. Ordinary variant and stereo

  Input:
  x[N*S]   input samples, Q31, Q15, floating point
  h[M]     filter coefficients in normal order, Q31, Q15, floating point
  hl[M]    for stereo filters: filter coefficients for left channel
  hr[M]    for stereo filters: filter coefficients for right channel
  N        length of sample block, should be a multiple of 4
  M        length of filter, should be a multiple of 4
  extIR    if zero, IR is copied from original location, otherwise not
           but user should keep alignment, order of coefficients 
           and zero padding requirements shown below
  S        1 for ordinary (single channel) filters, 2 - for stereo variant
  
  Output:
  y[N*S]   output samples, Q31, Q15, floating point

  Alignment, ordering and zero padding for external IR  (extIR!=0)
  ------------------------+----------+--------------+--------------+----------------
  Function                |Alignment,|Left zero     |   Coefficient| Right zero 
                          | bytes    |padding, bytes|   order      | padding, bytes
  ------------------------+----------+--------------+--------------+----------------
  bkfir16x16_init         |    16    |      2       |  inverted    |  6
  bkfir32x16_init         |    16    |      2       |  inverted    |  6
  bkfir32x32_init         |    16    |      4       |  inverted    |  12
  bkfir32x32ep_init       |    16    |      4       |  inverted    |  12
  bkfirf_init             |    16    |      4       |  inverted    |  12
  stereo_bkfir16x16_init  |    16    |      2       |  inverted    |  6
  stereo_bkfir32x32_init  |    16    |      4       |  inverted    |  12
  stereo_bkfirf_init      |    16    |      4       |  inverted    |  12
  ------------------------+----------+--------------+--------------+----------------

  Restrictions:
  x, y     should not be overlapping
  x, h     aligned on a 16-bytes boundary
  N, M     multiples of 4 
-------------------------------------------------------------------------*/

void bkfir32x32_process( bkfir32x32_handle_t handle,
                         int32_t * restrict  y,
                   const int32_t * restrict  x, int N )
{
    bkfir32x32_t * bkfir = (bkfir32x32_ptr_t)handle;

    const ae_int32x4 * restrict pX;
    const ae_int32x4 * restrict pDr;
          ae_int32x4 * restrict pDw;
    const ae_int32x4 * restrict S0;
    const ae_int32x4 * restrict S1;
    const ae_int32x4 * restrict S2;
    const ae_int32x2 * restrict pH;
          ae_int32x2 * restrict pY;
    ae_valign aY;
    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_int32x2 d0, d1, d2, d3, d4, d5;
    ae_int32x2 h0, h1;

    int M;
    int n, m;

    M = bkfir->M;
    NASSERT(bkfir && bkfir->magic == BKFIR32X32_MAGIC && y && x);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN((bkfir->coef), 16);
    NASSERT_ALIGN((bkfir->delayLine), 16);
    NASSERT_ALIGN((bkfir->delayLine + bkfir->delayLen), 16);
    NASSERT_ALIGN((bkfir->delayPos), 16);
    NASSERT(N % 4 == 0 && M % 4 == 0);
    if (N <= 0) return;

    // Setup pointers and circular delay line buffer.
    pX = (const ae_int32x4 *)x;
    pY = (      ae_int32x2 *)y;
    pDw= (      ae_int32x4 *)bkfir->delayPos;
    WUR_AE_CBEGIN0((uintptr_t)(bkfir->delayLine));
    WUR_AE_CEND0  ((uintptr_t)(bkfir->delayLine + bkfir->delayLen));
    aY = AE_ZALIGN64();

    for (n = 0; n < (N >> 3); n++)
    {
        AE_L32X2X2_IP(d0, d1, pX, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(d2, d3, pX, 4 * sizeof(int32_t));
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
    if (N & 4)
    {
        AE_L32X2X2_IP(d0, d1, pX, 4 * sizeof(int32_t));
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

    bkfir->delayPos = (int32_t*)pDw;
} // bkfir32x32_process()
