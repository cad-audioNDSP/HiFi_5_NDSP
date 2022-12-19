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
    Real block FIR filter, 32x16-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "bkfir32x16_common.h"

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

void bkfir32x16_process( bkfir32x16_handle_t handle,
                         int32_t * restrict  y,
                   const int32_t * restrict  x, int N)
{
    bkfir32x16_ptr_t bkfir = (bkfir32x16_ptr_t)handle;

    const ae_int32x4 * restrict pX;
    const ae_int32x4 * restrict pDr;
          ae_int32x4 * restrict pDw;
    const ae_int32x4 * restrict S0;
    const ae_int32x4 * restrict S1;
    const ae_int32x4 * restrict S2;
    const ae_int32x4 * restrict S3;
    const ae_int32x4 * restrict S4;
    const ae_int16x4 * restrict pH;
          ae_int32x4 * restrict pY;
    ae_valignx2 aY;
    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
    ae_int32x2 d0, d1, d2, d3, d4, d5;
    ae_int32x2 d6, d7, d8, d9;
    ae_int16x4 h0;

    int M;
    int n, m;

    M = bkfir->M;
    NASSERT(bkfir && bkfir->magic == BKFIR32X16_MAGIC && y && x);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN((bkfir->coef), 16);
    NASSERT_ALIGN((bkfir->delayLine), 16);
    NASSERT_ALIGN((bkfir->delayLine + bkfir->delayLen), 16);
    NASSERT_ALIGN((bkfir->delayPos), 16);
    NASSERT(N % 4 == 0 && M % 4 == 0);
    if (N <= 0) return;

    // Setup pointers and circular delay line buffer.
    pX = (const ae_int32x4 *)x;
    pY = (      ae_int32x4 *)y;
    pDw= (      ae_int32x4 *)bkfir->delayPos;
    WUR_AE_CBEGIN0((uintptr_t)(bkfir->delayLine));
    WUR_AE_CEND0  ((uintptr_t)(bkfir->delayLine + bkfir->delayLen));
    aY = AE_ZALIGN128();

    for (n = 0; n < (N >> 4); n++)
    {
        AE_L32X2X2_IP(d0, d1, pX, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(d2, d3, pX, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(d4, d5, pX, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(d6, d7, pX, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d0, d1, pDw, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d2, d3, pDw, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d4, d5, pDw, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d6, d7, pDw, 4 * sizeof(int32_t));

        pH = (const ae_int16x4 *)bkfir->coef;
        AE_L16X4_IP(h0, pH, 4 * sizeof(int16_t));
        pDr = pDw;
        AE_L32X2X2_XC(d0, d1, pDr, 4 * sizeof(int32_t)); S0 = pDr;
        AE_L32X2X2_XC(d2, d3, pDr, 4 * sizeof(int32_t)); S1 = pDr;
        AE_L32X2X2_XC(d4, d5, pDr, 4 * sizeof(int32_t)); S2 = pDr;
        AE_L32X2X2_XC(d6, d7, pDr, 4 * sizeof(int32_t)); S3 = pDr;
        AE_L32X2X2_XC(d8, d9, pDr, 4 * sizeof(int32_t)); S4 = pDr;
        AE_MUL2Q32X16_FIR_H(q0, q1, d0, d1, d2, h0);
        AE_MUL2Q32X16_FIR_H(q2, q3, d1, d2, d3, h0);
        AE_MUL2Q32X16_FIR_H(q4, q5, d2, d3, d4, h0);
        AE_MUL2Q32X16_FIR_H(q6, q7, d3, d4, d5, h0);
        AE_MUL2Q32X16_FIR_H(q8, q9, d4, d5, d6, h0);
        AE_MUL2Q32X16_FIR_H(qa, qb, d5, d6, d7, h0);
        AE_MUL2Q32X16_FIR_H(qc, qd, d6, d7, d8, h0);
        AE_MUL2Q32X16_FIR_H(qe, qf, d7, d8, d9, h0);
        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L16X4_IP(h0, pH, 4 * sizeof(int16_t));
            AE_L32X2X2_XC(d0, d1, S0, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d2, d3, S1, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d4, d5, S2, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d6, d7, S3, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d8, d9, S4, 4 * sizeof(int32_t));
            AE_MULA2Q32X16_FIR_H(q0, q1, d0, d1, d2, h0);
            AE_MULA2Q32X16_FIR_H(q2, q3, d1, d2, d3, h0);
            AE_MULA2Q32X16_FIR_H(q4, q5, d2, d3, d4, h0);
            AE_MULA2Q32X16_FIR_H(q6, q7, d3, d4, d5, h0);
            AE_MULA2Q32X16_FIR_H(q8, q9, d4, d5, d6, h0);
            AE_MULA2Q32X16_FIR_H(qa, qb, d5, d6, d7, h0);
            AE_MULA2Q32X16_FIR_H(qc, qd, d6, d7, d8, h0);
            AE_MULA2Q32X16_FIR_H(qe, qf, d7, d8, d9, h0);
        }

        q0 = AE_SLAI64(q0, 1); q1 = AE_SLAI64(q1, 1);
        d0 = AE_ROUND32X2F48SASYM(q0, q1);
        AE_PKSR32(d1, q2, 1); AE_PKSR32(d1, q3, 1);
        AE_SA32X2X2_IP(d0, d1, aY, pY);
        q4 = AE_SLAI64(q4, 1); q5 = AE_SLAI64(q5, 1);
        d0 = AE_ROUND32X2F48SASYM(q4, q5);
        AE_PKSR32(d1, q6, 1); AE_PKSR32(d1, q7, 1);
        AE_SA32X2X2_IP(d0, d1, aY, pY);
        q8 = AE_SLAI64(q8, 1); q9 = AE_SLAI64(q9, 1);
        d0 = AE_ROUND32X2F48SASYM(q8, q9);
        AE_PKSR32(d1, qa, 1); AE_PKSR32(d1, qb, 1);
        AE_SA32X2X2_IP(d0, d1, aY, pY);
        qc = AE_SLAI64(qc, 1); qd = AE_SLAI64(qd, 1);
        d0 = AE_ROUND32X2F48SASYM(qc, qd);
        AE_PKSR32(d1, qe, 1); AE_PKSR32(d1, qf, 1);
        AE_SA32X2X2_IP(d0, d1, aY, pY);
    }
    if (N & 8)
    {
        AE_L32X2X2_IP(d0, d1, pX, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(d2, d3, pX, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d0, d1, pDw, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d2, d3, pDw, 4 * sizeof(int32_t));

        pH = (const ae_int16x4 *)bkfir->coef;
        pDr = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 8 * sizeof(int32_t)); S0 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(int32_t)); S1 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(int32_t)); S2 = pDr;

        q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            AE_L16X4_IP(h0, pH, 4 * sizeof(int16_t));
            AE_L32X2X2_XC(d0, d1, S0, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d2, d3, S1, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(d4, d5, S2, 4 * sizeof(int32_t));
            AE_MULA2Q32X16_FIR_H(q0, q1, d0, d1, d2, h0);
            AE_MULA2Q32X16_FIR_H(q2, q3, d1, d2, d3, h0);
            AE_MULA2Q32X16_FIR_H(q4, q5, d2, d3, d4, h0);
            AE_MULA2Q32X16_FIR_H(q6, q7, d3, d4, d5, h0);
        }

        q0 = AE_SLAI64(q0, 1); q1 = AE_SLAI64(q1, 1);
        d0 = AE_ROUND32X2F48SASYM(q0, q1);
        AE_PKSR32(d1, q2, 1); AE_PKSR32(d1, q3, 1);
        AE_SA32X2X2_IP(d0, d1, aY, pY);
        q4 = AE_SLAI64(q4, 1); q5 = AE_SLAI64(q5, 1);
        d0 = AE_ROUND32X2F48SASYM(q4, q5);
        AE_PKSR32(d1, q6, 1); AE_PKSR32(d1, q7, 1);
        AE_SA32X2X2_IP(d0, d1, aY, pY);
    }
    if (N & 4)
    {
        AE_L32X2X2_IP(d0, d1, pX, 4 * sizeof(int32_t));
        AE_S32X2X2_XC(d0, d1, pDw, 4 * sizeof(int32_t));

        pH = (const ae_int16x4 *)bkfir->coef;
        pDr = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 12 * sizeof(int32_t));
        AE_L32X2X2_XC(d0, d1, pDr, 4 * sizeof(int32_t));
        q0 = q1 = q2 = q3 = AE_ZERO64();

        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            AE_L16X4_IP(h0, pH, 4 * sizeof(int16_t));
            AE_L32X2X2_XC(d2, d3, pDr, 4 * sizeof(int32_t));
            AE_MULA2Q32X16_FIR_H(q0, q1, d0, d1, d2, h0);
            AE_MULA2Q32X16_FIR_H(q2, q3, d1, d2, d3, h0);
            d0 = d2; d1 = d3;
        }

        q0 = AE_SLAI64(q0, 1); q1 = AE_SLAI64(q1, 1);
        d0 = AE_ROUND32X2F48SASYM(q0, q1);
        q2 = AE_SLAI64(q2, 1); q3 = AE_SLAI64(q3, 1);
        d1 = AE_ROUND32X2F48SASYM(q2, q3);
        AE_SA32X2X2_IP(d0, d1, aY, pY);
    }
    AE_SA128POS_FP(aY, pY);

    bkfir->delayPos = (int32_t*)pDw;
}
