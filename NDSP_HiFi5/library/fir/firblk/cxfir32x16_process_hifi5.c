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
    Complex block FIR filter, 32x16-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "cxfir32x16_common.h"

/*-------------------------------------------------------------------------
  Complex Block FIR Filter
  Computes a complex FIR filter (direct-form) using complex IR stored in 
  vector h. The complex data input is stored in vector x. The filter output
  result is stored in vector r. The filter calculates N output samples using 
  M coefficients and requires last M-1 samples in the delay line. Real and
  imaginary parts are interleaved and real parts go first (at even indexes).
  NOTE: 
  1. User application is not responsible for management of delay lines
  2. User has an option to set IR externally or copy from original location 
     (i.e. from the slower constant memory). In the first case, user is 
     responsible for right alignment, ordering and zero padding of filter 
     coefficients - usually array is composed from zeroes (left padding), 
     reverted IR and right zero padding.

  Precision: 
  16x16     16-bit data, 16-bit coefficients, 16-bit outputs
  32x16     32-bit data, 16-bit coefficients, 32-bit outputs
  32x32     32-bit data, 32-bit coefficients, 32-bit outputs
  32x32ep   the same as above but using 72-bit accumulator for intermediate 
            computations
  f         floating point

  Input:
  h[M]      complex filter coefficients; h[0] is to be multiplied with the 
            newest sample, Q15, Q31, floating point
  x[N]      input samples, Q15, Q31, floating point
  N         length of sample block (in complex samples) 
  M         length of filter 
  extIR     if zero, IR is copied from original location, otherwise not
            but user should keep alignment, order of coefficients 
            and zero padding requirements shown below
  Output:			
  y[N]      output samples, Q15, Q31, floating point

  Alignment, ordering and zero padding for external IR  (extIR!=0)
  -----------------+----------+--------------+--------------+----------------
  Function	       |Alignment,|Left zero     |   Coefficient| Right zero 
                   | bytes    |padding, bytes|   order      | padding, bytes
  -----------------+----------+--------------+--------------+----------------
  cxfir16x16_init, |    16    |  2 before    |  *           |  6 after
                   |          |  each copy   |              |  each copy
  cxfir32x16_init  |    16    |  2 before    |  *           |  6 after
                   |          |  each copy   |              |  each copy
  cxfir32x32_init  |    16    |    8         |  inverted    |  0
  cxfir32x32ep_init|    16    |    0         |  inv,conj    |  0
  cxfirf_init      |    16    |    0         |  direct      |  0
  -----------------+----------+--------------+--------------+----------------
  * inverted: conjugated copy and (imaginary; real) copy at 4*(M+4) bytes offset

  Restriction:
  x,y       should not overlap
  x,h       aligned on a 16-bytes boundary
  N,M       multiples of 4
-------------------------------------------------------------------------*/

void cxfir32x16_process( cxfir32x16_handle_t handle,
                         complex_fract32 * restrict  y,
                   const complex_fract32 * restrict  x, int N )
{
    cxfir32x16_ptr_t cxfir = (cxfir32x16_ptr_t)handle;

    const ae_int32x4 *          pX;
          ae_int32x4 * restrict pDw;
    const ae_int32x4 *          S0;
    const ae_int32x4 *          S1;
    const ae_int16x8 *          pH;
          ae_int32x2 * restrict pY;
    ae_f64 q0r, q1r, q2r, q3r;
    ae_f64 q0i, q1i, q2i, q3i;
    ae_int32x2 d0, d1, d2, d3, d4, d5, d6, d7;
    ae_int16x4 c0, c0i, c1, c1i;

    int M;
    int n, m;

    M = cxfir->M;
    NASSERT(cxfir && cxfir->magic == CXFIR32X16_MAGIC && y && x);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN((cxfir->coef), 16);
    NASSERT_ALIGN((cxfir->delayLine), 16);
    NASSERT_ALIGN((cxfir->delayLine + cxfir->delayLen), 16);
    NASSERT_ALIGN((cxfir->delayPos), 16);
    NASSERT(N % 4 == 0 && M % 4 == 0);
    if (N <= 0) return;

    // Setup pointers and circular delay line buffer.
    pX = (const ae_int32x4 *)x;
    pY = (      ae_int32x2 *)y;
    pDw= (      ae_int32x4 *)(cxfir->delayPos);
    WUR_AE_CBEGIN0((uintptr_t)(cxfir->delayLine));
    WUR_AE_CEND0  ((uintptr_t)(cxfir->delayLine + cxfir->delayLen));

    //
    // Break the input signal into 4-samples blocks. For each block, store 4
    // samples to the delay line and compute the filter response.
    //
    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 2); n++)
    {
        AE_L32X2X2_IP(d0, d1, pX, 2 * sizeof(complex_fract32));
        AE_L32X2X2_IP(d2, d3, pX, 2 * sizeof(complex_fract32));
        AE_S32X2X2_XC(d0, d1, pDw, 2 * sizeof(complex_fract32));
        AE_S32X2X2_XC(d2, d3, pDw, 2 * sizeof(complex_fract32));

        pH = (const ae_int16x8 *)cxfir->coef;
        S0 = pDw;
        S1 = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, S1), 4 * sizeof(complex_fract32));

        q0r = q1r = q2r = q3r = AE_ZERO64();
        q0i = q1i = q2i = q3i = AE_ZERO64();

        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            AE_L32X2X2_XC(d0, d1, S0, 2 * sizeof(complex_fract32));
            AE_L32X2X2_XC(d2, d3, S0, 2 * sizeof(complex_fract32));
            AE_L32X2X2_XC(d4, d5, S1, 2 * sizeof(complex_fract32));
            AE_L32X2X2_XC(d6, d7, S1, 2 * sizeof(complex_fract32));

            AE_L16X4X2_X(c0i, c1i, pH, (M + 4)*sizeof(complex_fract16));
            AE_L16X4X2_IP(c0, c1, pH, 4 * sizeof(complex_fract16));

            AE_MULAAAA2Q32X16(q0r, q0i, d0, d1, c0, c0i);
            AE_MULAAAA2Q32X16(q1r, q1i, d1, d2, c0, c0i);
            AE_MULAAAA2Q32X16(q2r, q2i, d2, d3, c0, c0i);
            AE_MULAAAA2Q32X16(q3r, q3i, d3, d4, c0, c0i);

            AE_MULAAAA2Q32X16(q0r, q0i, d2, d3, c1, c1i);
            AE_MULAAAA2Q32X16(q1r, q1i, d3, d4, c1, c1i);
            AE_MULAAAA2Q32X16(q2r, q2i, d4, d5, c1, c1i);
            AE_MULAAAA2Q32X16(q3r, q3i, d5, d6, c1, c1i);
        }

        q0r = AE_SLAI64(q0r, 1); q0i = AE_SLAI64(q0i, 1);
        q1r = AE_SLAI64(q1r, 1); q1i = AE_SLAI64(q1i, 1);
        q2r = AE_SLAI64(q2r, 1); q2i = AE_SLAI64(q2i, 1);
        q3r = AE_SLAI64(q3r, 1); q3i = AE_SLAI64(q3i, 1);
        AE_S32X2_IP(AE_ROUND32X2F48SASYM(q0r, q0i), pY, sizeof(complex_fract32));
        AE_S32X2_IP(AE_ROUND32X2F48SASYM(q1r, q1i), pY, sizeof(complex_fract32));
        AE_S32X2_IP(AE_ROUND32X2F48SASYM(q2r, q2i), pY, sizeof(complex_fract32));
        AE_S32X2_IP(AE_ROUND32X2F48SASYM(q3r, q3i), pY, sizeof(complex_fract32));
    }

    cxfir->delayPos = (complex_fract32*)pDw;
} // cxfir32x16_process()
