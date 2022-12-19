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
    Complex block FIR filter, 32x32-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "cxfir32x32_common.h"

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

void cxfir32x32_process( cxfir32x32_handle_t handle,
                         complex_fract32 * restrict y,
                   const complex_fract32 * restrict x, int N )
{
    cxfir32x32_ptr_t cxfir = (cxfir32x32_ptr_t)handle;

    const ae_int32x4 * restrict pX;
    const ae_int32x2 * restrict pDr;
          ae_int32x4 * restrict pDw;
    const ae_int32x2 * restrict pH;
          ae_int32x2 * restrict pY;
    ae_f64 q0_re, q1_re, q2_re, q3_re;
    ae_f64 q0_im, q1_im, q2_im, q3_im;
    ae_int32x2 X0, X1, X2, X3, Y;

    int M;
    int n, m;

    M = cxfir->M;
    NASSERT(cxfir && cxfir->magic == CXFIR32X32_MAGIC && y && x);
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
    pDw= (      ae_int32x4 *)cxfir->delayPos;
    WUR_AE_CBEGIN0((uintptr_t)(cxfir->delayLine));
    WUR_AE_CEND0  ((uintptr_t)(cxfir->delayLine + cxfir->delayLen));

#if SMALLER_CODESIZE
    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 2); n++)
    {
        ae_int32x2 X4;

        pDr = (ae_int32x2 *)pDw;
        AE_L32X2X2_IP(X0, X1, pX, 2 * sizeof(complex_fract32));
        AE_L32X2X2_IP(X2, X3, pX, 2 * sizeof(complex_fract32));
        AE_S32X2X2_XC(X0, X1, pDw, 2 * sizeof(complex_fract32));
        AE_S32X2X2_XC(X2, X3, pDw, 2 * sizeof(complex_fract32));

        pH = (const ae_int32x2 *)cxfir->coef;
        /* preload data from x */
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(complex_fract32));
        AE_L32X2X2_XC(X0, X1, castxcc(ae_int32x4, pDr), 2 * sizeof(complex_fract32));
        AE_L32X2X2_XC(X2, X3, castxcc(ae_int32x4, pDr), 2 * sizeof(complex_fract32));
        /* load data from y */
        AE_L32X2_IP(Y, pH, sizeof(complex_fract32));
        /* compute correlation of 8 values */
        AE_MULFC32RA(q0_re, q0_im, X0, Y);
        AE_MULFC32RA(q1_re, q1_im, X1, Y);
        AE_MULFC32RA(q2_re, q2_im, X2, Y);
        AE_MULFC32RA(q3_re, q3_im, X3, Y);
        /* shift input line for the next iteration */
        X0 = X1; X1 = X2; X2 = X3;

        __Pragma("loop_count min=2, factor=2");
        for (m = 0; m < (M >> 1); m++)
        {
            /* load data from x */
            AE_L32X2X2_XC(X3, X4, castxcc(ae_int32x4, pDr), 2 * sizeof(complex_fract32));
            /* load data from y */
            AE_L32X2_IP(Y, pH, sizeof(complex_fract32));
            /* compute correlation of 8 values */
            AE_MULAFC32RA(q0_re, q0_im, X0, Y);
            AE_MULAFC32RA(q1_re, q1_im, X1, Y);
            AE_MULAFC32RA(q2_re, q2_im, X2, Y);
            AE_MULAFC32RA(q3_re, q3_im, X3, Y);
            /* load data from y */
            AE_L32X2_IP(Y, pH, sizeof(complex_fract32));
            /* compute correlation of 8 values */
            AE_MULAFC32RA(q0_re, q0_im, X1, Y);
            AE_MULAFC32RA(q1_re, q1_im, X2, Y);
            AE_MULAFC32RA(q2_re, q2_im, X3, Y);
            AE_MULAFC32RA(q3_re, q3_im, X4, Y);
            /* shift input line for the next iteration */
            X0 = X2; X1 = X3; X2 = X4;
        }
        /* save computed samples */
        X0 = AE_ROUND32X2F48SASYM(q0_re, q0_im);
        X1 = AE_ROUND32X2F48SASYM(q1_re, q1_im);
        X2 = AE_ROUND32X2F48SASYM(q2_re, q2_im);
        X3 = AE_ROUND32X2F48SASYM(q3_re, q3_im);
        AE_S32X2_IP(X0, pY, sizeof(complex_fract32));
        AE_S32X2_IP(X1, pY, sizeof(complex_fract32));
        AE_S32X2_IP(X2, pY, sizeof(complex_fract32));
        AE_S32X2_IP(X3, pY, sizeof(complex_fract32));
    }
#else
    for (n = 0; n < (N >> 3); n++)
    {
        ae_f64 q4_re, q5_re, q6_re, q7_re;
        ae_f64 q4_im, q5_im, q6_im, q7_im;
        ae_int32x2 X4, X5, X6, X7;

        AE_L32X2X2_IP(X0, X1, pX, 2 * sizeof(complex_fract32));
        AE_L32X2X2_IP(X2, X3, pX, 2 * sizeof(complex_fract32));
        AE_L32X2X2_IP(X4, X5, pX, 2 * sizeof(complex_fract32));
        AE_L32X2X2_IP(X6, X7, pX, 2 * sizeof(complex_fract32));
        AE_S32X2X2_XC(X0, X1, pDw, 2 * sizeof(complex_fract32));
        AE_S32X2X2_XC(X2, X3, pDw, 2 * sizeof(complex_fract32));
        AE_S32X2X2_XC(X4, X5, pDw, 2 * sizeof(complex_fract32));
        AE_S32X2X2_XC(X6, X7, pDw, 2 * sizeof(complex_fract32));

        pH = (const ae_int32x2 *)cxfir->coef;
        /* preload data from x */
        pDr = (ae_int32x2 *)pDw;
        AE_L32X2X2_XC(X0, X1, castxcc(ae_int32x4, pDr), 2 * sizeof(complex_fract32));
        AE_L32X2X2_XC(X2, X3, castxcc(ae_int32x4, pDr), 2 * sizeof(complex_fract32));
        AE_L32X2X2_XC(X4, X5, castxcc(ae_int32x4, pDr), 2 * sizeof(complex_fract32));
        AE_L32X2X2_XC(X6, X7, castxcc(ae_int32x4, pDr), 2 * sizeof(complex_fract32));
        /* load data from y */
        AE_L32X2_IP(Y, pH, sizeof(complex_fract32));
        /* compute correlation of 8 values */
        AE_MULFC32RA(q0_re, q0_im, X0, Y);
        AE_MULFC32RA(q1_re, q1_im, X1, Y);
        AE_MULFC32RA(q2_re, q2_im, X2, Y);
        AE_MULFC32RA(q3_re, q3_im, X3, Y);
        AE_MULFC32RA(q4_re, q4_im, X4, Y);
        AE_MULFC32RA(q5_re, q5_im, X5, Y);
        AE_MULFC32RA(q6_re, q6_im, X6, Y);
        AE_MULFC32RA(q7_re, q7_im, X7, Y);
        /* shift input line for the next iteration */
        X0 = X1; X1 = X2; X2 = X3;
        X3 = X4; X4 = X5; X5 = X6; X6 = X7;

        __Pragma("loop_count min=4, factor=4");
        for (m = 0; m < M; m++)
        {
            /* load data from x */
            AE_L32X2_XC(X7, pDr, sizeof(complex_fract32));
            /* load data from y */
            AE_L32X2_IP(Y, pH, sizeof(complex_fract32));
            /* compute correlation of 8 values */
            AE_MULAFC32RA(q0_re, q0_im, X0, Y);
            AE_MULAFC32RA(q1_re, q1_im, X1, Y);
            AE_MULAFC32RA(q2_re, q2_im, X2, Y);
            AE_MULAFC32RA(q3_re, q3_im, X3, Y);
            AE_MULAFC32RA(q4_re, q4_im, X4, Y);
            AE_MULAFC32RA(q5_re, q5_im, X5, Y);
            AE_MULAFC32RA(q6_re, q6_im, X6, Y);
            AE_MULAFC32RA(q7_re, q7_im, X7, Y);
            /* shift input line for the next iteration */
            X0 = X1; X1 = X2; X2 = X3;
            X3 = X4; X4 = X5; X5 = X6; X6 = X7;
        }
        /* save computed samples */
        X0 = AE_ROUND32X2F48SASYM(q0_re, q0_im);
        X1 = AE_ROUND32X2F48SASYM(q1_re, q1_im);
        X2 = AE_ROUND32X2F48SASYM(q2_re, q2_im);
        X3 = AE_ROUND32X2F48SASYM(q3_re, q3_im);
        AE_S32X2_IP(X0, pY, sizeof(complex_fract32));
        AE_S32X2_IP(X1, pY, sizeof(complex_fract32));
        AE_S32X2_IP(X2, pY, sizeof(complex_fract32));
        AE_S32X2_IP(X3, pY, sizeof(complex_fract32));
        X0 = AE_ROUND32X2F48SASYM(q4_re, q4_im);
        X1 = AE_ROUND32X2F48SASYM(q5_re, q5_im);
        X2 = AE_ROUND32X2F48SASYM(q6_re, q6_im);
        X3 = AE_ROUND32X2F48SASYM(q7_re, q7_im);
        AE_S32X2_IP(X0, pY, sizeof(complex_fract32));
        AE_S32X2_IP(X1, pY, sizeof(complex_fract32));
        AE_S32X2_IP(X2, pY, sizeof(complex_fract32));
        AE_S32X2_IP(X3, pY, sizeof(complex_fract32));
    }
    if (N & 4)
    {
        AE_L32X2X2_IP(X0, X1, pX, 2 * sizeof(complex_fract32));
        AE_L32X2X2_IP(X2, X3, pX, 2 * sizeof(complex_fract32));
        AE_S32X2X2_XC(X0, X1, pDw, 2 * sizeof(complex_fract32));
        AE_S32X2X2_XC(X2, X3, pDw, 2 * sizeof(complex_fract32));

        pH = (const ae_int32x2 *)cxfir->coef;
        /* preload data from x */
        pDr = (ae_int32x2 *)pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(complex_fract32));
        AE_L32X2X2_XC(X0, X1, castxcc(ae_int32x4, pDr), 2 * sizeof(complex_fract32));
        AE_L32X2X2_XC(X2, X3, castxcc(ae_int32x4, pDr), 2 * sizeof(complex_fract32));
        /* load data from y */
        AE_L32X2_IP(Y, pH, sizeof(complex_fract32));
        /* compute correlation of 4 values */
        AE_MULFC32RA(q0_re, q0_im, X0, Y);
        AE_MULFC32RA(q1_re, q1_im, X1, Y);
        AE_MULFC32RA(q2_re, q2_im, X2, Y);
        AE_MULFC32RA(q3_re, q3_im, X3, Y);
        /* shift input line for the next iteration */
        X0 = X1; X1 = X2; X2 = X3;

        __Pragma("loop_count min=4, factor=4");
        for (m = 0; m < M; m++)
        {
            /* load data from x */
            AE_L32X2_XC(X3, pDr, sizeof(complex_fract32));
            /* load data from y */
            AE_L32X2_IP(Y, pH, sizeof(complex_fract32));
            /* compute correlation of 4 values */
            AE_MULAFC32RA(q0_re, q0_im, X0, Y);
            AE_MULAFC32RA(q1_re, q1_im, X1, Y);
            AE_MULAFC32RA(q2_re, q2_im, X2, Y);
            AE_MULAFC32RA(q3_re, q3_im, X3, Y);
            /* shift input line for the next iteration */
            X0 = X1; X1 = X2; X2 = X3;
        }
        /* save computed samples */
        X0 = AE_ROUND32X2F48SASYM(q0_re, q0_im);
        X1 = AE_ROUND32X2F48SASYM(q1_re, q1_im);
        X2 = AE_ROUND32X2F48SASYM(q2_re, q2_im);
        X3 = AE_ROUND32X2F48SASYM(q3_re, q3_im);
        AE_S32X2_IP(X0, pY, sizeof(complex_fract32));
        AE_S32X2_IP(X1, pY, sizeof(complex_fract32));
        AE_S32X2_IP(X2, pY, sizeof(complex_fract32));
        AE_S32X2_IP(X3, pY, sizeof(complex_fract32));
    }
#endif

    cxfir->delayPos = (complex_fract32*)pDw;
} // cxfir32x32_process()
