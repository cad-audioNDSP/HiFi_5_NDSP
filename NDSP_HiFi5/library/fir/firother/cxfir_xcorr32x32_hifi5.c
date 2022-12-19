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
    complex data circular cross-correlation, complex 32x32
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

/*-------------------------------------------------------------------------
  Circular Correlation
  Estimates the circular cross-correlation between vectors x (of length N) 
  and y (of length M)  resulting in vector r of length N. It is a similar 
  to correlation but x is read in opposite direction.

  Precision: 
  16x16     16x16-bit data, 16-bit outputs
  32x16     32x16-bit data, 32-bit outputs
  32x32     32x32-bit data, 32-bit outputs
  32x32ep   the same as above but using 72-bit accumulator for intermediate 
            computations
  f         floating point 


  Input:
  x[N]      input data Q15, Q31 or floating point
  y[M]      input data Q15, Q31 or floating point
  N         length of x
  M         length of y
  Output:
  r[N]      output data, Q15, Q31 or floating point

  Restriction:
  x,y,r     should not overlap
  x,y,r     aligned on an 16-bytes boundary
  N,M       multiples of 4 and >0
-------------------------------------------------------------------------*/

void cxfir_xcorr32x32( complex_fract32 * restrict r,
                 const complex_fract32 * restrict x,
                 const complex_fract32 * restrict y,
                 int N, int M )
{
    //
    // Circular cross-correlation algorithm:
    //
    //   r[n] = sum( x[mod(n+m,N)]*y[m] )
    //        m=0..M-1
    //
    //   where n = 0..N-1
    //
    const ae_int32x2 * restrict pX;
    const ae_int64   * restrict pY;
          ae_int32x2 * restrict pR;

    ae_f64 q0_re, q1_re, q2_re, q3_re;
    ae_f64 q0_im, q1_im, q2_im, q3_im;
    ae_int32x2 X0, X1, X2, X3, Y;
    ae_int64 Y_;

    int n, m;

    NASSERT(r);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(r, 16);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT(M > 0 && M % 4 == 0);
    NASSERT(N > 0 && N % 4 == 0);

    pR = (ae_int32x2*)r;
    /* set circular buffer boundaries */
    WUR_AE_CBEGIN0((uintptr_t)(x + 0));
    WUR_AE_CEND0  ((uintptr_t)(x + N));

    for (n = 0; n < (N >> 3); n++)
    {
        ae_f64 q4_re, q5_re, q6_re, q7_re;
        ae_f64 q4_im, q5_im, q6_im, q7_im;
        ae_int32x2 X4, X5, X6, X7;

        pX = (const ae_int32x2 *)(x + 8 * n);
        pY = (const ae_int64   *)y;
        q0_re = q1_re = q2_re = q3_re = AE_ZERO64();
        q0_im = q1_im = q2_im = q3_im = AE_ZERO64();
        q4_re = q5_re = q6_re = q7_re = AE_ZERO64();
        q4_im = q5_im = q6_im = q7_im = AE_ZERO64();
        /* preload data from x */
        AE_L32X2_XC(X0, pX, sizeof(complex_fract32));
        AE_L32X2_XC(X1, pX, sizeof(complex_fract32));
        AE_L32X2_XC(X2, pX, sizeof(complex_fract32));
        AE_L32X2_XC(X3, pX, sizeof(complex_fract32));
        AE_L32X2_XC(X4, pX, sizeof(complex_fract32));
        AE_L32X2_XC(X5, pX, sizeof(complex_fract32));
        AE_L32X2_XC(X6, pX, sizeof(complex_fract32));

        __Pragma("loop_count min=4, factor=4");
        for (m = 0; m < M; m++)
        {
            /* load data from x */
            AE_L32X2_XC(X7, pX, sizeof(complex_fract32));
            /* load data from y */
            AE_L64_IP(Y_, pY, sizeof(complex_fract32));
            Y = AE_MOVINT32X2_FROMINT64(Y_);
            /* compute correlation of 8 values */
            AE_MULAFC32RA(q0_im, q0_re, X0, Y);
            AE_MULAFC32RA(q1_im, q1_re, X1, Y);
            AE_MULAFC32RA(q2_im, q2_re, X2, Y);
            AE_MULAFC32RA(q3_im, q3_re, X3, Y);
            AE_MULAFC32RA(q4_im, q4_re, X4, Y);
            AE_MULAFC32RA(q5_im, q5_re, X5, Y);
            AE_MULAFC32RA(q6_im, q6_re, X6, Y);
            AE_MULAFC32RA(q7_im, q7_re, X7, Y);
            /* shift input line for the next iteration */
            X0 = X1; X1 = X2; X2 = X3;
            X3 = X4; X4 = X5; X5 = X6; X6 = X7;
        }
        /* save computed samples */
        X0 = AE_ROUND32X2F48SASYM(q0_re, q0_im);
        X1 = AE_ROUND32X2F48SASYM(q1_re, q1_im);
        X2 = AE_ROUND32X2F48SASYM(q2_re, q2_im);
        X3 = AE_ROUND32X2F48SASYM(q3_re, q3_im);
        AE_S32X2_IP(X0, pR, sizeof(complex_fract32));
        AE_S32X2_IP(X1, pR, sizeof(complex_fract32));
        AE_S32X2_IP(X2, pR, sizeof(complex_fract32));
        AE_S32X2_IP(X3, pR, sizeof(complex_fract32));
        X0 = AE_ROUND32X2F48SASYM(q4_re, q4_im);
        X1 = AE_ROUND32X2F48SASYM(q5_re, q5_im);
        X2 = AE_ROUND32X2F48SASYM(q6_re, q6_im);
        X3 = AE_ROUND32X2F48SASYM(q7_re, q7_im);
        AE_S32X2_IP(X0, pR, sizeof(complex_fract32));
        AE_S32X2_IP(X1, pR, sizeof(complex_fract32));
        AE_S32X2_IP(X2, pR, sizeof(complex_fract32));
        AE_S32X2_IP(X3, pR, sizeof(complex_fract32));
    }
    if (N & 4)
    {
        pX = (const ae_int32x2 *)(x + 8 * n);
        pY = (const ae_int64   *)y;
        q0_re = q1_re = q2_re = q3_re = AE_ZERO64();
        q0_im = q1_im = q2_im = q3_im = AE_ZERO64();
        /* preload data from x */
        AE_L32X2_XC(X0, pX, sizeof(complex_fract32));
        AE_L32X2_XC(X1, pX, sizeof(complex_fract32));
        AE_L32X2_XC(X2, pX, sizeof(complex_fract32));

        __Pragma("loop_count min=4, factor=4");
        for (m = 0; m < M; m++)
        {
            /* load data from x */
            AE_L32X2_XC(X3, pX, sizeof(complex_fract32));
            /* load data from y */
            AE_L64_IP(Y_, pY, sizeof(complex_fract32));
            Y = AE_MOVINT32X2_FROMINT64(Y_);
            /* compute correlation of 4 values */
            AE_MULAFC32RA(q0_im, q0_re, X0, Y);
            AE_MULAFC32RA(q1_im, q1_re, X1, Y);
            AE_MULAFC32RA(q2_im, q2_re, X2, Y);
            AE_MULAFC32RA(q3_im, q3_re, X3, Y);
            /* shift input line for the next iteration */
            X0 = X1; X1 = X2; X2 = X3;
        }
        /* save computed samples */
        X0 = AE_ROUND32X2F48SASYM(q0_re, q0_im);
        X1 = AE_ROUND32X2F48SASYM(q1_re, q1_im);
        X2 = AE_ROUND32X2F48SASYM(q2_re, q2_im);
        X3 = AE_ROUND32X2F48SASYM(q3_re, q3_im);
        AE_S32X2_IP(X0, pR, sizeof(complex_fract32));
        AE_S32X2_IP(X1, pR, sizeof(complex_fract32));
        AE_S32X2_IP(X2, pR, sizeof(complex_fract32));
        AE_S32X2_IP(X3, pR, sizeof(complex_fract32));
    }
} // cxfir_xcorr32x32()
