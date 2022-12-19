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
    Complex data circular convolution, 32x16-bit
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
  Circular Convolution
  Performs circular convolution between vectors x (of length N) and y (of 
  length M)  resulting in vector r of length N.

  Precision: 
  16x16     16x16-bit data, 16-bit outputs
  32x16     32x16-bit data, 32-bit outputs 
  32x32     32x32-bit data, 32-bit outputs
  32x32ep   the same as above but using 72-bit accumulator for intermediate 
            computations
  f         floating point

  Input:
  x[N]      input data, Q15, Q31 or floating point
  y[M]      input data, Q15, Q31 or floating point
  N         length of x
  M         length of y
  Output:
  r[N]      output data, Q15, Q31 or floating point

  Restriction:
  x,y,r     should not overlap
  x,y,r     aligned on an 16-bytes boundary
  N,M       multiples of 4 and >0
-------------------------------------------------------------------------*/

void cxfir_convol32x16(complex_fract32 * restrict r,
                 const complex_fract32 * restrict x,
                 const complex_fract16 * restrict y,
                 int N, int M)
{
    int n, m;
    const int32_t    *          xn;
    const ae_int32x2 *          X;
    const ae_int32x2 *          Y;
          ae_int32x2 * restrict R;

    ae_f64 q0r, q0i, q1r, q1i, q2r, q2i, q3r, q3i;
    ae_int32x2 x0, x1, x2, x3, x4, x5;
    ae_int32x2 x0_conj, x1_conj, x2_conj, x3_conj;
    ae_int32x2 y_tmp;
    ae_int16x4 y0, y1;

    NASSERT(r);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(r, 16);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT(M % 4 == 0 && N % 4 == 0);

    xn = (const int32_t    *)(x + 3);
    R  = (      ae_int32x2 *)r;

    WUR_AE_CBEGIN0((uintptr_t)(x + 0));
    WUR_AE_CEND0  ((uintptr_t)(x + N));

#if 0
    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 2); n++)
    {
        Y = (const ae_int32x2 *)y;
        X = (const ae_int32x2 *)xn;
        xn += 8;

        AE_L32X2_XC(x0, X, -(int)sizeof(complex_fract32));
        AE_L32X2_XC(x1, X, -(int)sizeof(complex_fract32));
        AE_L32X2_XC(x2, X, -(int)sizeof(complex_fract32));
        AE_L32X2_XC(x3, X, -(int)sizeof(complex_fract32));

        q0r = q1r = q2r = q3r = AE_ZERO64();
        q0i = q1i = q2i = q3i = AE_ZERO64();

        __Pragma("loop_count min=2, factor=2");
        for (m = 0; m < (M >> 1); m++)
        {
            AE_L16X4_IP(y0, castxcc(ae_int16x4, Y), 2 * sizeof(complex_fract16));

            AE_L32X2_XC(x4, X, -(int)sizeof(complex_fract32));
            AE_L32X2_XC(x5, X, -(int)sizeof(complex_fract32));

            AE_MULAFC32X16W_H(q3r, q3i, x0, y0);
            AE_MULAFC32X16W_L(q3r, q3i, x1, y0);
            AE_MULAFC32X16W_H(q2r, q2i, x1, y0);
            AE_MULAFC32X16W_L(q2r, q2i, x2, y0);
            AE_MULAFC32X16W_H(q1r, q1i, x2, y0);
            AE_MULAFC32X16W_L(q1r, q1i, x3, y0);
            AE_MULAFC32X16W_H(q0r, q0i, x3, y0);
            AE_MULAFC32X16W_L(q0r, q0i, x4, y0);

            x0 = x2; x1 = x3; x2 = x4; x3 = x5;
        }

        AE_S32X2RA64S_IP(q0r, q0i, R);
        AE_S32X2RA64S_IP(q1r, q1i, R);
        AE_S32X2RA64S_IP(q2r, q2i, R);
        AE_S32X2RA64S_IP(q3r, q3i, R);
    }
#else
    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 2); n++)
    {
        Y = (const ae_int32x2 *)y;
        X = (const ae_int32x2 *)xn;
        xn += 8;

        AE_L32X2_XC(x0, X, -(int)sizeof(complex_fract32));
        AE_L32X2_XC(x1, X, -(int)sizeof(complex_fract32));
        AE_L32X2_XC(x2, X, -(int)sizeof(complex_fract32));
        AE_L32X2_XC(x3, X, -(int)sizeof(complex_fract32));
        x0_conj = AE_MUL32JS(x0);
        x1_conj = AE_MUL32JS(x1);
        x2_conj = AE_MUL32JS(x2);
        x3_conj = AE_MUL32JS(x3);

        AE_L32X2_IP(y_tmp, Y, 2 * sizeof(complex_fract16));
        y0 = AE_MOVF16X4_FROMF32X2(y_tmp);

        AE_MULZAAAA2Q32X16(q1r, q3r, x0_conj, x1_conj, AE_ZERO16(), y0);
        AE_MULZAAAA2Q32X16(q1i, q3i, x0, x1, AE_ZERO16(), y0);
        AE_MULZAAAA2Q32X16(q0r, q2r, x1_conj, x2_conj, AE_ZERO16(), y0);
        AE_MULZAAAA2Q32X16(q0i, q2i, x1, x2, AE_ZERO16(), y0);

        AE_L32X2_XC(x4, X, -(int)sizeof(complex_fract32));
        AE_L32X2_XC(x5, X, -(int)sizeof(complex_fract32));

        AE_L32X2_XP(y_tmp, Y, 2 * sizeof(complex_fract16));
        y1 = AE_MOVF16X4_FROMF32X2(y_tmp);

        AE_MULAAAA2Q32X16(q1r, q3r, x2_conj, x3_conj, y0, y1);
        AE_MULAAAA2Q32X16(q1i, q3i, x2, x3, y0, y1);
        AE_MULAAAA2Q32X16(q0r, q2r, x3_conj, AE_MUL32JS(x4), y0, y1);
        AE_MULAAAA2Q32X16(q0i, q2i, x3, x4, y0, y1);

        x2_conj = AE_MUL32JS(x4);
        x3_conj = AE_MUL32JS(x5);
        x2 = x4; x3 = x5;

        y0 = y1;

        __Pragma("loop_count factor=2");
        for (m = 0; m < (M >> 1) - 2; m++)
        {
            AE_L32X2_XC(x4, X, -(int)sizeof(complex_fract32));
            AE_L32X2_XC(x5, X, -(int)sizeof(complex_fract32));

            AE_L32X2_XP(y_tmp, Y, 2 * sizeof(complex_fract16));
            y1 = AE_MOVF16X4_FROMF32X2(y_tmp);

            AE_MULAAAA2Q32X16(q1r, q3r, x2_conj, x3_conj, y0, y1);
            AE_MULAAAA2Q32X16(q1i, q3i, x2, x3, y0, y1);
            AE_MULAAAA2Q32X16(q0r, q2r, x3_conj, AE_MUL32JS(x4), y0, y1);
            AE_MULAAAA2Q32X16(q0i, q2i, x3, x4, y0, y1);

            x2_conj = AE_MUL32JS(x4);
            x3_conj = AE_MUL32JS(x5);
            x2 = x4; x3 = x5;

            y0 = y1;
        }

        AE_L32X2_XC(x4, X, -(int)sizeof(complex_fract32));
        AE_MULAAAAQ32X16(q1r, x2_conj, x3_conj, y0);
        AE_MULAAAAQ32X16(q1i, x2, x3, y0);
        AE_MULAAAAQ32X16(q0r, x3_conj, AE_MUL32JS(x4), y0);
        AE_MULAAAAQ32X16(q0i, x3, x4, y0);

        q0r = AE_SLAI64(q0r, 1); q0i = AE_SLAI64(q0i, 1);
        q1r = AE_SLAI64(q1r, 1); q1i = AE_SLAI64(q1i, 1);
        q2r = AE_SLAI64(q2r, 1); q2i = AE_SLAI64(q2i, 1);
        q3r = AE_SLAI64(q3r, 1); q3i = AE_SLAI64(q3i, 1);

        AE_S32X2RA64S_IP(q0r, q0i, R);
        AE_S32X2RA64S_IP(q1r, q1i, R);
        AE_S32X2RA64S_IP(q2r, q2i, R);
        AE_S32X2RA64S_IP(q3r, q3i, R);
    }
#endif
} /* cxfir_convol32x16() */
