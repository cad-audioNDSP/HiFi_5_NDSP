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
    Real data circular convolution, 32x32-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

#define AE_L32X2X2_RXC(x0, x1, pX, ofs)                    \
{                                                          \
    ae_int64 tmp0, tmp1;                                   \
    AE_L64X2_XC(tmp1, tmp0, castxcc(ae_int64x2, pX), ofs); \
    x0 = AE_MOVINT32X2_FROMINT64(tmp0);                    \
    x1 = AE_MOVINT32X2_FROMINT64(tmp1);                    \
}

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

void fir_convol32x32( int32_t * restrict r,
                const int32_t * restrict x,
                const int32_t * restrict y,
                int N, int M )
{
    // Circular convolution algorithm:
    //
    //   r[n] = sum( x[mod(n-m,N)]*y[m] )
    //        m=0..M-1
    //
    //   where n = 0..N-1
    const ae_int32x4 *          pX;
    const ae_int32x4 *          S0;
    const ae_int32x4 *          S1;
    const ae_int32x4 *          S2;
    const ae_int32x2 *          pY;
          ae_int32x2 * restrict pR;

    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_int32x2 x0, x1, x2, x3, x4, x5;
    ae_int32x2 y0, y1;

    int n, m;

    NASSERT(r);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(r, 16);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT(M > 0 && M % 4 == 0);
    NASSERT(N > 0 && N % 4 == 0);

    pX = (const ae_int32x4 *)x;
    pR = (      ae_int32x2 *)r;
    WUR_AE_CBEGIN0((uintptr_t)(x + 0));
    WUR_AE_CEND0  ((uintptr_t)(x + N));
    AE_ADDCIRC32X2_XC(castxcc(ae_int32x2, pX), -4 * (int)sizeof(int32_t));

    for (n = 0; n < (N >> 3); n++)
    {
        pY = (const ae_int32x2 *)y;
        S2 = pX;
        AE_ADDCIRC32X2_XC(castxcc(ae_int32x2, pX), 4 * sizeof(int32_t));
        S1 = pX;
        AE_ADDCIRC32X2_XC(castxcc(ae_int32x2, pX), 4 * sizeof(int32_t));
        S0 = pX;

        q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            AE_L32X2_IP(y1, pY, 2 * sizeof(int32_t));
            AE_L32X2X2_RXC(x0, x1, S0, -4 * (int)sizeof(int32_t));
            AE_L32X2X2_RXC(x2, x3, S1, -4 * (int)sizeof(int32_t));
            AE_L32X2X2_RXC(x4, x5, S2, -4 * (int)sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_H(q7, q6, x0, x1, y0);
            AE_MULAFD32X2RA_FIR_H(q5, q4, x1, x2, y0);
            AE_MULAFD32X2RA_FIR_H(q3, q2, x2, x3, y0);
            AE_MULAFD32X2RA_FIR_H(q1, q0, x3, x4, y0);
            AE_MULAFD32X2RA_FIR_H(q7, q6, x1, x2, y1);
            AE_MULAFD32X2RA_FIR_H(q5, q4, x2, x3, y1);
            AE_MULAFD32X2RA_FIR_H(q3, q2, x3, x4, y1);
            AE_MULAFD32X2RA_FIR_H(q1, q0, x4, x5, y1);
        }

        x0 = AE_ROUND32X2F48SASYM(q0, q1);
        x1 = AE_ROUND32X2F48SASYM(q2, q3);
        x2 = AE_ROUND32X2F48SASYM(q4, q5);
        x3 = AE_ROUND32X2F48SASYM(q6, q7);
        AE_S32X2_IP(x0, pR, 2 * sizeof(int32_t));
        AE_S32X2_IP(x1, pR, 2 * sizeof(int32_t));
        AE_S32X2_IP(x2, pR, 2 * sizeof(int32_t));
        AE_S32X2_IP(x3, pR, 2 * sizeof(int32_t));
    }
    if (N & 4)
    {
        pY = (const ae_int32x2 *)y;
        S1 = pX;
        AE_ADDCIRC32X2_XC(castxcc(ae_int32x2, pX), 4 * sizeof(int32_t));
        S0 = pX;
        q0 = q1 = q2 = q3 = AE_ZERO64();

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            AE_L32X2_IP(y1, pY, 2 * sizeof(int32_t));
            AE_L32X2X2_RXC(x0, x1, S0, -4 * (int)sizeof(int32_t));
            AE_L32X2X2_RXC(x2, x3, S1, -4 * (int)sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_H(q3, q2, x0, x1, y0);
            AE_MULAFD32X2RA_FIR_H(q1, q0, x1, x2, y0);
            AE_MULAFD32X2RA_FIR_H(q3, q2, x1, x2, y1);
            AE_MULAFD32X2RA_FIR_H(q1, q0, x2, x3, y1);
        }

        x0 = AE_ROUND32X2F48SASYM(q0, q1);
        x1 = AE_ROUND32X2F48SASYM(q2, q3);
        AE_S32X2_IP(x0, pR, 2 * sizeof(int32_t));
        AE_S32X2_IP(x1, pR, 2 * sizeof(int32_t));
    }
} // fir_convol32x32()
