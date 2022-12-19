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
    Real data circular convolution, 16x16-bit
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

void fir_convol16x16( int16_t * restrict r,
                const int16_t * restrict x,
                const int16_t * restrict y,
                int N, int M )
{
    // Circular convolution algorithm:
    //
    //   r[n] = sum( x[mod(n-m,N)]*y[m] )
    //        m=0..M-1
    //
    //   where n = 0..N-1
    //
    const ae_int16x4 *          pX;
    const ae_int16x4 *          pX0;
    const ae_int16x4 *          S0;
    const ae_int16x4 *          S1;
    const ae_int16x4 *          S2;
    const ae_int16x4 *          pY;
          ae_int16x4 * restrict pR;

    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_int16x4 y0, x0, x1, x2;
    ae_int32x2 t0, t1;

    int n, m;

    NASSERT(r);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(r, 8);
    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, 8);
    NASSERT((M > 0) && ((M % 4) == 0));
    NASSERT((N > 0) && ((N % 4) == 0));

    pX0= (const ae_int16x4 *)x;
    pR = (      ae_int16x4 *)r;
    WUR_AE_CBEGIN0((uintptr_t)(x + 0));
    WUR_AE_CEND0  ((uintptr_t)(x + N));

    for (n = 0; n < (N >> 4); n++)
    {
        ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
        ae_int16x4 x3, x4;

        pY = (const ae_int16x4 *)y;
        AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));

        pX = pX0 + 3;
        pX0 += 4;
        AE_L16X4_RIC(x0, pX);
        S0 = pX;
        AE_L16X4_RIC(x1, pX);
        AE_L16X4_RIC(x2, pX);
        S1 = pX;
        AE_L16X4_RIC(x3, pX);
        S2 = pX;
        AE_L16X4_RIC(x4, S2);

        AE_MULFQ16X2_FIR_3(q0, q1, x0, x1, y0);
        AE_MULFQ16X2_FIR_1(q2, q3, x0, x1, y0);
        AE_MULFQ16X2_FIR_3(q4, q5, x1, x2, y0);
        AE_MULFQ16X2_FIR_1(q6, q7, x1, x2, y0);
        AE_MULFQ16X2_FIR_3(q8, q9, x2, x3, y0);
        AE_MULFQ16X2_FIR_1(qa, qb, x2, x3, y0);
        AE_MULFQ16X2_FIR_3(qc, qd, x3, x4, y0);
        AE_MULFQ16X2_FIR_1(qe, qf, x3, x4, y0);

        for (m = 0; m < (M >> 2) - 1; m++)
        {
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));

            AE_L16X4_RIC(x0, S0);
            x1 = AE_L16X4_RI(S0, 0);
            AE_L16X4_RIC(x2, S1);
            x3 = AE_L16X4_RI(S1, 0);
            AE_L16X4_RIC(x4, S2);

            AE_MULAFQ16X2_FIR_3(q0, q1, x0, x1, y0);
            AE_MULAFQ16X2_FIR_1(q2, q3, x0, x1, y0);
            AE_MULAFQ16X2_FIR_3(q4, q5, x1, x2, y0);
            AE_MULAFQ16X2_FIR_1(q6, q7, x1, x2, y0);
            AE_MULAFQ16X2_FIR_3(q8, q9, x2, x3, y0);
            AE_MULAFQ16X2_FIR_1(qa, qb, x2, x3, y0);
            AE_MULAFQ16X2_FIR_3(qc, qd, x3, x4, y0);
            AE_MULAFQ16X2_FIR_1(qe, qf, x3, x4, y0);
        }

        t0 = AE_TRUNCA32X2F64S(qf, qe, 32);
        t1 = AE_SAT32X2(qd, qc);
        AE_S16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), pR, 4 * sizeof(int16_t));
        t0 = AE_TRUNCA32X2F64S(qb, qa, 32);
        t1 = AE_SAT32X2(q9, q8);
        AE_S16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), pR, 4 * sizeof(int16_t));
        t0 = AE_TRUNCA32X2F64S(q7, q6, 32);
        t1 = AE_SAT32X2(q5, q4);
        AE_S16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), pR, 4 * sizeof(int16_t));
        t0 = AE_TRUNCA32X2F64S(q3, q2, 32);
        t1 = AE_SAT32X2(q1, q0);
        AE_S16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), pR, 4 * sizeof(int16_t));
    }

    if (N & 8)
    {
        pY = (const ae_int16x4 *)y;
        pX = pX0 + 1;
        pX0 += 2;
        S0 = pX;
        AE_L16X4_RIC(x0, pX);
        S1 = pX;
        AE_L16X4_RIC(x0, pX);
        q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));

            AE_L16X4_RIC(x0, S0);
            AE_L16X4_RIC(x1, S1);
            x2 = AE_L16X4_RI(S1, 0);

            AE_MULAFQ16X2_FIR_3(q0, q1, x0, x1, y0);
            AE_MULAFQ16X2_FIR_1(q2, q3, x0, x1, y0);
            AE_MULAFQ16X2_FIR_3(q4, q5, x1, x2, y0);
            AE_MULAFQ16X2_FIR_1(q6, q7, x1, x2, y0);
        }

        t0 = AE_TRUNCA32X2F64S(q7, q6, 32);
        t1 = AE_SAT32X2(q5, q4);
        AE_S16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), pR, 4 * sizeof(int16_t));
        t0 = AE_TRUNCA32X2F64S(q3, q2, 32);
        t1 = AE_SAT32X2(q1, q0);
        AE_S16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), pR, 4 * sizeof(int16_t));
    }

    if (N & 4)
    {
        pY = (const ae_int16x4 *)y;
        pX = pX0;
        AE_L16X4_RIC(x0, pX);
        q0 = q1 = q2 = q3 = AE_ZERO64();
        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
            AE_L16X4_RIC(x1, pX);
            AE_MULAFQ16X2_FIR_3(q0, q1, x0, x1, y0);
            AE_MULAFQ16X2_FIR_1(q2, q3, x0, x1, y0);
            x0 = x1;
        }
        t0 = AE_TRUNCA32X2F64S(q3, q2, 32);
        t1 = AE_SAT32X2(q1, q0);
        AE_S16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), pR, 4 * sizeof(int16_t));
    }
}
