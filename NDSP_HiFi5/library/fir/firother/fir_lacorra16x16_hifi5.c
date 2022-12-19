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
    Real data linear auto-correlation, 16x16-bit, no requirements on vectors
    length and alignment.
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

/*-------------------------------------------------------------------------
  Linear Autocorrelation 
  Functions estimate the linear auto-correlation of vector x. Returns 
  autocorrelation of length N.

  Precision: 
  16x16   16-bit data, 16-bit outputs
  32x32   32-bit data, 32-bit outputs

  Input:
  s[]       scratch area of
            FIR_LACORRA16X16_SCRATCH_SIZE( N )
            FIR_LACORRA32X32_SCRATCH_SIZE( N ) bytes
            
  x[N]      input data Q15, Q31 
  N         length of x
  Output:
  r[N]      output data, Q15, Q31

  Restrictions:
  x,r,s   should not overlap
  N       >0
  s       aligned on an 16-bytes boundary
-------------------------------------------------------------------------*/

void fir_lacorra16x16 ( void         * restrict s,
                      int16_t        * restrict r,
                      const int16_t  * restrict x, int N )
{
    void       * s_ptr;
    int16_t    * x_buf;
    ae_int16x4 * restrict D;
    ae_valign S_va;

    const ae_int16x4 *          pX;
    const ae_int16x4 *          S0;
    const ae_int16x4 *          S1;
    const ae_int16x4 *          pY;
          ae_int16x4 * restrict pR;

    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_int16x4 y0, x0, x1, x2;
    ae_int32x2 t0, t1;
    ae_valign ar;

    int n, m, m_iter;

    NASSERT(s);
    NASSERT(r);
    NASSERT(x);
    NASSERT_ALIGN(s, 8);
    NASSERT(N>0);

    //----------------------------------------------------------------------------
    // Partition the scratch memory area.
    s_ptr = s;
    x_buf = (int16_t*)ALIGNED_ADDR(s_ptr, 8);
    s_ptr = x_buf + ((N + 4 + 3)&~3);
    ASSERT((int8_t *)s_ptr - (int8_t *)s <= (int)FIR_LACORRA16X16_SCRATCH_SIZE(N));

    //----------------------------------------------------------------------------
    // Copy x[N].
    S0 = (const ae_int16x4 *)x;
    D  = (      ae_int16x4 *)x_buf;
    S_va = AE_LA64_PP(S0);
    for (n = 0; n < (N >> 2); n++)
    {
        AE_LA16X4_IP(x0, S_va, S0);
        AE_S16X4_IP(x0, D, 4 * sizeof(int16_t));
    }
    AE_S16X4_I(AE_ZERO16(), D, 0);
    if (N & 3)
    {
        AE_S16X4_I(AE_ZERO16(), D, 4 * sizeof(int16_t));
        for (n = (N&~3); n < N; n++)
        {
            AE_L16_IP(x0, castxcc(ae_int16, S0), sizeof(int16_t));
            AE_S16_0_IP(x0, castxcc(ae_int16, D), sizeof(int16_t));
        }
    }

    m_iter = ((N + 3) >> 2) - 1;

    pX = (const ae_int16x4 *)x_buf;
    pR = (      ae_int16x4 *)r;
    ar = AE_ZALIGN64();

    for (n = 0; n < (N >> 4); n++, m_iter -= 4)
    {
        ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
        ae_int16x4 x3, x4;

        pY = (const ae_int16x4 *)x_buf;
        AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));

        AE_L16X4_IP(x0, pX, 4 * sizeof(int16_t));
        S0 = pX;
        AE_L16X4_IP(x1, pX, 4 * sizeof(int16_t));
        AE_L16X4_IP(x2, pX, 4 * sizeof(int16_t));
        S1 = pX;
        AE_L16X4_IP(x3, pX, 4 * sizeof(int16_t));
        x4 = AE_L16X4_I(pX, 0);

        AE_MULFQ16X2_FIR_3(q0, q1, x0, x1, y0);
        AE_MULFQ16X2_FIR_1(q2, q3, x0, x1, y0);
        AE_MULFQ16X2_FIR_3(q4, q5, x1, x2, y0);
        AE_MULFQ16X2_FIR_1(q6, q7, x1, x2, y0);
        AE_MULFQ16X2_FIR_3(q8, q9, x2, x3, y0);
        AE_MULFQ16X2_FIR_1(qa, qb, x2, x3, y0);
        AE_MULFQ16X2_FIR_3(qc, qd, x3, x4, y0);
        AE_MULFQ16X2_FIR_1(qe, qf, x3, x4, y0);

        for (m = 0; m < m_iter - 3; m++)
        {
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));

            AE_L16X4_IP(x0, S0, 4 * sizeof(int16_t));
            x1 = AE_L16X4_I(S0, 0);
            AE_L16X4_IP(x2, S1, 4 * sizeof(int16_t));
            x3 = AE_L16X4_I(S1, 0);
            x4 = AE_L16X4_I(S1, 4 * sizeof(int16_t));

            AE_MULAFQ16X2_FIR_3(q0, q1, x0, x1, y0);
            AE_MULAFQ16X2_FIR_1(q2, q3, x0, x1, y0);
            AE_MULAFQ16X2_FIR_3(q4, q5, x1, x2, y0);
            AE_MULAFQ16X2_FIR_1(q6, q7, x1, x2, y0);
            AE_MULAFQ16X2_FIR_3(q8, q9, x2, x3, y0);
            AE_MULAFQ16X2_FIR_1(qa, qb, x2, x3, y0);
            AE_MULAFQ16X2_FIR_3(qc, qd, x3, x4, y0);
            AE_MULAFQ16X2_FIR_1(qe, qf, x3, x4, y0);
        }

        AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
        AE_MULAFQ16X2_FIR_3(q0, q1, x1, x2, y0);
        AE_MULAFQ16X2_FIR_1(q2, q3, x1, x2, y0);
        AE_MULAFQ16X2_FIR_3(q4, q5, x2, x3, y0);
        AE_MULAFQ16X2_FIR_1(q6, q7, x2, x3, y0);
        AE_MULAFQ16X2_FIR_3(q8, q9, x3, x4, y0);
        AE_MULAFQ16X2_FIR_1(qa, qb, x3, x4, y0);

        AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
        AE_MULAFQ16X2_FIR_3(q0, q1, x2, x3, y0);
        AE_MULAFQ16X2_FIR_1(q2, q3, x2, x3, y0);
        AE_MULAFQ16X2_FIR_3(q4, q5, x3, x4, y0);
        AE_MULAFQ16X2_FIR_1(q6, q7, x3, x4, y0);

        AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
        AE_MULAFQ16X2_FIR_3(q0, q1, x3, x4, y0);
        AE_MULAFQ16X2_FIR_1(q2, q3, x3, x4, y0);

        t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
        t1 = AE_SAT32X2(q2, q3);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
        t0 = AE_TRUNCA32X2F64S(q4, q5, 32);
        t1 = AE_SAT32X2(q6, q7);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
        t0 = AE_TRUNCA32X2F64S(q8, q9, 32);
        t1 = AE_SAT32X2(qa, qb);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
        t0 = AE_TRUNCA32X2F64S(qc, qd, 32);
        t1 = AE_SAT32X2(qe, qf);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
    }

    if (N & 8)
    {
        pY = (const ae_int16x4 *)x_buf;
        AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));

        AE_L16X4_IP(x0, pX, 4 * sizeof(int16_t));
        S0 = pX;
        AE_L16X4_IP(x1, pX, 4 * sizeof(int16_t));
        S1 = pX;
        x2 = AE_L16X4_I(S1, 0);

        AE_MULFQ16X2_FIR_3(q0, q1, x0, x1, y0);
        AE_MULFQ16X2_FIR_1(q2, q3, x0, x1, y0);
        AE_MULFQ16X2_FIR_3(q4, q5, x1, x2, y0);
        AE_MULFQ16X2_FIR_1(q6, q7, x1, x2, y0);

        for (m = 0; m < m_iter - 1; m++)
        {
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));

            AE_L16X4_XP(x0, S0, 4 * sizeof(int16_t));
            AE_L16X4_XP(x1, S1, 4 * sizeof(int16_t));
            x2 = AE_L16X4_I(S1, 0);

            AE_MULAFQ16X2_FIR_3(q0, q1, x0, x1, y0);
            AE_MULAFQ16X2_FIR_1(q2, q3, x0, x1, y0);
            AE_MULAFQ16X2_FIR_3(q4, q5, x1, x2, y0);
            AE_MULAFQ16X2_FIR_1(q6, q7, x1, x2, y0);
        }

        AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
        AE_MULAFQ16X2_FIR_3(q0, q1, x1, x2, y0);
        AE_MULAFQ16X2_FIR_1(q2, q3, x1, x2, y0);

        t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
        t1 = AE_SAT32X2(q2, q3);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
        t0 = AE_TRUNCA32X2F64S(q4, q5, 32);
        t1 = AE_SAT32X2(q6, q7);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
        m_iter -= 2;
    }

    if (N & 4)
    {
        pY = (const ae_int16x4 *)x_buf;
        AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
        AE_L16X4_IP(x0, pX, 4 * sizeof(int16_t));
        x1 = AE_L16X4_I(pX, 0);
        AE_MULFQ16X2_FIR_3(q0, q1, x0, x1, y0);
        AE_MULFQ16X2_FIR_1(q2, q3, x0, x1, y0);
        if (m_iter)
        {
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
            x2 = AE_L16X4_I(pX, 4 * sizeof(int16_t));
            AE_MULAFQ16X2_FIR_3(q0, q1, x1, x2, y0);
            AE_MULAFQ16X2_FIR_1(q2, q3, x1, x2, y0);
        }
        t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
        t1 = AE_SAT32X2(q2, q3);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
    }
    AE_SA64POS_FP(ar, pR);

    N &= 3;
    if (N)
    {
        pY = (const ae_int16x4 *)x_buf;
        AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
        AE_L16X4_IP(x0, pX, 4 * sizeof(int16_t));
        AE_L16X4_IP(x1, pX, 4 * sizeof(int16_t));
        AE_MULFQ16X2_FIR_3(q0, q1, x0, x1, y0);
        AE_MULFQ16X2_FIR_1(q2, q3, x0, x1, y0);
        t0 = AE_TRUNCA32X2F64S(q1, q0, 32);
        t1 = AE_SAT32X2(q2, q2);
        x0 = AE_ROUND16X4F32SASYM(t1, t0);
        AE_S16_0_IP(x0, castxcc(ae_int16, pR), sizeof(int16_t));
        if (N > 1) { x0 = AE_SEL16_4321(x0, x0); AE_S16_0_IP(x0, castxcc(ae_int16, pR), sizeof(int16_t)); }
        if (N > 2) { x0 = AE_SEL16_4321(x0, x0); AE_S16_0_IP(x0, castxcc(ae_int16, pR), sizeof(int16_t)); }
    }
} /* fir_lacorra16x16() */
