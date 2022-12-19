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
    helper for correlation/convolution
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "raw_lxcorr16x16.h"

/*-----------------------------------------------------
raw linear correlation:
Input:
x[N] padded with extra 6 zeroes
y[M] padded with extra 4 zeroes
Output:
r[N+M-1]
restrictions:
M should be >0
x,y aligned on 16-byte boundary
-----------------------------------------------------*/
void raw_lxcorr16x16(int16_t * r, const int16_t * restrict x, const int16_t * restrict y, int N, int M)
{
    const ae_int16x4 *          pX;
    const ae_int16x4 *          S0;
    const ae_int16x4 *          S1;
    const ae_int16x4 *          pY;
          ae_int16x4 * restrict pR;

    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_int16x4 y0, x0, x1, x2;
    ae_int32x2 t0, t1;
    ae_valign ar;

    int n, m;

    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, 8);
    NASSERT(M > 0);
    NASSERT(N > 0);

    /*
     * Compute first M-1 entries.
     */
    {
        int n_iter = (M - 1);
        int m_iter = ((n_iter + 3) >> 2) - 1;

        pX = (const ae_int16x4 *)y;
        pR = (ae_int16x4 *)(r + M - 1 - 1);
        ar = AE_ZALIGN64();

        for (n = 0; n < (n_iter >> 4); n++, m_iter -= 4)
        {
            ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
            ae_int16x4 x3, x4;

            pY = (const ae_int16x4 *)x;
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));

            AE_L16X4_IP(x0, pX, 4 * sizeof(int16_t));
            S0 = pX;
            AE_L16X4_IP(x1, pX, 4 * sizeof(int16_t));
            AE_L16X4_IP(x2, pX, 4 * sizeof(int16_t));
            S1 = pX;
            AE_L16X4_IP(x3, pX, 4 * sizeof(int16_t));
            x4 = AE_L16X4_I(pX, 0);

            AE_MULFQ16X2_FIR_2(q0, q1, x0, x1, y0);
            AE_MULFQ16X2_FIR_0(q2, q3, x0, x1, y0);
            AE_MULFQ16X2_FIR_2(q4, q5, x1, x2, y0);
            AE_MULFQ16X2_FIR_0(q6, q7, x1, x2, y0);
            AE_MULFQ16X2_FIR_2(q8, q9, x2, x3, y0);
            AE_MULFQ16X2_FIR_0(qa, qb, x2, x3, y0);
            AE_MULFQ16X2_FIR_2(qc, qd, x3, x4, y0);
            AE_MULFQ16X2_FIR_0(qe, qf, x3, x4, y0);

            for (m = 0; m < m_iter - 3; m++)
            {
                AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));

                AE_L16X4_IP(x0, S0, 4 * sizeof(int16_t));
                x1 = AE_L16X4_I(S0, 0);
                AE_L16X4_IP(x2, S1, 4 * sizeof(int16_t));
                x3 = AE_L16X4_I(S1, 0);
                x4 = AE_L16X4_I(S1, 4 * sizeof(int16_t));

                AE_MULAFQ16X2_FIR_2(q0, q1, x0, x1, y0);
                AE_MULAFQ16X2_FIR_0(q2, q3, x0, x1, y0);
                AE_MULAFQ16X2_FIR_2(q4, q5, x1, x2, y0);
                AE_MULAFQ16X2_FIR_0(q6, q7, x1, x2, y0);
                AE_MULAFQ16X2_FIR_2(q8, q9, x2, x3, y0);
                AE_MULAFQ16X2_FIR_0(qa, qb, x2, x3, y0);
                AE_MULAFQ16X2_FIR_2(qc, qd, x3, x4, y0);
                AE_MULAFQ16X2_FIR_0(qe, qf, x3, x4, y0);
            }

            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
            AE_MULAFQ16X2_FIR_2(q0, q1, x1, x2, y0);
            AE_MULAFQ16X2_FIR_0(q2, q3, x1, x2, y0);
            AE_MULAFQ16X2_FIR_2(q4, q5, x2, x3, y0);
            AE_MULAFQ16X2_FIR_0(q6, q7, x2, x3, y0);
            AE_MULAFQ16X2_FIR_2(q8, q9, x3, x4, y0);
            AE_MULAFQ16X2_FIR_0(qa, qb, x3, x4, y0);

            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
            AE_MULAFQ16X2_FIR_2(q0, q1, x2, x3, y0);
            AE_MULAFQ16X2_FIR_0(q2, q3, x2, x3, y0);
            AE_MULAFQ16X2_FIR_2(q4, q5, x3, x4, y0);
            AE_MULAFQ16X2_FIR_0(q6, q7, x3, x4, y0);

            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
            AE_MULAFQ16X2_FIR_2(q0, q1, x3, x4, y0);
            AE_MULAFQ16X2_FIR_0(q2, q3, x3, x4, y0);

            t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
            t1 = AE_SAT32X2(q2, q3);
            AE_SA16X4_RIP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
            t0 = AE_TRUNCA32X2F64S(q4, q5, 32);
            t1 = AE_SAT32X2(q6, q7);
            AE_SA16X4_RIP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
            t0 = AE_TRUNCA32X2F64S(q8, q9, 32);
            t1 = AE_SAT32X2(qa, qb);
            AE_SA16X4_RIP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
            t0 = AE_TRUNCA32X2F64S(qc, qd, 32);
            t1 = AE_SAT32X2(qe, qf);
            AE_SA16X4_RIP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
        }

        if (n_iter & 8)
        {
            pY = (const ae_int16x4 *)x;
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));

            AE_L16X4_IP(x0, pX, 4 * sizeof(int16_t));
            S0 = pX;
            AE_L16X4_IP(x1, pX, 4 * sizeof(int16_t));
            S1 = pX;
            x2 = AE_L16X4_I(S1, 0);

            AE_MULFQ16X2_FIR_2(q0, q1, x0, x1, y0);
            AE_MULFQ16X2_FIR_0(q2, q3, x0, x1, y0);
            AE_MULFQ16X2_FIR_2(q4, q5, x1, x2, y0);
            AE_MULFQ16X2_FIR_0(q6, q7, x1, x2, y0);

            for (m = 0; m < m_iter - 1; m++)
            {
                AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));

                AE_L16X4_XP(x0, S0, 4 * sizeof(int16_t));
                AE_L16X4_XP(x1, S1, 4 * sizeof(int16_t));
                x2 = AE_L16X4_I(S1, 0);

                AE_MULAFQ16X2_FIR_2(q0, q1, x0, x1, y0);
                AE_MULAFQ16X2_FIR_0(q2, q3, x0, x1, y0);
                AE_MULAFQ16X2_FIR_2(q4, q5, x1, x2, y0);
                AE_MULAFQ16X2_FIR_0(q6, q7, x1, x2, y0);
            }

            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
            AE_MULAFQ16X2_FIR_2(q0, q1, x1, x2, y0);
            AE_MULAFQ16X2_FIR_0(q2, q3, x1, x2, y0);

            t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
            t1 = AE_SAT32X2(q2, q3);
            AE_SA16X4_RIP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
            t0 = AE_TRUNCA32X2F64S(q4, q5, 32);
            t1 = AE_SAT32X2(q6, q7);
            AE_SA16X4_RIP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
            m_iter -= 2;
        }

        if (n_iter & 4)
        {
            pY = (const ae_int16x4 *)x;
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
            AE_L16X4_IP(x0, pX, 4 * sizeof(int16_t));
            x1 = AE_L16X4_I(pX, 0);
            AE_MULFQ16X2_FIR_2(q0, q1, x0, x1, y0);
            AE_MULFQ16X2_FIR_0(q2, q3, x0, x1, y0);
            if (m_iter)
            {
                AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
                x2 = AE_L16X4_I(pX, 4 * sizeof(int16_t));
                AE_MULAFQ16X2_FIR_2(q0, q1, x1, x2, y0);
                AE_MULAFQ16X2_FIR_0(q2, q3, x1, x2, y0);
            }
            t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
            t1 = AE_SAT32X2(q2, q3);
            AE_SA16X4_RIP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
        }
        AE_SA64NEG_FP(ar, pR);

        n_iter &= 3;
        if (n_iter)
        {
            pY = (const ae_int16x4 *)x;
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
            AE_L16X4_IP(x0, pX, 4 * sizeof(int16_t));
            AE_L16X4_IP(x1, pX, 4 * sizeof(int16_t));
            AE_MULFQ16X2_FIR_2(q0, q1, x0, x1, y0);
            AE_MULFQ16X2_FIR_0(q2, q3, x0, x1, y0);
            t0 = AE_TRUNCA32X2F64S(q1, q0, 32);
            t1 = AE_SAT32X2(q2, q2);
            x0 = AE_ROUND16X4F32SASYM(t1, t0);
            AE_S16_0_IP(x0, castxcc(ae_int16, pR), -(int)sizeof(int16_t));
            if (n_iter > 1) { x0 = AE_SEL16_4321(x0, x0); AE_S16_0_IP(x0, castxcc(ae_int16, pR), -(int)sizeof(int16_t)); }
            if (n_iter > 2) { x0 = AE_SEL16_4321(x0, x0); AE_S16_0_IP(x0, castxcc(ae_int16, pR), -(int)sizeof(int16_t)); }
        }
    }

    /*
     * Compute r[M]...r[MIN(N + M - 1,(N + 1)&~1] entries.
     */
    {
        int n_iter = XT_MIN(N, ((N - M + 1) + 3)&~3);
        int m_iter = (M + 3) >> 2;

        pX = (const ae_int16x4 *)x;
        pR = (      ae_int16x4 *)(r + M - 1);
        ar = AE_ZALIGN64();

        for (n = 0; n < (n_iter >> 4); n++)
        {
            ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
            ae_int16x4 x3, x4;

            pY = (const ae_int16x4 *)y;
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

            for (m = 0; m < m_iter - 1; m++)
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

        if (n_iter & 8)
        {
            pY = (const ae_int16x4 *)y;
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

            t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
            t1 = AE_SAT32X2(q2, q3);
            AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
            t0 = AE_TRUNCA32X2F64S(q4, q5, 32);
            t1 = AE_SAT32X2(q6, q7);
            AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
        }

        if (n_iter & 4)
        {
            pY = (const ae_int16x4 *)y;
            AE_L16X4_XP(x0, pX, 4 * sizeof(int16_t));
            S0 = pX;
            q0 = q1 = q2 = q3 = AE_ZERO64();
            __Pragma("loop_count min=1");
            for (m = 0; m < m_iter; m++)
            {
                AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
                AE_L16X4_XP(x1, S0, 4 * sizeof(int16_t));
                AE_MULAFQ16X2_FIR_3(q0, q1, x0, x1, y0);
                AE_MULAFQ16X2_FIR_1(q2, q3, x0, x1, y0);
                x0 = x1;
            }
            t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
            t1 = AE_SAT32X2(q2, q3);
            AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), ar, pR);
        }
        AE_SA64POS_FP(ar, pR);

        n_iter &= 3;
        if (n_iter)
        {
            pY = (const ae_int16x4 *)y;
            AE_L16X4_IP(x0, pX, 4 * sizeof(int16_t));
            q0 = q1 = q2 = q3 = AE_ZERO64();
            __Pragma("loop_count min=1");
            for (m = 0; m < m_iter; m++)
            {
                AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
                AE_L16X4_XP(x1, pX, 4 * sizeof(int16_t));
                AE_MULAFQ16X2_FIR_3(q0, q1, x0, x1, y0);
                AE_MULAFQ16X2_FIR_1(q2, q3, x0, x1, y0);
                x0 = x1;
            }
            t0 = AE_TRUNCA32X2F64S(q1, q0, 32);
            t1 = AE_SAT32X2(q2, q2);
            x0 = AE_ROUND16X4F32SASYM(t1, t0);
            AE_S16_0_IP(x0, castxcc(ae_int16, pR), sizeof(int16_t));
            if (n_iter > 1) { x0 = AE_SEL16_4321(x0, x0); AE_S16_0_IP(x0, castxcc(ae_int16, pR), sizeof(int16_t)); }
            if (n_iter > 2) { x0 = AE_SEL16_4321(x0, x0); AE_S16_0_IP(x0, castxcc(ae_int16, pR), sizeof(int16_t)); }
        }
    }

    /*
     * Compute r[MIN(N + M - 1,(N + 1)&~1]....r[N+M-1] entries.
     */
    {
        int offset = (-((N - (M - 1))) & 3);
        int n_iter = XT_MAX(0, (M - 1) - offset);
        int m_iter = ((n_iter + 3) >> 2) - 1;

        pX = (const ae_int16x4 *)(x + N + offset - (M - 1));
        pR = (      ae_int16x4 *)(r + N + offset);
        ar = AE_ZALIGN64();

        for (n = 0; n < (n_iter >> 4); n++, m_iter -= 4)
        {
            ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
            ae_int16x4 x3, x4;

            pY = (const ae_int16x4 *)y;
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

        if (n_iter & 8)
        {
            pY = (const ae_int16x4 *)y;
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

        if (n_iter & 4)
        {
            pY = (const ae_int16x4 *)y;
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

        n_iter &= 3;
        if (n_iter)
        {
            pY = (const ae_int16x4 *)y;
            AE_L16X4_IP(y0, pY, 4 * sizeof(int16_t));
            AE_L16X4_IP(x0, pX, 4 * sizeof(int16_t));
            AE_L16X4_IP(x1, pX, 4 * sizeof(int16_t));
            AE_MULFQ16X2_FIR_3(q0, q1, x0, x1, y0);
            AE_MULFQ16X2_FIR_1(q2, q3, x0, x1, y0);
            t0 = AE_TRUNCA32X2F64S(q1, q0, 32);
            t1 = AE_SAT32X2(q2, q2);
            x0 = AE_ROUND16X4F32SASYM(t1, t0);
            AE_S16_0_IP(x0, castxcc(ae_int16, pR), sizeof(int16_t));
            if (n_iter > 1) { x0 = AE_SEL16_4321(x0, x0); AE_S16_0_IP(x0, castxcc(ae_int16, pR), sizeof(int16_t)); }
            if (n_iter > 2) { x0 = AE_SEL16_4321(x0, x0); AE_S16_0_IP(x0, castxcc(ae_int16, pR), sizeof(int16_t)); }
        }
    }
} /* raw_lxcorr16x16() */
