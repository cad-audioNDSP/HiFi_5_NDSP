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
#include "raw_lxcorr32x32.h"

/*-----------------------------------------------------
raw linear correlation:
Input:
x[N] padded with extra 3 zeroes
y[M] padded with extra 3 zeroes
Output:
r[N+M-1]
restriction:
M should be >0
-----------------------------------------------------*/
void fir_lxcorr32x32(int32_t * r, const int32_t * restrict x, const int32_t * restrict y, int N, int M)
{

    const ae_int32x2 *          pX;
    const ae_int32x2 *          S0;
    const ae_int32x2 *          S1;
    const ae_int32x2 *          S2;
    const ae_int32x2 *          pY;
          ae_int32x4 * restrict pR;
    ae_valignx2 ar2;
    ae_valign ar;

    int n, m;

    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 8);
    NASSERT(M > 0);
    NASSERT(N > 0);

    /*
     * Compute first M-1 entries.
     */
    {
        int n_iter = M - 1;
        int m_iter = (n_iter - 3) >> 2;

        pX = (const ae_int32x2 *)y;
        pR = (      ae_int32x4 *)(r + M - 1 - 1);
        ar = AE_ZALIGN64();

        for (n = 0; n < (n_iter >> 3); n++, m_iter -= 2)
        {
            ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
            ae_int32x2 x0, x1, x2, x3, x4, x5;
            ae_int32x2 y0, y1;

            pY = (const ae_int32x2 *)x;
            S0 = pX;
            S1 = pX + 2;
            pX += 4;
            S2 = pX;

            q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

            __Pragma("loop_count min=1");
            for (m = 0; m < m_iter; m++)
            {
                AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
                AE_L32X2_IP(y1, pY, 2 * sizeof(int32_t));

                AE_L32X2_IP(x0, S0, 2 * sizeof(int32_t));
                AE_L32X2_IP(x1, S0, 2 * sizeof(int32_t));
                AE_L32X2_IP(x2, S1, 2 * sizeof(int32_t));
                AE_L32X2_IP(x3, S1, 2 * sizeof(int32_t));
                AE_L32X2_IP(x4, S2, 2 * sizeof(int32_t));
                AE_L32X2_IP(x5, S2, 2 * sizeof(int32_t));
                //AE_L32X2X2_IP(x0, x1, castxcc(ae_int32x4, S0), 4 * sizeof(int32_t));
                //AE_L32X2X2_IP(x2, x3, castxcc(ae_int32x4, S1), 4 * sizeof(int32_t));
                //AE_L32X2X2_IP(x4, x5, castxcc(ae_int32x4, S2), 4 * sizeof(int32_t));

                AE_MULAFD32X2RA_FIR_L(q0, q1, x0, x1, y0);
                AE_MULAFD32X2RA_FIR_L(q2, q3, x1, x2, y0);
                AE_MULAFD32X2RA_FIR_L(q4, q5, x2, x3, y0);
                AE_MULAFD32X2RA_FIR_L(q6, q7, x3, x4, y0);

                AE_MULAFD32X2RA_FIR_L(q0, q1, x1, x2, y1);
                AE_MULAFD32X2RA_FIR_L(q2, q3, x2, x3, y1);
                AE_MULAFD32X2RA_FIR_L(q4, q5, x3, x4, y1);
                AE_MULAFD32X2RA_FIR_L(q6, q7, x4, x5, y1);
            }
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            AE_L32X2_IP(y1, pY, 2 * sizeof(int32_t));

            AE_L32X2_IP(x0, S0, 2 * sizeof(int32_t));
            AE_L32X2_IP(x1, S0, 2 * sizeof(int32_t));
            AE_L32X2_IP(x2, S1, 2 * sizeof(int32_t));
            AE_L32X2_IP(x3, S1, 2 * sizeof(int32_t));
            //AE_L32X2X2_IP(x0, x1, castxcc(ae_int32x4, S0), 4 * sizeof(int32_t));
            //AE_L32X2X2_IP(x2, x3, castxcc(ae_int32x4, S1), 4 * sizeof(int32_t));

            AE_MULAFD32X2RA_FIR_L(q0, q1, x0, x1, y0);
            AE_MULAFD32X2RA_FIR_L(q2, q3, x1, x2, y0);
            AE_MULAFD32X2RA_FIR_L(q4, q5, x2, x3, y0);

            AE_MULAFD32X2RA_FIR_L(q0, q1, x1, x2, y1);
            AE_MULAFD32X2RA_FIR_L(q2, q3, x2, x3, y1);

            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_L(q0, q1, x2, x3, y0);

            x0 = AE_ROUND32X2F48SASYM(q0, q1);
            x1 = AE_ROUND32X2F48SASYM(q2, q3);
            x2 = AE_ROUND32X2F48SASYM(q4, q5);
            x3 = AE_ROUND32X2F48SASYM(q6, q7);
            AE_SA32X2_RIP(x0, ar, castxcc(ae_int32x2, pR));
            AE_SA32X2_RIP(x1, ar, castxcc(ae_int32x2, pR));
            AE_SA32X2_RIP(x2, ar, castxcc(ae_int32x2, pR));
            AE_SA32X2_RIP(x3, ar, castxcc(ae_int32x2, pR));
        }

        //m_iter = ((n_iter & 7) - 1) >> 1;
        if (n_iter & 4)
        {
            ae_f64 q0, q1, q2, q3;
            ae_int32x2 x0, x1, x2, y0;

            pY = (const ae_int32x2 *)x;
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));

            AE_L32X2_IP(x0, pX, 2 * sizeof(int32_t));
            AE_L32X2_IP(x1, pX, 2 * sizeof(int32_t));
            S0 = pX;
            AE_L32X2_IP(x2, S0, 2 * sizeof(int32_t));

            AE_MULFD32X2RA_FIR_L(q0, q1, x0, x1, y0);
            AE_MULFD32X2RA_FIR_L(q2, q3, x1, x2, y0);

            if (n_iter & 3)//m_iter>1
            {
                AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
                x0 = x1; x1 = x2;
                AE_L32X2_IP(x2, S0, 2 * sizeof(int32_t));
                AE_MULAFD32X2RA_FIR_L(q0, q1, x0, x1, y0);
                AE_MULAFD32X2RA_FIR_L(q2, q3, x1, x2, y0);
            }
            if ((n_iter & 3) == 3)//m_iter>2
            {
                AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
                x0 = x1; x1 = x2;
                AE_L32X2_IP(x2, S0, 2 * sizeof(int32_t));
                AE_MULAFD32X2RA_FIR_L(q0, q1, x0, x1, y0);
                AE_MULAFD32X2RA_FIR_L(q2, q3, x1, x2, y0);
            }

            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_L(q0, q1, x1, x2, y0);

            x0 = AE_ROUND32X2F48SASYM(q0, q1);
            x1 = AE_ROUND32X2F48SASYM(q2, q3);
            AE_SA32X2_RIP(x0, ar, castxcc(ae_int32x2, pR));
            AE_SA32X2_RIP(x1, ar, castxcc(ae_int32x2, pR));
            //m_iter -= 2;
        }

        if (n_iter & 2)
        {
            ae_f64 q0, q1;
            ae_int32x2 x0, x1, x2, y0;
            pY = (const ae_int32x2 *)x;
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            AE_L32X2_IP(x0, pX, 2 * sizeof(int32_t));
            x1 = AE_L32X2_I(pX, 0);
            AE_MULFD32X2RA_FIR_L(q0, q1, x0, x1, y0);
            if (n_iter & 1)//m_iter==1
            {
                AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
                x2 = AE_L32X2_I(pX, 2 * sizeof(int32_t));
                AE_MULAFD32X2RA_FIR_L(q0, q1, x1, x2, y0);
            }
            x0 = AE_ROUND32X2F48SASYM(q0, q1);
            AE_SA32X2_RIP(x0, ar, castxcc(ae_int32x2, pR));
        }
        AE_SA64NEG_FP(ar, pR);

        if (n_iter & 1)
        {
            ae_int32x2 x0, y0;
            y0 = AE_L32_I((const ae_int32 *)x, 0);
            x0 = AE_L32_I((const ae_int32 *)pX, sizeof(int32_t));
            x0 = AE_MULFP32X2RAS(x0, y0);
            AE_S32_L_I(x0, (ae_int32 *)pR, 0);
        }
    }

    /*
     * Compute r[M]...r[MIN(N + M - 1,(N + 1)&~1] entries.
     */
    {
        int n_iter = XT_MIN(N, ((N - (M - 1)) + 1)&~1);
        int m_iter = (M + 3) >> 2;

        pX = (const ae_int32x2 *)x;
        pR = (      ae_int32x4 *)(r + M - 1);
        ar2 = AE_ZALIGN128();

        for (n = 0; n < (n_iter >> 3); n++)
        {
            ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
            ae_int32x2 x0, x1, x2, x3, x4, x5;
            ae_int32x2 y0, y1;

            pY = (const ae_int32x2 *)y;
            S0 = pX;
            S1 = pX + 2;
            pX += 4;
            S2 = pX;

            q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

            __Pragma("loop_count min=1");
            for (m = 0; m < m_iter; m++)
            {
                AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
                AE_L32X2_IP(y1, pY, 2 * sizeof(int32_t));

                AE_L32X2X2_IP(x0, x1, castxcc(ae_int32x4, S0), 4 * sizeof(int32_t));
                AE_L32X2X2_IP(x2, x3, castxcc(ae_int32x4, S1), 4 * sizeof(int32_t));
                AE_L32X2X2_IP(x4, x5, castxcc(ae_int32x4, S2), 4 * sizeof(int32_t));

                AE_MULAFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
                AE_MULAFD32X2RA_FIR_H(q2, q3, x1, x2, y0);
                AE_MULAFD32X2RA_FIR_H(q4, q5, x2, x3, y0);
                AE_MULAFD32X2RA_FIR_H(q6, q7, x3, x4, y0);

                AE_MULAFD32X2RA_FIR_H(q0, q1, x1, x2, y1);
                AE_MULAFD32X2RA_FIR_H(q2, q3, x2, x3, y1);
                AE_MULAFD32X2RA_FIR_H(q4, q5, x3, x4, y1);
                AE_MULAFD32X2RA_FIR_H(q6, q7, x4, x5, y1);
            }

            x0 = AE_ROUND32X2F48SASYM(q0, q1);
            x1 = AE_ROUND32X2F48SASYM(q2, q3);
            x2 = AE_ROUND32X2F48SASYM(q4, q5);
            x3 = AE_ROUND32X2F48SASYM(q6, q7);
            AE_SA32X2X2_IP(x0, x1, ar2, pR);
            AE_SA32X2X2_IP(x2, x3, ar2, pR);
        }

        if (n_iter & 4)
        {
            ae_f64 q0, q1, q2, q3;
            ae_int32x2 x0, x1, x2, x3;
            ae_int32x2 y0, y1;

            pY = (const ae_int32x2 *)y;
            S0 = pX;
            pX += 2;
            S1 = pX;

            q0 = q1 = q2 = q3 = AE_ZERO64();

            __Pragma("loop_count min=1");
            for (m = 0; m < m_iter; m++)
            {
                AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
                AE_L32X2_IP(y1, pY, 2 * sizeof(int32_t));
                AE_L32X2X2_IP(x0, x1, castxcc(ae_int32x4, S0), 4 * sizeof(int32_t));
                AE_L32X2X2_IP(x2, x3, castxcc(ae_int32x4, S1), 4 * sizeof(int32_t));
                AE_MULAFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
                AE_MULAFD32X2RA_FIR_H(q2, q3, x1, x2, y0);
                AE_MULAFD32X2RA_FIR_H(q0, q1, x1, x2, y1);
                AE_MULAFD32X2RA_FIR_H(q2, q3, x2, x3, y1);
            }

            x0 = AE_ROUND32X2F48SASYM(q0, q1);
            x1 = AE_ROUND32X2F48SASYM(q2, q3);
            AE_SA32X2X2_IP(x0, x1, ar2, pR);
        }
        AE_SA128POS_FP(ar2, pR);

        if (n_iter & 2)
        {
            ae_f64 q0, q1;
            ae_int32x2 x0, x1, x2, y0, y1;

            pY = (const ae_int32x2 *)y;
            S0 = pX;
            pX += 1;
            S1 = pX + 1;

            q0 = q1 = AE_ZERO64();

            __Pragma("loop_count min=1");
            for (m = 0; m < m_iter; m++)
            {
                AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
                AE_L32X2_IP(y1, pY, 2 * sizeof(int32_t));
                AE_L32X2X2_IP(x0, x1, castxcc(ae_int32x4, S0), 4 * sizeof(int32_t));
                AE_L32X2_IP(x2, S1, 4 * sizeof(int32_t));
                AE_MULAFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
                AE_MULAFD32X2RA_FIR_H(q0, q1, x1, x2, y1);
            }
            AE_S32RA64S_IP(q0, castxcc(ae_int32, pR), sizeof(int32_t));
            AE_S32RA64S_IP(q1, castxcc(ae_int32, pR), sizeof(int32_t));
        }

        if (n_iter & 1)
        {
            ae_int32x2 x0, y0;
            y0 = AE_L32_I((const ae_int32 *)y, 0);
            x0 = AE_L32_I((const ae_int32 *)pX, 0);
            x0 = AE_MULFP32X2RAS(x0, y0);
            AE_S32_L_I(x0, (ae_int32 *)pR, 0);
        }
    }

    /*
     * Compute r[MIN(N + M - 1,(N + 1)&~1]....r[N+M-1] entries.
     */
    {
        int offset = 1 - (((N - (M - 1)) + 1) & 1);
        int n_iter = XT_MAX(0, (M - 1) - offset);
        int m_iter = (n_iter - 3) >> 2;

        pX = (const ae_int32x2 *)(x + N - (M - 1) + offset);
        pR = (      ae_int32x4 *)(r + N + offset);
        ar2 = AE_ZALIGN128();

        for (n = 0; n < (n_iter >> 3); n++, m_iter -= 2)
        {
            ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
            ae_int32x2 x0, x1, x2, x3, x4, x5;
            ae_int32x2 y0, y1;

            pY = (const ae_int32x2 *)y;
            S0 = pX;
            S1 = pX + 2;
            pX += 4;
            S2 = pX;

            q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

            __Pragma("loop_count min=1");
            for (m = 0; m < m_iter; m++)
            {
                AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
                AE_L32X2_IP(y1, pY, 2 * sizeof(int32_t));

                AE_L32X2_IP(x0, S0, 2 * sizeof(int32_t));
                AE_L32X2_IP(x1, S0, 2 * sizeof(int32_t));
                AE_L32X2_IP(x2, S1, 2 * sizeof(int32_t));
                AE_L32X2_IP(x3, S1, 2 * sizeof(int32_t));
                AE_L32X2_IP(x4, S2, 2 * sizeof(int32_t));
                AE_L32X2_IP(x5, S2, 2 * sizeof(int32_t));
                //AE_L32X2X2_IP(x0, x1, castxcc(ae_int32x4, S0), 4 * sizeof(int32_t));
                //AE_L32X2X2_IP(x2, x3, castxcc(ae_int32x4, S1), 4 * sizeof(int32_t));
                //AE_L32X2X2_IP(x4, x5, castxcc(ae_int32x4, S2), 4 * sizeof(int32_t));

                AE_MULAFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
                AE_MULAFD32X2RA_FIR_H(q2, q3, x1, x2, y0);
                AE_MULAFD32X2RA_FIR_H(q4, q5, x2, x3, y0);
                AE_MULAFD32X2RA_FIR_H(q6, q7, x3, x4, y0);

                AE_MULAFD32X2RA_FIR_H(q0, q1, x1, x2, y1);
                AE_MULAFD32X2RA_FIR_H(q2, q3, x2, x3, y1);
                AE_MULAFD32X2RA_FIR_H(q4, q5, x3, x4, y1);
                AE_MULAFD32X2RA_FIR_H(q6, q7, x4, x5, y1);
            }
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            AE_L32X2_IP(y1, pY, 2 * sizeof(int32_t));

            AE_L32X2_IP(x0, S0, 2 * sizeof(int32_t));
            AE_L32X2_IP(x1, S0, 2 * sizeof(int32_t));
            AE_L32X2_IP(x2, S1, 2 * sizeof(int32_t));
            AE_L32X2_IP(x3, S1, 2 * sizeof(int32_t));
            //AE_L32X2X2_IP(x0, x1, castxcc(ae_int32x4, S0), 4 * sizeof(int32_t));
            //AE_L32X2X2_IP(x2, x3, castxcc(ae_int32x4, S1), 4 * sizeof(int32_t));

            AE_MULAFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
            AE_MULAFD32X2RA_FIR_H(q2, q3, x1, x2, y0);
            AE_MULAFD32X2RA_FIR_H(q4, q5, x2, x3, y0);

            AE_MULAFD32X2RA_FIR_H(q0, q1, x1, x2, y1);
            AE_MULAFD32X2RA_FIR_H(q2, q3, x2, x3, y1);

            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_H(q0, q1, x2, x3, y0);

            x0 = AE_ROUND32X2F48SASYM(q0, q1);
            x1 = AE_ROUND32X2F48SASYM(q2, q3);
            x2 = AE_ROUND32X2F48SASYM(q4, q5);
            x3 = AE_ROUND32X2F48SASYM(q6, q7);
            AE_SA32X2X2_IP(x0, x1, ar2, pR);
            AE_SA32X2X2_IP(x2, x3, ar2, pR);
        }

        //m_iter = ((n_iter & 7) - 1) >> 1;
        if (n_iter & 4)
        {
            ae_f64 q0, q1, q2, q3;
            ae_int32x2 x0, x1, x2, y0;

            pY = (const ae_int32x2 *)y;
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));

            AE_L32X2_IP(x0, pX, 2 * sizeof(int32_t));
            AE_L32X2_IP(x1, pX, 2 * sizeof(int32_t));
            S0 = pX;
            AE_L32X2_IP(x2, S0, 2 * sizeof(int32_t));

            AE_MULFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
            AE_MULFD32X2RA_FIR_H(q2, q3, x1, x2, y0);

            if (n_iter & 3)//m_iter>1
            {
                AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
                x0 = x1; x1 = x2;
                AE_L32X2_IP(x2, S0, 2 * sizeof(int32_t));
                AE_MULAFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
                AE_MULAFD32X2RA_FIR_H(q2, q3, x1, x2, y0);
            }
            if ((n_iter & 3) == 3)//m_iter>2
            {
                AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
                x0 = x1; x1 = x2;
                AE_L32X2_IP(x2, S0, 2 * sizeof(int32_t));
                AE_MULAFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
                AE_MULAFD32X2RA_FIR_H(q2, q3, x1, x2, y0);
            }

            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_H(q0, q1, x1, x2, y0);

            x0 = AE_ROUND32X2F48SASYM(q0, q1);
            x1 = AE_ROUND32X2F48SASYM(q2, q3);
            AE_SA32X2X2_IP(x0, x1, ar2, pR);
            //m_iter -= 2;
        }
        AE_SA128POS_FP(ar2, pR);

        if (n_iter & 2)
        {
            ae_f64 q0, q1;
            ae_int32x2 x0, x1, x2, y0;
            pY = (const ae_int32x2 *)y;
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            AE_L32X2_IP(x0, pX, 2 * sizeof(int32_t));
            x1 = AE_L32X2_I(pX, 0);
            AE_MULFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
            if (n_iter & 1)//m_iter==1
            {
                AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
                x2 = AE_L32X2_I(pX, 2 * sizeof(int32_t));
                AE_MULAFD32X2RA_FIR_H(q0, q1, x1, x2, y0);
            }
            AE_S32RA64S_IP(q0, castxcc(ae_int32, pR), sizeof(int32_t));
            AE_S32RA64S_IP(q1, castxcc(ae_int32, pR), sizeof(int32_t));
        }

        if (n_iter & 1)
        {
            ae_int32x2 x0, y0;
            y0 = AE_L32_I((const ae_int32 *)y, 0);
            x0 = AE_L32_I((const ae_int32 *)pX, 0);
            x0 = AE_MULFP32X2RAS(x0, y0);
            AE_S32_L_I(x0, (ae_int32 *)pR, 0);
        }
    }
} /* fir_lxcorr32x32() */
