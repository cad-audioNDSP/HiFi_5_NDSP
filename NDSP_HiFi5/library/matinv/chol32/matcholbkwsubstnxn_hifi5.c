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
#include "NatureDSP_types.h"
#include "common.h"
#include "NatureDSP_Signal_matinv.h"
#include "chol32x32_common.h"

/*-----------------------------------------------
backward recursion:
makes backward recursion from transposed matrix
yt[NxN]
Input:
yt[NxN]    - transposed right part
Rt         - transposed Cholesky decomposition
D[N]       - diagonal terms
Output:
x[NxN] 
-----------------------------------------------*/
void matcholbkwsubstnxn(int32_t * restrict x,
                     const int32_t * Rt,
                     const complex_fract32 * restrict D,
                     const int32_t * restrict yt, 
                     int N, int qXYR)
{
    int n, m, k;
    int d_exp;

    const ae_int32* restrict pY;
    const ae_int32* restrict pRt1;
    const ae_int32* restrict pRt2;
    const int32_t* restrict pD;
    ae_int32* restrict pX;
    ae_int32* restrict pX1;
    ae_int32* restrict pX2;
    ae_int32x2 r0, x0, x1, y0, y1;
    ae_int64  A_re0, A_re1;
    ae_int64 B_re0, B_re1;
    ae_ep ep_re0, ep_re1;
    ae_int32x2 d0;

    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(Rt, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(yt, HIFI_SIMD_WIDTH);

    pRt1 = (const ae_int32*)(Rt + N*(N+1)/2);
    //pRt1 = (const ae_int32x2*)Rt + N*(N+1)/2-2;
    pD = (const int32_t*)&D[N];
    for (k = N - 1; k >= 0; k--)
    {
        pY = (const ae_int32*)& yt[k];
        pX1 = (ae_int32*)& x[((N - 1) * N)];
        pRt1--;
        pD -= 2;
        for (n = 0; n < N; n+=2)
        {
            AE_L32_XP(y0, pY, N * sizeof(ae_int32));
            B_re0 = AE_SLAA64S(AE_MUL32_HH(y0, 1), qXYR);
            ep_re0 = AE_SEXT72(B_re0);

            AE_L32_XP(y1, pY, N * sizeof(ae_int32));
            B_re1 = AE_SLAA64S(AE_MUL32_HH(y1, 1), qXYR);
            ep_re1 = AE_SEXT72(B_re1);

            pX = pX1;
            pX2 = pX1 + 1;
            pRt2 = pRt1;
            for (m = N - k - 2; m >= 0; m--)
            {
                AE_L32_XP(x0, pX, -N * (int)sizeof(ae_int32));
                AE_L32_XP(r0, pRt2, -(k + m + 1) * (int)sizeof(ae_int32));
                AE_MULS32EP_HH(ep_re0, B_re0, r0, x0);

                AE_L32_XP(x1, pX2, -N * (int)sizeof(ae_int32));
                AE_MULS32EP_HH(ep_re1, B_re1, r0, x1);
            }
            d_exp = pD[1];
            AE_L32_IP(d0, castxcc(ae_int32, pD), 0);

            A_re0 = AE_SRAI72(ep_re0, B_re0, 4);
            AE_MUL32USEP_LH(ep_re0, B_re0, AE_MOVINT32X2_FROMINT64(A_re0), d0);
            B_re0 = AE_SRAI72(ep_re0, B_re0, 31);
            AE_MULAF32S_HH(B_re0, AE_MOVINT32X2_FROMINT64(A_re0), d0);
            AE_S32_L_IP(AE_TRUNCA32F64S(B_re0, d_exp + 5), pX, sizeof(ae_int32));


            A_re1 = AE_SRAI72(ep_re1, B_re1, 4);
            AE_MUL32USEP_LH(ep_re1, B_re1, AE_MOVINT32X2_FROMINT64(A_re1), d0);
            B_re1 = AE_SRAI72(ep_re1, B_re1, 31);
            AE_MULAF32S_HH(B_re1, AE_MOVINT32X2_FROMINT64(A_re1), d0);
            AE_S32_L_IP(AE_TRUNCA32F64S(B_re1, d_exp + 5), pX2, sizeof(ae_int32));
            pX1 += 2;
        }
    }
}
