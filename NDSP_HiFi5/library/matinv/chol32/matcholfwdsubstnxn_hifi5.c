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

/* ---------------------------------------------------------------------
make forward recursion for R'\y=A'
Input:
R,D       - Cholesky decomposition
A[NxN]    - original input maatrix
Output:
yt[NxN]     transposed output
---------------------------------------------------------------------*/
void matcholfwdsubstnxn(int32_t* yt,
                  const int32_t* R, 
                  const complex_fract32* D, 
                  const int32_t* A, int N, int qYBRA)
{
    int n, m, k;
    int d_exp, d_exp1;
    ae_int64 A_re0, A_re1;
    ae_int64 B_re0, B_re1;
    ae_int32x2 r0, r1;
    ae_int32x2 y0, y1;
    ae_int32x2 a0, a1;
    ae_int32x2 d0, d1;
    const ae_int32* restrict pR0;
    const ae_int32* restrict pR1;
    const int32_t* restrict pD;
    const ae_int32* restrict pA;
    ae_int32* restrict pY;
    ae_ep ep0, ep1;

    NASSERT_ALIGN(yt, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);

    for (m = 0; m < N; m++)
    {
        pA = (const ae_int32*)& A[(m * N)];
        pR0 = (const ae_int32*)(R);
        pR1 = (const ae_int32*)(R + 1);
        pD = (const int32_t*) D;
        for (n = 0; n < N; n+=2)
        {
            AE_L32_IP(a0, pA, sizeof(ae_int32));
            AE_L32_IP(a1, pA, sizeof(ae_int32));
            B_re0 = AE_SLAA64S(AE_MUL32_HH(a0, AE_MOVDA32X2(1, -1)), 31 + qYBRA);
            B_re1 = AE_SLAA64S(AE_MUL32_HH(a1, AE_MOVDA32X2(1, -1)), 31 + qYBRA);
            ep0 = AE_SEXT72(B_re0);
            ep1 = AE_SEXT72(B_re1);
            pY = (ae_int32*)& yt[(m * N)];
            ae_valign va0, va1, vay;
            va0 = AE_LA64_PP(pR0);
            va1 = AE_LA64_PP(pR1);
            vay = AE_LA64_PP(pY);
            for (k = 0; k < n; k+=2)
            {
                AE_LA32X2_IP(r0, va0, castxcc(ae_int32x2, pR0));
                AE_LA32X2_IP(r1, va1, castxcc(ae_int32x2, pR1));
                AE_LA32X2_IP(y0, vay, castxcc(ae_int32x2, pY));
                AE_MULSSD32EP_HH_LL(ep0, B_re0, r0, y0);
                AE_MULSSD32EP_HH_LL(ep1, B_re1, r1, y0);
            }
            d_exp = pD[1];
            AE_L32_XP(d0, castxcc(ae_int32, pD), sizeof(ae_int32x2));
            A_re0 = AE_SRAI72(ep0, B_re0, 4);
            AE_MUL32USEP_LH(ep0, B_re0, AE_MOVINT32X2_FROMINT64(A_re0), d0);
            B_re0 = AE_SRAI72(ep0, B_re0, 31);
            AE_MULAF32S_HH(B_re0, AE_MOVINT32X2_FROMINT64(A_re0), d0);
            y1 = AE_TRUNCA32F64S(B_re0, d_exp + 5);
            AE_S32_L_IP(y1, pY, sizeof(ae_int32));

            AE_L32_IP(r1, pR1, sizeof(ae_int32));
            AE_MULS32EP_HH(ep1, B_re1, r1, y1);
            d_exp1 = pD[1];
            AE_L32_XP(d1, castxcc(ae_int32, pD), sizeof(ae_int32x2));
            A_re1 = AE_SRAI72(ep1, B_re1, 4);
            AE_MUL32USEP_LH(ep1, B_re1, AE_MOVINT32X2_FROMINT64(A_re1), d1);
            B_re1 = AE_SRAI72(ep1, B_re1, 31);
            AE_MULAF32S_HH(B_re1, AE_MOVINT32X2_FROMINT64(A_re1), d1);
            AE_S32_L_IP(AE_TRUNCA32F64S(B_re1, d_exp1 + 5), pY, sizeof(ae_int32));
            pR0 += n + 3;
            pR1 += n + 4;
        }
    }
}


