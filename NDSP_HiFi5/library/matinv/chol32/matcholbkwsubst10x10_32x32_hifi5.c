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

/*-------------------------------------------------------------------------
Cholesky Backward Substitution for Pseudo-inversion
These functions make backward substitution stage of pseudo-inversion. They
use Cholesky decomposition of original matrices and results of forward
substitution.

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
R[((N+1)*N)/2]
            Cholesky upper-triangle matrix R. For fixed point API, the 
            representation is Q(qR)
y[N*P]      Results of forward substitution stage. For fixed point API, 
            the representation is Q(qY)
D[N]        sequence of reciprocals of main diagonal R. NOTE: for the 
            fixed point API, these data are stored internally in special 
            format with separate mantissa and exponent for better accuracy 
            and dynamic range control. So, even for the real data, they 
            stored as pairs of 2 integers and packed to the complex_fract32 
            format
qXYR        combination of fixed point representation (matrices R, x and y) 
            qXYR=qX-qY+qR (for fixed point API only)
Output:
x[N*P]      Decision matrix x. For fixed point API, the representation is Q(qX)

N = 4, 6, 8, 10
P = 1

Restrictions:
All matrices and scratch should not overlap and be aligned on
16-bytes boundary
---------------------------------------------------------------------------*/
void matcholbkwsubst10x10_32x32(void *pScr,
                 int32_t*         restrict x,
           const int32_t*         restrict R,
           const complex_fract32* restrict D,
           const int32_t*         restrict y,
                 int qXYR)
{
	const int N = 10;
	NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);

    int m, k;
    const int32_t* restrict D0 = (const int32_t*)& D[N - 1];
    ae_int32* restrict pX;
    const ae_int32* restrict pY = (const ae_int32*)& y[N - 1];
    ae_ep ep0_re, ep1_re;
    ae_int64 A0_re, A1_re;
    ae_int64 B0_re, B1_re;
    ae_int32x2 r0_re, r1_re;
    ae_int32x2 x0_re, x1_re;
    ae_int32x2 d0, d1;
    int d0_exp, d1_exp;
    ae_int32x2 t;

    const ae_int32* restrict pRt0 = (const ae_int32*)R + N * (N + 1) / 2 - 1;
    const ae_int32* restrict pRt1;
    const ae_int32* restrict pRt2;
    for (k = N - 1; k >= 0; k -= 2)
    {
        AE_L32_XP(t, pY, -(int)sizeof(ae_int32));
        B0_re = AE_SLAA64S(AE_MUL32_HH(t, 1), qXYR);
        ep0_re = AE_SEXT72(B0_re);
        AE_L32_XP(t, pY, -(int)sizeof(ae_int32));
        B1_re = AE_SLAA64S(AE_MUL32_HH(t, 1), qXYR);
        ep1_re = AE_SEXT72(B1_re);
        pRt1 = pRt0;
        pRt2 = pRt0 - 1;
        pRt0 -= 2;
        pX = (ae_int32*)& x[(N - 1)];
        for (m = 0; m < N - k - 1; m += 2)
        {
            AE_L32_XP(x0_re, pX, -(int)sizeof(ae_int32));
            AE_L32_XP(r0_re, pRt1, (-N + 1 + m) * (int)sizeof(ae_int32));
            AE_MULS32EP_HH(ep0_re, B0_re, x0_re, r0_re);

            AE_L32_XP(r1_re, pRt2, (-N + 1 + m) * (int)sizeof(ae_int32));
            AE_MULS32EP_HH(ep1_re, B1_re, x0_re, r1_re);

            AE_L32_XP(x0_re, pX, -(int)sizeof(ae_int32));
            AE_L32_XP(r0_re, pRt1, (-N + 2 + m) * (int)sizeof(ae_int32));
            AE_MULS32EP_HH(ep0_re, B0_re, x0_re, r0_re);

            AE_L32_XP(r1_re, pRt2, (-N + 2 + m) * (int)sizeof(ae_int32));
            AE_MULS32EP_HH(ep1_re, B1_re, x0_re, r1_re);

        }
        d0_exp = D0[1];
        AE_L32_XP(d0, castxcc(ae_int32, D0), -2 * (int)sizeof(ae_int32));
        A0_re = AE_SRAI72(ep0_re, B0_re, 4);
        AE_MUL32USEP_LH(ep0_re, B0_re, AE_MOVINT32X2_FROMINT64(A0_re), d0);
        B0_re = AE_SRAI72(ep0_re, B0_re, 31);
        AE_MULAF32S_HH(B0_re, AE_MOVINT32X2_FROMINT64(A0_re), d0);
        d0 = AE_TRUNCA32F64S(B0_re, d0_exp + 5);
        AE_S32_L_XP(d0, pX, -(int)sizeof(ae_int32));

        x1_re = d0;
        AE_L32_IP(r1_re, pRt2, 0);
        AE_MULS32EP_HH(ep1_re, B1_re, x1_re, r1_re);

        d1_exp = D0[1];
        AE_L32_XP(d1, castxcc(ae_int32, D0), -2 * (int)sizeof(ae_int32));
        A1_re = AE_SRAI72(ep1_re, B1_re, 4);
        AE_MUL32USEP_LH(ep1_re, B1_re, AE_MOVINT32X2_FROMINT64(A1_re), d1);
        B1_re = AE_SRAI72(ep1_re, B1_re, 31);
        AE_MULAF32S_HH(B1_re, AE_MOVINT32X2_FROMINT64(A1_re), d1);
        d1 = AE_TRUNCA32F64S(B1_re, d1_exp + 5);
        AE_S32_L_XP(d1, pX, -(int)sizeof(ae_int32));
    }
} /*  matcholbkwsubst10x10_32x32() */

size_t matcholbkwsubst10x10_32x32_getScratchSize()   { return 0;}
