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
#include "chol32x32_common.h"

/*--------------------------------------------------------
   make forward recursion to update n new column elements
   Input:
   Z[L][N+1][C]  convolutions in N-th column, Q(2*qR-4)
   D[L][SD]      reciprocals of main diagonal
   stride        width of matrix Z
   C             1 for real, 2 for complex
   Output:
   y             pointer to the N-th column in R[L][SR]
				 (only N elements filled), Q(qR)
--------------------------------------------------------*/
void matcholfwdrec(int32_t *restrict y, const int32_t * restrict R, const int32_t * restrict D, const int64_t *Z, int N, int stride)
{
	int m, n;
	const ae_int64* restrict pZ = (const ae_int64*)Z;
	const ae_int32* restrict pR0;
	const ae_int32* restrict pR1;
	ae_int32* restrict pY = (ae_int32*)y;
	ae_ep ep_re0, ep_re1;
	ae_int32x2 y0, r0, d0;
	ae_int32x2 y1, r1, d1;
	ae_int64 A_re0, A_re1;
	ae_int64 B_re0, B_re1;
	int d_exp0, d_exp1;

	//NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(Z, sizeof(ae_int64));

	for (n = 0; n < (N&~1); n += 2)
	{
		ep_re0 = ep_re1 = AE_MOVEA(0);
		B_re0 = B_re1 = 0;
		R += n;
		AE_L64_XP(B_re0, pZ, stride*sizeof(ae_int64));
		AE_L64_XP(B_re1, pZ, stride*sizeof(ae_int64));
        AE_SLAI72(ep_re0, B_re0, B_re0, 4);
        AE_SLAI72(ep_re1, B_re1, B_re1, 4);
		pY = (ae_int32*)y;
		ae_valign va0, va1, vay;
		pR0 = (const ae_int32*)R;
		pR1 = (const ae_int32*)R + (n + 1);
		va0 = AE_LA64_PP(pR0);
        va1 = AE_LA64_PP(pR1);
        vay = AE_LA64_PP(pY);
        for (m = 0; m < n; m += 2)
		{
			AE_LA32X2_IP(r0, va0, castxcc(ae_int32x2, pR0));
			AE_LA32X2_IP(r1, va1, castxcc(ae_int32x2, pR1));
			AE_LA32X2_IP(y0, vay, castxcc(ae_int32x2, pY));
			AE_MULSSD32EP_HH_LL(ep_re0, B_re0, r0, y0);
			AE_MULSSD32EP_HH_LL(ep_re1, B_re1, r1, y0);
		}
		d_exp0 = D[1];
		AE_L32_IP(d0, castxcc(ae_int32, D), sizeof(ae_int32x2));
		A_re0 = AE_SRAI72(ep_re0, B_re0, 4);
		AE_MUL32USEP_LH(ep_re0, B_re0, AE_MOVINT32X2_FROMINT64(A_re0), d0);
		B_re0 = AE_SRAI72(ep_re0, B_re0, 31);
		AE_MULAF32S_HH(B_re0, AE_MOVINT32X2_FROMINT64(A_re0), d0);
		y0 = AE_TRUNCA32F64S(B_re0, d_exp0 + 5);
		AE_S32_L_IP(y0, pY, sizeof(ae_int32));

		AE_L32_IP(r1, pR1, 0);
		AE_MULS32EP_HH(ep_re1, B_re1, r1, y0);
		d_exp1 = D[1];
		AE_L32_IP(d1, castxcc(ae_int32, D), sizeof(ae_int32x2));
		A_re1 = AE_SRAI72(ep_re1, B_re1, 4);
		AE_MUL32USEP_LH(ep_re1, B_re1, AE_MOVINT32X2_FROMINT64(A_re1), d1);
		B_re1 = AE_SRAI72(ep_re1, B_re1, 31);
		AE_MULAF32S_HH(B_re1, AE_MOVINT32X2_FROMINT64(A_re1), d1);
		y1 = AE_TRUNCA32F64S(B_re1, d_exp1 + 5);
		AE_S32_L_IP(y1, pY, sizeof(ae_int32));
		R += n+1;
	}

    if (N & 1)
    {
        ep_re0 = ep_re1 = AE_MOVEA(0);
        B_re0 = B_re1 = 0;
        R += n;
        AE_L64_XP(B_re0, pZ, stride * sizeof(ae_int64));
        AE_SLAI72(ep_re0, B_re0, B_re0, 4);
        pR0 = (const ae_int32*)R;
        pY = (ae_int32*)y;
        for (m = 0; m < n; m++)
        {
            AE_L32_IP(r0, pR0, sizeof(ae_int32));
            AE_L32_IP(y0, pY, sizeof(ae_int32));
            AE_MULS32EP_HH(ep_re0, B_re0, r0, y0);
        }
        d_exp0 = D[1];
        AE_L32_IP(d0, castxcc(ae_int32, D), 2 * sizeof(ae_int32));

        A_re0 = AE_SRAI72(ep_re0, B_re0, 4);
        AE_MUL32USEP_LH(ep_re0, B_re0, AE_MOVINT32X2_FROMINT64(A_re0), d0);
        B_re0 = AE_SRAI72(ep_re0, B_re0, 31);
        AE_MULAF32S_HH(B_re0, AE_MOVINT32X2_FROMINT64(A_re0), d0);

        AE_S32_L_IP(AE_TRUNCA32F64S(B_re0, d_exp0 + 5), pY, sizeof(ae_int32));
    }
}
