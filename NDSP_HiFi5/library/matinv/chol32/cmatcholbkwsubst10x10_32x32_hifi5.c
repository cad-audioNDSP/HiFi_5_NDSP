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
void cmatcholbkwsubst10x10_32x32(void *pScr,
                 complex_fract32* restrict x,
           const complex_fract32* restrict R,
           const complex_fract32* restrict D,
           const complex_fract32* restrict y,
                 int qXYR)
{
	const int N = 10;
	int m, k, p;
	const ae_int32x2 * restrict pR;
	ae_int32x2  r0reim, r1reim, xreim0, xreim1, t, t0, t1, d;
	ae_int64 A_re, A_im;
	ae_int64 B0_re, B0_im, B1_re, B1_im;
	ae_ep    ep0_re, ep0_im, ep1_re, ep1_im;
	int d_exp;
	NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
	D += N - 1;
	R = R + (((N)*(N + 1)) >> 1) - 2;
	y += N - 2;
	x += N - 2;
	for (k = N - 1; k >= 0; k -= 2)
	{
		pR = (const ae_int32x2*)R;
		R -= 2;
		AE_L32X2X2_XP(t1, t0, castxcc(ae_int32x4, y), -2 * (int)sizeof(ae_int32x2));
		B0_re = AE_SLAA64S(AE_MUL32_HH(t0, 1), qXYR); ep0_re = AE_SEXT72(B0_re);
		B0_im = AE_SLAA64S(AE_MUL32_LL(t0, 1), qXYR); ep0_im = AE_SEXT72(B0_im);
		B1_re = AE_SLAA64S(AE_MUL32_HH(t1, 1), qXYR); ep1_re = AE_SEXT72(B1_re);
		B1_im = AE_SLAA64S(AE_MUL32_LL(t1, 1), qXYR); ep1_im = AE_SEXT72(B1_im);
		p = -(N - 1)*(int)sizeof(ae_int32x2);
		for (m = N - 1; m > k; m -= 2)
		{
			AE_L32X2X2_XP(xreim1, xreim0, castxcc(ae_int32x4, x), -2 * (int)sizeof(ae_int32x2));
			r0reim = AE_L32X2_I(pR, sizeof(ae_int32x2));
			AE_L32X2_XP(r1reim, pR, p);
			p += (int)sizeof(ae_int32x2);
			t = AE_SEL32_LH(xreim0, xreim0);
			AE_MULSSD32EP_HH_LL(ep0_im, B0_im, t, r0reim);
			AE_MULSSD32EP_HH_LL(ep1_im, B1_im, t, r1reim);
			AE_MOVT32X2(xreim0, AE_NEG32S(xreim0), AE_MOVBA2(2));
			AE_MULAAD32EP_HH_LL(ep0_re, B0_re, xreim0, r0reim);
			AE_MULAAD32EP_HH_LL(ep1_re, B1_re, xreim0, r1reim);

			r0reim = AE_L32X2_I(pR, sizeof(ae_int32x2));
			AE_L32X2_XP(r1reim, pR, p);
			p += (int)sizeof(ae_int32x2);
			t = AE_SEL32_LH(xreim1, xreim1);
			AE_MULSSD32EP_HH_LL(ep0_im, B0_im, t, r0reim);
			AE_MULSSD32EP_HH_LL(ep1_im, B1_im, t, r1reim);
			AE_MOVT32X2(xreim1, AE_NEG32S(xreim1), AE_MOVBA2(2));
			AE_MULAAD32EP_HH_LL(ep0_re, B0_re, xreim1, r0reim);
			AE_MULAAD32EP_HH_LL(ep1_re, B1_re, xreim1, r1reim);
		}
		r1reim = AE_L32X2_I(pR, 0);
		d_exp = D[0].s.im;
		AE_L32_XP(d, castxcc(ae_int32, D), -2 * (int)sizeof(int32_t));
		A_re = AE_SRAI72(ep0_re, B0_re, 4);
		A_im = AE_SRAI72(ep0_im, B0_im, 4);
		AE_MUL32USEP_LH(ep0_re, B0_re, AE_MOVINT32X2_FROMINT64(A_re), d);
		AE_MUL32USEP_LH(ep0_im, B0_im, AE_MOVINT32X2_FROMINT64(A_im), d);
		B0_re = AE_SRAI72(ep0_re, B0_re, 31);
		B0_im = AE_SRAI72(ep0_im, B0_im, 31);
		AE_MULAF32S_HH(B0_re, AE_MOVINT32X2_FROMINT64(A_re), d);
		AE_MULAF32S_HH(B0_im, AE_MOVINT32X2_FROMINT64(A_im), d);
		d = AE_TRUNCA32X2F64S(B0_re, B0_im, d_exp + 5);
		xreim1 = xreim0 = d;

		AE_MULSSD32EP_HH_LL(ep1_im, B1_im, AE_SEL32_LH(xreim1, xreim1), r1reim);
		AE_MOVT32X2(xreim1, AE_NEG32S(xreim1), AE_MOVBA2(2));
		AE_MULAAD32EP_HH_LL(ep1_re, B1_re, xreim1, r1reim);
		d_exp = D[0].s.im;
		AE_L32_XP(d, castxcc(ae_int32, D), -2 * (int)sizeof(int32_t));
		A_re = AE_SRAI72(ep1_re, B1_re, 4);
		A_im = AE_SRAI72(ep1_im, B1_im, 4);
		AE_MUL32USEP_LH(ep1_re, B1_re, AE_MOVINT32X2_FROMINT64(A_re), d);
		AE_MUL32USEP_LH(ep1_im, B1_im, AE_MOVINT32X2_FROMINT64(A_im), d);
		B1_re = AE_SRAI72(ep1_re, B1_re, 31);
		B1_im = AE_SRAI72(ep1_im, B1_im, 31);
		AE_MULAF32S_HH(B1_re, AE_MOVINT32X2_FROMINT64(A_re), d);
		AE_MULAF32S_HH(B1_im, AE_MOVINT32X2_FROMINT64(A_im), d);
		d = AE_TRUNCA32X2F64S(B1_re, B1_im, d_exp + 5);
		xreim1 = d;
		AE_S32X2X2_XP(xreim1, xreim0, castxcc(ae_int32x4, x), (N - k - 1)*(int)sizeof(ae_int32x2));
	}
} /* cmatcholbkwsubst4x4_32x32() */

size_t cmatcholbkwsubst10x10_32x32_getScratchSize()   
{ 
    return 0;
}
