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
#include <math.h>
#include "common_fpu.h"
#include "cholnf_common.h"
/*
code optimized for HiFi4 with VFPU
*/

#if (HAVE_VFPU)
#define SZ_F32 (sizeof(float32_t))

/* Load 32bit and replicate to xtfloatx2*/
#define _L32_SX2_IP(a,b,c) \
{ \
ae_int32x2 tmp; \
AE_L32_IP(tmp, castxcc(ae_int32, b), c); \
a = XT_AE_MOVXTFLOATX2_FROMINT32X2(tmp); }
/*
compute matrix product Z[L][NxP]=A[L][MxN]'*B[L][MxP]
Input:
A[SA]    matrix MxN
B[SB]    matrix MxP
Output:
Z[N*P]   matrix NxP
*/
static void realcomputeAB10f(float32_t* Z, const float32_t* A, const float32_t* B)
{
	xtfloatx4 * restrict pZ = (xtfloatx4 *)Z;
	const xtfloatx4 * restrict pA = (const xtfloatx4 *)A;
	const xtfloatx2 * restrict pB = (const xtfloatx2 *)B;
	xtfloatx2 A0, A1, A2, A3, A4, A5, A6, A7, A8, A9;
	xtfloatx2 B0, B1;
	xtfloatx2 Z0, Z1, Z2, Z3, Z4;
	xtfloatx2 Z0t, Z1t, Z2t, Z3t, Z4t;
	//int m;
	
	//Z0 = Z1 = Z2 = Z3 = Z4 = XT_CONST_S(0);
	//Z0t = Z1t = Z2t = Z3t = Z4t = XT_CONST_S(0);
	//for (m = 0; m<5; m++)
	{
		AE_LSX2X2_IP(A0, A1, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A2, A3, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A4, A5, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A6, A7, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A8, A9, pA, 4 * SZ_F32);
		_L32_SX2_IP(B0, pB, SZ_F32);
		_L32_SX2_IP(B1, pB, SZ_F32);
		MUL_SX2X2(Z0, Z1, A0, A1, B0, B0);
		MUL_SX2X2(Z2, Z3, A2, A3, B0, B0);
		MUL_SX2X2(Z4, Z0t, A4, A5, B0, B1);
		MUL_SX2X2(Z1t, Z2t, A6, A7, B1, B1);
		MUL_SX2X2(Z3t, Z4t, A8, A9, B1, B1);

		AE_LSX2X2_IP(A0, A1, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A2, A3, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A4, A5, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A6, A7, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A8, A9, pA, 4 * SZ_F32);
		_L32_SX2_IP(B0, pB, SZ_F32);
		_L32_SX2_IP(B1, pB, SZ_F32);
		MADD_SX2X2(Z0, Z1, A0, A1, B0, B0);
		MADD_SX2X2(Z2, Z3, A2, A3, B0, B0);
		MADD_SX2X2(Z4, Z0t, A4, A5, B0, B1);
		MADD_SX2X2(Z1t, Z2t, A6, A7, B1, B1);
		MADD_SX2X2(Z3t, Z4t, A8, A9, B1, B1);

		AE_LSX2X2_IP(A0, A1, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A2, A3, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A4, A5, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A6, A7, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A8, A9, pA, 4 * SZ_F32);
		_L32_SX2_IP(B0, pB, SZ_F32);
		_L32_SX2_IP(B1, pB, SZ_F32);
		MADD_SX2X2(Z0, Z1, A0, A1, B0, B0);
		MADD_SX2X2(Z2, Z3, A2, A3, B0, B0);
		MADD_SX2X2(Z4, Z0t, A4, A5, B0, B1);
		MADD_SX2X2(Z1t, Z2t, A6, A7, B1, B1);
		MADD_SX2X2(Z3t, Z4t, A8, A9, B1, B1);

		AE_LSX2X2_IP(A0, A1, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A2, A3, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A4, A5, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A6, A7, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A8, A9, pA, 4 * SZ_F32);
		_L32_SX2_IP(B0, pB, SZ_F32);
		_L32_SX2_IP(B1, pB, SZ_F32);
		MADD_SX2X2(Z0, Z1, A0, A1, B0, B0);
		MADD_SX2X2(Z2, Z3, A2, A3, B0, B0);
		MADD_SX2X2(Z4, Z0t, A4, A5, B0, B1);
		MADD_SX2X2(Z1t, Z2t, A6, A7, B1, B1);
		MADD_SX2X2(Z3t, Z4t, A8, A9, B1, B1);

		AE_LSX2X2_IP(A0, A1, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A2, A3, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A4, A5, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A6, A7, pA, 4 * SZ_F32);
		AE_LSX2X2_IP(A8, A9, pA, 4 * SZ_F32);
		_L32_SX2_IP(B0, pB, SZ_F32);
		_L32_SX2_IP(B1, pB, SZ_F32);
		MADD_SX2X2(Z0, Z1, A0, A1, B0, B0);
		MADD_SX2X2(Z2, Z3, A2, A3, B0, B0);
		MADD_SX2X2(Z4, Z0t, A4, A5, B0, B1);
		MADD_SX2X2(Z1t, Z2t, A6, A7, B1, B1);
		MADD_SX2X2(Z3t, Z4t, A8, A9, B1, B1);
	}
	ADD_SX2X2(Z0, Z1, Z0, Z1, Z0t, Z1t);
	ADD_SX2X2(Z2, Z3, Z2, Z3, Z2t, Z3t);
	Z4 = XT_ADD_SX2(Z4, Z4t);
	AE_SSX2X2_IP(Z0, Z1, pZ, 4 * SZ_F32);
	AE_SSX2X2_IP(Z2, Z3, pZ, 4 * SZ_F32);
	XT_SSX2IP(Z4, castxcc(xtfloatx2, pZ), 2 * SZ_F32);
	

}
/*-------------------------------------------------------------------------
Cholesky Forward Substitution for Pseudo-inversion
These functions make forward recursion stage of pseudo-inversion. They use
Cholesky decomposition R[NxN] of original matrices A[MxN]

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
R[((N+1)*N)/2]
            upper triangular matrix R. For fixed point, representation is Q(qR)
A[M*N]      matrix A. For fixed point, representation is Q(qA)
B[M*P]      original right-side matrix B. For fixed point, representation is Q(qB)
D[N]        reciprocals of main diagonal. NOTE: for the fixed point API,
            these data are stored internally in special format with separate
            mantissa and exponent for better accuracy and dynamic range 
            control. So, even for the real data, they stored as pairs of 2
            integers and packed to the complex_fract32 format
qYBRA       qY-qB+qR-qA, combination of fixed point representations of 
            matrices y,B,R and A (for the fixed point API only)
Output:
y[N*P]      Decision matrix y. For fixed point, representation is Q(qY)
Temporary:
pScr        Scratch memory

N = M = 4, 6, 8, 10
P = 1

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
16-bytes boundary
---------------------------------------------------------------------------*/
void  matcholfwdsubst10x10f(void * pScr,
	float32_t * y,
	const float32_t * R,
	const float32_t * D,
	const float32_t * A,
	const float32_t * B)
{
	float32_t* _Z = (float32_t*)pScr;

	NASSERT(y);
	NASSERT(R);
	NASSERT(D);
	NASSERT(A);
	NASSERT(B);
	NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(B, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

	// compute A'*B
	realcomputeAB10f(_Z, A, B);

	// forward recursion
	realcholFwdrec10f(y, R, D, _Z, 1);
}


size_t  matcholfwdsubst10x10f_getScratchSize()
{
	return 10 * sizeof(float32_t);
}
#elif (HAVE_FPU)
#define SZ_F32 (sizeof(float32_t))
/*
compute matrix product Z[L][NxP]=A[L][MxN]'*B[L][MxP]
Input:
A[SA]    matrix MxN
B[SB]    matrix MxP
Output:
Z[N*P]   matrix NxP
*/
static void realcomputeAB10f(float32_t* Z, const float32_t* A, const float32_t* B)
{
    float32_t B_re, B_im;
    int n, m;
    for (n = 0; n<10; n++)
    {
        B_re = B_im = 0;
        for (m = 0; m<10; m++)
        {
            float32_t a_re, b_re;
            a_re = A[n + m * 10];
            b_re = B[m];
            B_re += (a_re*b_re);
        }
        Z[n] = B_re;
    }
}
/*-------------------------------------------------------------------------
Cholesky Forward Substitution for Pseudo-inversion
These functions make forward recursion stage of pseudo-inversion. They use
Cholesky decomposition R[NxN] of original matrices A[MxN]

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
R[((N+1)*N)/2]
upper triangular matrix R. For fixed point, representation is Q(qR)
A[M*N]      matrix A. For fixed point, representation is Q(qA)
B[M*P]      original right-side matrix B. For fixed point, representation is Q(qB)
D[N]        reciprocals of main diagonal. NOTE: for the fixed point API,
these data are stored internally in special format with separate
mantissa and exponent for better accuracy and dynamic range
control. So, even for the real data, they stored as pairs of 2
integers and packed to the complex_fract32 format
qYBRA       qY-qB+qR-qA, combination of fixed point representations of
matrices y,B,R and A (for the fixed point API only)
Output:
y[N*P]      Decision matrix y. For fixed point, representation is Q(qY)
Temporary:
pScr        Scratch memory

N = M = 4, 6, 8, 10
P = 1

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
16-bytes boundary
---------------------------------------------------------------------------*/
void  matcholfwdsubst10x10f(void * pScr,
    float32_t * y,
    const float32_t * R,
    const float32_t * D,
    const float32_t * A,
    const float32_t * B)
{
    float32_t* _Z = (float32_t*)pScr;

    NASSERT(y);
    NASSERT(R);
    NASSERT(D);
    NASSERT(A);
    NASSERT(B);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(B, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    // compute A'*B
    realcomputeAB10f(_Z, A, B);

    // forward recursion
    realcholFwdrec10f(y, R, D, _Z, 1);
}


size_t  matcholfwdsubst10x10f_getScratchSize()
{
    return 10 * sizeof(float32_t);
}
#else
DISCARD_FUN(void, matcholfwdsubst10x10f, (void * pScr,
	float32_t * y,
	const float32_t * R,
	const float32_t * D,
	const float32_t * A,
	const float32_t * B))

	size_t  matcholfwdsubst10x10f_getScratchSize()
{
	return 0;
}
#endif
