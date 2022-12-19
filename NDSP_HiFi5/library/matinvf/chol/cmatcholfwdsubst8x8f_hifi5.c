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
#define SZ_CF32 (2*SZ_F32)
/*
compute matrix product Z[Nx1]=A[MxN]'*B[Mx1]
Input:
A[SA]    complex matrix MxN
B[SB]    complex matrix Mx1
Output:
Z[N*P]   complex matrix Nx1
*/
void cplxcomputeAB8f(complex_float * Z,
	const complex_float * A,
	const complex_float * B)
{
    xtfloatx4 * restrict pZ = (xtfloatx4 *)Z;
    const xtfloatx4 * restrict pA = (const xtfloatx4 *)A;
    const xtfloatx4 * restrict pB = (const xtfloatx4 *)B;

    xtfloatx2 A0, A1, A2, A3, A4, A5, A6, A7;
	xtfloatx2 B0, B1;
    xtfloatx2 Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7;
    xtfloatx2 Z0t, Z1t, Z2t, Z3t, Z4t, Z5t, Z6t, Z7t;
    int m;
	
	
	
    Z0 = Z1 = Z2 = Z3 = Z4 = Z5 = Z6 = Z7 = XT_CONST_S(0);
	Z0t = Z1t = Z2t = Z3t = Z4t = Z5t = Z6t = Z7t = XT_CONST_S(0);
	for (m = 0; m<8; m+=2)
	{
		AE_LSX2X2_IP(A0, A1, pA, 2 * SZ_CF32);
		AE_LSX2X2_IP(A2, A3, pA, 2 * SZ_CF32);
		AE_LSX2X2_IP(A4, A5, pA, 2 * SZ_CF32);
		AE_LSX2X2_IP(A6, A7, pA, 2 * SZ_CF32);
		AE_LSX2X2_IP(B0, B1, pB, 2 * SZ_CF32);
		MADDMUX_SX2X2(Z0, Z1, A0, A1, B0, B0, 0);
		MADDMUX_SX2X2(Z0t, Z1t, A0, A1, B0, B0, 3);
		MADDMUX_SX2X2(Z2, Z3, A2, A3, B0, B0, 0);
		MADDMUX_SX2X2(Z2t, Z3t, A2, A3, B0, B0, 3);
		MADDMUX_SX2X2(Z4, Z5, A4, A5, B0, B0, 0);
		MADDMUX_SX2X2(Z4t, Z5t, A4, A5, B0, B0, 3);
		MADDMUX_SX2X2(Z6, Z7, A6, A7, B0, B0, 0);
		MADDMUX_SX2X2(Z6t, Z7t, A6, A7, B0, B0, 3);

		AE_LSX2X2_IP(A0, A1, pA, 2 * SZ_CF32);
		AE_LSX2X2_IP(A2, A3, pA, 2 * SZ_CF32);
		AE_LSX2X2_IP(A4, A5, pA, 2 * SZ_CF32);
		AE_LSX2X2_IP(A6, A7, pA, 2 * SZ_CF32);
		MADDMUX_SX2X2(Z0, Z1, A0, A1, B1, B1, 0);
		MADDMUX_SX2X2(Z0t, Z1t, A0, A1, B1, B1, 3);
		MADDMUX_SX2X2(Z2, Z3, A2, A3, B1, B1, 0);
		MADDMUX_SX2X2(Z2t, Z3t, A2, A3, B1, B1, 3);
		MADDMUX_SX2X2(Z4, Z5, A4, A5, B1, B1, 0);
		MADDMUX_SX2X2(Z4t, Z5t, A4, A5, B1, B1, 3);
		MADDMUX_SX2X2(Z6, Z7, A6, A7, B1, B1, 0);
		MADDMUX_SX2X2(Z6t, Z7t, A6, A7, B1, B1, 3);
	}
	ADD_SX2X2(Z0, Z1, Z0, Z1, Z0t, Z1t);
	ADD_SX2X2(Z2, Z3, Z2, Z3, Z2t, Z3t);
	ADD_SX2X2(Z4, Z5, Z4, Z5, Z4t, Z5t);
	ADD_SX2X2(Z6, Z7, Z6, Z7, Z6t, Z7t);
	AE_SSX2X2_IP(Z0, Z1, pZ, 2 * SZ_CF32);
	AE_SSX2X2_IP(Z2, Z3, pZ, 2 * SZ_CF32);
	AE_SSX2X2_IP(Z4, Z5, pZ, 2 * SZ_CF32);
	AE_SSX2X2_IP(Z6, Z7, pZ, 2 * SZ_CF32);
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
void  cmatcholfwdsubst8x8f(void * pScr,
    complex_float * y,
    const complex_float * R,
    const complex_float * D,
    const complex_float * A,
    const complex_float * B)
{
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

    complex_float * _Z = (complex_float *)pScr;

    // compute A'*B
    cplxcomputeAB8f(_Z, A, B);

    // forward recursion
    cplxcholFwdrec8f(y, R, D, _Z, 1);
}


size_t  cmatcholfwdsubst8x8f_getScratchSize()
{
	return 2 * 8 * sizeof(float32_t);
}

#elif (HAVE_FPU)
#define SZ_F32 (sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
/*
compute matrix product Z[Nx1]=A[MxN]'*B[Mx1]
Input:
A[SA]    complex matrix MxN
B[SB]    complex matrix Mx1
Output:
Z[N*1]   complex matrix Nx1
*/
void cplxcomputeAB8f(complex_float * Z,
    const complex_float * A,
    const complex_float * B)
{
    float32_t* restrict _Z = (float32_t*)Z;
    const float32_t* restrict _B = (const float32_t*)B;
    const float32_t* restrict _A = (const float32_t*)A;
    float32_t B_re, B_im;
    int n, m;

    for (n = 0; n<8; n++)
    {
        B_re = B_im = 0;
        for (m = 0; m<8; m++)
        {
            float32_t a_re, a_im, b_re, b_im;
            a_re = _A[2 * n + m * 8 * 2 + 0]; a_im = _A[2 * n + m * 8 * 2 + 1];
            b_re = _B[m * 1 * 2 + 0]; b_im = _B[m * 1 * 2 + 1];
            B_re += (a_re*b_re) + (a_im*b_im);
            B_im += (a_re*b_im) - (a_im*b_re);
        }
        _Z[2 * n * 1 + 0] = B_re;
        _Z[2 * n * 1 + 1] = B_im;
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
void  cmatcholfwdsubst8x8f(void * pScr,
    complex_float * y,
    const complex_float * R,
    const complex_float * D,
    const complex_float * A,
    const complex_float * B)
{
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

    complex_float * _Z = (complex_float *)pScr;

    // compute A'*B
    cplxcomputeAB8f(_Z, A, B);

    // forward recursion
    cplxcholFwdrec8f(y, R, D, _Z, 1);
}

size_t  cmatcholfwdsubst8x8f_getScratchSize()
{
    return 2 * 8 * sizeof(float32_t);
}
#else
DISCARD_FUN(void, cmatcholfwdsubst8x8f, (void * pScr,
	complex_float * y,
	const complex_float * R,
	const complex_float * D,
	const complex_float * A,
	const complex_float * B))

	size_t  cmatcholfwdsubst8x8f_getScratchSize()
{
	return 0;
}
#endif
