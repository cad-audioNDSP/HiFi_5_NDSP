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
/*
code optimized for HiFi4 with VFPU
*/

#if (HAVE_VFPU)
#define SZ_F32 (int)(sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
#define CF32_WIDTH (HIFI_SIMD_WIDTH/SZ_CF32)
/*
Forward recursion for pseudo inversion
Input:
R[((N+1)*N)/2]
upper triangular matrix R
Z[N*P]      matrix A
D[N]        reciprocals of main diagonal
Output:
y[N*P]		Decision matrix y
*/
static void fwd10x10f(
    complex_float* restrict y,
    const complex_float* restrict R,
    const complex_float* restrict D,
    const complex_float* restrict Z)
{
#if 0
    int n, m, p;
    xtfloatx2 Yw, Yr, R0, D0;
    xtfloatx2 * restrict pY;
    const xtfloatx2 * restrict pR;
    const xtfloatx2 * restrict pD;
    const xtfloatx2 * restrict pA;
    for (n = 0; n<10; n++)
        for (p = 0; p<10; p++)
        {
            pR = (const xtfloatx2 *)(R + (n*(n + 1)) / 2);
            pA = (const xtfloatx2 *)(Z + n + p * 10);
            XT_LSX2IP(Yw, pA, SZ_CF32);
            Yw = XT_CONJC_S(Yw);
            pY = (xtfloatx2*)(y + p);
            for (m = 0; m<n; m++)
            {
                XT_LSX2XP(Yr, pY, 10 * SZ_CF32);
                XT_LSX2IP(R0, pR, SZ_CF32);
                XT_MADDMUX_S(Yw, R0, Yr, 2);
                XT_MADDMUX_S(Yw, R0, Yr, 1);
            }
            pD = (const xtfloatx2 *)(D + n);
            XT_LSX2IP(D0, pD, SZ_CF32);
            Yw = XT_MUL_SX2(Yw, D0);
            XT_SSX2IP(Yw, pY, SZ_CF32);
        }
#else
	int n, m;
	xtfloatx2 R0, D0;
	xtfloatx2 Yr0, Yr1, Yr2, Yr3, Yr4, Yr5, Yr6, Yr7, Yr8, Yr9;
	xtfloatx2 Yw0, Yw1, Yw2, Yw3, Yw4, Yw5, Yw6, Yw7, Yw8, Yw9;
	xtfloatx2 Yw0t, Yw1t, Yw2t, Yw3t, Yw4t, Yw5t, Yw6t, Yw7t, Yw8t, Yw9t;
	xtfloatx4 * restrict pY;
	const xtfloatx2 * restrict pR = (const xtfloatx2 *)R;
	const xtfloatx2 * restrict pD = (const xtfloatx2 *)D;
	const xtfloatx2 * restrict pA = (const xtfloatx2 *)Z;
	for (n = 0; n < 10; n++)
		//for (p = 0; p<P; p++)
	{
		pY = (xtfloatx4*)y;
		XT_LSX2XP(Yw0, pA, 10 * SZ_CF32);
		XT_LSX2XP(Yw1, pA, 10 * SZ_CF32);
		XT_LSX2XP(Yw2, pA, 10 * SZ_CF32);
		XT_LSX2XP(Yw3, pA, 10 * SZ_CF32);
		XT_LSX2XP(Yw4, pA, 10 * SZ_CF32);
		XT_LSX2XP(Yw5, pA, 10 * SZ_CF32);
		XT_LSX2XP(Yw6, pA, 10 * SZ_CF32);
		XT_LSX2XP(Yw7, pA, 10 * SZ_CF32);
		XT_LSX2XP(Yw8, pA, 10 * SZ_CF32);
		XT_LSX2XP(Yw9, pA, (-10 * 9 + 1)*SZ_CF32);
		CONJC_SX2X2(Yw0, Yw1, Yw0, Yw1);
		CONJC_SX2X2(Yw2, Yw3, Yw2, Yw3);
		CONJC_SX2X2(Yw4, Yw5, Yw4, Yw5);
		CONJC_SX2X2(Yw6, Yw7, Yw6, Yw7);
		CONJC_SX2X2(Yw8, Yw9, Yw8, Yw9);

		CONST_SX2X2(Yw0t, Yw1t, 0);
		CONST_SX2X2(Yw2t, Yw3t, 0);
		CONST_SX2X2(Yw4t, Yw5t, 0);
		CONST_SX2X2(Yw6t, Yw7t, 0);
		CONST_SX2X2(Yw8t, Yw9t, 0);
		for (m = 0; m < n; m++)
		{
			AE_LSX2X2_IP(Yr0, Yr1, pY, 2 * SZ_CF32);
			AE_LSX2X2_IP(Yr2, Yr3, pY, 2 * SZ_CF32);
			AE_LSX2X2_IP(Yr4, Yr5, pY, 2 * SZ_CF32);
			AE_LSX2X2_IP(Yr6, Yr7, pY, 2 * SZ_CF32);
			AE_LSX2X2_IP(Yr8, Yr9, pY, 2 * SZ_CF32);
			XT_LSX2IP(R0, pR, SZ_CF32);

			MADDMUX_SX2X2(Yw0, Yw1, R0, R0, Yr0, Yr1, 2);
			MADDMUX_SX2X2(Yw0t, Yw1t, R0, R0, Yr0, Yr1, 1);
			MADDMUX_SX2X2(Yw2, Yw3, R0, R0, Yr2, Yr3, 2);
			MADDMUX_SX2X2(Yw2t, Yw3t, R0, R0, Yr2, Yr3, 1);
			MADDMUX_SX2X2(Yw4, Yw5, R0, R0, Yr4, Yr5, 2);
			MADDMUX_SX2X2(Yw4t, Yw5t, R0, R0, Yr4, Yr5, 1);
			MADDMUX_SX2X2(Yw6, Yw7, R0, R0, Yr6, Yr7, 2);
			MADDMUX_SX2X2(Yw6t, Yw7t, R0, R0, Yr6, Yr7, 1);
			MADDMUX_SX2X2(Yw8, Yw9, R0, R0, Yr8, Yr9, 2);
			MADDMUX_SX2X2(Yw8t, Yw9t, R0, R0, Yr8, Yr9, 1);
		}
		ADD_SX2X2(Yw0, Yw1, Yw0, Yw1, Yw0t, Yw1t);
		ADD_SX2X2(Yw2, Yw3, Yw2, Yw3, Yw2t, Yw3t);
		ADD_SX2X2(Yw4, Yw5, Yw4, Yw5, Yw4t, Yw5t);
		ADD_SX2X2(Yw6, Yw7, Yw6, Yw7, Yw6t, Yw7t);
		ADD_SX2X2(Yw8, Yw9, Yw8, Yw9, Yw8t, Yw9t);

		XT_LSX2IP(D0, pD, SZ_CF32);
		MUL_SX2X2(Yw0, Yw1, Yw0, Yw1, D0, D0);
		MUL_SX2X2(Yw2, Yw3, Yw2, Yw3, D0, D0);
		MUL_SX2X2(Yw4, Yw5, Yw4, Yw5, D0, D0);
		MUL_SX2X2(Yw6, Yw7, Yw6, Yw7, D0, D0);
		MUL_SX2X2(Yw8, Yw9, Yw8, Yw9, D0, D0);

		AE_SSX2X2_IP(Yw0, Yw1, pY, 2 * SZ_CF32);
		AE_SSX2X2_IP(Yw2, Yw3, pY, 2 * SZ_CF32);
		AE_SSX2X2_IP(Yw4, Yw5, pY, 2 * SZ_CF32);
		AE_SSX2X2_IP(Yw6, Yw7, pY, 2 * SZ_CF32);
		AE_SSX2X2_IP(Yw8, Yw9, pY, 2 * SZ_CF32);
		pR += 1;
	}
#endif
}
/*
backward recursion: P!=1
does not require tranformed R, reverse inner loop
*/
static void bkw10x10f(
    complex_float* restrict x,
    const complex_float* restrict R,
    const complex_float* restrict D,
    const complex_float* restrict y)
{
	int m, k;
	xtfloatx4 * restrict pX;
	const xtfloatx2 * restrict pR;
	const xtfloatx2 * restrict tpR = (xtfloatx2*)(R + 10 * 11 / 2 - 1); //last element in R
	const xtfloatx2 * restrict pD = (xtfloatx2*)(D + (10 - 1));
	const xtfloatx4 * restrict pY = (xtfloatx4*)(y + (10 * 10 - 2));
	xtfloatx2 X0, X1, X2, X3, X4, X5, X6, X7, X8, X9;
	xtfloatx2 Xw0, Xw1, Xw2, Xw3, Xw4, Xw5, Xw6, Xw7, Xw8, Xw9;
	xtfloatx2 Xw0t, Xw1t, Xw2t, Xw3t, Xw4t, Xw5t, Xw6t, Xw7t, Xw8t, Xw9t;
	xtfloatx2 R0, D0;

	for (k = 10 - 1; k >= 0; k--)
	{
		//for (p = 0; p < P; p++)
		{
			pX = (xtfloatx4*)(x + 10 * 10 - 2); // last element in X
			pR = (xtfloatx2*)(tpR--); //points to the end of row k in R
			AE_LSX2X2_XP(Xw8, Xw9, pY, -2 * SZ_CF32);
			AE_LSX2X2_XP(Xw6, Xw7, pY, -2 * SZ_CF32);
			AE_LSX2X2_XP(Xw4, Xw5, pY, -2 * SZ_CF32);
			AE_LSX2X2_XP(Xw2, Xw3, pY, -2 * SZ_CF32);
			AE_LSX2X2_XP(Xw0, Xw1, pY, -2 * SZ_CF32);

			CONST_SX2X2(Xw0t, Xw1t, 0);
			CONST_SX2X2(Xw2t, Xw3t, 0);
			CONST_SX2X2(Xw4t, Xw5t, 0);
			CONST_SX2X2(Xw6t, Xw7t, 0);
			CONST_SX2X2(Xw8t, Xw9t, 0);
			for (m = 0; m < 10 - k - 1; m++)
			{
				AE_LSX2X2_XP(X8, X9, pX, -2 * SZ_CF32);
				AE_LSX2X2_XP(X6, X7, pX, -2 * SZ_CF32);
				AE_LSX2X2_XP(X4, X5, pX, -2 * SZ_CF32);
				AE_LSX2X2_XP(X2, X3, pX, -2 * SZ_CF32);
				AE_LSX2X2_XP(X0, X1, pX, -2 * SZ_CF32);
				XT_LSX2XP(R0, pR, -(9 - m)*SZ_CF32);
				MADDMUX_SX2X2(Xw8, Xw9, X8, X9, R0, R0, 2);
				MADDMUX_SX2X2(Xw8t, Xw9t, X8, X9, R0, R0, 3);
				MADDMUX_SX2X2(Xw6, Xw7, X6, X7, R0, R0, 2);
				MADDMUX_SX2X2(Xw6t, Xw7t, X6, X7, R0, R0, 3);
				MADDMUX_SX2X2(Xw4, Xw5, X4, X5, R0, R0, 2);
				MADDMUX_SX2X2(Xw4t, Xw5t, X4, X5, R0, R0, 3);
				MADDMUX_SX2X2(Xw2, Xw3, X2, X3, R0, R0, 2);
				MADDMUX_SX2X2(Xw2t, Xw3t, X2, X3, R0, R0, 3);
				MADDMUX_SX2X2(Xw0, Xw1, X0, X1, R0, R0, 2);
				MADDMUX_SX2X2(Xw0t, Xw1t, X0, X1, R0, R0, 3);
			}
			ADD_SX2X2(Xw0, Xw1, Xw0, Xw1, Xw0t, Xw1t);
			ADD_SX2X2(Xw2, Xw3, Xw2, Xw3, Xw2t, Xw3t);
			ADD_SX2X2(Xw4, Xw5, Xw4, Xw5, Xw4t, Xw5t);
			ADD_SX2X2(Xw6, Xw7, Xw6, Xw7, Xw6t, Xw7t);
			ADD_SX2X2(Xw8, Xw9, Xw8, Xw9, Xw8t, Xw9t);

			XT_LSX2XP(D0, pD, -SZ_CF32);
			MUL_SX2X2(Xw0, Xw1, Xw0, Xw1, D0, D0);
			MUL_SX2X2(Xw2, Xw3, Xw2, Xw3, D0, D0);
			MUL_SX2X2(Xw4, Xw5, Xw4, Xw5, D0, D0);
			MUL_SX2X2(Xw6, Xw7, Xw6, Xw7, D0, D0);
			MUL_SX2X2(Xw8, Xw9, Xw8, Xw9, D0, D0);

			AE_SSX2X2_XP(Xw8, Xw9, pX, -2 * SZ_CF32);
			AE_SSX2X2_XP(Xw6, Xw7, pX, -2 * SZ_CF32);
			AE_SSX2X2_XP(Xw4, Xw5, pX, -2 * SZ_CF32);
			AE_SSX2X2_XP(Xw2, Xw3, pX, -2 * SZ_CF32);
			AE_SSX2X2_XP(Xw0, Xw1, pX, -2 * SZ_CF32);
		}
	}
}
/*-------------------------------------------------------------------------
Matrix (Pseudo) Inversion
Obtain Left Inverse of a matrix using Cholesky Decomposition
The result is matrix x = A^-1
Fixed point API requires explicit setting of fixed point representation of 
input/output matrices as well as for internal temporary matrices such as R 
(Cholesky decomposition) and y (decision of R'y=(A'*B))


Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
A[M*N]      matrix A, for fixed point API, the representation is Q(qA)
sigma2      regularization term, for fixed point API, the 
            representation is Q(2*qA-30)
qRA         qR-qA; difference between fixed point representations of R
            and A matrices (for the fixed point API only). Should be 
            equal or less than 0 (typically -2).
qYBRA       qY-qB+qR-qA, combination of fixed point representations of 
            matrices y, B, R and A (for the fixed point API only). Since 
            for matrix inversion we simply use identity matrix B, we may 
            always suppose qB=31 
qXYR        combination of fixed point representation (matrices R, x and y) 
            qXYR=qX-qY+qR (for the fixed point API only)
Output:
x[N*M]      Left Inverse of the matrix A, for fixed point API, the 
            representation is Q(qX)
Temporary:
pScr            Scratch memory

N = M = 4, 6, 8, 10

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
16-bytes boundary
---------------------------------------------------------------------------*/
void  cmatcholpseudoinv10x10f(void* pScr,
    complex_float *x,
    const complex_float * A,
    const float32_t  sigma2)
{
    complex_float * restrict D; int SD;
    complex_float * restrict R; int SR;
    complex_float * restrict y; int SY;

    NASSERT(x);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    SD = 10;
    SR = (((10 + 1) * 10) >> 1);
    SY = 10 * 10;
    SD = (SD + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SR = (SR + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SY = (SY + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SD = 2 * SD * sizeof(float32_t);
    SR = 2 * SR * sizeof(float32_t);
    SY = 2 * SY * sizeof(float32_t);
    R = (complex_float *)pScr;
    D = (complex_float *)(((uintptr_t)R) + SR);
    y = (complex_float *)(((uintptr_t)D) + SD);
    pScr = (complex_float *)(((uintptr_t)y) + SY);

    cmatcholdecomp10x10f(pScr, R, D, A, sigma2);
    fwd10x10f(y, R, D, A);
    bkw10x10f(x, R, D, y);
}

size_t  cmatcholpseudoinv10x10f_getScratchSize()
{
    size_t s_dc;
    int SD, SR, SY;

    SD = 10;
    SR = (((10 + 1) * 10) >> 1);
    SY = 10 * 10;
    SD = (SD + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SR = (SR + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SY = (SY + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SD = 2 * SD * sizeof(float32_t);
    SR = 2 * SR * sizeof(float32_t);
    SY = 2 * SY * sizeof(float32_t);

    s_dc = cmatcholdecomp10x10f_getScratchSize();
    return SD + SR + SY + +s_dc;
}
#elif (HAVE_FPU)
#define SZ_F32 (int)(sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
#define CF32_WIDTH (HIFI_SIMD_WIDTH/SZ_CF32)
/*
Forward recursion for pseudo inversion
Input:
R[((N+1)*N)/2]
upper triangular matrix R
Z[N*P]      matrix A
D[N]        reciprocals of main diagonal
Output:
y[N*P]		Decision matrix y
*/
static void fwd10x10f(
    complex_float* restrict y,
    const complex_float* restrict R,
    const complex_float* restrict D,
    const complex_float* restrict Z)
{
#if 0
    int n, m, p;
    xtfloat Yw_re, Yw_im, Yr_re, Yr_im, R_re, R_im, D_re, D_im;
    xtfloat * restrict pY;
    const xtfloat * restrict pR;
    const xtfloat * restrict pD;
    const xtfloat * restrict pA;
    for (n = 0; n<10; n++)
    {
        //pA = (const xtfloat *)(Z + n);
        for (p = 0; p < 10; p++)
        {
            pR = (const xtfloat *)(R + (n*(n + 1)) / 2);
            pA = (const xtfloat *)(Z + n + p * 10);
            XT_LSIP(Yw_re, pA, SZ_F32);
            XT_LSIP(Yw_im, pA, 10 * SZ_CF32 - SZ_F32);
            Yw_im = XT_NEG_S(Yw_im);
            pY = (xtfloat*)(y + p);
            for (m = 0; m < n; m++)
            {
                XT_LSXP(Yr_re, pY, SZ_F32);
                XT_LSXP(Yr_im, pY, 10 * SZ_CF32 - SZ_F32);
                XT_LSIP(R_re, pR, SZ_F32);
                XT_LSIP(R_im, pR, SZ_F32);
                XT_MSUB_S(Yw_re, R_re, Yr_re);
                XT_MSUB_S(Yw_re, R_im, Yr_im);
                XT_MSUB_S(Yw_im, R_re, Yr_im);
                XT_MADD_S(Yw_im, R_im, Yr_re);

                //XT_MADDMUX_S(Yw, R0, Yr, 2);
                //XT_MADDMUX_S(Yw, R0, Yr, 1);
            }
            pD = (const xtfloat *)(D + n);
            XT_LSIP(D_re, pD, SZ_F32);
            XT_LSIP(D_im, pD, SZ_F32);
            Yw_re = XT_MUL_S(Yw_re, D_re);
            Yw_im = XT_MUL_S(Yw_im, D_im);
            XT_SSIP(Yw_re, pY, SZ_F32);
            XT_SSIP(Yw_im, pY, SZ_F32);
        }
    }
#else
    int n, m;
    xtfloat R_re, R_im, D_re, D_im;
    xtfloat Yw0_re, Yw1_re, Yw2_re, Yw3_re, Yw4_re;
    xtfloat Yw0_im, Yw1_im, Yw2_im, Yw3_im, Yw4_im;
    xtfloat Yr0_re, Yr1_re, Yr2_re, Yr3_re, Yr4_re;
    xtfloat Yr0_im, Yr1_im, Yr2_im, Yr3_im, Yr4_im;

    xtfloat * restrict pY;
    const xtfloat * restrict pR = (const xtfloat *)(R);
    const xtfloat * restrict pRt;
    const xtfloat * restrict pD = (const xtfloat *)(D);
    const xtfloat * restrict pA;
    for (n = 0; n<10; n++)
    {
        pRt = pR;
        pA = (const xtfloat *)(Z + n);
        XT_LSIP(Yw0_re, pA, SZ_F32);
        XT_LSIP(Yw0_im, pA, 10 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw1_re, pA, SZ_F32);
        XT_LSIP(Yw1_im, pA, 10 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw2_re, pA, SZ_F32);
        XT_LSIP(Yw2_im, pA, 10 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw3_re, pA, SZ_F32);
        XT_LSIP(Yw3_im, pA, 10 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw4_re, pA, SZ_F32);
        XT_LSIP(Yw4_im, pA, 10 * SZ_CF32 - SZ_F32);
        Yw0_im = XT_NEG_S(Yw0_im);
        Yw1_im = XT_NEG_S(Yw1_im);
        Yw2_im = XT_NEG_S(Yw2_im);
        Yw3_im = XT_NEG_S(Yw3_im);
        Yw4_im = XT_NEG_S(Yw4_im);
        pY = (xtfloat*)(y);
        for (m = 0; m<n; m++)
        {
            XT_LSXP(Yr0_re, pY, SZ_F32);
            XT_LSXP(Yr0_im, pY, SZ_F32);
            XT_LSXP(Yr1_re, pY, SZ_F32);
            XT_LSXP(Yr1_im, pY, SZ_F32);
            XT_LSXP(Yr2_re, pY, SZ_F32);
            XT_LSXP(Yr2_im, pY, SZ_F32);
            XT_LSXP(Yr3_re, pY, SZ_F32);
            XT_LSXP(Yr3_im, pY, SZ_F32);
            XT_LSXP(Yr4_re, pY, SZ_F32);
            XT_LSXP(Yr4_im, pY, SZ_F32 + 5 * SZ_CF32);
            XT_LSIP(R_re, pR, SZ_F32);
            XT_LSIP(R_im, pR, SZ_F32);
            XT_MSUB_S(Yw0_re, R_re, Yr0_re);
            XT_MSUB_S(Yw0_re, R_im, Yr0_im);
            XT_MSUB_S(Yw0_im, R_re, Yr0_im);
            XT_MADD_S(Yw0_im, R_im, Yr0_re);
            XT_MSUB_S(Yw1_re, R_re, Yr1_re);
            XT_MSUB_S(Yw1_re, R_im, Yr1_im);
            XT_MSUB_S(Yw1_im, R_re, Yr1_im);
            XT_MADD_S(Yw1_im, R_im, Yr1_re);
            XT_MSUB_S(Yw2_re, R_re, Yr2_re);
            XT_MSUB_S(Yw2_re, R_im, Yr2_im);
            XT_MSUB_S(Yw2_im, R_re, Yr2_im);
            XT_MADD_S(Yw2_im, R_im, Yr2_re);
            XT_MSUB_S(Yw3_re, R_re, Yr3_re);
            XT_MSUB_S(Yw3_re, R_im, Yr3_im);
            XT_MSUB_S(Yw3_im, R_re, Yr3_im);
            XT_MADD_S(Yw3_im, R_im, Yr3_re);
            XT_MSUB_S(Yw4_re, R_re, Yr4_re);
            XT_MSUB_S(Yw4_re, R_im, Yr4_im);
            XT_MSUB_S(Yw4_im, R_re, Yr4_im);
            XT_MADD_S(Yw4_im, R_im, Yr4_re);
        }
        //pD = (const xtfloat *)(D + n);
        D_re = XT_LSI(pD, 0);
        D_im = XT_LSI(pD, SZ_F32);
        Yw0_re = XT_MUL_S(Yw0_re, D_re);
        Yw0_im = XT_MUL_S(Yw0_im, D_im);
        Yw1_re = XT_MUL_S(Yw1_re, D_re);
        Yw1_im = XT_MUL_S(Yw1_im, D_im);
        Yw2_re = XT_MUL_S(Yw2_re, D_re);
        Yw2_im = XT_MUL_S(Yw2_im, D_im);
        Yw3_re = XT_MUL_S(Yw3_re, D_re);
        Yw3_im = XT_MUL_S(Yw3_im, D_im);
        Yw4_re = XT_MUL_S(Yw4_re, D_re);
        Yw4_im = XT_MUL_S(Yw4_im, D_im);
        XT_SSIP(Yw0_re, pY, SZ_F32);
        XT_SSIP(Yw0_im, pY, SZ_F32);
        XT_SSIP(Yw1_re, pY, SZ_F32);
        XT_SSIP(Yw1_im, pY, SZ_F32);
        XT_SSIP(Yw2_re, pY, SZ_F32);
        XT_SSIP(Yw2_im, pY, SZ_F32);
        XT_SSIP(Yw3_re, pY, SZ_F32);
        XT_SSIP(Yw3_im, pY, SZ_F32);
        XT_SSIP(Yw4_re, pY, SZ_F32);
        XT_SSIP(Yw4_im, pY, SZ_F32);

        pR = pRt;
        XT_LSIP(Yw0_re, pA, SZ_F32);
        XT_LSIP(Yw0_im, pA, 10 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw1_re, pA, SZ_F32);
        XT_LSIP(Yw1_im, pA, 10 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw2_re, pA, SZ_F32);
        XT_LSIP(Yw2_im, pA, 10 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw3_re, pA, SZ_F32);
        XT_LSIP(Yw3_im, pA, 10 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw4_re, pA, SZ_F32);
        XT_LSIP(Yw4_im, pA, 10 * SZ_CF32 - SZ_F32);
        Yw0_im = XT_NEG_S(Yw0_im);
        Yw1_im = XT_NEG_S(Yw1_im);
        Yw2_im = XT_NEG_S(Yw2_im);
        Yw3_im = XT_NEG_S(Yw3_im);
        Yw4_im = XT_NEG_S(Yw4_im);
        pY = (xtfloat*)(y + 5);
        for (m = 0; m<n; m++)
        {
            XT_LSXP(Yr0_re, pY, SZ_F32);
            XT_LSXP(Yr0_im, pY, SZ_F32);
            XT_LSXP(Yr1_re, pY, SZ_F32);
            XT_LSXP(Yr1_im, pY, SZ_F32);
            XT_LSXP(Yr2_re, pY, SZ_F32);
            XT_LSXP(Yr2_im, pY, SZ_F32);
            XT_LSXP(Yr3_re, pY, SZ_F32);
            XT_LSXP(Yr3_im, pY, SZ_F32);
            XT_LSXP(Yr4_re, pY, SZ_F32);
            XT_LSXP(Yr4_im, pY, SZ_F32 + 5 * SZ_CF32);
            XT_LSIP(R_re, pR, SZ_F32);
            XT_LSIP(R_im, pR, SZ_F32);
            XT_MSUB_S(Yw0_re, R_re, Yr0_re);
            XT_MSUB_S(Yw0_re, R_im, Yr0_im);
            XT_MSUB_S(Yw0_im, R_re, Yr0_im);
            XT_MADD_S(Yw0_im, R_im, Yr0_re);
            XT_MSUB_S(Yw1_re, R_re, Yr1_re);
            XT_MSUB_S(Yw1_re, R_im, Yr1_im);
            XT_MSUB_S(Yw1_im, R_re, Yr1_im);
            XT_MADD_S(Yw1_im, R_im, Yr1_re);
            XT_MSUB_S(Yw2_re, R_re, Yr2_re);
            XT_MSUB_S(Yw2_re, R_im, Yr2_im);
            XT_MSUB_S(Yw2_im, R_re, Yr2_im);
            XT_MADD_S(Yw2_im, R_im, Yr2_re);
            XT_MSUB_S(Yw3_re, R_re, Yr3_re);
            XT_MSUB_S(Yw3_re, R_im, Yr3_im);
            XT_MSUB_S(Yw3_im, R_re, Yr3_im);
            XT_MADD_S(Yw3_im, R_im, Yr3_re);
            XT_MSUB_S(Yw4_re, R_re, Yr4_re);
            XT_MSUB_S(Yw4_re, R_im, Yr4_im);
            XT_MSUB_S(Yw4_im, R_re, Yr4_im);
            XT_MADD_S(Yw4_im, R_im, Yr4_re);
        }
        //pD = (const xtfloat *)(D + n);
        XT_LSIP(D_re, pD, SZ_F32);
        XT_LSIP(D_im, pD, SZ_F32);
        Yw0_re = XT_MUL_S(Yw0_re, D_re);
        Yw0_im = XT_MUL_S(Yw0_im, D_im);
        Yw1_re = XT_MUL_S(Yw1_re, D_re);
        Yw1_im = XT_MUL_S(Yw1_im, D_im);
        Yw2_re = XT_MUL_S(Yw2_re, D_re);
        Yw2_im = XT_MUL_S(Yw2_im, D_im);
        Yw3_re = XT_MUL_S(Yw3_re, D_re);
        Yw3_im = XT_MUL_S(Yw3_im, D_im);
        Yw4_re = XT_MUL_S(Yw4_re, D_re);
        Yw4_im = XT_MUL_S(Yw4_im, D_im);
        XT_SSIP(Yw0_re, pY, SZ_F32);
        XT_SSIP(Yw0_im, pY, SZ_F32);
        XT_SSIP(Yw1_re, pY, SZ_F32);
        XT_SSIP(Yw1_im, pY, SZ_F32);
        XT_SSIP(Yw2_re, pY, SZ_F32);
        XT_SSIP(Yw2_im, pY, SZ_F32);
        XT_SSIP(Yw3_re, pY, SZ_F32);
        XT_SSIP(Yw3_im, pY, SZ_F32);
        XT_SSIP(Yw4_re, pY, SZ_F32);
        XT_SSIP(Yw4_im, pY, SZ_F32);
        pR += 2;
    }
#endif
}
/*
backward recursion: P!=1
does not require tranformed R, reverse inner loop
*/
static void bkw10x10f(
    complex_float* restrict x,
    const complex_float* restrict R,
    const complex_float* restrict D,
    const complex_float* restrict y)
{
    int m, k;
    xtfloat * restrict pX;
    const xtfloat * restrict pR;
    const xtfloat * restrict tpR = (xtfloat*)(R + 10 * (10 + 1) / 2 - 1) + 1; //last element in R
    const xtfloat * restrict pD = (xtfloat*)(D + (10 - 1)) + 1;
    const xtfloat * restrict pY = (xtfloat*)(y + (10 * 10 - 1)) + 1;
    xtfloat X0_re, X1_re, X2_re, X3_re, X4_re;
    xtfloat X0_im, X1_im, X2_im, X3_im, X4_im;
    xtfloat Xw0_re, Xw1_re, Xw2_re, Xw3_re, Xw4_re;
    xtfloat Xw0_im, Xw1_im, Xw2_im, Xw3_im, Xw4_im;
    xtfloat R0_re, R0_im, D0_re, D0_im;

    for (k = 10 - 1; k >= 0; k--)
    {
        pX = (xtfloat*)(x + 10 * 10 - 1) + 1; // last element in X
        pR = (xtfloat*)tpR; //points to the end of row k in R
        XT_LSXP(Xw4_im, pY, -SZ_F32);
        XT_LSXP(Xw4_re, pY, -SZ_F32);
        XT_LSXP(Xw3_im, pY, -SZ_F32);
        XT_LSXP(Xw3_re, pY, -SZ_F32);
        XT_LSXP(Xw2_im, pY, -SZ_F32);
        XT_LSXP(Xw2_re, pY, -SZ_F32);
        XT_LSXP(Xw1_im, pY, -SZ_F32);
        XT_LSXP(Xw1_re, pY, -SZ_F32);
        XT_LSXP(Xw0_im, pY, -SZ_F32);
        XT_LSXP(Xw0_re, pY, -SZ_F32);

        for (m = 0; m < 10 - k - 1; m++)
        {
            XT_LSXP(X4_im, pX, -SZ_F32);
            XT_LSXP(X4_re, pX, -SZ_F32);
            XT_LSXP(X3_im, pX, -SZ_F32);
            XT_LSXP(X3_re, pX, -SZ_F32);
            XT_LSXP(X2_im, pX, -SZ_F32);
            XT_LSXP(X2_re, pX, -SZ_F32);
            XT_LSXP(X1_im, pX, -SZ_F32);
            XT_LSXP(X1_re, pX, -SZ_F32);
            XT_LSXP(X0_im, pX, -SZ_F32);
            XT_LSXP(X0_re, pX, -SZ_F32 - 5 * SZ_CF32);
            XT_LSXP(R0_im, pR, -SZ_F32);
            XT_LSXP(R0_re, pR, -(10 - 1 - m)*SZ_CF32 + SZ_F32);
            XT_MSUB_S(Xw0_re, X0_re, R0_re);
            XT_MADD_S(Xw0_re, X0_im, R0_im);
            XT_MSUB_S(Xw0_im, X0_re, R0_im);
            XT_MSUB_S(Xw0_im, X0_im, R0_re);
            XT_MSUB_S(Xw1_re, X1_re, R0_re);
            XT_MADD_S(Xw1_re, X1_im, R0_im);
            XT_MSUB_S(Xw1_im, X1_re, R0_im);
            XT_MSUB_S(Xw1_im, X1_im, R0_re);
            XT_MSUB_S(Xw2_re, X2_re, R0_re);
            XT_MADD_S(Xw2_re, X2_im, R0_im);
            XT_MSUB_S(Xw2_im, X2_re, R0_im);
            XT_MSUB_S(Xw2_im, X2_im, R0_re);
            XT_MSUB_S(Xw3_re, X3_re, R0_re);
            XT_MADD_S(Xw3_re, X3_im, R0_im);
            XT_MSUB_S(Xw3_im, X3_re, R0_im);
            XT_MSUB_S(Xw3_im, X3_im, R0_re);
            XT_MSUB_S(Xw4_re, X4_re, R0_re);
            XT_MADD_S(Xw4_re, X4_im, R0_im);
            XT_MSUB_S(Xw4_im, X4_re, R0_im);
            XT_MSUB_S(Xw4_im, X4_im, R0_re);
        }
        XT_LSXP(D0_im, pD, -SZ_F32);
        XT_LSXP(D0_re, pD, -SZ_F32);
        Xw0_re = XT_MUL_S(Xw0_re, D0_re);
        Xw0_im = XT_MUL_S(Xw0_im, D0_im);
        Xw1_re = XT_MUL_S(Xw1_re, D0_re);
        Xw1_im = XT_MUL_S(Xw1_im, D0_im);
        Xw2_re = XT_MUL_S(Xw2_re, D0_re);
        Xw2_im = XT_MUL_S(Xw2_im, D0_im);
        Xw3_re = XT_MUL_S(Xw3_re, D0_re);
        Xw3_im = XT_MUL_S(Xw3_im, D0_im);
        Xw4_re = XT_MUL_S(Xw4_re, D0_re);
        Xw4_im = XT_MUL_S(Xw4_im, D0_im);

        XT_SSXP(Xw4_im, pX, -SZ_F32);
        XT_SSXP(Xw4_re, pX, -SZ_F32);
        XT_SSXP(Xw3_im, pX, -SZ_F32);
        XT_SSXP(Xw3_re, pX, -SZ_F32);
        XT_SSXP(Xw2_im, pX, -SZ_F32);
        XT_SSXP(Xw2_re, pX, -SZ_F32);
        XT_SSXP(Xw1_im, pX, -SZ_F32);
        XT_SSXP(Xw1_re, pX, -SZ_F32);
        XT_SSXP(Xw0_im, pX, -SZ_F32);
        XT_SSXP(Xw0_re, pX, -SZ_F32);

        pR = (xtfloat*)tpR; //points to the end of row k in R
        XT_LSXP(Xw4_im, pY, -SZ_F32);
        XT_LSXP(Xw4_re, pY, -SZ_F32);
        XT_LSXP(Xw3_im, pY, -SZ_F32);
        XT_LSXP(Xw3_re, pY, -SZ_F32);
        XT_LSXP(Xw2_im, pY, -SZ_F32);
        XT_LSXP(Xw2_re, pY, -SZ_F32);
        XT_LSXP(Xw1_im, pY, -SZ_F32);
        XT_LSXP(Xw1_re, pY, -SZ_F32);
        XT_LSXP(Xw0_im, pY, -SZ_F32);
        XT_LSXP(Xw0_re, pY, -SZ_F32);

        pX = (xtfloat*)(x + 10 * 10 - 5 - 1) + 1;
        for (m = 0; m < 10 - k - 1; m++)
        {
            XT_LSXP(X4_im, pX, -SZ_F32);
            XT_LSXP(X4_re, pX, -SZ_F32);
            XT_LSXP(X3_im, pX, -SZ_F32);
            XT_LSXP(X3_re, pX, -SZ_F32);
            XT_LSXP(X2_im, pX, -SZ_F32);
            XT_LSXP(X2_re, pX, -SZ_F32);
            XT_LSXP(X1_im, pX, -SZ_F32);
            XT_LSXP(X1_re, pX, -SZ_F32);
            XT_LSXP(X0_im, pX, -SZ_F32);
            XT_LSXP(X0_re, pX, -SZ_F32 - 5 * SZ_CF32);
            XT_LSXP(R0_im, pR, -SZ_F32);
            XT_LSXP(R0_re, pR, -(10 - 1 - m)*SZ_CF32 + SZ_F32);
            XT_MSUB_S(Xw0_re, X0_re, R0_re);
            XT_MADD_S(Xw0_re, X0_im, R0_im);
            XT_MSUB_S(Xw0_im, X0_re, R0_im);
            XT_MSUB_S(Xw0_im, X0_im, R0_re);
            XT_MSUB_S(Xw1_re, X1_re, R0_re);
            XT_MADD_S(Xw1_re, X1_im, R0_im);
            XT_MSUB_S(Xw1_im, X1_re, R0_im);
            XT_MSUB_S(Xw1_im, X1_im, R0_re);
            XT_MSUB_S(Xw2_re, X2_re, R0_re);
            XT_MADD_S(Xw2_re, X2_im, R0_im);
            XT_MSUB_S(Xw2_im, X2_re, R0_im);
            XT_MSUB_S(Xw2_im, X2_im, R0_re);
            XT_MSUB_S(Xw3_re, X3_re, R0_re);
            XT_MADD_S(Xw3_re, X3_im, R0_im);
            XT_MSUB_S(Xw3_im, X3_re, R0_im);
            XT_MSUB_S(Xw3_im, X3_im, R0_re);
            XT_MSUB_S(Xw4_re, X4_re, R0_re);
            XT_MADD_S(Xw4_re, X4_im, R0_im);
            XT_MSUB_S(Xw4_im, X4_re, R0_im);
            XT_MSUB_S(Xw4_im, X4_im, R0_re);
        }
        Xw0_re = XT_MUL_S(Xw0_re, D0_re);
        Xw0_im = XT_MUL_S(Xw0_im, D0_im);
        Xw1_re = XT_MUL_S(Xw1_re, D0_re);
        Xw1_im = XT_MUL_S(Xw1_im, D0_im);
        Xw2_re = XT_MUL_S(Xw2_re, D0_re);
        Xw2_im = XT_MUL_S(Xw2_im, D0_im);
        Xw3_re = XT_MUL_S(Xw3_re, D0_re);
        Xw3_im = XT_MUL_S(Xw3_im, D0_im);
        Xw4_re = XT_MUL_S(Xw4_re, D0_re);
        Xw4_im = XT_MUL_S(Xw4_im, D0_im);

        XT_SSXP(Xw4_im, pX, -SZ_F32);
        XT_SSXP(Xw4_re, pX, -SZ_F32);
        XT_SSXP(Xw3_im, pX, -SZ_F32);
        XT_SSXP(Xw3_re, pX, -SZ_F32);
        XT_SSXP(Xw2_im, pX, -SZ_F32);
        XT_SSXP(Xw2_re, pX, -SZ_F32);
        XT_SSXP(Xw1_im, pX, -SZ_F32);
        XT_SSXP(Xw1_re, pX, -SZ_F32);
        XT_SSXP(Xw0_im, pX, -SZ_F32);
        XT_SSXP(Xw0_re, pX, -SZ_F32);
        tpR -= 2;
    }
}
/*-------------------------------------------------------------------------
Matrix (Pseudo) Inversion
Obtain Left Inverse of a matrix using Cholesky Decomposition
The result is matrix x = A^-1
Fixed point API requires explicit setting of fixed point representation of
input/output matrices as well as for internal temporary matrices such as R
(Cholesky decomposition) and y (decision of R'y=(A'*B))


Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
A[M*N]      matrix A, for fixed point API, the representation is Q(qA)
sigma2      regularization term, for fixed point API, the
representation is Q(2*qA-30)
qRA         qR-qA; difference between fixed point representations of R
and A matrices (for the fixed point API only). Should be
equal or less than 0 (typically -2).
qYBRA       qY-qB+qR-qA, combination of fixed point representations of
matrices y, B, R and A (for the fixed point API only). Since
for matrix inversion we simply use identity matrix B, we may
always suppose qB=31
qXYR        combination of fixed point representation (matrices R, x and y)
qXYR=qX-qY+qR (for the fixed point API only)
Output:
x[N*M]      Left Inverse of the matrix A, for fixed point API, the
representation is Q(qX)
Temporary:
pScr            Scratch memory

N = M = 4, 6, 8, 10

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
16-bytes boundary
---------------------------------------------------------------------------*/
void  cmatcholpseudoinv10x10f(void* pScr,
    complex_float *x,
    const complex_float * A,
    const float32_t  sigma2)
{
    complex_float * restrict D; int SD;
    complex_float * restrict R; int SR;
    complex_float * restrict y; int SY;

    NASSERT(x);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    SD = 10;
    SR = (((10 + 1) * 10) >> 1);
    SY = 10 * 10;
    SD = (SD + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SR = (SR + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SY = (SY + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SD = 2 * SD * sizeof(float32_t);
    SR = 2 * SR * sizeof(float32_t);
    SY = 2 * SY * sizeof(float32_t);
    R = (complex_float *)pScr;
    D = (complex_float *)(((uintptr_t)R) + SR);
    y = (complex_float *)(((uintptr_t)D) + SD);
    pScr = (complex_float *)(((uintptr_t)y) + SY);

    cmatcholdecomp10x10f(pScr, R, D, A, sigma2);
    fwd10x10f(y, R, D, A);
    bkw10x10f(x, R, D, y);
}

size_t  cmatcholpseudoinv10x10f_getScratchSize()
{
    size_t s_dc;
    int SD, SR, SY;

    SD = 10;
    SR = (((10 + 1) * 10) >> 1);
    SY = 10 * 10;
    SD = (SD + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SR = (SR + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SY = (SY + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SD = 2 * SD * sizeof(float32_t);
    SR = 2 * SR * sizeof(float32_t);
    SY = 2 * SY * sizeof(float32_t);

    s_dc = cmatcholdecomp10x10f_getScratchSize();
    return SD + SR + SY + s_dc;
}
#else
DISCARD_FUN(void, cmatcholpseudoinv10x10f, (void* pScr,
	complex_float *R,
	const complex_float * A,
	const float32_t sigma2))

	size_t  cmatcholpseudoinv10x10f_getScratchSize()
{
	return 0;
}
#endif
