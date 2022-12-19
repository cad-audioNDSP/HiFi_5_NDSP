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
#define F32_WIDTH (HIFI_SIMD_WIDTH/SZ_F32)
/*
Forward recursion
Input:
R[((N+1)*N)/2]
upper triangular matrix R
Z[N*P]      matrix A'*B
D[N]        reciprocals of main diagonal
Output:
y[N*P]		Decision matrix y
Temporary:
pScr        Scratch memory
*/
static void rfwd6x6f(
    float32_t* restrict y,
    const float32_t* restrict R,
    const float32_t* restrict D,
    const float32_t* restrict Z)
{
	int n, m;
	xtfloat R0, R0t, D0;
	xtfloatx2 Yr01, Yr23, Yr45;
	xtfloatx2 Yw0, Yw1, Yw2, Yw3, Yw4, Yw5;
	xtfloatx2 Yw01, Yw23, Yw45;
	xtfloatx2 Yw01t, Yw23t, Yw45t;
	xtfloatx2 * restrict pY;
	const xtfloat * restrict pR = (const xtfloat *)R;
	const xtfloat * restrict pRt = (const xtfloat *)R;

	const xtfloat * restrict pD = (const xtfloat *)D;
	const xtfloatx2 * restrict pA = (const xtfloatx2 *)Z;
	for (n = 0; n<6; n += 2)
	{
		pY = (xtfloatx2*)y;
		XT_LSX2IP(Yw0, pA, 6 * SZ_F32);
		XT_LSX2IP(Yw1, pA, 6 * SZ_F32);
		XT_LSX2IP(Yw2, pA, 6 * SZ_F32);
		XT_LSX2IP(Yw3, pA, 6 * SZ_F32);
		XT_LSX2IP(Yw4, pA, 6 * SZ_F32);
		XT_LSX2XP(Yw5, pA, (-6 * 5 + 2)*SZ_F32);

		Yw01 = XT_SEL32_HH_SX2(Yw0, Yw1);
		Yw23 = XT_SEL32_HH_SX2(Yw2, Yw3);
		Yw45 = XT_SEL32_HH_SX2(Yw4, Yw5);
		Yw01t = XT_SEL32_LL_SX2(Yw0, Yw1);
		Yw23t = XT_SEL32_LL_SX2(Yw2, Yw3);
		Yw45t = XT_SEL32_LL_SX2(Yw4, Yw5);
		pRt = pR + n + 1;
		for (m = 0; m<n; m++)
		{
			AE_LSX2IP(Yr01, pY, 2 * SZ_F32);
			AE_LSX2IP(Yr23, pY, 2 * SZ_F32);
			AE_LSX2IP(Yr45, pY, 2 * SZ_F32);

			XT_LSIP(R0, pR, SZ_F32);
			XT_LSIP(R0t, pRt, SZ_F32);
			MSUB_SX2X2(Yw01, Yw23, R0, R0, Yr01, Yr23);
			MSUB_SX2X2(Yw01t, Yw23t, R0t, R0t, Yr01, Yr23);
			MSUB_SX2X2(Yw45, Yw45t, R0, R0t, Yr45, Yr45);

		}

		XT_LSIP(D0, pD, SZ_F32);
		MUL_SX2X2(Yw01, Yw23, Yw01, Yw23, D0, D0);
		Yw45 = MUL_SX2(Yw45, D0);

		Yr01 = Yw01; Yr23 = Yw23; Yr45 = Yw45;

		XT_LSIP(R0t, pRt, SZ_F32);
		MSUB_SX2X2(Yw01t, Yw23t, R0t, R0t, Yr01, Yr23);
		MSUB_SX2(Yw45t, R0t, Yr45);

		XT_LSIP(D0, pD, SZ_F32);
		MUL_SX2X2(Yw01t, Yw23t, Yw01t, Yw23t, D0, D0);
		Yw45t = MUL_SX2(Yw45t, D0);

		AE_SSX2X2_IP(Yw01, Yw23, castxcc(xtfloatx4, pY), 4 * SZ_F32);
		AE_SSX2X2_IP(Yw45, Yw01t, castxcc(xtfloatx4, pY), 4 * SZ_F32);
		AE_SSX2X2_IP(Yw23t, Yw45t, castxcc(xtfloatx4, pY), 4 * SZ_F32);
		pR = pRt + 1;
	}
}
/*
backward recursion
does not require tranformed R, reverse inner loop
*/
static void rbkw6x6f(
    float32_t* restrict x,
    const float32_t* restrict R,
    const float32_t* restrict D,
    const float32_t* restrict y)
{
	int m, k;
	xtfloatx2 * restrict pX;
	const xtfloat * restrict pR0;
	const xtfloat * restrict pR1;
	const xtfloat * restrict tpR = (xtfloat*)(R + 6 * 7 / 2 - 1); //last element in R
	const xtfloat * restrict pD = (xtfloat*)(D + (6 - 1));
	const xtfloatx4 * restrict pY = (xtfloatx4*)(y + (6 * 6 - 4));
	xtfloatx2 Xw01, Xw23, Xw45;
	xtfloatx2 Xw01t, Xw23t, Xw45t;
	xtfloatx2 Xr01, Xr23, Xr45;
	xtfloat R0, R0t, D0;

	for (k = 6 - 1; k >= 0; k -= 2)
	{

		pX = (xtfloatx2*)(x + 6 * 6 - 2); // last 2 elements in X
		pR0 = (xtfloat*)(tpR--); //points to the end of row k in R
		pR1 = (xtfloat*)(tpR--); //points to the end of row k in R

		AE_LSX2X2_IP(Xw23, Xw45, pY, -4 * SZ_F32);
		AE_LSX2X2_IP(Xw45t, Xw01, pY, -4 * SZ_F32);
		AE_LSX2X2_IP(Xw01t, Xw23t, pY, -4 * SZ_F32);

		for (m = 0; m < 6 - k - 1; m++)
		{
			AE_LSX2XP(Xr45, pX, -2 * SZ_F32);
			AE_LSX2XP(Xr23, pX, -2 * SZ_F32);
			AE_LSX2XP(Xr01, pX, -2 * SZ_F32);

			XT_LSXP(R0, pR0, -(5 - m)*SZ_F32);
			XT_LSXP(R0t, pR1, -(5 - m)*SZ_F32);
			MSUB_SX2X2(Xw01, Xw23, Xr01, Xr23, R0, R0);
			MSUB_SX2X2(Xw01t, Xw23t, Xr01, Xr23, R0t, R0t);
			MSUB_SX2X2(Xw45, Xw45t, Xr45, Xr45, R0, R0t);
		}
		XT_LSXP(D0, pD, -SZ_F32);
		MUL_SX2X2(Xw01, Xw23, Xw01, Xw23, D0, D0);
		Xw45 = XT_MUL_SX2(Xw45, D0);

		Xr01 = Xw01; Xr23 = Xw23; Xr45 = Xw45;

		XT_LSXP(R0, pR1, -k * SZ_F32);
		MSUB_SX2X2(Xw01t, Xw23t, Xr01, Xr23, R0, R0);
		MSUB_SX2(Xw45t, Xr45, R0);

		XT_LSXP(D0, pD, -SZ_F32);
		MUL_SX2X2(Xw01t, Xw23t, Xw01t, Xw23t, D0, D0);
		Xw45t = XT_MUL_SX2(Xw45t, D0);

		pX--;
		AE_SSX2X2_IP(Xw23, Xw45, castxcc(xtfloatx4, pX), -4 * SZ_F32);
		AE_SSX2X2_IP(Xw45t, Xw01, castxcc(xtfloatx4, pX), -4 * SZ_F32);
		AE_SSX2X2_IP(Xw01t, Xw23t, castxcc(xtfloatx4, pX), -4 * SZ_F32);

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
void  matcholpseudoinv6x6f(void* pScr,
    float32_t *x,
    const float32_t * A,
    const float32_t  sigma2)
{
    float32_t * restrict D; int SD;
    float32_t * restrict R; int SR;
    float32_t * restrict y; int SY;

    NASSERT(x);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    SD = 6;
    SR = (((6 + 1) * 6) >> 1);
    SY = 6 * 6;
    SD = (SD + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SR = (SR + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SY = (SY + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SD = SD * sizeof(float32_t);
    SR = SR * sizeof(float32_t);
    SY = SY * sizeof(float32_t);
    R = (float32_t *)pScr;
    D = (float32_t *)(((uintptr_t)R) + SR);
    y = (float32_t *)(((uintptr_t)D) + SD);
    pScr = (float32_t *)(((uintptr_t)y) + SY);
    matcholdecomp6x6f(pScr, R, D, A, sigma2);
    rfwd6x6f(y, R, D, A);
    rbkw6x6f(x, R, D, y);
}

size_t   matcholpseudoinv6x6f_getScratchSize()
{
    size_t s_dc;
    int SD, SR, SY;

    SD = 6;
    SR = (((6 + 1) * 6) >> 1);
    SY = 6 * 6;
    SD = (SD + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SR = (SR + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SY = (SY + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SD = SD * sizeof(float32_t);
    SR = SR * sizeof(float32_t);
    SY = SY * sizeof(float32_t);

    s_dc = matcholdecomp6x6f_getScratchSize();

    return SD + SR + SY + s_dc;
}
#elif (HAVE_FPU)
#define SZ_F32 (int)(sizeof(float32_t))
#define F32_WIDTH (HIFI_SIMD_WIDTH/SZ_F32)
/*
Forward recursion
Input:
R[((N+1)*N)/2]
upper triangular matrix R
Z[N*P]      matrix A'*B
D[N]        reciprocals of main diagonal
Output:
y[N*P]		Decision matrix y
Temporary:
pScr        Scratch memory
*/
static void rfwd6x6f(
    float32_t* restrict y,
    const float32_t* restrict R,
    const float32_t* restrict D,
    const float32_t* restrict Z)
{
#if 0
    int n, m, p;
    xtfloat Yw, Yr, R0, D0;
    xtfloat * restrict pY;
    const xtfloat * restrict pR;
    const xtfloat * restrict pD;
    const xtfloat * restrict pA;
    for (n = 0; n<6; n++)
        for (p = 0; p<6; p++)
        {
            pR = (const xtfloat *)(R + (n*(n + 1)) / 2);
            pA = (const xtfloat *)(Z + n + p * 6);
            XT_LSIP(Yw, pA, SZ_F32);
            pY = (xtfloat*)(y + p);
            for (m = 0; m<n; m++)
            {
                XT_LSXP(Yr, pY, 6 * SZ_F32);
                XT_LSIP(R0, pR, SZ_F32);
                XT_MSUB_S(Yw, R0, Yr);
            }
            pD = (const xtfloat *)(D + n);
            XT_LSIP(D0, pD, SZ_F32);
            Yw = XT_MUL_S(Yw, D0);
            XT_SSIP(Yw, pY, SZ_F32);
        }
#else
    int n, m;
    xtfloat R0, D0;
    xtfloat Yr0, Yr1, Yr2, Yr3, Yr4, Yr5;
    xtfloat Yw0, Yw1, Yw2, Yw3, Yw4, Yw5;
    xtfloat * restrict pY;
    const xtfloat * restrict pR = (const xtfloat *)R;
    const xtfloat * restrict pD = (const xtfloat *)D;
    const xtfloat * restrict pA = (const xtfloat *)Z;
    for (n = 0; n<6; n++)
    {
        pY = (xtfloat*)y;
        XT_LSIP(Yw0, pA, 6 * SZ_F32);
        XT_LSIP(Yw1, pA, 6 * SZ_F32);
        XT_LSIP(Yw2, pA, 6 * SZ_F32);
        XT_LSIP(Yw3, pA, 6 * SZ_F32);
        XT_LSIP(Yw4, pA, 6 * SZ_F32);
        XT_LSXP(Yw5, pA, (-6 * (6 - 1) + 1)*SZ_F32);

        for (m = 0; m<n; m++)
        {
            XT_LSIP(Yr0, pY, SZ_F32);
            XT_LSIP(Yr1, pY, SZ_F32);
            XT_LSIP(Yr2, pY, SZ_F32);
            XT_LSIP(Yr3, pY, SZ_F32);
            XT_LSIP(Yr4, pY, SZ_F32);
            XT_LSIP(Yr5, pY, SZ_F32);
            XT_LSIP(R0, pR, SZ_F32);
            XT_MSUB_S(Yw0, R0, Yr0);
            XT_MSUB_S(Yw1, R0, Yr1);
            XT_MSUB_S(Yw2, R0, Yr2);
            XT_MSUB_S(Yw3, R0, Yr3);
            XT_MSUB_S(Yw4, R0, Yr4);
            XT_MSUB_S(Yw5, R0, Yr5);
        }

        XT_LSIP(D0, pD, SZ_F32);
        Yw0 = XT_MUL_S(Yw0, D0);
        Yw1 = XT_MUL_S(Yw1, D0);
        Yw2 = XT_MUL_S(Yw2, D0);
        Yw3 = XT_MUL_S(Yw3, D0);
        Yw4 = XT_MUL_S(Yw4, D0);
        Yw5 = XT_MUL_S(Yw5, D0);

        XT_SSIP(Yw0, pY, SZ_F32);
        XT_SSIP(Yw1, pY, SZ_F32);
        XT_SSIP(Yw2, pY, SZ_F32);
        XT_SSIP(Yw3, pY, SZ_F32);
        XT_SSIP(Yw4, pY, SZ_F32);
        XT_SSIP(Yw5, pY, SZ_F32);
        pR += 1;
    }
#endif
}
/*
backward recursion: P!=1
*/
static void rbkw6x6f(
    float32_t* restrict x,
    const float32_t* restrict R,
    const float32_t* restrict D,
    const float32_t* restrict y)
{
    int m, k;
    xtfloat * restrict pX;
    const xtfloat * restrict pR;
    const xtfloat * restrict tpR = (xtfloat*)(R + 6 * (6 + 1) / 2 - 1); //last element in R
    const xtfloat * restrict pD = (xtfloat*)(D + (6 - 1));
    const xtfloat * restrict pY = (xtfloat*)(y + (6 * 6 - 1));
    xtfloat Xw0, Xw1, Xw2, Xw3, Xw4, Xw5;
    xtfloat Xr0, Xr1, Xr2, Xr3, Xr4, Xr5;
    xtfloat R0, D0;

    for (k = 6 - 1; k >= 0; k--)
    {
        {
            pX = (xtfloat*)(x + 6 * 6 - 1); // last element in X
            pR = (xtfloat*)(tpR--); //points to the end of row k in R
            XT_LSXP(Xw5, pY, -SZ_F32);
            XT_LSXP(Xw4, pY, -SZ_F32);
            XT_LSXP(Xw3, pY, -SZ_F32);
            XT_LSXP(Xw2, pY, -SZ_F32);
            XT_LSXP(Xw1, pY, -SZ_F32);
            XT_LSXP(Xw0, pY, -SZ_F32);

            for (m = 0; m < 6 - k - 1; m++)
            {
                XT_LSXP(Xr5, pX, -SZ_F32);
                XT_LSXP(Xr4, pX, -SZ_F32);
                XT_LSXP(Xr3, pX, -SZ_F32);
                XT_LSXP(Xr2, pX, -SZ_F32);
                XT_LSXP(Xr1, pX, -SZ_F32);
                XT_LSXP(Xr0, pX, -SZ_F32);

                XT_LSXP(R0, pR, -(6 - 1 - m)*SZ_F32);
                XT_MSUB_S(Xw5, Xr5, R0);
                XT_MSUB_S(Xw4, Xr4, R0);
                XT_MSUB_S(Xw3, Xr3, R0);
                XT_MSUB_S(Xw2, Xr2, R0);
                XT_MSUB_S(Xw1, Xr1, R0);
                XT_MSUB_S(Xw0, Xr0, R0);
            }
            XT_LSXP(D0, pD, -SZ_F32);
            Xw5 = XT_MUL_S(Xw5, D0);
            Xw4 = XT_MUL_S(Xw4, D0);
            Xw3 = XT_MUL_S(Xw3, D0);
            Xw2 = XT_MUL_S(Xw2, D0);
            Xw1 = XT_MUL_S(Xw1, D0);
            Xw0 = XT_MUL_S(Xw0, D0);

            XT_SSXP(Xw5, pX, -SZ_F32);
            XT_SSXP(Xw4, pX, -SZ_F32);
            XT_SSXP(Xw3, pX, -SZ_F32);
            XT_SSXP(Xw2, pX, -SZ_F32);
            XT_SSXP(Xw1, pX, -SZ_F32);
            XT_SSXP(Xw0, pX, -SZ_F32);
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
void  matcholpseudoinv6x6f(void* pScr,
    float32_t *x,
    const float32_t * A,
    const float32_t  sigma2)
{
    float32_t * restrict D; int SD;
    float32_t * restrict R; int SR;
    float32_t * restrict y; int SY;

    NASSERT(x);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    SD = 6;
    SR = (((6 + 1) * 6) >> 1);
    SY = 6 * 6;
    SD = (SD + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SR = (SR + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SY = (SY + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SD = SD * sizeof(float32_t);
    SR = SR * sizeof(float32_t);
    SY = SY * sizeof(float32_t);
    R = (float32_t *)pScr;
    D = (float32_t *)(((uintptr_t)R) + SR);
    y = (float32_t *)(((uintptr_t)D) + SD);
    pScr = (float32_t *)(((uintptr_t)y) + SY);
    matcholdecomp6x6f(pScr, R, D, A, sigma2);
    rfwd6x6f(y, R, D, A);
    rbkw6x6f(x, R, D, y);
}

size_t   matcholpseudoinv6x6f_getScratchSize()
{
    size_t s_dc;
    int SD, SR, SY;

    SD = 6;
    SR = (((6 + 1) * 6) >> 1);
    SY = 6 * 6;
    SD = (SD + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SR = (SR + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SY = (SY + F32_WIDTH - 1)&~(F32_WIDTH - 1);
    SD = SD * sizeof(float32_t);
    SR = SR * sizeof(float32_t);
    SY = SY * sizeof(float32_t);

    s_dc = matcholdecomp6x6f_getScratchSize();

    return SD + SR + SY + s_dc;
}

#else
DISCARD_FUN(void, matcholpseudoinv6x6f, (void* pScr,
	float32_t *R,
	const float32_t * A,
	const float32_t  sigma2))

	size_t  matcholpseudoinv6x6f_getScratchSize()
{
	return 0;
}
#endif
