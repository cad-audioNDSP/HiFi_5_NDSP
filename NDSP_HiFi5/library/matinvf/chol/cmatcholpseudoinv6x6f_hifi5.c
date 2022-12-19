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
static void fwd6x6f(
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
    for (n = 0; n<6; n++)
        for (p = 0; p<6; p++)
        {
            pR = (const xtfloatx2 *)(R + (n*(n + 1))/2);
            pA = (const xtfloatx2 *)(Z + n + p * 6);
            XT_LSX2IP(Yw, pA, SZ_CF32);
            Yw = XT_CONJC_S(Yw);
            pY = (xtfloatx2*)(y + p);
            for (m = 0; m<n; m++)
            {

                XT_LSX2XP(Yr, pY, 6*SZ_CF32);
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
    xtfloatx2 Yr0, Yr1, Yr2, Yr3, Yr4, Yr5;
    xtfloatx2 Yw0, Yw1, Yw2, Yw3, Yw4, Yw5;
    xtfloatx2 Yw0t, Yw1t, Yw2t, Yw3t, Yw4t, Yw5t;
    xtfloatx4 * restrict pY;
    const xtfloatx2 * restrict pR = (const xtfloatx2 *)R;
    const xtfloatx2 * restrict pD = (const xtfloatx2 *)D;
    const xtfloatx2 * restrict pA = (const xtfloatx2 *)Z;
    for (n = 0; n < 6; n++)
        //for (p = 0; p<P; p++)
    {
        pY = (xtfloatx4*)y;
        XT_LSX2IP(Yw0, pA, 6 * SZ_CF32);
        XT_LSX2IP(Yw1, pA, 6 * SZ_CF32);
        XT_LSX2IP(Yw2, pA, 6 * SZ_CF32);
        XT_LSX2IP(Yw3, pA, 6 * SZ_CF32);
        XT_LSX2IP(Yw4, pA, 6 * SZ_CF32);
        XT_LSX2XP(Yw5, pA, (-6 * 5 + 1)*SZ_CF32);
		CONJC_SX2X2(Yw0, Yw1, Yw0, Yw1);
		CONJC_SX2X2(Yw2, Yw3, Yw2, Yw3);
		CONJC_SX2X2(Yw4, Yw5, Yw4, Yw5);

		CONST_SX2X2(Yw0t, Yw1t, 0);
		CONST_SX2X2(Yw2t, Yw3t, 0);
		CONST_SX2X2(Yw4t, Yw5t, 0);
		for (m = 0; m < n; m++)
        {
			AE_LSX2X2_IP(Yr0, Yr1, pY, 2 * SZ_CF32);
			AE_LSX2X2_IP(Yr2, Yr3, pY, 2 * SZ_CF32);
			AE_LSX2X2_IP(Yr4, Yr5, pY, 2 * SZ_CF32);
			XT_LSX2IP(R0, pR, SZ_CF32);

			MADDMUX_SX2X2(Yw0, Yw1, R0, R0, Yr0, Yr1, 2);
			MADDMUX_SX2X2(Yw0t, Yw1t, R0, R0, Yr0, Yr1, 1);
			MADDMUX_SX2X2(Yw2, Yw3, R0, R0, Yr2, Yr3, 2);
			MADDMUX_SX2X2(Yw2t, Yw3t, R0, R0, Yr2, Yr3, 1);
			MADDMUX_SX2X2(Yw4, Yw5, R0, R0, Yr4, Yr5, 2);
			MADDMUX_SX2X2(Yw4t, Yw5t, R0, R0, Yr4, Yr5, 1);
		}
		ADD_SX2X2(Yw0, Yw1, Yw0, Yw1, Yw0t, Yw1t);
		ADD_SX2X2(Yw2, Yw3, Yw2, Yw3, Yw2t, Yw3t);
		ADD_SX2X2(Yw4, Yw5, Yw4, Yw5, Yw4t, Yw5t);

        XT_LSX2IP(D0, pD, SZ_CF32);
		MUL_SX2X2(Yw0, Yw1, Yw0, Yw1, D0, D0);
		MUL_SX2X2(Yw2, Yw3, Yw2, Yw3, D0, D0);
		MUL_SX2X2(Yw4, Yw5, Yw4, Yw5, D0, D0);

		AE_SSX2X2_IP(Yw0, Yw1, pY, 2 * SZ_CF32);
		AE_SSX2X2_IP(Yw2, Yw3, pY, 2 * SZ_CF32);
		AE_SSX2X2_IP(Yw4, Yw5, pY, 2 * SZ_CF32);
		pR += 1;
    }
#endif
}
/*
backward recursion: P!=1
does not require tranformed R, reverse inner loop
*/
static void bkw6x6f(
    complex_float* restrict x,
    const complex_float* restrict R,
    const complex_float* restrict D,
    const complex_float* restrict y)
{
    int m, k;
    xtfloatx4 * restrict pX;
    const xtfloatx2 * restrict pR;
    const xtfloatx2 * restrict tpR = (xtfloatx2*)(R + 6 * 7 / 2 - 1); //last element in R
    const xtfloatx2 * restrict pD = (xtfloatx2*)(D + (6 - 1));
    const xtfloatx4 * restrict pY = (xtfloatx4*)(y + (6 * 6 - 2));
    xtfloatx2 X0, X1, X2, X3, X4, X5;
    xtfloatx2 Xw0, Xw1, Xw2, Xw3, Xw4, Xw5;
    xtfloatx2 Xw0t, Xw1t, Xw2t, Xw3t, Xw4t, Xw5t;
    xtfloatx2 R0, D0;

    for (k = 6 - 1; k >= 0; k--)
    {
        //for (p = 0; p < P; p++)
        {
            pX = (xtfloatx4*)(x + 6 * 6 - 2); // last element in X
            pR = (xtfloatx2*)(tpR--); //points to the end of row k in R
			AE_LSX2X2_XP(Xw4, Xw5, pY, -2 * SZ_CF32);
			AE_LSX2X2_XP(Xw2, Xw3, pY, -2 * SZ_CF32);
			AE_LSX2X2_XP(Xw0, Xw1, pY, -2 * SZ_CF32);

			CONST_SX2X2(Xw0t, Xw1t, 0);
			CONST_SX2X2(Xw2t, Xw3t, 0);
			CONST_SX2X2(Xw4t, Xw5t, 0);
			for (m = 0; m < 6 - k - 1; m++)
            {
				AE_LSX2X2_XP(X4, X5, pX, -2 * SZ_CF32);
				AE_LSX2X2_XP(X2, X3, pX, -2 * SZ_CF32);
				AE_LSX2X2_XP(X0, X1, pX, -2 * SZ_CF32);
				XT_LSX2XP(R0, pR, -(5 - m)*SZ_CF32);
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

            XT_LSX2XP(D0, pD, -SZ_CF32);
			MUL_SX2X2(Xw0, Xw1, Xw0, Xw1, D0, D0);
			MUL_SX2X2(Xw2, Xw3, Xw2, Xw3, D0, D0);
			MUL_SX2X2(Xw4, Xw5, Xw4, Xw5, D0, D0);

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
void  cmatcholpseudoinv6x6f(void* pScr,
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

    SD = 6;
    SR = (((6 + 1) * 6) >> 1);
    SY = 6 * 6;
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

    cmatcholdecomp6x6f(pScr, R, D, A, sigma2);
    fwd6x6f(y, R, D, A);
    bkw6x6f(x, R, D, y);
}

size_t  cmatcholpseudoinv6x6f_getScratchSize()
{
    size_t s_dc;
    int SD, SR, SY;

    SD = 6;
    SR = (((6 + 1) * 6) >> 1);
    SY = 6 * 6;
    SD = (SD + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SR = (SR + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SY = (SY + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SD = 2 * SD * sizeof(float32_t);
    SR = 2 * SR * sizeof(float32_t);
    SY = 2 * SY * sizeof(float32_t);

    s_dc = cmatcholdecomp6x6f_getScratchSize();
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
static void fwd6x6f(
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
    for (n = 0; n<6; n++)
    {
        //pA = (const xtfloat *)(Z + n);
        for (p = 0; p < 6; p++)
        {
            pR = (const xtfloat *)(R + (n*(n + 1)) / 2);
            pA = (const xtfloat *)(Z + n + p * 6);
            XT_LSIP(Yw_re, pA, SZ_F32);
            XT_LSIP(Yw_im, pA, 6 * SZ_CF32 - SZ_F32);
            Yw_im = XT_NEG_S(Yw_im);
            pY = (xtfloat*)(y + p);
            for (m = 0; m < n; m++)
            {
                XT_LSXP(Yr_re, pY, SZ_F32);
                XT_LSXP(Yr_im, pY, 6 * SZ_CF32 - SZ_F32);
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
#elif 1
    int n, m;
    xtfloat R_re, R_im, D_re, D_im;
    xtfloat Yw0_re, Yw1_re, Yw2_re;
    xtfloat Yw0_im, Yw1_im, Yw2_im;
    xtfloat Yr0_re, Yr1_re, Yr2_re;
    xtfloat Yr0_im, Yr1_im, Yr2_im;

    xtfloat * restrict pY;
    const xtfloat * restrict pR = (const xtfloat *)(R);
    const xtfloat * restrict pRt;
    const xtfloat * restrict pD = (const xtfloat *)(D);
    const xtfloat * restrict pA;
    for (n = 0; n<6; n++)
    {
        pRt = pR;
        pA = (const xtfloat *)(Z + n);
        XT_LSIP(Yw0_re, pA, SZ_F32);
        XT_LSIP(Yw0_im, pA, 6 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw1_re, pA, SZ_F32);
        XT_LSIP(Yw1_im, pA, 6 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw2_re, pA, SZ_F32);
        XT_LSIP(Yw2_im, pA, 6 * SZ_CF32 - SZ_F32);
        Yw0_im = XT_NEG_S(Yw0_im);
        Yw1_im = XT_NEG_S(Yw1_im);
        Yw2_im = XT_NEG_S(Yw2_im);
        pY = (xtfloat*)(y);
        for (m = 0; m<n; m++)
        {
            XT_LSXP(Yr0_re, pY, SZ_F32);
            XT_LSXP(Yr0_im, pY, SZ_F32);
            XT_LSXP(Yr1_re, pY, SZ_F32);
            XT_LSXP(Yr1_im, pY, SZ_F32);
            XT_LSXP(Yr2_re, pY, SZ_F32);
            XT_LSXP(Yr2_im, pY, SZ_F32 + 3 * SZ_CF32);
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
        XT_SSIP(Yw0_re, pY, SZ_F32);
        XT_SSIP(Yw0_im, pY, SZ_F32);
        XT_SSIP(Yw1_re, pY, SZ_F32);
        XT_SSIP(Yw1_im, pY, SZ_F32);
        XT_SSIP(Yw2_re, pY, SZ_F32);
        XT_SSIP(Yw2_im, pY, SZ_F32);
        
        pR = pRt;
        XT_LSIP(Yw0_re, pA, SZ_F32);
        XT_LSIP(Yw0_im, pA, 6 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw1_re, pA, SZ_F32);
        XT_LSIP(Yw1_im, pA, 6 * SZ_CF32 - SZ_F32);
        XT_LSIP(Yw2_re, pA, SZ_F32);
        XT_LSIP(Yw2_im, pA, 6 * SZ_CF32 - SZ_F32);
        Yw0_im = XT_NEG_S(Yw0_im);
        Yw1_im = XT_NEG_S(Yw1_im);
        Yw2_im = XT_NEG_S(Yw2_im);
        pY = (xtfloat*)(y + 3);
        for (m = 0; m<n; m++)
        {
            XT_LSXP(Yr0_re, pY, SZ_F32);
            XT_LSXP(Yr0_im, pY, SZ_F32);
            XT_LSXP(Yr1_re, pY, SZ_F32);
            XT_LSXP(Yr1_im, pY, SZ_F32);
            XT_LSXP(Yr2_re, pY, SZ_F32);
            XT_LSXP(Yr2_im, pY, SZ_F32 + 3*SZ_CF32);
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
        XT_SSIP(Yw0_re, pY, SZ_F32);
        XT_SSIP(Yw0_im, pY, SZ_F32);
        XT_SSIP(Yw1_re, pY, SZ_F32);
        XT_SSIP(Yw1_im, pY, SZ_F32);
        XT_SSIP(Yw2_re, pY, SZ_F32);
        XT_SSIP(Yw2_im, pY, SZ_F32);
        pR += 2;
    }
#endif
}
/*
backward recursion: P!=1
does not require tranformed R, reverse inner loop
*/
static void bkw6x6f(
    complex_float* restrict x,
    const complex_float* restrict R,
    const complex_float* restrict D,
    const complex_float* restrict y)
{
    int m, k;
    xtfloat * restrict pX;
    const xtfloat * restrict pR;
    const xtfloat * restrict tpR = (xtfloat*)(R + 6 * (6 + 1) / 2 - 1) + 1; //last element in R
    const xtfloat * restrict pD = (xtfloat*)(D + (6 - 1)) + 1;
    const xtfloat * restrict pY = (xtfloat*)(y + (6 * 6 - 1)) + 1;
    xtfloat X0_re, X1_re, X2_re;
    xtfloat X0_im, X1_im, X2_im;
    xtfloat Xw0_re, Xw1_re, Xw2_re;
    xtfloat Xw0_im, Xw1_im, Xw2_im;
    xtfloat R0_re, R0_im, D0_re, D0_im;

    for (k = 6 - 1; k >= 0; k--)
    {
        pX = (xtfloat*)(x + 6 * 6 - 1) + 1; // last element in X
        pR = (xtfloat*)tpR; //points to the end of row k in R
        XT_LSXP(Xw2_im, pY, -SZ_F32);
        XT_LSXP(Xw2_re, pY, -SZ_F32);
        XT_LSXP(Xw1_im, pY, -SZ_F32);
        XT_LSXP(Xw1_re, pY, -SZ_F32);
        XT_LSXP(Xw0_im, pY, -SZ_F32);
        XT_LSXP(Xw0_re, pY, -SZ_F32);

        for (m = 0; m < 6 - k - 1; m++)
        {
            XT_LSXP(X2_im, pX, -SZ_F32);
            XT_LSXP(X2_re, pX, -SZ_F32);
            XT_LSXP(X1_im, pX, -SZ_F32);
            XT_LSXP(X1_re, pX, -SZ_F32);
            XT_LSXP(X0_im, pX, -SZ_F32);
            XT_LSXP(X0_re, pX, -SZ_F32 - 3 * SZ_CF32);
            XT_LSXP(R0_im, pR, -SZ_F32);
            XT_LSXP(R0_re, pR, -(6 - 1 - m)*SZ_CF32 + SZ_F32);
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
        }
        XT_LSXP(D0_im, pD, -SZ_F32);
        XT_LSXP(D0_re, pD, -SZ_F32);
        Xw0_re = XT_MUL_S(Xw0_re, D0_re);
        Xw0_im = XT_MUL_S(Xw0_im, D0_im);
        Xw1_re = XT_MUL_S(Xw1_re, D0_re);
        Xw1_im = XT_MUL_S(Xw1_im, D0_im);
        Xw2_re = XT_MUL_S(Xw2_re, D0_re);
        Xw2_im = XT_MUL_S(Xw2_im, D0_im);

        XT_SSXP(Xw2_im, pX, -SZ_F32);
        XT_SSXP(Xw2_re, pX, -SZ_F32);
        XT_SSXP(Xw1_im, pX, -SZ_F32);
        XT_SSXP(Xw1_re, pX, -SZ_F32);
        XT_SSXP(Xw0_im, pX, -SZ_F32);
        XT_SSXP(Xw0_re, pX, -SZ_F32);

        pR = (xtfloat*)tpR; //points to the end of row k in R
        XT_LSXP(Xw2_im, pY, -SZ_F32);
        XT_LSXP(Xw2_re, pY, -SZ_F32);
        XT_LSXP(Xw1_im, pY, -SZ_F32);
        XT_LSXP(Xw1_re, pY, -SZ_F32);
        XT_LSXP(Xw0_im, pY, -SZ_F32);
        XT_LSXP(Xw0_re, pY, -SZ_F32);

        pX = (xtfloat*)(x + 6 * 6 - 4) + 1; 
        for (m = 0; m < 6 - k - 1; m++)
        {
            XT_LSXP(X2_im, pX, -SZ_F32);
            XT_LSXP(X2_re, pX, -SZ_F32);
            XT_LSXP(X1_im, pX, -SZ_F32);
            XT_LSXP(X1_re, pX, -SZ_F32);
            XT_LSXP(X0_im, pX, -SZ_F32);
            XT_LSXP(X0_re, pX, -SZ_F32 - 3 * SZ_CF32);
            XT_LSXP(R0_im, pR, -SZ_F32);
            XT_LSXP(R0_re, pR, -(6 - 1 - m)*SZ_CF32 + SZ_F32);
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
        }
        Xw0_re = XT_MUL_S(Xw0_re, D0_re);
        Xw0_im = XT_MUL_S(Xw0_im, D0_im);
        Xw1_re = XT_MUL_S(Xw1_re, D0_re);
        Xw1_im = XT_MUL_S(Xw1_im, D0_im);
        Xw2_re = XT_MUL_S(Xw2_re, D0_re);
        Xw2_im = XT_MUL_S(Xw2_im, D0_im);

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
void  cmatcholpseudoinv6x6f(void* pScr,
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

    SD = 6;
    SR = (((6 + 1) * 6) >> 1);
    SY = 6 * 6;
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

    cmatcholdecomp6x6f(pScr, R, D, A, sigma2);
    fwd6x6f(y, R, D, A);
    bkw6x6f(x, R, D, y);
}

size_t  cmatcholpseudoinv6x6f_getScratchSize()
{
    size_t s_dc;
    int SD, SR, SY;

    SD = 6;
    SR = (((6 + 1) * 6) >> 1);
    SY = 6 * 6;
    SD = (SD + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SR = (SR + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SY = (SY + CF32_WIDTH - 1)&~(CF32_WIDTH - 1);
    SD = 2 * SD * sizeof(float32_t);
    SR = 2 * SR * sizeof(float32_t);
    SY = 2 * SY * sizeof(float32_t);

    s_dc = cmatcholdecomp6x6f_getScratchSize();
    return SD + SR + SY + s_dc;
}
#else
DISCARD_FUN(void, cmatcholpseudoinv6x6f, (void* pScr,
	complex_float *R,
	const complex_float * A,
	const float32_t sigma2))

	size_t  cmatcholpseudoinv6x6f_getScratchSize()
{
	return 0;
}
#endif
