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
void  cmatcholbkwsubst4x4f(void * pScr,
    complex_float * x,
    const complex_float * R,
    const complex_float * D,
    const complex_float * y)
{
    NASSERT(x);
    NASSERT(R);
    NASSERT(D);
    NASSERT(y);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    //NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    const xtfloatx2 * restrict pR3 = (const xtfloatx2*)(R + (4 * 3) / 2 + 0); //N*(N-1)/2
    const xtfloatx2 * restrict pR2 = (const xtfloatx2*)(R + (4 * 3) / 2 + 1);
    const xtfloatx2 * restrict pR1 = (const xtfloatx2*)(R + (4 * 3) / 2 + 2);
    const xtfloatx4 * restrict pY = (const xtfloatx4*)(y + 2);
    const xtfloatx4 * restrict pD = (const xtfloatx4*)(D + 2);
    xtfloatx4 * restrict pX = (xtfloatx4*)(x + 2);
    xtfloatx2 D0, D1;
    xtfloatx2 R0, R1, R2;
    xtfloatx2 X0, X1, X2, X3, Xt;

	// X0
	AE_LSX2X2_IP(D1, D0, pD, -2*SZ_CF32);
	AE_LSX2X2_IP(X1, X0, pY, -2*SZ_CF32);
	X0 = XT_MUL_SX2(X0, D0);

	// X1
	XT_LSX2IP(R0, pR1, 0);
	XT_MADDMUX_S(X1, X0, R0, 2);
	XT_MADDMUX_S(X1, X0, R0, 3);
	X1 = XT_MUL_SX2(X1, D1);
	AE_SSX2X2_IP(X1, X0, pX, -2*SZ_CF32);

	// X2
	Xt = XT_CONST_S(0);
	AE_LSX2X2_IP(D1, D0, pD, -2*SZ_CF32);
	AE_LSX2X2_IP(X3, X2, pY, -2*SZ_CF32);
	XT_LSX2XP(R0, pR2, -SZ_CF32 * 3); // N-1
	XT_LSX2IP(R1, pR2, 0);
	MADDMUX_SX2X2(X2, Xt, X0, X1, R0, R1, 2);
	MADDMUX_SX2X2(X2, Xt, X0, X1, R0, R1, 3);
	X2 = XT_ADD_SX2(X2, Xt);
	X2 = XT_MUL_SX2(X2, D0);

	// X3
	Xt = XT_CONST_S(0);
	XT_LSX2XP(R0, pR3, -SZ_CF32 * 3);
	XT_LSX2XP(R1, pR3, -SZ_CF32 * 2);
	XT_LSX2IP(R2, pR3, 0);
	MADDMUX_SX2X2(X3, Xt, X0, X1, R0, R1, 2);
	MADDMUX_SX2X2(X3, Xt, X0, X1, R0, R1, 3);
	XT_MADDMUX_S(X3, X2, R2, 2);
	XT_MADDMUX_S(Xt, X2, R2, 3);
	X3 = XT_ADD_SX2(X3, Xt);
	X3 = XT_MUL_SX2(X3, D1);
	AE_SSX2X2_IP(X3, X2, pX, -2*SZ_CF32);
}

size_t  cmatcholbkwsubst4x4f_getScratchSize()
{
    return 0;
}
#elif (HAVE_FPU)
#define SZ_F32 (int)(sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
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
void  cmatcholbkwsubst4x4f(void * pScr,
    complex_float * x,
    const complex_float * R,
    const complex_float * D,
    const complex_float * y)
{
    NASSERT(x);
    NASSERT(R);
    NASSERT(D);
    NASSERT(y);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    //NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    const xtfloat * restrict pR3 = (const xtfloat*)(R + (4 * 3) / 2 + 0) + 1; //N*(N-1)/2
    const xtfloat * restrict pR2 = (const xtfloat*)(R + (4 * 3) / 2 + 1) + 1;
    const xtfloat * restrict pR1 = (const xtfloat*)(R + (4 * 3) / 2 + 2) + 1;
    const xtfloat * restrict pY = (const xtfloat*)(y + 3) + 1;
    const xtfloat * restrict pD = (const xtfloat*)(D + 3) + 1;
    xtfloat * restrict pX = (xtfloat*)(x + 3) + 1;
    xtfloat D0_re, D0_im;
    xtfloat R1_re, R2_re, R3_re;
    xtfloat R1_im, R2_im, R3_im;
    xtfloat X0_re, X1_re, X2_re, X3_re;
    xtfloat X0_im, X1_im, X2_im, X3_im;

    // X0
    XT_LSXP(D0_im, pD, -SZ_F32);
    XT_LSXP(D0_re, pD, -SZ_F32);
    XT_LSXP(X0_im, pY, -SZ_F32);
    XT_LSXP(X0_re, pY, -SZ_F32);
    X0_re = XT_MUL_S(X0_re, D0_re);
    X0_im = XT_MUL_S(X0_im, D0_im);
    XT_SSXP(X0_im, pX, -SZ_F32);
    XT_SSXP(X0_re, pX, -SZ_F32);

    // X1
    XT_LSXP(D0_im, pD, -SZ_F32);
    XT_LSXP(D0_re, pD, -SZ_F32);
    XT_LSXP(X1_im, pY, -SZ_F32);
    XT_LSXP(X1_re, pY, -SZ_F32);
    XT_LSXP(R1_im, pR1, -SZ_F32);
    XT_LSIP(R1_re, pR1, 0);
    XT_MSUB_S(X1_re, X0_re, R1_re);
    XT_MADD_S(X1_re, X0_im, R1_im);
    XT_MSUB_S(X1_im, X0_re, R1_im);
    XT_MSUB_S(X1_im, X0_im, R1_re);
    X1_re = XT_MUL_S(X1_re, D0_re);
    X1_im = XT_MUL_S(X1_im, D0_im);
    XT_SSXP(X1_im, pX, -SZ_F32);
    XT_SSXP(X1_re, pX, -SZ_F32);

    // X2
    XT_LSXP(D0_im, pD, -SZ_F32);
    XT_LSXP(D0_re, pD, -SZ_F32);
    XT_LSXP(X2_im, pY, -SZ_F32);
    XT_LSXP(X2_re, pY, -SZ_F32);
    XT_LSXP(R2_im, pR2, -SZ_F32);
    XT_LSXP(R2_re, pR2, -3 * SZ_CF32 + SZ_F32); //N-1
    XT_MSUB_S(X2_re, X0_re, R2_re);
    XT_MADD_S(X2_re, X0_im, R2_im);
    XT_MSUB_S(X2_im, X0_re, R2_im);
    XT_MSUB_S(X2_im, X0_im, R2_re);
    XT_LSXP(R2_im, pR2, -SZ_F32);
    XT_LSIP(R2_re, pR2, 0);
    XT_MSUB_S(X2_re, X1_re, R2_re);
    XT_MADD_S(X2_re, X1_im, R2_im);
    XT_MSUB_S(X2_im, X1_re, R2_im);
    XT_MSUB_S(X2_im, X1_im, R2_re);
    X2_re = XT_MUL_S(X2_re, D0_re);
    X2_im = XT_MUL_S(X2_im, D0_im);
    XT_SSXP(X2_im, pX, -SZ_F32);
    XT_SSXP(X2_re, pX, -SZ_F32);

    // X3
    XT_LSXP(D0_im, pD, -SZ_F32);
    XT_LSXP(D0_re, pD, -SZ_F32);
    XT_LSXP(X3_im, pY, -SZ_F32);
    XT_LSXP(X3_re, pY, -SZ_F32);
    XT_LSXP(R3_im, pR3, -SZ_F32);
    XT_LSXP(R3_re, pR3, -3 * SZ_CF32 + SZ_F32); //N-1
    XT_MSUB_S(X3_re, X0_re, R3_re);
    XT_MADD_S(X3_re, X0_im, R3_im);
    XT_MSUB_S(X3_im, X0_re, R3_im);
    XT_MSUB_S(X3_im, X0_im, R3_re);
    XT_LSXP(R3_im, pR3, -SZ_F32);
    XT_LSXP(R3_re, pR3, -2 * SZ_CF32 + SZ_F32); 
    XT_MSUB_S(X3_re, X1_re, R3_re);
    XT_MADD_S(X3_re, X1_im, R3_im);
    XT_MSUB_S(X3_im, X1_re, R3_im);
    XT_MSUB_S(X3_im, X1_im, R3_re);
    XT_LSXP(R3_im, pR3, -SZ_F32);
    XT_LSIP(R3_re, pR3, 0);
    XT_MSUB_S(X3_re, X2_re, R3_re);
    XT_MADD_S(X3_re, X2_im, R3_im);
    XT_MSUB_S(X3_im, X2_re, R3_im);
    XT_MSUB_S(X3_im, X2_im, R3_re);
    X3_re = XT_MUL_S(X3_re, D0_re);
    X3_im = XT_MUL_S(X3_im, D0_im);
    XT_SSXP(X3_im, pX, -SZ_F32);
    XT_SSXP(X3_re, pX, -SZ_F32);
}

size_t  cmatcholbkwsubst4x4f_getScratchSize()
{
    return 0;
}
#else
DISCARD_FUN(void, cmatcholbkwsubst4x4f, (void * pScr,
	complex_float * x,
	const complex_float * R,
	const complex_float * D,
	const complex_float * y))

	size_t  cmatcholbkwsubst4x4f_getScratchSize()
{
	return 0;
}
#endif
