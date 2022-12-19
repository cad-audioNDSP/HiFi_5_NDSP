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
#include "common_fpu.h"
#include "cholnf_common.h"
/*
code optimized for HiFi4 with VFPU
*/

#if (HAVE_VFPU)
#define SZ_F32 (sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
/*--------------------------------------------------------
Forward recursion P=1
Input:
R[((N+1)*N)/2]  upper triangular matrix R
stride          width of matrix A'*B
Z[N*1]          column in matrix A'*B
D[N]            reciprocals of main diagonal
Output:
y[N*1]		Decision matrix y
N = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
--------------------------------------------------------*/
void cplxcholFwdrec2f(complex_float * y,
	const complex_float * R,
	const complex_float * D,
	const complex_float * Z, int stride)
{
	xtfloatx2 * restrict pY = (xtfloatx2 *)y;
    const xtfloatx2 * restrict pZ = (const xtfloatx2 *)Z;
    const xtfloatx2 * restrict pR1 = (const xtfloatx2 *)(R + 1);/* + ((n*(n+1))/2) */
    const xtfloatx2 * restrict pD = (const xtfloatx2 *)D;
	xtfloatx2 Z0, Z1;
	xtfloatx2 Z1t;
	xtfloatx2 D0, D1;
	xtfloatx2 Y0, Y1;
	xtfloatx2 R0;

	Z1t = XT_CONST_S(0);
	// Y0
	XT_LSX2XP(Z0, pZ, stride*SZ_CF32);
	AE_LSX2X2_I(D0, D1, (xtfloatx4*)pD, 0 * SZ_CF32);
	Y0 = XT_MUL_SX2(Z0, D0);
	XT_SSX2I(Y0, pY, 0 * SZ_CF32);

	// Y1
	XT_LSX2XP(Z1, pZ, stride*SZ_CF32);
	R0 = XT_LSX2I(pR1, 0);
	XT_MADDMUX_S(Z1, R0, Y0, 2);
	XT_MADDMUX_S(Z1t, R0, Y0, 1);
	Z1 = XT_ADD_SX2(Z1, Z1t);
	Y1 = XT_MUL_SX2(Z1, D1);
	XT_SSX2I(Y1, pY, 1 * SZ_CF32);
}
#elif (HAVE_FPU)
#define SZ_F32 (sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
/*--------------------------------------------------------
Forward recursion P=1
Input:
R[((N+1)*N)/2]  upper triangular matrix R
stride          width of matrix A'*B
Z[N*1]          column in matrix A'*B
D[N]            reciprocals of main diagonal
Output:
y[N*1]		Decision matrix y
N = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
--------------------------------------------------------*/
void cplxcholFwdrec2f(complex_float * y,
    const complex_float * R,
    const complex_float * D,
    const complex_float * Z, int stride)
{
    xtfloat * restrict pY = (xtfloat*)y;
    const xtfloat * restrict pZ = (const xtfloat*)Z;
    const xtfloat * restrict pD = (const xtfloat*)D;
    const xtfloat * restrict pR = (const xtfloat*)(R + 1);
    xtfloat Y0_re, Y1_re;
    xtfloat Y0_im, Y1_im;
    xtfloat R_re, R_im;
    xtfloat D_re, D_im;

    //Y0
    XT_LSIP(Y0_re, pZ, SZ_F32);
    XT_LSXP(Y0_im, pZ, stride * SZ_CF32 - SZ_F32);
    XT_LSIP(D_re, pD, SZ_F32);
    XT_LSIP(D_im, pD, SZ_F32);

    Y0_re = XT_MUL_S(Y0_re, D_re);
    Y0_im = XT_MUL_S(Y0_im, D_im);
    XT_SSIP(Y0_re, pY, SZ_F32);
    XT_SSIP(Y0_im, pY, SZ_F32);

    //Y1
    XT_LSIP(Y1_re, pZ, SZ_F32);
    XT_LSXP(Y1_im, pZ, stride * SZ_CF32 - SZ_F32);
    XT_LSIP(D_re, pD, SZ_F32);
    XT_LSIP(D_im, pD, SZ_F32);

    XT_LSIP(R_re, pR, SZ_F32);
    XT_LSIP(R_im, pR, 3 * SZ_F32);
    XT_MSUB_S(Y1_re, Y0_re, R_re);
    XT_MSUB_S(Y1_re, Y0_im, R_im);
    XT_MSUB_S(Y1_im, Y0_im, R_re);
    XT_MADD_S(Y1_im, Y0_re, R_im);

    Y1_re = XT_MUL_S(Y1_re, D_re);
    Y1_im = XT_MUL_S(Y1_im, D_im);
    XT_SSIP(Y1_re, pY, SZ_F32);
    XT_SSIP(Y1_im, pY, SZ_F32);

}
#endif
