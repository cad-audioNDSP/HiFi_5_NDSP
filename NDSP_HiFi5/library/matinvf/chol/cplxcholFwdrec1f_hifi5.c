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
void cplxcholFwdrec1f(complex_float * y,
	const complex_float * R,
	const complex_float * D,
	const complex_float * Z, int stride)
{
	xtfloatx2 * restrict pY = (xtfloatx2 *)y;
	const xtfloatx2 * restrict pZ = (const xtfloatx2 *)Z;
	const xtfloatx2 * restrict pD = (const xtfloatx2 *)D;
	xtfloatx2 Z0;
	xtfloatx2 D0;
	xtfloatx2 Y0;

	// Y0
	Z0 = XT_LSX2I(pZ, 0);
	D0 = XT_LSX2I(pD, 0);
	Y0 = XT_MUL_SX2(Z0, D0);
	XT_SSX2I(Y0, pY, 0);
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
void cplxcholFwdrec1f(complex_float * y,
    const complex_float * R,
    const complex_float * D,
    const complex_float * Z, int stride)
{
    xtfloat * restrict pY = (xtfloat*)y;
    const xtfloat * restrict pZ = (const xtfloat*)Z;
    const xtfloat * restrict pD = (const xtfloat*)D;
    xtfloat Y0_re;
    xtfloat Y0_im;
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

}
#endif
