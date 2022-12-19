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
Preprocessing for Least Square Solutions
The result is matrix Z[NxN], such that
Z = A'*A + sigma2*I[NxN], where ' denotes the conjugate transpose of
a matrix, and sigma2*I[NxN] is the NxN identity matrix multiplied with
the regularization term.

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
A[M*N]          matrix A. For the fixed point, the representation is Q(qA)
sigma2          regularization term. For fixed point, the representation 
                should be Q(2*qA-30)
qRA             qR-qA; difference between fixed point representations of 
                decomposition matrix R and original matrix A (for the fixed 
                point API only). Should be equal or less than 0 (typically 
                -2).
Output:
Z[N*N]          matrix Z. For the fixed point, the representation is Q(2*qR-4)
Temporary:
pScr            Scratch memory

N = M = 4, 6, 8, 10

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
16-bytes boundary
---------------------------------------------------------------------------*/
/* single-matrix API */
void  cmatcholpreprocess10x10f(void* pScr,
	complex_float *R,
	const complex_float * A,
	const float32_t  sigma2)
{
    xtfloatx2 sigma = sigma2;
    xtfloatx2 sigma_w_zeros;

    NASSERT(R);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    sigma_w_zeros = XT_AE_MOVXTFLOATX2_FROMF32X2(AE_MOVF32X2_FROMINT64(AE_SLAI64(AE_MOVINT64_FROMF32X2(XT_AE_MOVF32X2_FROMXTFLOATX2(sigma)), 32)));

	const xtfloatx4 * restrict pApn = (xtfloatx4 *)A;
	int p;
	xtfloatx4 * restrict pRv = (xtfloatx4 *)R;
	xtfloatx4 * restrict pRh;
	xtfloatx2 Acc00, Acc01, Acc02, Acc03, Acc04, Acc05, Acc06, Acc07, Acc08, Acc09;
	xtfloatx2        Acc11, Acc12, Acc13, Acc14, Acc15, Acc16, Acc17, Acc18, Acc19;
	xtfloatx2               Acc22, Acc23, Acc24, Acc25, Acc26, Acc27, Acc28, Acc29;
	xtfloatx2                      Acc33, Acc34, Acc35, Acc36, Acc37, Acc38, Acc39;
	xtfloatx2                             Acc44, Acc45, Acc46, Acc47, Acc48, Acc49;
	xtfloatx2                                    Acc55, Acc56, Acc57, Acc58, Acc59;
	xtfloatx2                                           Acc66, Acc67, Acc68, Acc69;
	xtfloatx2                                                  Acc77, Acc78, Acc79;
	xtfloatx2                                                         Acc88, Acc89;
	xtfloatx2                                                                Acc99;
	xtfloatx2 Apn0, Apn1, Apn2, Apn3, Apn4, Apn5, Apn6, Apn7, Apn8, Apn9;
	Acc00 = Acc11 = sigma_w_zeros;
	Acc01 = XT_CONST_S(0);
	CONST_SX2X2(Acc02, Acc12, 0);
	CONST_SX2X2(Acc03, Acc13, 0);
	CONST_SX2X2(Acc04, Acc14, 0);
	CONST_SX2X2(Acc05, Acc15, 0);
	CONST_SX2X2(Acc06, Acc16, 0);
	CONST_SX2X2(Acc07, Acc17, 0);
	CONST_SX2X2(Acc08, Acc18, 0);
	CONST_SX2X2(Acc09, Acc19, 0);
	                       
	/* FIRST PART */
	for (p = 0; p < 10; p++)
	{
		AE_LSX2X2_IP(Apn0, Apn1, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn2, Apn3, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn4, Apn5, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn6, Apn7, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn8, Apn9, pApn, 2 * SZ_CF32);

		MADDMUX_SX2X2(Acc00, Acc01, Apn0, Apn0, Apn0, Apn1, 0);
		MADDMUX_SX2X2(Acc00, Acc01, Apn0, Apn0, Apn0, Apn1, 3);
		MADDMUX_SX2X2(Acc02, Acc03, Apn0, Apn0, Apn2, Apn3, 0);
		MADDMUX_SX2X2(Acc02, Acc03, Apn0, Apn0, Apn2, Apn3, 3);
		MADDMUX_SX2X2(Acc04, Acc05, Apn0, Apn0, Apn4, Apn5, 0);
		MADDMUX_SX2X2(Acc04, Acc05, Apn0, Apn0, Apn4, Apn5, 3);
		MADDMUX_SX2X2(Acc06, Acc07, Apn0, Apn0, Apn6, Apn7, 0);
		MADDMUX_SX2X2(Acc06, Acc07, Apn0, Apn0, Apn6, Apn7, 3);
		MADDMUX_SX2X2(Acc08, Acc09, Apn0, Apn0, Apn8, Apn9, 0);
		MADDMUX_SX2X2(Acc08, Acc09, Apn0, Apn0, Apn8, Apn9, 3);

		MADDMUX_SX2X2(Acc11, Acc12, Apn1, Apn1, Apn1, Apn2, 0);
		MADDMUX_SX2X2(Acc11, Acc12, Apn1, Apn1, Apn1, Apn2, 3);
		MADDMUX_SX2X2(Acc13, Acc14, Apn1, Apn1, Apn3, Apn4, 0);
		MADDMUX_SX2X2(Acc13, Acc14, Apn1, Apn1, Apn3, Apn4, 3);
		MADDMUX_SX2X2(Acc15, Acc16, Apn1, Apn1, Apn5, Apn6, 0);
		MADDMUX_SX2X2(Acc15, Acc16, Apn1, Apn1, Apn5, Apn6, 3);
		MADDMUX_SX2X2(Acc17, Acc18, Apn1, Apn1, Apn7, Apn8, 0);
		MADDMUX_SX2X2(Acc17, Acc18, Apn1, Apn1, Apn7, Apn8, 3);
		XT_MADDMUX_S(Acc19, Apn1, Apn9, 0);
		XT_MADDMUX_S(Acc19, Apn1, Apn9, 3);
	}
	AE_SSX2X2_IP(Acc00, Acc01, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc02, Acc03, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc04, Acc05, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc06, Acc07, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc08, Acc09, pRv, 2 * SZ_CF32);

	AE_SSX2X2_IP(XT_CONJC_S(Acc01), Acc11, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc12, Acc13, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc14, Acc15, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc16, Acc17, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc18, Acc19, pRv, 4 * SZ_CF32);

	pRh = (xtfloatx4 *)R + 10;

	CONJC_SX2X2(Acc02, Acc12, Acc02, Acc12);
	CONJC_SX2X2(Acc03, Acc13, Acc03, Acc13);
	CONJC_SX2X2(Acc04, Acc14, Acc04, Acc14);
	CONJC_SX2X2(Acc05, Acc15, Acc05, Acc15);
	CONJC_SX2X2(Acc06, Acc16, Acc06, Acc16);
	CONJC_SX2X2(Acc07, Acc17, Acc07, Acc17);
	CONJC_SX2X2(Acc08, Acc18, Acc08, Acc18);
	CONJC_SX2X2(Acc09, Acc19, Acc09, Acc19);

	AE_SSX2X2_IP(Acc02, Acc12, pRh, 10 * SZ_CF32);
	AE_SSX2X2_IP(Acc03, Acc13, pRh, 10 * SZ_CF32);
	AE_SSX2X2_IP(Acc04, Acc14, pRh, 10 * SZ_CF32);
	AE_SSX2X2_IP(Acc05, Acc15, pRh, 10 * SZ_CF32);
	AE_SSX2X2_IP(Acc06, Acc16, pRh, 10 * SZ_CF32);
	AE_SSX2X2_IP(Acc07, Acc17, pRh, 10 * SZ_CF32);
	AE_SSX2X2_IP(Acc08, Acc18, pRh, 10 * SZ_CF32);
	AE_SSX2X2_IP(Acc09, Acc19, pRh, 10 * SZ_CF32);

	__Pragma("flush_memory");

	/* SECOND PART */
	Acc22 = Acc33 = sigma_w_zeros;
	Acc23 = XT_CONST_S(0);
	CONST_SX2X2(Acc24, Acc34, 0);
	CONST_SX2X2(Acc25, Acc35, 0);
	CONST_SX2X2(Acc26, Acc36, 0);
	CONST_SX2X2(Acc27, Acc37, 0);
	CONST_SX2X2(Acc28, Acc38, 0);
	CONST_SX2X2(Acc29, Acc39, 0);

	
	pApn = (xtfloatx4 *)A;
	for (p = 0; p < 10; p++)
	{
		pApn++;
		AE_LSX2X2_IP(Apn2, Apn3, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn4, Apn5, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn6, Apn7, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn8, Apn9, pApn, 2 * SZ_CF32);

		XT_MADDMUX_S(Acc22, Apn2, Apn2, 0);
		XT_MADDMUX_S(Acc22, Apn2, Apn2, 3);
		
		MADDMUX_SX2X2(Acc23, Acc24, Apn2, Apn2, Apn3, Apn4, 0);
		MADDMUX_SX2X2(Acc23, Acc24, Apn2, Apn2, Apn3, Apn4, 3);
		MADDMUX_SX2X2(Acc25, Acc26, Apn2, Apn2, Apn5, Apn6, 0);
		MADDMUX_SX2X2(Acc25, Acc26, Apn2, Apn2, Apn5, Apn6, 3);
		MADDMUX_SX2X2(Acc27, Acc28, Apn2, Apn2, Apn7, Apn8, 0);
		MADDMUX_SX2X2(Acc27, Acc28, Apn2, Apn2, Apn7, Apn8, 3);

		MADDMUX_SX2X2(Acc29, Acc33, Apn2, Apn3, Apn9, Apn3, 0);
		MADDMUX_SX2X2(Acc29, Acc33, Apn2, Apn3, Apn9, Apn3, 3);
		MADDMUX_SX2X2(Acc34, Acc35, Apn3, Apn3, Apn4, Apn5, 0);
		MADDMUX_SX2X2(Acc34, Acc35, Apn3, Apn3, Apn4, Apn5, 3);
		MADDMUX_SX2X2(Acc36, Acc37, Apn3, Apn3, Apn6, Apn7, 0);
		MADDMUX_SX2X2(Acc36, Acc37, Apn3, Apn3, Apn6, Apn7, 3);
		MADDMUX_SX2X2(Acc38, Acc39, Apn3, Apn3, Apn8, Apn9, 0);
		MADDMUX_SX2X2(Acc38, Acc39, Apn3, Apn3, Apn8, Apn9, 3);
	}
	AE_SSX2X2_IP(Acc22, Acc23, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc24, Acc25, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc26, Acc27, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc28, Acc29, pRv, 4 * SZ_CF32);

	AE_SSX2X2_IP(XT_CONJC_S(Acc23), Acc33, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc34, Acc35, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc36, Acc37, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc38, Acc39, pRv, 6 * SZ_CF32);

	pRh = (xtfloatx4 *)R + 21;

	CONJC_SX2X2(Acc24, Acc34, Acc24, Acc34);
	CONJC_SX2X2(Acc25, Acc35, Acc25, Acc35);
	CONJC_SX2X2(Acc26, Acc36, Acc26, Acc36);
	CONJC_SX2X2(Acc27, Acc37, Acc27, Acc37);
	CONJC_SX2X2(Acc28, Acc38, Acc28, Acc38);
	CONJC_SX2X2(Acc29, Acc39, Acc29, Acc39);

	AE_SSX2X2_IP(Acc24, Acc34, pRh, 10 * SZ_CF32);
	AE_SSX2X2_IP(Acc25, Acc35, pRh, 10 * SZ_CF32);
	AE_SSX2X2_IP(Acc26, Acc36, pRh, 10 * SZ_CF32);
	AE_SSX2X2_IP(Acc27, Acc37, pRh, 10 * SZ_CF32);
	AE_SSX2X2_IP(Acc28, Acc38, pRh, 10 * SZ_CF32);
	AE_SSX2X2_IP(Acc29, Acc39, pRh, 10 * SZ_CF32);

	__Pragma("flush_memory");

	/* THIRD PART */
	Acc44 = Acc55 = Acc66 = Acc77 = Acc88 = Acc99 = sigma_w_zeros;
	Acc45 = XT_CONST_S(0);
	CONST_SX2X2(Acc46, Acc56, 0);
	CONST_SX2X2(Acc47, Acc57, 0);
	CONST_SX2X2(Acc48, Acc58, 0);
	CONST_SX2X2(Acc49, Acc59, 0);
	CONST_SX2X2(Acc67, Acc68, 0);
	CONST_SX2X2(Acc69, Acc78, 0);
	CONST_SX2X2(Acc79, Acc89, 0);

	pApn = (xtfloatx4 *)A;
	for (p = 0; p < 10; p++)
	{
		pApn += 2;
		AE_LSX2X2_IP(Apn4, Apn5, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn6, Apn7, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn8, Apn9, pApn, 2 * SZ_CF32);

		MADDMUX_SX2X2(Acc44, Acc45, Apn4, Apn4, Apn4, Apn5, 0);
		MADDMUX_SX2X2(Acc44, Acc45, Apn4, Apn4, Apn4, Apn5, 3);
		MADDMUX_SX2X2(Acc46, Acc47, Apn4, Apn4, Apn6, Apn7, 0);
		MADDMUX_SX2X2(Acc46, Acc47, Apn4, Apn4, Apn6, Apn7, 3);
		MADDMUX_SX2X2(Acc48, Acc49, Apn4, Apn4, Apn8, Apn9, 0);
		MADDMUX_SX2X2(Acc48, Acc49, Apn4, Apn4, Apn8, Apn9, 3);

		MADDMUX_SX2X2(Acc55, Acc56, Apn5, Apn5, Apn5, Apn6, 0);
		MADDMUX_SX2X2(Acc55, Acc56, Apn5, Apn5, Apn5, Apn6, 3);
		MADDMUX_SX2X2(Acc57, Acc58, Apn5, Apn5, Apn7, Apn8, 0);
		MADDMUX_SX2X2(Acc57, Acc58, Apn5, Apn5, Apn7, Apn8, 3);

		MADDMUX_SX2X2(Acc59, Acc66, Apn5, Apn6, Apn9, Apn6, 0);
		MADDMUX_SX2X2(Acc59, Acc66, Apn5, Apn6, Apn9, Apn6, 3);
		MADDMUX_SX2X2(Acc67, Acc68, Apn6, Apn6, Apn7, Apn8, 0);
		MADDMUX_SX2X2(Acc67, Acc68, Apn6, Apn6, Apn7, Apn8, 3);

		MADDMUX_SX2X2(Acc69, Acc77, Apn6, Apn7, Apn9, Apn7, 0);
		MADDMUX_SX2X2(Acc69, Acc77, Apn6, Apn7, Apn9, Apn7, 3);
		MADDMUX_SX2X2(Acc78, Acc79, Apn7, Apn7, Apn8, Apn9, 0);
		MADDMUX_SX2X2(Acc78, Acc79, Apn7, Apn7, Apn8, Apn9, 3);

		MADDMUX_SX2X2(Acc88, Acc89, Apn8, Apn8, Apn8, Apn9, 0);
		MADDMUX_SX2X2(Acc88, Acc89, Apn8, Apn8, Apn8, Apn9, 3);

		XT_MADDMUX_S(Acc99, Apn9, Apn9, 0);
		XT_MADDMUX_S(Acc99, Apn9, Apn9, 3);
	}
	
	AE_SSX2X2_IP(Acc44, Acc45, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc46, Acc47, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc48, Acc49, pRv, 6 * SZ_CF32);

	AE_SSX2X2_IP(XT_CONJC_S(Acc45), Acc55, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc56, Acc57, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc58, Acc59, pRv, 6 * SZ_CF32);

	AE_SSX2X2_IP(XT_CONJC_S(Acc46), XT_CONJC_S(Acc56), pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc66, Acc67, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc68, Acc69, pRv, 6 * SZ_CF32);

	AE_SSX2X2_IP(XT_CONJC_S(Acc47), XT_CONJC_S(Acc57), pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(XT_CONJC_S(Acc67), Acc77, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc78, Acc79, pRv, 6 * SZ_CF32);

	AE_SSX2X2_IP(XT_CONJC_S(Acc48), XT_CONJC_S(Acc58), pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(XT_CONJC_S(Acc68), XT_CONJC_S(Acc78), pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc88, Acc89, pRv, 6 * SZ_CF32);

	AE_SSX2X2_IP(XT_CONJC_S(Acc49), XT_CONJC_S(Acc59), pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(XT_CONJC_S(Acc69), XT_CONJC_S(Acc79), pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(XT_CONJC_S(Acc89), Acc99, pRv, 0);
}

size_t  cmatcholpreprocess10x10f_getScratchSize()
{
	return 0;
}
#elif(HAVE_FPU)
#define SZ_F32 (sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)
/*-------------------------------------------------------------------------
Preprocessing for Least Square Solutions
The result is matrix Z[NxN], such that
Z = A'*A + sigma2*I[NxN], where ' denotes the conjugate transpose of
a matrix, and sigma2*I[NxN] is the NxN identity matrix multiplied with
the regularization term.

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
A[M*N]          matrix A. For the fixed poiont, the representation is Q(qA)
sigma2          regularization term. For fixed point, the representation
should be Q(2*qA-30)
qRA             qR-qA; difference between fixed point representations of
decomposition matrix R and original matrix A (for the fixed
point API only). Should be equal or less than 0 (typically
-2).
Output:
Z[N*N]          matrix Z. For the fixed poiont, the representation is Q(2*qR-4)
Temporary:
pScr            Scratch memory

N = M = 4, 6, 8, 10

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
16-bytes boundary
---------------------------------------------------------------------------*/
void  cmatcholpreprocess10x10f(void* pScr,
    complex_float *R,
    const complex_float * A,
    const float32_t  sigma2)
{
#if 0
    const float32_t* restrict Af = (const float32_t*)A;
    float32_t* restrict Rf = (float32_t*)R;
    int m, n, p;

    NASSERT(R);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    float32_t re_Apm;
    float32_t im_Apm;
    float32_t re_Apn;
    float32_t im_Apn;
    float32_t re_acc;
    float32_t im_acc;
    for (m = 0; m < 10; m++)
    {
        re_acc = sigma2;
        for (p = 0; p < 10; p++)
        {
            re_Apm = Af[(p * 10 + m) * 2 + 0];
            im_Apm = Af[(p * 10 + m) * 2 + 1];

            re_acc += re_Apm*re_Apm + im_Apm*im_Apm;
        }
        Rf[(m * 10 + m) * 2 + 0] = re_acc;
        Rf[(m * 10 + m) * 2 + 1] = 0.f;

        for (n = m + 1; n < 10; n++)
        {
            re_acc = 0.f;//m == n ? sigma2 : 0.f;
            im_acc = 0.f;

            for (p = 0; p < 10; p++)
            {
                re_Apm = Af[(p * 10 + m) * 2 + 0];
                im_Apm = Af[(p * 10 + m) * 2 + 1];
                re_Apn = Af[(p * 10 + n) * 2 + 0];
                im_Apn = Af[(p * 10 + n) * 2 + 1];

                re_acc += re_Apm*re_Apn + im_Apm*im_Apn;
                im_acc += re_Apm*im_Apn - im_Apm*re_Apn;
            }
            Rf[(m * 10 + n) * 2 + 0] = re_acc;
            Rf[(m * 10 + n) * 2 + 1] = im_acc;
            Rf[(m + n * 10) * 2 + 0] = re_acc;
            Rf[(m + n * 10) * 2 + 1] = -im_acc;
        }
    }
#elif 1
    const float32_t* restrict Af = (const float32_t*)A;
    float32_t* restrict Rf = (float32_t*)R;
    int m, n, p;

    NASSERT(R);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    float32_t re_Apm0, re_Apm1;
    float32_t im_Apm0, im_Apm1;
    float32_t re_Apn0, re_Apn1;
    float32_t im_Apn0, im_Apn1;
    float32_t re_acc00, re_acc01, re_acc10, re_acc11;
    float32_t im_acc00, im_acc01, im_acc10, im_acc11;
    /* processing by squares 2x2 */
    for (m = 0; m < 10; m += 2)
    {
        /* diagonal squares */
        re_acc00 = re_acc11 = sigma2;
        re_acc01 = re_acc10 = 0.f;
        im_acc00 = im_acc01 = im_acc10 = im_acc11 = 0.f;

        for (p = 0; p < 10; p++)
        {
            re_Apm0 = Af[(p * 10 + (m + 0)) * 2 + 0];
            re_Apm1 = Af[(p * 10 + (m + 1)) * 2 + 0];
            im_Apm0 = Af[(p * 10 + (m + 0)) * 2 + 1];
            im_Apm1 = Af[(p * 10 + (m + 1)) * 2 + 1];

            re_acc00 += re_Apm0*re_Apm0 + im_Apm0*im_Apm0;
            re_acc01 += re_Apm0*re_Apm1 + im_Apm0*im_Apm1;
            re_acc11 += re_Apm1*re_Apm1 + im_Apm1*im_Apm1;
            im_acc00 += re_Apm0*im_Apm0 - im_Apm0*re_Apm0;
            im_acc01 += re_Apm0*im_Apm1 - im_Apm0*re_Apm1;
            im_acc11 += re_Apm1*im_Apm1 - im_Apm1*re_Apm1;
        }
        Rf[((m + 0) * 10 + (m + 0)) * 2 + 0] = re_acc00;
        Rf[((m + 0) * 10 + (m + 0)) * 2 + 1] = im_acc00;
        Rf[((m + 0) * 10 + (m + 1)) * 2 + 0] = re_acc01;
        Rf[((m + 0) * 10 + (m + 1)) * 2 + 1] = im_acc01;
        Rf[((m + 1) * 10 + (m + 0)) * 2 + 0] = re_acc01;
        Rf[((m + 1) * 10 + (m + 0)) * 2 + 1] = -im_acc01;
        Rf[((m + 1) * 10 + (m + 1)) * 2 + 0] = re_acc11;
        Rf[((m + 1) * 10 + (m + 1)) * 2 + 1] = im_acc11;

        /* rest */
        for (n = m + 2; n < 10; n += 2)
        {
            re_acc00 = re_acc11 = re_acc01 = re_acc10 = 0.f;
            im_acc00 = im_acc01 = im_acc10 = im_acc11 = 0.f;

            for (p = 0; p < 10; p++)
            {
                re_Apm0 = Af[(p * 10 + (m + 0)) * 2 + 0];
                re_Apm1 = Af[(p * 10 + (m + 1)) * 2 + 0];
                re_Apn0 = Af[(p * 10 + (n + 0)) * 2 + 0];
                re_Apn1 = Af[(p * 10 + (n + 1)) * 2 + 0];
                im_Apm0 = Af[(p * 10 + (m + 0)) * 2 + 1];
                im_Apm1 = Af[(p * 10 + (m + 1)) * 2 + 1];
                im_Apn0 = Af[(p * 10 + (n + 0)) * 2 + 1];
                im_Apn1 = Af[(p * 10 + (n + 1)) * 2 + 1];

                re_acc00 += re_Apm0*re_Apn0 + im_Apm0*im_Apn0;
                re_acc01 += re_Apm0*re_Apn1 + im_Apm0*im_Apn1;
                re_acc10 += re_Apm1*re_Apn0 + im_Apm1*im_Apn0;
                re_acc11 += re_Apm1*re_Apn1 + im_Apm1*im_Apn1;
                im_acc00 += re_Apm0*im_Apn0 - im_Apm0*re_Apn0;
                im_acc01 += re_Apm0*im_Apn1 - im_Apm0*re_Apn1;
                im_acc10 += re_Apm1*im_Apn0 - im_Apm1*re_Apn0;
                im_acc11 += re_Apm1*im_Apn1 - im_Apm1*re_Apn1;
            }
            Rf[((m + 0) * 10 + (n + 0)) * 2 + 0] = re_acc00;
            Rf[((m + 0) * 10 + (n + 0)) * 2 + 1] = im_acc00;
            Rf[((m + 0) * 10 + (n + 1)) * 2 + 0] = re_acc01;
            Rf[((m + 0) * 10 + (n + 1)) * 2 + 1] = im_acc01;
            Rf[((m + 1) * 10 + (n + 0)) * 2 + 0] = re_acc10;
            Rf[((m + 1) * 10 + (n + 0)) * 2 + 1] = im_acc10;
            Rf[((m + 1) * 10 + (n + 1)) * 2 + 0] = re_acc11;
            Rf[((m + 1) * 10 + (n + 1)) * 2 + 1] = im_acc11;

            Rf[((m + 0) + (n + 0) * 10) * 2 + 0] = re_acc00;
            Rf[((m + 0) + (n + 0) * 10) * 2 + 1] = -im_acc00;
            Rf[((m + 0) + (n + 1) * 10) * 2 + 0] = re_acc01;
            Rf[((m + 0) + (n + 1) * 10) * 2 + 1] = -im_acc01;
            Rf[((m + 1) + (n + 0) * 10) * 2 + 0] = re_acc10;
            Rf[((m + 1) + (n + 0) * 10) * 2 + 1] = -im_acc10;
            Rf[((m + 1) + (n + 1) * 10) * 2 + 0] = re_acc11;
            Rf[((m + 1) + (n + 1) * 10) * 2 + 1] = -im_acc11;

        }
    }
#else
    const float32_t* restrict Af = (const float32_t*)A;
    float32_t* restrict Rf = (float32_t*)R;
    int m, n, p;

    NASSERT(R);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    xtfloat * restrict pRv;
    xtfloat * restrict pRh;
    const xtfloat * restrict pApm;
    const xtfloat * restrict pApn;
    xtfloat R_re, R_im;
    xtfloat Apm_re, Apm_im;
    xtfloat Apn_re, Apn_im;

    for (m = 0; m < 10; m++)
    {
        R_re = sigma2;
        pApm = (xtfloat*)(Af + (0 * 10 + m) * 2);
        pRv = (xtfloat*)(Rf + (m * 10 + m) * 2);
        pRh = (xtfloat*)(Rf + (m + (m + 1) * 10) * 2);

        for (p = 0; p < 10; p++)
        {
            XT_LSIP(Apm_re, pApm, SZ_F32);
            XT_LSIP(Apm_im, pApm, (10 * 2 - 1) * SZ_F32);
            XT_MADD_S(R_re, Apm_re, Apm_re);
            XT_MADD_S(R_re, Apm_im, Apm_im);
        }
        XT_SSIP(R_re, pRv, SZ_F32);
        XT_SSIP(XT_CONST_S(0), pRv, SZ_F32);

        for (n = m + 1; n < 10; n++)
        {
            pApm = (xtfloat*)(Af + (0 * 10 + m) * 2);
            pApn = (xtfloat*)(Af + (0 * 10 + n) * 2);

            R_re = R_im = XT_CONST_S(0);
            for (p = 0; p < 10; p++)
            {
                XT_LSIP(Apm_re, pApm, SZ_F32);
                XT_LSIP(Apm_im, pApm, (10 * 2 - 1) * SZ_F32);
                XT_LSIP(Apn_re, pApn, SZ_F32);
                XT_LSIP(Apn_im, pApn, (10 * 2 - 1) * SZ_F32);
                XT_MADD_S(R_re, Apm_re, Apn_re);
                XT_MADD_S(R_re, Apm_im, Apn_im);
                XT_MADD_S(R_im, Apm_re, Apn_im);
                XT_MSUB_S(R_im, Apm_im, Apn_re);
            }
            XT_SSIP(R_re, pRv, SZ_F32);
            XT_SSIP(R_im, pRv, SZ_F32);
            XT_SSIP(R_re, pRh, SZ_F32);
            XT_SSIP(XT_NEG_S(R_im), pRh, (10 * 2 - 1) * SZ_F32);
        }
    }

#endif
}

size_t  cmatcholpreprocess10x10f_getScratchSize()
{
    return 0;
}
#else
DISCARD_FUN(void, cmatcholpreprocess10x10f, (void* pScr,
	complex_float *R,
	const complex_float * A,
	const float32_t  sigma2))

	size_t  cmatcholpreprocess10x10f_getScratchSize()
{
	return 0;
}
#endif
