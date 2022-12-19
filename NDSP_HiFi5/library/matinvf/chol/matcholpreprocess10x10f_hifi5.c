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

/* Load 32bit and replicate to xtfloatx2*/
#define _L32_SX2_IP(a,b,c) \
{ \
ae_int32x2 tmp; \
AE_L32_IP(tmp, castxcc(ae_int32, b), c); \
a = XT_AE_MOVXTFLOATX2_FROMINT32X2(tmp); }
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
void  matcholpreprocess10x10f(void* pScr,
    float32_t *R,
    const float32_t * A,
    const float32_t  sigma2)
{
	NASSERT(R);
	NASSERT(A);
	NASSERT(sigma2);
	NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
	NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);


	const xtfloat * restrict pApn = (xtfloat *)A;
	int p;
	xtfloatx2 * restrict pRv = (xtfloatx2 *)R;
	xtfloatx2 * restrict pRh;
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
	xtfloatx2 Accw0, Accw1;
	Acc00 = Acc11 = sigma2;
	Acc01 = Acc02 = Acc03 = Acc04 = Acc05 = Acc06 = Acc07 = Acc08 = Acc09 = XT_CONST_S(0);
	Acc12 = Acc13 = Acc14 = Acc15 = Acc16 = Acc17 = Acc18 = Acc19 = XT_CONST_S(0);

	/* FIRST PART */
	for (p = 0; p < 10; p++)
	{
		_L32_SX2_IP(Apn0, pApn, SZ_F32);
		_L32_SX2_IP(Apn1, pApn, SZ_F32);
		_L32_SX2_IP(Apn2, pApn, SZ_F32);
		_L32_SX2_IP(Apn3, pApn, SZ_F32);
		_L32_SX2_IP(Apn4, pApn, SZ_F32);
		_L32_SX2_IP(Apn5, pApn, SZ_F32);
		_L32_SX2_IP(Apn6, pApn, SZ_F32);
		_L32_SX2_IP(Apn7, pApn, SZ_F32);
		_L32_SX2_IP(Apn8, pApn, SZ_F32);
		_L32_SX2_IP(Apn9, pApn, SZ_F32);

		MADD_SX2X2(Acc00, Acc01, Apn0, Apn0, Apn0, Apn1);
		MADD_SX2X2(Acc02, Acc03, Apn0, Apn0, Apn2, Apn3);
		MADD_SX2X2(Acc04, Acc05, Apn0, Apn0, Apn4, Apn5);
		MADD_SX2X2(Acc06, Acc07, Apn0, Apn0, Apn6, Apn7);
		MADD_SX2X2(Acc08, Acc09, Apn0, Apn0, Apn8, Apn9);

		MADD_SX2X2(Acc11, Acc12, Apn1, Apn1, Apn1, Apn2);
		MADD_SX2X2(Acc13, Acc14, Apn1, Apn1, Apn3, Apn4);
		MADD_SX2X2(Acc15, Acc16, Apn1, Apn1, Apn5, Apn6);
		MADD_SX2X2(Acc17, Acc18, Apn1, Apn1, Apn7, Apn8);
		MADD_SX2(Acc19, Apn1, Apn9);
	}
	Accw0 = AE_SELSX2IR(Acc01, Acc00, 1);
	Accw1 = AE_SELSX2IR(Acc03, Acc02, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc05, Acc04, 1);
	Accw1 = AE_SELSX2IR(Acc07, Acc06, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc09, Acc08, 1);
	Accw1 = AE_SELSX2IR(Acc11, Acc01, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc13, Acc12, 1);
	Accw1 = AE_SELSX2IR(Acc15, Acc14, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc17, Acc16, 1);
	Accw1 = AE_SELSX2IR(Acc19, Acc18, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 4 * SZ_F32);


	pRh = (xtfloatx2 *)R + 10;
	Accw0 = AE_SELSX2IR(Acc12, Acc02, 1);
	AE_SSX2IP(Accw0, pRh, 10 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc13, Acc03, 1);
	AE_SSX2IP(Accw0, pRh, 10 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc14, Acc04, 1);
	AE_SSX2IP(Accw0, pRh, 10 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc15, Acc05, 1);
	AE_SSX2IP(Accw0, pRh, 10 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc16, Acc06, 1);
	AE_SSX2IP(Accw0, pRh, 10 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc17, Acc07, 1);
	AE_SSX2IP(Accw0, pRh, 10 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc18, Acc08, 1);
	AE_SSX2IP(Accw0, pRh, 10 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc19, Acc09, 1);
	AE_SSX2IP(Accw0, pRh, 10 * SZ_F32);

	__Pragma("flush_memory");

	/* SECOND PART */
	Acc22 = Acc33 = sigma2;
	Acc23 = Acc24 = Acc25 = Acc26 = Acc27 = Acc28 = Acc29 = XT_CONST_S(0);
	Acc34 = Acc35 = Acc36 = Acc37 = Acc38 = Acc39 = XT_CONST_S(0);

	pApn = (xtfloat *)A;
	for (p = 0; p < 10; p++)
	{
		pApn+=2;
		_L32_SX2_IP(Apn2, pApn, SZ_F32);
		_L32_SX2_IP(Apn3, pApn, SZ_F32);
		_L32_SX2_IP(Apn4, pApn, SZ_F32);
		_L32_SX2_IP(Apn5, pApn, SZ_F32);
		_L32_SX2_IP(Apn6, pApn, SZ_F32);
		_L32_SX2_IP(Apn7, pApn, SZ_F32);
		_L32_SX2_IP(Apn8, pApn, SZ_F32);
		_L32_SX2_IP(Apn9, pApn, SZ_F32);

		MADD_SX2(Acc22, Apn2, Apn2);

		MADD_SX2X2(Acc23, Acc24, Apn2, Apn2, Apn3, Apn4);
		MADD_SX2X2(Acc25, Acc26, Apn2, Apn2, Apn5, Apn6);
		MADD_SX2X2(Acc27, Acc28, Apn2, Apn2, Apn7, Apn8);
		
		MADD_SX2X2(Acc29, Acc33, Apn2, Apn3, Apn9, Apn3);
		MADD_SX2X2(Acc34, Acc35, Apn3, Apn3, Apn4, Apn5);
		MADD_SX2X2(Acc36, Acc37, Apn3, Apn3, Apn6, Apn7);
		MADD_SX2X2(Acc38, Acc39, Apn3, Apn3, Apn8, Apn9);
	}
	Accw0 = AE_SELSX2IR(Acc23, Acc22, 1);
	Accw1 = AE_SELSX2IR(Acc25, Acc24, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc27, Acc26, 1);
	Accw1 = AE_SELSX2IR(Acc29, Acc28, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 4 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc33, Acc23, 1);
	Accw1 = AE_SELSX2IR(Acc35, Acc34, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc37, Acc36, 1);
	Accw1 = AE_SELSX2IR(Acc39, Acc38, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 6 * SZ_F32);

	pRh = (xtfloatx2 *)R + 21;
	Accw0 = AE_SELSX2IR(Acc34, Acc24, 1);
	AE_SSX2IP(Accw0, pRh, 10 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc35, Acc25, 1);
	AE_SSX2IP(Accw0, pRh, 10 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc36, Acc26, 1);
	AE_SSX2IP(Accw0, pRh, 10 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc37, Acc27, 1);
	AE_SSX2IP(Accw0, pRh, 10 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc38, Acc28, 1);
	AE_SSX2IP(Accw0, pRh, 10 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc39, Acc29, 1);
	AE_SSX2IP(Accw0, pRh, 10 * SZ_F32);

	__Pragma("flush_memory");

	/* THIRD PART */
	Acc44 = Acc55 = Acc66 = Acc77 = Acc88 = Acc99 = sigma2;
	Acc45 = Acc46 = Acc47 = Acc48 = Acc49 = XT_CONST_S(0);
	Acc56 = Acc57 = Acc58 = Acc59 = XT_CONST_S(0);
	Acc67 = Acc68 = Acc69 = XT_CONST_S(0);
	Acc78 = Acc79 = XT_CONST_S(0);
	Acc89 = XT_CONST_S(0);

	pApn = (xtfloat *)A;
	for (p = 0; p < 10; p++)
	{
		pApn += 4;
		_L32_SX2_IP(Apn4, pApn, SZ_F32);
		_L32_SX2_IP(Apn5, pApn, SZ_F32);
		_L32_SX2_IP(Apn6, pApn, SZ_F32);
		_L32_SX2_IP(Apn7, pApn, SZ_F32);
		_L32_SX2_IP(Apn8, pApn, SZ_F32);
		_L32_SX2_IP(Apn9, pApn, SZ_F32);

		MADD_SX2X2(Acc44, Acc45, Apn4, Apn4, Apn4, Apn5);
		MADD_SX2X2(Acc46, Acc47, Apn4, Apn4, Apn6, Apn7);
		MADD_SX2X2(Acc48, Acc49, Apn4, Apn4, Apn8, Apn9);

		MADD_SX2X2(Acc55, Acc56, Apn5, Apn5, Apn5, Apn6);
		MADD_SX2X2(Acc57, Acc58, Apn5, Apn5, Apn7, Apn8);
		
		MADD_SX2X2(Acc59, Acc66, Apn5, Apn6, Apn9, Apn6);
		MADD_SX2X2(Acc67, Acc68, Apn6, Apn6, Apn7, Apn8);
		
		MADD_SX2X2(Acc69, Acc77, Apn6, Apn7, Apn9, Apn7);
		MADD_SX2X2(Acc78, Acc79, Apn7, Apn7, Apn8, Apn9);
		
		MADD_SX2X2(Acc88, Acc89, Apn8, Apn8, Apn8, Apn9);

		MADD_SX2(Acc99, Apn9, Apn9);
	}
	Accw0 = AE_SELSX2IR(Acc45, Acc44, 1);
	Accw1 = AE_SELSX2IR(Acc47, Acc46, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc49, Acc48, 1);
	Accw1 = AE_SELSX2IR(Acc55, Acc45, 1);
	AE_SSX2IP(Accw0, pRv, 6 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc57, Acc56, 1);
	Accw1 = AE_SELSX2IR(Acc59, Acc58, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 6 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc56, Acc46, 1);
	Accw1 = AE_SELSX2IR(Acc67, Acc66, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc69, Acc68, 1);
	Accw1 = AE_SELSX2IR(Acc57, Acc47, 1);
	AE_SSX2IP(Accw0, pRv, 6 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc77, Acc67, 1);
	Accw1 = AE_SELSX2IR(Acc79, Acc78, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 6 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc58, Acc48, 1);
	Accw1 = AE_SELSX2IR(Acc78, Acc68, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc89, Acc88, 1);
	Accw1 = AE_SELSX2IR(Acc59, Acc49, 1);
	AE_SSX2IP(Accw0, pRv, 6 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc79, Acc69, 1);
	Accw1 = AE_SELSX2IR(Acc99, Acc89, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 6 * SZ_F32);

}
/* scratch allocation functions */
size_t   matcholpreprocess10x10f_getScratchSize()
{
	return 0;
}
#elif(HAVE_FPU)
#define SZ_F32 (sizeof(float32_t))
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
/* single-matrix API */
void  matcholpreprocess10x10f(void* restrict pScr,
    float32_t * restrict R,
    const float32_t * restrict A,
    const float32_t  sigma2)
{
#if 0
    int m, n, p;

    NASSERT(R);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);
    float32_t acc;
    float32_t Apm;
    float32_t Apn;
    for (m = 0; m < 10; m++)
    {
        acc = sigma2;
        for (p = 0; p < 10; p++)
        {
            Apm = A[p * 10 + m];

            acc += Apm*Apm;
        }
        R[m * 10 + m] = acc;

        for (n = m + 1; n < 10; n++)
        {
            acc = 0.f;
            for (p = 0; p < 10; p++)
            {
                Apm = A[p * 10 + m];
                Apn = A[p * 10 + n];

                acc += Apm*Apn;
            }
            R[m * 10 + n] = acc;
            R[m + n * 10] = acc;
        }
    }
#else
    int m, n, p;

    NASSERT(R);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);
    float32_t acc00, acc01, acc10, acc11;
    float32_t Apm0, Apm1;
    float32_t Apn0, Apn1;
    /* processing by squares 2x2 */
    for (m = 0; m < 10; m+=2)
    {
        /* diagonal squares */
        acc00 = acc11 = sigma2;
        acc01 = 0.f;
        for (p = 0; p < 10; p++)
        {
            Apm0 = A[p * 10 + m + 0];
            Apm1 = A[p * 10 + m + 1];
            acc00 += Apm0*Apm0;
            acc01 += Apm0*Apm1;
            acc11 += Apm1*Apm1;
        }
        R[(m + 0) * 10 + (m + 0)] = acc00;
        R[(m + 0) * 10 + (m + 1)] = acc01;
        R[(m + 1) * 10 + (m + 0)] = acc01;
        R[(m + 1) * 10 + (m + 1)] = acc11;

        /* rest */
        for (n = m+2; n < 10; n+=2)
        {
            acc00 = acc11 = 0.f;
            acc01 = acc10 = 0.f;
            for (p = 0; p < 10; p++)
            {
                Apm0 = A[p * 10 + m + 0];
                Apm1 = A[p * 10 + m + 1];
                Apn0 = A[p * 10 + n + 0];
                Apn1 = A[p * 10 + n + 1];
                acc00 += Apm0*Apn0;
                acc01 += Apm0*Apn1;
                acc10 += Apm1*Apn0;
                acc11 += Apm1*Apn1;
            }
            R[(m + 0) * 10 + (n + 0)] = acc00;
            R[(m + 0) * 10 + (n + 1)] = acc01;
            R[(m + 1) * 10 + (n + 0)] = acc10;
            R[(m + 1) * 10 + (n + 1)] = acc11;

            R[(m + 0) + (n + 0) * 10] = acc00;
            R[(m + 0) + (n + 1) * 10] = acc01;
            R[(m + 1) + (n + 0) * 10] = acc10;
            R[(m + 1) + (n + 1) * 10] = acc11;
        }
    }
#endif
}

size_t  matcholpreprocess10x10f_getScratchSize()
{
    return 0;
}
#else
DISCARD_FUN(void, matcholpreprocess10x10f, (void* pScr,
	float32_t *R,
	const float32_t * A,
	const float32_t  sigma2))

	size_t  matcholpreprocess10x10f_getScratchSize()
{
	return 0;
}
#endif
