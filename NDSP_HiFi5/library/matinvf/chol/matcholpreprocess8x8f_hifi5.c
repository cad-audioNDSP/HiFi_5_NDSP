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
void  matcholpreprocess8x8f(void* pScr,
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
	xtfloatx2 Acc00, Acc01, Acc02, Acc03, Acc04, Acc05, Acc06, Acc07;
	xtfloatx2        Acc11, Acc12, Acc13, Acc14, Acc15, Acc16, Acc17;
	xtfloatx2               Acc22, Acc23, Acc24, Acc25, Acc26, Acc27;
	xtfloatx2                      Acc33, Acc34, Acc35, Acc36, Acc37;
	xtfloatx2                             Acc44, Acc45, Acc46, Acc47;
	xtfloatx2                                    Acc55, Acc56, Acc57;
	xtfloatx2                                           Acc66, Acc67;
	xtfloatx2                                                  Acc77;
	xtfloatx2 Apn0, Apn1, Apn2, Apn3, Apn4, Apn5, Apn6, Apn7;
	xtfloatx2 Accw0, Accw1;


	/* FIRST PART */
	Acc00 = Acc11 = sigma2;
	Acc01 = Acc02 = Acc03 = Acc04 = Acc05 = Acc06 = Acc07 = XT_CONST_S(0);
	Acc12 = Acc13 = Acc14 = Acc15 = Acc16 = Acc17 = XT_CONST_S(0);
	for (p = 0; p < 8; p++)
	{
		_L32_SX2_IP(Apn0, pApn, SZ_F32);
		_L32_SX2_IP(Apn1, pApn, SZ_F32);
		_L32_SX2_IP(Apn2, pApn, SZ_F32);
		_L32_SX2_IP(Apn3, pApn, SZ_F32);
		_L32_SX2_IP(Apn4, pApn, SZ_F32);
		_L32_SX2_IP(Apn5, pApn, SZ_F32);
		_L32_SX2_IP(Apn6, pApn, SZ_F32);
		_L32_SX2_IP(Apn7, pApn, SZ_F32);

		MADD_SX2X2(Acc00, Acc01, Apn0, Apn0, Apn0, Apn1);
		MADD_SX2X2(Acc02, Acc03, Apn0, Apn0, Apn2, Apn3);
		MADD_SX2X2(Acc04, Acc05, Apn0, Apn0, Apn4, Apn5);
		MADD_SX2X2(Acc06, Acc07, Apn0, Apn0, Apn6, Apn7);

		MADD_SX2X2(Acc11, Acc12, Apn1, Apn1, Apn1, Apn2);
		MADD_SX2X2(Acc13, Acc14, Apn1, Apn1, Apn3, Apn4);
		MADD_SX2X2(Acc15, Acc16, Apn1, Apn1, Apn5, Apn6);

		MADD_SX2(Acc17, Apn1, Apn7);
	}

	Accw0 = AE_SELSX2IR(Acc01, Acc00, 1);
	Accw1 = AE_SELSX2IR(Acc03, Acc02, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);
	
	Accw0 = AE_SELSX2IR(Acc05, Acc04, 1);
	Accw1 = AE_SELSX2IR(Acc07, Acc06, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc11, Acc01, 1);
	Accw1 = AE_SELSX2IR(Acc13, Acc12, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);
	
	Accw0 = AE_SELSX2IR(Acc15, Acc14, 1);
	Accw1 = AE_SELSX2IR(Acc17, Acc16, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 4 * SZ_F32);

	pRh = (xtfloatx2 *)R + 8;
	Accw0 = AE_SELSX2IR(Acc12, Acc02, 1);
	AE_SSX2IP(Accw0, pRh, 8 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc13, Acc03, 1);
	AE_SSX2IP(Accw0, pRh, 8 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc14, Acc04, 1);
	AE_SSX2IP(Accw0, pRh, 8 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc15, Acc05, 1);
	AE_SSX2IP(Accw0, pRh, 8 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc16, Acc06, 1);
	AE_SSX2IP(Accw0, pRh, 8 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc17, Acc07, 1);
	AE_SSX2IP(Accw0, pRh, 8 * SZ_F32);

	__Pragma("flush_memory");

	/* SECOND PART */
	Acc22 = Acc33 = sigma2;
	Acc23 = Acc24 = Acc25 = Acc26 = Acc27 = XT_CONST_S(0);
	Acc34 = Acc35 = Acc36 = Acc37 = XT_CONST_S(0);

	pApn = (xtfloat *)A;
	for (p = 0; p < 8; p++)
	{
		pApn+=2;
		_L32_SX2_IP(Apn2, pApn, SZ_F32);
		_L32_SX2_IP(Apn3, pApn, SZ_F32);
		_L32_SX2_IP(Apn4, pApn, SZ_F32);
		_L32_SX2_IP(Apn5, pApn, SZ_F32);
		_L32_SX2_IP(Apn6, pApn, SZ_F32);
		_L32_SX2_IP(Apn7, pApn, SZ_F32);

		MADD_SX2(Acc22, Apn2, Apn2);

		MADD_SX2X2(Acc23, Acc24, Apn2, Apn2, Apn3, Apn4);
		MADD_SX2X2(Acc25, Acc26, Apn2, Apn2, Apn5, Apn6);
		
		MADD_SX2X2(Acc27, Acc33, Apn2, Apn3, Apn7, Apn3);
		MADD_SX2X2(Acc34, Acc35, Apn3, Apn3, Apn4, Apn5);
		MADD_SX2X2(Acc36, Acc37, Apn3, Apn3, Apn6, Apn7);
	}
	Accw0 = AE_SELSX2IR(Acc23, Acc22, 1);
	Accw1 = AE_SELSX2IR(Acc25, Acc24, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc27, Acc26, 1);
	Accw1 = AE_SELSX2IR(Acc33, Acc23, 1);
	AE_SSX2IP(Accw0, pRv, 4 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc35, Acc34, 1);
	Accw1 = AE_SELSX2IR(Acc37, Acc36, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 6 * SZ_F32);

	pRh = (xtfloatx2 *)R + 17;
	Accw0 = AE_SELSX2IR(Acc34, Acc24, 1);
	AE_SSX2IP(Accw0, pRh, 8 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc35, Acc25, 1);
	AE_SSX2IP(Accw0, pRh, 8 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc36, Acc26, 1);
	AE_SSX2IP(Accw0, pRh, 8 * SZ_F32);
	Accw0 = AE_SELSX2IR(Acc37, Acc27, 1);
	AE_SSX2IP(Accw0, pRh, 8 * SZ_F32);

	__Pragma("flush_memory");
	/* THIRD PART */
	Acc44 = Acc55 = Acc66 = Acc77 = sigma2;
	Acc45 = Acc46 = Acc47 = XT_CONST_S(0);
	Acc56 = Acc57 = XT_CONST_S(0);
	Acc67 = XT_CONST_S(0);

	pApn = (xtfloat *)A;
	for (p = 0; p < 8; p++)
	{
		pApn += 4;
		_L32_SX2_IP(Apn4, pApn, SZ_F32);
		_L32_SX2_IP(Apn5, pApn, SZ_F32);
		_L32_SX2_IP(Apn6, pApn, SZ_F32);
		_L32_SX2_IP(Apn7, pApn, SZ_F32);

		MADD_SX2X2(Acc44, Acc45, Apn4, Apn4, Apn4, Apn5);
		MADD_SX2X2(Acc46, Acc47, Apn4, Apn4, Apn6, Apn7);
		
		MADD_SX2X2(Acc55, Acc56, Apn5, Apn5, Apn5, Apn6);
		
		MADD_SX2X2(Acc57, Acc66, Apn5, Apn6, Apn7, Apn6);
		
		MADD_SX2X2(Acc67, Acc77, Apn6, Apn7, Apn7, Apn7);

	}
	Accw0 = AE_SELSX2IR(Acc45, Acc44, 1);
	Accw1 = AE_SELSX2IR(Acc47, Acc46, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 6 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc55, Acc45, 1);
	Accw1 = AE_SELSX2IR(Acc57, Acc56, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 6 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc56, Acc46, 1);
	Accw1 = AE_SELSX2IR(Acc67, Acc66, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 6 * SZ_F32);

	Accw0 = AE_SELSX2IR(Acc57, Acc47, 1);
	Accw1 = AE_SELSX2IR(Acc77, Acc67, 1);
	AE_SSX2IP(Accw0, pRv, 2 * SZ_F32);
	AE_SSX2IP(Accw1, pRv, 2 * SZ_F32);
}

/* scratch allocation functions */
size_t   matcholpreprocess8x8f_getScratchSize()
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
void  matcholpreprocess8x8f(void* pScr,
    float32_t * R,
    const float32_t * A,
    const float32_t  sigma2)
{
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
    for (m = 0; m < 8; m += 2)
    {
        acc00 = acc11 = sigma2;
        acc01 = 0.f;
        /* diagonal squares */
        for (p = 0; p < 8; p++)
        {
            Apm0 = A[p * 8 + m + 0];
            Apm1 = A[p * 8 + m + 1];
            acc00 += Apm0*Apm0;
            acc01 += Apm0*Apm1;
            acc11 += Apm1*Apm1;
        }
        R[(m + 0) * 8 + (m + 0)] = acc00;
        R[(m + 0) * 8 + (m + 1)] = acc01;
        R[(m + 1) * 8 + (m + 0)] = acc01;
        R[(m + 1) * 8 + (m + 1)] = acc11;

        /* rest */
        for (n = m + 2; n < 8; n += 2)
        {
            acc00 = acc11 = 0.f;
            acc01 = acc10 = 0.f;
            for (p = 0; p < 8; p++)
            {
                Apm0 = A[p * 8 + m + 0];
                Apm1 = A[p * 8 + m + 1];
                Apn0 = A[p * 8 + n + 0];
                Apn1 = A[p * 8 + n + 1];
                acc00 += Apm0*Apn0;
                acc01 += Apm0*Apn1;
                acc10 += Apm1*Apn0;
                acc11 += Apm1*Apn1;
            }
            R[(m + 0) * 8 + (n + 0)] = acc00;
            R[(m + 0) * 8 + (n + 1)] = acc01;
            R[(m + 1) * 8 + (n + 0)] = acc10;
            R[(m + 1) * 8 + (n + 1)] = acc11;

            R[(m + 0) + (n + 0) * 8] = acc00;
            R[(m + 0) + (n + 1) * 8] = acc01;
            R[(m + 1) + (n + 0) * 8] = acc10;
            R[(m + 1) + (n + 1) * 8] = acc11;
        }
    }
}

size_t  matcholpreprocess8x8f_getScratchSize()
{
    return 0;
}
#else
DISCARD_FUN(void, matcholpreprocess8x8f, (void* pScr,
	float32_t *R,
	const float32_t * A,
	const float32_t  sigma2))

	size_t  matcholpreprocess8x8f_getScratchSize()
{
	return 0;
}
#endif
