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
void  cmatcholpreprocess8x8f(void* pScr,
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
	xtfloatx4 * restrict pR = (xtfloatx4 *)R;
	xtfloatx4 * restrict pRv = pR;
	xtfloatx4 * restrict pRh = pR;
	xtfloatx2 Acc00, Acc01, Acc02, Acc03, Acc04, Acc05, Acc06, Acc07;
	xtfloatx2        Acc11, Acc12, Acc13, Acc14, Acc15, Acc16, Acc17;
	xtfloatx2               Acc22, Acc23, Acc24, Acc25, Acc26, Acc27;
	xtfloatx2                      Acc33, Acc34, Acc35, Acc36, Acc37;
	xtfloatx2                             Acc44, Acc45, Acc46, Acc47;
	xtfloatx2                                    Acc55, Acc56, Acc57;
	xtfloatx2                                           Acc66, Acc67;
	xtfloatx2                                                  Acc77;
	xtfloatx2 Apn0, Apn1, Apn2, Apn3, Apn4, Apn5, Apn6, Apn7;
	

	/* FIRST PART */
	Acc00 = Acc11 = sigma_w_zeros;
	CONST_SX2X2(Acc01, Acc02, 0);
	CONST_SX2X2(Acc03, Acc04, 0);
	CONST_SX2X2(Acc05, Acc06, 0);
	CONST_SX2X2(Acc07, Acc12, 0);
	CONST_SX2X2(Acc13, Acc14, 0);
	CONST_SX2X2(Acc15, Acc16, 0);
	Acc17 = XT_CONST_S(0);

	for (p = 0; p < 8; p++)
	{
		AE_LSX2X2_IP(Apn0, Apn1, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn2, Apn3, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn4, Apn5, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn6, Apn7, pApn, 2 * SZ_CF32);

		MADDMUX_SX2X2(Acc00, Acc01, Apn0, Apn0, Apn0, Apn1, 0);
		MADDMUX_SX2X2(Acc00, Acc01, Apn0, Apn0, Apn0, Apn1, 3);
		MADDMUX_SX2X2(Acc02, Acc03, Apn0, Apn0, Apn2, Apn3, 0);
		MADDMUX_SX2X2(Acc02, Acc03, Apn0, Apn0, Apn2, Apn3, 3);
		MADDMUX_SX2X2(Acc04, Acc05, Apn0, Apn0, Apn4, Apn5, 0);
		MADDMUX_SX2X2(Acc04, Acc05, Apn0, Apn0, Apn4, Apn5, 3);
		MADDMUX_SX2X2(Acc06, Acc07, Apn0, Apn0, Apn6, Apn7, 0);
		MADDMUX_SX2X2(Acc06, Acc07, Apn0, Apn0, Apn6, Apn7, 3);

		MADDMUX_SX2X2(Acc11, Acc12, Apn1, Apn1, Apn1, Apn2, 0);
		MADDMUX_SX2X2(Acc11, Acc12, Apn1, Apn1, Apn1, Apn2, 3);
		MADDMUX_SX2X2(Acc13, Acc14, Apn1, Apn1, Apn3, Apn4, 0);
		MADDMUX_SX2X2(Acc13, Acc14, Apn1, Apn1, Apn3, Apn4, 3);
		MADDMUX_SX2X2(Acc15, Acc16, Apn1, Apn1, Apn5, Apn6, 0);
		MADDMUX_SX2X2(Acc15, Acc16, Apn1, Apn1, Apn5, Apn6, 3);

		XT_MADDMUX_S(Acc17, Apn1, Apn7, 0);
		XT_MADDMUX_S(Acc17, Apn1, Apn7, 3);
	}
	AE_SSX2X2_IP(Acc00, Acc01, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc02, Acc03, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc04, Acc05, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc06, Acc07, pRv, 2 * SZ_CF32);

	AE_SSX2X2_IP(XT_CONJC_S(Acc01), Acc11, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc12, Acc13, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc14, Acc15, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc16, Acc17, pRv, 4 * SZ_CF32);

	pRh = (xtfloatx4 *)R + 8;
	CONJC_SX2X2(Acc02, Acc12, Acc02, Acc12);
	CONJC_SX2X2(Acc03, Acc13, Acc03, Acc13);
	CONJC_SX2X2(Acc04, Acc14, Acc04, Acc14);
	CONJC_SX2X2(Acc05, Acc15, Acc05, Acc15);
	CONJC_SX2X2(Acc06, Acc16, Acc06, Acc16);
	CONJC_SX2X2(Acc07, Acc17, Acc07, Acc17);
	AE_SSX2X2_XP(Acc02, Acc12, pRh, 8 * SZ_CF32);
	AE_SSX2X2_XP(Acc03, Acc13, pRh, 8 * SZ_CF32);
	AE_SSX2X2_XP(Acc04, Acc14, pRh, 8 * SZ_CF32);
	AE_SSX2X2_XP(Acc05, Acc15, pRh, 8 * SZ_CF32);
	AE_SSX2X2_XP(Acc06, Acc16, pRh, 8 * SZ_CF32);
	AE_SSX2X2_IP(Acc07, Acc17, pRh, 0);

	__Pragma("flush_memory");

	/* SECOND PART */
	Acc22 = Acc33 = sigma_w_zeros;
	CONST_SX2X2(Acc23, Acc24, 0);
	CONST_SX2X2(Acc25, Acc26, 0);
	CONST_SX2X2(Acc27, Acc34, 0);
	CONST_SX2X2(Acc35, Acc36, 0);

	Acc37 = XT_CONST_S(0);

	pApn = (xtfloatx4 *)A;
	for (p = 0; p < 8; p++)
	{
		pApn++;
		AE_LSX2X2_IP(Apn2, Apn3, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn4, Apn5, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn6, Apn7, pApn, 2 * SZ_CF32);

		XT_MADDMUX_S(Acc22, Apn2, Apn2, 0);
		XT_MADDMUX_S(Acc22, Apn2, Apn2, 3);

		MADDMUX_SX2X2(Acc23, Acc24, Apn2, Apn2, Apn3, Apn4, 0);
		MADDMUX_SX2X2(Acc23, Acc24, Apn2, Apn2, Apn3, Apn4, 3);
		MADDMUX_SX2X2(Acc25, Acc26, Apn2, Apn2, Apn5, Apn6, 0);
		MADDMUX_SX2X2(Acc25, Acc26, Apn2, Apn2, Apn5, Apn6, 3);

		MADDMUX_SX2X2(Acc27, Acc33, Apn2, Apn3, Apn7, Apn3, 0);
		MADDMUX_SX2X2(Acc27, Acc33, Apn2, Apn3, Apn7, Apn3, 3);
		MADDMUX_SX2X2(Acc34, Acc35, Apn3, Apn3, Apn4, Apn5, 0);
		MADDMUX_SX2X2(Acc34, Acc35, Apn3, Apn3, Apn4, Apn5, 3);
		MADDMUX_SX2X2(Acc36, Acc37, Apn3, Apn3, Apn6, Apn7, 0);
		MADDMUX_SX2X2(Acc36, Acc37, Apn3, Apn3, Apn6, Apn7, 3);
	}

	AE_SSX2X2_IP(Acc22, Acc23, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc24, Acc25, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc26, Acc27, pRv, 4 * SZ_CF32);

	AE_SSX2X2_IP(XT_CONJC_S(Acc23), Acc33, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc34, Acc35, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc36, Acc37, pRv, 6 * SZ_CF32);

	pRh = (xtfloatx4 *)R + 17;
	AE_SSX2X2_XP(XT_CONJC_S(Acc24), XT_CONJC_S(Acc34), pRh, 8 * SZ_CF32);
	AE_SSX2X2_XP(XT_CONJC_S(Acc25), XT_CONJC_S(Acc35), pRh, 8 * SZ_CF32);
	AE_SSX2X2_XP(XT_CONJC_S(Acc26), XT_CONJC_S(Acc36), pRh, 8 * SZ_CF32);
	AE_SSX2X2_IP(XT_CONJC_S(Acc27), XT_CONJC_S(Acc37), pRh, 0);

	__Pragma("flush_memory");
	/* THIRD PART */
	Acc44 = Acc55 = Acc66 = Acc77 = sigma_w_zeros;
	CONST_SX2X2(Acc45, Acc46, 0);
	CONST_SX2X2(Acc47, Acc56, 0);
	CONST_SX2X2(Acc57, Acc67, 0);

	pApn = (xtfloatx4 *)A;
	for (p = 0; p < 8; p++)
	{
		pApn += 2;
		AE_LSX2X2_IP(Apn4, Apn5, pApn, 2 * SZ_CF32);
		AE_LSX2X2_IP(Apn6, Apn7, pApn, 2 * SZ_CF32);

		MADDMUX_SX2X2(Acc44, Acc45, Apn4, Apn4, Apn4, Apn5, 0);
		MADDMUX_SX2X2(Acc44, Acc45, Apn4, Apn4, Apn4, Apn5, 3);
		MADDMUX_SX2X2(Acc46, Acc47, Apn4, Apn4, Apn6, Apn7, 0);
		MADDMUX_SX2X2(Acc46, Acc47, Apn4, Apn4, Apn6, Apn7, 3);

		MADDMUX_SX2X2(Acc55, Acc56, Apn5, Apn5, Apn5, Apn6, 0);
		MADDMUX_SX2X2(Acc55, Acc56, Apn5, Apn5, Apn5, Apn6, 3);

		MADDMUX_SX2X2(Acc57, Acc66, Apn5, Apn6, Apn7, Apn6, 0);
		MADDMUX_SX2X2(Acc57, Acc66, Apn5, Apn6, Apn7, Apn6, 3);

		MADDMUX_SX2X2(Acc67, Acc77, Apn6, Apn7, Apn7, Apn7, 0);
		MADDMUX_SX2X2(Acc67, Acc77, Apn6, Apn7, Apn7, Apn7, 3);

	}

	AE_SSX2X2_IP(Acc44, Acc45, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc46, Acc47, pRv, 6 * SZ_CF32);

	CONJC_SX2X2(Acc45, Acc46, Acc45, Acc46);

	AE_SSX2X2_IP(Acc45, Acc55, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc56, Acc57, pRv, 6 * SZ_CF32);

	CONJC_SX2X2(Acc56, Acc47, Acc56, Acc47);

	AE_SSX2X2_IP(Acc46, Acc56, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc66, Acc67, pRv, 6 * SZ_CF32);

	CONJC_SX2X2(Acc57, Acc67, Acc57, Acc67);

	AE_SSX2X2_IP(Acc47, Acc57, pRv, 2 * SZ_CF32);
	AE_SSX2X2_IP(Acc67, Acc77, pRv, 2 * SZ_CF32);
}

size_t  cmatcholpreprocess8x8f_getScratchSize()
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
void  cmatcholpreprocess8x8f(void* pScr,
    complex_float *R,
    const complex_float * A,
    const float32_t  sigma2)
{
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
    for (m = 0; m < 8; m += 2)
    {
        /* diagonal squares */
        re_acc00 = re_acc11 = sigma2;
        re_acc01 = re_acc10 = 0.f;
        im_acc00 = im_acc01 = im_acc10 = im_acc11 = 0.f;

        for (p = 0; p < 8; p++)
        {
            re_Apm0 = Af[(p * 8 + (m + 0)) * 2 + 0];
            re_Apm1 = Af[(p * 8 + (m + 1)) * 2 + 0];
            im_Apm0 = Af[(p * 8 + (m + 0)) * 2 + 1];
            im_Apm1 = Af[(p * 8 + (m + 1)) * 2 + 1];

            re_acc00 += re_Apm0*re_Apm0 + im_Apm0*im_Apm0;
            re_acc01 += re_Apm0*re_Apm1 + im_Apm0*im_Apm1;
            re_acc11 += re_Apm1*re_Apm1 + im_Apm1*im_Apm1;
            im_acc00 += re_Apm0*im_Apm0 - im_Apm0*re_Apm0;
            im_acc01 += re_Apm0*im_Apm1 - im_Apm0*re_Apm1;
            im_acc11 += re_Apm1*im_Apm1 - im_Apm1*re_Apm1;
        }
        Rf[((m + 0) * 8 + (m + 0)) * 2 + 0] = re_acc00;
        Rf[((m + 0) * 8 + (m + 0)) * 2 + 1] = im_acc00;
        Rf[((m + 0) * 8 + (m + 1)) * 2 + 0] = re_acc01;
        Rf[((m + 0) * 8 + (m + 1)) * 2 + 1] = im_acc01;
        Rf[((m + 1) * 8 + (m + 0)) * 2 + 0] = re_acc01;
        Rf[((m + 1) * 8 + (m + 0)) * 2 + 1] = -im_acc01;
        Rf[((m + 1) * 8 + (m + 1)) * 2 + 0] = re_acc11;
        Rf[((m + 1) * 8 + (m + 1)) * 2 + 1] = im_acc11;

        /* rest */
        for (n = m + 2; n < 8; n += 2)
        {
            re_acc00 = re_acc11 = re_acc01 = re_acc10 = 0.f;
            im_acc00 = im_acc01 = im_acc10 = im_acc11 = 0.f;

            for (p = 0; p < 8; p++)
            {
                re_Apm0 = Af[(p * 8 + (m + 0)) * 2 + 0];
                re_Apm1 = Af[(p * 8 + (m + 1)) * 2 + 0];
                re_Apn0 = Af[(p * 8 + (n + 0)) * 2 + 0];
                re_Apn1 = Af[(p * 8 + (n + 1)) * 2 + 0];
                im_Apm0 = Af[(p * 8 + (m + 0)) * 2 + 1];
                im_Apm1 = Af[(p * 8 + (m + 1)) * 2 + 1];
                im_Apn0 = Af[(p * 8 + (n + 0)) * 2 + 1];
                im_Apn1 = Af[(p * 8 + (n + 1)) * 2 + 1];

                re_acc00 += re_Apm0*re_Apn0 + im_Apm0*im_Apn0;
                re_acc01 += re_Apm0*re_Apn1 + im_Apm0*im_Apn1;
                re_acc10 += re_Apm1*re_Apn0 + im_Apm1*im_Apn0;
                re_acc11 += re_Apm1*re_Apn1 + im_Apm1*im_Apn1;
                im_acc00 += re_Apm0*im_Apn0 - im_Apm0*re_Apn0;
                im_acc01 += re_Apm0*im_Apn1 - im_Apm0*re_Apn1;
                im_acc10 += re_Apm1*im_Apn0 - im_Apm1*re_Apn0;
                im_acc11 += re_Apm1*im_Apn1 - im_Apm1*re_Apn1;
            }
            Rf[((m + 0) * 8 + (n + 0)) * 2 + 0] = re_acc00;
            Rf[((m + 0) * 8 + (n + 0)) * 2 + 1] = im_acc00;
            Rf[((m + 0) * 8 + (n + 1)) * 2 + 0] = re_acc01;
            Rf[((m + 0) * 8 + (n + 1)) * 2 + 1] = im_acc01;
            Rf[((m + 1) * 8 + (n + 0)) * 2 + 0] = re_acc10;
            Rf[((m + 1) * 8 + (n + 0)) * 2 + 1] = im_acc10;
            Rf[((m + 1) * 8 + (n + 1)) * 2 + 0] = re_acc11;
            Rf[((m + 1) * 8 + (n + 1)) * 2 + 1] = im_acc11;

            Rf[((m + 0) + (n + 0) * 8) * 2 + 0] = re_acc00;
            Rf[((m + 0) + (n + 0) * 8) * 2 + 1] = -im_acc00;
            Rf[((m + 0) + (n + 1) * 8) * 2 + 0] = re_acc01;
            Rf[((m + 0) + (n + 1) * 8) * 2 + 1] = -im_acc01;
            Rf[((m + 1) + (n + 0) * 8) * 2 + 0] = re_acc10;
            Rf[((m + 1) + (n + 0) * 8) * 2 + 1] = -im_acc10;
            Rf[((m + 1) + (n + 1) * 8) * 2 + 0] = re_acc11;
            Rf[((m + 1) + (n + 1) * 8) * 2 + 1] = -im_acc11;

        }
    }
}

size_t  cmatcholpreprocess8x8f_getScratchSize()
{
    return 0;
}
#else
DISCARD_FUN(void, cmatcholpreprocess8x8f, (void* pScr,
	complex_float *R,
	const complex_float * A,
	const float32_t  sigma2))

	size_t  cmatcholpreprocess8x8f_getScratchSize()
{
	return 0;
}
#endif
