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
#include "cholnf_common.h"
/*
code optimized for HiFi4 with VFPU
*/

#if (HAVE_VFPU || HAVE_FPU)
#define SZ_F32 (int)(sizeof(float32_t))
/*-------------------------------------------------------------------------
Cholesky decomposition for Pseudo-Inversion
Apply the Cholesky decomposition to the matrix of normal equations system
associated with a complex-valued least squares problem:
A[MxN]*X[NxP]=B[MxP], MxN.
The decomposition results in an upper triangular complex matrix R[NxN]
with real and positive numbers on the main diagonal, such that
R'*R = A'*A + sigma2*I[NxN], where ' denotes the conjugate transpose of
a matrix, and sigma2*I[NxN] is the NxN identity matrix multiplied with
the regularization term.

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
A[M*N]          matrix A. For fixed point API, the representation is Q(qA)
sigma2          regularization term. For fixed point, the representation 
                should be Q(2*qA-30)
qRA             qR-qA; difference between fixed point representations of R
                and A matrices (for the fixed point API only). Should be 
                equal or less than 0 (typically -2).
Output:
R[((N+1)*N)/2]  upper triangular matrix R
D[N]            reciprocals of the main diagonal. NOTE: for the fixed point API,
                these data are stored internally in special format with separate
                mantissa and exponent for better accuracy and dynamic range 
                control. So, even for the real data, they stored as pairs of 2
                integers and packed to the complex_fract32 format
Temporary:
pScr            Scratch memory

N = M = 4, 6, 8, 10

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
16-bytes boundary
---------------------------------------------------------------------------*/
void  matcholdecomp4x4f(void* pScr,
    float32_t * R,
    float32_t * D,
    const float32_t * A,
    const float32_t   sigma2)
{
    static const f_realcholFwdrecf f_Fwdrec[4] = { realcholFwdrec1f, realcholFwdrec2f, realcholFwdrec3f, realcholFwdrec4f };
    int n;
    float32_t *y; // pointer to the new column in R
    float32_t* prep = (float32_t*)pScr;

    NASSERT(R);
    NASSERT(D);
    NASSERT(A);
    NASSERT(sigma2);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(sigma2, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    y = R;

    // compute  A'*A 
    matcholpreprocess4x4f(NULL, prep, A, sigma2);

#if (HAVE_VFPU)
    /* to avoid ferret warnings */
    {
      xtfloatx4 *T0 = (xtfloatx4 *)R;
      xtfloatx4 *T1 = (xtfloatx4 *)D;

      AE_SSX2X2_IP(0, 0, T1, 4 * (sizeof(float32_t)));

      for (n = 0; n<((4 * 5) >> 3); n++)
      {
        AE_SSX2X2_IP(0, 0, T0, 4 * (sizeof(float32_t)));
      }
      XT_SSX2IP(0, castxcc(xtfloatx2, T0), 2 * (sizeof(float32_t)));

    }
#else
    {
      float32_t *T0 = (float32_t *)R;
      float32_t *T1 = (float32_t *)D;
      for (n = 0; n<(4); n++)
      {
        T1[n] = 0.0f;
      }
      for (n = 0; n<(5 * 4 / 2); n++)
      {
        T0[n] = 0.0f;
      }
    }
#endif

    for (n = 0; n<4; n++)
    {
        y += n; // go to the next column
        // make forward recursion to update n new column elements
        f_Fwdrec[n](y, R, D, prep + n, 4);
        // update n-th diagonal element
        realcholDiagUpdf(y, D + n, prep + n * 4 + n, n);
    }
}

/* scratch allocation functions */
size_t  matcholdecomp4x4f_getScratchSize()
{
    return 4 * 4 * SZ_F32;
}
#else
DISCARD_FUN(void, matcholdecomp4x4f, (void * pScr,
	float32_t * R,
	float32_t * D,
	const float32_t * A,
	const float32_t sigma2))

	size_t  matcholdecomp4x4f_getScratchSize()
{
	return 0;
}
#endif
