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
#include "chol32x32_common.h"

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
void matcholdecomp10x10_32x32  (void * restrict pScr,
                    int32_t * restrict R,
                    complex_fract32 * restrict D,
              const int32_t * restrict A,
                    int32_t sigma2,
                    int qRA)
{
    const int N = 10;
    int64_t* Z = (int64_t*)pScr; // L columns of A'*A
    int64_t* Zdiag = Z;
    int32_t* y; // pointer to the new column in R
    int n;

    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A, HIFI_SIMD_WIDTH);

    matcholpreprocess10x10_32x32(pScr, Z, A, sigma2, qRA);

    /* to avoid ferret warnings */
    {
    
      ae_int32x4 *T0 = (ae_int32x4 *)R;
      ae_int32x4 *T1 = (ae_int32x4 *)D;

      AE_S32X2X2_IP(AE_ZERO32(), AE_ZERO32(), T1, sizeof(ae_int32x4));
      AE_S32X2X2_IP(AE_ZERO32(), AE_ZERO32(), T1, sizeof(ae_int32x4));
      AE_S32X2_IP(AE_ZERO32(), castxcc(ae_int32x2, T1), sizeof(ae_int32x2));
      
      for (n = 0; n<((10 * 11) >> 3); n++)
      {
        AE_S32X2X2_IP(AE_ZERO32(), AE_ZERO32(), T0, sizeof(ae_int32x4));
      }
      AE_S32X2_IP(AE_ZERO32(), castxcc(ae_int32x2, T0), sizeof(ae_int32x2));
      AE_S32_H_IP(AE_ZERO32(), castxcc(ae_int32, T0), 0);
    }

    y = R;
    for (n = 0; n < N; n++)
    {
        y += n; // go to the next column
        // make forward recursion to update n new column elements
        matcholfwdrec(y, R, (int32_t*)D, Z, n, N);
        // update n-th diagonal element
        matcholdiagUpd(y, (int32_t*)(D + n), Zdiag, n);
        Z += 1;
        Zdiag += N + 1;
    }
}

size_t    matcholdecomp10x10_32x32_getScratchSize()   { return (10 * 10 * sizeof(int64_t)); }
