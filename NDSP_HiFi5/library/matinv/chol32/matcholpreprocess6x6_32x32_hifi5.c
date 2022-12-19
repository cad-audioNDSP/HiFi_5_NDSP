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
void  matcholpreprocess6x6_32x32(void* pScr,
            int64_t * Z,
    const   int32_t * A,
            int32_t  sigma2, int qRA)
{
    ae_int64x2 * restrict pZ0;
    ae_int64x2 * restrict pZ1;
    const ae_int32 * restrict pAkm;
    const ae_int32 * restrict pAkn;
    ae_int32x2 akm0, akm1;
    ae_int32x2 akn0, akn1;
    ae_int64 B_re00, B_re01, B_re10, B_re11;
    ae_ep ep_re00, ep_re01, ep_re10, ep_re11;
    int m, k, n;
#define N 6
    (void)pScr;
    for (m = 0; m<N; m += 2)
    {
        B_re00 = B_re11 = AE_MUL32_HH(sigma2, 1 << 30);
        pZ0 = (ae_int64x2*)&Z[(m*N + m)];
        pZ1 = (ae_int64x2*)&Z[(m*N + m)];
        for (n = m; n < N; n += 2)
        {
            B_re01 = B_re10 = 0;
            ep_re00 = AE_SEXT72(B_re00);
            ep_re11 = AE_SEXT72(B_re11);
            ep_re01 = ep_re10 = AE_MOVEA(0);
            pAkm = (ae_int32 *)&A[m];
            pAkn = (ae_int32 *)&A[n];
            for (k = 0; k < N; k++)
            {
                akm1 = AE_L32_I(pAkm, sizeof(ae_int32));
                AE_L32_IP(akm0, pAkm, N*sizeof(ae_int32));
                akn1 = AE_L32_I(pAkn, sizeof(ae_int32));
                AE_L32_IP(akn0, pAkn, N*sizeof(ae_int32));

                AE_MULA32EP_HH(ep_re00, B_re00, akm0, akn0);
                AE_MULA32EP_HH(ep_re01, B_re01, akm0, akn1);
                AE_MULA32EP_HH(ep_re10, B_re10, akm1, akn0);
                AE_MULA32EP_HH(ep_re11, B_re11, akm1, akn1);
            }
            B_re00 = AE_SLAA64S(AE_SRAI72(ep_re00, B_re00, 4), 2 * qRA);
            B_re01 = AE_SLAA64S(AE_SRAI72(ep_re01, B_re01, 4), 2 * qRA);
            B_re10 = AE_SLAA64S(AE_SRAI72(ep_re10, B_re10, 4), 2 * qRA);
            B_re11 = AE_SLAA64S(AE_SRAI72(ep_re11, B_re11, 4), 2 * qRA);

			AE_S64X2_X(B_re10, B_re11, pZ0, N*sizeof(ae_int64));
			AE_S64X2_IP(B_re00, B_re01, pZ0, 2 * sizeof(ae_int64));
			AE_S64X2_X(B_re01, B_re11, pZ1, N*sizeof(ae_int64));
			AE_S64X2_XP(B_re00, B_re10, pZ1, 2 * N*sizeof(ae_int64));
			B_re00 = B_re11 = 0;
        }
    }
#undef N
}
/* scratch allocation function */
size_t  matcholpreprocess6x6_32x32_getScratchSize()   { return 0; }
