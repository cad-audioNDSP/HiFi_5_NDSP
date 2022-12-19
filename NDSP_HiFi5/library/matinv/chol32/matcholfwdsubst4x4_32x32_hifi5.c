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

/*---------------------------------------------------------------
    compute L of matrix product Z[Nx1]=A[MxN]'*B[Mx1]
    Input:
    A[SA][C]    matrix MxN
    B[SB][C]    matrix Mx1
    C           1 for real, 2 for complex
    Output:
    Z[N*1][C]   matrix Nx1
---------------------------------------------------------------*/
static void rcomputeAB(int64_t* Z, const int32_t* A, const int32_t* B, int qYBRA)
{
#define N 4
    ae_int64x2 * restrict pZ = (ae_int64x2 *)Z;
    const ae_int32 * restrict pA;
    const ae_int32 * restrict pB;
    ae_ep ep_re0, ep_re1;
    ae_int32x2 a_re0, a_re1, b_re0;
    ae_int64 B_re0, B_re1;

    int n;
    for(n=0; n<N; n+=2)
    {
        ep_re0 = ep_re1 = AE_MOVEA(0);
        B_re0 = B_re1 = 0;
        pB = (const ae_int32 *)B;
        pA = (const ae_int32 *)&A[n];
        //for (m=0; m<M; m++) 
        {
            a_re1 = AE_L32_I(pA, sizeof(ae_int32));
            AE_L32_IP(a_re0, pA, N * sizeof(ae_int32));
            AE_L32_IP(b_re0, pB, sizeof(ae_int32));
            AE_MULA32EP_HH(ep_re0, B_re0, a_re0, b_re0);
            AE_MULA32EP_HH(ep_re1, B_re1, a_re1, b_re0);

            a_re1 = AE_L32_I(pA, sizeof(ae_int32));
            AE_L32_IP(a_re0, pA, N * sizeof(ae_int32));
            AE_L32_IP(b_re0, pB, sizeof(ae_int32));
            AE_MULA32EP_HH(ep_re0, B_re0, a_re0, b_re0);
            AE_MULA32EP_HH(ep_re1, B_re1, a_re1, b_re0);

            a_re1 = AE_L32_I(pA, sizeof(ae_int32));
            AE_L32_IP(a_re0, pA, N * sizeof(ae_int32));
            AE_L32_IP(b_re0, pB, sizeof(ae_int32));
            AE_MULA32EP_HH(ep_re0, B_re0, a_re0, b_re0);
            AE_MULA32EP_HH(ep_re1, B_re1, a_re1, b_re0);

            a_re1 = AE_L32_I(pA, sizeof(ae_int32));
            AE_L32_IP(a_re0, pA, N * sizeof(ae_int32));
            AE_L32_IP(b_re0, pB, sizeof(ae_int32));
            AE_MULA32EP_HH(ep_re0, B_re0, a_re0, b_re0);
            AE_MULA32EP_HH(ep_re1, B_re1, a_re1, b_re0);
        }
        B_re0 = AE_SLAA64S(AE_SRAI72(ep_re0, B_re0, 4), qYBRA + 4);
        B_re1 = AE_SLAA64S(AE_SRAI72(ep_re1, B_re1, 4), qYBRA + 4);

        AE_S64X2_IP(B_re0, B_re1, pZ, 2*sizeof(ae_int64));
    }
#undef N
}

/*-------------------------------------------------------
    forward with P==1
-------------------------------------------------------*/
static void rfwdnx1(
                  int32_t* restrict y,
            const int32_t* restrict R, 
            const int32_t* restrict D,
            const int64_t* restrict Z)
{
#define N 4
    int m, n;
    const ae_int64x2 * restrict pZ = (const ae_int64x2 *)Z;
    const ae_int32 * restrict pR0;
    const ae_int32 * restrict pR1;
    ae_int32 * restrict pY = (ae_int32 *)y;
    ae_ep ep_re0, ep_re1;
    ae_int32x2 y0, r0, d0;
    ae_int32x2 y1, r1, d1;
    ae_int64 A_re0, A_re1;
    ae_int64 B_re0, B_re1;
    int d_exp0, d_exp1;

    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(Z, HIFI_SIMD_WIDTH);

    for (n = 0; n<N; n += 2)
    {
        ep_re0 = ep_re1 = AE_MOVEA(0);
        B_re0 = B_re1 = 0;
        R += n;
        AE_L64X2_IP(B_re0, B_re1, pZ, 2*sizeof(ae_int64));
        ep_re0 = AE_SEXT72(B_re0);
        ep_re1 = AE_SEXT72(B_re1);
        pY = (ae_int32 *)y;
        ae_valign va0, va1;
        pR0 = (const ae_int32 *)R;
        pR1 = (const ae_int32 *)R + (n + 1);
        va0 = AE_LA64_PP(pR0);
        va1 = AE_LA64_PP(pR1);
        for (m = 0; m<n; m+=2)
		{
            AE_LA32X2_IP(r0, va0, castxcc(ae_int32x2, pR0));
            AE_LA32X2_IP(r1, va1, castxcc(ae_int32x2, pR1));
            AE_L32X2_IP(y0, castxcc(ae_int32x2, pY), sizeof(ae_int32x2));
            AE_MULSSD32EP_HH_LL(ep_re0, B_re0, r0, y0);
            AE_MULSSD32EP_HH_LL(ep_re1, B_re1, r1, y0);
        }
        d_exp0 = D[1];
        AE_L32_IP(d0, castxcc(ae_int32, D), sizeof(ae_int32x2));
        A_re0 = AE_SRAI72(ep_re0, B_re0, 4);
        AE_MUL32USEP_LH(ep_re0, B_re0, AE_MOVINT32X2_FROMINT64(A_re0), d0);
        B_re0 = AE_SRAI72(ep_re0, B_re0, 31);
        AE_MULAF32S_HH(B_re0, AE_MOVINT32X2_FROMINT64(A_re0), d0);
        y0 = AE_TRUNCA32F64S(B_re0, d_exp0 + 5);
        AE_S32_L_IP(y0, pY, sizeof(ae_int32));

        AE_L32_IP(r1, pR1, 0);
        AE_MULS32EP_HH(ep_re1, B_re1, r1, y0);
        d_exp1 = D[1];
        AE_L32_IP(d1, castxcc(ae_int32, D), sizeof(ae_int32x2));
        A_re1 = AE_SRAI72(ep_re1, B_re1, 4);
        AE_MUL32USEP_LH(ep_re1, B_re1, AE_MOVINT32X2_FROMINT64(A_re1), d1);
        B_re1 = AE_SRAI72(ep_re1, B_re1, 31);
        AE_MULAF32S_HH(B_re1, AE_MOVINT32X2_FROMINT64(A_re1), d1);
        y1 = AE_TRUNCA32F64S(B_re1, d_exp1 + 5);
        AE_S32_L_IP(y1, pY, sizeof(ae_int32));
        R += 1;
    }
#undef N
}

/*-------------------------------------------------------------------------
Cholesky Forward Substitution for Pseudo-inversion
These functions make forward recursion stage of pseudo-inversion. They use
Cholesky decomposition R[NxN] of original matrices A[MxN]

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
R[((N+1)*N)/2]
            upper triangular matrix R. For fixed point, representation is Q(qR)
A[M*N]      matrix A. For fixed point, representation is Q(qA)
B[M*P]      original right-side matrix B. For fixed point, representation is Q(qB)
D[N]        reciprocals of main diagonal. NOTE: for the fixed point API,
            these data are stored internally in special format with separate
            mantissa and exponent for better accuracy and dynamic range 
            control. So, even for the real data, they stored as pairs of 2
            integers and packed to the complex_fract32 format
qYBRA       qY-qB+qR-qA, combination of fixed point representations of 
            matrices y,B,R and A (for the fixed point API only)
Output:
y[N*P]      Decision matrix y. For fixed point, representation is Q(qY)
Temporary:
pScr        Scratch memory

N = M = 4, 6, 8, 10
P = 1

Restrictions:
All matrices and scratch memory should not overlap and be aligned on
16-bytes boundary
---------------------------------------------------------------------------*/
void matcholfwdsubst4x4_32x32(void * restrict pScr,
                         int32_t* restrict y,
                   const int32_t* restrict R,
                   const complex_fract32* restrict D,
                   const int32_t* restrict A,
                   const int32_t* restrict B,
                         int qYBRA)
{
    int64_t *Z=(int64_t *)pScr;

    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(B,HIFI_SIMD_WIDTH);
    rcomputeAB(Z, A, B, qYBRA);
    rfwdnx1(y,R,(int32_t*)D,Z);
}


/* scratch allocation functions */
size_t  matcholfwdsubst4x4_32x32_getScratchSize()  { return (size_t)( 4* 4*sizeof(int64_t));}
