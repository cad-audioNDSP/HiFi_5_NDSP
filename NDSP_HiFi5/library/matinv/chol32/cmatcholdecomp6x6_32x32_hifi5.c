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

/*----------------------------------------------------------------------------------------
the same as cmatcholpreprocess10x10_32x32, but computes only lower triangle and
conjugate result
The result is matrix Z[NxN], such that
Z = (A'*A)' + sigma2*I[NxN], where ' denotes the conjugate transpose of
a matrix, and sigma2*I[NxN] is the NxN identity matrix multiplied with
the regularization term.

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

----------------------------------------------------------------------------------------*/
static void  preprocess_conj(
            complex_fract64 *Z,
    const   complex_fract32 * A,
            int32_t  sigma2, int qRA)
#if 0
{
    int m,k,n;
    const int N=6;
    for (m=0; m<N; m++)
    for (n=0; n<N; n++)
    {
        int72_t B_re,B_im;
        int32_t akm_re,akn_re, akm_im,akn_im;
        B_re = sla72((m==n)?sigma2:0,30);
        B_im.hi=0; B_im.lo = 0;
        for (k=0; k<N; k++)
        {
            akm_re=A[(k*N+m)].s.re; akm_im=A[(k*N+m)].s.im;
            akn_re=A[(k*N+n)].s.re; akn_im=A[(k*N+n)].s.im;
            B_re = add72x64(B_re, (int64_t)akm_re*akn_re); B_re = add72x64(B_re, (int64_t)akm_im*akn_im);
            B_im = add72x64(B_im,-(int64_t)akm_re*akn_im); B_im = add72x64(B_im, (int64_t)akm_im*akn_re);
        }
        Z[m*N+n].s.re = LL_shl_ll(sra72(B_re,4),2*qRA); /* -> Q(2*qA)->Q(2*qR-4) */
        Z[m*N+n].s.im = LL_shl_ll(sra72(B_im,4),2*qRA);
    }
}
#else
{
    ae_int64x2 * restrict pZ1;
    const ae_int32x4 * restrict pAm;
    const ae_int32x2 * restrict pAn;
    int m,n;
#define N 6
    for (m=0; m<N; m+=2)
    {
        int cond=2;
        pZ1=(ae_int64x2*)&Z[(m*N+m)];
        pAn=(const ae_int32x2 *)&A[m];
		pAm=(const ae_int32x4 *)&A[m];
        for (n=m; n<N; n++)
        {
            ae_int64 B_re0,B_im0,B_re1,B_im1;
            ae_ep ep_re0,ep_re1,ep_im0,ep_im1;
            ae_int32x2 akmreim0,akmreim1,aknreim,t;
            ae_int32x2 s=sigma2;
            AE_MOVF32X2(s,0,AE_MOVBA2(cond)); cond>>=1;
            B_re0=AE_MUL32_HH(s,1<<30);
            B_re1=AE_MUL32_LH(s,1<<30);
            ep_re0=AE_SEXT72(B_re0);
            ep_re1=AE_SEXT72(B_re1);

            AE_L32X2X2_IP(akmreim0,akmreim1,pAm,N*sizeof(ae_int32x2)); AE_L32X2_IP(aknreim ,pAn,N*sizeof(ae_int32x2));
            AE_MULAAD32EP_HH_LL(ep_re0,B_re0,akmreim0,aknreim);  AE_MULAAD32EP_HH_LL(ep_re1,B_re1,akmreim1,aknreim);
            t=AE_SEL32_LH(aknreim,AE_NEG32S(aknreim));
            AE_MULZAAD32EP_HH_LL(ep_im0,B_im0,akmreim0,t);        AE_MULZAAD32EP_HH_LL(ep_im1,B_im1,akmreim1,t);

            AE_L32X2X2_IP(akmreim0,akmreim1,pAm,N*sizeof(ae_int32x2)); AE_L32X2_IP(aknreim ,pAn,N*sizeof(ae_int32x2));
            AE_MULAAD32EP_HH_LL(ep_re0,B_re0,akmreim0,aknreim);  AE_MULAAD32EP_HH_LL(ep_re1,B_re1,akmreim1,aknreim);
            t=AE_SEL32_LH(aknreim,AE_NEG32S(aknreim));
            AE_MULAAD32EP_HH_LL(ep_im0,B_im0,akmreim0,t);        AE_MULAAD32EP_HH_LL(ep_im1,B_im1,akmreim1,t);

            AE_L32X2X2_IP(akmreim0,akmreim1,pAm,N*sizeof(ae_int32x2)); AE_L32X2_IP(aknreim ,pAn,N*sizeof(ae_int32x2));
            AE_MULAAD32EP_HH_LL(ep_re0,B_re0,akmreim0,aknreim);  AE_MULAAD32EP_HH_LL(ep_re1,B_re1,akmreim1,aknreim);
            t=AE_SEL32_LH(aknreim,AE_NEG32S(aknreim));
            AE_MULAAD32EP_HH_LL(ep_im0,B_im0,akmreim0,t);        AE_MULAAD32EP_HH_LL(ep_im1,B_im1,akmreim1,t);

            AE_L32X2X2_IP(akmreim0,akmreim1,pAm,N*sizeof(ae_int32x2)); AE_L32X2_IP(aknreim ,pAn,N*sizeof(ae_int32x2));
            AE_MULAAD32EP_HH_LL(ep_re0,B_re0,akmreim0,aknreim);  AE_MULAAD32EP_HH_LL(ep_re1,B_re1,akmreim1,aknreim);
            t=AE_SEL32_LH(aknreim,AE_NEG32S(aknreim));
            AE_MULAAD32EP_HH_LL(ep_im0,B_im0,akmreim0,t);        AE_MULAAD32EP_HH_LL(ep_im1,B_im1,akmreim1,t);

            AE_L32X2X2_IP(akmreim0,akmreim1,pAm,N*sizeof(ae_int32x2)); AE_L32X2_IP(aknreim ,pAn,N*sizeof(ae_int32x2));
            AE_MULAAD32EP_HH_LL(ep_re0,B_re0,akmreim0,aknreim);  AE_MULAAD32EP_HH_LL(ep_re1,B_re1,akmreim1,aknreim);
            t=AE_SEL32_LH(aknreim,AE_NEG32S(aknreim));
            AE_MULAAD32EP_HH_LL(ep_im0,B_im0,akmreim0,t);        AE_MULAAD32EP_HH_LL(ep_im1,B_im1,akmreim1,t);

            AE_L32X2X2_XP(akmreim0,akmreim1,pAm,-(N-1)*N*(int)sizeof(ae_int32x2)); AE_L32X2_XP(aknreim ,pAn,-(N-1)*N*(int)sizeof(ae_int32x2));
            AE_MULAAD32EP_HH_LL(ep_re0,B_re0,akmreim0,aknreim);  AE_MULAAD32EP_HH_LL(ep_re1,B_re1,akmreim1,aknreim);
            t=AE_SEL32_LH(aknreim,AE_NEG32S(aknreim));
            AE_MULAAD32EP_HH_LL(ep_im0,B_im0,akmreim0,t);        AE_MULAAD32EP_HH_LL(ep_im1,B_im1,akmreim1,t);
            pAn++;

            B_re0=AE_SLAA64S(AE_SRAI72(ep_re0,B_re0,4),2*qRA);
            B_re1=AE_SLAA64S(AE_SRAI72(ep_re1,B_re1,4),2*qRA);
            B_im0=AE_SLAA64S(AE_SRAI72(ep_im0,B_im0,4),2*qRA);
            B_im1=AE_SLAA64S(AE_SRAI72(ep_im1,B_im1,4),2*qRA);
            /* -> Q(2*qA)->Q(2*qR-4) */
			AE_S64X2_I(B_re1, B_im1, pZ1, 2 * sizeof(ae_int64));
			AE_S64X2_XP(B_re0, B_im0, pZ1, 2 * N*sizeof(ae_int64));
        }
    }
#undef N
}
#endif

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
void cmatcholdecomp6x6_32x32  (void * restrict pScr,
                      complex_fract32 * restrict R,
                      complex_fract32 * restrict D,
                const complex_fract32 * restrict A,
                      int32_t  sigma2,
                      int qRA)
{
    const int N=6 ;
    int64_t *Z = (int64_t *)pScr; // L columns of A'*A
    int32_t *y; // pointer to the new column in R
    int n;
    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A,HIFI_SIMD_WIDTH);
    preprocess_conj((complex_fract64*)Z,A,sigma2, qRA);

    /* to avoid ferret warnings */
    {

      ae_int32x4 *T0 = (ae_int32x4 *)R;
      ae_int32x4 *T1 = (ae_int32x4 *)D;

      AE_S32X2X2_IP(AE_ZERO32(), AE_ZERO32(), T1, sizeof(ae_int32x4));
      AE_S32X2X2_IP(AE_ZERO32(), AE_ZERO32(), T1, sizeof(ae_int32x4));
      AE_S32X2X2_IP(AE_ZERO32(), AE_ZERO32(), T1, sizeof(ae_int32x4));
    
      for (n = 0; n<((6 * 7) >> 2); n++)
      {
        AE_S32X2X2_IP(AE_ZERO32(), AE_ZERO32(), T0, sizeof(ae_int32x4));
      }
      AE_S32X2_IP(AE_ZERO32(), castxcc(ae_int32x2, T0), sizeof(ae_int32x2));
    }


    y=(int32_t *)R;
    for (n=0; n<N; n++)
    {
        y+=2*n; // go to the next column
        // make forward recursion to update n new column elements
        cmatcholfwdrec(y,(int32_t *)R,(int32_t *)D,Z,n);
        // update n-th diagonal element
        cmatcholdiagUpd(y,(int32_t *)(D+n),Z,n);
        Z+=2*N;
    }
}
size_t   cmatcholdecomp6x6_32x32_getScratchSize()   { return (2* 6*6*sizeof(int64_t));}
