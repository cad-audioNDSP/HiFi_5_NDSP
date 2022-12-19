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
static void computeAB(int64_t* restrict Z, const int32_t* restrict A, const int32_t* restrict B, int qYBRA)
#if 0
{
    int72_t B_re,B_im;
    int n,m;

    for(n=0; n<N; n++)
    {
        B_re = B_im = sext72(0);
        for (m=0; m<M; m++) 
        {
            int32_t a_re, a_im, b_re, b_im;

            a_re=A[2*n+m*N*2+0];
            a_im=A[2*n+m*N*2+1];
            b_re=B[m*2+0];
            b_im=B[m*2+1];

            B_re = add72x64(B_re, (int64_t)a_re*b_re);
            B_re = add72x64(B_re, (int64_t)a_im*b_im);
            B_im = add72x64(B_im, (int64_t)a_re*b_im);
            B_im = add72x64(B_im,-(int64_t)a_im*b_re);
        }
        Z[2*n+0] = LL_shl_ll(sra72(B_re,4),qYBRA+4);
        Z[2*n+1] = LL_shl_ll(sra72(B_im,4),qYBRA+4);
    }
}
#elif 0
{
    const int M=10,N=10;
    ae_int64 B_re0,B_im0;
    ae_int64 B_re1,B_im1;
    ae_ep ep_re0,ep_im0;
    ae_ep ep_re1,ep_im1;
    ae_int32x2 areim0,areim1,breim,t;
    int n,m;
    const int32_t * restrict A0=A;
    for(n=0; n<N; n+=2)
    {
        A=A0;
        A0+=4;
        areim1=AE_L32X2_I((const ae_int32x2*)A,sizeof(ae_int32x2));
        AE_L32X2_XP(areim0,castxcc(ae_int32x2,A),N*sizeof(ae_int32x2));
        AE_L32X2_IP(breim,castxcc(ae_int32x2,B),sizeof(ae_int32x2));
        AE_MULZAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
        AE_MULZAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
        t=AE_SEL32_LH(breim,AE_NEG32S(breim));
        AE_MULZAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
        AE_MULZAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);
        for (m=1; m<M; m++) 
        {
            areim1=AE_L32X2_I((const ae_int32x2*)A,sizeof(ae_int32x2));
            AE_L32X2_XP(areim0,castxcc(ae_int32x2,A),N*sizeof(ae_int32x2));
            AE_L32X2_IP(breim,castxcc(ae_int32x2,B),sizeof(ae_int32x2));
            AE_MULAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
            AE_MULAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
            t=AE_SEL32_LH(breim,AE_NEG32S(breim));
            AE_MULAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
            AE_MULAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);
        }
        B-=2*M;
        AE_S64_IP(AE_SLAA64S(AE_SRAI72(ep_re0,B_re0,4),qYBRA+4),castxcc(ae_int64,Z),sizeof(ae_int64));
        AE_S64_IP(AE_SLAA64S(AE_SRAI72(ep_im0,B_im0,4),qYBRA+4),castxcc(ae_int64,Z),sizeof(ae_int64));
        AE_S64_IP(AE_SLAA64S(AE_SRAI72(ep_re1,B_re1,4),qYBRA+4),castxcc(ae_int64,Z),sizeof(ae_int64));
        AE_S64_IP(AE_SLAA64S(AE_SRAI72(ep_im1,B_im1,4),qYBRA+4),castxcc(ae_int64,Z),sizeof(ae_int64));
    }
}
#else
{
    const int M=10,N=10;
    ae_int64 B_re0,B_im0;
    ae_int64 B_re1,B_im1;
    ae_ep ep_re0,ep_im0;
    ae_ep ep_re1,ep_im1;
    ae_int32x2 areim0,areim1,breim,t;
    int n;
    const int32_t * restrict A0=A;
    for(n=0; n<N; n+=2)
    {
        A=A0;
        A0+=4;
        AE_L32X2X2_XP(areim0,areim1,castxcc(ae_int32x4,A),N*sizeof(ae_int32x2));
        AE_L32X2_IP(breim,castxcc(ae_int32x2,B),sizeof(ae_int32x2));
        AE_MULZAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
        AE_MULZAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
        t=AE_SEL32_LH(breim,AE_NEG32S(breim));
        AE_MULZAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
        AE_MULZAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);

        AE_L32X2X2_XP(areim0,areim1,castxcc(ae_int32x4,A),N*sizeof(ae_int32x2));
        AE_L32X2_IP(breim,castxcc(ae_int32x2,B),sizeof(ae_int32x2));
        AE_MULAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
        AE_MULAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
        t=AE_SEL32_LH(breim,AE_NEG32S(breim));
        AE_MULAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
        AE_MULAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);

        AE_L32X2X2_XP(areim0,areim1,castxcc(ae_int32x4,A),N*sizeof(ae_int32x2));
        AE_L32X2_IP(breim,castxcc(ae_int32x2,B),sizeof(ae_int32x2));
        AE_MULAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
        AE_MULAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
        t=AE_SEL32_LH(breim,AE_NEG32S(breim));
        AE_MULAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
        AE_MULAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);

        AE_L32X2X2_XP(areim0,areim1,castxcc(ae_int32x4,A),N*sizeof(ae_int32x2));
        AE_L32X2_IP(breim,castxcc(ae_int32x2,B),sizeof(ae_int32x2));
        AE_MULAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
        AE_MULAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
        t=AE_SEL32_LH(breim,AE_NEG32S(breim));
        AE_MULAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
        AE_MULAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);

        AE_L32X2X2_XP(areim0,areim1,castxcc(ae_int32x4,A),N*sizeof(ae_int32x2));
        AE_L32X2_IP(breim,castxcc(ae_int32x2,B),sizeof(ae_int32x2));
        AE_MULAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
        AE_MULAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
        t=AE_SEL32_LH(breim,AE_NEG32S(breim));
        AE_MULAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
        AE_MULAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);

        AE_L32X2X2_XP(areim0,areim1,castxcc(ae_int32x4,A),N*sizeof(ae_int32x2));
        AE_L32X2_IP(breim,castxcc(ae_int32x2,B),sizeof(ae_int32x2));
        AE_MULAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
        AE_MULAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
        t=AE_SEL32_LH(breim,AE_NEG32S(breim));
        AE_MULAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
        AE_MULAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);

        AE_L32X2X2_XP(areim0,areim1,castxcc(ae_int32x4,A),N*sizeof(ae_int32x2));
        AE_L32X2_IP(breim,castxcc(ae_int32x2,B),sizeof(ae_int32x2));
        AE_MULAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
        AE_MULAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
        t=AE_SEL32_LH(breim,AE_NEG32S(breim));
        AE_MULAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
        AE_MULAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);

        AE_L32X2X2_XP(areim0,areim1,castxcc(ae_int32x4,A),N*sizeof(ae_int32x2));
        AE_L32X2_IP(breim,castxcc(ae_int32x2,B),sizeof(ae_int32x2));
        AE_MULAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
        AE_MULAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
        t=AE_SEL32_LH(breim,AE_NEG32S(breim));
        AE_MULAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
        AE_MULAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);

        AE_L32X2X2_XP(areim0,areim1,castxcc(ae_int32x4,A),N*sizeof(ae_int32x2));
        AE_L32X2_IP(breim,castxcc(ae_int32x2,B),sizeof(ae_int32x2));
        AE_MULAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
        AE_MULAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
        t=AE_SEL32_LH(breim,AE_NEG32S(breim));
        AE_MULAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
        AE_MULAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);

        AE_L32X2X2_XP(areim0,areim1,castxcc(ae_int32x4,A),N*sizeof(ae_int32x2));
        AE_L32X2_XP(breim,castxcc(ae_int32x2,B),-(M-1)*(int)sizeof(ae_int32x2));
        AE_MULAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
        AE_MULAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
        t=AE_SEL32_LH(breim,AE_NEG32S(breim));
        AE_MULAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
        AE_MULAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);
		AE_S64X2_IP(AE_SLAA64S(AE_SRAI72(ep_re0, B_re0, 4), qYBRA + 4), AE_SLAA64S(AE_SRAI72(ep_im0, B_im0, 4), qYBRA + 4), castxcc(ae_int64x2, Z), 2 * sizeof(ae_int64));
		AE_S64X2_IP(AE_SLAA64S(AE_SRAI72(ep_re1, B_re1, 4), qYBRA + 4), AE_SLAA64S(AE_SRAI72(ep_im1, B_im1, 4), qYBRA + 4), castxcc(ae_int64x2, Z), 2 * sizeof(ae_int64));
	}
}
#endif

/*-------------------------------------------------------
    forward with P==1
-------------------------------------------------------*/
static void fwdnx1(
                  int32_t* restrict y,
            const int32_t* restrict R, 
            const int32_t* restrict D,
            const int64_t* restrict Z)
#if 0
{
    int n, m;
    const int32_t *pR;
    int64_t A_re, A_im;
    int32_t d;
    int d_exp;

    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(Z,HIFI_SIMD_WIDTH);

    for(n=0; n<N; n++)
    {
        int72_t B_re, B_im;
        pR=R+(n*(n+1));
        // calculate A(:,n)'*B-Rn'*Y, 1xP
        B_re = sext72(Z[2*n+0]);
        B_im = sext72(Z[2*n+1]);
        for (m=0; m<n; m++)   
        {
            int32_t r_re, r_im;
            int32_t y_re, y_im;
            r_re = pR[m*2+0];
            r_im = pR[m*2+1];
            y_re = y[m*2+0];
            y_im = y[m*2+1];
            B_re = add72x64(B_re,-(int64_t)y_re*r_re);
            B_re = add72x64(B_re,-(int64_t)y_im*r_im);
            B_im = add72x64(B_im,-(int64_t)y_im*r_re);
            B_im = add72x64(B_im, (int64_t)y_re*r_im);
        }
        d     = D[2*n+0];
        d_exp = D[2*n+1];
        A_re = sra72(B_re,4);
        A_im = sra72(B_im,4);
        A_re = mulQ63Q31(A_re,d);
        A_im = mulQ63Q31(A_im,d);
        y[m*2+0] = satQ31(LL_shl_ll(A_re,d_exp-31+4));
        y[m*2+1] = satQ31(LL_shl_ll(A_im,d_exp-31+4));
    }
}
#elif 0
{
    const int N=10;
    int n, m;
    const ae_int32x2 * restrict pR;
          ae_int32x2 * restrict pY;
    ae_int64 A_re, A_im;
    ae_int32x2 t,d,rreim,yreim;
    ae_int64 B_re, B_im;
    ae_ep ep_re,ep_im;
    int d_exp;

    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(Z,HIFI_SIMD_WIDTH);

    pY=(ae_int32x2 *)y;
    for(n=0; n<N; n++)
    {
        pR=(const ae_int32x2 *)(R+(n*(n+1)));
        // calculate A(:,n)'*B-Rn'*Y, 1xP
        AE_L64_IP(B_re,castxcc(ae_int64,Z),sizeof(int64_t));
        AE_L64_IP(B_im,castxcc(ae_int64,Z),sizeof(int64_t));
        ep_re=AE_SEXT72(B_re); ep_im=AE_SEXT72(B_im);
        for (m=0; m<n; m++)   
        {
            AE_L32X2_IP(rreim,pR,sizeof(ae_int32x2));
            AE_L32X2_IP(yreim,pY,sizeof(ae_int32x2));
            AE_MULSSD32EP_HH_LL(ep_re,B_re,yreim,rreim);
            t=AE_SEL32_LH(rreim,AE_NEG32S(rreim));
            AE_MULAAD32EP_HH_LL(ep_im,B_im,yreim,t);
        }
        d_exp = D[1];
        AE_L32_IP(d,castxcc(ae_int32,D),sizeof(ae_int32x2));
        A_re = AE_SRAI72(ep_re,B_re, 4);
        A_im = AE_SRAI72(ep_im,B_im, 4);
        AE_MUL32USEP_LH(ep_re,B_re,AE_MOVINT32X2_FROMINT64(A_re),d);
        AE_MUL32USEP_LH(ep_im,B_im,AE_MOVINT32X2_FROMINT64(A_im),d);
        B_re=AE_SRAI72(ep_re,B_re,31);
        B_im=AE_SRAI72(ep_im,B_im,31);
        AE_MULAF32S_HH(B_re,AE_MOVINT32X2_FROMINT64(A_re),d);
        AE_MULAF32S_HH(B_im,AE_MOVINT32X2_FROMINT64(A_im),d);
        AE_S32X2_XP(AE_TRUNCA32X2F64S(B_re,B_im,d_exp+5),pY,-n*(int)sizeof(ae_int32x2));
    }
}
#else
{
    const int N=10;
    int m,n;
    ae_ep ep0_re,ep0_im;
    ae_ep ep1_re,ep1_im;
    ae_int64 B0_re, B0_im;
    ae_int64 B1_re, B1_im;
    ae_int64 A_re,A_im;
    ae_int32x2 d0,d1;
    int d_exp;

    for(n=0; n<N; n+=2)
    {
        ae_int32x2 rreim0,rreim1,yreim,t;
        R+=2*n;
        AE_L64X2_IP(B0_re,B0_im,castxcc(ae_int64x2,Z),2*sizeof(int64_t));
        AE_L64X2_IP(B1_re,B1_im,castxcc(ae_int64x2,Z),2*sizeof(int64_t));
        ep0_re=AE_SEXT72(B0_re);
        ep0_im=AE_SEXT72(B0_im);
        ep1_re=AE_SEXT72(B1_re);
        ep1_im=AE_SEXT72(B1_im);
        for (m=0; m<n; m++)   
        {
            rreim1=AE_L32X2_X ((ae_int32x2*)R,(n+1)*sizeof(ae_int32x2));
            AE_L32X2_IP(rreim0,castxcc(ae_int32x2,R),sizeof(ae_int32x2));
            AE_L32X2_IP(yreim,castxcc(ae_int32x2,y),sizeof(ae_int32x2));
            AE_MULSSD32EP_HH_LL(ep0_re,B0_re,rreim0,yreim);
            AE_MULSSD32EP_HH_LL(ep1_re,B1_re,rreim1,yreim);
            t=AE_SEL32_LH(yreim,AE_NEG32S(yreim));
            AE_MULSSD32EP_HH_LL(ep0_im,B0_im,rreim0,t);
            AE_MULSSD32EP_HH_LL(ep1_im,B1_im,rreim1,t);
        }
        d_exp=D[1];
        AE_L32_IP(d0,castxcc(ae_int32,D),2*sizeof(int32_t));
        A_re = AE_SRAI72(ep0_re,B0_re, 4);
        A_im = AE_SRAI72(ep0_im,B0_im, 4);
        AE_MUL32USEP_LH(ep0_re,B0_re,AE_MOVINT32X2_FROMINT64(A_re),d0);
        AE_MUL32USEP_LH(ep0_im,B0_im,AE_MOVINT32X2_FROMINT64(A_im),d0);
        B0_re=AE_SRAI72(ep0_re,B0_re,31);
        B0_im=AE_SRAI72(ep0_im,B0_im,31);
        AE_MULAF32S_HH(B0_re,AE_MOVINT32X2_FROMINT64(A_re),d0);
        AE_MULAF32S_HH(B0_im,AE_MOVINT32X2_FROMINT64(A_im),d0);
        d0=AE_TRUNCA32X2F64S(B0_re,B0_im,d_exp+5);
        yreim=d0;
        rreim1=AE_L32X2_X ((ae_int32x2*)R,(n+1)*sizeof(ae_int32x2));
        AE_MULSSD32EP_HH_LL(ep1_re,B1_re,rreim1,yreim);
        t=AE_SEL32_LH(yreim,AE_NEG32S(yreim));
        AE_MULSSD32EP_HH_LL(ep1_im,B1_im,rreim1,t);
        d_exp=D[1];
        AE_L32_IP(d1,castxcc(ae_int32,D),2*sizeof(int32_t));
        A_re = AE_SRAI72(ep1_re,B1_re, 4);
        A_im = AE_SRAI72(ep1_im,B1_im, 4);
        AE_MUL32USEP_LH(ep1_re,B1_re,AE_MOVINT32X2_FROMINT64(A_re),d1);
        AE_MUL32USEP_LH(ep1_im,B1_im,AE_MOVINT32X2_FROMINT64(A_im),d1);
        B1_re=AE_SRAI72(ep1_re,B1_re,31);
        B1_im=AE_SRAI72(ep1_im,B1_im,31);
        AE_MULAF32S_HH(B1_re,AE_MOVINT32X2_FROMINT64(A_re),d1);
        AE_MULAF32S_HH(B1_im,AE_MOVINT32X2_FROMINT64(A_im),d1);
        d1=AE_TRUNCA32X2F64S(B1_re,B1_im,d_exp+5);
		AE_S32X2X2_XP(d0, d1, castxcc(ae_int32x4, y), -(n)*(int)sizeof(ae_int32x2));
        R+=2;
    }
} 
#endif

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
void cmatcholfwdsubst10x10_32x32(void * restrict pScr,
                         complex_fract32* restrict y,
                   const complex_fract32* restrict R,
                   const complex_fract32* restrict D,
                   const complex_fract32* restrict A,
                   const complex_fract32* restrict B,
                         int qYBRA)
{
    int64_t *Z=(int64_t *)pScr;

    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(B,HIFI_SIMD_WIDTH);
    // compute A'*B
    computeAB(Z, (int32_t*)A, (int32_t*)B, qYBRA);
    fwdnx1((int32_t*)y,(int32_t*)R,(int32_t*)D,Z);
}

/* scratch allocation functions */
size_t  cmatcholfwdsubst10x10_32x32_getScratchSize()  { return (size_t)(2* 10* 10*sizeof(int64_t));}
