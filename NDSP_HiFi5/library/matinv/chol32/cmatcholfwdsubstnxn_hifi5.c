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

/* ---------------------------------------------------------------------
make forward recursion for R'\y=A'
Input:
R,D       - Cholesky decomposition
A[NxN]    - original input maatrix
Output:
yt[NxN]     transposed output
---------------------------------------------------------------------*/
void cmatcholfwdsubstnxn(   complex_fract32 * restrict yt,
                const complex_fract32 * restrict R, 
                const complex_fract32 * restrict D, 
                const complex_fract32 * restrict A, int N, int qYBRA)
#if 0
{
    int n, m,k;
    const int32_t *pR;
    int64_t A_re, A_im;
    int32_t d;
    int d_exp;

    NASSERT_ALIGN(yt,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A,HIFI_SIMD_WIDTH);

    for(m=0; m<N; m++)
    for(n=0; n<N; n++)
    {
        int72_t B_re, B_im;
        pR=(const int32_t*)(R+((n*(n+1))>>1));
        // calculate A(:,n)'*B-Rn'*Y, 1xP
        B_re = sext72(LL_shl_ll(        A[(n+m*N)].s.re ,31+qYBRA));
        B_im = sext72(LL_shl_ll(L_neg_l(A[(n+m*N)].s.im),31+qYBRA));
        for (k=0; k<n; k++)   
        {
            int32_t r_re, r_im;
            int32_t y_re, y_im;
            r_re = pR[k*2+0];
            r_im = pR[k*2+1];
            y_re = yt[(k+m*N)].s.re;
            y_im = yt[(k+m*N)].s.im;
            B_re = add72x64(B_re,-(int64_t)y_re*r_re);
            B_re = add72x64(B_re,-(int64_t)y_im*r_im);
            B_im = add72x64(B_im,-(int64_t)y_im*r_re);
            B_im = add72x64(B_im, (int64_t)y_re*r_im);
        }
        d     = D[n].s.re;
        d_exp = D[n].s.im;
        A_re = sra72(B_re,4);
        A_im = sra72(B_im,4);
        A_re = mulQ63Q31(A_re,d);
        A_im = mulQ63Q31(A_im,d);
        yt[(k+m*N)].s.re = satQ31(LL_shl_ll(A_re,d_exp-31+4));
        yt[(k+m*N)].s.im = satQ31(LL_shl_ll(A_im,d_exp-31+4));
    }
}
#else
{
    int n, m,k;
    const ae_int32x2 * restrict pR;
          ae_int32x2 * restrict pY;
    ae_int64 A_re0, A_im0;
    ae_int64 A_re1, A_im1;
    ae_int32x2 d;
    int d_exp;

    NASSERT_ALIGN(yt,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A,HIFI_SIMD_WIDTH);

    pY=(      ae_int32x2*)yt;
    for(m=0; m<N; m+=2)
    {
        for(n=0; n<N; n++)
        {
            ae_int64 B_re0, B_im0;
            ae_int64 B_re1, B_im1;
            ae_ep ep_re0,ep_im0;
            ae_ep ep_re1,ep_im1;
            ae_int32x2 a0,a1;
            pR=(const ae_int32x2*)(R+((n*(n+1))>>1));
            // calculate A(:,n)'*B-Rn'*Y, 1xP
            a1=AE_L32X2_X((const ae_int32x2*)A,N*sizeof(ae_int32x2));
            AE_L32X2_IP(a0,castxcc(ae_int32x2,A),sizeof(ae_int32x2));
            B_re0=AE_SLAA64S(AE_MUL32_HH(a0,AE_MOVDA32X2(1,-1)),31+qYBRA);
            B_im0=AE_SLAA64S(AE_MUL32_LL(a0,AE_MOVDA32X2(1,-1)),31+qYBRA);
            B_re1=AE_SLAA64S(AE_MUL32_HH(a1,AE_MOVDA32X2(1,-1)),31+qYBRA);
            B_im1=AE_SLAA64S(AE_MUL32_LL(a1,AE_MOVDA32X2(1,-1)),31+qYBRA);
            ep_re0=AE_SEXT72(B_re0); ep_im0=AE_SEXT72(B_im0);
            ep_re1=AE_SEXT72(B_re1); ep_im1=AE_SEXT72(B_im1);
            for (k=0; k<n; k++)   
            {
                ae_int32x2 rreim,yreim0,yreim1,t;
                AE_L32X2_IP(rreim,pR,sizeof(ae_int32x2));
                yreim1=AE_L32X2_X (pY,N*sizeof(ae_int32x2));
                AE_L32X2_IP(yreim0,pY,sizeof(ae_int32x2));
                AE_MULSSD32EP_HH_LL(ep_re0,B_re0,yreim0,rreim);
                AE_MULSSD32EP_HH_LL(ep_re1,B_re1,yreim1,rreim);
                t=AE_SEL32_LH(rreim,AE_NEG32S(rreim));
                AE_MULAAD32EP_HH_LL(ep_im0,B_im0,yreim0,t);
                AE_MULAAD32EP_HH_LL(ep_im1,B_im1,yreim1,t);
            }
            d_exp = D[0].s.im;
            AE_L32_XP(d,castxcc(ae_int32,D),sizeof(ae_int32x2));
            A_re0 = AE_SRAI72(ep_re0,B_re0, 4);
            A_re1 = AE_SRAI72(ep_re1,B_re1, 4);
            A_im0 = AE_SRAI72(ep_im0,B_im0, 4);
            A_im1 = AE_SRAI72(ep_im1,B_im1, 4);
            AE_MUL32USEP_LH(ep_re0,B_re0,AE_MOVINT32X2_FROMINT64(A_re0),d);
            AE_MUL32USEP_LH(ep_re1,B_re1,AE_MOVINT32X2_FROMINT64(A_re1),d);
            AE_MUL32USEP_LH(ep_im0,B_im0,AE_MOVINT32X2_FROMINT64(A_im0),d);
            AE_MUL32USEP_LH(ep_im1,B_im1,AE_MOVINT32X2_FROMINT64(A_im1),d);
            B_re0=AE_SRAI72(ep_re0,B_re0,31);
            B_re1=AE_SRAI72(ep_re1,B_re1,31);
            B_im0=AE_SRAI72(ep_im0,B_im0,31);
            B_im1=AE_SRAI72(ep_im1,B_im1,31);
            AE_MULAF32S_HH(B_re0,AE_MOVINT32X2_FROMINT64(A_re0),d);
            AE_MULAF32S_HH(B_re1,AE_MOVINT32X2_FROMINT64(A_re1),d);
            AE_MULAF32S_HH(B_im0,AE_MOVINT32X2_FROMINT64(A_im0),d);
            AE_MULAF32S_HH(B_im1,AE_MOVINT32X2_FROMINT64(A_im1),d);
            AE_S32X2_X (AE_TRUNCA32X2F64S(B_re1,B_im1,d_exp+5),pY,N*sizeof(ae_int32x2));
            AE_S32X2_XP(AE_TRUNCA32X2F64S(B_re0,B_im0,d_exp+5),pY,-n*(int)sizeof(ae_int32x2));
        }
        D-=N;
        A+=N;
        pY+=2*N;
    }
}
#endif
