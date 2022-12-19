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

/*-----------------------------------------------
backward recursion:
makes backward recursion from transposed matrix
yt[NxN]
Input:
yt[NxN]    - transposed right part
Rt         - transposed Cholesky decomposition
D[N]       - diagonal terms
Output:
x[NxN] 
-----------------------------------------------*/
void cmatcholbkwsubstnxn(      complex_fract32* restrict x ,
                   const complex_fract32* Rt,
                   const complex_fract32 * restrict D,
                   const complex_fract32 * restrict yt, 
                   int N, int qXYR)
#if 0
{
    int n, m, k;
    const int32_t *pRt;
    int32_t  r_re, r_im;
    int32_t  x_re, x_im;
    int64_t  A_re, A_im;
    int72_t B_re, B_im;
    int32_t d;
    int d_exp;

    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(Rt, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(yt, HIFI_SIMD_WIDTH);
    for (k = N - 1; k >= 0; k--)
    {
        for (n=0; n<N; n++)
        {
            pRt = (const int32_t*)(Rt + (((N - k + 1)*(N - k - 2))>>1));
            // calculate y(m,:)-R(m,:)*X, 1xP
            B_re = sext72(LL_shl_ll(yt[k+n*N].s.re,qXYR)); // qY->qX+qR
            B_im = sext72(LL_shl_ll(yt[k+n*N].s.im,qXYR));
            for (m = N - k - 2; m>=0; m--)
            {
                x_re = x[((k + 1 + m)*N+n) ].s.re;
                x_im = x[((k + 1 + m)*N+n) ].s.im;
                r_re = pRt[0];
                r_im = pRt[1];
                pRt-=2;
                B_re = add72x64(B_re, -(int64_t)x_re*r_re);
                B_re = add72x64(B_re,  (int64_t)x_im*r_im);
                B_im = add72x64(B_im, -(int64_t)x_re*r_im);
                B_im = add72x64(B_im, -(int64_t)x_im*r_re);
            }
            d     = D[k].s.re;
            d_exp = D[k].s.im;
            A_re = sra72(B_re,4);
            A_im = sra72(B_im,4);
            A_re = mulQ63Q31(A_re,d);
            A_im = mulQ63Q31(A_im,d);
            x[(k*N+n)].s.re = satQ31(LL_shl_ll(A_re,d_exp-31+4));
            x[(k*N+n)].s.im = satQ31(LL_shl_ll(A_im,d_exp-31+4));
        }
    }
}
#else
{
    int n, m, k;
    const ae_int32x2 * restrict pY;
    const ae_int32x2 * restrict pRt;
    const ae_int32x2 * restrict pRt0;
          ae_int32x4 * restrict pX;
    ae_int32x2 rreim,xreim0,xreim1,t,y0,y1;
    ae_int64  A_re0, A_im0;
    ae_int64  A_re1, A_im1;
    ae_int64 B_re0, B_im0;
    ae_int64 B_re1, B_im1;
    ae_ep ep_re0,ep_im0;
    ae_ep ep_re1,ep_im1;
    ae_int32x2 d;
    int d_exp;

    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(Rt, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(yt, HIFI_SIMD_WIDTH);
    D+=N-1;
    x+=(N-1)*N;
    for (k = N - 1; k >= 0; k--)
    {
        d_exp = D[0].s.im;
        AE_L32_XP(d,castxcc(ae_int32,D),-(int)sizeof(ae_int32x2));
        pRt0 = (const ae_int32x2*)(Rt + (((N - k + 1)*(N - k - 2))>>1));
        pY=(const ae_int32x2*)&yt[k];
        for (n=0; n<N; n+=2)
        {
            pRt = pRt0;
            pX=(ae_int32x4*)&x[n];
            // calculate y(m,:)-R(m,:)*X, 1xP
            AE_L32X2_XP(y0,pY,N*sizeof(ae_int32x2));
            AE_L32X2_XP(y1,pY,N*sizeof(ae_int32x2));
            B_re0=AE_SLAA64S(AE_MUL32_HH(y0,1),qXYR);
            B_im0=AE_SLAA64S(AE_MUL32_LL(y0,1),qXYR);
            B_re1=AE_SLAA64S(AE_MUL32_HH(y1,1),qXYR);
            B_im1=AE_SLAA64S(AE_MUL32_LL(y1,1),qXYR);
            ep_re0=AE_SEXT72(B_re0); ep_im0=AE_SEXT72(B_im0);
            ep_re1=AE_SEXT72(B_re1); ep_im1=AE_SEXT72(B_im1);
            for (m = N - k - 2; m>=0; m--)
            {
                AE_L32X2X2_XP(xreim0,xreim1,pX ,-N*(int)sizeof(ae_int32x2));
                AE_L32X2_XP(rreim,pRt,-1*(int)sizeof(ae_int32x2));
                t=rreim; AE_MOVT32X2(t,AE_NEG32S(t),AE_MOVBA2(2));
                AE_MULAAD32EP_HH_LL(ep_re0,B_re0,xreim0,t);
                AE_MULAAD32EP_HH_LL(ep_re1,B_re1,xreim1,t);
                AE_MULSSD32EP_HH_LL(ep_im0,B_im0,xreim0,AE_SEL32_LH(rreim,rreim));
                AE_MULSSD32EP_HH_LL(ep_im1,B_im1,xreim1,AE_SEL32_LH(rreim,rreim));
            }
            A_re0 = AE_SRAI72(ep_re0,B_re0, 4);
            A_re1 = AE_SRAI72(ep_re1,B_re1, 4);
            A_im0 = AE_SRAI72(ep_im0,B_im0, 4);
            A_im1 = AE_SRAI72(ep_im1,B_im1, 4);
            AE_MUL32USEP_LH(ep_re0,B_re0,AE_MOVINT32X2_FROMINT64(A_re0),d);
            AE_MUL32USEP_LH(ep_im1,B_im1,AE_MOVINT32X2_FROMINT64(A_im1),d);
            AE_MUL32USEP_LH(ep_re0,B_re0,AE_MOVINT32X2_FROMINT64(A_re0),d);
            AE_MUL32USEP_LH(ep_re1,B_re1,AE_MOVINT32X2_FROMINT64(A_re1),d);
            B_re0=AE_SRAI72(ep_re0,B_re0,31);
            B_re1=AE_SRAI72(ep_re1,B_re1,31);
            B_im0=AE_SRAI72(ep_im0,B_im0,31);
            B_im1=AE_SRAI72(ep_im1,B_im1,31);
            AE_MULAF32S_HH(B_re0,AE_MOVINT32X2_FROMINT64(A_re0),d);
            AE_MULAF32S_HH(B_re1,AE_MOVINT32X2_FROMINT64(A_re1),d);
            AE_MULAF32S_HH(B_im0,AE_MOVINT32X2_FROMINT64(A_im0),d);
            AE_MULAF32S_HH(B_im1,AE_MOVINT32X2_FROMINT64(A_im1),d);
            AE_S32X2X2_XP(AE_TRUNCA32X2F64S(B_re0,B_im0,d_exp+5),AE_TRUNCA32X2F64S(B_re1,B_im1,d_exp+5),pX,0);
        }
        __Pragma("no_reorder")
    }
}
#endif
