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
#include "chol32x32_common.h"

/*--------------------------------------------------------
   make forward recursion to update n new column elements
   Input:
   Z[L][N+1][C]  convolutions in N-th column, Q(2*qR-4)
   D[L][SD]      reciprocals of main diagonal
   C             1 for real, 2 for complex
   Output:
   y             pointer to the N-th column in R[L][SR] 
                 (only N elements filled), Q(qR)
--------------------------------------------------------*/
void cmatcholfwdrec(int32_t *restrict y, const int32_t * restrict R, const int32_t * restrict D, const int64_t *Z, int N)
#if 0
{
    int m,n;
    int72_t B_re, B_im;
    int64_t A_re,A_im;
    int32_t d;
    int d_exp;

    for(n=0; n<N; n++)
    {
        R+=2*n;
        B_re = sla72(Z[2*n+0],4); // Q(2*qR-4)
        B_im = sla72(Z[2*n+1],4);
        for (m=0; m<n; m++)   
        {
            int32_t r_re,r_im;
            int32_t y_re,y_im;
            r_re=R[2*m+0]; r_im=R[2*m+1];
            y_re=y[2*m+0]; y_im=y[2*m+1];

            B_re = add72x64(B_re,-(int64_t)y_re*r_re);   
            B_re = add72x64(B_re,-(int64_t)y_im*r_im);
            B_im = add72x64(B_im,-(int64_t)y_im*r_re);
            B_im = add72x64(B_im, (int64_t)y_re*r_im);
        }
        d = D[2*n+0];
        d_exp = D[2*n+1];
        A_re = sra72(B_re, 4);
        A_im = sra72(B_im, 4);
        A_re = mulQ63Q31(A_re, d);
        A_im = mulQ63Q31(A_im, d);
        A_re = LL_shl_ll(A_re, d_exp-27);
        A_im = LL_shl_ll(A_im, d_exp-27);
        y[2*m+0] = satQ31(A_re);
        y[2*m+1] = satQ31(A_im);
    }
} /* fwdrec() */
#else
{
    int m,n;
    ae_ep ep0_re,ep0_im;
    ae_ep ep1_re,ep1_im;
    ae_int64 B0_re, B0_im;
    ae_int64 B1_re, B1_im;
    ae_int64 A_re,A_im;
    ae_int32x2 d;
    int d_exp;

    if (N&1)
    {
        AE_L64X2_IP(A_re,A_im,castxcc(ae_int64x2,Z),2*sizeof(int64_t));
        d_exp=D[1];
        AE_L32_IP(d,castxcc(ae_int32,D),2*sizeof(int32_t));
        AE_MUL32USEP_LH(ep0_re,B0_re,AE_MOVINT32X2_FROMINT64(A_re),d);
        AE_MUL32USEP_LH(ep0_im,B0_im,AE_MOVINT32X2_FROMINT64(A_im),d);
        B0_re=AE_SRAI72(ep0_re,B0_re,31);
        B0_im=AE_SRAI72(ep0_im,B0_im,31);
        AE_MULAF32S_HH(B0_re,AE_MOVINT32X2_FROMINT64(A_re),d);
        AE_MULAF32S_HH(B0_im,AE_MOVINT32X2_FROMINT64(A_im),d);
        d=AE_TRUNCA32X2F64S(B0_re,B0_im,d_exp+5);
        AE_S32X2_I(d,(ae_int32x2*)y,0);
    }
    for(n=(N&1); n<N; n+=2)
    {
        ae_int32x2 rreim0,rreim1,yreim,t;
        R+=2*n;
        AE_L64X2_IP(B0_re,B0_im,castxcc(ae_int64x2,Z),2*sizeof(int64_t));
        AE_L64X2_IP(B1_re,B1_im,castxcc(ae_int64x2,Z),2*sizeof(int64_t));
        AE_SLAI72(ep0_re,B0_re,B0_re,4);// Q(2*qR-4)
        AE_SLAI72(ep0_im,B0_im,B0_im,4);
        AE_SLAI72(ep1_re,B1_re,B1_re,4);
        AE_SLAI72(ep1_im,B1_im,B1_im,4);

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
        AE_L32_IP(d,castxcc(ae_int32,D),2*sizeof(int32_t));
        A_re = AE_SRAI72(ep0_re,B0_re, 4);
        A_im = AE_SRAI72(ep0_im,B0_im, 4);
        AE_MUL32USEP_LH(ep0_re,B0_re,AE_MOVINT32X2_FROMINT64(A_re),d);
        AE_MUL32USEP_LH(ep0_im,B0_im,AE_MOVINT32X2_FROMINT64(A_im),d);
        B0_re=AE_SRAI72(ep0_re,B0_re,31);
        B0_im=AE_SRAI72(ep0_im,B0_im,31);
        AE_MULAF32S_HH(B0_re,AE_MOVINT32X2_FROMINT64(A_re),d);
        AE_MULAF32S_HH(B0_im,AE_MOVINT32X2_FROMINT64(A_im),d);
        d=AE_TRUNCA32X2F64S(B0_re,B0_im,d_exp+5);
        AE_S32X2_IP(d,castxcc(ae_int32x2,y),sizeof(ae_int32x2));
        yreim=d;
        rreim1=AE_L32X2_X ((ae_int32x2*)R,(n+1)*sizeof(ae_int32x2));
        AE_MULSSD32EP_HH_LL(ep1_re,B1_re,rreim1,yreim);
        t=AE_SEL32_LH(yreim,AE_NEG32S(yreim));
        AE_MULSSD32EP_HH_LL(ep1_im,B1_im,rreim1,t);
        d_exp=D[1];
        AE_L32_IP(d,castxcc(ae_int32,D),2*sizeof(int32_t));
        A_re = AE_SRAI72(ep1_re,B1_re, 4);
        A_im = AE_SRAI72(ep1_im,B1_im, 4);
        AE_MUL32USEP_LH(ep1_re,B1_re,AE_MOVINT32X2_FROMINT64(A_re),d);
        AE_MUL32USEP_LH(ep1_im,B1_im,AE_MOVINT32X2_FROMINT64(A_im),d);
        B1_re=AE_SRAI72(ep1_re,B1_re,31);
        B1_im=AE_SRAI72(ep1_im,B1_im,31);
        AE_MULAF32S_HH(B1_re,AE_MOVINT32X2_FROMINT64(A_re),d);
        AE_MULAF32S_HH(B1_im,AE_MOVINT32X2_FROMINT64(A_im),d);
        d=AE_TRUNCA32X2F64S(B1_re,B1_im,d_exp+5);
        AE_S32X2_XP(d,castxcc(ae_int32x2,y),-(n+1)*(int)sizeof(ae_int32x2));

        R+=2;
    }
} 
#endif
