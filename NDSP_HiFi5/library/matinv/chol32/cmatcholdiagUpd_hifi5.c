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
   update N-th diagonal element
   Input:
   Z[L][N+1][C]  convolutions in N-th column
   C             1 for real, 2 for complex
   Input/output:
   y             pointer to the begining of column in matrix 
                 R[L][SR] (N+1 elements is written)
   Output:
   D[L][SD]      reciprocals of main diagonal (pointer to 
                 the N-th element
--------------------------------------------------------*/
void cmatcholdiagUpd(int32_t* y, int32_t* D, const int64_t* Z, int N)
#if 0
{
    int64_t A_re;
    int32_t d;
    int d_exp;
    int m;
    int72_t a;

    Z+=2*N; // points to diagonal element
    a=sla72(Z[0],4);
    for (m=0; m<N; m++)   
    {
        int32_t r_re, r_im;
        r_re=y[2*m+0]; r_im=y[2*m+1];
        a = add72x64(a,-(int64_t)r_re*r_re);
        a = add72x64(a,-(int64_t)r_im*r_im);
    }
    A_re=sra72(a,4);
    invSqrt_hifi(D, A_re);
    d     = D[0];
    d_exp = D[1];
    A_re = mulQ63Q31(A_re,d);
    A_re = LL_shl_ll(A_re,d_exp-27);
    y[2*N+0] = satQ31(A_re);
    y[2*N+1] = 0;
} 
#else
{
    ae_int64 a,A_re,t;
    ae_int32x2 d;
    int d_exp;
    int m;
    ae_ep ep;

    t=AE_L64_X((const ae_int64*)Z,2*N*sizeof(int64_t));
    AE_SLAI72(ep,a,t,4);
    for (m=0; m<N; m++)   
    {
        ae_int32x2 r;
        AE_L32X2_IP(r,castxcc(ae_int32x2,y),sizeof(ae_int32x2));
        AE_MULSSD32EP_HH_LL(ep ,a, r,r);
    }
    A_re=AE_SRAI72(ep,a,4);
    invSqrt_hifi(D, (int64_t)A_re);
    d     = D[0];
    d_exp = D[1];
    AE_MUL32USEP_LH(ep,a,AE_MOVINT32X2_FROMINT64(A_re),d);
    t=AE_SRAI72(ep,a,31);
    AE_MULAF32S_HH(t,AE_MOVINT32X2_FROMINT64(A_re),d);
    A_re=t;
    d=AE_TRUNCA32X2F64S(A_re,0,d_exp+5);
    AE_S32X2_I(d,(ae_int32x2*)y,0);
} 
#endif
/* diagUpd() */
