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
void cmatcholfwdsubst4x4_32x32(void * restrict pScr,
                         complex_fract32* restrict y,
                   const complex_fract32* restrict R,
                   const complex_fract32* restrict D,
                   const complex_fract32* restrict A,
                   const complex_fract32* restrict B,
                         int qYBRA)
{
    const complex_fract32 * restrict A0=A;
    #define N 4
    #define M 4
    ae_int64x2 *Z=(ae_int64x2 *)pScr;

    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(A,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(B,HIFI_SIMD_WIDTH);

    // compute A'*B
    {
        ae_int64 B_re0,B_im0;
        ae_int64 B_re1,B_im1;
        ae_ep ep_re0,ep_im0;
        ae_ep ep_re1,ep_im1;
        ae_int32x2 areim0,areim1,breim,t;
        int n;
        for(n=0; n<N; n+=2)
        {
            A=A0;
            A0+=2;
            //areim1=AE_L32X2_I((const ae_int32x2*)A,sizeof(ae_int32x2));
            AE_L32X2X2_IP(areim0,areim1,castxcc(ae_int32x4,A),N*sizeof(ae_int32x2));
            AE_L32X2_IP(breim,castxcc(ae_int32x2,B),sizeof(ae_int32x2));
            AE_MULZAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
            AE_MULZAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
            t=AE_SEL32_LH(breim,AE_NEG32S(breim));
            AE_MULZAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
            AE_MULZAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);

            AE_L32X2X2_IP(areim0,areim1,castxcc(ae_int32x4,A),N*sizeof(ae_int32x2));
            AE_L32X2_IP(breim,castxcc(ae_int32x2,B),sizeof(ae_int32x2));
            AE_MULAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
            AE_MULAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
            t=AE_SEL32_LH(breim,AE_NEG32S(breim));
            AE_MULAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
            AE_MULAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);

            AE_L32X2X2_IP(areim0,areim1,castxcc(ae_int32x4,A),N*sizeof(ae_int32x2));
            AE_L32X2_IP(breim,castxcc(ae_int32x2,B),sizeof(ae_int32x2));
            AE_MULAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
            AE_MULAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
            t=AE_SEL32_LH(breim,AE_NEG32S(breim));
            AE_MULAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
            AE_MULAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);

            AE_L32X2X2_IP(areim0,areim1,castxcc(ae_int32x4,A),N*sizeof(ae_int32x2));
            AE_L32X2_XP(breim,castxcc(ae_int32x2,B),-(M-1)*(int)sizeof(ae_int32x2));
            AE_MULAAD32EP_HH_LL(ep_re0,B_re0,areim0,breim);
            AE_MULAAD32EP_HH_LL(ep_re1,B_re1,areim1,breim);
            t=AE_SEL32_LH(breim,AE_NEG32S(breim));
            AE_MULAAD32EP_HH_LL(ep_im0,B_im0,areim0,t);
            AE_MULAAD32EP_HH_LL(ep_im1,B_im1,areim1,t);

			AE_S64X2_IP(AE_SLAA64S(AE_SRAI72(ep_re0, B_re0, 4), qYBRA + 4), AE_SLAA64S(AE_SRAI72(ep_im0, B_im0, 4), qYBRA + 4), castxcc(ae_int64x2, Z), 2*sizeof(ae_int64));
			AE_S64X2_IP(AE_SLAA64S(AE_SRAI72(ep_re1, B_re1, 4), qYBRA + 4), AE_SLAA64S(AE_SRAI72(ep_im1, B_im1, 4), qYBRA + 4), castxcc(ae_int64x2, Z), 2*sizeof(ae_int64));
        }
    }
    Z=(ae_int64x2*)pScr;
    // forward recursion
    {
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
            R+=n;
            AE_L64X2_IP(B0_re,B0_im,Z,2*sizeof(int64_t));
            AE_L64X2_IP(B1_re,B1_im,Z,2*sizeof(int64_t));
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
            d_exp=D[0].s.im;
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
            d_exp=D[0].s.im;
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
            R+=1;
        }
    } 
    #undef N
    #undef M

}
/* scratch allocation functions */
size_t  cmatcholfwdsubst4x4_32x32_getScratchSize()  { return (size_t)(2* 4* 4*sizeof(int64_t));}
