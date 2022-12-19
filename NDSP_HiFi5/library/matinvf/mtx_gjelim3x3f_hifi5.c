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
/*
   Real Matrix Gauss-Jordan Elimination for linear equation problem, 
   floating point API
   code optimized for HiFi5
 */

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Matrix functions */
#include "NatureDSP_Signal_matinv.h"
/* Common helper macros. */
#include "common.h"
#include "common_fpu.h"


#if (HAVE_VFPU)

/*-------------------------------------------------------------------------
  These functions implement Gauss elimination elimination process with full 
  pivoting to find solution of linear equations A*y=x
  
  Fixed point version takes representation of input matrix/vector and forms 
  representation of output vector with proper scaling.

  Precision: 
  f     floating point
  32x32 32-bit input, 32-bit output

  Input:
  A[N*N]      input matrix, representation is defined by parameter qA
  x[N]        input rigth side of equation. For fixed point API, 
              representation is defined by parameter qX
  qA          input matrix A representation (for fixed point API only)
  qX          input vector x representation (for fixed point API only)
  Output:
  y[N]        output vector
  Temporary:
  pScr        scratch memory. Size in bytes is defined by corresponding 
              scratch allocation function 
  N is 2,3,4,6,8,10

  Returned value: fixed point functions return fixed-point representation 
                  of resulted vector
  Restrictions:
  none
-------------------------------------------------------------------------*/
void  mtx_gjelim3x3f  (void* pScr, float32_t *y, const float32_t* A,const float32_t * x) 
{
    xtfloatx4 * restrict pk ; 
    xtfloatx4 * restrict pA  ;
    xtfloatx4 * restrict pA0 ;
    xtfloatx4 * restrict pA1 ;
    xtfloatx4 * restrict pmax;
    const int N=3;
    int k,n;
    float32_t *B; // [N+1][N]

    // allocate on scratch
    B=(float32_t *)pScr;
    /* Set bounds of the buffer */
    NASSERT_ALIGN16(B);
    WUR_AE_CBEGIN0((uintptr_t)(B));
    WUR_AE_CEND0  ((uintptr_t)(B+N*(N+1)));

    // copy inputs
    {
        ae_valignx2 aA;
        ae_valign   aX;
        xtfloatx2 a01,a23,a45,a67,a8,x01,x2;
        aA=AE_LA128_PP(A);
        AE_LASX2X2_IP(a01,a23,aA,castxcc(xtfloatx4,A));
        AE_LASX2X2_IP(a45,a67,aA,castxcc(xtfloatx4,A));
        a8=XT_LSI((const xtfloat*)A,0);
        aX=AE_LA64_PP(x);
        AE_LASX2IP(x01,aX,castxcc(xtfloatx2,x));
        x2=XT_LSI((const xtfloat*)x,0);
        pA  = (xtfloatx4 *)(B);
        AE_SSX2X2_IP(a01                     ,AE_SEL32_HH_SX2(a23,x01),pA,sizeof(xtfloatx4));
        AE_SSX2X2_IP(AE_SEL32_LH_SX2(a23,a45),AE_SEL32_LL_SX2(a45,x01),pA,sizeof(xtfloatx4));
        AE_SSX2X2_IP(a67                     ,AE_SEL32_HH_SX2(a8 ,x2 ),pA,sizeof(xtfloatx4));
    }

    pk  = (xtfloatx4 *)(B);
    pA0 = (xtfloatx4 *)(B+4);
    pA1 = (xtfloatx4 *)(B+8);
    for (k=0; k<N-1; k++)
    {
        xtfloatx2 Ak0,Ak1,Amax0,Amax1;        
        xtfloatx2 C0,C1,A00,A01,A10,A11;
        xtfloatx2 amax,R;
        int imax;

        /* pivoting */
        imax=k;
        amax=XT_CONST_S(0);
        /* find absolute max value in the k-th column */
        pA=(xtfloatx4*)&B[(N+1)*k+k];
        __Pragma("loop_count min=1,max=3")
        __Pragma("no_unroll")
        for(n=k; n<3; n++)
        {
            xtbool2 cond;
            xtfloatx2 t;
            {
                ae_int32x2 tt;
                AE_L32_XP(tt,castxcc(ae_int32,pA),(N+1)*sizeof(xtfloat));
                t=XT_AE_MOVXTFLOATX2_FROMINT32X2(tt);
            }
            BMAXNUMABS_SX2(cond,amax,amax,t);
            XT_MOVF_SX2(R  ,  t,cond);
            XT_MOVEQZ(imax, n, AE_MOVAB2(cond));
        }
        R = XT_RECIP_SX2(R);

        /* swap k-th and imax-th rows and scale them */
        pmax = (xtfloatx4 *)(B+imax*(N+1));
        AE_LSX2X2_I (Ak0,Ak1    ,pmax, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (Amax0,Amax1,pk  , 0*sizeof(xtfloatx4));
        AE_SSX2X2_IP(Amax0,Amax1,pmax, sizeof(xtfloatx4));
        __Pragma("no_reorder")
        /* multiply k-th row by the reciprocal pivot element during swapping rows  */
        MUL_SX2X2(Ak0,Ak1,Ak0,Ak1,R,R);
        AE_SSX2X2_IP(Ak0,Ak1,pk, sizeof(xtfloatx4));

        /* elimination */
        /* join forward and back substitution */
        C0 = XT_LSX((xtfloat *)pA0, k*sizeof(float32_t));
        C1 = XT_LSX((xtfloat *)pA1, k*sizeof(float32_t));
        AE_LSX2X2_I(A00,A01,(const xtfloatx4*)pA0,0*sizeof(xtfloatx4));
        AE_LSX2X2_I(A10,A11,(const xtfloatx4*)pA1,0*sizeof(xtfloatx4));
        MSUB_SX2X2(A00,A01, Ak0,Ak1, C0,C0);
        MSUB_SX2X2(A10,A11, Ak0,Ak1, C1,C1);
        AE_SSX2X2_XC(A00,A01,castxcc(xtfloatx4,pA0),sizeof(xtfloatx4));
        AE_SSX2X2_XC(A10,A11,castxcc(xtfloatx4,pA1),sizeof(xtfloatx4));
        __Pragma("no_reorder")
    }
    // last unrolled iteration with saving on the output
    {
        xtfloatx2 Ak0,Ak1;        
        xtfloatx2 C0,C1,A00,A01,A10,A11;
        xtfloatx2 R;

        {
            ae_int32x2 tt;
            tt=AE_L32_I((const ae_int32*)pA,0);
            tt=AE_L32_X((const ae_int32*)pk,k*sizeof(float32_t));
            R=XT_AE_MOVXTFLOATX2_FROMINT32X2(tt);
        }
        R = XT_RECIP_SX2(R);
        /* swap k-th and imax-th rows and scale them */
        AE_LSX2X2_I (Ak0,Ak1 ,pk, 0*sizeof(xtfloatx4));
        /* multiply k-th row by the reciprocal pivot element during swapping rows  */
        MUL_SX2X2(Ak0,Ak1,Ak0,Ak1,R,R);
        XT_SSI(Ak1,(xtfloat*)y,2*sizeof(float32_t));

        /* elimination */
        /* join forward and back substitution */
        C0 = XT_LSX((xtfloat *)pA0, k*sizeof(float32_t));
        C1 = XT_LSX((xtfloat *)pA1, k*sizeof(float32_t));
        AE_LSX2X2_I(A00,A01,(const xtfloatx4*)pA0,0*sizeof(xtfloatx4));
        AE_LSX2X2_I(A10,A11,(const xtfloatx4*)pA1,0*sizeof(xtfloatx4));
        MSUB_SX2X2(A01,A11, Ak1,Ak1, C0,C1);
        XT_SSI(A01,(xtfloat*)y,0*sizeof(float32_t));
        XT_SSI(A11,(xtfloat*)y,1*sizeof(float32_t));
    }
}

size_t mtx_gjelim3x3f_getScratchSize   () { const int N=3; return  N*(N+1)*sizeof(float32_t);  }
#elif (HAVE_FPU)
// code for scalar FPU

void  mtx_gjelim3x3f  (void* pScr, float32_t * restrict y, const float32_t* restrict A,const float32_t * restrict x) 
{
    const int N=3;
    int k,n;
    float32_t *restrict B; // [N][N]
    float32_t *restrict C; // [N]
    float32_t *restrict T; // [N]

    // allocate on scratch
    B=(float32_t *)pScr;
    C=B+N*N;
    T=C+N;
    // copy input
    for (k=0; k<N*N; k++) B[k]=A[k];
    for (k=0; k<N  ; k++) C[k]=x[k];

    for (k=0; k<N; k++)
    {
        float32_t bmax;
        int i,imax;

        // pivoting
        imax=k; bmax=0;
        for (n=k; n<N; n++)
        {
            xtbool cond;
            float32_t b=XT_ABS_S(B[n*N+k]);
            cond=XT_OLE_S(bmax,b);
            XT_MOVT_S(bmax,  b, cond);
            XT_MOVT  (imax,  n, cond);
        }
        for (n=0; n<N; n++) 
        {
            float32_t t;
            t=B[imax*N+n]; B[imax*N+n]=B[k*N+n]; B[k*N+n]=t;
        }
        {
            float32_t t;
            t=C[imax    ]; C[imax    ]=C[k    ]; C[k]=t;
        }
        // find normalization factor
        {
            float32_t rden,tk;
            for (n=0; n<N; n++) T[n]=B[n*N+k];
            tk=T[k]; T[k]=XT_CONST_S(1);
            rden=XT_RECIP_S(tk); // reciprocal of tk
            for (n=0; n<N; n++)
            {
                T[n]=T[n]*rden;
            }
        }
        for (i=0; i<N; i++)
        {
            xtfloat t;
            if(i==k) continue; 
            for (n=0; n<N; n++)
            {
                t=B[i*N+n];
                XT_MSUB_S(t,B[k*N+n],T[i]);
                B[i*N+n]=t;
            }
            t=C[i];
            XT_MSUB_S(t,C[k],T[i]);
            C[i]=t;
        }
        for (n=0; n<N; n++)
        {
            B[k*N+n]=XT_MUL_S(B[k*N+n],T[k]);
        }
        C[k]=XT_MUL_S(C[k],T[k]);
    }
    // copy back to the output
    for (k=0; k<N; k++) y[k]=C[k];
}

size_t mtx_gjelim3x3f_getScratchSize   () { const int N=3; return  N*(N+2)*sizeof(float32_t);  }
#else
DISCARD_FUN(void, mtx_gjelim3x3f, (void* pScr, float32_t *y, const float32_t* A, const float32_t * x))
size_t mtx_gjelim3x3f_getScratchSize()
{
  return 0;
}
#endif
