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
void  mtx_gjelim4x4f  (void* pScr, float32_t *y, const float32_t* A,const float32_t * x) 
{
    xtfloatx2 *restrict pk;
    xtfloatx2 *restrict pmax;
    xtfloatx2 *restrict pA;
    xtfloatx2 *restrict pA0;
    xtfloatx2 *restrict pA1;
    xtfloatx2 *restrict pA2;
    xtfloatx2 R, C0, C1, C2;
    xtfloatx2 A00, A01, A02;
    xtfloatx2 A10, A11, A12;
    xtfloatx2 A20, A21, A22;
    xtfloatx2 Amax0, Amax1, Amax2,
                Ak0,   Ak1,   Ak2;
    const int N=4;
    int k,n;
    float32_t *B; // [N+4][N]

    // allocate on scratch
    B=(float32_t *)pScr;
    /* Set bounds of the buffer */
    NASSERT_ALIGN16(B);
    WUR_AE_CBEGIN0((uintptr_t)(B));
    WUR_AE_CEND0  ((uintptr_t)(B+N*(N+4)));

    // copy input
    {
        xtfloatx2 x0,x1,x2,x3,x4,x5,x6,x7,y0,y1;
        ae_valignx2 aA,aX;
        pA =(xtfloatx2*)B;
        aA=AE_LA128_PP(A);
        aX=AE_LA128_PP(x);
        AE_LASX2X2_IP(x0,x1,aA,castxcc(xtfloatx4,A));
        AE_LASX2X2_IP(x2,x3,aA,castxcc(xtfloatx4,A));
        AE_LASX2X2_IP(x4,x5,aA,castxcc(xtfloatx4,A));
        AE_LASX2X2_IP(x6,x7,aA,castxcc(xtfloatx4,A));
        AE_LASX2X2_IP(y0,y1,aX,castxcc(xtfloatx4,x));
        AE_SSX2X2_IP(x0,x1,castxcc(xtfloatx4,pA) ,sizeof(xtfloatx4));
        AE_SSX2X2_IP(y0,XT_CONST_S(0),castxcc(xtfloatx4,pA) ,sizeof(xtfloatx4));
        AE_SSX2X2_IP(x2,x3,castxcc(xtfloatx4,pA) ,sizeof(xtfloatx4));
        AE_SSX2X2_IP(AE_SEL32_LL_SX2(y0,y0),XT_CONST_S(0),castxcc(xtfloatx4,pA) ,sizeof(xtfloatx4));
        AE_SSX2X2_IP(x4,x5,castxcc(xtfloatx4,pA) ,sizeof(xtfloatx4));
        AE_SSX2X2_IP(y1,XT_CONST_S(0),castxcc(xtfloatx4,pA) ,sizeof(xtfloatx4));
        AE_SSX2X2_IP(x6,x7,castxcc(xtfloatx4,pA) ,sizeof(xtfloatx4));
        AE_SSX2X2_IP(AE_SEL32_LL_SX2(y1,y1),XT_CONST_S(0),castxcc(xtfloatx4,pA) ,sizeof(xtfloatx4));
    }
    pk  = (xtfloatx2 *)(B);
    pA0 = (xtfloatx2 *)(B+8);
    pA1 = (xtfloatx2 *)(B+16);
    pA2 = (xtfloatx2 *)(B+24);
    /* Gauss elimination */
    for(k=0; k<N-1; k++)
    {
        xtfloatx2 amax;
        unsigned int imax;
        /* pivoting */
        imax=k;
        amax=XT_CONST_S(0);
        /* find absolute max value in the k-th column */
        pA=(xtfloatx2*)&B[k*(N+4)+k];
        for(n=k; n<4; n++)
        {
            xtbool2 cond;
            xtfloatx2 t;
            //XT_LSXP(t,castxcc(xtfloat,pA),2*N*sizeof(xtfloat));
            {
                ae_int32x2 tt;
                AE_L32_XP(tt,castxcc(ae_int32,pA),(N+4)*sizeof(xtfloat));
                t=XT_AE_MOVXTFLOATX2_FROMINT32X2(tt);
            }
            BMAXNUMABS_SX2(cond,amax,amax,t);
            XT_MOVF_SX2(R  ,  t,cond);
            XT_MOVEQZ(imax, n, AE_MOVAB2(cond));
        }
        R = XT_RECIP_SX2(R);
    
        /* swap k-th and imax-th rows */
        pmax = (xtfloatx2 *)(B+imax*(N+4));
        AE_LSX2X2_I (Ak0,Ak1,(const xtfloatx4*)pmax, 0*sizeof(xtfloatx4));
        Ak2=XT_LSX2I (                         pmax, 1*sizeof(xtfloatx4));
        AE_LSX2X2_I (Amax0,Amax1,(const xtfloatx4*)pk, 0*sizeof(xtfloatx4));
        Amax2=XT_LSX2I (                           pk, 1*sizeof(xtfloatx4));
        AE_SSX2X2_IP(Amax0,Amax1,castxcc(xtfloatx4,pmax), sizeof(xtfloatx4));
        XT_SSX2IP   (Amax2,                        pmax,  sizeof(xtfloatx4));
        __Pragma("no_reorder")
        /* multiply k-th row by the reciprocal *
            * pivot element during swapping rows  */
        MUL_SX2X2(Ak0,Ak1,Ak0,Ak1,R,R);
        Ak2=MUL_SX2(Ak2,R);
        AE_SSX2X2_IP(Ak0,Ak1,castxcc(xtfloatx4,pk), sizeof(xtfloatx4));
        XT_SSX2IP   (Ak2,                      pk , sizeof(xtfloatx4));

        /* elimination */
        /* join forward and back substitution */

        C0 = XT_LSX((xtfloat *)pA0, k*sizeof(float32_t));
        C1 = XT_LSX((xtfloat *)pA1, k*sizeof(float32_t));
        C2 = XT_LSX((xtfloat *)pA2, k*sizeof(float32_t));

        AE_LSX2X2_I (A00,A01,(const xtfloatx4*)pA0, 0*sizeof(xtfloatx4));
        A02=XT_LSX2I    (                      pA0, 1*sizeof(xtfloatx4));
        AE_LSX2X2_I (A10,A11,(const xtfloatx4*)pA1, 0*sizeof(xtfloatx4));
        A12=XT_LSX2I    (                      pA1, 1*sizeof(xtfloatx4));
        AE_LSX2X2_I (A20,A21,(const xtfloatx4*)pA2, 0*sizeof(xtfloatx4));
        A22=XT_LSX2I    (                      pA2, 1*sizeof(xtfloatx4));

        MSUB_SX2X2(A00,A01, Ak0, Ak1, C0,C0);
        MSUB_SX2X2(A02,A12,Ak2,Ak2,C0,C1);
        MSUB_SX2X2(A10,A11, Ak0, Ak1, C1,C1);
        MSUB_SX2X2(A20,A21, Ak0, Ak1, C2,C2);
        MSUB_SX2(A22,Ak2, C2);

        AE_SSX2X2_IP (A00,A01,castxcc(xtfloatx4,pA0), sizeof(xtfloatx4));
        XT_SSX2XC    (A02,                      pA0 , sizeof(xtfloatx4));
        AE_SSX2X2_IP (A10,A11,castxcc(xtfloatx4,pA1), sizeof(xtfloatx4));
        XT_SSX2XC    (A12,                      pA1 , sizeof(xtfloatx4));
        AE_SSX2X2_IP (A20,A21,castxcc(xtfloatx4,pA2), sizeof(xtfloatx4));
        XT_SSX2XC    (A22,                      pA2 , sizeof(xtfloatx4));
    }
    // last unrolled iteration
    {
        ae_valign aY=AE_ZALIGN64();
        pmax = (xtfloatx2 *)(B+(N-1)*(N+4));
        {
            ae_int32x2 tt;
            tt=AE_L32_X((const ae_int32*)pmax,(N-1)*sizeof(float32_t));
            R=XT_AE_MOVXTFLOATX2_FROMINT32X2(tt);
        }
        R = XT_RECIP_SX2(R);
        Ak2=XT_LSX2I (pmax, 1*sizeof(xtfloatx4));
        Ak2=MUL_SX2(Ak2,R);
        C0 = XT_LSX((xtfloat *)pA0, k*sizeof(float32_t));
        C1 = XT_LSX((xtfloat *)pA1, k*sizeof(float32_t));
        C2 = XT_LSX((xtfloat *)pA2, k*sizeof(float32_t));
        A02=XT_LSX2I(pA0, 1*sizeof(xtfloatx4));
        A12=XT_LSX2I(pA1, 1*sizeof(xtfloatx4));
        A22=XT_LSX2I(pA2, 1*sizeof(xtfloatx4));
        MSUB_SX2X2(A02,A12,Ak2,Ak2,C0,C1);
        MSUB_SX2(A22,Ak2, C2);
        // copy back to the output
        XT_SASX2IP(AE_SEL32_HH_SX2(A02,A12),aY,castxcc(xtfloatx2,y));
        XT_SASX2IP(AE_SEL32_HH_SX2(A22,Ak2),aY,castxcc(xtfloatx2,y));
        AE_SA64POS_FP(aY,y);
    }
}

size_t mtx_gjelim4x4f_getScratchSize   () { const int N=4; return(N+4)*N*sizeof(float32_t);  }
#elif (HAVE_FPU)
// code for scalar FPU

void  mtx_gjelim4x4f  (void* pScr, float32_t * restrict y, const float32_t* restrict A,const float32_t * restrict x) 
{
    const int N=4;
    int k,n;
    float32_t *restrict B; // [N][N]
    float32_t *restrict C; // [N]
    float32_t *restrict T; // [N]
    const float32_t *restrict pBrd; 
          float32_t *restrict pBwr; 
    const float32_t *restrict pTrd; 
          float32_t *restrict pTwr; 

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
            __Pragma("no_reorder")
            pTrd=T;
            pTwr=T;
            for (n=0; n<N; n++)
            {
                xtfloat t;
                XT_LSIP(t,pTrd,sizeof(float32_t));
                XT_SSIP(XT_MUL_S(t,rden),pTwr,sizeof(float32_t));
            }
        }
        __Pragma("no_reorder")
        for (i=0; i<N; i++)
        {
            xtfloat t,Ti;
            Ti=T[i];
            pBwr=&B[i*N];
            pBrd=&B[i*N];
            if(i==k) continue; 
            for (n=0; n<N; n++)
            {
                XT_LSIP(t,pBrd,sizeof(xtfloat));
                XT_MSUB_S(t,B[k*N+n],Ti);
                XT_SSIP(t,pBwr,sizeof(xtfloat));
            }
            t=C[i];
            XT_MSUB_S(t,C[k],Ti);
            C[i]=t;
        }
        // pivot row
        {
            xtfloat t,Tk=T[k];
            pBwr=&B[k*N];
            pBrd=&B[k*N];
            for (n=0; n<N; n++)
            {
                XT_LSIP(t,pBrd,sizeof(xtfloat));
                t=XT_MUL_S(t,Tk);
                XT_SSIP(t,pBwr,sizeof(xtfloat));
            }
            C[k]=XT_MUL_S(C[k],Tk);
        }
    }
    // copy back to the output
    for (k=0; k<N; k++) y[k]=C[k];
}

size_t mtx_gjelim4x4f_getScratchSize   () { const int N=4; return  N*(N+2)*sizeof(float32_t);  }
#else
DISCARD_FUN(void, mtx_gjelim4x4f, (void* pScr, float32_t *y, const float32_t* A, const float32_t * x))
size_t mtx_gjelim4x4f_getScratchSize()
{
  return 0;
}
#endif
