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
void  mtx_gjelim8x8f  (void* pScr, float32_t *y, const float32_t* A,const float32_t * x)
{
    const int N=8;
    int k,n;
    float32_t *B; // [N+2][N]
    xtfloatx4 * restrict pA;
    xtfloatx4 * restrict pAwr;

    // allocate on scratch
    B=(float32_t *)pScr;
    /* Set bounds of the buffer */
    NASSERT_ALIGN16(B);
    WUR_AE_CBEGIN0((uintptr_t)(B));
    WUR_AE_CEND0  ((uintptr_t)(B+N*(N+4)));

    // copy input
    {
        ae_valignx2 aA;
        pAwr=(xtfloatx4*)B;
        aA=AE_LA128_PP(A);
        for (k=0; k<N; k++) 
        {
            xtfloatx2 x0,x1,x2,x3,x4;
            AE_LASX2X2_IP(x0,x1,aA,castxcc(xtfloatx4,A));
            AE_LASX2X2_IP(x2,x3,aA,castxcc(xtfloatx4,A));
            {
                ae_int32x2 tt;
                AE_L32_IP(tt,castxcc(ae_int32,x),sizeof(xtfloat));
                x4=XT_AE_MOVXTFLOATX2_FROMINT32X2(tt);
            }
            AE_SSX2X2_IP(x0,x1,pAwr,sizeof(xtfloatx4));
            AE_SSX2X2_IP(x2,x3,pAwr,sizeof(xtfloatx4));
            AE_SSX2X2_IP(x4,XT_CONST_S(0),pAwr,sizeof(xtfloatx4));
        }
    }
    /* Gauss elimination */
    for(k=0; k<N; k++)
    {
        xtfloatx2 Ak0,Ak1,Ak2,Ak3,Ak4;
        xtfloatx2 Amax0,Amax1,Amax2,Amax3,Amax4;
        xtfloatx2 amax;
        xtfloatx2 R=XT_CONST_S(0);
        unsigned int imax;
        /* pivoting */
        imax=k;
        amax=XT_CONST_S(0);
        /* find absolute max value in the k-th column */
        pA=(xtfloatx4*)&B[k*(N+4)+k];
        for(n=k; n<N; n++)
        {
          xtbool2 cond;
          xtfloatx2  t;
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

        pAwr=(xtfloatx4*)&B[k*(N+4)];
        pA  =(xtfloatx4*)&B[imax*(N+4)];
        AE_LSX2X2_I (Ak0,Ak1,(const xtfloatx4*)pA, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (Ak2,Ak3,(const xtfloatx4*)pA, 1*sizeof(xtfloatx4));
        Ak4=XT_LSX2I (       (const xtfloatx2*)pA, 2*sizeof(xtfloatx4));
        AE_LSX2X2_I (Amax0,Amax1,(const xtfloatx4*)pAwr, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (Amax2,Amax3,(const xtfloatx4*)pAwr, 1*sizeof(xtfloatx4));
        Amax4=XT_LSX2I(          (const xtfloatx2*)pAwr, 2*sizeof(xtfloatx4));
        AE_SSX2X2_IP(Amax0,Amax1,castxcc(xtfloatx4,pA), sizeof(xtfloatx4));
        AE_SSX2X2_IP(Amax2,Amax3,castxcc(xtfloatx4,pA), sizeof(xtfloatx4));
        XT_SSX2IP   (Amax4      ,castxcc(xtfloatx2,pA), sizeof(xtfloatx4));
        __Pragma("no_reorder")
        /* multiply k-th row by the reciprocal *
        * pivot element during swapping rows  */
        MUL_SX2X2(Ak0,Ak1,Ak0,Ak1,R,R);
        MUL_SX2X2(Ak2,Ak3,Ak2,Ak3,R,R);
        Ak4=MUL_SX2(Ak4,R);
        AE_SSX2X2_I (Ak2,Ak3,pAwr,1*sizeof(xtfloatx4));
        XT_SSX2I    (Ak4,(xtfloatx2*)pAwr,2*sizeof(xtfloatx4));
        AE_SSX2X2_XC(Ak0,Ak1,pAwr,(N+4)*sizeof(float32_t));
        pA=pAwr;
        for (n=0; n<N-1; n++)
        {
            xtfloatx2 c,An0,An1,An2,An3,An4;
            c=XT_LSX((const xtfloat*)pA,k*sizeof(float32_t));

            AE_LSX2X2_I (An2,An3,pA,1*sizeof(xtfloatx4));
            An4=XT_LSX2I ((xtfloatx2*)pA,2*sizeof(xtfloatx4));
            AE_LSX2X2_XC(An0,An1,pA,(N+4)*sizeof(float32_t));

            MSUB_SX2X2(An0,An1,Ak0,Ak1,c,c);
            MSUB_SX2X2(An2,An3,Ak2,Ak3,c,c);
            MSUB_SX2(An4,Ak4,c);

            AE_SSX2X2_I (An2,An3,pAwr,1*sizeof(xtfloatx4));
            XT_SSX2I    (An4,(xtfloatx2*)pAwr,2*sizeof(xtfloatx4));
            AE_SSX2X2_XC(An0,An1,pAwr,(N+4)*sizeof(float32_t));
        }
        __Pragma("no_reorder")
    }
    // copy back to the output
    {
        ae_valign aY=AE_ZALIGN64();
        xtfloatx2 x0,x1;
        pA=(xtfloatx4*)(B+N);
        for (k=0; k<N/2; k++)
        {
            XT_LSX2XP(x0,castxcc(xtfloatx2,pA),(N+4)*sizeof(float32_t));
            XT_LSX2XP(x1,castxcc(xtfloatx2,pA),(N+4)*sizeof(float32_t));
            XT_SASX2IP(AE_SEL32_HH_SX2(x0,x1),aY,castxcc(xtfloatx2,y));
        }
        AE_SA64POS_FP(aY,y);
    }
}
size_t mtx_gjelim8x8f_getScratchSize   () { const int N=8; return(N+4)*N*sizeof(float32_t);  }
#elif (HAVE_FPU)
// code for scalar FPU

void  mtx_gjelim8x8f  (void* pScr, float32_t * restrict y, const float32_t* restrict A,const float32_t * restrict x) 
{
    const int N=8;
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

size_t mtx_gjelim8x8f_getScratchSize   () { const int N=8; return  N*(N+2)*sizeof(float32_t);  }
#else
DISCARD_FUN(void, mtx_gjelim8x8f, (void* pScr, float32_t *y, const float32_t* A, const float32_t * x))
size_t mtx_gjelim8x8f_getScratchSize()
{
  return 0;
}
#endif
