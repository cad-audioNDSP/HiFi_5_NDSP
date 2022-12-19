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
   Complex Matrix Gauss-Jordan Elimination for linear equation problem
   floating point 
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
void  cmtx_gjelim3x3f  (void* pScr, complex_float *y, const complex_float* A,const complex_float * x)
#if 0
{
    const int N=3;
    int k,n;
    complex_float *B; // [N][N+1]

    // allocate on scratch
    B=(complex_float *)pScr;
    // copy input
    for (k=0; k<N; k++) 
    {
        for (n=0; n<N; n++) B[k*(N+1)+n]=A[k*N+n];
        B[k*(N+1)+N]=x[k];
    }

    for (k=0; k<N; k++)
    {
        complex_float r;
        float32_t bmax;
        int p,imax,m;

        // pivoting
        imax=k; bmax=0;
        for (n=k; n<N; n++)
        {
            float32_t b=abs_fl32c(B[n*(N+1)+k]);
            if (b>=bmax) { bmax=b; imax=n; }
        }

        r=recip_fl32c(B[imax*(N+1)+k]);
        for(m=0; m<(N+1); m++) { complex_float t; t=B[imax*(N+1)+m]; B[imax*(N+1)+m]=B[k*(N+1)+m]; B[k*(N+1)+m]=mul_fl32c(t,r); }
        /* elimination */
        {
            complex_float c;
            for (p=1; p<N; p++)
            {
                n = (k+p) % N;
                c=B[n*(N+1)+k];
                for (m=0; m<(N+1); m++)
                {
                    B[n*(N+1)+m]=sub_fl32c(B[n*(N+1)+m],mul_fl32c(B[k*(N+1)+m],c));
                }
            }
        }
    }
    // copy back to the output
    for (k=0; k<N; k++) y[k]=B[k*(N+1)+N];
}
#else
{
    const int N=3;
    int n,k;
    complex_float *B; // [N][N+1]

    ae_valignx2 ax;
    xtfloatx2 *restrict pX;
    xtfloatx4 *restrict pk;
    xtfloatx4 *restrict pA;
    xtfloatx4 *restrict pA0;
    xtfloatx4 *restrict pA1;

    NASSERT_ALIGN16(pScr);
    // allocate on scratch
    B=(complex_float *)pScr;
    /* Set bounds of the buffer */
    WUR_AE_CBEGIN0((uintptr_t)(B));
    WUR_AE_CEND0  ((uintptr_t)(B+(N+1)*N));
    // copy inputs
    pX=(xtfloatx2*)x;
    pk=(xtfloatx4*)B;
    pA=(xtfloatx4*)A;
    for(n=0; n<N; n++)
    {
        xtfloatx2 An0,An1,An2,Bn;
        ax=AE_LA128_PP(pA);
        AE_LASX2X2_IP (An0,An1,ax,pA);
        XT_LSX2IP (An2,castxcc(xtfloatx2,pA),sizeof(xtfloatx2));
        XT_LSX2IP (Bn, pX, sizeof(xtfloatx2));
        AE_SSX2X2_IP(An0,An1,pk,sizeof(xtfloatx4));
        AE_SSX2X2_IP(An2,Bn ,pk,sizeof(xtfloatx4));
    }
    __Pragma("no_reorder")

    pk  = (xtfloatx4 *)(B);
    pA0 = (xtfloatx4 *)(B+1*(N+1));
    pA1 = (xtfloatx4 *)(B+2*(N+1));
    /* Gauss elimination */
    for(k=0; k<2; k++)
    {
        xtfloatx2 rden,Amax;
        xtfloatx2 Ak0,Ak1,Ak2,Ak3;
        xtfloatx2 An0,An1,An2,An3;
        xtfloatx2 Ank,R;

        int imax;
        /* pivoting */
        pA=(xtfloatx4*)(B+k*(N+1)+k);
        Amax=XT_CONST_S(0);
        imax=k;
        __Pragma("loop_count max=3")
        __Pragma("no_unroll")   // 2 or 3 iterations only
        for(n=k; n<N; n++)
        {
            xtbool2 cond;
            xtfloatx2 a;
            XT_LSX2XP(a,castxcc(xtfloatx2,pA),(N+1)*sizeof(complex_float));
            BMAXNUMABS_SX2(cond,Amax,Amax,MAXNUMABS_SX2(a,AE_SEL32_LH_SX2(a,a)));
            XT_MOVF_SX2(R  ,  a,cond);
            XT_MOVEQZ(imax, n, AE_MOVAB2(cond));
        }
        // permute and process pivot row first
        pA=(xtfloatx4*)(B+imax*(N+1));
        Amax=R;
        rden=XT_MULCCONJ_S(Amax,Amax);
        rden=MULMUX_S(rden,XT_CONST_S(1),1);    // make (re,-re)
        rden=XT_RECIP_SX2(rden);
        rden=XT_MUL_SX2(Amax,rden);
        AE_LSX2X2_I (Ak0,Ak1,pk, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (Ak2,Ak3,pk, 1*sizeof(xtfloatx4));
        AE_LSX2X2_I (An0,An1,pA, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (An2,An3,pA, 1*sizeof(xtfloatx4));
        __Pragma("no_reorder")
        AE_SSX2X2_I (Ak0,Ak1,pA, 0*sizeof(xtfloatx4));
        AE_SSX2X2_I (Ak2,Ak3,pA, 1*sizeof(xtfloatx4));
        __Pragma("no_reorder")
        MULMUX_SX2X2(Ak0,Ak1,An0,An1,rden,rden, 0); MADDMUXQ_S (Ak0,Ak1,An0,An1,rden, 1);
        MULMUX_SX2X2(Ak2,Ak3,An2,An3,rden,rden, 0); MADDMUXQ_S (Ak2,Ak3,An2,An3,rden, 1);
        AE_SSX2X2_I (Ak2,Ak3,pk, 1*sizeof(xtfloatx4));
        AE_SSX2X2_XC(Ak0,Ak1,pk, (N+1)*sizeof(xtfloatx2));
        /* elimination in other rows */
        Ank=XT_LSX2X((const xtfloatx2*)pA0,k*sizeof(complex_float));
        AE_LSX2X2_I (An0,An1,pA0, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (An2,An3,pA0, 1*sizeof(xtfloatx4));
        MADDMUXQ_S(An0,An1,Ak0,Ak1,Ank,2); MADDMUXQ_S(An0,An1,Ak0,Ak1,Ank,3);
        MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,2); MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,3);
        AE_SSX2X2_I (An2,An3,pA0, 1*sizeof(xtfloatx4));
        AE_SSX2X2_XC(An0,An1,pA0, (N+1)*sizeof(xtfloatx2));
        Ank=XT_LSX2X((const xtfloatx2*)pA1,k*sizeof(complex_float));
        AE_LSX2X2_I (An0,An1,pA1, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (An2,An3,pA1, 1*sizeof(xtfloatx4));
        MADDMUXQ_S(An0,An1,Ak0,Ak1,Ank,2); MADDMUXQ_S(An0,An1,Ak0,Ak1,Ank,3);
        MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,2); MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,3);
        AE_SSX2X2_I (An2,An3,pA1, 1*sizeof(xtfloatx4));
        AE_SSX2X2_XC(An0,An1,pA1, (N+1)*sizeof(xtfloatx2));
        __Pragma("no_reorder")
    }
    // last simpler iteration together with saving the output
    pX=(xtfloatx2*)y;
    {
        xtfloatx2 rden,Amax;
        xtfloatx2 Ak2,Ak3;
        xtfloatx2 An0,An1,An2,An3;
        xtfloatx2 Ank;

        AE_LSX2X2_I (An0,An1,pk, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (An2,An3,pk, 1*sizeof(xtfloatx4));
        // permute and process pivot row first
        Amax=XT_LSX2I((const xtfloatx2*)pk,2*sizeof(xtfloatx2));
        rden=XT_MULCCONJ_S(Amax,Amax);
        rden=MULMUX_S(rden,XT_CONST_S(1),1);    // make (re,-re)
        rden=XT_RECIP_SX2(rden);
        rden=XT_MUL_SX2(Amax,rden);
        MULMUX_SX2X2(Ak2,Ak3,An2,An3,rden,rden, 0); MADDMUXQ_S (Ak2,Ak3,An2,An3,rden, 1);
        XT_SSX2I(Ak3,pX,2*sizeof(complex_float));
        /* elimination in other rows */
        Ank=XT_LSX2X((const xtfloatx2*)pA0,k*sizeof(complex_float));
        AE_LSX2X2_I (An2,An3,pA0, 1*sizeof(xtfloatx4));
        MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,2); MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,3);
        XT_SSX2I(An3,pX,0*sizeof(complex_float));
        Ank=XT_LSX2X((const xtfloatx2*)pA1,k*sizeof(complex_float));
        AE_LSX2X2_I (An2,An3,pA1, 1*sizeof(xtfloatx4));
        MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,2); MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,3);
        XT_SSX2I(An3,pX,1*sizeof(complex_float));
    }
}
#endif

size_t cmtx_gjelim3x3f_getScratchSize   () { return ((3+1)*3)*sizeof(complex_float);  }
#elif (HAVE_FPU)
// code for scalar FPU

/* pseudoabs of complex number, single precision. */
#define absc(z,x) __absc((float32_t*)&z,(const float32_t*)&x)
inline_ void __absc(float32_t* z, const float32_t* x)
{
    float32_t re=XT_ABS_S(x[0]), im=XT_ABS_S(x[1]);
    z[0]=XT_MAX_S(re,im);
}

/* Complex floating-point multiplication, single precision. */
#define mulc(z,x,y) __mulc((float32_t*)&z,(const float32_t*)&x,(const float32_t*)&y)
inline_ void __mulc(float32_t* z, const float32_t* x, const float32_t* y )
{
    xtfloat re,im;
    re=XT_MUL_S(x[0],y[0]); XT_MSUB_S(re,x[1],y[1]);
    im=XT_MUL_S(x[1],y[0]); XT_MADD_S(im,x[0],y[1]);
    z[0]=re; 
    z[1]=im;
}

#define msubc(z,x,y) __msubc((float32_t*)&z,(const float32_t*)&x,(const float32_t*)&y)
inline_ void __msubc(float32_t* z, const float32_t* x, const float32_t* y )
{
    float32_t re=z[0];
    float32_t im=z[1];
    XT_MSUB_S(re,x[0],y[0]); XT_MADD_S(re,x[1],y[1]);
    XT_MSUB_S(im,x[0],y[1]); XT_MSUB_S(im,x[1],y[0]); 
    z[0]=re;
    z[1]=im;
}

#define recipc(y,x) __recipc((float32_t*)&y,(const float32_t*)&x);
inline_ void __recipc( float32_t* y, const float32_t* x)
{
    xtfloat c,re,im;
    c=XT_MUL_S(x[0],x[0]); XT_MADD_S(c,x[1],x[1]);
    c=XT_RECIP_S(c);
    re=XT_MUL_S(x[0],c);
    im=XT_MUL_S(x[1],XT_NEG_S(c));
    y[0]=re;
    y[1]=im;
}

void  cmtx_gjelim3x3f(void* pScr, complex_float *restrict y, const complex_float* restrict A,const complex_float * restrict x) 
{
    const int N=3;
    int k,n;
    complex_float *restrict B; // [N][N]
    complex_float *restrict C; // [N]
    complex_float *restrict T; // [N]
    const complex_float * restrict pBk ; 
    const complex_float * restrict pBrd; 
          complex_float * restrict pBwr; 

    // allocate on scratch
    B=(complex_float *)pScr;
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
        imax=k; bmax=XT_CONST_S(0);
        for (n=k; n<N; n++)
        {
            xtbool cond;
            float32_t b;
            absc(b,B[n*N+k]);
            cond=XT_OLE_S(bmax,b);
            XT_MOVT_S(bmax,  b, cond);
            XT_MOVT  (imax,  n, cond);
        }
        for (n=0; n<N; n++) 
        {
            complex_float t;
            t=B[imax*N+n]; B[imax*N+n]=B[k*N+n]; B[k*N+n]=t;
        }
        {
            complex_float t;
            t=C[imax    ]; C[imax    ]=C[k    ]; C[k]=t;
        }
        // find normalization factor
        {
            complex_float rden,tk;
            for (n=0; n<N; n++) T[n]=B[n*N+k];
            tk=T[k]; 
            ((xtfloat*)(T+k))[0]=XT_CONST_S(1); 
            ((xtfloat*)(T+k))[1]=XT_CONST_S(0); 
            recipc(rden,tk); // reciprocal of tk
            for (n=0; n<N; n++)
            {
                mulc(T[n],T[n],rden);
            }
        }
        __Pragma("no_reorder")
        for (i=0; i<N; i++)
        {
            complex_float Ti=T[i];
            if(i==k) continue;
            pBrd=B+i*N;
            pBwr=B+i*N;
            pBk =B+k*N;
            for (n=0; n<N; n++)
            {
                complex_float t=*pBrd++;
                complex_float bk=*pBk++;
                msubc(t,bk,Ti);
                *pBwr++=t;
            }
            msubc(C[i],C[k],Ti);
        }
        {
            complex_float t,Tk=T[k];
            pBrd=B+k*N;
            pBwr=B+k*N;
            for (n=0; n<N; n++)
            {
                t=*pBrd++;
                mulc(t,B[k*N+n],Tk);
                *pBwr++=t;
            }
            mulc(C[k],C[k],Tk);
        }
    }
    // copy back to the output
    for (k=0; k<N; k++) y[k]=C[k];
}
size_t cmtx_gjelim3x3f_getScratchSize   () { const int N=3; return (N+2)*N*sizeof(complex_float);  }

#else
DISCARD_FUN(void, cmtx_gjelim3x3f, (void* pScr, complex_float *y, const complex_float* A, const complex_float * x))
size_t cmtx_gjelim3x3f_getScratchSize()
{
  return 0;
}
#endif
