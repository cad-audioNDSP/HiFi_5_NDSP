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
 * Complex Matrix Inversion
 * floating point API
 * code optimized for HiFi5
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
  These functions implement in-place matrix inversion by Gauss elimination 
  with full pivoting
  NOTE: user may detect "invalid" or "divide-by-zero" exception in the CPU 
  flags which MAY indicate that inversion results are not accurate. Also 
  it's responsibility of the user to provide valid input matrix for 
  inversion.
  Fixed point version takes representation of input matrix and forms 
  representation of output matrix with proper scaling.

  Precision: 
  f     floating point
  32x32 32-bit input, 32-bit output

  Input:
  x[N*N]      input matrix
  qX          input matrix representation (for fixed point API only)
  Output:
  x[N*N]      result
  Temporary:
  pScr        scratch memory. Size in bytes is defined by corresponding 
              scratch allocation function 
  N is 2,3,4,6,8,10

  Returned value: floating functions return none, fixed point functions 
                  return fixed-point representation of inverted matrix
  Restrictions:
  none
-------------------------------------------------------------------------*/
void cmtx_inv3x3f(void * pScr,complex_float* x)
{
    const int N=3;
    int n,k;
    complex_float* A=(complex_float* )pScr;

    ae_valignx2 ax;
    xtfloatx4 *restrict pk;
    xtfloatx4 *restrict pA;
    xtfloatx4 *restrict pA0;
    xtfloatx4 *restrict pA1;

    NASSERT_ALIGN16(pScr);

    /* Set bounds of the buffer */
    WUR_AE_CBEGIN0((uintptr_t)(A));
    WUR_AE_CEND0  ((uintptr_t)(A+2*N*N));
    /* write orignial matrix together with eye diagonal matrix */
    pk=(xtfloatx4*)x;
    pA=(xtfloatx4*)A;
    for(n=0; n<N; n++)
    {
        xtfloatx2 An0,An1,An2;
        ax=AE_LA128_PP(pk);
        AE_LASX2X2_IP (An0,An1,ax,pk);
        XT_LSX2IP (An2,castxcc(xtfloatx2,pk),sizeof(xtfloatx2));
        AE_SSX2X2_IP(An0,An1,pA,sizeof(xtfloatx4));
        AE_SSX2X2_IP(An2,XT_CONST_S(0),pA,sizeof(xtfloatx4));
        AE_SSX2X2_IP(XT_CONST_S(0),XT_CONST_S(0),pA,sizeof(xtfloatx4));
        XT_SSX      (XT_CONST_S(1),(xtfloat*)pA,(n-N)*(int)sizeof(complex_float));
    }
    __Pragma("no_reorder")

    pk  = (xtfloatx4 *)(A);
    pA0 = (xtfloatx4 *)(A+2*1*N);
    pA1 = (xtfloatx4 *)(A+2*2*N);
    /* Gauss elimination */
    for(k=0; k<2; k++)
    {
        xtfloatx2 rden,Amax;
        xtfloatx2 Ak0,Ak1,Ak2,Ak3,Ak4,Ak5;
        xtfloatx2 An0,An1,An2,An3,An4,An5;
        xtfloatx2 Ank,R;

        int imax;
        /* pivoting */
        pA=(xtfloatx4*)(A+k*2*N+k);
        Amax=XT_CONST_S(0);
        imax=k;
        __Pragma("loop_count max=3")
        __Pragma("no_unroll")   // 2 or 3 iterations only
        for(n=k; n<N; n++)
        {
            xtbool2 cond;
            xtfloatx2 a;
            XT_LSX2XP(a,castxcc(xtfloatx2,pA),2*N*sizeof(complex_float));
            BMAXNUMABS_SX2(cond,Amax,Amax,MAXNUMABS_SX2(a,AE_SEL32_LH_SX2(a,a)));
            XT_MOVF_SX2(R  ,  a,cond);
            XT_MOVEQZ(imax, n, AE_MOVAB2(cond));
        }
        // permute and process pivot row first
        pA=(xtfloatx4*)(A+imax*2*N);
        Amax=R;
        rden=XT_MULCCONJ_S(Amax,Amax);
        rden=MULMUX_S(rden,XT_CONST_S(1),1);    // make (re,-re)
        rden=XT_RECIP_SX2(rden);
        rden=XT_MUL_SX2(Amax,rden);
        AE_LSX2X2_I (Ak0,Ak1,pk, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (Ak2,Ak3,pk, 1*sizeof(xtfloatx4));
        AE_LSX2X2_I (Ak4,Ak5,pk, 2*sizeof(xtfloatx4));
        AE_LSX2X2_I (An0,An1,pA, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (An2,An3,pA, 1*sizeof(xtfloatx4));
        AE_LSX2X2_I (An4,An5,pA, 2*sizeof(xtfloatx4));
        __Pragma("no_reorder")
        AE_SSX2X2_I (Ak0,Ak1,pA, 0*sizeof(xtfloatx4));
        AE_SSX2X2_I (Ak2,Ak3,pA, 1*sizeof(xtfloatx4));
        AE_SSX2X2_I (Ak4,Ak5,pA, 2*sizeof(xtfloatx4));
        __Pragma("no_reorder")
        MULMUX_SX2X2(Ak0,Ak1,An0,An1,rden,rden, 0); MADDMUXQ_S (Ak0,Ak1,An0,An1,rden, 1);
        MULMUX_SX2X2(Ak2,Ak3,An2,An3,rden,rden, 0); MADDMUXQ_S (Ak2,Ak3,An2,An3,rden, 1);
        MULMUX_SX2X2(Ak4,Ak5,An4,An5,rden,rden, 0); MADDMUXQ_S (Ak4,Ak5,An4,An5,rden, 1);
        AE_SSX2X2_I (Ak2,Ak3,pk, 1*sizeof(xtfloatx4));
        AE_SSX2X2_I (Ak4,Ak5,pk, 2*sizeof(xtfloatx4));
        AE_SSX2X2_XC(Ak0,Ak1,pk, N*sizeof(xtfloatx4));
        /* elimination in other rows */
        Ank=XT_LSX2X((const xtfloatx2*)pA0,k*sizeof(complex_float));
        AE_LSX2X2_I (An0,An1,pA0, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (An2,An3,pA0, 1*sizeof(xtfloatx4));
        AE_LSX2X2_I (An4,An5,pA0, 2*sizeof(xtfloatx4));
        MADDMUXQ_S(An0,An1,Ak0,Ak1,Ank,2); MADDMUXQ_S(An0,An1,Ak0,Ak1,Ank,3);
        MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,2); MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,3);
        MADDMUXQ_S(An4,An5,Ak4,Ak5,Ank,2); MADDMUXQ_S(An4,An5,Ak4,Ak5,Ank,3);
        AE_SSX2X2_I (An2,An3,pA0, 1*sizeof(xtfloatx4));
        AE_SSX2X2_I (An4,An5,pA0, 2*sizeof(xtfloatx4));
        AE_SSX2X2_XC(An0,An1,pA0, N*sizeof(xtfloatx4));
        Ank=XT_LSX2X((const xtfloatx2*)pA1,k*sizeof(complex_float));
        AE_LSX2X2_I (An0,An1,pA1, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (An2,An3,pA1, 1*sizeof(xtfloatx4));
        AE_LSX2X2_I (An4,An5,pA1, 2*sizeof(xtfloatx4));
        MADDMUXQ_S(An0,An1,Ak0,Ak1,Ank,2); MADDMUXQ_S(An0,An1,Ak0,Ak1,Ank,3);
        MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,2); MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,3);
        MADDMUXQ_S(An4,An5,Ak4,Ak5,Ank,2); MADDMUXQ_S(An4,An5,Ak4,Ak5,Ank,3);
        AE_SSX2X2_I (An2,An3,pA1, 1*sizeof(xtfloatx4));
        AE_SSX2X2_I (An4,An5,pA1, 2*sizeof(xtfloatx4));
        AE_SSX2X2_XC(An0,An1,pA1, N*sizeof(xtfloatx4));
        __Pragma("no_reorder")
    }
    // last simpler iteration: no pivoting there
    {
        xtfloatx2 rden,Amax;
        xtfloatx2 Ak2,Ak3,Ak4,Ak5;
        xtfloatx2 An0,An1,An2,An3,An4,An5;
        xtfloatx2 Ank;

        AE_LSX2X2_I (An0,An1,pk, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (An2,An3,pk, 1*sizeof(xtfloatx4));
        AE_LSX2X2_I (An4,An5,pk, 2*sizeof(xtfloatx4));
        // permute and process pivot row first
        Amax=XT_LSX2I((const xtfloatx2*)pk,2*sizeof(xtfloatx2));
        rden=XT_MULCCONJ_S(Amax,Amax);
        rden=MULMUX_S(rden,XT_CONST_S(1),1);    // make (re,-re)
        rden=XT_RECIP_SX2(rden);
        rden=XT_MUL_SX2(Amax,rden);
        MULMUX_SX2X2(Ak2,Ak3,An2,An3,rden,rden, 0); MADDMUXQ_S (Ak2,Ak3,An2,An3,rden, 1);
        MULMUX_SX2X2(Ak4,Ak5,An4,An5,rden,rden, 0); MADDMUXQ_S (Ak4,Ak5,An4,An5,rden, 1);
        AE_SSX2X2_I (Ak2,Ak3,pk, 1*sizeof(xtfloatx4));
        AE_SSX2X2_I (Ak4,Ak5,pk, 2*sizeof(xtfloatx4));
        /* elimination in other rows */
        Ank=XT_LSX2X((const xtfloatx2*)pA0,k*sizeof(complex_float));
        AE_LSX2X2_I (An2,An3,pA0, 1*sizeof(xtfloatx4));
        AE_LSX2X2_I (An4,An5,pA0, 2*sizeof(xtfloatx4));
        MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,2); MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,3);
        MADDMUXQ_S(An4,An5,Ak4,Ak5,Ank,2); MADDMUXQ_S(An4,An5,Ak4,Ak5,Ank,3);
        AE_SSX2X2_I (An2,An3,pA0, 1*sizeof(xtfloatx4));
        AE_SSX2X2_I (An4,An5,pA0, 2*sizeof(xtfloatx4));
        Ank=XT_LSX2X((const xtfloatx2*)pA1,k*sizeof(complex_float));
        AE_LSX2X2_I (An2,An3,pA1, 1*sizeof(xtfloatx4));
        AE_LSX2X2_I (An4,An5,pA1, 2*sizeof(xtfloatx4));
        MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,2); MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,3);
        MADDMUXQ_S(An4,An5,Ak4,Ak5,Ank,2); MADDMUXQ_S(An4,An5,Ak4,Ak5,Ank,3);
        AE_SSX2X2_I (An2,An3,pA1, 1*sizeof(xtfloatx4));
        AE_SSX2X2_I (An4,An5,pA1, 2*sizeof(xtfloatx4));
        __Pragma("no_reorder")
    }
    /* copy to x */
    ax=AE_ZALIGN128();
    pA=(xtfloatx4*)(A+4);
    pk=(xtfloatx4*)x;
    {
        xtfloatx2 An0,An1,An2,An3,An4,An5,An6,An7,An8;
        An0=XT_LSX2I((const xtfloatx2*)pA,-1*(int)sizeof(xtfloatx2));
        AE_LSX2X2_IP(An1,An2,pA,  3*sizeof(xtfloatx4));
        An3=XT_LSX2I((const xtfloatx2*)pA,-1*(int)sizeof(xtfloatx2));
        AE_LSX2X2_IP(An4,An5,pA,  3*sizeof(xtfloatx4));
        An6=XT_LSX2I((const xtfloatx2*)pA,-1*(int)sizeof(xtfloatx2));
        AE_LSX2X2_IP(An7,An8,pA,  3*sizeof(xtfloatx4));

        AE_SASX2X2_IP(An0,An1,ax,pk);
        AE_SASX2X2_IP(An2,An3,ax,pk);
        AE_SASX2X2_IP(An4,An5,ax,pk);
        AE_SASX2X2_IP(An6,An7,ax,pk);
        AE_SA128POS_FP(ax,pk);
        XT_SSX2I(An8,(xtfloatx2*)pk,0);
    }
}
size_t cmtx_inv3x3f_getScratchSize() { return 18*sizeof(complex_float); }
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

#define makec(z,re,im) { ((float32_t*)&(z))[0]=re;((float32_t*)&(z))[1]=im; }

void cmtx_inv3x3f(void * pScr,complex_float* restrict x) 
{
    const int N=3;
    int p,n,m,k;
    complex_float r;
    complex_float  *restrict A=(complex_float* )pScr;
          complex_float *restrict  pAwr;
    const complex_float *restrict  pArd;
          complex_float *restrict  pAk;
          complex_float *restrict  pAnext;
    WUR_AE_CBEGIN0((uintptr_t)(A));
    WUR_AE_CEND0  ((uintptr_t)(A+2*N*N));

    for(n=0; n<N; n++)
    {
        for (m=0; m<N; m++)
        {
            A[n*2*N+m]=x[n*N+m];
            makec(A[n*2*N+m+N],0.f,0.f);
        }
        makec(A[n*2*N+n+N],1.f,0.f);
    }
    __Pragma("no_reorder")
    /* Gauss elimination */
    for(k=0; k<N; k++)
    {
        float32_t bmax;
        int imax;
        /* pivoting */
        imax=k; bmax=0.f;
        pArd=&A[k];
        for(n=k; n<N; n++)
        {
            xtbool cond;
            float32_t b;
            absc(b,pArd[n*2*N]);
            cond=XT_OLE_S(bmax,b);
            XT_MOVT_S(bmax,  b, cond);
            XT_MOVT  (imax,  n, cond);
        }
        __Pragma("no_reorder")
        pAk=&A[k*2*N];
        pAwr=&A[imax*2*N];
        recipc(r,pAwr[k]);
        for(m=0; m<N; m++) 
        {
            complex_float t0,t1,s0,s1; 
            t0=pAk[2*(imax-k)*N+0]; 
            t1=pAk[2*(imax-k)*N+1]; 
            s0=pAk[0]; 
            s1=pAk[1]; 
            mulc(t0,t0,r); 
            mulc(t1,t1,r); 
            pAk[2*(imax-k)*N+0]=s0;
            pAk[2*(imax-k)*N+1]=s1;
            *pAk++=t0;
            *pAk++=t1;
        }
        __Pragma("no_reorder")
        /* elimination */
        pAnext=A+k*2*N;
        for (p=1; p<N; p++)
        {
            complex_float c,t0,t1,ak0,ak1;
            AE_ADDCIRC_XC(castxcc(ae_int64,pAnext),2*N*sizeof(complex_float));
            pAwr=pAnext;
            pArd=pAwr;
            pAk=A+k*2*N;
            c=pArd[k];
            for (m=0; m<N; m++)
            {
                t0=*pArd++; t1=*pArd++;
                ak0=*pAk++; ak1=*pAk++;
                msubc(t0,ak0,c);
                msubc(t1,ak1,c);
                *pAwr++=t0;
                *pAwr++=t1;
            }
        }
        __Pragma("no_reorder")
    }
    /* copy to x */
    for(n=0; n<N; n++)
    {
        ((ae_int64*)x)[n*N+0]=((const ae_int64*)A)[n*2*N+0+N];
        ((ae_int64*)x)[n*N+1]=((const ae_int64*)A)[n*2*N+1+N];
        ((ae_int64*)x)[n*N+2]=((const ae_int64*)A)[n*2*N+2+N];
    }
}
size_t cmtx_inv3x3f_getScratchSize() { const int N=3; return 2*N*N*sizeof(complex_float); }
#else
DISCARD_FUN(void, cmtx_inv3x3f, (void* pScr, complex_float* x))
size_t cmtx_inv3x3f_getScratchSize()
{
  return 0;
}
#endif
