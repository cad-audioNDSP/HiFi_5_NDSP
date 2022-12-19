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
void cmtx_inv8x8f(void * pScr,complex_float* x)
{
    const int N=8;
    int n,k;
    complex_float* A=(complex_float* )pScr;
    ae_valignx2 ax;
    xtfloatx4 * restrict pA  ;
    xtfloatx4 * restrict pAwr;

    NASSERT_ALIGN16(pScr);
    /* Set bounds of the buffer */
    WUR_AE_CBEGIN0((uintptr_t)(A));
    WUR_AE_CEND0  ((uintptr_t)(A+2*N*N));

    /* write orignial matrix together with eye diagonal matrix */
    pA  =(xtfloatx4*)x;
    pAwr=(xtfloatx4*)A;
    ax=AE_LA128_PP(pA);
    for(n=0; n<N; n++)
    {
        xtfloatx2 An0,An1,An2,An3,An4,An5,An6,An7;
        AE_LASX2X2_IP (An0,An1,ax,pA);
        AE_LASX2X2_IP (An2,An3,ax,pA);
        AE_LASX2X2_IP (An4,An5,ax,pA);
        AE_LASX2X2_IP (An6,An7,ax,pA);
        AE_SSX2X2_IP(An0,An1,pAwr,sizeof(xtfloatx4));
        AE_SSX2X2_IP(An2,An3,pAwr,sizeof(xtfloatx4));
        AE_SSX2X2_IP(An4,An5,pAwr,sizeof(xtfloatx4));
        AE_SSX2X2_IP(An6,An7,pAwr,sizeof(xtfloatx4));
        AE_SSX2X2_IP(XT_CONST_S(0),XT_CONST_S(0),pAwr,sizeof(xtfloatx4));
        AE_SSX2X2_IP(XT_CONST_S(0),XT_CONST_S(0),pAwr,sizeof(xtfloatx4));
        AE_SSX2X2_IP(XT_CONST_S(0),XT_CONST_S(0),pAwr,sizeof(xtfloatx4));
        AE_SSX2X2_IP(XT_CONST_S(0),XT_CONST_S(0),pAwr,sizeof(xtfloatx4));
        XT_SSX      (XT_CONST_S(1),(xtfloat*)pAwr,(n-N)*(int)sizeof(complex_float));
    }
    __Pragma("no_reorder")

    /* Gauss elimination */
    for(k=0; k<N; k++)
    {
        xtfloatx2 An0,An1,An2,An3,An4,An5,An6,An7;
        xtfloatx2 Ak0,Ak1,Ak2,Ak3,Ak4,Ak5,Ak6,Ak7;
        xtfloatx2 Bn0,Bn1,Bn2,Bn3,Bn4,Bn5,Bn6,Bn7;
        xtfloatx2 Bk0,Bk1,Bk2,Bk3,Bk4,Bk5,Bk6,Bk7;
        xtfloatx2 Ank,Amax,R,rden;
        int imax;
        pA=(xtfloatx4*)(A+k*2*N+k);
        Amax=XT_CONST_S(0);
        imax=k;
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
        pA  =(xtfloatx4*)(A+imax*2*N);
        pAwr=(xtfloatx4*)(A+   k*2*N);
        Amax=R;
        rden=XT_MULCCONJ_S(Amax,Amax);
        rden=MULMUX_S(rden,XT_CONST_S(1),1);    // make (re,-re)
        rden=XT_RECIP_SX2(rden);
        rden=XT_MUL_SX2(Amax,rden);
        AE_LSX2X2_I (Ak0,Ak1,pAwr, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (Ak2,Ak3,pAwr, 1*sizeof(xtfloatx4));
        AE_LSX2X2_I (Ak4,Ak5,pAwr, 2*sizeof(xtfloatx4));
        AE_LSX2X2_I (Ak6,Ak7,pAwr, 3*sizeof(xtfloatx4));
        AE_LSX2X2_I (Bk0,Bk1,pAwr, 4*sizeof(xtfloatx4));
        AE_LSX2X2_I (Bk2,Bk3,pAwr, 5*sizeof(xtfloatx4));
        AE_LSX2X2_I (Bk4,Bk5,pAwr, 6*sizeof(xtfloatx4));
        AE_LSX2X2_I (Bk6,Bk7,pAwr, 7*sizeof(xtfloatx4));
        AE_LSX2X2_I (An0,An1,pA, 0*sizeof(xtfloatx4));
        AE_LSX2X2_I (An2,An3,pA, 1*sizeof(xtfloatx4));
        AE_LSX2X2_I (An4,An5,pA, 2*sizeof(xtfloatx4));
        AE_LSX2X2_I (An6,An7,pA, 3*sizeof(xtfloatx4));
        AE_LSX2X2_I (Bn0,Bn1,pA, 4*sizeof(xtfloatx4));
        AE_LSX2X2_I (Bn2,Bn3,pA, 5*sizeof(xtfloatx4));
        AE_LSX2X2_I (Bn4,Bn5,pA, 6*sizeof(xtfloatx4));
        AE_LSX2X2_I (Bn6,Bn7,pA, 7*sizeof(xtfloatx4));
        __Pragma("no_reorder")
        AE_SSX2X2_I (Ak0,Ak1,pA, 0*sizeof(xtfloatx4));
        AE_SSX2X2_I (Ak2,Ak3,pA, 1*sizeof(xtfloatx4));
        AE_SSX2X2_I (Ak4,Ak5,pA, 2*sizeof(xtfloatx4));
        AE_SSX2X2_I (Ak6,Ak7,pA, 3*sizeof(xtfloatx4));
        AE_SSX2X2_I (Bk0,Bk1,pA, 4*sizeof(xtfloatx4));
        AE_SSX2X2_I (Bk2,Bk3,pA, 5*sizeof(xtfloatx4));
        AE_SSX2X2_I (Bk4,Bk5,pA, 6*sizeof(xtfloatx4));
        AE_SSX2X2_I (Bk6,Bk7,pA, 7*sizeof(xtfloatx4));
        __Pragma("no_reorder")
        MULMUX_SX2X2(Ak0,Ak1,An0,An1,rden,rden, 0); MADDMUXQ_S (Ak0,Ak1,An0,An1,rden, 1);
        MULMUX_SX2X2(Ak2,Ak3,An2,An3,rden,rden, 0); MADDMUXQ_S (Ak2,Ak3,An2,An3,rden, 1);
        MULMUX_SX2X2(Ak4,Ak5,An4,An5,rden,rden, 0); MADDMUXQ_S (Ak4,Ak5,An4,An5,rden, 1);
        MULMUX_SX2X2(Ak6,Ak7,An6,An7,rden,rden, 0); MADDMUXQ_S (Ak6,Ak7,An6,An7,rden, 1);
        MULMUX_SX2X2(Bk0,Bk1,Bn0,Bn1,rden,rden, 0); MADDMUXQ_S (Bk0,Bk1,Bn0,Bn1,rden, 1);
        MULMUX_SX2X2(Bk2,Bk3,Bn2,Bn3,rden,rden, 0); MADDMUXQ_S (Bk2,Bk3,Bn2,Bn3,rden, 1);
        MULMUX_SX2X2(Bk4,Bk5,Bn4,Bn5,rden,rden, 0); MADDMUXQ_S (Bk4,Bk5,Bn4,Bn5,rden, 1);
        MULMUX_SX2X2(Bk6,Bk7,Bn6,Bn7,rden,rden, 0); MADDMUXQ_S (Bk6,Bk7,Bn6,Bn7,rden, 1);
        AE_SSX2X2_I (Ak2,Ak3,pAwr, 1*sizeof(xtfloatx4));
        AE_SSX2X2_I (Ak4,Ak5,pAwr, 2*sizeof(xtfloatx4));
        AE_SSX2X2_I (Ak6,Ak7,pAwr, 3*sizeof(xtfloatx4));
        AE_SSX2X2_I (Bk0,Bk1,pAwr, 4*sizeof(xtfloatx4));
        AE_SSX2X2_I (Bk2,Bk3,pAwr, 5*sizeof(xtfloatx4));
        AE_SSX2X2_I (Bk4,Bk5,pAwr, 6*sizeof(xtfloatx4));
        AE_SSX2X2_I (Bk6,Bk7,pAwr, 7*sizeof(xtfloatx4));
        AE_SSX2X2_XC(Ak0,Ak1,pAwr, N*sizeof(xtfloatx4));
        pA=pAwr;
        // elimination
        for (n=0; n<N-1; n++)
        {
            Ank=XT_LSX2X((const xtfloatx2*)pA,k*sizeof(complex_float));
            AE_LSX2X2_I (An2,An3,pA, 1*sizeof(xtfloatx4));
            AE_LSX2X2_I (An4,An5,pA, 2*sizeof(xtfloatx4));
            AE_LSX2X2_I (An6,An7,pA, 3*sizeof(xtfloatx4));
            AE_LSX2X2_I (Bn0,Bn1,pA, 4*sizeof(xtfloatx4));
            AE_LSX2X2_I (Bn2,Bn3,pA, 5*sizeof(xtfloatx4));
            AE_LSX2X2_I (Bn4,Bn5,pA, 6*sizeof(xtfloatx4));
            AE_LSX2X2_I (Bn6,Bn7,pA, 7*sizeof(xtfloatx4));
            AE_LSX2X2_XC(An0,An1,pA, N*sizeof(xtfloatx4));
            MADDMUXQ_S(An0,An1,Ak0,Ak1,Ank,2); MADDMUXQ_S(An0,An1,Ak0,Ak1,Ank,3);
            MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,2); MADDMUXQ_S(An2,An3,Ak2,Ak3,Ank,3);
            MADDMUXQ_S(An4,An5,Ak4,Ak5,Ank,2); MADDMUXQ_S(An4,An5,Ak4,Ak5,Ank,3);
            MADDMUXQ_S(An6,An7,Ak6,Ak7,Ank,2); MADDMUXQ_S(An6,An7,Ak6,Ak7,Ank,3);
            MADDMUXQ_S(Bn0,Bn1,Bk0,Bk1,Ank,2); MADDMUXQ_S(Bn0,Bn1,Bk0,Bk1,Ank,3);
            MADDMUXQ_S(Bn2,Bn3,Bk2,Bk3,Ank,2); MADDMUXQ_S(Bn2,Bn3,Bk2,Bk3,Ank,3);
            MADDMUXQ_S(Bn4,Bn5,Bk4,Bk5,Ank,2); MADDMUXQ_S(Bn4,Bn5,Bk4,Bk5,Ank,3);
            MADDMUXQ_S(Bn6,Bn7,Bk6,Bk7,Ank,2); MADDMUXQ_S(Bn6,Bn7,Bk6,Bk7,Ank,3);
            AE_SSX2X2_I (An2,An3,pAwr, 1*sizeof(xtfloatx4));
            AE_SSX2X2_I (An4,An5,pAwr, 2*sizeof(xtfloatx4));
            AE_SSX2X2_I (An6,An7,pAwr, 3*sizeof(xtfloatx4));
            AE_SSX2X2_I (Bn0,Bn1,pAwr, 4*sizeof(xtfloatx4));
            AE_SSX2X2_I (Bn2,Bn3,pAwr, 5*sizeof(xtfloatx4));
            AE_SSX2X2_I (Bn4,Bn5,pAwr, 6*sizeof(xtfloatx4));
            AE_SSX2X2_I (Bn6,Bn7,pAwr, 7*sizeof(xtfloatx4));
            AE_SSX2X2_XC(An0,An1,pAwr, N*sizeof(xtfloatx4));
        }
        __Pragma("no_reorder")
    }
    /* copy to x */
    ax=AE_ZALIGN128();
    pA=(xtfloatx4*)(A+N);
    pAwr=(xtfloatx4*)x;
    for(n=0; n<N; n++)
    {
        xtfloatx2 An0,An1,An2,An3,An4,An5,An6,An7;
        AE_LSX2X2_IP(An0,An1,pA, 1*sizeof(xtfloatx4));
        AE_LSX2X2_IP(An2,An3,pA, 1*sizeof(xtfloatx4));
        AE_LSX2X2_IP(An4,An5,pA, 1*sizeof(xtfloatx4));
        AE_LSX2X2_XP(An6,An7,pA, (1+N/2)*sizeof(xtfloatx4));
        AE_SASX2X2_IP(An0,An1,ax,pAwr);
        AE_SASX2X2_IP(An2,An3,ax,pAwr);
        AE_SASX2X2_IP(An4,An5,ax,pAwr);
        AE_SASX2X2_IP(An6,An7,ax,pAwr);
    }
    AE_SA128POS_FP(ax,pAwr);
}
size_t cmtx_inv8x8f_getScratchSize()   { return 2*8*8*sizeof(complex_float); }
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

void cmtx_inv8x8f(void * pScr,complex_float* restrict x) 
{
    const int N=8;
    int p,n,m,k;
    complex_float r;
    complex_float  *restrict A=(complex_float* )pScr;
          complex_float *restrict  pAwr;
    const complex_float *restrict  pArd;
          complex_float *restrict  pAk;
          complex_float *restrict  pAnext;
    ae_valignx2 aX;
    WUR_AE_CBEGIN0((uintptr_t)(A));
    WUR_AE_CEND0  ((uintptr_t)(A+2*N*N));

    pArd=x;
    pAwr=A;
    aX=AE_LA128_PP(pArd);
    for(n=0; n<N; n++)
    {
        ae_int32x2 x0,x1,x2,x3,x4,x5,x6,x7;
        AE_LA32X2X2_IP(x0,x1,aX,castxcc(ae_int32x4,pArd));
        AE_LA32X2X2_IP(x2,x3,aX,castxcc(ae_int32x4,pArd));
        AE_LA32X2X2_IP(x4,x5,aX,castxcc(ae_int32x4,pArd));
        AE_LA32X2X2_IP(x6,x7,aX,castxcc(ae_int32x4,pArd));
        AE_S32X2X2_IP(x0,x1,castxcc(ae_int32x4,pAwr),sizeof(ae_int32x4));
        AE_S32X2X2_IP(x2,x3,castxcc(ae_int32x4,pAwr),sizeof(ae_int32x4));
        AE_S32X2X2_IP(x4,x5,castxcc(ae_int32x4,pAwr),sizeof(ae_int32x4));
        AE_S32X2X2_IP(x6,x7,castxcc(ae_int32x4,pAwr),sizeof(ae_int32x4));
        AE_S32X2X2_IP(0 ,0 ,castxcc(ae_int32x4,pAwr),sizeof(ae_int32x4));
        AE_S32X2X2_IP(0 ,0 ,castxcc(ae_int32x4,pAwr),sizeof(ae_int32x4));
        AE_S32X2X2_IP(0 ,0 ,castxcc(ae_int32x4,pAwr),sizeof(ae_int32x4));
        AE_S32X2X2_IP(0 ,0 ,castxcc(ae_int32x4,pAwr),sizeof(ae_int32x4));
        ((float32_t*)&pAwr[n-N])[0]=XT_CONST_S(1);
    }
    __Pragma("no_reorder")
    /* Gauss elimination */
    for(k=0; k<N; k++)
    {
        float32_t bmax;
        int imax;
        /* pivoting */
        imax=k; bmax=XT_CONST_S(0);
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
        pAk=&A[k*2*N];
        recipc(r,pAk[2*N*(imax-k)+k]);
        for(m=0; m<N; m++) 
        {
            complex_float t0,t1;
            ae_int64 s0,s1; 
            t0=pAk[2*N*(imax-k)+0]; 
            t1=pAk[2*N*(imax-k)+1]; 
            AE_L64X2_I(s0,s1,(const ae_int64x2*)pAk,0);
            mulc(t0,t0,r); 
            mulc(t1,t1,r); 
            AE_S64X2_X(s0,s1,(ae_int64x2*)pAk,2*N*(imax-k)*sizeof(complex_float));
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
    pArd=A+N;
    pAwr=x;
    aX=AE_ZALIGN128();
    for(n=0; n<N; n++)
    {
        ae_int32x2 x0,x1,x2,x3,x4,x5,x6,x7;
        AE_L32X2X2_IP(x0,x1,castxcc(ae_int32x4,pArd),sizeof(ae_int32x4));
        AE_L32X2X2_IP(x2,x3,castxcc(ae_int32x4,pArd),sizeof(ae_int32x4));
        AE_L32X2X2_IP(x4,x5,castxcc(ae_int32x4,pArd),sizeof(ae_int32x4));
        AE_L32X2X2_XP(x6,x7,castxcc(ae_int32x4,pArd),(N/2+1)*sizeof(ae_int32x4));
        AE_SA32X2X2_IP(x0,x1,aX,castxcc(ae_int32x4,pAwr));
        AE_SA32X2X2_IP(x2,x3,aX,castxcc(ae_int32x4,pAwr));
        AE_SA32X2X2_IP(x4,x5,aX,castxcc(ae_int32x4,pAwr));
        AE_SA32X2X2_IP(x6,x7,aX,castxcc(ae_int32x4,pAwr));
    }
    AE_SA128POS_FP(aX,pAwr);
}
size_t cmtx_inv8x8f_getScratchSize() { const int N=8; return 2*N*N*sizeof(complex_float); }
#else
DISCARD_FUN(void, cmtx_inv8x8f, (void* pScr, complex_float* x))
size_t cmtx_inv8x8f_getScratchSize()
{
  return 0;
}
#endif
