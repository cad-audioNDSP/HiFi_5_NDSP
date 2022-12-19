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
void cmtx_inv2x2f(void * pScr,complex_float* x)
{
    xtfloatx4* pXwr=(xtfloatx4*)x;
    ae_valignx2 ax;
    (void)pScr;
    xtfloatx2 a,b,c,d,r,ad,bc,rden,mr;
    ax=AE_LA128_PP(x);
    AE_LASX2X2_IP(a,b,ax,castxcc(xtfloatx4,x));
    AE_LASX2X2_IP(c,d,ax,castxcc(xtfloatx4,x));
    // determinant
    MULC_SX2(ad,bc,a,b,d,c);
    r=XT_SUB_SX2(ad,bc);
    // complex reciprocal
    rden=XT_MULCCONJ_S(r,r);
    rden=MULMUX_S(rden,XT_CONST_S(1),1);    // make (re,-re)
    rden=XT_RECIP_SX2(rden);
    r=XT_MUL_SX2(r,rden);
    mr=XT_NEG_SX2(r);
    // scaling
    MULC_SX2(d,b,d,b,r,mr);
    MULC_SX2(c,a,c,a,mr,r);
    ax=AE_ZALIGN128();
    AE_SASX2X2_IP(d,b,ax,pXwr);
    AE_SASX2X2_IP(c,a,ax,pXwr);
    AE_SA128POS_FP(ax,pXwr);
}

size_t cmtx_inv2x2f_getScratchSize() { return 0; }
#elif (HAVE_FPU)
// code for scalar FPU

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

void cmtx_inv2x2f(void * pScr,complex_float* x)
{
    complex_float a,b,c,d,r,t;
    a=x[0]; b=x[1];
    c=x[2]; d=x[3];
    mulc(t,a,d); msubc(t,b,c);
    recipc(r,t);
    ((float32_t*)&t)[0]=XT_NEG_S(((float32_t*)&r)[0]);
    ((float32_t*)&t)[1]=XT_NEG_S(((float32_t*)&r)[1]);
    mulc(x[0],d,r); mulc(x[1],b,t);
    mulc(x[2],c,t); mulc(x[3],a,r) ;
}

size_t cmtx_inv2x2f_getScratchSize() { return 0; }
#else
DISCARD_FUN(void, cmtx_inv2x2f, (void* pScr, complex_float* x))
size_t cmtx_inv2x2f_getScratchSize()
{
  return 0;
}
#endif
