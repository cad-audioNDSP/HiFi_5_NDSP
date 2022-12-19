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
// specialized variant for 2x2
void  cmtx_gjelim2x2f  (void* pScr, complex_float *y, const complex_float* A,const complex_float * x) 
{
    xtfloatx2 a,b,c,d,r,x0,x1,y0,y1,ad,bc,rden;
    ae_valignx2 aA,aX,aY;
    aA=AE_LA128_PP(A);
    AE_LASX2X2_IP(a,b,aA,castxcc(xtfloatx4,A));
    AE_LASX2X2_IP(c,d,aA,castxcc(xtfloatx4,A));
    aX=AE_LA128_PP(x);
    AE_LASX2X2_IP(x0,x1,aX,castxcc(xtfloatx4,x));

    MULC_SX2(ad,bc,a,b,d,c);
    r=XT_SUB_SX2(ad,bc);
    // complex reciprocal
    rden=XT_MULCCONJ_S(r,r);
    rden=MULMUX_S(rden,XT_CONST_S(1),1);    // make (re,-re)
    rden=XT_RECIP_SX2(rden);
    r=XT_MUL_SX2(r,rden);
    // compute result via determinant
    MULC_SX2(y0,y1,d,a,x0,x1);
    MSUBC_SX2(y0,y1,b,c,x1,x0);
    MULC_SX2(y0,y1,y0,y1,r,r);
    aY=AE_ZALIGN128();
    AE_SASX2X2_IP(y0,y1,aY,castxcc(xtfloatx4,y));
    AE_SA128POS_FP(aY,y);
}
size_t cmtx_gjelim2x2f_getScratchSize   () { return 0; }

#elif (HAVE_FPU)
// code for scalar FPU

/* Complex floating-point multiplication, single precision. */
#define mulc(z,x,y) __mulc((float32_t*)&z,(const float32_t*)&x,(const float32_t*)&y)
inline_ void __mulc(float32_t* z, const float32_t* x, const float32_t* y )
{
  xtfloat re,im;
  re=XT_MUL_S(x[0],y[0]); XT_MSUB_S(re,x[1],y[1]);
  im=XT_MUL_S(x[1],y[0]); XT_MADD_S(im,x[0],y[1]);
  z[0]=re; z[1]=im;
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

#define subc(z,x,y) __subc((float32_t*)&z,(const float32_t*)&x,(const float32_t*)&y)
inline_ void __subc(float32_t* z, const float32_t* x, const float32_t* y )
{
    xtfloat re,im;
    re=XT_SUB_S(x[0],y[0]);
    im=XT_SUB_S(x[1],y[1]);
    z[0]=re;
    z[1]=im;
}

void  cmtx_gjelim2x2f  (void* pScr, complex_float *y, const complex_float* A,const complex_float * x) 
{
    complex_float a,b,c,d,r,x0,x1,y0,y1,t0,t1;
    complex_float det;
    a=A[0]; b=A[1]; c=A[2]; d=A[3];
    x0=x[0]; x1=x[1];
    mulc(t0,a,d);mulc(t1,c,b); subc(det,t0,t1);
    recipc(r,det);

    mulc(t0,d,x0);mulc(t1,b,x1); subc(y0,t0,t1); 
    mulc(t0,a,x1);mulc(t1,c,x0); subc(y1,t0,t1);
    mulc(y[0],y0,r);
    mulc(y[1],y1,r);
}
size_t cmtx_gjelim2x2f_getScratchSize   () { return 0; }

#else
DISCARD_FUN(void,cmtx_gjelim2x2f,(void* pScr, complex_float *y, const complex_float* A,const complex_float * x))
size_t cmtx_gjelim2x2f_getScratchSize        () 
{
    return 0;
}
#endif
