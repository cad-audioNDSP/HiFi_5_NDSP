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
// specialized variant for 2x2
void  mtx_gjelim2x2f  (void* pScr, float32_t *y, const float32_t* A,const float32_t * x) 
#if 1
{
    xtfloat a,b,c,d,r,x0,x1,y0,y1;
    xtfloat det;
    a=XT_LSI((const xtfloat*)A,0*sizeof(float32_t)); 
    b=XT_LSI((const xtfloat*)A,1*sizeof(float32_t)); 
    c=XT_LSI((const xtfloat*)A,2*sizeof(float32_t)); 
    d=XT_LSI((const xtfloat*)A,3*sizeof(float32_t));
    x0=XT_LSI((const xtfloat*)x,0*sizeof(float32_t)); 
    x1=XT_LSI((const xtfloat*)x,1*sizeof(float32_t)); 
    det=XT_MUL_S(a,d);XT_MSUB_S(det,c,b);
    r=XT_RECIP_S(det);
    y0=XT_MUL_S(d,x0); 
    y1=XT_MUL_S(a,x1); 
    XT_MSUB_S(y0,b,x1);
    XT_MSUB_S(y1,c,x0);
    y0=XT_MUL_S(y0,r);
    y1=XT_MUL_S(y1,r);
    XT_SSI(y0,(xtfloat*)y,0*sizeof(float32_t));          
    XT_SSI(y1,(xtfloat*)y,1*sizeof(float32_t));
}
#else
{
    ae_valignx2 aA;
    ae_valign aX,aY;
    xtfloatx2 ab,cd,x01,x10,r,da,bc,y01;
    aA=AE_LA128_PP(A);
    aX=AE_LA64_PP(x);
    aY=AE_ZALIGN64();
    AE_LASX2X2_IP(ab,cd,aA,castxcc(xtfloatx4,A));
    AE_LASX2IP(x01,aX,castxcc(xtfloatx2,x));
    x10=AE_SEL32_LH_SX2(x01,x01);
    da=AE_SEL32_LH_SX2(cd,ab);
    bc=AE_SEL32_LH_SX2(ab,cd);

    r=XT_MULMUX_S(ab,AE_SEL32_LL_SX2(cd,cd),0); // ad,ad
    XT_MADDMUX_S(r,cd,AE_SEL32_LL_SX2(ab,ab),2); // ad-cb,ad-cb
    r=XT_RECIP_SX2(r);
    y01=XT_MUL_SX2(da,x01);
    XT_MSUB_SX2(y01,bc,x10);
    y01=XT_MUL_SX2(y01,r);
    XT_SASX2IP(y01,aY,castxcc(xtfloatx2,y));
    AE_SA64POS_FP(aY,y);
}
#endif

size_t mtx_gjelim2x2f_getScratchSize   () { return 0; }

#elif (HAVE_FPU)
// for scalar FPU
void  mtx_gjelim2x2f  (void* pScr, float32_t *y, const float32_t* A,const float32_t * x) 
{
    xtfloat a,b,c,d,r,x0,x1,y0,y1;
    xtfloat det;
    a=XT_LSI((const xtfloat*)A,0*sizeof(float32_t)); 
    b=XT_LSI((const xtfloat*)A,1*sizeof(float32_t)); 
    c=XT_LSI((const xtfloat*)A,2*sizeof(float32_t)); 
    d=XT_LSI((const xtfloat*)A,3*sizeof(float32_t));
    x0=XT_LSI((const xtfloat*)x,0*sizeof(float32_t)); 
    x1=XT_LSI((const xtfloat*)x,1*sizeof(float32_t)); 
    det=XT_MUL_S(a,d);XT_MSUB_S(det,c,b);
    r=XT_RECIP_S(det);
    y0=XT_MUL_S(d,x0); 
    y1=XT_MUL_S(a,x1); 
    XT_MSUB_S(y0,b,x1);
    XT_MSUB_S(y1,c,x0);
    y0=XT_MUL_S(y0,r);
    y1=XT_MUL_S(y1,r);
    XT_SSI(y0,(xtfloat*)y,0*sizeof(float32_t));          
    XT_SSI(y1,(xtfloat*)y,1*sizeof(float32_t));
}

size_t mtx_gjelim2x2f_getScratchSize   () { return 0; }
#else
DISCARD_FUN(void, mtx_gjelim2x2f, (void* pScr, float32_t *y, const float32_t* A, const float32_t * x))
size_t mtx_gjelim2x2f_getScratchSize()
{
  return 0;
}
#endif
