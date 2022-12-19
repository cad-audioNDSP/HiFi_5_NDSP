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
  NatureDSP Signal Processing Library. Vector Operations
    Vector Dot product
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_vector.h"
#include "common.h"
#include "common_fpu.h"

#if (HAVE_VFPU==0 && HAVE_FPU==0)
DISCARD_FUN(float32_t,vec_dotf, (const float32_t * restrict x,const float32_t * restrict y,int N))
#elif (HAVE_VFPU)

/*===========================================================================
  Vector matematics:
  vec_dot              Vector Dot Product
===========================================================================*/
/*-------------------------------------------------------------------------
  Vector Dot product
  These routines take two vectors and calculates their dot product.
  Two versions of routines are available: regular versions (vec_dot64x32,
  vec_dot64x64, vec_dot64x64i, vec_dot32x16, vec_dot32x32,vec_dot16x16, 
  vec_dotf) work with arbitrary arguments, faster versions (vec_dot64x32_fast, 
  vec_dot64x64_fast, vec_dot64x64i_fast, vec_dot32x16_fast, vec_dot32x32_fast,
  vec_dot16x16_fast) apply some restrictions.  
  NOTE:
  vec_dot16x16_fast utilizes 32-bit saturating accumulator, so input data 
  should be scaled properly to avoid erroneous results.

  Precision: 
  64x32  64x32-bit data, 64-bit output (fractional multiply Q63xQ31->Q63)
  64x64  64x64-bit data, 64-bit output (fractional multiply Q63xQ63->Q63)
  64x64i 64x64-bit data, 64-bit output (low 64 bit of integer multiply)
  32x32  32x32-bit data, 64-bit output
  32x16  32x16-bit data, 64-bit output
  16x16  16x16-bit data, 64-bit output for regular version and 32-bit for 
                        fast version
  f      single precision floating point

  Input:
  x[N]  input data, Q15, Q31, Q63 or floating point
  y[N]  input data, Q15, Q31, Q63 or floating point
  N	    length of vectors
  Returns:
  dot product of all data pairs, Q31, Q63 or floating point

  Restrictions:
  Regular versions:
    none
  Faster versions:
    x,y - aligned on 16-byte boundary
    N   - multiple of 4
-------------------------------------------------------------------------*/
#define sz_f32    (int)sizeof(float32_t)
float32_t vec_dotf   (const float32_t * restrict x,const float32_t * restrict y,int N)
{
    static const int32_t ALIGN(16) seq[]={0,1,2,3,4,5,6,7};
    xtbool2 b01,b23,b45,b67;
    ae_int32x2 s01,s23,s45,s67;
    const xtfloatx4* restrict pX0;
    const xtfloatx4* restrict pY0;
    xtfloatx2 x0,x1,x2,x3;
    xtfloatx2 y0,y1,y2,y3;
    xtfloatx2 a0,a1,a2,a3;
    xtfloatx2 a4,a5,a6,a7;
    ae_valignx2 ax0,ay0;
    int n,N0;
    if(N<=0) return XT_CONST_S(0);
    if (N<=7)
    {
        xtfloat a0;
        a0=XT_CONST_S(0);
        __Pragma("no_unroll")
        for (n=0; n<N; n++) 
        {
            xtfloat x0,y0;
            XT_LSIP(x0,castxcc(xtfloat,x),sizeof(float32_t));
            XT_LSIP(y0,castxcc(xtfloat,y),sizeof(float32_t));
            XT_MADD_S(a0,x0,y0);
        }
        return (float32_t)a0;
    }
    AE_L32X2X2_I(s01,s23,(const ae_int32x4*)seq,0*sizeof(ae_int32x4));
    AE_L32X2X2_I(s45,s67,(const ae_int32x4*)seq,1*sizeof(ae_int32x4));
    N0=((N-1)&7)+1;
    b01=AE_LT32(s01,N0);    // mask unnessesary elements on the first iteration
    b23=AE_LT32(s23,N0);
    b45=AE_LT32(s45,N0);
    b67=AE_LT32(s67,N0);
    pX0=(const xtfloatx4 *)x;
    pY0=(const xtfloatx4 *)y;
    ax0=AE_LA128_PP(pX0);
    ay0=AE_LA128_PP(pY0);
    AE_LASX2X2_IP(x0,x1,ax0,pX0);
    AE_LASX2X2_IP(x2,x3,ax0,pX0);
    AE_LASX2X2_IP(y0,y1,ay0,pY0);
    AE_LASX2X2_IP(y2,y3,ay0,pY0);
    {
        ae_int32x2 t;
        t=AE_MOVINT32X2_FROMXTFLOATX2(x0);  AE_MOVF32X2(t,0,b01); x0=AE_MOVXTFLOATX2_FROMINT32X2(t);   
        t=AE_MOVINT32X2_FROMXTFLOATX2(x1);  AE_MOVF32X2(t,0,b23); x1=AE_MOVXTFLOATX2_FROMINT32X2(t); 
        t=AE_MOVINT32X2_FROMXTFLOATX2(x2);  AE_MOVF32X2(t,0,b45); x2=AE_MOVXTFLOATX2_FROMINT32X2(t); 
        t=AE_MOVINT32X2_FROMXTFLOATX2(x3);  AE_MOVF32X2(t,0,b67); x3=AE_MOVXTFLOATX2_FROMINT32X2(t); 
        t=AE_MOVINT32X2_FROMXTFLOATX2(y0);  AE_MOVF32X2(t,0,b01); y0=AE_MOVXTFLOATX2_FROMINT32X2(t); 
        t=AE_MOVINT32X2_FROMXTFLOATX2(y1);  AE_MOVF32X2(t,0,b23); y1=AE_MOVXTFLOATX2_FROMINT32X2(t); 
        t=AE_MOVINT32X2_FROMXTFLOATX2(y2);  AE_MOVF32X2(t,0,b45); y2=AE_MOVXTFLOATX2_FROMINT32X2(t); 
        t=AE_MOVINT32X2_FROMXTFLOATX2(y3);  AE_MOVF32X2(t,0,b67); y3=AE_MOVXTFLOATX2_FROMINT32X2(t); 
    }
    MUL_SX2X2(a0,a1,x0,x1,y0,y1);
    MUL_SX2X2(a2,a3,x2,x3,y2,y3);
    x+=N0;
    y+=N0;
    N-=N0;
    pX0=(const xtfloatx4 *)(x);
    pY0=(const xtfloatx4 *)(y);
    ax0=AE_LA128_PP(pX0);
    ay0=AE_LA128_PP(pY0);
    a4=a5=a6=a7=CONST_SX2(0);
    if(N&8)
    {
        AE_LASX2X2_IP(x0,x1,ax0,pX0);
        AE_LASX2X2_IP(x2,x3,ax0,pX0);
        AE_LASX2X2_IP(y0,y1,ay0,pY0);
        AE_LASX2X2_IP(y2,y3,ay0,pY0);
        MUL_SX2X2(a4,a5,x0,x1,y0,y1);
        MUL_SX2X2(a6,a7,x2,x3,y2,y3);
    }
    for (n=0; n<(N>>4); n++) 
    {
        AE_LASX2X2_IP(x0,x1,ax0,pX0);
        AE_LASX2X2_IP(x2,x3,ax0,pX0);
        AE_LASX2X2_IP(y0,y1,ay0,pY0);
        AE_LASX2X2_IP(y2,y3,ay0,pY0);
        MADD_SX2X2(a0,a1,x0,x1,y0,y1);
        MADD_SX2X2(a2,a3,x2,x3,y2,y3);
        AE_LASX2X2_IP(x0,x1,ax0,pX0);
        AE_LASX2X2_IP(x2,x3,ax0,pX0);
        AE_LASX2X2_IP(y0,y1,ay0,pY0);
        AE_LASX2X2_IP(y2,y3,ay0,pY0);
        MADD_SX2X2(a4,a5,x0,x1,y0,y1);
        MADD_SX2X2(a6,a7,x2,x3,y2,y3);
    }
    a0=XT_ADD_SX2(a0,a4);
    a1=XT_ADD_SX2(a1,a5);
    a2=XT_ADD_SX2(a2,a6);
    a3=XT_ADD_SX2(a3,a7);
    a0=XT_ADD_SX2(a0,a1);
    a2=XT_ADD_SX2(a2,a3);
    a0=XT_ADD_SX2(a0,a2);
    a0=ADD_HL_LH_S(a0,a0);
        return (float32_t)XT_LOW_S(a0);
}
#elif (HAVE_FPU)
float32_t vec_dotf   (const float32_t * restrict x,const float32_t * restrict y,int N)
{
  xtfloat acc0, acc1,x0,y0;
  int n;
  const xtfloat  * restrict pX = (const xtfloat  *)x;
  const xtfloat  * restrict pY = (const xtfloat  *)y;
  if (N <= 0) return 0.f;
  acc0 = acc1 = XT_CONST_S(0);
  for (n = 0; n<(N&~1); n+=2)
  {
    XT_LSIP(x0, pX, sizeof(xtfloat));
    XT_LSIP(y0, pY, sizeof(xtfloat));
    XT_MADD_S(acc0,x0,y0);
    XT_LSIP(x0, pX, sizeof(xtfloat));
    XT_LSIP(y0, pY, sizeof(xtfloat));
    XT_MADD_S(acc1,x0,y0);
  }
  if (N&1)
  {
    XT_LSIP(x0, pX, sizeof(xtfloat));
    XT_LSIP(y0, pY, sizeof(xtfloat));
    XT_MADD_S(acc0,x0,y0);
  }
  return XT_ADD_S(acc0 , acc1 );
}
#endif
