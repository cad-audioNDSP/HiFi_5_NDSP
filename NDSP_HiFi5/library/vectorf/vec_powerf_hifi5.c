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
    Power of a Vector
    C code optimized for HiFi4
  IntegrIT, 2006-2019
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_vector.h"
#include "NatureDSP_types.h"
#include "common.h"
#include "common_fpu.h"

#if (HAVE_VFPU==0 && HAVE_FPU==0)
DISCARD_FUN(float32_t ,vec_powerf,( const float32_t * restrict x,int N))
#elif (HAVE_VFPU)

/*===========================================================================
  Vector matematics:
  vec_power            Power of a Vector
===========================================================================*/
/*-------------------------------------------------------------------------
  Power of a Vector
  These routines compute power of vector with scaling output result by rsh 
  bits. Fixed point rountines make accumulation in the 64-bit wide 
  accumulator and output may scaled down with saturation by rsh bits. 
  So, if representation of x input is Qx, result will be represented in 
  Q(2x-rsh) format.
  Two versions of routines are available: regular versions (vec_power32x32, 
  vec_power16x16, vec_powerf) work with arbitrary arguments, faster versions 
  (vec_power32x32_fast, vec_power16x16_fast) apply some restrictions.

  Precision: 
  32x32 32x32-bit data, 64-bit output
  16x16 16x16-bit data, 64-bit output
  f     single precision floating point

  Input:
  x[N]  input data, Q31, Q15 or floating point
  rsh   right shift of result
  N     length of vector
  Returns: 
  Sum of squares of a vector, Q(2x-rsh)

  Restrictions:
  for vec_power32x32(): rsh in range 31...62
  for vec_power16x16(): rsh in range 0...31
  For regular versions (vec_power32x32, vec_power16x16, vec_powerf):
  none
  For faster versions (vec_power32x32_fast, vec_power16x16_fast ):
  x - aligned on 16-byte boundary
  N - multiple of 4
-------------------------------------------------------------------------*/
#define sz_f32    (int)sizeof(float32_t)
float32_t   vec_powerf     ( const float32_t * restrict x,int N)
#if 0
{
  int n;

  xtfloatx2 vacc0, vacc1, vacc2, vacc3;
  xtfloatx2 vacc4, vacc5, vacc6, vacc7;
  xtfloat xf, zf, acc;
  xtfloatx2 x0, x1;
  const xtfloatx2 * restrict px0 = (const xtfloatx2 *)x;

  NASSERT(x);
  if (N <= 0) return 0;

  vacc0 = XT_MOV_SX2(0.f);
  vacc1 = XT_MOV_SX2(0.f);
  vacc2 = XT_MOV_SX2(0.f);
  vacc3 = XT_MOV_SX2(0.f);
  vacc4 = XT_MOV_SX2(0.f);
  vacc5 = XT_MOV_SX2(0.f);
  vacc6 = XT_MOV_SX2(0.f);
  vacc7 = XT_MOV_SX2(0.f);
  zf = XT_MOV_S(0.f);

  if ((((uintptr_t)(x)) & 7))
  {
    ae_int32x2 tmp;
    xf = XT_LSI((const xtfloat *)px0, 0);
    AE_L32_IP(tmp, castxcc(ae_int32, px0), sz_f32);
    zf = XT_MUL_S(xf, xf);

    N--;
  }
  for (n = 0; n<(N>>4); n ++)
  {
    XT_LSX2IP(x0, px0, 2 * sz_f32);
    XT_LSX2IP(x1, px0, 2 * sz_f32);
    XT_MADD_SX2(vacc0, x0, x0);
    XT_MADD_SX2(vacc1, x1, x1);
    XT_LSX2IP(x0, px0, 2 * sz_f32);
    XT_LSX2IP(x1, px0, 2 * sz_f32);
    XT_MADD_SX2(vacc2, x0, x0);
    XT_MADD_SX2(vacc3, x1, x1);
    XT_LSX2IP(x0, px0, 2 * sz_f32);
    XT_LSX2IP(x1, px0, 2 * sz_f32);
    XT_MADD_SX2(vacc4, x0, x0);
    XT_MADD_SX2(vacc5, x1, x1);
    XT_LSX2IP(x0, px0, 2 * sz_f32);
    XT_LSX2IP(x1, px0, 2 * sz_f32);
    XT_MADD_SX2(vacc6, x0, x0);
    XT_MADD_SX2(vacc7, x1, x1); 
  }
  if (N & 8)
  {
    XT_LSX2IP(x0, px0, 2 * sz_f32);
    XT_MADD_SX2(vacc0, x0, x0);
    XT_LSX2IP(x0, px0, 2 * sz_f32);
    XT_MADD_SX2(vacc0, x0, x0);
    XT_LSX2IP(x0, px0, 2 * sz_f32);
    XT_MADD_SX2(vacc0, x0, x0);
    XT_LSX2IP(x0, px0, 2 * sz_f32);
    XT_MADD_SX2(vacc0, x0, x0);
  }
  if (N & 4)
  {
    XT_LSX2IP(x0, px0, 2 * sz_f32);
    XT_MADD_SX2(vacc0, x0, x0);
    XT_LSX2IP(x0, px0, 2 * sz_f32);
    XT_MADD_SX2(vacc0, x0, x0);
  }
  if (N & 2)
  {
    XT_LSX2IP(x0, px0, 2 * sz_f32);
    XT_MADD_SX2(vacc0, x0, x0);
  }
  if (N & 1)
  {
    xf = XT_LSI((const xtfloat *)px0, 0);
    XT_MADD_S(zf, xf, xf);
  }
  vacc0 = XT_ADD_SX2(vacc0, vacc1);
  vacc2 = XT_ADD_SX2(vacc3, vacc2);
  vacc4 = XT_ADD_SX2(vacc4, vacc5);
  vacc6 = XT_ADD_SX2(vacc6, vacc7);
  vacc0 = XT_ADD_SX2(vacc0, vacc2);
  vacc4 = XT_ADD_SX2(vacc4, vacc6);
  vacc0 = XT_ADD_SX2(vacc0, vacc4);

  acc = XT_RADD_SX2(vacc0);
  acc = XT_ADD_S(acc, zf);

  return acc;
} /* vec_powerf() */
#else
{
    static const int32_t ALIGN(16) seq[]={0,1,2,3,4,5,6,7};
    xtbool2 b01,b23,b45,b67;
    ae_int32x2 s01,s23,s45,s67;
    const xtfloatx4* restrict pX0;
    const xtfloatx4* restrict pX1;
    xtfloatx2 x0,x1,x2,x3;
    xtfloatx2 a0,a1,a2,a3;
    xtfloatx2 a4,a5,a6,a7;
    xtfloatx2 a8,a9,aa,ab;
    xtfloatx2 ac,ad,ae,af;
    ae_valignx2 ax0,ax1;
    int n,N0;
    if(N<=0) return XT_CONST_S(0);
    if (N<=7)
    {
        xtfloat a0;
        a0=XT_CONST_S(0);
        __Pragma("no_unroll")
        for (n=0; n<N; n++) 
        {
            xtfloat x0;
            XT_LSIP(x0,castxcc(xtfloat,x),sizeof(float32_t));
            XT_MADD_S(a0,x0,x0);
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
    ax0=AE_LA128_PP(pX0);
    AE_LASX2X2_IP(x0,x1,ax0,pX0);
    AE_LASX2X2_IP(x2,x3,ax0,pX0);
    {
        ae_int32x2 t;
        t=AE_MOVINT32X2_FROMXTFLOATX2(x0);  AE_MOVF32X2(t,0,b01); x0=AE_MOVXTFLOATX2_FROMINT32X2(t);   
        t=AE_MOVINT32X2_FROMXTFLOATX2(x1);  AE_MOVF32X2(t,0,b23); x1=AE_MOVXTFLOATX2_FROMINT32X2(t); 
        t=AE_MOVINT32X2_FROMXTFLOATX2(x2);  AE_MOVF32X2(t,0,b45); x2=AE_MOVXTFLOATX2_FROMINT32X2(t); 
        t=AE_MOVINT32X2_FROMXTFLOATX2(x3);  AE_MOVF32X2(t,0,b67); x3=AE_MOVXTFLOATX2_FROMINT32X2(t); 
    }
    MUL_SX2X2(a0,a1,x0,x1,x0,x1);
    MUL_SX2X2(a2,a3,x2,x3,x2,x3);
    x+=N0;
    N-=N0;
    pX0=(const xtfloatx4 *)(x);
    pX1=(const xtfloatx4 *)(x+(N>>1));
    ax0=AE_LA128_PP(pX0);
    ax1=AE_LA128_PP(pX1);
    CONST_SX2X2(a4,a5,0);CONST_SX2X2(a6,a7,0);
    CONST_SX2X2(a8,a9,0);CONST_SX2X2(aa,ab,0);
    CONST_SX2X2(ac,ad,0);CONST_SX2X2(ae,af,0);
    if (N&8)
    {
        AE_LASX2X2_IP(x0,x1,ax0,pX0);
        AE_LASX2X2_IP(x2,x3,ax1,pX1);
        MADD_SX2X2(a0,a1,x0,x1,x0,x1);
        MADD_SX2X2(a2,a3,x2,x3,x2,x3);
    }
    if (N&16)
    {
        AE_LASX2X2_IP(x0,x1,ax0,pX0);
        AE_LASX2X2_IP(x2,x3,ax1,pX1);
        MADD_SX2X2(a0,a1,x0,x1,x0,x1);
        MADD_SX2X2(a2,a3,x2,x3,x2,x3);
        AE_LASX2X2_IP(x0,x1,ax0,pX0);
        AE_LASX2X2_IP(x2,x3,ax1,pX1);
        MADD_SX2X2(a4,a5,x0,x1,x0,x1);
        MADD_SX2X2(a6,a7,x2,x3,x2,x3);
    }
    for (n=0; n<(N>>5); n++) 
    {
        AE_LASX2X2_IP(x0,x1,ax0,pX0);
        AE_LASX2X2_IP(x2,x3,ax1,pX1);
        MADD_SX2X2(a0,a1,x0,x1,x0,x1);
        MADD_SX2X2(a2,a3,x2,x3,x2,x3);
        AE_LASX2X2_IP(x0,x1,ax0,pX0);
        AE_LASX2X2_IP(x2,x3,ax1,pX1);
        MADD_SX2X2(a4,a5,x0,x1,x0,x1);
        MADD_SX2X2(a6,a7,x2,x3,x2,x3);
        AE_LASX2X2_IP(x0,x1,ax0,pX0);
        AE_LASX2X2_IP(x2,x3,ax1,pX1);
        MADD_SX2X2(a8,a9,x0,x1,x0,x1);
        MADD_SX2X2(aa,ab,x2,x3,x2,x3);
        AE_LASX2X2_IP(x0,x1,ax0,pX0);
        AE_LASX2X2_IP(x2,x3,ax1,pX1);
        MADD_SX2X2(ac,ad,x0,x1,x0,x1);
        MADD_SX2X2(ae,af,x2,x3,x2,x3);
    }
    a8=XT_ADD_SX2(a8,ac); a9=XT_ADD_SX2(a9,ad);
    aa=XT_ADD_SX2(aa,ae); ab=XT_ADD_SX2(ab,af);
    a0=XT_ADD_SX2(a0,a4); a1=XT_ADD_SX2(a1,a5);
    a2=XT_ADD_SX2(a2,a6); a3=XT_ADD_SX2(a3,a7);
    a0=XT_ADD_SX2(a0,a8); a1=XT_ADD_SX2(a1,a9);
    a2=XT_ADD_SX2(a2,aa); a3=XT_ADD_SX2(a3,ab);
    a0=XT_ADD_SX2(a0,a1);
    a2=XT_ADD_SX2(a2,a3);
    a0=XT_ADD_SX2(a0,a2);
    a0=ADD_HL_LH_S(a0,a0);
    return (float32_t)XT_LOW_S(a0);
}
#endif
#else
float32_t   vec_powerf     ( const float32_t * restrict x,int N)
{
  xtfloat acc0, acc1,x0;
  int n;
  const xtfloat  * restrict pX = (const xtfloat  *)x;
  if (N <= 0) return 0.f;
  acc0 = acc1 = XT_CONST_S(0);
  for (n = 0; n<(N&~1); n+=2)
  {
    XT_LSIP(x0, pX, sizeof(xtfloat));
    XT_MADD_S(acc0,x0,x0);
    XT_LSIP(x0, pX, sizeof(xtfloat));
    XT_MADD_S(acc1,x0,x0);
  }
  if (N&1)
  {
    XT_LSIP(x0, pX, sizeof(xtfloat));
    XT_MADD_S(acc0,x0,x0);
  }
  return XT_ADD_S(acc0 , acc1 );
}
#endif
