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
  NatureDSP Signal Processing Library. Fitting and Interpolation Routines
    Polynomial approximation
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/
#include "NatureDSP_types.h"
/* Library API */
#include "NatureDSP_Signal_fit.h"
#include "common.h"
#include "common_fpu.h"

#if (HAVE_VFPU==0 && HAVE_FPU==0)
DISCARD_FUN(void,vec_poly4f,(float32_t * restrict z, const float32_t * restrict x, const float32_t * restrict c, int N ))
#elif (HAVE_VFPU)

/*-------------------------------------------------------------------------
  Polynomial approximation
  Functions calculate polynomial approximation for all values from given 
  vector. Fixed point functions take polynomial coefficients in Q31 precision. 
  NOTE:
  approximation is calculated like Taylor series that is why overflow may 
  potentially occur if cumulative sum of coefficients given from the last 
  to the first coefficient is bigger that 1. To avoid this negative effect,
  all the coefficients may be scaled down and result will be shifted left 
  after all intermediate computations.

  Precision: 
  32x32  32-bit inputs, 32-bit coefficients, 32-bit output.
  f      floating point

  Input:
  x[N]    input data, Q31 or floating point
  N       length of vector
  lsh     additional left shift for result
  c[M+1]  coefficients (M=4 coefficients for vec_poly4_xxx 
          and M=8 coefficients for vec_poly8_xxx), Q31 or floating point
  Output:			
  z[N]    result, Q31 or floating point

  Restriction:
  x,c,z should not overlap
  lsh   should be in range 0...31

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  x aligned on 16-byte boundary
  N   - multiple of 2
-------------------------------------------------------------------------*/
#define sz_f32    (int)sizeof(float32_t)
void vec_poly4f      (float32_t * restrict z, const float32_t * restrict x, const float32_t * restrict c, int N )
{
  /*
  * float32_t X,t;
  *  X=(x[n]);
  *  t=X*c4+c3;
  *  t=X*t +c2;
  *  t=X*t +c1;
  *  t=X*t +c0;
  *  z[n]=(t);
  */
    int n,Nhalf;

    xtfloatx2 c0f, c1f, c2f, c3f, c4f;
    ae_int32x2 d0,d1,d2,d3,d4;

    const xtfloatx4* restrict pX = (const xtfloatx4 *)x;
    const xtfloatx4* restrict pY;
    const ae_int32 *  restrict pc = (const ae_int32 *)c;
          xtfloatx4* restrict pZ = (xtfloatx4 *)z;
          xtfloatx4* restrict pW;
          ae_int32x4* restrict pT;
    ae_valignx2 aX,aY,aZ,aW;
    NASSERT(x);
    NASSERT(c);
    NASSERT(z);
    if (N <= 0) return;
    pc+=3;
    d4=AE_L32_I ( pc, 1 * sz_f32);   c4f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d4);
    if (N>7)
    {
        Nhalf=(N&~7)>>1;
        pY=(const xtfloatx4*)(x+Nhalf);
        pW=(      xtfloatx4*)(z+Nhalf);
        aX= AE_LA128_PP(pX);
        aY= AE_LA128_PP(pY);
        aW=aZ= AE_ZALIGN128();
        for (n = 0; n < (N>>3); n++)
        {
            xtfloatx2 x0,x1,x2,x3;
            xtfloatx2  t01,t00,t11,t10,t21,t20,t31,t30;
            xtfloatx2  t03,t02,t13,t12,t23,t22,t33,t32;

            AE_L32_IP(d3, pc, -1 * sz_f32);   c3f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d3);
            AE_L32_IP(d2, pc, -1 * sz_f32);   c2f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d2);
            AE_L32_IP(d1, pc, -1 * sz_f32);   c1f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d1);
            AE_L32_XP(d0, pc,  3 * sz_f32);   c0f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d0);
    
            AE_LASX2X2_IP(x0,x1,aX,pX);
            AE_LASX2X2_IP(x2,x3,aY,pY);

            MOV_SX2X2(t33,t32,c3f,c3f); MOV_SX2X2(t31,t20,c3f,c2f);
            MOV_SX2X2(t23,t22,c2f,c2f); MOV_SX2X2(t21,t30,c2f,c3f);
            MOV_SX2X2(t13,t12,c1f,c1f); MOV_SX2X2(t11,t00,c1f,c0f);
            MOV_SX2X2(t03,t02,c0f,c0f); MOV_SX2X2(t01,t10,c0f,c1f);
            MADDQ_S(t30,t31, x0,x1, c4f);
            MADDQ_S(t32,t33, x2,x3, c4f);
            MADD_SX2X2(t20,t21, x0,x1, t30,t31);
            MADD_SX2X2(t22,t23, x2,x3, t32,t33);
            MADD_SX2X2(t10,t11, x0,x1, t20,t21);
            MADD_SX2X2(t12,t13, x2,x3, t22,t23);
            MADD_SX2X2(t00,t01, x0,x1, t10,t11);
            MADD_SX2X2(t02,t03, x2,x3, t12,t13);
            AE_SASX2X2_IP(t00,t01, aZ,pZ);
            AE_SASX2X2_IP(t02,t03, aW,pW);
        }
        AE_SA128POS_FP(aZ, pZ);
        AE_SA128POS_FP(aW, pW);
        pX=(const xtfloatx4*)(x+(N&~7));
        pZ=(      xtfloatx4*)(z+(N&~7));
        N&=7;
    }
    if (N)
    {
        xtfloatx2  t01,t00,t11,t10,t21,t20,t31,t30;
        xtfloatx2  t03,t02,t13,t12,t23,t22,t33,t32;
        xtfloatx2 x0,x1,x2,x3;
        float32_t ALIGN(16) buf[8];
        pT=(ae_int32x4*)buf; AE_S32X2X2_I(0,0,pT,sizeof(ae_int32x4));
        __Pragma("loop_count min=1,max=7")
        for (n=0; n<N; n++) 
        {
            ae_int32x2 x0;
            AE_L32_IP(x0,castxcc(ae_int32,pX),sizeof(ae_int32));
            AE_S32_L_IP(x0,castxcc(ae_int32,pT),sizeof(ae_int32));
        }
        pT=(ae_int32x4*)buf;
        AE_LSX2X2_I (x0,x1,(const xtfloatx4*)pT,0*sizeof(ae_int32x4)); 
        AE_LSX2X2_I (x2,x3,(const xtfloatx4*)pT,1*sizeof(ae_int32x4)); 

        AE_L32_IP(d3, pc, -1 * sz_f32);   c3f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d3);
        AE_L32_IP(d2, pc, -1 * sz_f32);   c2f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d2);
        AE_L32_IP(d1, pc, -1 * sz_f32);   c1f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d1);
        AE_L32_XP(d0, pc,  3 * sz_f32);   c0f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d0);
     
        MOV_SX2X2(t33,t32,c3f,c3f); MOV_SX2X2(t31,t20,c3f,c2f);
        MOV_SX2X2(t23,t22,c2f,c2f); MOV_SX2X2(t21,t30,c2f,c3f);
        MOV_SX2X2(t13,t12,c1f,c1f); MOV_SX2X2(t11,t00,c1f,c0f);
        MOV_SX2X2(t03,t02,c0f,c0f); MOV_SX2X2(t01,t10,c0f,c1f);
        MADDQ_S(t30,t31, x0,x1, c4f);
        MADDQ_S(t32,t33, x2,x3, c4f);
        MADD_SX2X2(t20,t21, x0,x1, t30,t31);
        MADD_SX2X2(t22,t23, x2,x3, t32,t33);
        MADD_SX2X2(t10,t11, x0,x1, t20,t21);
        MADD_SX2X2(t12,t13, x2,x3, t22,t23);
        MADD_SX2X2(t00,t01, x0,x1, t10,t11);
        MADD_SX2X2(t02,t03, x2,x3, t12,t13);

        AE_SSX2X2_I (t00,t01, (xtfloatx4*)pT,0*sizeof(ae_int32x4)); 
        AE_SSX2X2_I (t02,t03, (xtfloatx4*)pT,1*sizeof(ae_int32x4)); 
        __Pragma("loop_count min=1,max=7")
        for (n=0; n<N; n++) 
        {
            ae_int32x2 x0;
            AE_L32_IP(x0,castxcc(ae_int32,pT),sizeof(ae_int32));
            AE_S32_L_IP(x0,castxcc(ae_int32,pZ),sizeof(ae_int32));
        }
    }
} /* vec_poly4f() */
#else
void vec_poly4f      (float32_t * restrict z, const float32_t * restrict x, const float32_t * restrict c, int N )
{
  int n;
  const xtfloat * restrict pX=(const xtfloat *)x;
        xtfloat * restrict pZ=(      xtfloat *)z;
  xtfloat c0,c1,c2,c3,c4;
  c0=XT_LSI((const xtfloat*)c,sizeof(xtfloat)*0);
  c1=XT_LSI((const xtfloat*)c,sizeof(xtfloat)*1);
  c2=XT_LSI((const xtfloat*)c,sizeof(xtfloat)*2);
  c3=XT_LSI((const xtfloat*)c,sizeof(xtfloat)*3);
  c4=XT_LSI((const xtfloat*)c,sizeof(xtfloat)*4);
  /* NOTE: both orders of summation are possible */
  for (n=0; n<N; n++)
  {
    xtfloat t,z0,z1,x2,xn;
    XT_LSIP(xn,pX,sizeof(xtfloat));
    x2=XT_MUL_S(xn,xn);
    z0=c2; XT_MADD_S(z0,x2,c4);
    z1=c1; XT_MADD_S(z1,x2,c3);
    t=c0; XT_MADD_S(t,x2,z0);z0=t;
    t=z0; XT_MADD_S(t,xn,z1);z1=t;
    XT_SSIP(z1,pZ,sizeof(xtfloat));
  }
}
#endif
