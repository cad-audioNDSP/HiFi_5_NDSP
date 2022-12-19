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
DISCARD_FUN(void,vec_poly8f ,(float32_t * restrict z, const float32_t * restrict x, const float32_t * restrict c, int N ))
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
#define sz_f32  (int)sizeof(float32_t)
void vec_poly8f (float32_t * restrict z, const float32_t * restrict x, const float32_t * restrict c, int N )
{
    int n,Nhalf;
    xtfloatx2 c0f, c1f, c2f, c3f, c4f, c5f,c6f,c7f,c8f;
    ae_int32x2 d0,d1,d2,d3,d4,d5,d6,d7,d8;

    const xtfloatx4* restrict pX = (const xtfloatx4 *)x;
    const xtfloatx4* restrict pY;
    const ae_int32 *  restrict pc = (const ae_int32 *)c;
          xtfloatx4* restrict pZ = (xtfloatx4 *)z;
          ae_int32x4* restrict pT;
          xtfloatx4* restrict pW;
    ae_valignx2 aX,aY,aZ,aW;
    NASSERT(x);
    NASSERT(c);
    NASSERT(z);
    if (N <= 0) return;
    pc+=8;
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
            xtfloatx2 y0,y1,y2,y3;
            xtfloatx2 z0,z1,z2,z3;
            xtfloatx2  t01,t00,t11,t10,t21,t20,t31,t30,t41,t40;
            xtfloatx2  t03,t02,t13,t12,t23,t22,t33,t32,t43,t42;

            AE_L32_IP(d8, pc, -1 * sz_f32);   c8f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d8);
            AE_L32_XP(d7, pc, -1 *sz_f32);    c7f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d7);
            AE_L32_IP(d6, pc, -1 *sz_f32);    c6f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d6);
            AE_L32_XP(d5, pc, -1 *sz_f32);    c5f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d5);
            AE_L32_IP(d4, pc, -1 *sz_f32);    c4f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d4);
            AE_L32_XP(d3, pc, -1 *sz_f32);    c3f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d3);
            AE_L32_IP(d2, pc, -1 *sz_f32);    c2f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d2);
            AE_L32_XP(d1, pc, -1 *sz_f32);    c1f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d1);
            AE_L32_XP(d0, pc,  8 *sz_f32);    c0f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d0);
    
            AE_LASX2X2_IP(x0,x1,aX,pX);
            AE_LASX2X2_IP(x2,x3,aY,pY);
            MUL_SX2X2(y0,y1,x0,x1,x0,x1);
            MUL_SX2X2(y2,y3,x2,x3,x2,x3);
            MUL_SX2X2(z0,z1,y0,y1,y0,y1);
            MUL_SX2X2(z2,z3,y2,y3,y2,y3);
            MOV_SX2X2(t43,t42,c4f,c4f);MOV_SX2X2(t41,t30,c4f,c3f);
            MOV_SX2X2(t33,t32,c3f,c3f);MOV_SX2X2(t31,t40,c3f,c4f);
            MOV_SX2X2(t23,t22,c2f,c2f);MOV_SX2X2(t21,t10,c2f,c1f);
            MOV_SX2X2(t13,t12,c1f,c1f);MOV_SX2X2(t11,t20,c1f,c2f);
            MOV_SX2X2(t03,t02,c0f,c0f);MOV_SX2X2(t01,t00,c0f,c0f);
            MADDQ_S(t40,t41, z0,z1, c8f);            MADDQ_S(t42,t43, z2,z3, c8f);
            MADDQ_S(t30,t31, z0,z1, c7f);            MADDQ_S(t32,t33, z2,z3, c7f);
            MADDQ_S(t20,t21, z0,z1, c6f);            MADDQ_S(t22,t23, z2,z3, c6f);           
            MADDQ_S(t10,t11, z0,z1, c5f);            MADDQ_S(t12,t13, z2,z3, c5f);
            MADD_SX2X2(t00,t01, z0,z1, t40,t41);     MADD_SX2X2(t02,t03, z2,z3, t42,t43);
            MADD_SX2X2(t10,t11, y0,y1, t30,t31);     MADD_SX2X2(t12,t13, y2,y3, t32,t33);
            MADD_SX2X2(t00,t01, y0,y1, t20,t21);     MADD_SX2X2(t02,t03, y2,y3, t22,t23);
            MADD_SX2X2(t00,t01, x0,x1, t10,t11);     MADD_SX2X2(t02,t03, x2,x3, t12,t13);

            AE_SASX2X2_IP(t00,t01, aZ,pZ);
            AE_SASX2X2_IP(t02,t03, aW,pW);
        }
        AE_SA128POS_FP(aZ, pZ);
        AE_SA128POS_FP(aW, pW);
        pX=(const xtfloatx4*)(x+(N&~7));
        pZ=(      xtfloatx4*)(z+(N&~7));
        N&=7;
    }
#if 0
    if (N)
    {       // tail
        pc = (const ae_int32 *)c;

        for (n=0;n<N;n++)
        {
            xtfloat  x0,y0,z0;
            xtfloat  t00,t10,t20,t30,t40,t50,t60,t70,t80;

            AE_L32_IP(d0, pc, sz_f32);        c0f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d0);
            AE_L32_IP(d1, pc, sz_f32);        c1f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d1);
            AE_L32_IP(d2, pc, sz_f32);        c2f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d2);
            AE_L32_IP(d3, pc, sz_f32);        c3f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d3);
            AE_L32_IP(d4, pc, sz_f32);        c4f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d4);
            AE_L32_IP(d5, pc, sz_f32);        c5f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d5);
            AE_L32_IP(d6, pc, sz_f32);        c6f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d6);
            AE_L32_IP(d7, pc, sz_f32);        c7f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d7);
            AE_L32_IP(d8, pc, -8 * sz_f32);   c8f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d8);
    
            XT_LSIP(x0,castxcc(xtfloat,pX),sizeof(float32_t));

            y0=XT_MUL_S(x0,x0);
            z0=XT_MUL_S(y0,y0);
            t00=c0f;
            t10=c1f;
            t20=c2f;
            t30=c3f;
            t40=c4f;
            t50=c5f;
            t60=c6f;
            t70=c7f;
            t80=c8f;
            XT_MADD_S(t40,z0, t80); 
            XT_MADD_S(t30,z0, t70); 
            XT_MADD_S(t20,z0, t60);       
            XT_MADD_S(t10,z0, t50); 
            XT_MADD_S(t00,z0, t40); 
            XT_MADD_S(t10,y0, t30); 
            XT_MADD_S(t00,y0, t20); 
            XT_MADD_S(t00,x0, t10); 

            XT_SSIP(t00, castxcc(xtfloat,pZ), sizeof(float32_t));
        }
    }
#else
    if (N)
    {
        xtfloatx2 x0,x1,x2,x3;
        xtfloatx2 y0,y1,y2,y3;
        xtfloatx2 z0,z1,z2,z3;
        xtfloatx2  t01,t00,t11,t10,t21,t20,t31,t30,t41,t40;
        xtfloatx2  t03,t02,t13,t12,t23,t22,t33,t32,t43,t42;
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

        AE_L32_IP(d8, pc, -1 * sz_f32);   c8f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d8);
        AE_L32_XP(d7, pc, -1 *sz_f32);    c7f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d7);
        AE_L32_IP(d6, pc, -1 *sz_f32);    c6f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d6);
        AE_L32_XP(d5, pc, -1 *sz_f32);    c5f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d5);
        AE_L32_IP(d4, pc, -1 *sz_f32);    c4f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d4);
        AE_L32_XP(d3, pc, -1 *sz_f32);    c3f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d3);
        AE_L32_IP(d2, pc, -1 *sz_f32);    c2f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d2);
        AE_L32_XP(d1, pc, -1 *sz_f32);    c1f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d1);
        AE_L32_XP(d0, pc,  8 *sz_f32);    c0f =XT_AE_MOVXTFLOATX2_FROMINT32X2 (d0);
    
        MUL_SX2X2(y0,y1,x0,x1,x0,x1);
        MUL_SX2X2(y2,y3,x2,x3,x2,x3);
        MUL_SX2X2(z0,z1,y0,y1,y0,y1);
        MUL_SX2X2(z2,z3,y2,y3,y2,y3);
        MOV_SX2X2(t43,t42,c4f,c4f);MOV_SX2X2(t41,t30,c4f,c3f);
        MOV_SX2X2(t33,t32,c3f,c3f);MOV_SX2X2(t31,t40,c3f,c4f);
        MOV_SX2X2(t23,t22,c2f,c2f);MOV_SX2X2(t21,t10,c2f,c1f);
        MOV_SX2X2(t13,t12,c1f,c1f);MOV_SX2X2(t11,t20,c1f,c2f);
        MOV_SX2X2(t03,t02,c0f,c0f);MOV_SX2X2(t01,t00,c0f,c0f);
        MADDQ_S(t40,t41, z0,z1, c8f);            MADDQ_S(t42,t43, z2,z3, c8f);
        MADDQ_S(t30,t31, z0,z1, c7f);            MADDQ_S(t32,t33, z2,z3, c7f);
        MADDQ_S(t20,t21, z0,z1, c6f);            MADDQ_S(t22,t23, z2,z3, c6f);           
        MADDQ_S(t10,t11, z0,z1, c5f);            MADDQ_S(t12,t13, z2,z3, c5f);
        MADD_SX2X2(t00,t01, z0,z1, t40,t41);     MADD_SX2X2(t02,t03, z2,z3, t42,t43);
        MADD_SX2X2(t10,t11, y0,y1, t30,t31);     MADD_SX2X2(t12,t13, y2,y3, t32,t33);
        MADD_SX2X2(t00,t01, y0,y1, t20,t21);     MADD_SX2X2(t02,t03, y2,y3, t22,t23);
        MADD_SX2X2(t00,t01, x0,x1, t10,t11);     MADD_SX2X2(t02,t03, x2,x3, t12,t13);

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
#endif
} /* vec_poly8f() */
#else
void vec_poly8f (float32_t * restrict z, const float32_t * restrict x, const float32_t * restrict c, int N )
{
  int n;
  const xtfloat * restrict pX=(const xtfloat *)x;
        xtfloat * restrict pZ=(      xtfloat *)z;
  xtfloat c0,c1,c2,c3,c4,c5,c6,c7,c8;
  c0=XT_LSI((const xtfloat*)c,sizeof(xtfloat)*0);
  c1=XT_LSI((const xtfloat*)c,sizeof(xtfloat)*1);
  c2=XT_LSI((const xtfloat*)c,sizeof(xtfloat)*2);
  c3=XT_LSI((const xtfloat*)c,sizeof(xtfloat)*3);
  c4=XT_LSI((const xtfloat*)c,sizeof(xtfloat)*4);
  c5=XT_LSI((const xtfloat*)c,sizeof(xtfloat)*5);
  c6=XT_LSI((const xtfloat*)c,sizeof(xtfloat)*6);
  c7=XT_LSI((const xtfloat*)c,sizeof(xtfloat)*7);
  c8=XT_LSI((const xtfloat*)c,sizeof(xtfloat)*8);
  /* NOTE: both orders of summation are possible */
  for (n=0; n<N; n++)
  {
    xtfloat z0,z1,x2,xn,t;
    XT_LSIP(xn,pX,sizeof(xtfloat));
    x2=XT_MUL_S(xn,xn);
    z0=c6; XT_MADD_S(z0,x2,c8);
    z1=c5; XT_MADD_S(z1,x2,c7);
    t=c4; XT_MADD_S(t,x2,z0); z0=t;
    t=c3; XT_MADD_S(t,x2,z1); z1=t;
    t=c2; XT_MADD_S(t,x2,z0); z0=t;
    t=c1; XT_MADD_S(t,x2,z1); z1=t;
    t=c0; XT_MADD_S(t,x2,z0);z0=t;
    t=z0; XT_MADD_S(t,xn,z1);z1=t;
    XT_SSIP(z1,pZ,sizeof(xtfloat));
  }
}
#endif
