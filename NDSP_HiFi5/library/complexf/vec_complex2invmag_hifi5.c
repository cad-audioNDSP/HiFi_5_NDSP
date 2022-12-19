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
#include "NatureDSP_Signal_complex.h"
#include "common.h"
#include "common_fpu.h"
#include "inff_tbl.h"
// optimized for HiFi5

#if (HAVE_VFPU)
inline_ void rsqrt_sx2x2(xtfloatx2* pd0,xtfloatx2* pd1,xtfloatx2* ps0,xtfloatx2* ps1)
{
    xtfloatx2 s0,s1,t10,t11,t20,t21,t30,t40,t41,t50,t51,t60,t61,t70,t71,t80,t81;
    s0=*ps0; s1=*ps1;
    t10=XT_RSQRT0_SX2(s0); t11=XT_RSQRT0_SX2(s1);
    MUL_SX2X2 (t20,t21,s0,s1,t10,t11);
    t30=XT_CONST_S(3);
    MULQ_S(t40,t41,t10,t11,t30);
    CONST_SX2X2(t50,t51, 1);
    MSUB_SX2X2(t50,t51,t20,t21,t10,t11);
    MADD_SX2X2(t10,t11,t40,t41,t50,t51);
    MUL_SX2X2 (t60,t61,s0,s1,t10,t11);
    MULQ_S(t70,t71,t10,t11,t30);
    CONST_SX2X2(t80,t81, 1);
    MSUB_SX2X2(t80,t81,t60,t61,t10,t11);
    MADD_SX2X2(t10,t11,t70,t71,t80,t81);
    *pd0=t10; *pd1=t11;
}

#define RSQRT_SX2X2(vd0,vd1,vs0,vs1) rsqrt_sx2x2(&vd0,&vd1,&vs0,&vs1)
#endif

#if (HAVE_VFPU==0 && HAVE_FPU==0)
DISCARD_FUN(void,vec_complex2invmag,(float32_t  * restrict y, const complex_float  * restrict x, int N))
#elif (HAVE_VFPU)
/*===========================================================================
  Vector matematics:
  vec_complex2invmag     complex magnitude (reciprocal)
===========================================================================*/
/*-------------------------------------------------------------------------
  Complex magnitude
  Routines compute complex magnitude or its reciprocal

  Precision: 
  f     single precision floating point

  Input:
  x[N]  input complex data
  N     length of vector
  Output:
  y[N]  output data

  Restriction:
  none
-------------------------------------------------------------------------*/
#define sz_i32  (int)sizeof(int32_t)
#define sz_f32  (int)sizeof(float32_t)
#include <math.h>
void       vec_complex2invmag (float32_t  * restrict y, const complex_float  * restrict x, int N)
{
    ae_valignx2 aX;
    ae_valign aY;
    int n;
    static const union { uint32_t u; float32_t f; } 
                 xmin={0x20000000}, // sqrt(realmin)
                 xmax={0x5f000000}, // sqrt(realmax)/2
                 k   ={0x6a800000}, // direct and inverse scale factors 
                 kinv={0x14800000}; //
    if (N<=0) return;
    aY=AE_ZALIGN64();
    aX=AE_LA128_PP(x);
    for ( n=0; n<(N>>2); n++)
    {
        xtbool2 bsmall,bbig;
        xtfloatx2 reim0,reim1,rere,imim,maxreim,r0,r1,t0,t1;
        xtfloatx2 s0,s1;
        AE_LASX2X2_IP(reim0,reim1,aX,castxcc(xtfloatx4,x));
        ABS_SX2X2(reim0,reim1,reim0,reim1);
        rere=AE_SEL32_HH_SX2(reim0,reim1);
        imim=AE_SEL32_LL_SX2(reim0,reim1);
        maxreim=MAXNUM_SX2(rere,imim);
        bsmall=XT_OLE_SX2(maxreim,xmin.f);
        bbig  =XT_OLE_SX2(xmax.f,maxreim);
        CONST_SX2X2(s0,s1,1);
        XT_MOVT_SX2(s0,k.f,bsmall);
        XT_MOVT_SX2(s0,kinv.f,bbig);
        MUL_SX2X2(rere,imim,rere,imim,s0,s0);
        MUL_SX2X2(r0,t0,rere,imim,rere,imim);

        AE_LASX2X2_IP(reim0,reim1,aX,castxcc(xtfloatx4,x));
        ABS_SX2X2(reim0,reim1,reim0,reim1);
        rere=AE_SEL32_HH_SX2(reim0,reim1);
        imim=AE_SEL32_LL_SX2(reim0,reim1);
        maxreim=MAXNUM_SX2(rere,imim);
        bsmall=XT_OLE_SX2(maxreim,xmin.f);
        bbig  =XT_OLE_SX2(xmax.f,maxreim);
        XT_MOVT_SX2(s1,k.f,bsmall);
        XT_MOVT_SX2(s1,kinv.f,bbig);
        MUL_SX2X2(rere,imim,rere,imim,s1,s1);
        MUL_SX2X2(r1,t1,rere,imim,rere,imim);

        ADD_SX2X2(r0,r1,r0,r1,t0,t1);
        {
            xtbool2 bzero0,bzero1,binf0,binf1;
            bzero0=XT_OEQ_SX2(r0,0.f);
            bzero1=XT_OEQ_SX2(r1,0.f);
            binf0 =XT_OEQ_SX2(r0,plusInff.f);
            binf1 =XT_OEQ_SX2(r1,plusInff.f);
            RSQRT_SX2X2(r0,r1,r0,r1);
            XT_MOVT_SX2(r0,plusInff.f,bzero0);
            XT_MOVT_SX2(r1,plusInff.f,bzero1);
            XT_MOVT_SX2(r0,0.f,binf0);
            XT_MOVT_SX2(r1,0.f,binf1);
        }
        MUL_SX2X2(r0,r1,r0,r1,s0,s1);

        XT_SASX2IP(r0,aY,castxcc(xtfloatx2,y));
        XT_SASX2IP(r1,aY,castxcc(xtfloatx2,y));
    }
    AE_SA64POS_FP(aY,castxcc(xtfloatx2,y));
    // tail
    for (n=0; n<(N&3); n++)
    {
        xtbool2 bsmall,bbig;
        xtbool2 bzero, binf;
        xtfloatx2 reim;
        xtfloatx2 r, maxreim;
        xtfloatx2 s;
        XT_LSX2IP(reim,castxcc(xtfloatx2,x),sizeof(xtfloatx2));
        reim=XT_ABS_SX2(reim);
        maxreim=RMAXNUM_S(reim);
        bsmall=XT_OLE_SX2(maxreim,xmin.f);
        bbig  =XT_OLE_SX2(xmax.f,maxreim);
        s=XT_CONST_S(1); 
        XT_MOVT_SX2(s,k.f,bsmall);
        XT_MOVT_SX2(s,kinv.f,bbig);
        reim=XT_MUL_SX2(reim,s);
        reim=XT_MUL_SX2(reim,reim);
        r=XT_RADD_SX2(reim);
        bzero=XT_OEQ_S(r,0.f);
        binf =XT_OEQ_S(r,plusInff.f);
        r=XT_RSQRT_S(r);
        XT_MOVT_SX2(r,plusInff.f,bzero);
        XT_MOVT_SX2(r,0.f,binf);
        r=XT_MUL_S(r,s);
        XT_SSIP(r,castxcc(xtfloat,y),sizeof(float32_t));
    }
}
#else
// for scalar FPU
void       vec_complex2invmag (float32_t  * restrict y, const complex_float  * restrict x, int N)
{
    int n;
    const xtfloat* restrict pX=(const xtfloat*)x;
          xtfloat* restrict pY=(      xtfloat*)y;
    if (N<=0) return;
    for (n=0; n<N; n++)
    {
        xtfloat x0, x1, y0,coef;
        xtfloat xre, xim;
        int t0, t1, nsa;
        int e0;
        int nsa0;
        xtbool b0, b1;

        XT_LSIP(xre,pX,1*sizeof(xtfloat));
        XT_LSIP(xim,pX,1*sizeof(xtfloat));
        xre = XT_ABS_S(xre);
        xim = XT_ABS_S(xim);

        t0 = XT_RFR(xre);
        t1 = XT_RFR(xim);

        nsa = XT_MAX(t0, t1);
        nsa = ((uint32_t)nsa)>>23;
        nsa = nsa-127;
        nsa = XT_MIN(nsa, 127);
        e0 = (127-nsa);
        nsa0 = e0<<23;
        XT_MOVEQZ(nsa0,0x00400000,e0);
        coef = XT_WFR(nsa0);

        xre = XT_MUL_S(xre, coef);
        xim = XT_MUL_S(xim, coef);

        x0 = XT_MUL_S(xre, xre);
        x1 = XT_MUL_S(xim, xim);

        x0 = XT_ADD_S(x0, x1);

        y0 = XT_RSQRT_S(x0);
        b0 = XT_OEQ_S(x0, XT_CONST_S(0));
        b1 = XT_OEQ_S(x0, plusInff.f);
        XT_MOVT_S(y0,plusInff.f,b0);
        XT_MOVT_S(y0,XT_CONST_S(0),b1);
        y0 = XT_MUL_S(y0, coef);
        XT_SSIP(y0,pY,sizeof(xtfloat));
    }
}
#endif
