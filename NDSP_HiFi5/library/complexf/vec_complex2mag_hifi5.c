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
// optimized for HiFi5

#if HAVE_VFPU
inline_ void sqrt_sx2x2(xtfloatx2 *pvd0,xtfloatx2 *pvd1,const xtfloatx2 *pvs0,const xtfloatx2 *pvs1) 
{
        xtfloatx2 vs0=pvs0[0],vs1=pvs1[0];
        xtfloatx2 f20,f21,f30,f31,f40,f41,f00,f01,f60,f61,f70,f71,f10,f11; 
        f20=SQRT0_SX2(vs0); f21=SQRT0_SX2(vs1);
        MUL_SX2X2 (f30,f31,f20,f21,f20,f21);
        f40=NEXP01_SX2(vs0); f41=NEXP01_SX2(vs1);
        CONST_SX2X2(f00,f01, 3);
        ADDEXP_SX2(f40, f00); ADDEXP_SX2(f41, f01);
        MADD_SX2X2 (f00,f01, f30,f31,f40,f41);
        f30=NEXP01_SX2(vs0); f31=NEXP01_SX2(vs1);
        MADD_SX2X2 (f20,f21,f00,f01, f20,f21);
        MUL_SX2X2(f00,f01, f30,f31, f20,f21); NEG_SX2X2(f00,f01,f00,f01);
        MUL_SX2X2 (f60,f61, f20,f21, f40,f41);
        CONST_SX2X2(f40,f41, 3);
        MUL_SX2X2  (f70,f71, f40,f41, f20,f21);
        MADD_SX2X2 (f30,f31, f00,f01, f00,f01);
        MADD_SX2X2 (f40,f41, f60,f61, f20,f21);
        NEG_SX2X2 (f70,f71, f70,f71);
        MADD_SX2X2 (f00,f01, f30,f31, f70,f71);
        MSUB_SX2X2 (f70,f71, f40,f41, f70,f71);
        f30=f70; f31=f71;
        f20=MKSADJ_SX2(vs0); f21=MKSADJ_SX2(vs1);
        f10=NEXP01_SX2(vs0); f11=NEXP01_SX2(vs1);
        MADD_SX2X2(f10,f11, f00,f01, f00,f01);
        ADDEXPM_SX2(f00, f20); ADDEXPM_SX2(f01, f21);
        ADDEXP_SX2 (f30, f20); ADDEXP_SX2 (f31, f21);
        DIVN_SX2(f00,f10, f30); DIVN_SX2(f01,f11, f31);
        *pvd0=f00; *pvd1=f01;
}
#define SQRT_SX2X2(vd0,vd1,vs0,vs1) sqrt_sx2x2(&vd0,&vd1,&vs0,&vs1);
#endif

#if (HAVE_VFPU==0 && HAVE_FPU==0)
DISCARD_FUN(void ,vec_complex2mag,(float32_t  * restrict y, const complex_float  * restrict x, int N))
#elif (HAVE_VFPU)
/*===========================================================================
  Vector matematics:
  vec_complex2invmag     complex magnitude  
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
void       vec_complex2mag (float32_t  * restrict y, const complex_float  * restrict x, int N)
{
    int n;
    ae_valign aY;
    ae_valignx2 aX;
    static const union {  uint32_t  u;  float32_t f; }
                 xmin={0x20000000}, // sqrt(realmin)
                 xmax={0x5f000000}, // sqrt(realmax)/2
                 k   ={0x6a800000}, // direct and inverse scale factors 
                 kinv={0x14800000}; //
    NASSERT( x );
    NASSERT( y );
    if (N <= 0) return;

    aY=AE_ZALIGN64();
    aX=AE_LA128_PP(x);
    for ( n=0; n<(N>>2); n++)
    {
        xtbool2 bsmall,bbig;
        xtfloatx2 reim0,reim1,rere,imim,maxreim,r0,r1,t0,t1;
        xtfloatx2 s,sinv0,sinv1;
        AE_LASX2X2_IP(reim0,reim1,aX,castxcc(xtfloatx4,x));
        ABS_SX2X2(reim0,reim1,reim0,reim1);
        rere=AE_SEL32_HH_SX2(reim0,reim1);
        imim=AE_SEL32_LL_SX2(reim0,reim1);
        maxreim=MAXNUM_SX2(rere,imim);
        bsmall=XT_OLE_SX2(maxreim,xmin.f);
        bbig  =XT_OLE_SX2(xmax.f,maxreim);
        CONST_SX2X2(sinv0,s,1);
        XT_MOVT_SX2(s,k.f,bsmall);
        XT_MOVT_SX2(sinv0,kinv.f,bsmall);
        XT_MOVT_SX2(s,kinv.f,bbig);
        XT_MOVT_SX2(sinv0,k.f,bbig);
        MUL_SX2X2(rere,imim,rere,imim,s,s);
        MUL_SX2X2(r0,t0,rere,imim,rere,imim);

        AE_LASX2X2_IP(reim0,reim1,aX,castxcc(xtfloatx4,x));
        ABS_SX2X2(reim0,reim1,reim0,reim1);
        rere=AE_SEL32_HH_SX2(reim0,reim1);
        imim=AE_SEL32_LL_SX2(reim0,reim1);
        maxreim=MAXNUM_SX2(rere,imim);
        bsmall=XT_OLE_SX2(maxreim,xmin.f);
        bbig  =XT_OLE_SX2(xmax.f,maxreim);
        CONST_SX2X2(sinv1,s,1); 
        XT_MOVT_SX2(s,k.f,bsmall);
        XT_MOVT_SX2(sinv1,kinv.f,bsmall);
        XT_MOVT_SX2(s,kinv.f,bbig);
        XT_MOVT_SX2(sinv1,k.f,bbig);
        MUL_SX2X2(rere,imim,rere,imim,s,s);
        MUL_SX2X2(r1,t1,rere,imim,rere,imim);

        ADD_SX2X2(r0,r1,r0,r1,t0,t1);
        SQRT_SX2X2(r0,r1,r0,r1);
        MUL_SX2X2(r0,r1,r0,r1,sinv0,sinv1);

        XT_SASX2IP(r0,aY,castxcc(xtfloatx2,y));
        XT_SASX2IP(r1,aY,castxcc(xtfloatx2,y));
    }
    AE_SA64POS_FP(aY,castxcc(xtfloatx2,y));
    // tail
    __Pragma("loop_count max=3")
    for ( n=0; n<(N&3); n++ )
    {
        xtbool2 bsmall,bbig;
        xtfloatx2 reim;
        xtfloat r, maxreim;
        xtfloatx2 s,sinv;
        XT_LSX2IP(reim,castxcc(xtfloatx2,x),sizeof(xtfloatx2));
        reim=XT_ABS_SX2(reim);
        maxreim=RMAXNUM_S(reim);
        bsmall=XT_OLE_SX2(maxreim,xmin.f);
        bbig  =XT_OLE_SX2(xmax.f,maxreim);
        sinv=s=XT_CONST_S(1); 
        XT_MOVT_SX2(s,k.f,bsmall);
        XT_MOVT_SX2(sinv,kinv.f,bsmall);
        XT_MOVT_SX2(s,kinv.f,bbig);
        XT_MOVT_SX2(sinv,k.f,bbig);
        reim=XT_MUL_SX2(reim,s);
        reim=XT_MUL_SX2(reim,reim);
        r=XT_RADD_SX2(reim);
        r=XT_SQRT_S(r);
        r=XT_MUL_S(r,sinv);
        XT_SSIP(r,castxcc(xtfloat,y),sizeof(float32_t));
    }
}
#else
// for scalar FPU
void       vec_complex2mag (float32_t  * restrict y, const complex_float  * restrict x, int N)
{
    const int blkSize = (MAX_ALLOCA_SZ/(sizeof(xtfloat)))/20;
    xtfloat scr[blkSize];
    int n,NN;
    const xtfloat* restrict pX=(const xtfloat*)x;
    const xtfloat* restrict pZrd;
          xtfloat* restrict pZwr;
    const xtfloat* restrict pYrd=(const xtfloat*)y;
          xtfloat* restrict pYwr=(      xtfloat*)y;
    NN=N;
    while(NN>0)
    {
        __Pragma("no_reorder")
        N=XT_MIN(NN,blkSize);
        pZwr=scr;
        pYwr=(      xtfloat*)y;
        for (n=0; n<N; n++)
        {
            xtfloat x0, x1, y0;
            xtfloat xre, xim;
            int t0, t1, nsa;
            int e0;
            int nsa0;

            XT_LSIP(xre,pX,1*sizeof(xtfloat));
            XT_LSIP(xim,pX,1*sizeof(xtfloat));
            xre = XT_ABS_S(xre);
            xim = XT_ABS_S(xim);
            t0 = XT_RFR(xre);
            t1 = XT_RFR(xim);
            nsa = XT_MAX(t0, t1);
            nsa = ((uint32_t)nsa)>> 23;
            nsa = (nsa-127);
            nsa = XT_MIN(nsa, 127);
            e0 = (127-nsa);
            nsa0 = (e0<<23);
            XT_MOVEQZ(nsa0,0x00400000,e0);
            y0 = XT_WFR(nsa0);

            xre = XT_MUL_S(xre, y0);
            xim = XT_MUL_S(xim, y0);

            x0 = XT_MUL_S(xre, xre);
            x1 = XT_MUL_S(xim, xim);

            x0 = XT_ADD_S(x0, x1);
            XT_SSIP(x0,pYwr,sizeof(xtfloat));

            e0 = (127+nsa);
            nsa0 = (e0<<23);
            XT_MOVEQZ(nsa0, 0x00400000, e0);
            x0 = XT_WFR(nsa0);
            XT_SSIP(x0,pZwr,sizeof(xtfloat));
        }
        __Pragma("no_reorder")
        pZrd=scr;
        pYrd=(const xtfloat*)y;
        pYwr=(      xtfloat*)y;
        for (n=0; n<N; n++)
        {
            xtfloat y0,x0;
            XT_LSIP(y0,pYrd,sizeof(xtfloat));
            y0 = XT_SQRT_S(y0);
            XT_LSIP(x0,pZrd,sizeof(xtfloat));
            y0 = XT_MUL_S(y0, x0);    
            XT_SSIP(y0,pYwr,sizeof(xtfloat));
        }
        NN-=N;
        y+=N;
    }
}
#endif

