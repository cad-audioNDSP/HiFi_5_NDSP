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
    Vector multiply-accumulate for Emulated Floating Point
    optimized code for HiFi5 core
*/
/* Library API */
#include "NatureDSP_Signal_vector.h"
#include "common.h"

/* 
ALG==1 uses single normalization but bigger amount of 64-bit operations
ALG==2 uses 2 normalizations but utilized HiFi5 AE_ADDCEXP32_H/L instructions
for acceleration of additions
*/
#define ALG 2 
#define CONTROL_UNDERFLOW 0
#define CONTROL_OVERFLOW  1
/*-------------------------------------------------------------------------
  Vector Multiply-Accumulate for Emulated Floating Point
  routines multiply-accumulate vectors by scalar represented in emulated 
  floating point format

  Input:
  xmant[N]  mantissa of input data
  xexp[N]   exponent of input data
  ymant     mantissa of scalar
  yexp      exponent of scalar
  N         length of vectors
  Output:
  zmant[N]  mantissa of output data
  zexp[N]   exponent of output data

  Restriction:
  xmant,xexp,zmant,zexp should not overlap
-------------------------------------------------------------------------*/
void   vec_mac_32x16ef (      int32_t  *  restrict zmant,       int16_t  * restrict zexp, 
                        const int32_t  *  restrict xmant, const int16_t  * restrict xexp, 
                              int32_t              ymant,       int16_t             yexp, 
                        int N)
#if ALG==1
{
    const ae_int32x4 *restrict pZmantRd;
    const ae_int32x4 *restrict pXmant;
          ae_int32x4 *restrict pZmantWr;
    const ae_int16x4 *restrict pZexpRd;
    const ae_int16x4 *restrict pXexp;
          ae_int16x4 *restrict pZexpWr;
    ae_valignx2 aZmantRd,aXmant,aZmantWr;
    ae_valign   aZexpRd,aXexp,aZexpWr;
    int n;
    if(N<=0) return;
    pZmantRd=(const ae_int32x4 *)zmant;
    pXmant  =(const ae_int32x4 *)xmant;
    pZmantWr=(      ae_int32x4 *)zmant;
    pZexpRd =(const ae_int16x4 *)zexp ;
    pXexp   =(const ae_int16x4 *)xexp ;
    pZexpWr =(      ae_int16x4 *)zexp ;

    for (n=0; n<(N&3); n++)
    {
        ae_int32x2 aexp0,x0,cmant0,cexp0,exp0,nsa0,tmp;
        ae_int16x4 xexp ,cexp;
        ae_int64 a0;
        ae_int64 b0;
        int nsa;
        AE_L32_IP(x0    ,castxcc(ae_int32,pXmant  ),sizeof(int32_t));
        AE_L32_IP(cmant0,castxcc(ae_int32,pZmantRd),sizeof(int32_t));
        AE_L16_IP(xexp  ,castxcc(ae_int16,pXexp   ),sizeof(int16_t));
        AE_L16_IP(cexp  ,castxcc(ae_int16,pZexpRd ),sizeof(int16_t));
        AE_ADDW16(aexp0,tmp,xexp,yexp);
        a0=AE_MUL32_HH(x0,ymant);
        // set output exponent to the minimum value if resulted mantissa is zero
        AE_MOVT32X2(aexp0,MIN_INT16,AE_EQ32(AE_ZERO32(),AE_SAT32X2(a0,a0)));
        // the same with accumulator
        AE_MOVT16X4(cexp,MIN_INT16,AE_EQ16(AE_ZERO16(),AE_SAT16X4(cmant0,cmant0)));
        AE_ADDW16(cexp0,tmp,cexp,31);
        exp0=AE_MAX32(cexp0,aexp0);
        aexp0=AE_SUB32(aexp0,exp0);
        cexp0=AE_SUB32(cexp0,exp0);
        a0 = AE_SLAA64S(a0,AE_MOVAD32_H(aexp0));
        b0=AE_MUL32_HH(cmant0,1);
        b0 = AE_SLAA64S(b0,AE_MOVAD32_H(cexp0));
        a0 = AE_ADD64S(a0,b0);
        /* normalize output */
        nsa=AE_NSA64(a0);
        x0=AE_SEL32_HH(AE_TRUNCA32F64S(a0,nsa),AE_TRUNCA32F64S(a0,nsa));
        nsa0=AE_ADD32(1,AE_SUB32(exp0,AE_MOVDA32X2(nsa,nsa)));
        /* underflow processing */
        AE_MOVT32X2(x0,AE_ZERO32(),AE_LT32(nsa0,MIN_INT16));
        /* overflow processing */
        AE_MOVT32X2(x0,AE_SLAI32S(x0,31),AE_LT32(MAX_INT16,nsa0));
        AE_S32_L_IP(x0,castxcc(ae_int32,pZmantWr),sizeof(int32_t));
        AE_S16_0_IP(AE_SAT16X4(nsa0,nsa0),castxcc(ae_int16,pZexpWr),sizeof(int16_t));
    }
    N&=~3;
    aZexpWr=AE_ZALIGN64();
    aZmantWr=AE_ZALIGN128();
    for (n=0; n<N; n+=4)
    {
        ae_int32x2 aexp01,aexp23,x01,x23,cmant01,cmant23,cexp01,cexp23,exp01,exp23,nsa01,nsa23;
        ae_int16x4 xexp0123,cexp0123;
        ae_int64 a0,a1,a2,a3;
        ae_int64 b0,b1,b2,b3;
        int nsa0,nsa1,nsa2,nsa3;
        // load inputs and multiply mantissas
        aXexp   =AE_LA64_PP(pXexp  );
        aZexpRd =AE_LA64_PP(pZexpRd);
        aXmant  =AE_LA128_PP(pXmant  );
        aZmantRd=AE_LA128_PP(pZmantRd);
        AE_LA32X2X2_IP(x01    ,x23    ,aXmant  ,pXmant  );
        AE_LA32X2X2_IP(cmant01,cmant23,aZmantRd,pZmantRd);
        AE_LA16X4_IP(xexp0123,aXexp,pXexp);
        AE_LA16X4_IP(cexp0123,aZexpRd,pZexpRd);
        AE_ADDW16(aexp01,aexp23,xexp0123,yexp);
        AE_MUL32X2S_HH_LL(a0,a1,x01,ymant);
        AE_MUL32X2S_HH_LL(a2,a3,x23,ymant);
        // set output exponent to the minimum value if resulted mantissa is zero
        AE_MOVT32X2(aexp01,MIN_INT16,AE_EQ32(AE_ZERO32(),AE_SAT32X2(a0,a1)));
        AE_MOVT32X2(aexp23,MIN_INT16,AE_EQ32(AE_ZERO32(),AE_SAT32X2(a2,a3)));
        // the same with accumulator
        AE_MOVT16X4(cexp0123,MIN_INT16,AE_EQ16(AE_ZERO16(),AE_SAT16X4(cmant01,cmant23)));
        AE_ADDW16(cexp01,cexp23,cexp0123,31);
        exp01=AE_MAX32(cexp01,aexp01);
        exp23=AE_MAX32(cexp23,aexp23);
        aexp01=AE_SUB32(aexp01,exp01);
        aexp23=AE_SUB32(aexp23,exp23);
        cexp01=AE_SUB32(cexp01,exp01);
        cexp23=AE_SUB32(cexp23,exp23);
        a0 = AE_SLAA64S(a0,AE_MOVAD32_H(aexp01));
        a1 = AE_SLAA64S(a1,AE_MOVAD32_L(aexp01));
        a2 = AE_SLAA64S(a2,AE_MOVAD32_H(aexp23));
        a3 = AE_SLAA64S(a3,AE_MOVAD32_L(aexp23));
        AE_MUL32X2S_HH_LL(b0,b1,cmant01,1);
        AE_MUL32X2S_HH_LL(b2,b3,cmant23,1);
        b0 = AE_SLAA64S(b0,AE_MOVAD32_H(cexp01));
        b1 = AE_SLAA64S(b1,AE_MOVAD32_L(cexp01));
        b2 = AE_SLAA64S(b2,AE_MOVAD32_H(cexp23));
        b3 = AE_SLAA64S(b3,AE_MOVAD32_L(cexp23));
        a0 = AE_ADD64S(a0,b0);
        a1 = AE_ADD64S(a1,b1);
        a2 = AE_ADD64S(a2,b2);
        a3 = AE_ADD64S(a3,b3);
        /* normalize output */
        nsa0=AE_NSA64(a0);
        nsa1=AE_NSA64(a1);
        nsa2=AE_NSA64(a2);
        nsa3=AE_NSA64(a3);
        x01=AE_SEL32_HH(AE_TRUNCA32F64S(a0,nsa0),AE_TRUNCA32F64S(a1,nsa1));
        x23=AE_SEL32_HH(AE_TRUNCA32F64S(a2,nsa2),AE_TRUNCA32F64S(a3,nsa3));
        nsa01=AE_ADD32(1,AE_SUB32(exp01,AE_MOVDA32X2(nsa0,nsa1)));
        nsa23=AE_ADD32(1,AE_SUB32(exp23,AE_MOVDA32X2(nsa2,nsa3)));
        /* underflow processing */
        AE_MOVT32X2(x01,AE_ZERO32(),AE_LT32(nsa01,MIN_INT16));
        AE_MOVT32X2(x23,AE_ZERO32(),AE_LT32(nsa23,MIN_INT16));
        /* overflow processing */
        AE_MOVT32X2(x01,AE_SLAI32S(x01,31),AE_LT32(MAX_INT16,nsa01));
        AE_MOVT32X2(x23,AE_SLAI32S(x23,31),AE_LT32(MAX_INT16,nsa23));
        AE_SA32X2X2_IP(x01,x23,aZmantWr,pZmantWr);
        AE_SA16X4_IP(AE_SAT16X4(nsa01,nsa23),aZexpWr,pZexpWr);
    }
    AE_SA128POS_FP(aZmantWr,pZmantWr);
    AE_SA64POS_FP (aZexpWr ,pZexpWr );
}
#elif ALG==2
{
    const ae_int32x4 *restrict pZmantRd;
    const ae_int32x4 *restrict pXmant;
          ae_int32x4 *restrict pZmantWr;
    const ae_int16x4 *restrict pZexpRd;
    const ae_int16x4 *restrict pXexp;
          ae_int16x4 *restrict pZexpWr;
    ae_valignx2 aZmantRd,aXmant,aZmantWr;
    ae_valign   aZexpRd,aXexp,aZexpWr;
    int n;
    if (ymant==0) yexp=MIN_INT16;
    if (N<=0) return;
    pZmantRd=(const ae_int32x4 *)zmant;
    pXmant  =(const ae_int32x4 *)xmant;
    pZmantWr=(      ae_int32x4 *)zmant;
    pZexpRd =(const ae_int16x4 *)zexp ;
    pXexp   =(const ae_int16x4 *)xexp ;
    pZexpWr =(      ae_int16x4 *)zexp ;
    for (n=0; n<(N&3); n++)
    {
        ae_int16x4 vexp,vnsa;
        ae_int32x2 xm0,am0,temp;
        ae_int32x2 aexp0,nsa0;
        ae_int16x4 xexp0,ae0,exp0;
        int nsa;
        ae_int64 a0;
        AE_L32_IP(xm0  ,castxcc(ae_int32,pXmant),sizeof(int32_t));
        AE_L16_IP(xexp0,castxcc(ae_int16,pXexp ),sizeof(int16_t));
        // multiply and normalize
        a0=AE_MUL32_HH(xm0,ymant);
        nsa=AE_NSA64(a0);
        am0=AE_TRUNCA32F64S(a0,nsa);
        nsa0=AE_MOVDA32X2(nsa,nsa);
        AE_ADDW16(aexp0,temp,xexp0,yexp);
        nsa0=AE_ADD32(1,AE_SUB32(aexp0,nsa0));
        ae0 =AE_SAT16X4(nsa0,nsa0);
        AE_EXPADD16_L(ae0, 0, am0);
        AE_L32_IP (xm0  ,castxcc(ae_int32,pZmantRd),sizeof(int32_t));
        AE_L16_IP (xexp0,castxcc(ae_int16,pZexpRd ),sizeof(int16_t));
        AE_EXPADD16_L(xexp0, 0, xm0);
        exp0=AE_MAX16(xexp0,ae0);
        xexp0=AE_SUB16S(exp0,xexp0);
        ae0  =AE_SUB16S(exp0,ae0);
        xm0=AE_SRAV32RS(xm0,AE_SEXT32X2D16_32(xexp0));
        am0=AE_SRAV32RS(am0,AE_SEXT32X2D16_32(ae0));
        // add vectors with exponent biasing
        vexp=AE_ZERO16();
        AE_ADDCEXP32_H(am0,vexp,xm0,am0);
        vnsa=AE_NEG16S(AE_NSA32X4(am0,am0));
        am0=AE_SRAV32RS(am0,AE_SEXT32X2D16_32(vnsa));
        vexp=AE_ADD16(vexp,vnsa);
        AE_ADDW16(nsa0,temp,exp0,vexp);
#if CONTROL_UNDERFLOW
        // underflow processing
        AE_MOVT32X2(am0,0,AE_LT32(nsa0,MIN_INT16));
#endif
#if CONTROL_OVERFLOW
        // overflow processing
        AE_MOVT32X2(am0,AE_SLAI32S(am0,31),AE_LT32(MAX_INT16,nsa0));
#endif
        AE_S32_L_IP(am0,castxcc(ae_int32,pZmantWr),sizeof(int32_t));
        AE_S16_0_IP(AE_SAT16X4(nsa0,nsa0),castxcc(ae_int16,pZexpWr),sizeof(int16_t));
    }
    N&=~3;
    aZmantWr =AE_ZALIGN128();
    aZexpWr  =AE_ZALIGN64 ();
    for (n=0; n<(N>>2); n++)
    {
        ae_int16x4 vexp,vnsa;
        ae_int32x2 xm01,xm23,am01,am23;
        ae_int32x2 aexp01,aexp23,nsa01,nsa23;
        ae_int16x4 xexp0123,ae0123,exp0123;
        int nsa0,nsa1,nsa2,nsa3;
        ae_int64 a0,a1,a2,a3;
        aXexp   =AE_LA64_PP(pXexp  );
        aZexpRd =AE_LA64_PP(pZexpRd);
        aXmant  =AE_LA128_PP(pXmant  );
        aZmantRd=AE_LA128_PP(pZmantRd);
        AE_LA32X2X2_IP(xm01    ,xm23    ,aXmant  ,pXmant  );
        AE_LA16X4_IP(xexp0123,aXexp,pXexp);
        // multiply and normalize
        AE_MUL32X2S_HH_LL(a0,a1,xm01,ymant);
        AE_MUL32X2S_HH_LL(a2,a3,xm23,ymant);
        nsa0=AE_NSA64(a0);
        nsa1=AE_NSA64(a1);
        nsa2=AE_NSA64(a2);
        nsa3=AE_NSA64(a3);
        am01=AE_SEL32_HH(AE_TRUNCA32F64S(a0,nsa0),AE_TRUNCA32F64S(a1,nsa1));
        am23=AE_SEL32_HH(AE_TRUNCA32F64S(a2,nsa2),AE_TRUNCA32F64S(a3,nsa3));
        nsa01=AE_MOVDA32X2(nsa0,nsa1);
        nsa23=AE_MOVDA32X2(nsa2,nsa3);
        AE_ADDW16(aexp01,aexp23,xexp0123,yexp);
        nsa01=AE_ADD32(1,AE_SUB32(aexp01,nsa01));
        nsa23=AE_ADD32(1,AE_SUB32(aexp23,nsa23));
        ae0123 =AE_SAT16X4(nsa01,nsa23);
        AE_EXPADD16_H(ae0123, 0, am01);
        AE_EXPADD16_L(ae0123, 0, am23);
        AE_LA32X2X2_IP(xm01,xm23,aZmantRd,pZmantRd);
        AE_LA16X4_IP  (xexp0123,aZexpRd,pZexpRd);
        // denormalize
        AE_EXPADD16_H(xexp0123, 0, xm01);
        AE_EXPADD16_L(xexp0123, 0, xm23);
        exp0123=AE_MAX16(xexp0123,ae0123);
        xexp0123=AE_SUB16S(exp0123,xexp0123);
        ae0123  =AE_SUB16S(exp0123,ae0123  );
        xm01=AE_SRAV32RS(xm01,AE_SEXT32X2D16_32(xexp0123));
        xm23=AE_SRAV32RS(xm23,AE_SEXT32X2D16_10(xexp0123));
        am01=AE_SRAV32RS(am01,AE_SEXT32X2D16_32(ae0123));
        am23=AE_SRAV32RS(am23,AE_SEXT32X2D16_10(ae0123));
        // add vectors with exponent biasing
        vexp=AE_ZERO16();
        AE_ADDCEXP32_H(am01,vexp,xm01,am01);
        AE_ADDCEXP32_L(am23,vexp,xm23,am23);
        vnsa=AE_NEG16S(AE_NSA32X4(am01,am23));
        am01=AE_SRAV32RS(am01,AE_SEXT32X2D16_32(vnsa));
        am23=AE_SRAV32RS(am23,AE_SEXT32X2D16_10(vnsa));
        vexp=AE_ADD16(vexp,vnsa);
        AE_ADDW16(nsa01,nsa23,exp0123,vexp);
#if CONTROL_UNDERFLOW
        // underflow processing
        AE_MOVT32X2(am01,0,AE_LT32(nsa01,MIN_INT16));
        AE_MOVT32X2(am23,0,AE_LT32(nsa23,MIN_INT16));
#endif
#if CONTROL_OVERFLOW
        // overflow processing
        AE_MOVT32X2(am01,AE_SLAI32S(am01,31),AE_LT32(MAX_INT16,nsa01));
        AE_MOVT32X2(am23,AE_SLAI32S(am23,31),AE_LT32(MAX_INT16,nsa23));
#endif
        AE_SA32X2X2_IP(am01,am23,aZmantWr,pZmantWr);
        AE_SA16X4_IP(AE_SAT16X4(nsa01,nsa23),aZexpWr,pZexpWr);
    }
    AE_SA128POS_FP(aZmantWr,pZmantWr);
    AE_SA64POS_FP (aZexpWr ,pZexpWr );
}
#endif
