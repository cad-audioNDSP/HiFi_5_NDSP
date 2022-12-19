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
    Vector addition for Emulated Floating Point
    optimized code for HiFi5 core
*/
/* Library API */
#include "NatureDSP_Signal_vector.h"
#include "NatureDSP_types.h"
#include "common.h"
// code optimized for HiFi5 core
#define CONTROL_UNDERFLOW 0
#define CONTROL_OVERFLOW  1
/*-------------------------------------------------------------------------
  Vector Addition for Emulated Floating Point
  routines add two vectors represented in emulated floating point format

  Input:
  xmant[N]  mantissa of input data
  ymant[N]  mantissa of input data
  xexp[N]   exponent of input data
  yexp[N]   exponent of input data
  N         length of vectors
  Output:
  zmant[N]  mantissa of output data
  zexp[N]   exponent of output data

  Restriction:
  xmant,ymant,xexp,yexp,zmant,zexp should not overlap
-------------------------------------------------------------------------*/
void   vec_add_32x16ef (      int32_t  *  restrict zmant,       int16_t  *  restrict zexp, 
                        const int32_t  *  restrict xmant, const int16_t  *  restrict xexp, 
                        const int32_t  *  restrict ymant, const int16_t  *  restrict yexp, 
                        int N)
#if 0
{
    ae_valign aXmant,aYmant,aZmant;
    ae_int16x4 nsa16;
    xtbool2 ufl;
    ae_int64 a0,a1;
    ae_int32x2 x_exp,y_exp,exp,x_mant,y_mant,nsa,x;
    int x_exp0,y_exp0,nsa0;
    int x_exp1,y_exp1,nsa1;
    int n;

    if (N<=0) return;
    /* first, process unaligned odd element in scalar way */
    if(N&1)
    {
        x_exp0=*xexp++;
        y_exp0=*yexp++;
        x_exp=AE_MOVDA32(x_exp0);
        y_exp=AE_MOVDA32(y_exp0);

        AE_L32_IP(x_mant,castxcc(ae_int32,xmant),sizeof(int32_t));
        AE_L32_IP(y_mant,castxcc(ae_int32,ymant),sizeof(int32_t));
        /* set exponents to zeroes if mantissa is zero */
        AE_MOVT32X2(x_exp,MIN_INT16,AE_EQ32(x_mant,0));
        AE_MOVT32X2(y_exp,MIN_INT16,AE_EQ32(y_mant,0));
        exp=AE_MAX32(x_exp,y_exp);
        /* add with exponents */
        x_exp=AE_SUB32(x_exp,exp);
        y_exp=AE_SUB32(y_exp,exp);
        AE_MOVSARD7(AE_MAX32(-31,x_exp)); x_mant=AE_SLAS32S(x_mant);
        AE_MOVSARD7(AE_MAX32(-31,y_exp)); y_mant=AE_SLAS32S(y_mant);
        a0=AE_MUL32_HH(x_mant,AE_MOVDA32(1));
        AE_MULA32_HH(a0,y_mant,AE_MOVDA32(1));
        /* normalization */
        nsa0=AE_NSA64(a0);
        nsa=AE_MOVDA32(nsa0);
        x=AE_TRUNCA32X2F64S(a0,a0,nsa0);
        nsa=AE_SUB32(exp,AE_SUB32(nsa,32));
        /* underflow processing */
        ufl=AE_LT32(nsa,MIN_INT16);
        AE_MOVT32X2(x,0,ufl);
        /* overflow processing */
        AE_MOVT32X2(x,AE_SLLI32S(x,31),AE_LT32(MAX_INT16,nsa));
        AE_S32_L_IP(x,castxcc(ae_int32,zmant),sizeof(int32_t));
        nsa16=AE_SAT16X4(nsa,nsa);
        AE_S16_0_IP((nsa16),castxcc(ae_int16,zexp),sizeof(int16_t));
    }
    N&=~1;
    // process by 2 numbers in parallel
    aXmant=AE_LA64_PP(castxcc(ae_int32x2,xmant));
    aYmant=AE_LA64_PP(castxcc(ae_int32x2,ymant));
    aZmant=AE_ZALIGN64();
    for (n=0; n<N; n+=2) 
    {
        x_exp0=*xexp++,y_exp0=*yexp++;
        x_exp1=*xexp++,y_exp1=*yexp++;
        x_exp=AE_MOVDA32X2(x_exp0,x_exp1);
        y_exp=AE_MOVDA32X2(y_exp0,y_exp1);

        AE_LA32X2_IP(x_mant,aXmant,castxcc(ae_int32x2,xmant));
        AE_LA32X2_IP(y_mant,aYmant,castxcc(ae_int32x2,ymant));
        /* set exponents to zeroes if mantissa is zero */
        AE_MOVT32X2(x_exp,MIN_INT16,AE_EQ32(x_mant,0));
        AE_MOVT32X2(y_exp,MIN_INT16,AE_EQ32(y_mant,0));
        exp=AE_MAX32(x_exp,y_exp);
        /* add with exponents */
        x_exp=AE_SUB32(x_exp,exp);
        y_exp=AE_SUB32(y_exp,exp);
        AE_MOVSARD7(AE_MAX32(-31,x_exp)); x_mant=AE_SLAS32S(x_mant);
        AE_MOVSARD7(AE_MAX32(-31,y_exp)); y_mant=AE_SLAS32S(y_mant);
        a0=AE_MUL32_HH(x_mant,AE_MOVDA32(1));
        a1=AE_MUL32_LH(x_mant,AE_MOVDA32(1));
        AE_MULA32_HH(a0,y_mant,AE_MOVDA32(1));
        AE_MULA32_LH(a1,y_mant,AE_MOVDA32(1));
        /* normalization */
        nsa0=AE_NSA64(a0);
        nsa1=AE_NSA64(a1);
        nsa=AE_MOVDA32X2(nsa0,nsa1);
        x=AE_SEL32_HH(AE_TRUNCA32F64S(a0,nsa0),AE_TRUNCA32F64S(a1,nsa1));
        nsa=AE_SUB32(exp,AE_SUB32(nsa,32));
        /* underflow processing */
        ufl=AE_LT32(nsa,MIN_INT16);
        AE_MOVT32X2(x,0,ufl);
        /* overflow processing */
        AE_MOVT32X2(x,AE_SLLI32S(x,31),AE_LT32(MAX_INT16,nsa));
        AE_SA32X2_IP(x,aZmant,castxcc(ae_int32x2,zmant));
        nsa16=AE_SAT16X4(nsa,nsa);
        AE_S16_0_IP(AE_SHORTSWAP(nsa16),castxcc(ae_int16,zexp),sizeof(int16_t));
        AE_S16_0_IP((nsa16),castxcc(ae_int16,zexp),sizeof(int16_t));
    }
    AE_SA64POS_FP(aZmant,castxcc(ae_int32x2,zmant));
}
#else
{
          ae_int32x4 * restrict pZmant;
    const ae_int32x4 * restrict pXmant;
    const ae_int32x4 * restrict pYmant;
    const ae_int16x4 * restrict pXexp;
    const ae_int16x4 * restrict pYexp;
          ae_int16x4 * restrict pZexp;
    ae_valignx2 aZmant,aXmant,aYmant;
    ae_valign aXexp,aYexp,aZexp;
    int n;
    if (N<=0) return;
    pZmant=(      ae_int32x4 *)zmant;
    pXmant=(const ae_int32x4 *)xmant;
    pYmant=(const ae_int32x4 *)ymant;
    pXexp =(const ae_int16x4 *)xexp ;
    pYexp =(const ae_int16x4 *)yexp ;
    pZexp =(      ae_int16x4 *)zexp ;
    __Pragma("no_unroll")
    for (n=0; n<(N&3); n++)
    {
        ae_int32x2 xm0,ym0;
        ae_int16x4 xexp,yexp,nsa0123,exp0123;
        ae_int32x2 nsa0,temp;
        ae_int32x2 xexp0;
        ae_int32x2 yexp0;
        ae_int16x4 vexp,vnsa;
        AE_L32_IP(xm0 ,castxcc(ae_int32,pXmant),sizeof(int32_t));
        AE_L32_IP(ym0 ,castxcc(ae_int32,pYmant),sizeof(int32_t));
        AE_L16_IP(xexp,castxcc(ae_int16,pXexp ),sizeof(int16_t));
        AE_L16_IP(yexp,castxcc(ae_int16,pYexp ),sizeof(int16_t));
#if 0
        AE_MOVT16X4(xexp,MIN_INT16,AE_EQ16(AE_SAT16X4(xm0,xm0),AE_ZERO16()));
        AE_MOVT16X4(yexp,MIN_INT16,AE_EQ16(AE_SAT16X4(ym0,ym0),AE_ZERO16()));
#else
        AE_EXPADD16_L(xexp, 0, xm0);
        AE_EXPADD16_L(yexp, 0, ym0);
#endif
        exp0123=AE_MAX16(xexp,yexp);
        xexp=AE_SUB16S(exp0123,xexp);
        yexp=AE_SUB16S(exp0123,yexp);
        xexp0=AE_SEXT32X2D16_32(xexp);
        yexp0=AE_SEXT32X2D16_32(yexp);
        xm0=AE_SRAV32RS(xm0,xexp0);
        ym0=AE_SRAV32RS(ym0,yexp0);

        vexp=0;
        AE_ADDCEXP32_H(xm0,vexp,xm0,ym0);
        vnsa=AE_NEG16S(AE_NSA32X4(xm0,xm0));
        xexp0=AE_SEXT32X2D16_32(vnsa);
        xm0=AE_SRAV32RS(xm0,xexp0);
        vexp=AE_ADD16(vexp,vnsa);
        AE_ADDW16(nsa0,temp,exp0123,vexp);
#if CONTROL_UNDERFLOW
        // underflow processing
        AE_MOVT32X2(xm0,AE_ZERO32(),AE_LT32(nsa0,MIN_INT16));
#endif
#if CONTROL_OVERFLOW
        /* overflow processing */
        ym0=AE_SLAI32S(xm0,31);
        AE_MOVT32X2(xm0,ym0,AE_LT32(MAX_INT16,nsa0));
#endif
        nsa0123=AE_SAT16X4(nsa0,nsa0);
        AE_S32_L_IP(xm0    ,castxcc(ae_int32,pZmant),sizeof(int32_t));
        AE_S16_0_IP(nsa0123,castxcc(ae_int16,pZexp ),sizeof(int16_t));
    }
    N&=~3;
    aZmant=AE_ZALIGN128();
    aZexp =AE_ZALIGN64();
    for (n=0; n<(N>>2); n++)
    {
        ae_int32x2 xm01,xm23,ym01,ym23;
        ae_int16x4 xexp,yexp,nsa0123,exp0123;
        ae_int32x2 nsa01,nsa23;
        ae_int32x2 xexp01,xexp23;
        ae_int32x2 yexp01,yexp23;
        ae_int16x4 vexp,vnsa;
        aXmant=AE_LA128_PP(pXmant);
        aYmant=AE_LA128_PP(pYmant);
        aXexp =AE_LA64_PP (pXexp); 
        aYexp =AE_LA64_PP (pYexp);
        AE_LA32X2X2_IP(xm01,xm23,aXmant,pXmant);
        AE_LA32X2X2_IP(ym01,ym23,aYmant,pYmant);
        AE_LA16X4_IP(xexp,aXexp,pXexp);
        AE_LA16X4_IP(yexp,aYexp,pYexp);
        // set min exponent for all numbers with zero mantissa to avoid 
        // bad scaling of second argument
#if 0
        AE_MOVT16X4(xexp,MIN_INT16,AE_EQ16(AE_SAT16X4(xm01,xm23),AE_ZERO16()));
        AE_MOVT16X4(yexp,MIN_INT16,AE_EQ16(AE_SAT16X4(ym01,ym23),AE_ZERO16()));
#else
        AE_EXPADD16_H(xexp, 0, xm01);
        AE_EXPADD16_L(xexp, 0, xm23);
        AE_EXPADD16_H(yexp, 0, ym01);
        AE_EXPADD16_L(yexp, 0, ym23);
#endif
        // denormalization
        exp0123=AE_MAX16(xexp,yexp);
        xexp=AE_SUB16S(exp0123,xexp);
        yexp=AE_SUB16S(exp0123,yexp);
        xexp01=AE_SEXT32X2D16_32(xexp);
        xexp23=AE_SEXT32X2D16_10(xexp);
        yexp01=AE_SEXT32X2D16_32(yexp);
        yexp23=AE_SEXT32X2D16_10(yexp);
        xm01=AE_SRAV32RS(xm01,xexp01);
        xm23=AE_SRAV32RS(xm23,xexp23);
        ym01=AE_SRAV32RS(ym01,yexp01);
        ym23=AE_SRAV32RS(ym23,yexp23);
        // add with exponent biasing and normalization
        vexp=0;
        AE_ADDCEXP32_H(xm01,vexp,xm01,ym01);
        AE_ADDCEXP32_L(xm23,vexp,xm23,ym23);
        vnsa=AE_NEG16S(AE_NSA32X4(xm01,xm23));
        xexp01=AE_SEXT32X2D16_32(vnsa);
        xexp23=AE_SEXT32X2D16_10(vnsa);
        xm01=AE_SRAV32RS(xm01,xexp01);
        xm23=AE_SRAV32RS(xm23,xexp23);
        vexp=AE_ADD16(vexp,vnsa);
        AE_ADDW16(nsa01,nsa23,exp0123,vexp);
#if CONTROL_UNDERFLOW
        // underflow processing
        AE_MOVT32X2(xm01,AE_ZERO32(),AE_LT32(nsa01,MIN_INT16));
        AE_MOVT32X2(xm23,AE_ZERO32(),AE_LT32(nsa23,MIN_INT16));
#endif
#if CONTROL_OVERFLOW
        /* overflow processing */
        ym01=AE_SLAI32S(xm01,31);
        ym23=AE_SLAI32S(xm23,31);
        AE_MOVT32X2(xm01,ym01,AE_LT32(MAX_INT16,nsa01));
        AE_MOVT32X2(xm23,ym23,AE_LT32(MAX_INT16,nsa23));
#endif
        nsa0123=AE_SAT16X4(nsa01,nsa23);
        AE_SA32X2X2_IP(xm01,xm23,aZmant,pZmant);
        AE_SA16X4_IP(nsa0123,aZexp,pZexp);
    }
    AE_SA128POS_FP(aZmant,pZmant);
    AE_SA64POS_FP (aZexp ,pZexp );
}
#endif
