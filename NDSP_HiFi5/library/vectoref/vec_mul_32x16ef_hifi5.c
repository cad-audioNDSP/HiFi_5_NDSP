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
    Vector multiply for Emulated Floating Point
    optimized code for HiFi5 core
*/
/* Library API */
#include "NatureDSP_Signal_vector.h"
#include "NatureDSP_types.h"
#include "common.h"

#define CONTROL_UNDERFLOW 0
#define CONTROL_OVERFLOW  1

// code optimized for HiFi5 core
/*-------------------------------------------------------------------------
  Vector Multiply for Emulated Floating Point
  routines multiply two vectors represented in emulated floating point format

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
void   vec_mul_32x16ef (      int32_t  * restrict zmant,       int16_t  * restrict zexp, 
                        const int32_t  * restrict xmant, const int16_t  * restrict xexp, 
                        const int32_t  * restrict ymant, const int16_t  * restrict yexp, 
                        int N)
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
        ae_int32x2 xm,ym;
        ae_int32x2 exp,temp;
        ae_int32x2 nsa;
        ae_int64 zm;
        ae_int16x4 xe,ye;
        int nsa0;
        AE_L32_IP(xm,castxcc(ae_int32,pXmant),sizeof(int32_t));
        AE_L32_IP(ym,castxcc(ae_int32,pYmant),sizeof(int32_t));
        AE_L16_IP(xe,castxcc(ae_int16,pXexp ),sizeof(int16_t));
        AE_L16_IP(ye,castxcc(ae_int16,pYexp ),sizeof(int16_t));
        zm=AE_MUL32S_HH(xm,ym);
        AE_ADDW16(exp,temp,xe,ye); 
        nsa0=AE_NSA64(zm);
        xm=AE_TRUNCA32X2F64S(zm,zm,nsa0); 
        nsa=AE_MOVDA32X2(nsa0,nsa0);
        nsa=AE_ADD32(AE_SUB32(exp,nsa),1);
        // underflow processing
        AE_MOVT32X2(xm,AE_ZERO32(),AE_LT32(nsa,MIN_INT16));
        /* overflow processing */
        AE_MOVT32X2(xm,AE_SLLI32S(xm,31),AE_LT32(MAX_INT16,nsa));
        AE_S32_L_IP(xm,castxcc(ae_int32,pZmant),sizeof(int32_t));
        AE_S16_0_IP(AE_SAT16X4(nsa,nsa),castxcc(ae_int16,pZexp),sizeof(int16_t));
    }
    N&=~3;
    aZmant=AE_ZALIGN128();
    aZexp =AE_ZALIGN64();
    for (n=0; n<N; n+=4)
    {
        ae_int32x2 xm01,xm23,ym01,ym23;
        ae_int32x2 exp0,exp1;
        ae_int32x2 nsa01,nsa23;
        ae_int64 zm0,zm1,zm2,zm3;
        ae_int16x4 xe,ye,nsa;
        int nsa0,nsa1,nsa2,nsa3;
        aXmant=AE_LA128_PP(pXmant);
        aYmant=AE_LA128_PP(pYmant);
        aXexp =AE_LA64_PP (pXexp); 
        aYexp =AE_LA64_PP (pYexp);
        AE_LA32X2X2_IP(xm01,xm23,aXmant,pXmant);
        AE_LA32X2X2_IP(ym01,ym23,aYmant,pYmant);
        AE_LA16X4_IP(xe,aXexp,pXexp);
        AE_LA16X4_IP(ye,aYexp,pYexp);
        AE_MUL32X2S_HH_LL(zm0,zm1,xm01,ym01);
        AE_MUL32X2S_HH_LL(zm2,zm3,xm23,ym23);
        AE_ADDW16(exp0,exp1,xe,ye); 
        nsa0=AE_NSA64(zm0);
        nsa1=AE_NSA64(zm1);
        nsa2=AE_NSA64(zm2);
        nsa3=AE_NSA64(zm3);
        xm01=AE_SEL32_HH(AE_TRUNCA32X2F64S(zm0,zm0,nsa0),AE_TRUNCA32X2F64S(zm1,zm1,nsa1));
        xm23=AE_SEL32_HH(AE_TRUNCA32X2F64S(zm2,zm2,nsa2),AE_TRUNCA32X2F64S(zm3,zm3,nsa3));
        nsa01=AE_MOVDA32X2(nsa0,nsa1);
        nsa23=AE_MOVDA32X2(nsa2,nsa3);
        nsa01=AE_ADD32(AE_SUB32(exp0,nsa01),1);
        nsa23=AE_ADD32(AE_SUB32(exp1,nsa23),1);
#if CONTROL_UNDERFLOW
        // underflow processing
        AE_MOVT32X2(xm01,AE_ZERO32(),AE_LT32(nsa01,MIN_INT16));
        AE_MOVT32X2(xm23,AE_ZERO32(),AE_LT32(nsa23,MIN_INT16));
#endif
#if CONTROL_OVERFLOW
        /* overflow processing */
        AE_MOVT32X2(xm01,AE_SLLI32S(xm01,31),AE_LT32(MAX_INT16,nsa01));
        AE_MOVT32X2(xm23,AE_SLLI32S(xm23,31),AE_LT32(MAX_INT16,nsa23));
#endif
        nsa=AE_SAT16X4(nsa01,nsa23);
        AE_SA32X2X2_IP(xm01,xm23,aZmant,pZmant);
        AE_SA16X4_IP(nsa,aZexp,pZexp);
    }
    AE_SA128POS_FP(aZmant,pZmant);
    AE_SA64POS_FP (aZexp ,pZexp );
}


