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

/* DSP Library API */
#include "NatureDSP_Signal_math.h"
#include "common.h"
//    C code optimized for HiFi5

/*===========================================================================
  Vector matematics:
  vec_recip            Reciprocal on Q31/Q15 Numbers
===========================================================================*/
#define ALG 1 // 0 - 32-bit multiplies, 1 - 16-bit multiplies

/*-------------------------------------------------------------------------
  Reciprocal on Q63/Q31/Q15 Numbers
  These routines return the fractional and exponential portion of the 
  reciprocal of a vector x of Q31 or Q15 numbers. Since the reciprocal is 
  always greater than 1, it returns fractional portion frac in Q(31-exp) 
  or Q(15-exp) format and exponent exp so true reciprocal value in the 
  Q0.31/Q0.15 may be found by shifting fractional part left by exponent 
  value.

  Mantissa accuracy is 1 LSB, so relative accuracy is:
  vec_recip16x16, scl_recip16x16                   6.2e-5 
  scl_recip32x32                                   2.4e-7 
  vec_recip32x32                                   9.2e-10
  vec_recip64x64                                   2.2e-19

  Precision: 
  64x64  64-bit input, 64-bit output. 
  32x32  32-bit input, 32-bit output. 
  16x16  16-bit input, 16-bit output. 

  Input:
  x[N]    input data, Q63, Q31 or Q15
  N       length of vectors

  Output:
  frac[N] fractional part of result, Q(63-exp), Q(31-exp) or Q(15-exp)
  exp[N]  exponent of result 

  Restriction:
  x,frac,exp should not overlap

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  frac,x - aligned on 16-byte boundary
  N      - multiple of 4 and >4

  Scalar versions:
  ----------------
  Return packed value: 
  scl_recip64x64():
  bits 55:0 fractional part
  bits 63:56 exponent
  scl_recip32x32():
  bits 23:0 fractional part
  bits 31:24 exponent
  scl_recip16x16():
  bits 15:0 fractional part
  bits 31:16 exponent
-------------------------------------------------------------------------*/
void vec_recip16x16
(
  int16_t * restrict frac, 
  int16_t * exp, 
  const int16_t * restrict x, 
  int N
)
#if ALG==0  
// code with 32-bit multiplies
{
          ae_int16x8* restrict pFrac=(      ae_int16x8*)frac;
          ae_int16x8* restrict pExp =(      ae_int16x8*)exp;
    const ae_int16x8* restrict pX   =(const ae_int16x8*)x;
    ae_valignx2 ae,ax,ay;
    ae_int16x4 x0,x1;
    ae_int16x4 e0,e1;
    ae_int16x4 y0,y1;
    xtbool4 sx0,sx1;
    int n,N0;
    if (N<=0) return;
    N0=N&7;
    pFrac=(      ae_int16x8*)frac;
    pExp =(      ae_int16x8*)exp;
    pX   =(const ae_int16x8*)x;
    for (n=0; n<N0; n++)
    {
        AE_L16_IP(x0,castxcc(ae_int16,pX),sizeof(int16_t));
        e0=AE_NSA16X4(x0);
        AE_S16_0_IP(AE_ADD16(e0,1),castxcc(ae_int16,pExp),sizeof(int16_t) );
        x0=AE_SRAV16RS(x0,AE_NEG16S(e0));
        AE_MOVT16X4(x0,0x4000,AE_EQ16(x0,0));
        sx0=AE_LT16(x0,0);
        x0=AE_ABS16S(x0);
        y0=AE_SUB16((int16_t)47852,x0);
        e0=AE_SUB16(0x4000,AE_MULFP16X4S(x0,y0));
        e0=AE_SLLI16S(e0,1);
        y0=AE_ADD16S(y0,AE_MULFP16X4S(e0,y0));
        e0=AE_SUB16(0x4000,AE_MULFP16X4S(x0,y0));
        e0=AE_SLLI16S(e0,1);
        y0=AE_ADD16S(y0,AE_MULFP16X4S(e0,y0));
        e0=AE_SUB16(0x4000,AE_MULFP16X4S(x0,y0));
        e0=AE_SLLI16S(e0,1);
        y0=AE_ADD16S(y0,AE_MULFP16X4S(e0,y0));
        AE_MOVT16X4(y0,AE_NEG16S(y0),sx0);
        AE_S16_0_IP(y0,castxcc(ae_int16,pFrac),sizeof(int16_t));
    }

    N-=N0;

    pX   =(const ae_int16x8 *)(x+N0);
    pFrac=(      ae_int16x8 *)(frac+N0);
    pExp =(      ae_int16x8 *)(exp+N0);
    ax=AE_LA128_PP(pX);
    ay=AE_ZALIGN128();
    ae=AE_ZALIGN128();
    for (n=0; n<(N>>3); n++) 
    {
        ae_int32x2 z0,z1,z2,z3;
        ae_int32x2 w0,w1,w2,w3;
        AE_LA16X4X2_IP(x0,x1,ax,pX);
        e0=AE_NSA16X4(x0);
        e1=AE_NSA16X4(x1);
        AE_SA16X4X2_IP(AE_ADD16(e0,1),AE_ADD16(e1,1),ae,pExp );
        x0=AE_SRAV16RS(x0,AE_NEG16S(e0));
        x1=AE_SRAV16RS(x1,AE_NEG16S(e1));
        AE_MOVT16X4(x0,0x4000,AE_EQ16(x0,0));
        AE_MOVT16X4(x1,0x4000,AE_EQ16(x1,0));
                /* first approximation */
        sx0=AE_LT16(x0,0);
        sx1=AE_LT16(x1,0);
        x0=AE_ABS16S(x0);
        x1=AE_ABS16S(x1);
        y0=AE_SUB16((int16_t)47852,x0);
        y1=AE_SUB16((int16_t)47852,x1);
        AE_CVTI32X4F16S(z0,z1,y0,16);
        AE_CVTI32X4F16S(z2,z3,y1,16);
        w0=w1=w2=w3=0x40000000; 
        AE_MULSF2P32X16X4S(w0,w1,z0,z1,x0);
        AE_MULSF2P32X16X4S(w2,w3,z2,z3,x1);
        AE_MUL2P32X4(w0,w1,w0,w1,2,2);
        AE_MUL2P32X4(w2,w3,w2,w3,2,2);
        AE_MULAF2P32X4RS(z0,z1,z0,z1,w0,w1);
        AE_MULAF2P32X4RS(z2,z3,z2,z3,w2,w3);
        w0=w1=w2=w3=0x40000000; 
        AE_MULSF2P32X16X4S(w0,w1,z0,z1,x0);
        AE_MULSF2P32X16X4S(w2,w3,z2,z3,x1);
        AE_MUL2P32X4(w0,w1,w0,w1,2,2);
        AE_MUL2P32X4(w2,w3,w2,w3,2,2);
        AE_MULAF2P32X4RS(z0,z1,z0,z1,w0,w1);
        AE_MULAF2P32X4RS(z2,z3,z2,z3,w2,w3);
        w0=w1=w2=w3=0x40000000; 
        AE_MULSF2P32X16X4S(w0,w1,z0,z1,x0);
        AE_MULSF2P32X16X4S(w2,w3,z2,z3,x1);
        AE_MUL2P32X4(w0,w1,w0,w1,2,2);
        AE_MUL2P32X4(w2,w3,w2,w3,2,2);
        AE_MULAF2P32X4RS(z0,z1,z0,z1,w0,w1);
        AE_MULAF2P32X4RS(z2,z3,z2,z3,w2,w3);
        y0=AE_TRUNCA16X4F32S(z0,z1,0);
        y1=AE_TRUNCA16X4F32S(z2,z3,0);

        AE_MOVT16X4(y0,AE_NEG16S(y0),sx0);
        AE_MOVT16X4(y1,AE_NEG16S(y1),sx1);
        AE_SA16X4X2_IP(y0,y1,ay,pFrac);
    }
    AE_SA128POS_FP(ay,pFrac);
    AE_SA128POS_FP(ae,pExp );
}
#elif ALG==1  
// code with 16-bit multiplies
{
          ae_int16x8* restrict pFrac=(      ae_int16x8*)frac;
          ae_int16x8* restrict pExp =(      ae_int16x8*)exp;
    const ae_int16x8* restrict pX   =(const ae_int16x8*)x;
    ae_valignx2 ae,ax,ay;
    ae_int16x4 x0,x1;
    ae_int16x4 mx0,mx1;
    ae_int16x4 e0,e1;
    ae_int16x4 y0,y1;
    xtbool4 sx0,sx1;
    int n,N0;
    if (N<=0) return;
    N0=N&7;
    pFrac=(      ae_int16x8*)frac;
    pExp =(      ae_int16x8*)exp;
    pX   =(const ae_int16x8*)x;
    for (n=0; n<N0; n++)
    {
        AE_L16_IP(x0,castxcc(ae_int16,pX),sizeof(int16_t));
        e0=AE_NSA16X4(x0);
        AE_S16_0_IP(AE_ADD16(e0,1),castxcc(ae_int16,pExp),sizeof(int16_t) );
        x0=AE_SRAV16RS(x0,AE_NEG16S(e0));
        AE_MOVT16X4(x0,0x4000,AE_EQ16(x0,0));
        sx0=AE_LT16(x0,0);
        x0=AE_ABS16S(x0);
        y0=AE_SUB16((int16_t)47852,x0);
        mx0=AE_NEG16S(x0);
        e0=AE_MULFD16X16X4RAS(AE_MOVDA16(0x7fff),mx0,AE_MOVDA16(0x4000),y0);
        e0=AE_ADD16(e0,e0);
        y0=AE_MULFD16X16X4RAS(y0,e0,AE_MOVDA16(0x7fff),y0);

        e0=AE_MULFD16X16X4RAS(AE_MOVDA16(0x7fff),mx0,AE_MOVDA16(0x4000),y0);
        e0=AE_ADD16(e0,e0);
        y0=AE_MULFD16X16X4RAS(y0,e0,AE_MOVDA16(0x7fff),y0);

        e0=AE_SUB16(0x4000,AE_MULFP16X4S(x0,y0));
        e0=AE_ADD16(e0,e0);
        y0=AE_MULFD16X16X4RAS(y0,e0,AE_MOVDA16(0x7fff),y0);
        AE_MOVT16X4(y0,AE_NEG16S(y0),sx0);
        AE_S16_0_IP(y0,castxcc(ae_int16,pFrac),sizeof(int16_t));
    }

    N-=N0;

    pX   =(const ae_int16x8 *)(x+N0);
    pFrac=(      ae_int16x8 *)(frac+N0);
    pExp =(      ae_int16x8 *)(exp+N0);
    ax=AE_LA128_PP(pX);
    ay=AE_ZALIGN128();
    ae=AE_ZALIGN128();
    for (n=0; n<(N>>3); n++) 
    {
        AE_LA16X4X2_IP(x0,x1,ax,pX);
        e0=AE_NSA16X4(x0);
        e1=AE_NSA16X4(x1);
        AE_SA16X4X2_IP(AE_ADD16(e0,1),AE_ADD16(e1,1),ae,pExp );
        x0=AE_SRAV16RS(x0,AE_NEG16S(e0));
        x1=AE_SRAV16RS(x1,AE_NEG16S(e1));
        AE_MOVT16X4(x0,0x4000,AE_EQ16(x0,0));
        AE_MOVT16X4(x1,0x4000,AE_EQ16(x1,0));
                /* first approximation */
        sx0=AE_LT16(x0,0);
        sx1=AE_LT16(x1,0);
        x0=AE_ABS16S(x0);
        x1=AE_ABS16S(x1);
        y0=AE_SUB16((int16_t)47852,x0);
        y1=AE_SUB16((int16_t)47852,x1);
        mx0=AE_NEG16S(x0);
        mx1=AE_NEG16S(x1);
        e0=AE_MULFD16X16X4RAS(AE_MOVDA16(0x7fff),mx0,AE_MOVDA16(0x4000),y0);
        e1=AE_MULFD16X16X4RAS(AE_MOVDA16(0x7fff),mx1,AE_MOVDA16(0x4000),y1);
        e0=AE_ADD16(e0,e0);
        e1=AE_ADD16(e1,e1);
        y0=AE_MULFD16X16X4RAS(y0,e0,AE_MOVDA16(0x7fff),y0);
        y1=AE_MULFD16X16X4RAS(y1,e1,AE_MOVDA16(0x7fff),y1);

        e0=AE_MULFD16X16X4RAS(AE_MOVDA16(0x7fff),mx0,AE_MOVDA16(0x4000),y0);
        e1=AE_MULFD16X16X4RAS(AE_MOVDA16(0x7fff),mx1,AE_MOVDA16(0x4000),y1);
        e0=AE_ADD16(e0,e0);
        e1=AE_ADD16(e1,e1);
        y0=AE_MULFD16X16X4RAS(y0,e0,AE_MOVDA16(0x7fff),y0);
        y1=AE_MULFD16X16X4RAS(y1,e1,AE_MOVDA16(0x7fff),y1);

        e0=AE_SUB16(0x4000,AE_MULFP16X4S(x0,y0));
        e1=AE_SUB16(0x4000,AE_MULFP16X4S(x1,y1));
        e0=AE_ADD16(e0,e0);
        e1=AE_ADD16(e1,e1);
        y0=AE_MULFD16X16X4RAS(y0,e0,AE_MOVDA16(0x7fff),y0);
        y1=AE_MULFD16X16X4RAS(y1,e1,AE_MOVDA16(0x7fff),y1);

        AE_MOVT16X4(y0,AE_NEG16S(y0),sx0);
        AE_MOVT16X4(y1,AE_NEG16S(y1),sx1);
        AE_SA16X4X2_IP(y0,y1,ay,pFrac);
    }
    AE_SA128POS_FP(ay,pFrac);
    AE_SA128POS_FP(ae,pExp );
} /* vec_recip16x16() */
#endif
