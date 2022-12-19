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
  NatureDSP Signal Processing Library. Mathematics
    Vector operations
    Optimized code for HiFi5
  IntegrIT, 2006-2019
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Signal_math.h"
#include "common.h"

/*-------------------------------------------------------------------------
  Reciprocal Square Root
  These routines return the fractional and exponential portion of the 
  reciprocal square root of a vector x of Q31 or Q15 numbers. Since the 
  reciprocal square root is always greater than 1, they return fractional 
  portion frac in Q(31-exp) or Q(15-exp) format and exponent exp so true 
  reciprocal value in the Q0.31/Q0.15 may be found by shifting fractional 
  part left by exponent value.

  Mantissa accuracy is 1 LSB, so relative accuracy is:
  vec_rsqrt16x16, scl_rsqrt16x16	6.2e-5
  scl_rsqrt32x32	                2.4e-7
  vec_rsqrt32x32	                9.2e-10

  Precision: 
  16x16  16-bit inputs, 16-bit output. Accuracy: 2LSB
  32x32  32-bit inputs, 32-bit output. Accuracy: (2.6e-7*y+1LSB)

  Input:
  x[N]     input data, Q15, Q31 
  N        length of vectors
  Output:
  frac[N]  fractional part of result, Q(31-exp) or Q(15-exp)
  exp[N]   exponent of result 

  Restriction:
  x, fract, exp - should not overlap

  Scalar versions:
  ----------------
  Returned packed value: 
  scl_rsqrt32x32():
  bits 23…0 fractional part
  bits 31…24 exponent
  scl_rsqrt16x16():
  bits 15…0 fractional part
  bits 31…16 exponent

-------------------------------------------------------------------------*/
void vec_rsqrt16x16 (
                  int16_t * frac, 
                  int16_t * exp,
            const int16_t *  x, 
            int N)
{
    const ae_int16x8 * restrict pX   ;
          ae_int16x8 * restrict pFrac;
          ae_int16x8 * restrict pExp ;
    ae_valignx2 aX,aF,aE;
    int n;
    if (N<=0) return;
    pX   =(const ae_int16x8 *)x   ;
    pFrac=(      ae_int16x8 *)frac;
    pExp =(      ae_int16x8 *)exp ;
    for (n=0; n<(N&7); n++)
    {
        xtbool4 bltzero0;
        ae_int16x4 x0,y0,e0,nsa0;
        ae_int16x4 mx0;
        ae_int32x2 w0,w1;
        ae_int32x2 t0,t1;
        AE_L16_IP(x0,castxcc(ae_int16,pX),sizeof(int16_t));
        bltzero0=AE_LT16(x0,0);
        nsa0=AE_NSA16X4(x0);
        AE_MOVT16X4(nsa0,28,AE_LE16(x0,0));
        // normalize
        x0=AE_SRAV16RS(x0,AE_NEG16S(AE_AND16(nsa0,~1)));
        nsa0=AE_ADD16(1,AE_SRAI16(nsa0,1));
        AE_S16_0_IP(nsa0,castxcc(ae_int16,pExp),sizeof(int16_t));
        // first approximation
        y0=AE_SUB16(0x7fff,AE_SRAI16(x0,1));
        AE_MOVT16X4(x0,0,bltzero0);
        AE_MOVT16X4(y0,0x8000,bltzero0);
        mx0=AE_NEG16S(x0);
        // first 2 iterations
        e0=AE_MULFD16X16X4RAS(AE_MOVDA16(0x7fff),mx0,AE_MOVDA16(8192),AE_MULFP16X4RAS(y0,y0));
        e0=AE_ADD16(e0,e0);
        y0=AE_MULFD16X16X4RAS(y0,e0,AE_MOVDA16(0x7fff),y0);
        e0=AE_MULFD16X16X4RAS(AE_MOVDA16(0x7fff),mx0,AE_MOVDA16(8192),AE_MULFP16X4RAS(y0,y0));
        e0=AE_ADD16(e0,e0);
        y0=AE_MULFD16X16X4RAS(y0,e0,AE_MOVDA16(0x7fff),y0);
        // last iteration is done in better accuracy
        AE_MULF16X4SS(w0,w1,y0,y0);
        t0=t1=(8192<<16); AE_MULSF2P32X16X4RAS(t0,t1,w0,w1,x0); w0=t0; w1=t1;
        AE_CVTI32X4F16(t0,t1,y0,15);
        AE_MULAF2P32X16X4RAS(t0,t1,w0,w1,y0);
        y0=AE_TRUNCI16X4F32S(t0,t1,1);
        AE_S16_0_IP(y0,castxcc(ae_int16,pFrac),sizeof(int16_t));
    }
    N&=~7;
    aX=AE_LA128_PP(pX);
    aF=AE_ZALIGN128();
    aE=AE_ZALIGN128();
    for (n=0; n<(N>>3); n++)
    {
        xtbool4 bltzero0,bltzero1;
        ae_int16x4 x0,x1,y0,y1,e0,e1,nsa0,nsa1;
        ae_int16x4 mx0,mx1;
        ae_int32x2 w0,w1,w2,w3;
        ae_int32x2 t0,t1,t2,t3;
        AE_LA16X4X2_IP(x0,x1,aX,pX);
        bltzero0=AE_LT16(x0,0);
        bltzero1=AE_LT16(x1,0);
        nsa0=AE_NSA16X4(x0);
        nsa1=AE_NSA16X4(x1);
        AE_MOVT16X4(nsa0,28,AE_LE16(x0,0));
        AE_MOVT16X4(nsa1,28,AE_LE16(x1,0));
        // normalize
        x0=AE_SRAV16RS(x0,AE_NEG16S(AE_AND16(nsa0,~1)));
        x1=AE_SRAV16RS(x1,AE_NEG16S(AE_AND16(nsa1,~1)));
        nsa0=AE_ADD16(1,AE_SRAI16(nsa0,1));
        nsa1=AE_ADD16(1,AE_SRAI16(nsa1,1));
        AE_SA16X4X2_IP(nsa0,nsa1,aE,pExp);
        // first approximation
        y0=AE_SUB16(0x7fff,AE_SRAI16(x0,1));
        y1=AE_SUB16(0x7fff,AE_SRAI16(x1,1));
        AE_MOVT16X4(x0,0,bltzero0);
        AE_MOVT16X4(x1,0,bltzero1);
        AE_MOVT16X4(y0,0x8000,bltzero0);
        AE_MOVT16X4(y1,0x8000,bltzero1);
        mx0=AE_NEG16S(x0);
        mx1=AE_NEG16S(x1);
        // first 2 iterations
        e0=AE_MULFD16X16X4RAS(AE_MOVDA16(0x7fff),mx0,AE_MOVDA16(8192),AE_MULFP16X4RAS(y0,y0));
        e1=AE_MULFD16X16X4RAS(AE_MOVDA16(0x7fff),mx1,AE_MOVDA16(8192),AE_MULFP16X4RAS(y1,y1));
        e0=AE_ADD16(e0,e0);
        e1=AE_ADD16(e1,e1);
        y0=AE_MULFD16X16X4RAS(y0,e0,AE_MOVDA16(0x7fff),y0);
        y1=AE_MULFD16X16X4RAS(y1,e1,AE_MOVDA16(0x7fff),y1);
        e0=AE_MULFD16X16X4RAS(AE_MOVDA16(0x7fff),mx0,AE_MOVDA16(8192),AE_MULFP16X4RAS(y0,y0));
        e1=AE_MULFD16X16X4RAS(AE_MOVDA16(0x7fff),mx1,AE_MOVDA16(8192),AE_MULFP16X4RAS(y1,y1));
        e0=AE_ADD16(e0,e0);
        e1=AE_ADD16(e1,e1);
        y0=AE_MULFD16X16X4RAS(y0,e0,AE_MOVDA16(0x7fff),y0);
        y1=AE_MULFD16X16X4RAS(y1,e1,AE_MOVDA16(0x7fff),y1);
        // last iteration is done in better accuracy
        AE_MULF16X4SS(w0,w1,y0,y0);
        AE_MULF16X4SS(w2,w3,y1,y1);
        t0=t1=(8192<<16); AE_MULSF2P32X16X4RAS(t0,t1,w0,w1,x0); w0=t0; w1=t1;
        t2=t3=(8192<<16); AE_MULSF2P32X16X4RAS(t2,t3,w2,w3,x1); w2=t2; w3=t3;
        AE_CVTI32X4F16(t0,t1,y0,15);
        AE_CVTI32X4F16(t2,t3,y1,15);
        AE_MULAF2P32X16X4RAS(t0,t1,w0,w1,y0);
        AE_MULAF2P32X16X4RAS(t2,t3,w2,w3,y1);
        y0=AE_TRUNCI16X4F32S(t0,t1,1);
        y1=AE_TRUNCI16X4F32S(t2,t3,1);
        AE_SA16X4X2_IP(y0,y1,aF,pFrac);
    }
    AE_SA128POS_FP(aF,pFrac);
    AE_SA128POS_FP(aE,pExp);
} /* vec_rsqrt16x16() */
