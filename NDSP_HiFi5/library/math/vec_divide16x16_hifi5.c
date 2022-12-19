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
    Division
    Optimized code for HiFi5
  IntegrIT, 2006-2019
*/
#include "NatureDSP_Signal_math.h"
#include "common.h"


/*===========================================================================
  Vector matematics:
  vec_divide           Division of Q31/Q15 Numbers
===========================================================================*/

/*-------------------------------------------------------------------------
  Division
  These routines perform pair-wise division of vectors written in Q63, Q31 or 
  Q15 format. They return the fractional and exponential portion of the division 
  result. Since the division may generate result greater than 1, it returns 
  fractional portion frac in Q(63-exp), Q(31-exp) or Q(15-exp) format and 
  exponent exp so true division result in the Q0.31 may be found by shifting 
  fractional part left by exponent value.
  Additional routine makes integer division of 64-bit number to 32-bit 
  denominator forming 32-bit result. If result is overflown, 0x7fffffff 
  or 0x80000000 is returned depending on the signs of inputs.
  For division to 0, the result is not defined.

  Two versions of routines are available: regular versions (vec_divide64x32i,
  vec_divide64x64, vec_divide32x32, vec_divide16x16) work 
  with arbitrary arguments, faster versions (vec_divide32x32_fast, 
  vec_divide16x16_fast) apply some restrictions.

  Accuracy is measured as accuracy of fractional part (mantissa):
  vec_divide64x32i, scl_divide64x32                      :  1 LSB   
  vec_divide64x64                                        :  2 LSB 
  vec_divide32x32, vec_divide32x32_fast                  :  2 LSB (1.8e-9) 
  scl_divide32x32                                        :  2 LSB (4.8e-7) 
  vec_divide16x16, scl_divide16x16, vec_divide16x16_fast :  2 LSB (1.2e-4)

  Precision: 
  64x32i integer division, 64-bit nominator, 32-bit denominator, 32-bit output. 
  64x64  fractional division, 64-bit inputs, 64-bit output. 
  32x32  fractional division, 32-bit inputs, 32-bit output. 
  16x16  fractional division, 16-bit inputs, 16-bit output. 

  Input:
  x[N]    nominator, 64-bit integer, Q63, Q31 or Q15
  y[N]    denominator, 32-bit integer, Q63, Q31 or Q15
  N       length of vectors
  Output:
  frac[N] fractional parts of result, Q(63-exp), Q(31-exp) or Q(15-exp)
  exp[N]  exponents of result 

  Restriction:
  For regular versions (vec_divide64x32i, vec_divide64x64, vec_divide32x32,
  vec_divide16x16) :
  x,y,frac,exp should not overlap

  For faster versions (vec_divide32x32_fast, vec_divide16x16_fast) :
  x,y,frac,exp  should not overlap
  x,y,frac      to be aligned by 16-byte boundary, N - multiple of 4.

  Scalar versions:
  ----------------
  scl_divide64x32(): integer remainder
  Return packed value: 
  scl_divide64x64():
  bits 55:0 fractional part
  bits 63:56 exponent
  scl_divide32x32():
  bits 23:0 fractional part
  bits 31:24 exponent
  scl_divide16x16():
  bits 15:0 fractional part
  bits 31:16 exponent
-------------------------------------------------------------------------*/
void vec_divide16x16 
(
  int16_t *       restrict  frac,
  int16_t *       restrict  exp,
  const int16_t * restrict  x,
  const int16_t * restrict  y,
  int M)
{
    const ae_int16x8 *restrict pX;
    const ae_int16x8 *restrict pY;
          ae_int16x8 *restrict pF;
          ae_int16x8 *restrict pE;
    ae_valignx2 aX,aY,aF,aE;
    int n;
    if (M<=0) return;
    pX=(const ae_int16x8*)x;
    pY=(const ae_int16x8*)y;
    pF=(ae_int16x8*)frac;
    pE=(ae_int16x8*)exp ;
    if (M&7)
    {
        for (n=0; n<(M&7); n++)
        {
            ae_int16x4 x0,y0,z0,e0,expx0,expy0;
            ae_int16x4 my0;
            xtbool4 sy0;
            AE_L16_IP(x0,castxcc(ae_int16,pX),sizeof(int16_t));
            AE_L16_IP(y0,castxcc(ae_int16,pY),sizeof(int16_t));
            // normalization
            expx0 = AE_NSA16X4(x0);
            expy0 = AE_NSA16X4(y0);
            x0=AE_SRAV16RS(x0,AE_NEG16S(expx0));
            y0=AE_SRAV16RS(y0,AE_NEG16S(expy0));
            // save resulted exponent
            expy0=AE_ADD16(1,AE_SUB16(expy0,expx0));
            AE_S16_0_IP(expy0,castxcc(ae_int16,pE),sizeof(int16_t));
            sy0=AE_LT16(y0,AE_ZERO16());
            y0=AE_ABS16S(y0);
            my0=AE_NEG16S(y0);
            /* first approximation */
            z0=AE_SUB16(AE_MOVDA16((int16_t)47852),y0); 
            /* 3 iterations to achieve 1 LSB accuracy in mantissa */
            e0=AE_MULFD16X16X4RAS(AE_MOVDA16(0x4000),my0,AE_MOVDA16(0x7fff),z0);
            e0=AE_ADD16(e0,e0);
            z0=AE_MULFD16X16X4RAS(z0,z0,AE_MOVDA16(0x7fff),e0);
            e0=AE_MULFD16X16X4RAS(AE_MOVDA16(0x4000),my0,AE_MOVDA16(0x7fff),z0);
            e0=AE_ADD16(e0,e0);
            z0=AE_MULFD16X16X4RAS(z0,z0,AE_MOVDA16(0x7fff),e0);
            e0=AE_SUB16(0x4000,AE_MULFP16X4S(y0,z0)); 
            e0=AE_ADD16(e0,e0);
            z0=AE_MULFD16X16X4RAS(z0,z0,AE_MOVDA16(0x7fff),e0);
            /* restore original sign */
            y0=AE_NEG16S(z0);
            AE_MOVT16X4(z0,y0,sy0);
            /* multiply by X */
            z0=AE_MULFP16X4RAS(x0,z0);
            AE_S16_0_IP(z0,castxcc(ae_int16,pF),sizeof(int16_t));
        }    
    }
    M&=~7;
    aX=AE_LA128_PP(pX);
    aY=AE_LA128_PP(pY);
    aF=aE=AE_ZALIGN128();
    for (n=0; n<(M>>3); n++)
    {
        ae_int16x4 x0,x1,y0,y1,z0,z1,e0,e1,expx0,expx1,expy0,expy1;
        ae_int16x4 my0,my1;
        xtbool4 sy0,sy1;
        AE_LA16X4X2_IP(x0,x1,aX,pX);
        AE_LA16X4X2_IP(y0,y1,aY,pY);
        // normalization
        expx0 = AE_NSA16X4(x0);
        expx1 = AE_NSA16X4(x1);
        expy0 = AE_NSA16X4(y0);
        expy1 = AE_NSA16X4(y1);
        x0=AE_SRAV16RS(x0,AE_NEG16S(expx0));
        x1=AE_SRAV16RS(x1,AE_NEG16S(expx1));
        y0=AE_SRAV16RS(y0,AE_NEG16S(expy0));
        y1=AE_SRAV16RS(y1,AE_NEG16S(expy1));
        // save resulted exponent
        expy0=AE_ADD16(1,AE_SUB16(expy0,expx0));
        expy1=AE_ADD16(1,AE_SUB16(expy1,expx1));
        AE_SA16X4X2_IP(expy0,expy1,aE,pE);
        sy0=AE_LT16(y0,AE_ZERO16());
        sy1=AE_LT16(y1,AE_ZERO16());
        y0=AE_ABS16S(y0);
        y1=AE_ABS16S(y1);
        my0=AE_NEG16S(y0);
        my1=AE_NEG16S(y1);
        /* first approximation */
        z0=AE_SUB16(AE_MOVDA16((int16_t)47852),y0); 
        z1=AE_SUB16(AE_MOVDA16((int16_t)47852),y1); 
        /* 3 iterations to achieve 1 LSB accuracy in mantissa */
        e0=AE_MULFD16X16X4RAS(AE_MOVDA16(0x4000),my0,AE_MOVDA16(0x7fff),z0);
        e1=AE_MULFD16X16X4RAS(AE_MOVDA16(0x4000),my1,AE_MOVDA16(0x7fff),z1);
        e0=AE_ADD16(e0,e0);
        e1=AE_ADD16(e1,e1);
        z0=AE_MULFD16X16X4RAS(z0,z0,AE_MOVDA16(0x7fff),e0);
        z1=AE_MULFD16X16X4RAS(z1,z1,AE_MOVDA16(0x7fff),e1);
        e0=AE_MULFD16X16X4RAS(AE_MOVDA16(0x4000),my0,AE_MOVDA16(0x7fff),z0);
        e1=AE_MULFD16X16X4RAS(AE_MOVDA16(0x4000),my1,AE_MOVDA16(0x7fff),z1);
        e0=AE_ADD16(e0,e0);
        e1=AE_ADD16(e1,e1);
        z0=AE_MULFD16X16X4RAS(z0,z0,AE_MOVDA16(0x7fff),e0);
        z1=AE_MULFD16X16X4RAS(z1,z1,AE_MOVDA16(0x7fff),e1);
        e0=AE_SUB16(0x4000,AE_MULFP16X4S(y0,z0)); 
        e1=AE_SUB16(0x4000,AE_MULFP16X4S(y1,z1)); 
        e0=AE_ADD16(e0,e0);
        e1=AE_ADD16(e1,e1);
        z0=AE_MULFD16X16X4RAS(z0,z0,AE_MOVDA16(0x7fff),e0);
        z1=AE_MULFD16X16X4RAS(z1,z1,AE_MOVDA16(0x7fff),e1);
        /* restore original sign */
        y0=AE_NEG16S(z0);
        y1=AE_NEG16S(z1);
        AE_MOVT16X4(z0,y0,sy0);
        AE_MOVT16X4(z1,y1,sy1);
        /* multiply by X */
        z0=AE_MULFP16X4RAS(x0,z0);
        z1=AE_MULFP16X4RAS(x1,z1);
        AE_SA16X4X2_IP(z0,z1,aF,pF);
    }
    AE_SA128POS_FP(aF,pF);
    AE_SA128POS_FP(aE,pE);
} /* vec_divide16x16() */
