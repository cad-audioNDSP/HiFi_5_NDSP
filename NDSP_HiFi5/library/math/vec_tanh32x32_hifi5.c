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
  NatureDSP Signal Processing Library. Vector matematics
    Hyperbolic Tangent
    Optimized code for HiFi5
  IntegrIT, 2006-2019
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Signal_math.h"
#include "common.h"
#include "vec_tail32x32.h"

#define SMALLER_CODESIZE 1 // 1 - simpler variant, 2 - better codesize, but some performance loss

/*-------------------------------------------------------------------------
  Hyperbolic Tangent
  The functions compute the hyperbolic tangent of input argument. 32-bit
  fixed-point functions accept inputs in Q6.25 and form outputs in Q16.15
  format.

  Precision:
  32x32  32-bit inputs, 32-bit output. Accuracy: 2 LSB.
  f      single precision floating-point. Accuracy 2 ULP
  fp16   half precision floating-point. Accuracy 2 ULP
  Input:
  x[N]   input data, Q6.25 or floating point  
  N      length of vectors
  Output:
  y[N]   result, Q16.15 or floating point

  Restriction:
  x,y should not overlap

  Scalar versions:
  ----------------
  return result, Q16.15 or floating point
-------------------------------------------------------------------------*/
void vec_tanh32x32(int32_t* restrict y, const int32_t* restrict x, int N)
#if SMALLER_CODESIZE==1
{
  int n;
  static const int32_t ALIGN(32) polypow2[] = { 14685057, -114217091, 514075394, -1488269031, 2147475316 };

        ae_int32x4 * restrict pY = (      ae_int32x4 *)y;
  const ae_int32x4 * restrict pX = (const ae_int32x4 *)x;
  const ae_int32x4 * restrict pX1 = (const ae_int32x4 *)x;

  ae_int32x2 y0, y1, t0, t1, z0, z1;
  ae_int32x2 x0, x1, e0, e1, d0, d1;
  ae_valignx2 aX, aX1, aY;
  xtbool2 sign0, sign1;

  if (N<=0) return;
  if (N >= 8)
  {
		aX = AE_LA128_PP(pX);
		aX1 = AE_LA128_PP(pX1);

		aY = AE_ZALIGN128();
	  __Pragma("loop_count min=2,factor=2")
		  for (n = 0; n < ((N&~7) >> 2); n++)
		  {
			  AE_LA32X2X2_IP(x0, x1, aX, pX);
			  sign0 = AE_LT32(x0, 0);
			  sign1 = AE_LT32(x1, 0);
			  AE_MULF2P32X4RAS(z0, z1, x0, x1, 1549082005, 1549082005);
			  x0 = AE_ABS32S(z0);
			  x1 = AE_ABS32S(z1);
			  e0 = AE_SRAI32(x0, 23);
			  e1 = AE_SRAI32(x1, 23);
			  x0 = AE_MOVDEXT(x0, 23, 8);
			  x1 = AE_MOVDEXT(x1, 23, 8);
			  y0 = y1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 0);
			  t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 1);
			  AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
			  t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 2);
			  AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
			  t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 3);
			  AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
			  t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 4);
			  AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
			  x0 = AE_SRAV32RS(y0, e0);//Q31
			  x1 = AE_SRAV32RS(y1, e1);

			  /* 0.96-x/2 */
			  z0 = AE_SUB32(1030792151, AE_SRAI32(x0, 2));//Q30
			  z1 = AE_SUB32(1030792151, AE_SRAI32(x1, 2));
			  t0 = z0; t1 = z1;

			  AE_MULAF2P32X4RAS(t0, t1, z0, z1, x0, x1);  //Q30+Q(30+31-31)=Q30
			  d0 = AE_SUB32(1073741824, t0);//Q30
			  d1 = AE_SUB32(1073741824, t1);
			  t0 = AE_SRAI32(z0, 1); t1 = AE_SRAI32(z1, 1);//Q29
			  AE_MULAF2P32X4RAS(t0, t1, z0, z1, d0, d1); z0 = t0; z1 = t1; //Q29+Q(30+30-31)=Q29

			  AE_MULAF2P32X4RAS(t0, t1, z0, z1, x0, x1);//Q29+Q(29+31-31)
			  d0 = AE_SUB32(536870912, t0);//Q29
			  d1 = AE_SUB32(536870912, t1);
			  //    t0 = AE_SRAI32(z0, 2); t1 = AE_SRAI32(z1, 2);//Q27
			  AE_MULF2P32X4RAS(t0, t1, z0, z1, 0x20000000, 0x20000000);//Q27
			  AE_MULAF2P32X4RAS(t0, t1, z0, z1, d0, d1); z0 = t0; z1 = t1;//Q27+Q(29+29-31)=Q27

			  x0 = AE_SRAV32RS(x0, 12);//Q19
			  x1 = AE_SRAV32RS(x1, 12);
			  y0 = AE_SUB32(524288, x0);
			  y1 = AE_SUB32(524288, x1);
			  AE_MULF2P32X4RAS(z0, z1, z0, z1, y0, y1);//Q(27+19-31)=15

			  AE_LA32X2X2_IP(x0, x1, aX1, pX1);
			  sign0 = AE_LT32(x0, 0);
			  sign1 = AE_LT32(x1, 0);
			  z0 = AE_MOVNEG32S_T(z0, x0);
			  z1 = AE_MOVNEG32S_T(z1, x1);
			  AE_SA32X2X2_IP(z0, z1, aY, pY);
		  }
	  AE_SA128POS_FP(aY, pY);
  }
  N&=7;
  if (N) vec_tail32x32((int32_t*)pY,(const int32_t*)pX,vec_tanh32x32,N);
}
#elif SMALLER_CODESIZE==2 && defined (AE_LAV16X4X2_XP) && defined (AE_SAV16X4X2_XP)
{
  int n;
  static const int32_t ALIGN(32) polypow2[] = { 14685057, -114217091, 514075394, -1488269031, 2147475316 };
  ae_int16x4 tmp0,tmp1;
  static const int16_t ALIGN(16) change1632_tbl[]={6|(2<<8), 7|(3<<8), 4|(0<<8), 5|(1<<8)};
  ae_int16x4 change1632;

        ae_int32x4 * restrict pY = (      ae_int32x4 *)y;
  const ae_int32x4 * restrict pX = (const ae_int32x4 *)x;
  const ae_int32x4 * restrict pX1 = (const ae_int32x4 *)x;
  int nbytesRd,nbytesRd1,nbytesWr;
  ae_int32x2 y0, y1, t0, t1, z0, z1;
  ae_int32x2 x0, x1, e0, e1, d0, d1;

  ae_int32x2 y2, y3, t2, t3, z2, z3;
  ae_int32x2 x2, x3, e2, e3, d2, d3;
  ae_valignx2 aX, aX1, aY;
  xtbool2 sign0, sign1;
  xtbool2 sign2, sign3;
  if (N<=0) return;

  aX = AE_LA128_PP(pX);
  aX1 = AE_LA128_PP(pX1);
  aY = AE_ZALIGN128();
  change1632=AE_L16X4_I((const ae_int16x4*)change1632_tbl,0);

  nbytesRd=nbytesRd1=nbytesWr=N<<2;
  for (n=0; n<((N+7)>>3); n++)
  {
    AE_LAV16X4X2_XP(tmp0, tmp1,aX,castxcc(ae_int16x8,pX), nbytesRd); nbytesRd-=16;
    AE_DSEL16X4(tmp0,tmp1,tmp1,tmp0,change1632);
    x0=AE_MOVINT32X2_FROMINT16X4(tmp0);
    x1=AE_MOVINT32X2_FROMINT16X4(tmp1);
    AE_LAV16X4X2_XP(tmp0, tmp1,aX,castxcc(ae_int16x8,pX), nbytesRd); nbytesRd-=16;
    AE_DSEL16X4(tmp0,tmp1,tmp1,tmp0,change1632);
    x2=AE_MOVINT32X2_FROMINT16X4(tmp0);
    x3=AE_MOVINT32X2_FROMINT16X4(tmp1);

    sign0 = AE_LT32(x0, 0);
    sign1 = AE_LT32(x1, 0);
    sign2 = AE_LT32(x2, 0);
    sign3 = AE_LT32(x3, 0);
    AE_MULF2P32X4RAS(z0, z1, x0, x1, 1549082005, 1549082005);
    AE_MULF2P32X4RAS(z2, z3, x2, x3, 1549082005, 1549082005);
    x0 = AE_ABS32S(z0);
    x1 = AE_ABS32S(z1);
    x2 = AE_ABS32S(z2);
    x3 = AE_ABS32S(z3);
    e0 = AE_SRAI32(x0, 23);
    e1 = AE_SRAI32(x1, 23);
    e2 = AE_SRAI32(x2, 23);
    e3 = AE_SRAI32(x3, 23);
    x0 = AE_MOVDEXT(x0, 23,8);
    x1 = AE_MOVDEXT(x1, 23,8);
    x2 = AE_MOVDEXT(x2, 23,8);
    x3 = AE_MOVDEXT(x3, 23,8);
    y2 = y3 = y0 = y1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 0);
    t2 = t3 = t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 1);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    AE_MULAF2P32X4RAS(t2, t3, x2, x3, y2, y3); y2 = t2; y3 = t3;
    t2 = t3 = t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 2);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    AE_MULAF2P32X4RAS(t2, t3, x2, x3, y2, y3); y2 = t2; y3 = t3;
    t2 = t3 = t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 3);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    AE_MULAF2P32X4RAS(t2, t3, x2, x3, y2, y3); y2 = t2; y3 = t3;
    t2 = t3 =t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 4);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    AE_MULAF2P32X4RAS(t2, t3, x2, x3, y2, y3); y2 = t2; y3 = t3;
    x0 = AE_SRAV32RS(y0, e0);//Q31
    x1 = AE_SRAV32RS(y1, e1);
    x2 = AE_SRAV32RS(y2, e2);//Q31
    x3 = AE_SRAV32RS(y3, e3);

    /* 0.96-x/2 */
    z0 = AE_SUB32(1030792151, AE_SRAI32(x0, 2));//Q30
    z1 = AE_SUB32(1030792151, AE_SRAI32(x1, 2));
    z2 = AE_SUB32(1030792151, AE_SRAI32(x2, 2));//Q30
    z3 = AE_SUB32(1030792151, AE_SRAI32(x3, 2));
    t0 = z0; t1 = z1;
    t2 = z2; t3 = z3;

    AE_MULAF2P32X4RAS(t0, t1, z0, z1, x0, x1);  //Q30+Q(30+31-31)=Q30
    AE_MULAF2P32X4RAS(t2, t3, z2, z3, x2, x3);  //Q30+Q(30+31-31)=Q30
    d0 = AE_SUB32(1073741824, t0);//Q30
    d1 = AE_SUB32(1073741824, t1);
    d2 = AE_SUB32(1073741824, t2);//Q30
    d3 = AE_SUB32(1073741824, t3);
    t0 = AE_SRAI32(z0, 1); t1 = AE_SRAI32(z1, 1);//Q29
    t2 = AE_SRAI32(z2, 1); t3 = AE_SRAI32(z3, 1);//Q29
    AE_MULAF2P32X4RAS(t0, t1, z0, z1, d0, d1);z0 = t0; z1 = t1; //Q29+Q(30+30-31)=Q29
    AE_MULAF2P32X4RAS(t2, t3, z2, z3, d2, d3);z2 = t2; z3 = t3; //Q29+Q(30+30-31)=Q29

    AE_MULAF2P32X4RAS(t0, t1, z0, z1, x0, x1);//Q29+Q(29+31-31)
    AE_MULAF2P32X4RAS(t2, t3, z2, z3, x2, x3);//Q29+Q(29+31-31)
    d0 = AE_SUB32(536870912, t0);//Q29
    d1 = AE_SUB32(536870912, t1);
    d2 = AE_SUB32(536870912, t2);//Q29
    d3 = AE_SUB32(536870912, t3);
    AE_MULF2P32X4RAS(t0,t1,z0,z1,0x20000000,0x20000000);//Q27
    AE_MULF2P32X4RAS(t2,t3,z2,z3,0x20000000,0x20000000);//Q27
    AE_MULAF2P32X4RAS(t0, t1, z0, z1, d0, d1); z0 = t0; z1 = t1;//Q27+Q(29+29-31)=Q27
    AE_MULAF2P32X4RAS(t2, t3, z2, z3, d2, d3); z2 = t2; z3 = t3;//Q27+Q(29+29-31)=Q27

    x0 = AE_SRAV32RS(x0, 12);//Q19
    x1 = AE_SRAV32RS(x1, 12);
    x2 = AE_SRAV32RS(x2, 12);//Q19
    x3 = AE_SRAV32RS(x3, 12);
    y0 = AE_SUB32(524288, x0);
    y1 = AE_SUB32(524288, x1);
    y2 = AE_SUB32(524288, x2);
    y3 = AE_SUB32(524288, x3);
    AE_MULF2P32X4RAS(z0, z1, z0, z1, y0, y1);//Q(27+19-31)=15
    AE_MULF2P32X4RAS(z2, z3, z2, z3, y2, y3);//Q(27+19-31)=15
#if 1
    AE_LAV16X4X2_XP(tmp0, tmp1,aX1,castxcc(ae_int16x8,pX1), nbytesRd1); nbytesRd1-=16;
    AE_DSEL16X4(tmp0,tmp1,tmp1,tmp0,change1632);
    x0=AE_MOVINT32X2_FROMINT16X4(tmp0);
    x1=AE_MOVINT32X2_FROMINT16X4(tmp1);

    AE_LAV16X4X2_XP(tmp0, tmp1,aX1,castxcc(ae_int16x8,pX1), nbytesRd1); nbytesRd1-=16;
    AE_DSEL16X4(tmp0,tmp1,tmp1,tmp0,change1632);
    x2=AE_MOVINT32X2_FROMINT16X4(tmp0);
    x3=AE_MOVINT32X2_FROMINT16X4(tmp1);

    sign0 = AE_LT32(x0, 0);
    sign1 = AE_LT32(x1, 0);
    sign2 = AE_LT32(x2, 0);
    sign3 = AE_LT32(x3, 0);
    z0=AE_MOVNEG32S_T(z0,x0);
    z1=AE_MOVNEG32S_T(z1,x1);
    z2=AE_MOVNEG32S_T(z2,x2);
    z3=AE_MOVNEG32S_T(z3,x3);
#else
    x0 = AE_NEG32S(z0);
    x1 = AE_NEG32S(z1);
    x2 = AE_NEG32S(z2);
    x3 = AE_NEG32S(z3);

    AE_MOVT32X2(z0, x0, sign0);
    AE_MOVT32X2(z1, x1, sign1);
    AE_MOVT32X2(z2, x2, sign2);
    AE_MOVT32X2(z3, x3, sign3);

#endif

    tmp0=AE_MOVINT16X4_FROMINT32X2(z0);
    tmp1=AE_MOVINT16X4_FROMINT32X2(z1);
    AE_DSEL16X4(tmp0,tmp1,tmp1,tmp0,change1632);
    AE_SAV16X4X2_XP(tmp0, tmp1, aY,castxcc(ae_int16x8,pY), nbytesWr); nbytesWr-=16;
    tmp0=AE_MOVINT16X4_FROMINT32X2(z2);
    tmp1=AE_MOVINT16X4_FROMINT32X2(z3);
    AE_DSEL16X4(tmp0,tmp1,tmp1,tmp0,change1632);
    AE_SAV16X4X2_XP(tmp0, tmp1, aY,castxcc(ae_int16x8,pY), nbytesWr); nbytesWr-=16;

  }
  AE_SA128POS_FP(aY, pY);
} /* vec_tanh32x32() */
#else
{
  int n;
  static const int32_t ALIGN(32) polypow2[] = { 14685057, -114217091, 514075394, -1488269031, 2147475316 };

        ae_int32x4 * restrict pY = (      ae_int32x4 *)y;
  const ae_int32x4 * restrict pX = (const ae_int32x4 *)x;
  const ae_int32x4 * restrict pX1 = (const ae_int32x4 *)x;

  ae_int32x2 y0, y1, t0, t1, z0, z1;
  ae_int32x2 x0, x1, e0, e1, d0, d1;
  ae_valignx2 aX, aX1, aY;
  xtbool2 sign0, sign1;
  aX = AE_LA128_PP(pX);
  aX1 = AE_LA128_PP(pX1);

  aY = AE_ZALIGN128();

  if (N<=0) return;

  for (n=0; n<(N>>2); n++)
  {
    AE_LA32X2X2_IP(x0, x1, aX, pX);
    sign0 = AE_LT32(x0, 0);
    sign1 = AE_LT32(x1, 0);
    AE_MULF2P32X4RAS(z0, z1, x0, x1, 1549082005, 1549082005);
    x0 = AE_ABS32S(z0);
    x1 = AE_ABS32S(z1);
    e0 = AE_SRAI32(x0, 23);
    e1 = AE_SRAI32(x1, 23);
#if 0
    x0 = AE_AND32(x0, AE_MOVDA32X2(0x007fffff, 0x007fffff));
    x1 = AE_AND32(x1, AE_MOVDA32X2(0x007fffff, 0x007fffff));
    x0 = AE_SLAI32(x0, 8);//Q31
    x1 = AE_SLAI32(x1, 8);
#else
    x0 = AE_MOVDEXT(x0, 23,8);
    x1 = AE_MOVDEXT(x1, 23,8);
#endif
    y0 = y1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 0);
    t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 1);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 2);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 3);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 4);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    x0 = AE_SRAV32RS(y0, e0);//Q31
    x1 = AE_SRAV32RS(y1, e1);

    /* 0.96-x/2 */
    z0 = AE_SUB32(1030792151, AE_SRAI32(x0, 2));//Q30
    z1 = AE_SUB32(1030792151, AE_SRAI32(x1, 2));
    t0 = z0; t1 = z1;

    AE_MULAF2P32X4RAS(t0, t1, z0, z1, x0, x1);  //Q30+Q(30+31-31)=Q30
    d0 = AE_SUB32(1073741824, t0);//Q30
    d1 = AE_SUB32(1073741824, t1);
    t0 = AE_SRAI32(z0, 1); t1 = AE_SRAI32(z1, 1);//Q29
    AE_MULAF2P32X4RAS(t0, t1, z0, z1, d0, d1);z0 = t0; z1 = t1; //Q29+Q(30+30-31)=Q29

    AE_MULAF2P32X4RAS(t0, t1, z0, z1, x0, x1);//Q29+Q(29+31-31)
    d0 = AE_SUB32(536870912, t0);//Q29
    d1 = AE_SUB32(536870912, t1);
//    t0 = AE_SRAI32(z0, 2); t1 = AE_SRAI32(z1, 2);//Q27
    AE_MULF2P32X4RAS(t0,t1,z0,z1,0x20000000,0x20000000);//Q27
    AE_MULAF2P32X4RAS(t0, t1, z0, z1, d0, d1); z0 = t0; z1 = t1;//Q27+Q(29+29-31)=Q27

    x0 = AE_SRAV32RS(x0, 12);//Q19
    x1 = AE_SRAV32RS(x1, 12);
    y0 = AE_SUB32(524288, x0);
    y1 = AE_SUB32(524288, x1);
    AE_MULF2P32X4RAS(z0, z1, z0, z1, y0, y1);//Q(27+19-31)=15

    AE_LA32X2X2_IP(x0, x1, aX1, pX1);
    sign0 = AE_LT32(x0, 0);
    sign1 = AE_LT32(x1, 0);
#if 0
    x0 = AE_NEG32S(z0);
    x1 = AE_NEG32S(z1);

    AE_MOVT32X2(z0, x0, sign0);
    AE_MOVT32X2(z1, x1, sign1);
#else
    z0=AE_MOVNEG32S_T(z0,x0);
    z1=AE_MOVNEG32S_T(z1,x1);
#endif
    AE_SA32X2X2_IP(z0, z1, aY, pY);
  }
  AE_SA128POS_FP(aY, pY);
  x += (N&~3);
  y += (N&~3);
  N &= 3;
  if (N>0)
  {
    int32_t ALIGN(32) scratch[4];
    ae_int32x4 *pScr;
    pScr = (ae_int32x4*)scratch;
    pX = (const ae_int32x4*)x;
    pY = (      ae_int32x4*)y;
    AE_S32X2X2_I(0, 0, pScr, 0);
    __Pragma("no_unroll")
    for (n = 0; n<N; n++)
    {
      ae_int32x2 t;
      AE_L32_IP(t, castxcc(ae_int32, pX), sizeof(int32_t));
      AE_S32_L_IP(t, castxcc(ae_int32, pScr), sizeof(int32_t));
    }
    pScr = (ae_int32x4*)scratch;
    AE_L32X2X2_I(x0, x1, pScr, 0 * sizeof(ae_int32x4));
    sign0 = AE_LT32(x0, 0);
    sign1 = AE_LT32(x1, 0);
    AE_MULF2P32X4RAS(z0, z1, x0, x1, 1549082005, 1549082005);
    x0 = AE_ABS32S(z0);
    x1 = AE_ABS32S(z1);
    e0 = AE_SRAI32(x0, 23);
    e1 = AE_SRAI32(x1, 23);
#if 0
    x0 = AE_AND32(x0, AE_MOVDA32X2(0x007fffff, 0x007fffff));
    x1 = AE_AND32(x1, AE_MOVDA32X2(0x007fffff, 0x007fffff));
    x0 = AE_SLAI32(x0, 8);//Q31
    x1 = AE_SLAI32(x1, 8);
#else
    x0 = AE_MOVDEXT(x0, 23,8);
    x1 = AE_MOVDEXT(x1, 23,8);
#endif
    y0 = y1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 0);
    t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 1);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 2);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 3);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = AE_L32_I((const ae_int32 *)polypow2, 4 * 4);
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    x0 = AE_SRAV32RS(y0, e0);//Q31
    x1 = AE_SRAV32RS(y1, e1);

    /* 0.96-x/2 */
    z0 = AE_SUB32(1030792151, AE_SRAI32(x0, 2));//Q30
    z1 = AE_SUB32(1030792151, AE_SRAI32(x1, 2));
    t0 = z0; t1 = z1;

    AE_MULAF2P32X4RAS(t0, t1, z0, z1, x0, x1);  //Q30+Q(30+31-31)=Q30
    d0 = AE_SUB32(1073741824, t0);//Q30
    d1 = AE_SUB32(1073741824, t1);
    t0 = AE_SRAI32(z0, 1); t1 = AE_SRAI32(z1, 1);//Q29
    AE_MULAF2P32X4RAS(t0, t1, z0, z1, d0, d1); z0 = t0; z1 = t1; //Q29+Q(30+30-31)=Q29

    AE_MULAF2P32X4RAS(t0, t1, z0, z1, x0, x1);//Q29+Q(29+31-31)
    d0 = AE_SUB32(536870912, t0);//Q29
    d1 = AE_SUB32(536870912, t1);
    t0 = AE_SRAI32(z0, 2); t1 = AE_SRAI32(z1, 2);//Q27
    AE_MULAF2P32X4RAS(t0, t1, z0, z1, d0, d1); z0 = t0; z1 = t1;//Q27+Q(29+29-31)=Q27

    x0 = AE_SRAV32RS(x0, 12);//Q19
    x1 = AE_SRAV32RS(x1, 12);
    y0 = AE_SUB32(524288, x0);
    y1 = AE_SUB32(524288, x1);
    AE_MULF2P32X4RAS(z0, z1, z0, z1, y0, y1);//Q(27+19-31)=15

    AE_LA32X2X2_IP(x0, x1, aX1, pX1);
    sign0 = AE_LT32(x0, 0);
    sign1 = AE_LT32(x1, 0);
    x0 = AE_NEG32S(z0);
    x1 = AE_NEG32S(z1);

    AE_MOVT32X2(z0, x0, sign0);
    AE_MOVT32X2(z1, x1, sign1);
    AE_S32X2X2_I(z0, z1, pScr, 0 * sizeof(ae_int32x4));
    __Pragma("no_unroll")
    for (n = 0; n<N; n++)
    {
      ae_int32x2 t;
      AE_L32_IP(t, castxcc(ae_int32, pScr), sizeof(int32_t));
      AE_S32_L_IP(t, castxcc(ae_int32, pY), sizeof(int32_t));
    }
  }
} /* vec_tanh32x32() */
#endif
