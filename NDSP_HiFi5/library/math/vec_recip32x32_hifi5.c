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
#include "NatureDSP_Signal_math.h"
#include "common.h"

#define SMALLER_CODESIZE 1 // 1 no performance loss, 2 - some performance loss, but less code size even more

/*===========================================================================
  Vector matematics:
  vec_recip            Reciprocal on Q31/Q15 Numbers
===========================================================================*/

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
void vec_recip32x32
(
  int32_t * restrict        frac,
  int16_t *                 exp,
  const int32_t * restrict  x,
  int                       N
)
#if SMALLER_CODESIZE==1
{
  int n;
  ae_int32x4 * restrict pF = (ae_int32x4 *)frac;
  ae_int16x4 * restrict pE = (      ae_int16x4 *)exp;
  const ae_int32x4 * restrict pX = (const ae_int32x4 *)x;

  ae_valignx2 aX, aF;
  ae_valign aE;
  if (N <= 0) return;
  if (N >= 8)
  {
	  aX = AE_LA128_PP(pX);
	  aE = AE_ZALIGN64();
	  aF = AE_ZALIGN128();

	  /* compute exponent and normalize inputs */
	  pX = (const ae_int32x4 *)x;
	  pF = (ae_int32x4*)frac;
	  __Pragma("loop_count min=2, factor=2")
	  for (n = 0; n < ((N&~7) >> 2); n++)
	  {
		  xtbool2 isZero0, isZero1;
		  xtbool4 isZero;
		  ae_int32x2 x0, x1, t0, t1, y0, y1, e0, e1;
		  ae_int32x2 sx0, sx1;
		  ae_int32x2 _0x40000000;
		  ae_int16x4 nsay, expxy, t;
		  _0x40000000 = AE_MOVDA32(0x40000000);
		  AE_LA32X2X2_IP(y0, y1, aX, pX);
		  sx0 = y0; sx1 = y1;
		  nsay = AE_NSA32X4(y0, y1);
		  y0 = AE_SRAV32RS(y0, AE_SEXT32X2D16_32(AE_NEG16S(nsay)));
		  y1 = AE_SRAV32RS(y1, AE_SEXT32X2D16_10(AE_NEG16S(nsay)));

		  expxy = AE_ADD16(nsay, 1);
		  x0 = AE_ABS32S(y0);
		  x1 = AE_ABS32S(y1);
		  t = AE_TRUNC16X4F32(x0, x1);
		  isZero = AE_EQ16(t, AE_ZERO16());
		  isZero0 = AE_EQ32(x0, AE_ZERO32());
		  isZero1 = AE_EQ32(x1, AE_ZERO32());

		  AE_MOVT32X2(x0, AE_MOVDA32(0x7fffffff), isZero0);
		  AE_MOVT32X2(x1, AE_MOVDA32(0x7fffffff), isZero1);

		  AE_MOVT16X4(expxy, 0x20, isZero);

		  AE_SA16X4_IP(expxy, aE, pE);

		  y0 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), x0);
		  y1 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), x1);
		  /* 4 iterations to achieve 1 LSB accuracy in mantissa */
		  t0 = t1 = _0x40000000;
		  AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
		  e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
		  t0 = y0; t1 = y1;
		  AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
		  t0 = t1 = _0x40000000;
		  AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
		  e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
		  t0 = y0; t1 = y1;
		  AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
		  t0 = t1 = _0x40000000;
		  AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
		  e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
		  t0 = y0; t1 = y1;
		  AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
		  t0 = t1 = _0x40000000;
		  AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
		  e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
		  t0 = y0; t1 = y1;
		  AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;

		  /* restore original sign */
		  y0 = AE_MOVNEG32S_T(y0, sx0);
		  y1 = AE_MOVNEG32S_T(y1, sx1);

		  AE_SA32X2X2_IP(y0, y1, aF, pF);
	  }
	  AE_SA128POS_FP(aF, pF);
	  AE_SA64POS_FP(aE, pE);
  }
  N&=7;
  if (N)
  {
       int32_t ALIGN(16) buff[8];
       int16_t ALIGN(16) ee[8];
       int n;
       AE_S32X2X2_I(0,0,(ae_int32x4*)buff,0);
       AE_S32X2X2_I(0,0,(ae_int32x4*)buff,sizeof(ae_int32x4));
       __Pragma("loop_count min=1,max=7")
       __Pragma("no_unroll")
       __Pragma("no_simd")
       for (n=0; n<N; n++) buff[n]=((const int32_t*)pX)[n];
       vec_recip32x32(buff,ee,buff,8);
       __Pragma("loop_count min=1,max=7")
       __Pragma("no_unroll")
	   __Pragma("no_simd")
		   for (n=0; n<N; n++)
       {
           ((int32_t*)pF)[n]=buff[n];
           ((int16_t*)pE)[n]=ee[n];
       }
  }
}
#elif SMALLER_CODESIZE==2 && defined (AE_LAV16X4X2_XP) && defined (AE_SAV16X4X2_XP)
{
    ae_int16x4 tmp0,tmp1;
  static const int16_t ALIGN(16) change1632_tbl[]={6|(2<<8), 7|(3<<8), 4|(0<<8), 5|(1<<8)};
  ae_int16x4 change1632;
  int n,nbytesRd,nbytesWr,nbytesE;
  ae_int32x4 * restrict pF = (ae_int32x4 *)frac;
  ae_int16x4 * restrict pE = (      ae_int16x4 *)exp;
  const ae_int32x4 * restrict pX = (const ae_int32x4 *)x;

  ae_valignx2 aX, aF;
  ae_valignx2 aE;
  if (N<=0) return;
  aX = AE_LA128_PP(pX);
  aE = AE_ZALIGN128();
  aF = AE_ZALIGN128();

  change1632=AE_L16X4_I((const ae_int16x4*)change1632_tbl,0);
  /* compute exponent and normalize inputs */
  pX = (const ae_int32x4 *)x;
  pF = (      ae_int32x4*)frac;
  nbytesWr = nbytesRd=N<<2;
  nbytesE  = N<<1;
  for (n = 0; n<((N+7) >> 3); n++)
  {
    ae_int32x2 x0, x1, t0, t1, y0, y1, e0, e1, sx0, sx1;
    ae_int32x2 x2, x3, t2, t3, y2, y3, e2, e3, sx2, sx3;
    ae_int32x2 _0x40000000;
    ae_int16x4 nsay0, expxy0;
    ae_int16x4 nsay1, expxy1;
    _0x40000000 = AE_MOVDA32(0x40000000); 
    AE_LAV16X4X2_XP(tmp0, tmp1,aX,castxcc(ae_int16x8,pX), nbytesRd); nbytesRd-=16;
    AE_DSEL16X4(tmp0,tmp1,tmp1,tmp0,change1632);
    y0=AE_MOVINT32X2_FROMINT16X4(tmp0);
    y1=AE_MOVINT32X2_FROMINT16X4(tmp1);
    AE_LAV16X4X2_XP(tmp0, tmp1,aX,castxcc(ae_int16x8,pX), nbytesRd); nbytesRd-=16;
    AE_DSEL16X4(tmp0,tmp1,tmp1,tmp0,change1632);
    y2=AE_MOVINT32X2_FROMINT16X4(tmp0);
    y3=AE_MOVINT32X2_FROMINT16X4(tmp1);

    sx0=y0; sx1=y1;
    nsay0 = AE_NSA32X4(y0, y1);
    y0 = AE_SRAV32RS(y0, AE_SEXT32X2D16_32(AE_NEG16S(nsay0)));
    y1 = AE_SRAV32RS(y1, AE_SEXT32X2D16_10(AE_NEG16S(nsay0)));
    sx2=y2; sx3=y3;
    nsay1 = AE_NSA32X4(y2, y3);
    y2 = AE_SRAV32RS(y2, AE_SEXT32X2D16_32(AE_NEG16S(nsay1)));
    y3 = AE_SRAV32RS(y3, AE_SEXT32X2D16_10(AE_NEG16S(nsay1)));

    expxy0 = AE_ADD16(nsay0, 1);    
    expxy1 = AE_ADD16(nsay1, 1);    
    x0 = AE_ABS32S(y0);
    x1 = AE_ABS32S(y1);
    x2 = AE_ABS32S(y2);
    x3 = AE_ABS32S(y3);
    AE_MOVT16X4(expxy0, 0x20, AE_EQ16(AE_TRUNC16X4F32(x0, x1), AE_ZERO16()));
    AE_MOVT16X4(expxy1, 0x20, AE_EQ16(AE_TRUNC16X4F32(x2, x3), AE_ZERO16()));
    AE_MOVT32X2(x0, AE_MOVDA32(0x7fffffff), AE_EQ32(x0, AE_ZERO32()));
    AE_MOVT32X2(x1, AE_MOVDA32(0x7fffffff), AE_EQ32(x1, AE_ZERO32()));
    AE_MOVT32X2(x2, AE_MOVDA32(0x7fffffff), AE_EQ32(x2, AE_ZERO32()));
    AE_MOVT32X2(x3, AE_MOVDA32(0x7fffffff), AE_EQ32(x3, AE_ZERO32()));

    AE_SAV16X4X2_XP(expxy0,expxy1,  aE, castxcc(ae_int16x8,pE),nbytesE); nbytesE-=16;

    y0 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), x0);
    y1 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), x1);
    y2 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), x2);
    y3 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), x3);
    /* 4 iterations to achieve 1 LSB accuracy in mantissa */
    t0 = t1 = _0x40000000;
    t2 = t3 = _0x40000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    AE_MULSF2P32X4RAS(t2, t3, x2, x3, y2, y3); e2 = t2; e3 = t3;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    e2 = AE_ADD32(e2, e2); e3 = AE_ADD32(e3, e3);
    t0 = y0; t1 = y1;
    t2 = y2; t3 = y3;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
    AE_MULAF2P32X4RAS(t2, t3, e2, e3, y2, y3); y2 = t2; y3 = t3;
    t0 = t1 = _0x40000000;
    t2 = t3 = _0x40000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    AE_MULSF2P32X4RAS(t2, t3, x2, x3, y2, y3); e2 = t2; e3 = t3;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    e2 = AE_ADD32(e2, e2); e3 = AE_ADD32(e3, e3);
    t0 = y0; t1 = y1;
    t2 = y2; t3 = y3;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
    AE_MULAF2P32X4RAS(t2, t3, e2, e3, y2, y3); y2 = t2; y3 = t3;
    t0 = t1 = _0x40000000;
    t2 = t3 = _0x40000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    AE_MULSF2P32X4RAS(t2, t3, x2, x3, y2, y3); e2 = t2; e3 = t3;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    e2 = AE_ADD32(e2, e2); e3 = AE_ADD32(e3, e3);
    t0 = y0; t1 = y1;
    t2 = y2; t3 = y3;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
    AE_MULAF2P32X4RAS(t2, t3, e2, e3, y2, y3); y2 = t2; y3 = t3;
    t0 = t1 = _0x40000000;
    t2 = t3 = _0x40000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    AE_MULSF2P32X4RAS(t2, t3, x2, x3, y2, y3); e2 = t2; e3 = t3;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    e2 = AE_ADD32(e2, e2); e3 = AE_ADD32(e3, e3);
    t0 = y0; t1 = y1;
    t2 = y2; t3 = y3;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
    AE_MULAF2P32X4RAS(t2, t3, e2, e3, y2, y3); y2 = t2; y3 = t3;

    /* restore original sign */
    y0=AE_MOVNEG32S_T(y0,sx0);
    y1=AE_MOVNEG32S_T(y1,sx1);
    y2=AE_MOVNEG32S_T(y2,sx2);
    y3=AE_MOVNEG32S_T(y3,sx3);
    tmp0=AE_MOVINT16X4_FROMINT32X2(y0);
    tmp1=AE_MOVINT16X4_FROMINT32X2(y1);
    AE_DSEL16X4(tmp0,tmp1,tmp1,tmp0,change1632);
    AE_SAV16X4X2_XP(tmp0, tmp1, aF,castxcc(ae_int16x8,pF), nbytesWr); nbytesWr-=16;
    tmp0=AE_MOVINT16X4_FROMINT32X2(y2);
    tmp1=AE_MOVINT16X4_FROMINT32X2(y3);
    AE_DSEL16X4(tmp0,tmp1,tmp1,tmp0,change1632);
    AE_SAV16X4X2_XP(tmp0, tmp1, aF,castxcc(ae_int16x8,pF), nbytesWr); nbytesWr-=16;
  }
  AE_SA128POS_FP(aF, pF);
  AE_SA128POS_FP(aE, pE);
}
#else
{
  int n;
  ae_int32x4 * restrict pF = (ae_int32x4 *)frac;
  ae_int16x4 * restrict pE = (      ae_int16x4 *)exp;
  const ae_int32x4 * restrict pX = (const ae_int32x4 *)x;

  ae_valignx2 aX, aF;
  ae_valign aE;
  aX = AE_LA128_PP(pX);
  aE = AE_ZALIGN64();
  aF = AE_ZALIGN128();

  if (N<=0) return;
  /* compute exponent and normalize inputs */
  pX = (const ae_int32x4 *)x;
  pF = (      ae_int32x4*)frac;
  for (n = 0; n<(N >> 2); n++)
  {
    xtbool2 isZero0, isZero1;
    xtbool4 isZero;
    ae_int32x2 x0, x1, t0, t1, y0, y1, e0, e1;
    ae_int32x2 sx0, sx1;
    ae_int32x2 _0x40000000;
    ae_int16x4 nsay, expxy, t;
    _0x40000000 = AE_MOVDA32(0x40000000); 
    AE_LA32X2X2_IP(y0, y1, aX, pX);
    sx0=y0; sx1=y1;
    nsay = AE_NSA32X4(y0, y1);
    y0 = AE_SRAV32RS(y0, AE_SEXT32X2D16_32(AE_NEG16S(nsay)));
    y1 = AE_SRAV32RS(y1, AE_SEXT32X2D16_10(AE_NEG16S(nsay)));

    expxy = AE_ADD16(nsay, 1);    
    x0 = AE_ABS32S(y0);
    x1 = AE_ABS32S(y1);
    t = AE_TRUNC16X4F32(x0, x1);
    isZero = AE_EQ16(t, AE_ZERO16());
    isZero0 = AE_EQ32(x0, AE_ZERO32());
    isZero1 = AE_EQ32(x1, AE_ZERO32());

    AE_MOVT32X2(x0, AE_MOVDA32(0x7fffffff), isZero0);
    AE_MOVT32X2(x1, AE_MOVDA32(0x7fffffff), isZero1);

    AE_MOVT16X4(expxy, 0x20, isZero);

    AE_SA16X4_IP(expxy, aE, pE);

    y0 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), x0);
    y1 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), x1);
    /* 4 iterations to achieve 1 LSB accuracy in mantissa */
    t0 = t1 = _0x40000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    t0 = y0; t1 = y1;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = _0x40000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    t0 = y0; t1 = y1;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = _0x40000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    t0 = y0; t1 = y1;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = _0x40000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    t0 = y0; t1 = y1;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;

    /* restore original sign */
    y0=AE_MOVNEG32S_T(y0,sx0);
    y1=AE_MOVNEG32S_T(y1,sx1);

    AE_SA32X2X2_IP(y0, y1, aF, pF);
  }
  AE_SA128POS_FP(aF, pF);
  AE_SA64POS_FP(aE, pE);
  x += (N&~3);
  frac += (N&~3);
  exp += (N&~3);
  N &= 3;
  if (N>0)
  {
    xtbool2 isZero0, isZero1;
    xtbool4 isZero;
    ae_int32x2 sx0, sx1,x0, x1, t0, t1, y0, y1, e0, e1;
    ae_int32x2 _0x40000000;
    ae_int16x4 nsay, expxy, t;
    int32_t ALIGN(32) scratch[4+4];
    ae_int32x4 *pScr;
    _0x40000000 = AE_MOVDA32(0x40000000);
    pScr = (ae_int32x4*)scratch;
    pX = (const ae_int32x4*)x;
    pF = (      ae_int32x4*)frac;
    pE = (      ae_int16x4*)exp;
    AE_S32X2X2_I(0, 0, pScr, 0);
    AE_S32X2X2_I(0, 0, pScr, sizeof(ae_int32x4));
    __Pragma("no_unroll")
    for (n = 0; n<N; n++)
    {
      ae_int32x2 t;
      AE_L32_IP(t, castxcc(ae_int32, pX), sizeof(int32_t));
      AE_S32_L_IP(t, castxcc(ae_int32, pScr), sizeof(int32_t));
    }
    pScr = (ae_int32x4*)scratch;
    AE_L32X2X2_I(y0, y1, pScr, 0 * sizeof(ae_int32x4));
    sx0=y0;sx1=y1;
    nsay = AE_NSA32X4(y0, y1);
    y0 = AE_SRAV32RS(y0, AE_SEXT32X2D16_32(AE_NEG16S(nsay)));
    y1 = AE_SRAV32RS(y1, AE_SEXT32X2D16_10(AE_NEG16S(nsay)));

    expxy = AE_ADD16(nsay, 1);

    x0 = AE_ABS32S(y0);
    x1 = AE_ABS32S(y1);
    t = AE_TRUNC16X4F32(x0, x1);
    isZero = AE_EQ16(t, AE_ZERO16());
    isZero0 = AE_EQ32(x0, AE_ZERO32());
    isZero1 = AE_EQ32(x1, AE_ZERO32());

    AE_MOVT32X2(x0, AE_MOVDA32(0x7fffffff), isZero0);
    AE_MOVT32X2(x1, AE_MOVDA32(0x7fffffff), isZero1);

    AE_MOVT16X4(expxy, 0x20, isZero);

    y0 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), x0);
    y1 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), x1);
    /* 4 iterations to achieve 1 LSB accuracy in mantissa */
    t0 = t1 = _0x40000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    t0 = y0; t1 = y1;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = _0x40000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    t0 = y0; t1 = y1;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = _0x40000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    t0 = y0; t1 = y1;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = _0x40000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, y0, y1); e0 = t0; e1 = t1;
    e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
    t0 = y0; t1 = y1;
    AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;

    /* restore original sign */

    y0=AE_MOVNEG32S_T(y0,sx0);
    y1=AE_MOVNEG32S_T(y1,sx1);

    AE_S32X2X2_I(y0, y1, pScr, 0 * sizeof(ae_int32x4));
    AE_S16X4_I(expxy,(ae_int16x4 *) pScr, 1 * sizeof(ae_int32x4));
    __Pragma("no_unroll")
    for (n = 0; n<N; n++)
    {
      ae_int32x2 t;
      AE_L32_IP(t, castxcc(ae_int32, pScr), sizeof(int32_t));
      AE_S32_L_IP(t, castxcc(ae_int32, pF), sizeof(int32_t));
    }
    pScr = (ae_int32x4*)scratch+1;
    __Pragma("no_unroll")
    for (n = 0; n<N; n++)
    {
      ae_int16x4 h;
      AE_L16_IP(h, castxcc(ae_int16, pScr), sizeof(int16_t));
      AE_S16_0_IP(h, castxcc(ae_int16, pE), sizeof(int16_t));
    }
  }
}
#endif
