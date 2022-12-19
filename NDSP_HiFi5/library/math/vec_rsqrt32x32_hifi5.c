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
    Reciprocal Square Root
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
void vec_rsqrt32x32(int32_t * frac,
                    int16_t * exp,
                    const int32_t * x,
                    int N)
{
  int n;
        ae_int32x4 * restrict pF = (      ae_int32x4 *)frac;
        ae_int16x4 * restrict pE = (      ae_int16x4 *)exp;
  const ae_int32x4 * restrict pX = (const ae_int32x4 *)x;

  ae_int32x2 x0, x1, y0, y1, t0, t1;
  ae_int32x2 e0, e1, z0, z1;
  ae_int16x4 nsax;
  ae_valignx2 aX, aF;
  ae_valign aE;
  ae_int32x2 _0x80000000, _0x20000000;
  aX = AE_LA128_PP(pX);
  aE = AE_ZALIGN64();
  aF = AE_ZALIGN128();

  if (N<=0) return;
  /* compute exponent and normalize inputs */
  pX = (const ae_int32x4 *)x;
  pF = (      ae_int32x4*)frac;
  _0x80000000 = AE_MOVDA32(0x80000000);
  _0x20000000 = AE_MOVDA32(0x20000000);

  for (n = 0; n<(N >> 2); n++)
  {
    xtbool2 iszero0, iszero1, sx0, sx1;
    ae_int64 U064_l, U064_h, U164_l, U164_h;
    ae_int32x2 u0, u1;
    xtbool4 sx;
    ae_int16x4 t;
    AE_LA32X2X2_IP(x0, x1, aX, pX);
    nsax = AE_NSA32X4(x0, x1);   

    x0 = AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(AE_AND16(nsax, ~1))));
    x1 = AE_SRAV32RS(x1, AE_SEXT32X2D16_10(AE_NEG16S(AE_AND16(nsax, ~1))));
    nsax = AE_SRAI16(nsax, 1);
    nsax = AE_ADD16(nsax, 1);
    t = AE_TRUNC16X4F32(x0, x1);
    sx = AE_LT16(t, AE_ZERO16());
 
    AE_MOVT16X4(nsax, 0x1f, sx);
    AE_SA16X4_IP(nsax, aE, pE);
    
    sx0 = AE_LT32(x0, AE_ZERO32()); 
    sx1 = AE_LT32(x1, AE_ZERO32()); 
    iszero0 = AE_EQ32(x0, AE_ZERO32());
    iszero1 = AE_EQ32(x1, AE_ZERO32());
    
    y0 = AE_SUB32(_0x80000000, AE_SRAI32(x0, 1)); /* Q30 */
    y1 = AE_SUB32(_0x80000000, AE_SRAI32(x1, 1)); /* Q30 */
    
    AE_MOVT32X2(y0, AE_MOVDA32(0x7fffffff), iszero0);
    AE_MOVT32X2(y1, AE_MOVDA32(0x7fffffff), iszero1);
    
    /*1-st iteration*/
    AE_MULF2P32X4RAS(z0, z1, y0, y1, y0, y1);
    t0 = t1 = _0x20000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, z0, z1); e0 = t0; e1 = t1;
    AE_MULF2P32X4RAS(z0, z1, y0, y1, e0, e1);
    AE_MULA2P32X4(y0, y1, AE_MOVDA32(2), AE_MOVDA32(2), z0, z1); 
    
    AE_MULF2P32X4RAS(z0, z1, y0, y1, y0, y1);
    t0 = t1 = _0x20000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, z0, z1); e0 = t0; e1 = t1;
    AE_MULF2P32X4RAS(z0, z1, y0, y1, e0, e1);
    AE_MULA2P32X4(y0, y1, AE_MOVDA32(2), AE_MOVDA32(2), z0, z1);
    
    AE_MULF2P32X4RAS(z0, z1, y0, y1, y0, y1);
    t0 = t1 = _0x20000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, z0, z1); e0 = t0; e1 = t1;
    AE_MULF2P32X4RAS(z0, z1, y0, y1, e0, e1);
    AE_MULA2P32X4(y0, y1, AE_MOVDA32(2), AE_MOVDA32(2), z0, z1);
 
    /*last iteration will be done in bit higher precision*/
    AE_MUL32X2S_HH_LL(U064_h, U064_l, y0, y0);
    AE_MUL32X2S_HH_LL(U164_h, U164_l, y1, y1);
    U064_l = AE_ADD64(U064_l, 0x20000000);
    U164_l = AE_ADD64(U164_l, 0x20000000);
    U064_h = AE_ADD64(U064_h, 0x20000000);
    U164_h = AE_ADD64(U164_h, 0x20000000);
    U064_l = AE_SLLI64(U064_l, 2);
    U164_l = AE_SLLI64(U164_l, 2);
    U064_h = AE_SLLI64(U064_h, 2);
    U164_h = AE_SLLI64(U164_h, 2);
    u0 = AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(U064_h), AE_MOVINT32X2_FROMINT64(U064_l));
    u1 = AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(U164_h), AE_MOVINT32X2_FROMINT64(U164_l));
    z0 = AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(U064_l), AE_MOVINT32X2_FROMINT64(U064_h));
    z1 = AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(U164_l), AE_MOVINT32X2_FROMINT64(U164_h));
    U064_l = AE_MUL32U_LL(x0, u0);
    U164_l = AE_MUL32U_LL(x1, u1);
    x0 = AE_SEL32_HH(x0, x0);
    x1 = AE_SEL32_HH(x1, x1);
 
    U064_h = AE_MUL32U_LL(x0, z0);
    U164_h = AE_MUL32U_LL(x1, z1);
 
    U064_l = AE_SRLI64(U064_l, 29);
    U164_l = AE_SRLI64(U164_l, 29);
    U064_h = AE_SRLI64(U064_h, 29);
    U164_h = AE_SRLI64(U164_h, 29);
    e0 = AE_SEL32_LL(AE_MOVINT32X2_FROMINT64(U064_h), AE_MOVINT32X2_FROMINT64(U064_l));
    e1 = AE_SEL32_LL(AE_MOVINT32X2_FROMINT64(U164_h), AE_MOVINT32X2_FROMINT64(U164_l));
 
    AE_MULF2P32X4RAS(z0, z1, y0, y1, e0, e1);
    z0 = AE_ADD32(z0, AE_MOVDA32(2));
    z1 = AE_ADD32(z1, AE_MOVDA32(2));
    y0 = AE_SUB32S(y0, AE_SRAI32(z0, 2));
    y1 = AE_SUB32S(y1, AE_SRAI32(z1, 2));
    AE_MOVT32X2(y0, AE_MOVDA32(0x80000000), sx0);
    AE_MOVT32X2(y1, AE_MOVDA32(0x80000000), sx1);
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
    xtbool2 iszero0, iszero1, sx0, sx1;
    ae_int64 U064_l, U064_h, U164_l, U164_h;
    ae_int32x2 u0, u1;
    xtbool4 sx;
    ae_int16x4 t;
    int32_t ALIGN(32) scratch[4 + 4];
    ae_int32x4 *pScr;
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
    AE_L32X2X2_I(x0, x1, pScr, 0 * sizeof(ae_int32x4));
    nsax = AE_NSA32X4(x0, x1);   

    x0 = AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(AE_AND16(nsax, ~1))));
    x1 = AE_SRAV32RS(x1, AE_SEXT32X2D16_10(AE_NEG16S(AE_AND16(nsax, ~1))));
    nsax = AE_SRAI16(nsax, 1);
    nsax = AE_ADD16(nsax, 1);
    t = AE_TRUNC16X4F32(x0, x1);
    sx = AE_LT16(t, AE_ZERO16());
 
    AE_MOVT16X4(nsax, 0x1f, sx);
   
    sx0 = AE_LT32(x0, AE_ZERO32()); 
    sx1 = AE_LT32(x1, AE_ZERO32()); 
    iszero0 = AE_EQ32(x0, AE_ZERO32());
    iszero1 = AE_EQ32(x1, AE_ZERO32());
    
    y0 = AE_SUB32(_0x80000000, AE_SRAI32(x0, 1)); /* Q30 */
    y1 = AE_SUB32(_0x80000000, AE_SRAI32(x1, 1)); /* Q30 */
    
    AE_MOVT32X2(y0, AE_MOVDA32(0x7fffffff), iszero0);
    AE_MOVT32X2(y1, AE_MOVDA32(0x7fffffff), iszero1);
    
    /*1-st iteration*/
    AE_MULF2P32X4RAS(z0, z1, y0, y1, y0, y1);
    t0 = t1 = _0x20000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, z0, z1); e0 = t0; e1 = t1;
    AE_MULF2P32X4RAS(z0, z1, y0, y1, e0, e1);
    AE_MULA2P32X4(y0, y1, AE_MOVDA32(2), AE_MOVDA32(2), z0, z1); 
    
    AE_MULF2P32X4RAS(z0, z1, y0, y1, y0, y1);
    t0 = t1 = _0x20000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, z0, z1); e0 = t0; e1 = t1;
    AE_MULF2P32X4RAS(z0, z1, y0, y1, e0, e1);
    AE_MULA2P32X4(y0, y1, AE_MOVDA32(2), AE_MOVDA32(2), z0, z1);
    
    AE_MULF2P32X4RAS(z0, z1, y0, y1, y0, y1);
    t0 = t1 = _0x20000000;
    AE_MULSF2P32X4RAS(t0, t1, x0, x1, z0, z1); e0 = t0; e1 = t1;
    AE_MULF2P32X4RAS(z0, z1, y0, y1, e0, e1);
    AE_MULA2P32X4(y0, y1, AE_MOVDA32(2), AE_MOVDA32(2), z0, z1);
 
    /*last iteration will be done in bit higher precision*/
    AE_MUL32X2S_HH_LL(U064_h, U064_l, y0, y0);
    AE_MUL32X2S_HH_LL(U164_h, U164_l, y1, y1);
    U064_l = AE_ADD64(U064_l, 0x20000000);
    U164_l = AE_ADD64(U164_l, 0x20000000);
    U064_h = AE_ADD64(U064_h, 0x20000000);
    U164_h = AE_ADD64(U164_h, 0x20000000);
    U064_l = AE_SLLI64(U064_l, 2);
    U164_l = AE_SLLI64(U164_l, 2);
    U064_h = AE_SLLI64(U064_h, 2);
    U164_h = AE_SLLI64(U164_h, 2);
    u0 = AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(U064_h), AE_MOVINT32X2_FROMINT64(U064_l));
    u1 = AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(U164_h), AE_MOVINT32X2_FROMINT64(U164_l));
    z0 = AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(U064_l), AE_MOVINT32X2_FROMINT64(U064_h));
    z1 = AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(U164_l), AE_MOVINT32X2_FROMINT64(U164_h));
    U064_l = AE_MUL32U_LL(x0, u0);
    U164_l = AE_MUL32U_LL(x1, u1);
    x0 = AE_SEL32_HH(x0, x0);
    x1 = AE_SEL32_HH(x1, x1);
 
    U064_h = AE_MUL32U_LL(x0, z0);
    U164_h = AE_MUL32U_LL(x1, z1);
 
    U064_l = AE_SRLI64(U064_l, 29);
    U164_l = AE_SRLI64(U164_l, 29);
    U064_h = AE_SRLI64(U064_h, 29);
    U164_h = AE_SRLI64(U164_h, 29);
    e0 = AE_SEL32_LL(AE_MOVINT32X2_FROMINT64(U064_h), AE_MOVINT32X2_FROMINT64(U064_l));
    e1 = AE_SEL32_LL(AE_MOVINT32X2_FROMINT64(U164_h), AE_MOVINT32X2_FROMINT64(U164_l));
 
    AE_MULF2P32X4RAS(z0, z1, y0, y1, e0, e1);
    z0 = AE_ADD32(z0, AE_MOVDA32(2));
    z1 = AE_ADD32(z1, AE_MOVDA32(2));
    y0 = AE_SUB32S(y0, AE_SRAI32(z0, 2));
    y1 = AE_SUB32S(y1, AE_SRAI32(z1, 2));
    AE_MOVT32X2(y0, AE_MOVDA32(0x80000000), sx0);
    AE_MOVT32X2(y1, AE_MOVDA32(0x80000000), sx1);
    AE_S32X2X2_I(y0, y1, pScr, 0 * sizeof(ae_int32x4));
    AE_S16X4_I(nsax, (ae_int16x4 *)pScr, 1 * sizeof(ae_int32x4));
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
  
} /* vec_rsqrt32x32() */
