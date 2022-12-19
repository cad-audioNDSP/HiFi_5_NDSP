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


/*
  NatureDSP Signal Processing Library. Vector matematics
    Division
    Optimized code for HiFi5
  IntegrIT, 2006-2019
*/

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
void vec_divide64x32i
                (int32_t * restrict frac, 
                 const int64_t * restrict x, 
                 const int32_t * restrict y, int N)
{
/*algorithm:
  f(x,y) = div(x,y) = x/y = x*recip(y)

1)  f(y) = recip(y) = 1/y
Rough calculation:
  f(y) = f(y0) + f'(y0)*(y1-y0)
  1/y0 = f(y0)
  1/(y0^2) = f'(y0)
  y1 = y*2^nsa, y1 in 0.5..1
  y1-y0 = dy
Refiniment:
  err = f(y)*y - 1 ,error of calculation
  recip(y) = f(y) - err*f(y), refined value
2)f(x,y) = x*recip(y)
  */
  int n;
  const ae_int64    *          px   = (const ae_int64   *)x;
  const ae_int32x4  *          pY   = (const ae_int32x4 *)y;
        ae_int32x4  *          pZ   = (      ae_int32x4 *)frac;

  unsigned char sg0, sg1;
  xtbool2 b0, ovl;  
  xtbool  bl0, bl1;
  ae_int64 x2, x3;
  ae_int32x2 y1, z2, z3;
  ae_int64 x0, x1, s1, s2;
  ae_int32x2 y0, z0, s0, z1;
  ae_int32x2 maxint = AE_MOVDA32X2( 0x7FFFFFFF,0x7FFFFFFF );
  ae_valign y_align, z_align;
  ae_valignx2 aY, aZ;
  if(N<=0) return;

  aZ = AE_ZALIGN128();
  aY = AE_LA128_PP(pY);
  
  s0 = AE_MOVDA32X2( 0x00000000,0xFFFFFFFF );
  s2 = AE_MOVINT64_FROMINT32X2(s0);//;AE_MOV64( 0x00000000FFFFFFFF);
  
  s0 = AE_ZERO32(); 
  s1 = AE_ZERO64();  
 
  for (n=0;n<(N>>2);n++)
  {
    AE_L64_IP(x0, px, sizeof(ae_int64));
    AE_L64_IP(x1, px, sizeof(ae_int64));
    AE_L64_IP(x2, px, sizeof(ae_int64));
    AE_L64_IP(x3, px, sizeof(ae_int64));
    AE_LA32X2X2_IP(y0, y1, aY, pY);

    b0 = AE_LT32(y0,s0);
    bl0 = AE_LT64(x0,s1);
    bl1 = AE_LT64(x1,s1);
    sg1 = AE_MOVAB2(b0);
    sg0 = ((AE_MOVAB(bl0)*2)|(AE_MOVAB(bl1)))&0x3;
    sg0 = sg0^sg1;
    x0 = AE_ABS64S(x0);
    x1 = AE_ABS64S(x1);
    y0 = AE_ABS32S(y0);

    z0 = AE_TRUNCA32X2F64S(x0,x1,0);
    ovl = AE_LE32(y0,z0);

    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);

    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);
    AE_DIV64D32_L(x1,y0);

    /////////////////////////////

    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);

    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);
    AE_DIV64D32_H(x0,y0);

    x0 = AE_AND64(x0,s2);
    x1 = AE_AND64(x1,s2);
    x0 = AE_SLAI64S(x0,32);
    x1 = AE_SLAI64S(x1,32);
    z0 = AE_TRUNCA32X2F64S(x0,x1,0);

    AE_MOVT32X2(z0,maxint,ovl);

    z1 = AE_NEG32S(z0);
    AE_MOVT32X2(z0,z1,sg0);

    /**************************************/
    b0 = AE_LT32(y1, s0);
    bl0 = AE_LT64(x2, s1);
    bl1 = AE_LT64(x3, s1);
    sg1 = AE_MOVAB2(b0);
    sg0 = ((AE_MOVAB(bl0) * 2) | (AE_MOVAB(bl1))) & 0x3;
    sg0 = sg0^sg1;
    x2 = AE_ABS64S(x2);
    x3 = AE_ABS64S(x3);
    y1 = AE_ABS32S(y1);

    z2 = AE_TRUNCA32X2F64S(x2, x3, 0);
    ovl = AE_LE32(y1, z2);

    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);

    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);
    AE_DIV64D32_L(x3, y1);

    /////////////////////////////

    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);

    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);
    AE_DIV64D32_H(x2, y1);

    x2 = AE_AND64(x2, s2);
    x3 = AE_AND64(x3, s2);
    x2 = AE_SLAI64S(x2, 32);
    x3 = AE_SLAI64S(x3, 32);
    z2 = AE_TRUNCA32X2F64S(x2, x3, 0);

    AE_MOVT32X2(z2, maxint, ovl);

    z3 = AE_NEG32S(z2);
    AE_MOVT32X2(z2, z3, sg0);

    AE_SA32X2X2_IP(z0, z2, aZ, pZ);
  } 
  AE_SA128POS_FP(aZ, pZ);
  if ((N & 3)>1)
  {
    z_align = AE_ZALIGN64();
    y_align = AE_LA64_PP(pY);
    AE_L64_IP(x0, px, sizeof(ae_int64));
    AE_L64_IP(x1, px, sizeof(ae_int64));
    AE_LA32X2_IP(y0, y_align, castxcc(ae_int32x2, pY));

    b0 = AE_LT32(y0, s0);
    bl0 = AE_LT64(x0, s1);
    bl1 = AE_LT64(x1, s1);
    sg1 = AE_MOVAB2(b0);
    sg0 = ((AE_MOVAB(bl0) * 2) | (AE_MOVAB(bl1))) & 0x3;
    sg0 = sg0^sg1;
    x0 = AE_ABS64S(x0);
    x1 = AE_ABS64S(x1);
    y0 = AE_ABS32S(y0);

    z0 = AE_TRUNCA32X2F64S(x0, x1, 0);
    ovl = AE_LE32(y0, z0);

    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);

    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);
    AE_DIV64D32_L(x1, y0);

    /////////////////////////////

    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);

    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);

    x0 = AE_AND64(x0, s2);
    x1 = AE_AND64(x1, s2);
    x0 = AE_SLAI64S(x0, 32);
    x1 = AE_SLAI64S(x1, 32);
    z0 = AE_TRUNCA32X2F64S(x0, x1, 0);

    AE_MOVT32X2(z0, maxint, ovl);

    z1 = AE_NEG32S(z0);
    AE_MOVT32X2(z0, z1, sg0);
    AE_SA32X2_IP(z0, z_align, castxcc(ae_int32x2, pZ));
    AE_SA64POS_FP(z_align, castxcc(ae_int32x2, pZ));
  }

  if (N & 1)
  {
    AE_L64_IP(x0, px, sizeof(ae_int64));
    y0 = AE_L32_I((const ae_int32 *)pY, 0);

    b0 = AE_LT32(y0, s0);
    bl0 = AE_LT64(x0, s1);

    sg1 = AE_MOVAB2(b0);
    sg0 = ((AE_MOVAB(bl0) * 2) | (AE_MOVAB(bl0)));
    sg0 = sg0^sg1;
    x0 = AE_ABS64S(x0);
    y0 = AE_ABS32S(y0);

    z0 = AE_TRUNCA32X2F64S(x0, x0, 0);
    ovl = AE_LE32(y0, z0);

    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);

    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);

    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);

    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);

    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);

    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);

    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);

    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);
    AE_DIV64D32_H(x0, y0);

    x0 = AE_AND64(x0, s2);
    x0 = AE_SLAI64S(x0, 32);
    z0 = AE_TRUNCA32X2F64S(x0, x0, 0);

    AE_MOVT32X2(z0, maxint, ovl);

    z1 = AE_NEG32S(z0);
    AE_MOVT32X2(z0, z1, sg0);
    AE_S32_L_I(z0, (ae_int32 *)pZ, 0);
  }
} /* vec_divide64x32i() */
