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
  NatureDSP Signal Processing Library. Vector mathematics
   Power function
    code optimized for HiFi5 core
  IntegrIT, 2006-2019
*/

#include "NatureDSP_Signal_math.h"
#include "NatureDSP_types.h"
#include "common.h"
#include "vec_pow_32x32_table.h"

void scl_pow_32x32(int32_t *  m, int16_t *e, int32_t   x, int32_t  y);

#define sz_i32    (int)sizeof(int32_t)
/*-------------------------------------------------------------------------
  Power function
  This routine calculates power function for 32-bit fixed-point numbers or 
  floating point numbers. 
  For the fixed point API, The  base is represented in Q31, the exponent 
  is represented in Q6.25. Results are represented as normalized fixed point
  number with separate mantissa in Q31 and exponent.

  Precision:
  32x32  32-bit inputs, 32-bit outputs
  f      floating point input, floating point output

  Accuracy: 
  2 LSB for fixed point API
  2 ULP under condition that |y|<=100

  Notes:
1. Scalar floating point raise  to a power functions conform to ANSI C requirements on 
   standard math library functions in respect to treatment of errno and floating-
   point exceptions. Vectorized function do not touch errno and may raise or not raise 
   floating point exceptions.
2. For floating point API, If x<0 is finite, y is finite and not an integer value, 
   then the respective result z is set to NaN
3. For fixed point API, function returns zero for all non-positive x. Fixed point 
   functions never touch errno

    Special cases:
          x   |   y    | Result |  Extra Conditions    
      --------+--------+--------+---------------------
      floating point API
      --------+--------+--------+---------------------
        +/-0  | y      | +/-inf | odd y<0
        +/-0  | y      | +inf   | even y<0
        +/-0  | y      | +/-0   | odd y>0
        +/-0  | y      | 0      | even y>0
        +/-1  | +/-inf | 1      | 
        1     | y      | 1      | any y including NaN
        x     | +/-0   | 1      | any x including NaN
        x     | y      | NaN    | finite x<0 and finite 
              |        |        | non-integer y (see 
              |        |        | note 2)
        x     | -inf   | +inf   | |x|<1
        x     | -inf   | 0      | |x|>1
        x     | +inf   | 0      | |x|<1
        x     | +inf   | +inf   | |x|>1
        -inf  | y      | -0     | y an odd integer <0
        -inf  | y      | 0      | y<0 and not an odd 
              |        |        | integer
        -inf  | y      | -inf   | y an odd integer >0
        -inf  | y      | +inf   | y>0 and not an odd 
              |        |        | integer
        +inf  | y      | 0      | y<0
        +inf  | y      | +inf   | y>0
      --------+--------+--------+---------------------
      fixed point API
      --------+--------+--------+---------------------
         x    | y      | 0      | x<=0
      --------+--------+--------+---------------------

  Input:
  x[N]  input data,Q0.31 or floating point
  y[N]  input data,Q6.25 or floating point
  N     length of vectors
  Output (fixed point API):
  m[N]  mantissa of output, Q31 
  e[N]  exponent of output  
  Output (floating point API):
  z[N]  results: floating point

  Restriction:
  z,x,y,m should not overlap
-------------------------------------------------------------------------*/
static void mypow32(int32_t * scr, int32_t *  m, int16_t *e,
               const int32_t *  x, const int32_t*  y, int N)
{
  int n;
  const ae_int32*   restrict pPoly;
  const ae_int32x4* restrict pX;
        ae_int32x4* restrict pY;
        ae_int32x4* restrict pM;
        ae_int16x4* restrict pE;
  const ae_int32x4 * restrict S_rd;
        ae_int32x4 * restrict S_wr;
  ae_valignx2 aX, aY, aM;
  ae_valign aE; 
  ae_int32x2 x0, x1, y0, y1;
  /* Current block index; overall number of blocks; number of values in the current block */
  int blkLen;
  /* Block size, blkLen <= blkSize */
  const int blkSize = (MAX_ALLOCA_SZ / (4*(sizeof(int64_t)))&~3);
  /* Allocate a fixed-size scratch area on the stack. */

  NASSERT_ALIGN16(scr);
  /*
  * Data are processed in blocks of scratch area size. Further, the algorithm
  * implementation is splitted in order to feed the optimizing compiler with a
  * few loops of managable size.
  */
  blkLen = 0;
  for (; N>0; N -= blkLen, x += blkSize, y += blkSize, m += blkSize, e += blkSize)
  {
    blkLen = XT_MIN(N, blkSize);
    pX = (const ae_int32x4*)x;
    pY = (      ae_int32x4*)y;
    S_wr = (ae_int32x4*)scr;
    aX = AE_LA128_PP(pX);
    aY = AE_LA128_PP(pY);

    for (n = 0; n<(blkLen >> 2); n++)
    {
      ae_int16x4 nsa01;
      ae_int32x2 e0, e1;
      xtbool2 bsmall0, bsmall1;

      /* normalization */
      AE_LA32X2X2_IP(x0, x1, aX, pX);
      nsa01 = AE_NSAZ32X4(x0, x1);
      e0 = AE_SEXT32X2D16_32((nsa01));//!!!!!!!!!!!!!!!!!!!
      e1 = AE_SEXT32X2D16_10((nsa01));//!!!!!!!!!!!!!!!!!!!
      x0 = AE_SRAV32RS(x0, AE_NEG32S(e0));
      x1 = AE_SRAV32RS(x1, AE_NEG32S(e1));
      bsmall0 = AE_LT32(x0, 1518500250);
      bsmall1 = AE_LT32(x1, 1518500250);
      AE_MOVT32X2(x0, AE_SLLI32(x0, 1), bsmall0);
      AE_MOVT32X2(x1, AE_SLLI32(x1, 1), bsmall1);
      
      AE_MOVT32X2(e0, AE_ADD32(e0, 1), bsmall0);
      AE_MOVT32X2(e1, AE_ADD32(e1, 1), bsmall1);
      x0 = AE_SUB32(0x80000000, x0);
      x1 = AE_SUB32(0x80000000, x1);
      AE_S32X2X2_IP(x0, x1, S_wr, 2 * sizeof(ae_int32x2));
      AE_S32X2X2_IP(e0, e1, S_wr, 2 * 3 * sizeof(ae_int32x2));
    }
    __Pragma("no_reorder");
    /* compute polynomial */
    pPoly = ((ae_int32*)vec_pow_32x32_polylogQ63) + 1;
    S_rd = (ae_int32x4*)scr;
    S_wr = (ae_int32x4*)scr + 2;
    for (n = 0; n<(blkLen >> 2); n++)
    {
      ae_int32x2 p0, p1, p2, p3, p4, p5, p6, p7, p8;
      ae_int32x2 t0, t1, y0, y1;
      ae_int64 p9, v0, v1, v2, v3;
      AE_L32X2X2_IP(x0, x1, S_rd, 4 * 2 * sizeof(ae_int32x2));

      AE_L32_IP(p0, pPoly, 1 * sizeof(int64_t));
      AE_L32_IP(p1, pPoly, 1 * sizeof(int64_t));
      AE_L32_IP(p2, pPoly, 1 * sizeof(int64_t));
      AE_L32_IP(p3, pPoly, 1 * sizeof(int64_t));
      AE_L32_IP(p4, pPoly, 1 * sizeof(int64_t));
      AE_L32_IP(p5, pPoly, 1 * sizeof(int64_t));
      AE_L32_IP(p6, pPoly, 1 * sizeof(int64_t));
      AE_L32_IP(p7, pPoly, 1 * sizeof(int64_t));
      AE_L32_IP(p8, pPoly, 1 * sizeof(int32_t));
      AE_L64_XP(p9, castxcc(ae_int64, pPoly), -(int)(17 * sizeof(int32_t)));
      /* a portion of polynomial is computed by 32x32 multiplies, but remainder via 64x32 multiplies */
      t0 = t1 = p1; y0 = y1 = p0;      
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
      t0 = t1 = p2;
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
      t0 = t1 = p3;
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
      t0 = t1 = p4;
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
      t0 = t1 = p5;
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
      t0 = t1 = p6;
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
      t0 = t1 = p7;
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
      t0 = t1 = p8;
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;

      v0 = v1 = v2 = v3 = p9;
      //AE_MUL32X2S_HH_LL(v0, v1, y0, x0);
      //AE_MUL32X2S_HH_LL(v2, v3, y1, x1);
      AE_MULAF32S_HH(v0, y0, x0); AE_MULAF32S_LL(v1, y0, x0);
      AE_MULAF32S_HH(v2, y1, x1); AE_MULAF32S_LL(v3, y1, x1);
      AE_S64X2_IP(v0, v1, castxcc(ae_int64x2, S_wr), 2 * sizeof(ae_int64));
      AE_S64X2_IP(v2, v3, castxcc(ae_int64x2, S_wr), 3 * 2 * sizeof(ae_int64));
    }
    __Pragma("no_reorder");
    pPoly = ((const ae_int32*)(vec_pow_32x32_polylogQ63 + 10));
    S_rd = (ae_int32x4*)scr;
    S_wr = (ae_int32x4*)scr + 2;

    for (n = 0; n<(blkLen >> 2); n++)
    {
      ae_int64 w, u;
      ae_int64 v0, v1, v2, v3;
      ae_int32x2 t, e0, e1;
      ae_ep acc_ep;
      AE_L32X2X2_IP(x0, x1, S_rd, 2 * sizeof(ae_int32x2));
      AE_L32X2X2_IP(e0, e1, S_rd, 2 * sizeof(ae_int32x2));
      AE_L64X2_IP(v0, v1, castxcc(ae_int64x2, S_rd), 2 * sizeof(ae_int64));
      AE_L64X2_IP(v2, v3, castxcc(ae_int64x2, S_rd), 2 * sizeof(ae_int64));
     
      AE_L64_IP(w, castxcc(ae_int64, pPoly), sizeof(int64_t));
      t = AE_MOVINT32X2_FROMINT64(v0); AE_MUL32USEP_LH(acc_ep, u, t, x0); u = AE_SRAI64(u, 31); AE_MULAF32S_HH(u, t, x0); v0 = AE_ADD64(u, w);
      t = AE_MOVINT32X2_FROMINT64(v1); AE_MUL32USEP_LL(acc_ep, u, t, x0); u = AE_SRAI64(u, 31); AE_MULAF32S_HL(u, t, x0); v1 = AE_ADD64(u, w);
      t = AE_MOVINT32X2_FROMINT64(v2); AE_MUL32USEP_LH(acc_ep, u, t, x1); u = AE_SRAI64(u, 31); AE_MULAF32S_HH(u, t, x1); v2 = AE_ADD64(u, w);
      t = AE_MOVINT32X2_FROMINT64(v3); AE_MUL32USEP_LL(acc_ep, u, t, x1); u = AE_SRAI64(u, 31); AE_MULAF32S_HL(u, t, x1); v3 = AE_ADD64(u, w);
      AE_L64_IP(w, castxcc(ae_int64, pPoly), sizeof(int64_t));
      t = AE_MOVINT32X2_FROMINT64(v0); AE_MUL32USEP_LH(acc_ep, u, t, x0); u = AE_SRAI64(u, 31); AE_MULAF32S_HH(u, t, x0); v0 = AE_ADD64(u, w);
      t = AE_MOVINT32X2_FROMINT64(v1); AE_MUL32USEP_LL(acc_ep, u, t, x0); u = AE_SRAI64(u, 31); AE_MULAF32S_HL(u, t, x0); v1 = AE_ADD64(u, w);
      t = AE_MOVINT32X2_FROMINT64(v2); AE_MUL32USEP_LH(acc_ep, u, t, x1); u = AE_SRAI64(u, 31); AE_MULAF32S_HH(u, t, x1); v2 = AE_ADD64(u, w);
      t = AE_MOVINT32X2_FROMINT64(v3); AE_MUL32USEP_LL(acc_ep, u, t, x1); u = AE_SRAI64(u, 31); AE_MULAF32S_HL(u, t, x1); v3 = AE_ADD64(u, w);
      AE_L64_IP(w, castxcc(ae_int64, pPoly), sizeof(int64_t));
      t = AE_MOVINT32X2_FROMINT64(v0); AE_MUL32USEP_LH(acc_ep, u, t, x0); u = AE_SRAI64(u, 31); AE_MULAF32S_HH(u, t, x0); v0 = AE_ADD64(u, w);
      t = AE_MOVINT32X2_FROMINT64(v1); AE_MUL32USEP_LL(acc_ep, u, t, x0); u = AE_SRAI64(u, 31); AE_MULAF32S_HL(u, t, x0); v1 = AE_ADD64(u, w);
      t = AE_MOVINT32X2_FROMINT64(v2); AE_MUL32USEP_LH(acc_ep, u, t, x1); u = AE_SRAI64(u, 31); AE_MULAF32S_HH(u, t, x1); v2 = AE_ADD64(u, w);
      t = AE_MOVINT32X2_FROMINT64(v3); AE_MUL32USEP_LL(acc_ep, u, t, x1); u = AE_SRAI64(u, 31); AE_MULAF32S_HL(u, t, x1); v3 = AE_ADD64(u, w);
      AE_L64_IP(w, castxcc(ae_int64, pPoly), sizeof(int64_t));
      t = AE_MOVINT32X2_FROMINT64(v0); AE_MUL32USEP_LH(acc_ep, u, t, x0); u = AE_SRAI64(u, 31); AE_MULAF32S_HH(u, t, x0); v0 = AE_ADD64(u, w);
      t = AE_MOVINT32X2_FROMINT64(v1); AE_MUL32USEP_LL(acc_ep, u, t, x0); u = AE_SRAI64(u, 31); AE_MULAF32S_HL(u, t, x0); v1 = AE_ADD64(u, w);
      t = AE_MOVINT32X2_FROMINT64(v2); AE_MUL32USEP_LH(acc_ep, u, t, x1); u = AE_SRAI64(u, 31); AE_MULAF32S_HH(u, t, x1); v2 = AE_ADD64(u, w);
      t = AE_MOVINT32X2_FROMINT64(v3); AE_MUL32USEP_LL(acc_ep, u, t, x1); u = AE_SRAI64(u, 31); AE_MULAF32S_HL(u, t, x1); v3 = AE_ADD64(u, w);
     
      t = AE_MOVINT32X2_FROMINT64(v0); AE_MUL32USEP_LH(acc_ep, u, t, x0); u = AE_SRAI64(u, 31); AE_MULAF32S_HH(u, t, x0); v0 = u;
      t = AE_MOVINT32X2_FROMINT64(v0); AE_MUL32USEP_LH(acc_ep, u, t, x0); u = AE_SRAI64(u, 31); AE_MULAF32S_HH(u, t, x0); v0 = u;
      t = AE_MOVINT32X2_FROMINT64(v1); AE_MUL32USEP_LL(acc_ep, u, t, x0); u = AE_SRAI64(u, 31); AE_MULAF32S_HL(u, t, x0); v1 = u;
      t = AE_MOVINT32X2_FROMINT64(v1); AE_MUL32USEP_LL(acc_ep, u, t, x0); u = AE_SRAI64(u, 31); AE_MULAF32S_HL(u, t, x0); v1 = u;
      t = AE_MOVINT32X2_FROMINT64(v2); AE_MUL32USEP_LH(acc_ep, u, t, x1); u = AE_SRAI64(u, 31); AE_MULAF32S_HH(u, t, x1); v2 = u;
      t = AE_MOVINT32X2_FROMINT64(v2); AE_MUL32USEP_LH(acc_ep, u, t, x1); u = AE_SRAI64(u, 31); AE_MULAF32S_HH(u, t, x1); v2 = u;
      t = AE_MOVINT32X2_FROMINT64(v3); AE_MUL32USEP_LL(acc_ep, u, t, x1); u = AE_SRAI64(u, 31); AE_MULAF32S_HL(u, t, x1); v3 = u;
      t = AE_MOVINT32X2_FROMINT64(v3); AE_MUL32USEP_LL(acc_ep, u, t, x1); u = AE_SRAI64(u, 31); AE_MULAF32S_HL(u, t, x1); v3 = u;
     
      // final scaling by log2(e)
      u = AE_MOVINT64_FROMINT32X2(AE_SEL32_HH(x0, AE_ZERO32())); v0 = AE_SUB64(v0, u);
      u = AE_MOVINT64_FROMINT32X2(AE_SEL32_LL(x0, AE_ZERO32())); v1 = AE_SUB64(v1, u);
      u = AE_MOVINT64_FROMINT32X2(AE_SEL32_HH(x1, AE_ZERO32())); v2 = AE_SUB64(v2, u);
      u = AE_MOVINT64_FROMINT32X2(AE_SEL32_LL(x1, AE_ZERO32())); v3 = AE_SUB64(v3, u);
      AE_L64_XP(w, castxcc(ae_int64, pPoly), -(int)(4 * sizeof(int64_t)));
      t = AE_MOVINT32X2_FROMINT64(v0); AE_MULZAAD32USEP_HL_LH(acc_ep, v0, t, AE_MOVINT32X2_FROMINT64(w)); v0 = AE_SRAI72(acc_ep, v0, 31);  AE_MULAF32S_HH(v0, t, AE_MOVINT32X2_FROMINT64(w));
      t = AE_MOVINT32X2_FROMINT64(v1); AE_MULZAAD32USEP_HL_LH(acc_ep, v1, t, AE_MOVINT32X2_FROMINT64(w)); v1 = AE_SRAI72(acc_ep, v1, 31);  AE_MULAF32S_HH(v1, t, AE_MOVINT32X2_FROMINT64(w));
      t = AE_MOVINT32X2_FROMINT64(v2); AE_MULZAAD32USEP_HL_LH(acc_ep, v2, t, AE_MOVINT32X2_FROMINT64(w)); v2 = AE_SRAI72(acc_ep, v2, 31);  AE_MULAF32S_HH(v2, t, AE_MOVINT32X2_FROMINT64(w));
      t = AE_MOVINT32X2_FROMINT64(v3); AE_MULZAAD32USEP_HL_LH(acc_ep, v3, t, AE_MOVINT32X2_FROMINT64(w)); v3 = AE_SRAI72(acc_ep, v3, 31);  AE_MULAF32S_HH(v3, t, AE_MOVINT32X2_FROMINT64(w));
     
     
      e0 = AE_NEG32(e0); e1 = AE_NEG32(e1);
      v0 = AE_SRAI64(v0, 5); t = AE_SLAI32(AE_SEL32_HH(e0, AE_ZERO32()), 25); v0 = AE_ADD64(AE_MOVINT64_FROMINT32X2(t), v0);
      v1 = AE_SRAI64(v1, 5); t = AE_SLAI32(AE_SEL32_LL(e0, AE_ZERO32()), 25); v1 = AE_ADD64(AE_MOVINT64_FROMINT32X2(t), v1);
      v2 = AE_SRAI64(v2, 5); t = AE_SLAI32(AE_SEL32_HH(e1, AE_ZERO32()), 25); v2 = AE_ADD64(AE_MOVINT64_FROMINT32X2(t), v2);
      v3 = AE_SRAI64(v3, 5); t = AE_SLAI32(AE_SEL32_LL(e1, AE_ZERO32()), 25); v3 = AE_ADD64(AE_MOVINT64_FROMINT32X2(t), v3);
      /* here we have f[0] in range -0.5...0.5 and exponent in range -31...0 */
      AE_S64X2_IP(v0, v1, castxcc(ae_int64x2, S_wr), 2 * sizeof(ae_int64));
      AE_S64X2_IP(v2, v3, castxcc(ae_int64x2, S_wr), 3 * 2 * sizeof(ae_int64));
    }

    __Pragma("no_reorder");
    /* scale, extract exponential part */
    S_rd = (ae_int32x4*)scr + 2;
    S_wr = (ae_int32x4*)scr;
    pE = (ae_int16x4*)e;
    aE = AE_ZALIGN64();
    pX = (const ae_int32x4*)x;
    aX = AE_LA128_PP(pX);
    for (n = 0; n<(blkLen >> 2); n++)
    {
      ae_int16x4 h0;
      ae_int32x2 e0, e1;
      ae_int32x2 s0, s1, s2, s3;
      ae_int64 v0, v1, v2, v3;
      xtbool2 zerox0, zerox1;
      ae_ep acc_ep;
      AE_LA32X2X2_IP(x0, x1, aX, pX);
      AE_LA32X2X2_IP(y0, y1, aY, pY);
      AE_L64X2_IP(v0, v1, castxcc(ae_int64x2, S_rd), 2 * sizeof(ae_int64));
      AE_L64X2_IP(v2, v3, castxcc(ae_int64x2, S_rd), 3 * 2 * sizeof(ae_int64));
      
      zerox0 = AE_LE32(x0, AE_ZERO32());
      zerox1 = AE_LE32(x1, AE_ZERO32());
      /* multiply log2 to form Q51 result */
      s0 = AE_MOVINT32X2_FROMINT64(v0);
      s1 = AE_MOVINT32X2_FROMINT64(v1);
      s2 = AE_MOVINT32X2_FROMINT64(v2);
      s3 = AE_MOVINT32X2_FROMINT64(v3);
      
      AE_MUL32USEP_LH(acc_ep, v0, s0, y0);
      AE_MUL32USEP_LL(acc_ep, v1, s1, y0);
      AE_MUL32USEP_LH(acc_ep, v2, s2, y1);
      AE_MUL32USEP_LL(acc_ep, v3, s3, y1);
      v0 = AE_SRAI64(v0, 31);
      v1 = AE_SRAI64(v1, 31);
      v2 = AE_SRAI64(v2, 31);
      v3 = AE_SRAI64(v3, 31);
      
      AE_MULAF32S_HH(v0, s0, y0);
      AE_MULAF32S_HL(v1, s1, y0);
      AE_MULAF32S_HH(v2, s2, y1);
      AE_MULAF32S_HL(v3, s3, y1);
      /* exponent is found via ceil, mantissa is computed from fractional part */
      v0 = AE_NEG64(v0);
      v1 = AE_NEG64(v1);
      v2 = AE_NEG64(v2);
      v3 = AE_NEG64(v3);
      e0 = AE_TRUNCA32X2F64S(v0, v1, -19);
      e1 = AE_TRUNCA32X2F64S(v2, v3, -19);
      e0 = AE_NEG32(e0);
      e1 = AE_NEG32(e1);
      v0 = AE_AND64(v0, 0x7ffffffffffffULL);
      v1 = AE_AND64(v1, 0x7ffffffffffffULL);
      v2 = AE_AND64(v2, 0x7ffffffffffffULL);
      v3 = AE_AND64(v3, 0x7ffffffffffffULL);
      v0 = AE_NEG64(v0);
      v1 = AE_NEG64(v1);
      v2 = AE_NEG64(v2);
      v3 = AE_NEG64(v3);
      y0 = AE_TRUNCA32X2F64S(v0, v1, 12);
      y1 = AE_TRUNCA32X2F64S(v2, v3, 12);
      AE_S32X2X2_IP(y0, y1, S_wr, 4 * 2 * sizeof(ae_int32x2));
      
      AE_MOVT32X2(e0, AE_ZERO32(), zerox0);
      AE_MOVT32X2(e1, AE_ZERO32(), zerox1);
      
      h0 = AE_TRUNCI16X4F32S(e0, e1, 16);
      AE_SA16X4_IP(h0, aE, pE);
    }
    AE_SA64POS_FP(aE, pE);
    __Pragma("no_reorder");
    /* compute mantissa as 2^x from fractional part */
    S_rd = (ae_int32x4*)scr;
    pM = (ae_int32x4*)m;
    pX = (const ae_int32x4*)x;
    aX = AE_LA128_PP(pX);
    aM = AE_ZALIGN128();
    pPoly = (const ae_int32*)vec_pow_32x32_polypow2;
    for (n = 0; n<(blkLen >> 2); n++)
    {
      ae_int32x2 t0, t1;
      ae_int32x2 p0, p1, p2, p3, p4, p5, p6, p7, p8;
      xtbool2 zerox0, zerox1;
      AE_LA32X2X2_IP(x0, x1, aX, pX);
      zerox0 = AE_LE32(x0, AE_ZERO32());
      zerox1 = AE_LE32(x1, AE_ZERO32());
      AE_L32X2X2_IP(x0, x1, S_rd, 4 * 2 * sizeof(ae_int32x2));
      x0 = AE_ADD32(x0, (1 << 30));
      x1 = AE_ADD32(x1, (1 << 30));
      p0 = AE_L32_I(pPoly, 0 * sizeof(int32_t));
      p1 = AE_L32_I(pPoly, 1 * sizeof(int32_t));
      p2 = AE_L32_I(pPoly, 2 * sizeof(int32_t));
      p3 = AE_L32_I(pPoly, 3 * sizeof(int32_t));
      p4 = AE_L32_I(pPoly, 4 * sizeof(int32_t));
      p5 = AE_L32_I(pPoly, 5 * sizeof(int32_t));
      p6 = AE_L32_I(pPoly, 6 * sizeof(int32_t));
      p7 = AE_L32_I(pPoly, 7 * sizeof(int32_t));
      p8 = AE_L32_X(pPoly, 8 * sizeof(int32_t));

      t0 = t1 = p1; y0 = y1 = p0;
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
      t0 = t1 = p2;
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
      t0 = t1 = p3;
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
      t0 = t1 = p4;
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
      t0 = t1 = p5;
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
      t0 = t1 = p6;
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
      t0 = t1 = p7;
      AE_MULF2P32X4RAS(y0, y1, x0, x1, y0, y1);
      y0 = AE_SRAI32(y0, 1); y1 = AE_SRAI32(y1, 1);
      y0 = AE_ADD32S(y0, t0); y1 = AE_ADD32S(y1, t1);
      t0 = t1 = p8;
      AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
      AE_MOVT32X2(y0, AE_ZERO32(), zerox0);
      AE_MOVT32X2(y1, AE_ZERO32(), zerox1);
      AE_SA32X2X2_IP(y0, y1, aM, pM);
    }
    AE_SA128POS_FP(aM, pM);
  }
}
void vec_pow_32x32(int32_t *  m, int16_t *e,
                   const int32_t *  x, const int32_t*  y, int N)
{
  const int blkSize = (MAX_ALLOCA_SZ /( (sizeof(int64_t)))&~3);
  /* Allocate a fixed-size scratch area on the stack. */
  int32_t ALIGN(32) scr[blkSize];
  int32_t ALIGN(32) tmpInX[4], tmpInY[4], tmpOutM[4];
  int16_t ALIGN(32) tmpOutE[4];

  int M;
  if ( N<=0 ) return;
  M=N&~3;
  if ( M )
  {
    mypow32(scr,m,e,x,y,M); 
    y += M;
    x += M;
    m += M;
    e += M;
    N&=3;
  }
  if (N)
  {
    // processing the tail
    int off1,off2;
    ae_int32x2 x0, x1, x2;
    ae_int32x2 y0, y1, y2;
    ae_int16x4 h0, h1, h2;
    off1=XT_MIN(N-1,1)<<2;
    off2=XT_MIN(N-1,2)<<2;
    x0 = AE_L32_I((const ae_int32*)x, 0);
    x1 = AE_L32_X((const ae_int32*)x, off1);
    x2 = AE_L32_X((const ae_int32*)x, off2);
    AE_S32_L_I(x0, (ae_int32*)tmpInX, 0 * sizeof(int32_t));
    AE_S32_L_I(x1, (ae_int32*)tmpInX, 1 * sizeof(int32_t));
    AE_S32_L_I(x2, (ae_int32*)tmpInX, 2 * sizeof(int32_t));
    y0 = AE_L32_I((const ae_int32*)y, 0);
    y1 = AE_L32_X((const ae_int32*)y, off1);
    y2 = AE_L32_X((const ae_int32*)y, off2);
    AE_S32_L_I(y0, (ae_int32*)tmpInY, 0 * sizeof(int32_t));
    AE_S32_L_I(y1, (ae_int32*)tmpInY, 1 * sizeof(int32_t));
    AE_S32_L_I(y2, (ae_int32*)tmpInY, 2 * sizeof(int32_t));


    mypow32(scr, tmpOutM, tmpOutE, tmpInX, tmpInY, 4);

    x0 = AE_L32_I((const ae_int32*)tmpOutM, 0 * sizeof(int32_t));
    x1 = AE_L32_I((const ae_int32*)tmpOutM, 1 * sizeof(int32_t));
    x2 = AE_L32_I((const ae_int32*)tmpOutM,2*sizeof(int32_t));
    AE_S32_L_X(x2,(ae_int32*)m,off2);
    AE_S32_L_X(x1,(ae_int32*)m,off1);
    AE_S32_L_I(x0,(ae_int32*)m,0);

    h0 = AE_L16_I((const ae_int16*)tmpOutE, 0 * sizeof(int16_t));
    h1 = AE_L16_I((const ae_int16*)tmpOutE, 1 * sizeof(int16_t));
    h2 = AE_L16_I((const ae_int16*)tmpOutE, 2 * sizeof(int16_t));
    AE_S16_0_X(h2,(ae_int16*)e,(off2>>1));
    AE_S16_0_X(h1,(ae_int16*)e,(off1>>1));
    AE_S16_0_I(h0,(ae_int16*)e,0);
  }
} /* vec_pow_32x32() */
