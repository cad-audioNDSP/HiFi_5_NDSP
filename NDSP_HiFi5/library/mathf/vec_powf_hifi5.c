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
    Vector operations
    code optimized for HiFi5 core
  IntegrIT, 2006-2019
*/

#include "NatureDSP_Signal_math.h"
#include "NatureDSP_types.h"
#include "common.h"
/* Common helper macros. */
#include "common_fpu.h"
/* Constant tables. */
#include "pow2f_tbl.h"
/* +/-Infinity, single/double precision */
#include "inff_tbl.h"
#include "nanf_tbl.h"
#include "sqrt2f_tbl.h"
#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )


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

#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(void, vec_powf, (float32_t * restrict z, const float32_t * restrict y, const float32_t * restrict x, int N))
#elif HAVE_VFPU
#define sz_f32    (int)sizeof(float32_t)
static void mypowf(float32_t * scr,
                  float32_t * restrict z, 
            const float32_t * restrict x, 
            const float32_t * restrict y, 
            int N )
{
  /* Table of different constants used in computations */
  static const int32_t c_tbl[] =
  {
    -126,
    -150,
    (int32_t)0x007FFFFF,/* max denormalized floating-point number / mantissa mask */
    (int32_t)0x4B800000,/* 2^24 */
    (int32_t)0x3F3504F3,/* sqrt(0.5) */
    (int32_t)0x3F000000,/*  0.5 */
    (int32_t)0xBF000000,/* -0.5 */
    -252,
    254
  };
  int n;
  const xtfloatx4     *          pX;
  const xtfloatx4     *          pY;

  const xtfloatx4     * restrict S_rd;
        xtfloatx4     * restrict S_wr;
        xtfloatx4     * restrict pZ;
  const ae_int32      * restrict TBL;
  const  xtfloat      * restrict TBL_LOG2;
  const  xtfloat      * restrict TBL_POW2;
  xtfloatx2 x0, x1, y0, z0, z1, t0, t1;
  xtfloatx2 c2f, c3f, c4f;
  xtfloatx2 _0, _1, half;
  ae_int32x2 c0i, c1i, c5i, c7i, c8i;
  ae_int32x2 e0, e1, xi0, yi0, ex0, ex1;
  xtbool2 bsx, bsy;
  ae_valignx2 aX, aY;

  /* overall number of blocks; number of values in the current block */
  int blkLen;
  /* Block size, blkLen <= blkSize */
  const int blkSize = (MAX_ALLOCA_SZ / (3 * sz_f32))&~3;



  if (N <= 0) return;

  NASSERT(N % 4 == 0);
  NASSERT_ALIGN16(scr);

  /*
  * Data are processed in blocks of scratch area size. Further, the algorithm
  * implementation is splitted in order to feed the optimizing compiler with a
  * few loops of managable size.
  */


  blkLen = 0;
  TBL = (const ae_int32 *)c_tbl;
  for (; N>0; N -= blkLen, x += blkSize, y += blkSize, z += blkSize)
  {
    blkLen = XT_MIN(N, blkSize);
    _0 = 0.0f;
    _1 = (1.0f);
    half = (0.5f);
    {
      pX   = (const xtfloatx4*)x;
      S_wr = (      xtfloatx4*)scr;
      aX = AE_LA128_PP(pX);
      /*13*/
      for (n = 0; n<(blkLen >> 2); n++)
      {
        ae_int32x2 xi0, xi1, ex0, ex1;
        xtfloatx2 ef0, ef1;
        xtbool2 bdenorm0, bdenorm1, bsmall0, bsmall1;
        AE_LASX2X2_IP(x0, x1, aX, pX);

        ABS_SX2X2(x0, x1, x0, x1);
        c0i = AE_L32_I(TBL, 0 * 4); /*-126*/
        c1i = AE_L32_I(TBL, 1 * 4); /*-150*/
        c2f = XT_LSI((xtfloat*)TBL, 2 * 4);
        c3f = XT_LSI((xtfloat*)TBL, 3 * 4);
        /* process denormalized values */
        
        bdenorm0 = XT_OLE_SX2(x0, c2f);
        bdenorm1 = XT_OLE_SX2(x1, c2f);
        MULQ_S(t0, t1, x0, x1, c3f);
        XT_MOVT_SX2(x0, t0, bdenorm0);
        XT_MOVT_SX2(x1, t1, bdenorm1);
        e0 = e1 = c0i;
        AE_MOVT32X2(e0, c1i, bdenorm0);
        AE_MOVT32X2(e1, c1i, bdenorm1);
        /* extract exponent */
        xi0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(x0);
        xi1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(x1);

        ex0 = AE_SRLI32(xi0, 23);
        ex1 = AE_SRLI32(xi1, 23);
        e0 = AE_ADD32(e0, ex0);
        e1 = AE_ADD32(e1, ex1);
        /* extract mantissa */
        ex0 = ex1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(c2f);/* load mantissa mask */ 
        c5i = AE_L32_I(TBL, 5 * 4);/*  0.5 */
        xi0 = AE_AND32(xi0, ex0);
        xi1= AE_AND32(xi1, ex1);
        xi0 = AE_OR32(xi0, c5i);
        xi1 = AE_OR32(xi1, c5i);
        x0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(xi0);
        x1 = XT_AE_MOVXTFLOATX2_FROMINT32X2(xi1);
        /* adjust the mantissa to range [ sqrt(0.5) ; sqrt(2.0) ) */
        c4f = XT_LSI((xtfloat*)TBL, 4 * 4);
        bsmall0 = XT_OLT_SX2(x0, c4f);
        bsmall1 = XT_OLT_SX2(x1, c4f);
        ADD_SX2X2(t0, t1, x0, x1, x0, x1);
        ex0 = AE_SUB32(e0, 1);
        ex1 = AE_SUB32(e1, 1);
        XT_MOVT_SX2(x0, t0, bsmall0);
        XT_MOVT_SX2(x1, t1, bsmall1);
        AE_MOVT32X2(e0, ex0, bsmall0);
        AE_MOVT32X2(e1, ex1, bsmall1);
        SUB_SX2X2(x0, x1, _1, _1, x0, x1); 
        ef0 = XT_FLOAT_SX2(e0, 0);
        ef1 = XT_FLOAT_SX2(e1, 0);
        AE_SSX2X2_IP(x0, x1, S_wr, 4 * sz_f32);
        AE_SSX2X2_IP(ef0, ef1, S_wr, 2 * 4 * sz_f32);
      }
    }
    __Pragma("no_reorder");
    /*****************************************/
    {
      xtfloatx2 p0, p1, p2, p3, p4, p5, p6, p7, p8, p9;
      xtfloatx2 p10, p11, p12, p13;
      xtfloatx2 w0, w1, w2, w3;
      S_wr = (      xtfloatx4*)scr+2;
      S_rd = (const xtfloatx4*)scr;
      TBL_LOG2 = (const xtfloat *)log2f_coef;
      /*38/4*/
      for (n = 0; n<(blkLen >> 2); n++)
      {
        xtfloatx2 y0, y1;
        AE_LSX2X2_IP(x0, x1, S_rd, 3 * 4 * sz_f32);

        /* evaluate polynomial approximation */
        /* Load table of coefficients        */

        p0 = XT_LSI(TBL_LOG2, 0 * 4);
        p1 = XT_LSI(TBL_LOG2, 1 * 4);
        p2 = XT_LSI(TBL_LOG2, 2 * 4);
        p3 = XT_LSI(TBL_LOG2, 3 * 4);
        p4 = XT_LSI(TBL_LOG2, 4 * 4);
        p5 = XT_LSI(TBL_LOG2, 5 * 4);
        p6 = XT_LSI(TBL_LOG2, 6 * 4);
        p7 = XT_LSI(TBL_LOG2, 7 * 4);
        p8 = XT_LSX(TBL_LOG2, 8 * 4);
        p9 = XT_LSX(TBL_LOG2, 9 * 4);

        t0 = t1 = p1;
        MADDQ_S(t0, t1, x0, x1, p0); 
        y0 = y1 = p2;
        MADD_SX2X2(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
        y0 = y1 = p3;
        MADD_SX2X2(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
        y0 = y1 = p4;
        MADD_SX2X2(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
        y0 = y1 = p5;
        MADD_SX2X2(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
        y0 = y1 = p6;
        MADD_SX2X2(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
        y0 = y1 = p7;
        MADD_SX2X2(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
        y0 = y1 = p8;
        MADD_SX2X2(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
        y0 = y1 = p9;
        MADD_SX2X2(y0, y1, x0, x1, t0, t1); 
        
        AE_SSX2X2_IP(y0, y1, S_wr, 3 * 4 * sz_f32);   
      }
      __Pragma("no_reorder");
      S_wr = (      xtfloatx4*)scr;
      S_rd = (const xtfloatx4*)scr;
      /*29*/
      for (n = 0; n<(blkLen >> 2); n++)
      {
        xtfloatx2 y0, y1, t2, t3;
        xtfloatx2 ef0, ef1;
        p10 = XT_LSX(TBL_LOG2, 10 * 4);
        p11 = XT_LSX(TBL_LOG2, 11 * 4);
        p12 = XT_LSX(TBL_LOG2, 12 * 4);
        p13 = XT_LSX(TBL_LOG2, 13 * 4);

        AE_LSX2X2_IP(x0, x1, S_rd, 4 * sz_f32);
        AE_LSX2X2_IP(ef0, ef1, S_rd, 4 * sz_f32);
        AE_LSX2X2_IP(y0, y1, S_rd, 4 * sz_f32);

        /* next coefficients are computed in extended precision */
        MUL_SX2X2(t0, t1, x0, x1, y0, y1); t2 = t0; t3 = t1;
        MSUB_SX2X2(t2, t3, x0, x1, y0, y1);
        ADD_SX2X2(w0, w1, t0, t1, p10, p10);
        SUB_SX2X2(w2, w3, w0, w1, p10, p10);
        SUB_SX2X2(w2, w3, t0, t1, w2, w3);
        SUB_SX2X2(w2, w3, w2, w3, t2, t3);
        t0 = w0; t1 = w1;
        t2 = w2; t3 = w3;
        MUL_SX2X2(w0, w1, x0, x1, t0, t1); w2 = w0; w3 = w1;
        MSUB_SX2X2(w2, w3, x0, x1, t0, t1); t0 = w0; t1 = w1;
        MSUB_SX2X2(w2, w3, x0, x1, t2, t3); t2 = w2; t3 = w3;
        ADD_SX2X2(w0, w1, t0, t1, p11, p11);
        SUB_SX2X2(w2, w3, w0, w1, p11, p11);
        SUB_SX2X2(w2, w3, t0, t1, w2, w3);
        SUB_SX2X2(w2, w3, w2, w3, t2, t3);
        t0 = w0; t1 = w1;
        t2 = w2; t3 = w3;

        NEG_SX2X2(x0, x1, x0, x1);
        MUL_SX2X2(w0, w1, x0, x1, t0, t1); w2 = w0; w3 = w1;
        MSUB_SX2X2(w2, w3, x0, x1, t0, t1); t0 = w0; t1 = w1;
        MSUB_SX2X2(w2, w3, x0, x1, t2, t3); t2 = w2; t3 = w3;

       /* multiply by log2(e) */
        MULQ_S(w0, w1, t0, t1, p12); w2 = w0; w3 = w1;
        MSUBQ_S(w2, w3, t0, t1, p12);
        MADDQ_S(w2, w3, t2, t3, p12);
        MSUBQ_S(w2, w3, t0, t1, p13);
        t0 = w0; t1 = w1;
        t2 = w2; t3 = w3;
        /* add exponent */
        ADD_SX2X2(w0, w1, t0, t1, ef0, ef1);
        SUB_SX2X2(w2, w3, w0, w1, ef0, ef1);
        SUB_SX2X2(w2, w3, t0, t1, w2, w3);
        SUB_SX2X2(t2, t3, w2, w3, t2, t3);
        t0 = w0; t1 = w1;
        AE_SSX2X2_IP(t0, t1, S_wr, 4 * sz_f32);
        AE_SSX2X2_IP(t2, t3, S_wr, 2 * 4 * sz_f32);
      }    
    }
    __Pragma("no_reorder");
    /*******************************************/
    {
      xtfloatx2 xy0, xy1, dxy0, dxy1;
      xtfloatx2 c0_0, c1_0, c0_1, c1_1;
      xtfloatx2 t2, t3, y0, y1, y2, y3;
      xtfloatx2 p0, p1, p2, p3, p4, p5, p6;
      S_wr = (      xtfloatx4*)scr + 2;
      S_rd = (const xtfloatx4*)scr;
      TBL_POW2 = (const xtfloat *)pow2f_coef;
      pY   = (const xtfloatx4*)y;
      aY = AE_LA128_PP(pY);
      /*34/2*/
      for (n = 0; n<(blkLen >> 2); n++)
      {
        AE_LSX2X2_IP(t0, t1, S_rd, 4 * sz_f32);
        AE_LSX2X2_IP(t2, t3, S_rd, 2 * 4 * sz_f32);
        AE_LASX2X2_IP(y0, y1, aY, pY);

        /* compute y*log2(x) and separate result into integer and fractional parts */
        MUL_SX2X2(y2, y3, y0, y1, t0, t1);
        xy0 = XT_FIROUND_SX2(y2);
        xy1 = XT_FIROUND_SX2(y3);
        NEG_SX2X2(dxy0, dxy1, xy0, xy1);
        MADD_SX2X2(dxy0, dxy1, y0, y1, t0, t1);
        MADD_SX2X2(dxy0, dxy1, y0, y1, t2, t3);
        dxy0 = XT_MIN_SX2(dxy0, (xtfloatx2)1.0f);
        dxy0 = XT_MAX_SX2(dxy0, (xtfloatx2)-1.0f);
        dxy1 = XT_MIN_SX2(dxy1, (xtfloatx2)1.0f);
        dxy1 = XT_MAX_SX2(dxy1, (xtfloatx2)-1.0f);

        /* compute 2^fract */
        p0 = XT_LSI(TBL_POW2, 0 * 4);
        p1 = XT_LSI(TBL_POW2, 1 * 4);
        p2 = XT_LSI(TBL_POW2, 2 * 4);
        p3 = XT_LSI(TBL_POW2, 3 * 4);
        p4 = XT_LSI(TBL_POW2, 4 * 4);
        p5 = XT_LSI(TBL_POW2, 5 * 4);
        p6 = XT_LSI(TBL_POW2, 6 * 4);
        /* NOTE: do not change the order of computations and way of polynomial decomposition ! */
        t0 = t1 = p1;
        MADDQ_S(t0, t1, dxy0, dxy1, p0);
        y0 = y1 = p2;
        MADD_SX2X2(y0, y1, dxy0, dxy1, t0, t1); t0 = y0; t1 = y1;
        y0 = y1 = p3;
        MADD_SX2X2(y0, y1, dxy0, dxy1, t0, t1); t0 = y0; t1 = y1;
        y0 = y1 = p4;
        MADD_SX2X2(y0, y1, dxy0, dxy1, t0, t1); t0 = y0; t1 = y1;
        y0 = y1 = p5;
        MADD_SX2X2(y0, y1, dxy0, dxy1, t0, t1); t0 = y0; t1 = y1;
        y0 = y1 = p6;
        MADD_SX2X2(y0, y1, dxy0, dxy1, t0, t1);
        AE_SSX2X2_IP(y0, y1, S_wr, 3 * 4 * sz_f32);
      }
      __Pragma("no_reorder");
      S_wr = (      xtfloatx4*)scr;
      S_rd = (const xtfloatx4*)scr;
      TBL_POW2 = (const xtfloat *)pow2f_coef;
      pY = (const xtfloatx4*)y;
      aY = AE_LA128_PP(pY);
      /*14/2*/
      for (n = 0; n<(blkLen >> 2); n++)
      {
        AE_LSX2X2_IP(t0, t1, S_rd, 4 * sz_f32);
        AE_LSX2X2_IP(t2, t3, S_rd, 4 * sz_f32);
        AE_LASX2X2_IP(y0, y1, aY, pY);
        /* compute y*log2(x) and separate result into integer and fractional parts */
        MUL_SX2X2(y2, y3, y0, y1, t0, t1);
        xy0 = XT_FIROUND_SX2(y2);
        xy1 = XT_FIROUND_SX2(y3);
     
        AE_LSX2X2_IP(z0, z1, S_rd, 4 * sz_f32);
        /* apply integer part */
        e0 = XT_TRUNC_SX2(xy0, 0);
        e1 = XT_TRUNC_SX2(xy1, 0);
        c7i = AE_L32_I(TBL, 7 * 4);/* -252 */
        c8i = AE_L32_X(TBL, 8 * 4);/* 254 */
        AE_MINMAX32(e0, c7i, c8i);
        AE_MINMAX32(e1, c7i, c8i);

        ex0 = AE_SRAI32(e0, 1);
        ex1 = AE_SRAI32(e1, 1);

        e0 = AE_SUB32(e0, ex0);
        e1 = AE_SUB32(e1, ex1);
        
        c0_0 = FLOATEXP_SX2(AE_MOVINT8X8_FROMINT32X2(e0));
        c0_1 = FLOATEXP_SX2(AE_MOVINT8X8_FROMINT32X2(ex0));
        c1_0 = FLOATEXP_SX2(AE_MOVINT8X8_FROMINT32X2(e1));
        c1_1 = FLOATEXP_SX2(AE_MOVINT8X8_FROMINT32X2(ex1));
        
        MUL_SX2X2(z0, z1, z0, z1, c0_1, c1_1);
        MUL_SX2X2(z0, z1, z0, z1, c0_0, c1_0);
        AE_SSX2X2_IP(z0, z1, S_wr, 4 * sz_f32);
      }
    }
    __Pragma("no_reorder");
    /******************************************/
    {
      xtbool2 b_yint, b_e0, b0, b_notspec;
      xtbool2 b_yeqz, b_yinf, b_xeqz, b_xeq1, b_xinf;
      xtbool2 b_NaN1, b_NaN2, b_one, b_Inf, b_zero;
      uint32_t b0i, b1i;
      uint32_t yeqz, yinf, xeqz, xeq1, xinf, sx, sy, yint;
      uint32_t one, NaN1, Inf, zero;
      xtfloatx2 xabs, spec;
      ae_int32x2 sgn, zi0;
      ae_valign aX, aY, aZ;

      S_rd = (const xtfloatx4*)scr;
      pY   = (const xtfloatx4*)y;
      pX   = (const xtfloatx4*)x;
      pZ   = (      xtfloatx4*)z;
      aY = AE_LA64_PP(pY);
      aX = AE_LA64_PP(pX);
      aZ = AE_ZALIGN64();
      /*33*/
      /*31*/
      for (n = 0; n<(blkLen >> 1); n++)
      {
        XT_LSX2IP(z0, castxcc(xtfloatx2, S_rd), 2 * sz_f32);
        XT_LASX2IP(x0, aX, castxcc(xtfloatx2, pX));
        XT_LASX2IP(y0, aY, castxcc(xtfloatx2, pY));
        /* Take sign of x and y */
        xi0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(x0);
        yi0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(y0);
        bsx = XT_OLT_SX2(xi0, (xtfloatx2)0.0f);
        bsy = XT_OLT_SX2(yi0, (xtfloatx2)0.0f);

        xabs = XT_ABS_SX2(x0);
        /* check if y is integer */
        t0 = XT_FITRUNC_SX2(y0);
        b_yint = XT_OEQ_SX2(t0, y0);

        /* check if y is odd */
        e0 = XT_TRUNC_SX2(y0, 0); //temp0
        b_e0 = AE_EQ32(e0, MAX_INT32);//~b_tmp0
        b0i = AE_MOVAB2(b_e0);
        b1i = AE_MOVAB2(b_yint);
        b0i = b1i&(~b0i);
        b0 = AE_MOVBA2(b0i);
        AE_MOVF32X2(e0, AE_ZERO32(), b0);
        e0 = AE_SLLI32(e0, 31);
        sgn = AE_AND32(e0, xi0);
        /* process special numbers */
        b_yeqz = XT_OEQ_SX2((xtfloatx2)0.0f, y0);            /*  y ==0      */
        b_yinf = XT_OEQ_SX2(XT_ABS_SX2(y0), plusInff.f);     /* |y|==Inf    */
        b_xeqz = XT_OEQ_SX2(x0, (xtfloatx2)0.0f);            /*  x ==0      */
        b_xeq1 = XT_OEQ_SX2(xabs, (xtfloatx2)1.0f);          /* |x|==1      */
        b_xinf = XT_OEQ_SX2(xabs, plusInff.f);               /* |x|==INF    */

        yint = AE_MOVAB2(b_yint);
        yeqz = AE_MOVAB2(b_yeqz);
        yinf = AE_MOVAB2(b_yinf);
        xeqz = AE_MOVAB2(b_xeqz);
        xeq1 = AE_MOVAB2(b_xeq1);
        xinf = AE_MOVAB2(b_xinf);
        sx = AE_MOVAB2(bsx);
        sy = AE_MOVAB2(bsy);
        one = xeq1 & (yinf | (~sx));  /* |x|==1 && ( |y|==Inf || x>0 )                       */
        one = one | yeqz;           /* ( |x|==1 && ( |y|==Inf || x>0 ) ) || y==0 --> z=1.0 */
        NaN1 = sx&(~yint);          /* x<0 && y is not an integer --> z=NaN                */
        Inf = xinf&(~sy);          /* x==INF && y>0 --> z=INF */
        Inf = Inf | (xeqz & sy);    /* x==0   && y<0 --> z=INF */
        zero = xeqz &(~sy);         /* x==0   && y>0 --> z=0.0 */
        zero = zero | (xinf & sy);  /* x==INF && y<0 --> z=0.0 */

        b_NaN1 = AE_MOVBA2(NaN1);
        b_NaN2 = XT_UN_SX2(x0, y0);         /* isnan(x) || isnan(y) --> z=NaN                      */
        b_one = AE_MOVBA2(one);
        b_Inf = AE_MOVBA2(Inf);
        b_zero = AE_MOVBA2(zero);

        /* Save special numbers and mask for special numbers */
        spec = (xtfloatx2)qNaNf.f;
        XT_MOVF_SX2(spec, half, b_NaN1);
        XT_MOVT_SX2(spec, _0, b_zero);
        XT_MOVT_SX2(spec, plusInff.f, b_Inf);
        XT_MOVT_SX2(spec, qNaNf.f, b_NaN2);
        XT_MOVT_SX2(spec, _1, b_one);

        b_notspec = XT_OEQ_SX2(spec, half);
        /* Replace result with special numbers if needed */
        XT_MOVF_SX2(z0, spec, b_notspec);
        /* Restore sign and store result */
        zi0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(z0);
        zi0 = AE_XOR32(zi0, sgn);
        z0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(zi0);
        XT_SASX2IP(z0, aZ, castxcc(xtfloatx2, pZ));
      }
      AE_SA64POS_FP(aZ, castxcc(xtfloatx2, pZ));
    }
  }
} /* mypowf() */
void vec_powf (   float32_t * restrict z, 
            const float32_t * restrict x, 
            const float32_t * restrict y, 
            int N )
{
  const int blkSize = MAX_ALLOCA_SZ / sz_f32;
  /* Allocate a fixed-size scratch area on the stack. */
  float32_t ALIGN(32) scr[blkSize];
  float32_t ALIGN(32) tmpInx[4], tmpIny[4], tmpOut[4];
  int M;
  if (N <= 0) return;
  M = N&~3;
  if (M)
  {
    mypowf(scr, z, x, y, M);
    y += M;
    x += M;
    z += M;
    N &= 3;
  }
  if (N)
  {
    // processing the tail
    int off1, off2;
    xtfloat x0, x1, x2, y0, y1, y2;
    off1 = XT_MIN(N - 1, 1) << 2;
    off2 = XT_MIN(N - 1, 2) << 2;

    x0 = XT_LSI((const xtfloat*)x, 0);
    x1 = XT_LSX((const xtfloat*)x, off1);
    x2 = XT_LSX((const xtfloat*)x, off2);
    y0 = XT_LSI((const xtfloat*)y, 0);
    y1 = XT_LSX((const xtfloat*)y, off1);
    y2 = XT_LSX((const xtfloat*)y, off2);

    XT_SSI(x0, (xtfloat*)tmpInx, 0 * sizeof(xtfloat));
    XT_SSI(x1, (xtfloat*)tmpInx, 1 * sizeof(xtfloat));
    XT_SSI(x2, (xtfloat*)tmpInx, 2 * sizeof(xtfloat));
    XT_SSI(y0, (xtfloat*)tmpIny, 0 * sizeof(xtfloat));
    XT_SSI(y1, (xtfloat*)tmpIny, 1 * sizeof(xtfloat));
    XT_SSI(y2, (xtfloat*)tmpIny, 2 * sizeof(xtfloat));

    XT_SSI(XT_CONST_S(0), (xtfloat*)tmpInx, 3 * sizeof(xtfloat));
    XT_SSI(XT_CONST_S(0), (xtfloat*)tmpIny, 3 * sizeof(xtfloat));
    mypowf(scr, tmpOut, tmpInx, tmpIny, 4);
    x0 = XT_LSI((const xtfloat*)tmpOut, 0 * sizeof(xtfloat));
    x1 = XT_LSI((const xtfloat*)tmpOut, 1 * sizeof(xtfloat));
    x2 = XT_LSI((const xtfloat*)tmpOut, 2 * sizeof(xtfloat));

    XT_SSX(x2, (xtfloat*)z, off2);
    XT_SSX(x1, (xtfloat*)z, off1);
    XT_SSI(x0, (xtfloat*)z, 0);
  }
} /* vec_powf() */
#else
#define sz_f32    (int)sizeof(float32_t)
void vec_powf(float32_t * restrict z,
  const float32_t * restrict x,
  const float32_t * restrict y,
  int N)
{

  const int blkSizef = MAX_ALLOCA_SZ / sz_f32;
  /* Allocate a fixed-size scratch area on the stack. */
  float32_t ALIGN(16) scr[blkSizef];
  /* Table of different constants used in computations */
  static const int32_t c_tbl[] =
  {
    -126,
    -150,
    (int32_t)0x007FFFFF,/* max denormalized floating-point number / mantissa mask */
    (int32_t)0x4B800000,/* 2^24 */
    (int32_t)0x3F3504F3,/* sqrt(0.5) */
    (int32_t)0x3F000000,/*  0.5 */
    (int32_t)0xBF000000,/* -0.5 */
    -252,
    254
  };
  int n;
  const xtfloat     *          pX;
  const xtfloat     *          pY;

  const xtfloat     * restrict S_rd;
  xtfloat     * restrict S_wr;
  xtfloat     * restrict pZ;
  const ae_int32      * restrict TBL;
  const  xtfloat      * restrict TBL_LOG2;
  const  xtfloat      * restrict TBL_POW2;
  xtfloat x0, y0, z0, t0, t1, ef0;
  xtfloat c2f, c3f, c4f;
  xtfloat _0, _1, half;
  ae_int32x2 c0i, c1i, c5i, c6i, c7i, c8i;
  ae_int32 e0, xi0, yi0, ex0;
  xtbool bsx, bsy, bdenorm, bsmall;

  /* overall number of blocks; number of values in the current block */
  int blkLen;
  /* Block size, blkLen <= blkSize */
  const int blkSize = MAX_ALLOCA_SZ / (3 * sz_f32);


  if (N <= 0) return;

  NASSERT_ALIGN16(scr);

  /*
  * Data are processed in blocks of scratch area size. Further, the algorithm
  * implementation is splitted in order to feed the optimizing compiler with a
  * few loops of managable size.
  */

  blkLen = 0;
  TBL = (const ae_int32 *)c_tbl;
  for (; N>0; N -= blkLen, x += blkSize, y += blkSize, z += blkSize)
  {
    blkLen = XT_MIN(N, blkSize);
    _0 = 0.0f;
    _1 = (1.0f);
    half = (0.5f);
    {
      pX = (const xtfloat*)x;
      S_wr = (xtfloat*)scr;

      for (n = 0; n<(blkLen); n++)
      {
        XT_LSIP(x0, pX, sz_f32);

        x0 = XT_ABS_S(x0);
        c0i = AE_L32_I(TBL, 0 * 4); /*-126*/
        c1i = AE_L32_I(TBL, 1 * 4); /*-150*/
        c2f = XT_LSI((xtfloat*)TBL, 2 * 4);
        c3f = XT_LSI((xtfloat*)TBL, 3 * 4);
        /* process denormalized values */
        bdenorm = XT_OLE_S(x0, c2f);
        t0 = XT_MUL_S(x0, c3f);
        XT_MOVT_S(x0, t0, bdenorm);
        e0 = c0i;

        AE_MOVT_32(e0, c1i, bdenorm);
        /* extract exponent */
        xi0 = XT_RFR(x0);
        ex0 = AE_SRLI32(xi0, 23);
        e0 = AE_ADD32(e0, ex0);
        /* extract mantissa */
        ex0 = XT_RFR(c2f);/* load mantissa mask */ //!!!!!!!!!!!!!
        c5i = AE_L32_I(TBL, 5 * 4);/*  0.5 */
        xi0 = AE_AND32(xi0, ex0);
        xi0 = AE_OR32(xi0, c5i);
        x0 = XT_WFR(xi0);
        /* adjust the mantissa to range [ sqrt(0.5) ; sqrt(2.0) ) */
        c4f = XT_LSI((xtfloat*)TBL, 4 * 4);
        bsmall = XT_OLT_S(x0, c4f);
        t0 = XT_ADD_S(x0, x0);
        ex0 = AE_SUB32(e0, 1);
        XT_MOVT_S(x0, t0, bsmall);
        AE_MOVT_32(e0, ex0, bsmall);
        x0 = XT_SUB_S(_1, x0); //!!!
        ef0 = XT_FLOAT_S(e0, 0); //!!!
        XT_SSIP(x0, S_wr, sz_f32);
        XT_SSIP(ef0, S_wr, 2 * sz_f32);

      }
    }
    __Pragma("no_reorder");
    /*****************************************/
    {
      xtfloat p0, p1, p2, p3, p4, p5, p6, p7, p8, p9;
      xtfloat p10, p11, p12, p13;
      xtfloat t2, w0, w1;
      S_wr = (xtfloat*)scr + 2;
      S_rd = (const xtfloat*)scr;
      TBL_LOG2 = (const xtfloat *)log2f_coef;

      for (n = 0; n<(blkLen); n++)
      {
        XT_LSIP(x0, S_rd, 3 * sz_f32);

        /* evaluate polynomial approximation */
        /* Load table of coefficients */

        p0 = XT_LSI(TBL_LOG2, 0 * 4);
        p1 = XT_LSI(TBL_LOG2, 1 * 4);
        p2 = XT_LSI(TBL_LOG2, 2 * 4);
        p3 = XT_LSI(TBL_LOG2, 3 * 4);
        p4 = XT_LSI(TBL_LOG2, 4 * 4);
        p5 = XT_LSI(TBL_LOG2, 5 * 4);
        p6 = XT_LSI(TBL_LOG2, 6 * 4);
        p7 = XT_LSI(TBL_LOG2, 7 * 4);
        p8 = XT_LSX(TBL_LOG2, 8 * 4);
        p9 = XT_LSX(TBL_LOG2, 9 * 4);

        XT_MADD_S(p1, x0, p0);
        XT_MADD_S(p2, x0, p1);
        XT_MADD_S(p3, x0, p2);
        XT_MADD_S(p4, x0, p3);
        XT_MADD_S(p5, x0, p4);
        XT_MADD_S(p6, x0, p5);
        XT_MADD_S(p7, x0, p6);
        XT_MADD_S(p8, x0, p7);
        XT_MADD_S(p9, x0, p8);
        t2 = p9;
        XT_SSIP(t2, S_wr, 3 * sz_f32);
      }
      S_wr = (xtfloat*)scr;
      S_rd = (const xtfloat*)scr;

      for (n = 0; n<(blkLen); n++)
      {
        p10 = XT_LSX(TBL_LOG2, 10 * 4);
        p11 = XT_LSX(TBL_LOG2, 11 * 4);
        p12 = XT_LSX(TBL_LOG2, 12 * 4);
        p13 = XT_LSX(TBL_LOG2, 13 * 4);

        XT_LSIP(x0, S_rd, sz_f32);
        XT_LSIP(ef0, S_rd, sz_f32);
        XT_LSIP(t2, S_rd, sz_f32);

        /* next coefficients are computed in extended precision */
        t0 = XT_MUL_S(x0, t2); t1 = t0;
        XT_MSUB_S(t1, x0, t2);
        w0 = XT_ADD_S(t0, p10);
        w1 = XT_SUB_S(w0, p10);
        w1 = XT_SUB_S(t0, w1);
        w1 = XT_SUB_S(w1, t1);
        t0 = w0; t1 = w1;
        w0 = XT_MUL_S(x0, t0); w1 = w0;
        XT_MSUB_S(w1, x0, t0); t0 = w0;
        XT_MSUB_S(w1, x0, t1); t1 = w1;
        w0 = XT_ADD_S(t0, p11);
        w1 = XT_SUB_S(w0, p11);
        w1 = XT_SUB_S(t0, w1);
        w1 = XT_SUB_S(w1, t1);
        t0 = w0; t1 = w1;
        x0 = XT_NEG_S(x0);
        w0 = XT_MUL_S(x0, t0); w1 = w0;
        XT_MSUB_S(w1, x0, t0); t0 = w0;
        XT_MSUB_S(w1, x0, t1); t1 = w1;
        /* multiply by log2(e) */
        w0 = XT_MUL_S(t0, p12); w1 = w0;
        XT_MSUB_S(w1, t0, p12);
        XT_MADD_S(w1, t1, p12);
        XT_MSUB_S(w1, t0, p13);
        t0 = w0; t1 = w1;
        /* add exponent */
        w0 = XT_ADD_S(t0, ef0);
        w1 = XT_SUB_S(w0, ef0);
        w1 = XT_SUB_S(t0, w1);
        t1 = XT_SUB_S(w1, t1);//!!!!
        t0 = w0; ///!!!!!
        XT_SSIP(t0, S_wr, sz_f32);
        XT_SSIP(t1, S_wr, sz_f32);
      }
    }
    __Pragma("no_reorder");
    /*******************************************/
    {
      xtfloat xy, dxy, c0, c1, _m1;;
      xtfloat p0, p1, p2, p3, p4, p5, p6;
      S_wr = (xtfloat*)scr;
      S_rd = (const xtfloat*)scr;
      TBL_POW2 = (const xtfloat *)pow2f_coef;
      pY = (const xtfloat*)y;
      _m1 = -1.0f;
      for (n = 0; n<(blkLen); n++)
      {
        XT_LSIP(t0, S_rd, sz_f32);
        XT_LSIP(t1, S_rd, sz_f32);
        XT_LSIP(y0, pY, sz_f32);
        /* compute y*log2(x) and separate result into integer and fractional parts */
        xy = XT_FLOAT_S(XT_ROUND_S(XT_MUL_S(y0, t0), 0), 0);
        dxy = XT_NEG_S(xy);
        XT_MADD_S(dxy, y0, t0);
        XT_MADD_S(dxy, y0, t1);
        c5i = AE_L32_I(TBL, 5 * 4);/*  0.5 */
        c6i = AE_L32_I(TBL, 6 * 4);/*  -0.5 */
        dxy = XT_MIN_S(dxy, _1);
        dxy = XT_MAX_S(dxy, _m1);
        /* compute 2^fract */
        p0 = XT_LSI(TBL_POW2, 0 * 4);
        p1 = XT_LSI(TBL_POW2, 1 * 4);
        p2 = XT_LSI(TBL_POW2, 2 * 4);
        p3 = XT_LSI(TBL_POW2, 3 * 4);
        p4 = XT_LSI(TBL_POW2, 4 * 4);
        p5 = XT_LSI(TBL_POW2, 5 * 4);
        p6 = XT_LSI(TBL_POW2, 6 * 4);
        /* NOTE: do not change the order of computations and way of polynomial decomposition ! *///!!!!!!!!!!!!!!
        XT_MADD_S(p1, dxy, p0);
        XT_MADD_S(p2, dxy, p1);
        XT_MADD_S(p3, dxy, p2);
        XT_MADD_S(p4, dxy, p3);
        XT_MADD_S(p5, dxy, p4);
        XT_MADD_S(p6, dxy, p5);
        z0 = p6;
        /* apply integer part */
        e0 = XT_TRUNC_S(xy, 0);
        c7i = AE_L32_I(TBL, 7 * 4);/* -252 */
        c8i = AE_L32_X(TBL, 8 * 4);/* 254 */
        e0 = AE_MAX32(e0, c7i);
        e0 = AE_MIN32(e0, c8i);
        e0 = AE_ADD32(e0, c8i);
        ex0 = AE_SRAI32(e0, 1);
        e0 = AE_SUB32(e0, ex0);
        ex0 = AE_SLLI32(ex0, 23);
        e0 = AE_SLLI32(e0, 23);

        c0 = XT_WFR(e0);
        c1 = XT_WFR(ex0);
        z0 = XT_MUL_S(z0, c1);
        z0 = XT_MUL_S(z0, c0); //!!!!!!!!!!!!
        XT_SSIP(z0, S_wr, sz_f32);

      }
    }
    __Pragma("no_reorder");
    /******************************************/
    {
      xtbool b_yint, b_e0, b0, b_notspec;
      xtbool b_yeqz, b_yinf, b_xeqz, b_xeq1, b_xinf;
      xtbool b_NaN1, b_NaN2, b_one, b_Inf, b_zero;
      uint32_t b0i, b1i;
      uint32_t yeqz, yinf, xeqz, xeq1, xinf, sx, sy, yint;
      uint32_t one, NaN1, Inf, zero;
      xtfloat xabs, spec;
      ae_int32x2 sgn, zi0;

      S_rd = (const xtfloat*)scr;
      pY = (const xtfloat*)y;
      pX = (const xtfloat*)x;
      pZ = (xtfloat*)z;

      for (n = 0; n<(blkLen); n++)
      {
        XT_LSIP(z0, S_rd, sz_f32);
        XT_LSIP(x0, pX, sz_f32);
        XT_LSIP(y0, pY, sz_f32);

        /* Take sign of x and y */
        xi0 = XT_RFR(x0);
        yi0 = XT_RFR(y0);
        bsx = XT_OLT_S(x0, (xtfloat)0.0f);
        bsy = XT_OLT_S(y0, (xtfloat)0.0f);

        xabs = XT_ABS_S(x0);
        /* check if y is integer */
        {   /* validate if y is integral - all numbers bigger than 2^23 are assumed as integral */
          xtfloat t, c;
          t = XT_ABS_S((xtfloat)y0);
          c = 8388608.f;
          XT_MOVT_S(c, t, XT_ULT_S(t, 8388608.f));
          t = c;
          t0 = XT_FLOAT_S(XT_TRUNC_S(t, 0), 0);
          b_yint = XT_OEQ_S(XT_FLOAT_S(XT_TRUNC_S(t, 0), 0), t);
        }

        /* check if y is odd */
        e0 = XT_TRUNC_S(y0, 0); //temp0
        b_e0 = xtbool2_extract_0(AE_EQ32(e0, MAX_INT32));//~b_tmp0
        b0i = AE_MOVAB(b_e0);
        b1i = AE_MOVAB(b_yint);
        b0i = b1i&(~b0i);
        b0 = AE_MOVBA(b0i);
        AE_MOVF_32(e0, AE_ZERO32(), b0);
        e0 = AE_SLLI32(e0, 31);
        sgn = AE_AND32(e0, xi0);
        /* process special numbers */
        b_yeqz = XT_OEQ_S((xtfloat)0.0f, y0);            /*  y ==0      */
        b_yinf = XT_OEQ_S(XT_ABS_S(y0), plusInff.f);     /* |y|==Inf    */
        b_xeqz = XT_OEQ_S(x0, (xtfloat)0.0f);            /*  x ==0      */
        b_xeq1 = XT_OEQ_S(xabs, (xtfloat)1.0f);          /* |x|==1      */
        b_xinf = XT_OEQ_S(xabs, plusInff.f);               /* |x|==INF    */

        yint = AE_MOVAB(b_yint);
        yeqz = AE_MOVAB(b_yeqz);
        yinf = AE_MOVAB(b_yinf);
        xeqz = AE_MOVAB(b_xeqz);
        xeq1 = AE_MOVAB(b_xeq1);
        xinf = AE_MOVAB(b_xinf);
        sx = AE_MOVAB(bsx);
        sy = AE_MOVAB(bsy);
        one = xeq1 & (yinf | (~sx));  /* |x|==1 && ( |y|==Inf || x>0 )                       */
        one = one | yeqz;           /* ( |x|==1 && ( |y|==Inf || x>0 ) ) || y==0 --> z=1.0 */
        NaN1 = sx&(~yint);          /* x<0 && y is not an integer --> z=NaN                */
        Inf = xinf&(~sy);          /* x==INF && y>0 --> z=INF */
        Inf = Inf | (xeqz & sy);    /* x==0   && y<0 --> z=INF */
        zero = xeqz &(~sy);         /* x==0   && y>0 --> z=0.0 */
        zero = zero | (xinf & sy);  /* x==INF && y<0 --> z=0.0 */

        b_NaN1 = AE_MOVBA(NaN1);
        b_NaN2 = XT_UN_S(x0, y0);         /* isnan(x) || isnan(y) --> z=NaN                      */
        b_one = AE_MOVBA(one);
        b_Inf = AE_MOVBA(Inf);
        b_zero = AE_MOVBA(zero);

        /* Save special numbers and mask for special numbers */
        spec = (xtfloat)qNaNf.f;
        XT_MOVF_S(spec, half, b_NaN1);
        XT_MOVT_S(spec, _0, b_zero);
        XT_MOVT_S(spec, plusInff.f, b_Inf);
        XT_MOVT_S(spec, qNaNf.f, b_NaN2);
        XT_MOVT_S(spec, _1, b_one);

        b_notspec = XT_OEQ_S(spec, half);
        /* Replace result with special numbers if needed */
        XT_MOVF_S(z0, spec, b_notspec);
        /* Restore sign and store result */
        zi0 = XT_RFR(z0);
        zi0 = AE_XOR32(zi0, sgn);
        z0 = XT_WFR(zi0);
        XT_SSIP(z0, pZ, sz_f32);
      }
    }
  }
}
#endif
