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
    Code optimized for HiFi5 core
  IntegrIT, 2006-2019
*/
#include <errno.h>
#include <fenv.h>
 
#include "NatureDSP_Signal_math.h"
#include "NatureDSP_types.h"
#include "common.h"
#include "common_fpu.h"

/* Tables and constants. */
#include "tanh_fp16_tbl.h"


/* If non-zero, set errno and raise floating-point exceptions on errors. */

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
#if !HAVE_HPFPU
DISCARD_FUN(float16_t,scl_tanh_fp16,(float16_t x))
#else
float16_t scl_tanh_fp16( float16_t x )
{
#if 1
  int32_t SCF; /* Floating-point Status and Control Register values. */

  const xthalf * restrict pT = (const xthalf  *)tanh_fp16_tbl;
  xthalfx4 halfln3_fp16 = AE_MOVXTHALFX4_FROMINT16X4(AE_L16_I((ae_int16 *)pT, 0 * 2));
  xthalfx4 tanhf16_maxval = AE_MOVXTHALFX4_FROMINT16X4(AE_L16_I((ae_int16 *)pT, 1 * 2));
  ae_int16x4 _23637 = AE_L16_X((ae_int16 *)pT, 8 * 2);
  ae_int16x4 _7 = AE_L16_X((ae_int16 *)pT, 9 * 2);

  xthalfx4 x0, z0, zz0;
  ae_int16x4 a0, n0, a2;
  ae_int16x4 y0, y1, y2, y3, t0;
  ae_int32x2 d0, d1;

   x0 = AE_MOVXTHALFX4_FROMINT16X4(AE_L16_I((ae_int16 *)(&x), 0));
   if (AE_MOVAB4(UN_HX4(x0, x0)))
   {
     __Pragma("frequency_hint never");
     errno = EDOM;
     return 0x7e00;
   }
   SCF = XT_RUR_FSR(); /* Sample floating-point exception flags. */
   
   x0= ABS_HX4(x0);
   /* exact formula for x>0.54 */
   x0 = MINNUM_HX4(x0, tanhf16_maxval);
   a0 = TRUNC16_HX4(x0, 12);
   /* multiply by 1/ln(2) and convert to Q15 */
   AE_MUL16X4(d0, d1, a0, _23637);
   d0 = AE_SRAI32R(d0, 10);
   d1 = AE_SRAI32R(d1, 10);
   /* exponential part  */
   n0 = AE_TRUNCI16X4F32S(d0, d1, 1);
   d0 = AE_AND32(d0, 0x7fff);
   d1 = AE_AND32(d1, 0x7fff);
   /* mantissa          */
   a0 = AE_TRUNCI16X4F32S(d0, d1, 16);
   /* compute 2^a, 0..1 in Q14 */
   a2 = AE_MULFP16X4S(a0, a0);
   y1 = AE_L16_I((ae_int16 *)pT, 4 * 2);
   y2 = AE_L16_I((ae_int16 *)pT, 5 * 2);
   y3 = AE_L16_I((ae_int16 *)pT, 6 * 2);
   t0 = AE_MULFP16X4S(a2, y1);
   y0 = AE_ADD16S(y3, t0); 
  
   y3 = AE_L16_I((ae_int16 *)pT, 7 * 2);
   t0 = AE_MULFP16X4S(a2, y2);
   y2 = AE_ADD16S(y3, t0); 
  
   t0 = AE_MULFP16X4S(a0, y0);    y0 = AE_ADD16S(y2, t0);
  
   z0 = FLOAT16_HX4(y0, 7);
   n0 = AE_SLLI16S(AE_ADD16S(n0, _7), 10);
   z0 = MUL_HX4(z0, AE_MOVXTHALFX4_FROMINT16X4(n0));
   /* convert exp(2x)/2 to tanh */
   z0 = ADD_HX4( z0, CONST_HX4(3));
   z0 = RECIP_HX4(z0);
   z0 = SUB_HX4( CONST_HX4(1), z0);
   zz0=z0; 
  {
    xthalfx4 x0, x1, y0, z0, z1;
    xthalfx4 x2, t0;
    xtbool4 b0;
    ae_int16x4 a0;

    x0 = x1 = AE_MOVXTHALFX4_FROMINT16X4(AE_L16_I((ae_int16 *)(&x), 0));
    y0 = zz0; 
    a0 = AE_MOVINT16X4_FROMXTHALFX4(x0);
    ABS_HX4X2(x0, x1, x0, x1); 
    
    b0 = ULT_HX4(x0, halfln3_fp16);
    /* polynomial for small x */
    x2 = MUL_HX4( x0, x0);
    z1 = AE_MOVXTHALFX4_FROMINT16X4(AE_L16_I((ae_int16 *)pT, 2 * 2));
    t0 =  AE_MOVXTHALFX4_FROMINT16X4(AE_L16_I((ae_int16 *)pT, 3 * 2));
    MADD_HX4(t0,  x2,  z1); z0 = t0; 
    /// z0 = MUL_HX4( x2, z0);
    MUL_HX4X2(z0, z0, x2, x2, z0, z0);
     t0 = x0; MADD_HX4(t0, z0, x0); z0 = t0; 
    /* select variant and apply sign */
    MOVT_HX4(y0, z0, b0);  
    MOVT_HX4(y0, NEG_HX4(y0), AE_LT16(a0, AE_ZERO16()));
    XT_WUR_FSR(SCF);

    return AE_MOVAD16_0(AE_MOVINT16X4_FROMXTHALFX4(y0));
  }
#else
  const xthalf * restrict pT = (const xthalf  *)tanh_fp16_tbl;
  int32_t SCF; /* Floating-point Status and Control Register values. */
  xthalf halfln3_fp16   = AE_LHI(pT, 0*2);
  xthalf tanhf16_maxval = AE_LHI(pT, 1*2);
  int16_t _23637        = AE_L16_X((ae_int16 *)pT, 8 * 2);
  int16_t _7            = AE_L16_X((ae_int16 *)pT, 9 * 2);
  int16_t sx, a, n;
  ae_int16x4 y;
  xthalf x0, x2;
  xthalfx4 zbig, z;
  xthalf t, zsmall;
  xtbool4 bsmall;
  x0 = AE_LHI((xthalf *)(&x), 0);
  __Pragma("no_reorder");
  if (AE_MOVAB4(UN_H(x0, x0)))
  {
    __Pragma("frequency_hint never");
    errno = EDOM;
    return 0x7e00;
  }
  __Pragma("no_reorder");
  SCF = XT_RUR_FSR(); /* Sample floating-point exception flags. */
  sx = (x);
  x0 = ABS_H(x0);
  bsmall = ULT_H(x0, halfln3_fp16);
  /* polynomial for small x */
  x2 = MUL_H(x0, x0);
  zsmall = AE_LHX(pT, 2 * 2); 
  t = AE_LHX(pT, 3 * 2);  MADD_H(t, x2, zsmall); zsmall = t;
  __Pragma("no_reorder");
  zsmall = MUL_H(x2, zsmall);
  {
    xthalfx4 h0,h1;
    h0 = zsmall; h1 = x2;
    MUL_HX4X2(h0, h0, h1, h1, h0, h0);
    zsmall = h0;
  }
  __Pragma("no_reorder");
 // MUL_HX4X2(zsmall, zsmall, x2, x2, zsmall, zsmall);
  t = x0; MADD_H(t, zsmall, x0); zsmall = t;
  __Pragma("no_reorder");
  /* exact formula for x>0.54 */
  /* exp(2x)/2, x>=0 */
  x0 = MINNUM_H(x0, tanhf16_maxval);
  a = TRUNC16_H(x0, 12);
  /* multiply by 1/ln(2) and convert to Q15 */
  {
    ae_int32x2 tmp;
    tmp = AE_MUL16S_00(a, _23637);
    __Pragma("no_reorder");
    tmp = AE_SRAI32R(tmp, 10); 
    n = AE_TRUNCI16X4F32S(tmp, tmp, 1);  /* exponential part  */
    tmp = AE_AND32(tmp, 0x7fff);
    a = AE_TRUNCI16X4F32S(tmp, tmp, 16);  /* mantissa          */
  }
  /* compute 2^a, 0..1 in Q14 */
  { 
    ae_f16x4 y1, y2, c1, c2, a2, t;
    a2 = AE_MULFP16X4S(a, a);
    __Pragma("no_reorder");
    y1 = AE_L16_X((ae_int16 *)pT, 4 * 2);
    y2 = AE_L16_X((ae_int16 *)pT, 5 * 2);
    c1 = AE_L16_X((ae_int16 *)pT, 6 * 2); 
    t = AE_MULFP16X4S(a2, y1);  y1 = AE_ADD16S(c1, t);
    c2 = AE_L16_X((ae_int16 *)pT, 7 * 2); t = AE_MULFP16X4S(a2, y2);  y2 = AE_ADD16S(c2, t);
    t = AE_MULFP16X4S(a, y1);  y = AE_ADD16S(y2, t);
  }
  zbig = FLOAT16_HX4(y, 7);
  __Pragma("no_reorder");
  n = AE_SLLI16S(AE_ADD16S(n, _7), 10);
  __Pragma("no_reorder");
  zbig = MUL_HX4(zbig, AE_MOVXTHALFX4_FROMINT16X4(AE_MOVDA16(n)));
  /* convert exp(2x)/2 to tanh */
  zbig = ADD_HX4(zbig, CONST_HX4(3));
  zbig = RECIP_HX4(zbig);
  zbig = SUB_HX4(CONST_HX4(1), zbig);
  /* select small or big variant and apply original sign */
  z = zbig;
  zbig = AE_MOVXTHALFX4_FROMINT16X4(AE_L16_I((ae_int16 *)(&zsmall), 0));
  MOVT_HX4(z, zbig, bsmall);
  __Pragma("no_reorder");
  y = AE_MOVINT16X4_FROMXTHALFX4(z);
  AE_MOVT16X4(y, AE_MOVINT16X4_FROMXTHALFX4(NEG_HX4(z)), AE_LT16(sx, AE_ZERO16()));
  XT_WUR_FSR(SCF);
  return AE_MOVAD16_0(y);
  #endif
} /* scl_tanh_fp16() */
#endif

