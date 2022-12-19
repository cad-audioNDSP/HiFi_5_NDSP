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
    Square Root
    Optimized code for HiFi5
  IntegrIT, 2006-2019
*/

#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_math.h"
#include "common.h"
#include "scl_sqrt_table.h"

/*-------------------------------------------------------------------------
  Square Root
  These routines calculate square root.
  NOTE: functions return 0x80000000 on negative argument for 32-bit outputs
  or 0x8000 for 16-bit outputs.
  Two versions of functions available: regular version (vec_sqrt16x16, vec_sqrt32x32, 
  vec_sqrt32x16, vec_sqrt64x32) with arbitrary 
  arguments and faster version (vec_sqrt32x32_fast) that 
  apply some restrictions.

  Precision: 
  16x16  16-bit inputs, 16-bit output. Accuracy: 2LSB
  32x32  32-bit inputs, 32-bit output. Accuracy: (2.6e-7*y+1LSB)
  32x16  32-bit input, 16-bit output.  Accuracy: 2 LSB
  64x32  64-bit inputs, 32-bit output. Accuracy: 2LSB

  Input:
  x[N]  input data, Q15, Q31, Q63 
  N     length of vectors
  Output:
  y[N]  output data, Q15, Q31

  Restriction:
  Regular versions (vec_sqrt16x16, vec_sqrt32x32, vec_sqrt32x16, vec_sqrt64x32):
  x,y - should not overlap

  Faster versions (vec_sqrt32x32_fast):
  x,y - should not overlap
  x,y - aligned on 16-byte boundary
  N   - multiple of 2

  Scalar versions:
  ----------------
  return result, Q15, Q31
-------------------------------------------------------------------------*/
int32_t scl_sqrt64x32(int64_t x)
{
  int         nsa, shft_a;
  int32_t     y;
  ae_int64    vxw;
  ae_int32x2  vdw, viw, vxx, vzw;
  ae_int32x2  p0, p1, p2, p3, p4, p5;
  ae_int32x2  t0, y0, d0, x0;
  xtbool2     inf, zero;

  viw = AE_MOVDA32X2(MIN_INT32, MIN_INT32);
  vzw = AE_MOVDA32X2(0, 0);

  {
    ae_valign   movx;
    const ae_int32x2 * px = (const ae_int32x2 *)&x;
    movx = AE_LA64_PP(px);
    AE_LA32X2_IP(vxx, movx, px);
    vxx = AE_SEL32_LH(vxx, vxx);
    vxw = AE_MOVINT64_FROMINT32X2(vxx);
  }
  nsa = AE_NSA64(vxw);
  shft_a = (nsa & ~1);
  nsa = nsa >> 1;
  vxw = AE_SLAA64(vxw, shft_a);
  x0 = AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(vxw), AE_MOVINT32X2_FROMINT64(vxw));
  zero = AE_EQ32(x0, 0);
  inf = AE_LT32(x0, 0);

  p0 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 0);
  p1 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 1);
  p2 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 2);
  p3 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 3);
  p4 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 4);
  p5 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 5);

  t0 = p0; y0 = p1;
  AE_MULAFP32X2RAS(y0, x0, t0); t0 = y0;
  y0 = p2;
  AE_MULAFP32X2RAS(y0, x0, t0); t0 = y0;
  y0 = p3;
  AE_MULAFP32X2RAS(y0, x0, t0); t0 = y0;
  y0 = p4;
  AE_MULAFP32X2RAS(y0, x0, t0); t0 = y0;
  y0 = p5;
  AE_MULAFP32X2RAS(y0, x0, t0); t0 = y0;
  t0 = AE_SLAI32S(t0, 4);

  y0 = 0x08000000;
  d0 = AE_MULFP32X2RAS(t0, t0);
  AE_MULSFP32X2RAS(y0, x0, d0); d0 = y0;
  d0 = AE_MULFP32X2RAS(d0, t0);

  d0 = AE_SLAI32S(d0, 3);
  t0 = AE_ADD32S(d0, t0);
  
  /* compute sqrt and rescale back */
  y0 = AE_MULFP32X2RAS(t0, x0);
  y0 = AE_SLAI32S(y0, 2);

  d0 = x0; 
  AE_MULSFP32X2RAS(d0, y0, y0);
  d0 = AE_ADD32S(d0, d0);
  AE_MULAFP32X2RAS(y0, d0, t0);

  vdw = AE_SRAV32RS(y0, nsa);
  AE_MOVT32X2(vdw, viw, inf);
  AE_MOVT32X2(vdw, vzw, zero);
  y = AE_MOVAD32_H(vdw);
  return y;
} /* scl_sqrt64x32() */
