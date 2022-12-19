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
  NatureDSP Signal Processing Library. Matrix Operations
    Square root 32x16   
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/
#include "NatureDSP_Signal_math.h"
#include "NatureDSP_types.h"
#include "polyrsqrtq23_tbl.h"
#include "common.h"

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
int16_t scl_sqrt32x16(int32_t x)
{
  ae_int32x2 y0, t0;
  ae_int32x2 x0, d0;
  ae_int32x2 c0, c1, c2, c3, c4;
  ae_int16x4 nsa01;
  xtbool2 lezero0;

  x0 =x;
  /* Normalize x*/
  nsa01 = AE_NSAZ32X4(x0, x0);
  nsa01 = AE_AND16(nsa01, 0xFFFE);
  x0 = AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(nsa01)));

  lezero0 = AE_LT32(x0, 0);

  c0 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 0);
  c1 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 1);
  c2 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 2);
  c3 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 3);
  c4 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 4);

  t0 = c1;
  AE_MULAFP32X2RAS(t0, x0, c0); y0 = t0;
  t0 = c2;
  AE_MULAFP32X2RAS(t0, x0, y0); y0 = t0;
  t0 = c3;
  AE_MULAFP32X2RAS(t0, x0, y0); y0 = t0;
  t0 = c4;
  AE_MULAFP32X2RAS(t0, x0, y0); y0 = t0;
  /* reiterate rsqrt */
  y0 = AE_SLAI32S(y0, 11);
  d0 = AE_MULFP32X2RAS(y0, y0);
  t0 = 0x08000000;
  AE_MULSFP32X2RAS(t0, d0, x0); d0 = t0;
  d0 = AE_MULFP32X2RAS(d0, y0);
  /* compute sqrt and rescale back */
  d0 = AE_SLAI32(d0, 3);
  y0 = AE_ADD32S(y0, d0);
  /* compute sqrt and rescale back */
  y0 = AE_MULFP32X2RAS( x0, y0);
  nsa01 = AE_SRAI16(nsa01, 1);
  nsa01 = AE_SUB16(nsa01, 2);
  y0 = AE_SRAV32RS(y0, AE_SEXT32X2D16_32(nsa01));
  AE_MOVT32X2(y0, 0x80000000, lezero0);
  return (int16_t)(AE_MOVAD32_L(y0)>>16);
} /* scl_sqrt32x16() */
