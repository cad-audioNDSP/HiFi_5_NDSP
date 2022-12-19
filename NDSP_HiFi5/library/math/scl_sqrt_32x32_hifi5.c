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
#include "sqrt_table.h"
#include "common.h"
#include "NatureDSP_Signal_math.h"

/*
  NatureDSP Signal Processing Library. Scalar Mathematics
    Square Root
    Optimized code for HiFi5
  IntegrIT, 2006-2019
*/

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
int32_t scl_sqrt32x32(int32_t x) 
#if 0
{
  int         nsa, idx, shft_a;
  int32_t     y;
  ae_int32x2  vbw, vcw, vxw, zero, vdw, min_int32;
  ae_f32x2    vxf, vzf, vcf;
  xtbool2     inf, is_zero;
       
  zero = AE_MOVI(0);
  vxw = AE_MOVDA32X2(x,x);
  min_int32 = AE_MOVDA32X2(MIN_INT32, MIN_INT32);
  is_zero = AE_EQ32(vxw, zero);
  nsa = AE_NSAZ32_L(vxw);
  shft_a = ( nsa & ~1 );
  nsa   = ( nsa )/2;
  vxw = AE_SLAA32( vxw, shft_a );
  inf = AE_LT32(vxw, zero);

  /* get table indexes */
  vdw = AE_SRAI32(vxw, 24);
  vdw = AE_SLAI32(vdw, 3);
  idx=AE_MOVAD32_H(vdw); //
  vcw = AE_L32X2_X((const ae_int32x2 *)sqrt_table, idx); 
  vbw = AE_SEL32_LL(vcw, vcw); //f'(x0)
  vzf = (vcw);
  vxf = (vxw);
  AE_MULSFP32X2RAS(vxf, vzf, vzf); //dx
  vcf = (vbw);
  AE_MULAFP32X2RAS(vzf, vcf, vxf); //y 

  /*second iteration*/
  vxf = (vxw);
  AE_MULSFP32X2RAS(vxf, vzf, vzf);
  vcf = (vbw);
  AE_MULAFP32X2RAS(vzf, vcf, vxf); //y
  vdw =  (vzf);
  vdw = AE_SRAA32(vdw, nsa);
  AE_MOVT32X2(vdw, min_int32, inf);
  AE_MOVT32X2(vdw, zero, is_zero);
  y = AE_MOVAD32_H(vdw);
  return y;
}
#else
{
    ae_int32x2 x0;
    ae_int16x4 nsa01;
    ae_int16x4 y0,z0,w0,d0;
    ae_int32x2 t0;
    xtbool2 xnotpos0;
    ae_int32x2 c0;

    x0=AE_MOVDA32(x);
    // set masks for outbound conditions
    xnotpos0=AE_LE32(x0,0);
    c0=MIN_INT32; AE_MOVT32X2(c0,0,AE_EQ32(x0,0));
    // normalization
    nsa01=AE_AND16(AE_MOVDA16(~1),AE_NSA32X4(x0,x0));
    x0=AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(nsa01)));
    y0=AE_TRUNC16X4F32(x0,x0);
    // first approximation
    w0=AE_ADD16(AE_MOVDA16(0x8000),AE_SRAI16(y0,1));
    z0=AE_MULFD16X16X4RAS(AE_MOVDA16(21561),AE_MOVDA16(11534),y0,AE_MOVDA16(0x7fff));
    // first iteration
    d0=AE_MULFD16X16X4RAS(w0,AE_MOVDA16(0x4000),z0,AE_MOVDA16(0x7fff));
    d0=AE_ADD16(d0,d0); 
    w0=AE_MULFD16X16X4RAS(w0,w0,AE_MOVDA16(0x7fff),d0);
    d0=AE_MULFD16X16X4RAS(z0,y0,z0,AE_MOVDA16(MIN_INT16));
    z0=AE_MULFD16X16X4RAS(z0,w0,AE_MOVDA16(0x7fff),d0);
    // second iteration
    d0=AE_MULFD16X16X4RAS(w0,AE_MOVDA16(0x4000),z0,AE_MOVDA16(0x7fff));
    d0=AE_ADD16(d0,d0); 
    w0=AE_MULFD16X16X4RAS(w0,w0,AE_MOVDA16(0x7fff),d0);
    d0=AE_MULFD16X16X4RAS(z0,y0,z0,AE_MOVDA16(MIN_INT16));
    z0=AE_MULFD16X16X4RAS(z0,w0,AE_MOVDA16(0x7fff),d0);
    // third short iteration
    AE_CVTI32X4F16(t0,t0,z0,16);
    AE_MULSF16X4SS(x0,x0,z0,z0);
    AE_MULSFP32X16X2RAS_H(t0,x0,w0);
    // output scaling
    nsa01=AE_SRAI16(nsa01,1);
    t0=AE_SRAV32RS(t0, AE_SEXT32X2D16_32(nsa01));
    // masking
    AE_MOVT32X2(t0,c0,xnotpos0);
    return AE_MOVAD32_H(t0);
} /* scl_sqrt32x32() */
#endif
