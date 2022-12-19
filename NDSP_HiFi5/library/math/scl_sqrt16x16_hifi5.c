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
  NatureDSP Signal Processing Library. Mathematics
    Square Root
    Optimized code for HiFi5
  IntegrIT, 2006-2019
*/

#include "NatureDSP_types.h"
#include "NatureDSP_Signal_math.h"
#include "common.h"
#include "sqrtq15_tbl.h"

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
int16_t scl_sqrt16x16(int16_t x)
{
    ae_int16x4 t0,t1,t2,t3,t4,t5,t6,t7;
    ae_int16x4 x0,sh0,ix0,y0,z0;
    xtbool4 blezero0;
    AE_L16X4X2_I(t0,t1,(const ae_int16x8*)sqrtq15_tbl,0*sizeof(ae_int16x8));
    AE_L16X4X2_I(t2,t3,(const ae_int16x8*)sqrtq15_tbl,1*sizeof(ae_int16x8));
    AE_L16X4X2_I(t4,t5,(const ae_int16x8*)sqrtq15_tbl,2*sizeof(ae_int16x8));
    AE_L16X4X2_I(t6,t7,(const ae_int16x8*)sqrtq15_tbl,3*sizeof(ae_int16x8));

    x0=AE_MOVDA16(x);
    blezero0=AE_LE16(x0,0);
    // normalization
    sh0=AE_AND16(~1,AE_NSA16X4(x0));
    x0=AE_SRAV16RS(x0,AE_NEG16S(sh0));
    // split onto 8 segments
    ix0=AE_SRAI16(x0,12);
    x0 = AE_AND16(x0,4095);
    // polynomial approximation. First coefficient goes in Q14 so at the 
    // first step we have to double input
    y0=AE_SEL16X4(t0,t1, ix0);
    z0=AE_SEL16X4(t2,t3, ix0);
    y0=AE_MULFD16X16X4RAS(z0, AE_ADD16(x0,x0), AE_MOVDA16(0x7fff), y0);
    z0=AE_SEL16X4(t4,t5, ix0);
    y0=AE_MULFD16X16X4RAS(z0, x0, AE_MOVDA16(0x7fff), y0);
    z0=AE_SEL16X4(t6,t7, ix0);
    y0=AE_MULFD16X16X4RAS(z0, x0, AE_MOVDA16(0x7fff), y0);
    // final denormalization
    y0=AE_SRAV16RS(y0,AE_SRAI16(sh0,1));
    // range control
    AE_MOVT16X4(y0,0,blezero0);
    return AE_MOVAD16_0(y0);    
} /* scl_sqrt16x16() */
