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

/* DSP Library API */
/*    Code optimized for HiFi5 core */
#include "NatureDSP_Signal_math.h"
/* Common helper macros. */
#include "common_fpu.h"

#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(float32_t,scl_int2float,( int32_t x, int t ))
#else
/*
  NatureDSP Signal Processing Library. Scalar Mathematics
  Integer to float conversion
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/
/*-------------------------------------------------------------------------
  Integer to float conversion
  Routines convert integer to float and scale result up by 2^t.

  Precision: 
  f     32-bit input, floating point output

  Input:
  x[N]  input data, integer
  t     scale factor
  N     length of vector
  Output:
  y[N]  conversion results, floating point

  Restriction:
  t should be in range -126...126
-------------------------------------------------------------------------*/
float32_t scl_int2float( int32_t x, int t )
{
  /*
  * float32_t y;
  * ufloat32uint32 s;
  * ASSERT(t>=-126 && t<=126);
  * s.u = ((uint32_t)(t + 127)) << 23;
  * y=((float32_t)x);
  * y=y*s.f;
  * return y;
  */
  xtfloat f0, y0;
  ae_int32 x0;
  NASSERT(t >= -126 && t <= 126);
#if defined FLOATEXP_S
  y0 = FLOATEXP_S(AE_MOVINT8_FROMINT8X8(AE_MOVDA8(t)));
#else
  y0 = XT_WFR((uint32_t)(t + 127) << 23);
#endif
  x0 = x;
  f0 = XT_FLOAT_S(x0, 0);
  f0 = XT_MUL_S(f0, y0);
  return f0;
} /* scl_int2float() */
#endif
