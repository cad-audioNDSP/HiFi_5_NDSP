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

/* Library API */
#include "NatureDSP_Signal_vector.h"
#include "common.h"
// code optimized for HiFi5 core

/*===========================================================================
  Vector matematics:
  vec_shift,vec_scale  Vector Scaling with Saturation
===========================================================================*/

/*-------------------------------------------------------------------------
  Vector Scaling with Saturation
  These routines make shift with saturation of data values in the vector 
  by given scale factor (degree of 2).
  Functions vec_scale() make multiplication of a vector to a coefficient 
  which is not a power of 2 forming Q31, Q15 or floating-point result.
  Two versions of routines are available: regular versions (vec_shift32x32, 
  vec_shift16x16, vec_shiftf, vec_scale32x32, vec_scale16x16, vec_scalef, 
  vec_scale_sf) work with arbitrary arguments, faster versions 
  (vec_shift32x32_fast, vec_shift16x16_fast, vec_scale32x32_fast, 
  vec_scale16x16_fast) apply some restrictions

  For floating point:
  Fuction vec_shiftf() makes scaling without saturation of data values by given 
  scale factor (degree of 2). 
  Functions vec_scalef() and vec_scale_sf() make multiplication of input vector
  to coefficient which is not a power of 2.
  Two versions of routines are available: 
    without saturation - vec_scalef;
    with saturation - vec_scale_sf; 

  Precision:
  32x32 32-bit input, 32-bit output
  16x16 16-bit input, 16-bit output
  f     single precision floating point

  Input:
  x[N]    input data, Q31, Q15 or floating point
  t       shift count. If positive, it shifts left with saturation, if
          negative it shifts right
  s       scale factor, Q31, Q15 or floating point
  N       length of vector
  fmin    minimum output value (only for vec_scale_sf)
  fmax    maximum output value (only for vec_scale_sf)
  Output:
  y[N]    output data, Q31, Q15 or floating point

  Restrictions:
  For regular versions (vec_shift32x32, vec_shift16x16, vec_shiftf, 
  vec_scale32x32, vec_scale16x16, vec_scalef, vec_scale_sf):
  x,y should not overlap
  t   should be in range -31...31 for fixed-point functions and -129...146 
      for floating point
  For vec_scale_sf:
  fmin<=fmax;

  For faster versions (vec_shift32x32_fast, vec_shift16x16_fast, 
  vec_scale32x32_fast,vec_scale16x16_fast):
  x,y should not overlap
  t should be in range -31...31 
  x,y - aligned on 16-byte boundary
  N   - multiple of 4 
-------------------------------------------------------------------------*/
void vec_shift16x16_fast
(
        int16_t * restrict  y,
  const int16_t * restrict  x,
  int                       t,
  int                       N
)
{
    const ae_int16x8* restrict pX;
          ae_int16x8* restrict pZ;
    ae_int16x4 x0,x1,x2,x3;
    ae_int16x4 scale;
    int n;
    if(N<=0) return;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    pX=(const ae_int16x8 *)x;
    pZ=(      ae_int16x8 *)y;
    if (t<0)
    {   // negative shifts: execute via fractional multiplies
        t = -t;
        scale = AE_MOVDA16 ( (t>15) ? 1: 0x8000>>t);

        for (n=0; n<(N>>4); n++) 
        {
            AE_L16X4X2_IP(x0,x1,pX,sizeof(ae_int16x8));
            AE_L16X4X2_IP(x2,x3,pX,sizeof(ae_int16x8));
            x0=AE_MULFP16X4RS(x0,scale);
            x1=AE_MULFP16X4RS(x1,scale);
            x2=AE_MULFP16X4RS(x2,scale);
            x3=AE_MULFP16X4RS(x3,scale);
            AE_S16X4X2_IP(x0,x1,pZ,sizeof(ae_int16x8));
            AE_S16X4X2_IP(x2,x3,pZ,sizeof(ae_int16x8));
        }
        N&=15;
        if (N)
        {
            __Pragma("no_unroll")
            for (n=0; n<(N>>2); n++) 
            {
                AE_L16X4_IP(x0,castxcc(ae_int16x4,pX),sizeof(ae_int16x4));
                x0=AE_MULFP16X4RS(x0,scale);
                AE_S16X4_IP(x0,castxcc(ae_int16x4,pZ),sizeof(ae_int16x4));
            }
        }
        return;
    }
    // positive shifts
    scale = AE_SLAA16S(AE_MOVDA16(1), t);
    for (n=0; n<(N>>4); n++) 
    {
        AE_L16X4X2_IP(x0,x1,pX,sizeof(ae_int16x8));
        AE_L16X4X2_IP(x2,x3,pX,sizeof(ae_int16x8));
        x0=AE_MULP16X16X4S(x0,scale);
        x1=AE_MULP16X16X4S(x1,scale);
        x2=AE_MULP16X16X4S(x2,scale);
        x3=AE_MULP16X16X4S(x3,scale);
        AE_S16X4X2_IP(x0,x1,pZ,sizeof(ae_int16x8));
        AE_S16X4X2_IP(x2,x3,pZ,sizeof(ae_int16x8));
    }
    N&=15;
    if (N)
    {
        __Pragma("no_unroll")
        for (n=0; n<(N>>2); n++) 
        {
            AE_L16X4_IP(x0,castxcc(ae_int16x4,pX),sizeof(ae_int16x4));
            x0=AE_MULP16X16X4S(x0,scale);
            AE_S16X4_IP(x0,castxcc(ae_int16x4,pZ),sizeof(ae_int16x4));
        }
    }
}
