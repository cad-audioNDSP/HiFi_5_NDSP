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
void vec_scale32x32_fast
(
        int32_t * restrict y,
  const int32_t * restrict x,
        int32_t            s,
            int            N
)
{
    const ae_int32x4* restrict pX;
          ae_int32x4* restrict pY;
    ae_int32x2 x0,x1,x2,x3,vcf;
    int n;
    if(N<=0) return;
    vcf = AE_MOVDA32X2(s, s);
    pX=(const ae_int32x4 *)x;
    pY=(      ae_int32x4 *)y;
    for (n=0; n<(N>>3); n++) 
    {
        AE_L32X2X2_IP(x0,x1,pX,sizeof(ae_int32x4));
        AE_L32X2X2_IP(x2,x3,pX,sizeof(ae_int32x4));
        AE_MULF2P32X4RAS(x0,x1,x0,x1,vcf,vcf);
        AE_MULF2P32X4RAS(x2,x3,x2,x3,vcf,vcf);
        AE_S32X2X2_IP(x0,x1,pY,sizeof(ae_int32x4));
        AE_S32X2X2_IP(x2,x3,pY,sizeof(ae_int32x4));
    }
    if (N&4)
    {
        AE_L32X2X2_IP(x0,x1,pX,sizeof(ae_int32x4));
        AE_MULF2P32X4RAS(x0,x1,x0,x1,vcf,vcf);
        AE_S32X2X2_IP(x0,x1,pY,sizeof(ae_int32x4));
    }
}
