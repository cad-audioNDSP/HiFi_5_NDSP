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
  NatureDSP Signal Processing Library. Vector Operations
    Vector Scaling with Saturation
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/
/* Signal Processing Library API. */
#include "NatureDSP_Signal_vector.h"
#include "common.h"
#include "common_fpu.h"

#if (HAVE_VFPU==0 && HAVE_FPU==0)
DISCARD_FUN(void,vec_scalef,(     float32_t * restrict y,
                    const float32_t * restrict x,
                    float32_t s,
                    int N))
#elif (HAVE_VFPU)
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
#define sz_f32    (int)sizeof(float32_t)
void vec_scalef     (     float32_t * restrict y,
                    const float32_t * restrict x,
                    float32_t s,
                    int N)
{
    const xtfloatx4* restrict pX0;
    const xtfloatx4* restrict pX1;
          xtfloatx4* restrict pY0;
          xtfloatx4* restrict pY1;
    xtfloatx2 x0,x1,x2,x3,vcf;
    ae_valignx2 ax0,ay0;
    ae_valignx2 ax1,ay1;
    int n,N0;
    if(N<=0) return;
    vcf = s;
    pX0=(const xtfloatx4 *)x;
    pY0=(      xtfloatx4 *)y;
    if (N<=7)
    {
        __Pragma("no_unroll")
        for (n=0; n<N; n++) 
        {
            xtfloat x0;
            XT_LSIP(x0,castxcc(xtfloat,pX0),sizeof(float32_t));
            x0=XT_MUL_S(x0,vcf);
            XT_SSIP(x0,castxcc(xtfloat,pY0),sizeof(float32_t));
        }
        return ;
    }
    N0=((N-1)&7)+1;
    ax0=AE_LA128_PP(pX0);
    ay0=AE_ZALIGN128();
    ay1=AE_ZALIGN128();
    AE_LASX2X2_IP(x0,x1,ax0,pX0);
    AE_LASX2X2_IP(x2,x3,ax0,pX0);
    MULQ_S(x0,x1,x0,x1,vcf);
    MULQ_S(x2,x3,x2,x3,vcf);
    AE_SASX2X2_IP(x0,x1,ay0,pY0);
    AE_SASX2X2_IP(x2,x3,ay0,pY0);
    AE_SA128POS_FP(ay0,pY0);
    x+=N0;
    y+=N0;
    N-=N0;
    pX0=(const xtfloatx4 *)(x);
    pY0=(      xtfloatx4 *)(y);
    pX1=(const xtfloatx4 *)(x+(N>>1));
    pY1=(      xtfloatx4 *)(y+(N>>1));
    ax0=AE_LA128_PP(pX0);
    ax1=AE_LA128_PP(pX1);
    for (n=0; n<(N>>3); n++) 
    {
        AE_LASX2X2_IP(x0,x1,ax0,pX0);
        AE_LASX2X2_IP(x2,x3,ax1,pX1);
        MULQ_S(x0,x1,x0,x1,vcf);
        MULQ_S(x2,x3,x2,x3,vcf);
        AE_SASX2X2_IP(x0,x1,ay0,pY0);
        AE_SASX2X2_IP(x2,x3,ay1,pY1);
    }
    AE_SA128POS_FP(ay0,pY0);
    AE_SA128POS_FP(ay1,pY1);
}
#else
void vec_scalef     (     float32_t * restrict y,
                    const float32_t * restrict x,
                    float32_t s,
                    int N)
{
          xtfloat  * restrict pY = (      xtfloat  *)y;
    const xtfloat  * restrict pX = (const xtfloat  *)x;
    int n;
    for (n=0; n<N; n++)
    {
        xtfloat x0;
        XT_LSIP(x0, pX, sizeof(xtfloat));
        x0=XT_MUL_S(x0,s);
        XT_SSIP(x0, pY, sizeof(xtfloat));
    }
}
#endif
