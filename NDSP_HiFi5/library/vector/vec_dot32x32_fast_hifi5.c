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
  vec_dot              Vector Dot Product
===========================================================================*/

/*-------------------------------------------------------------------------
  Vector Dot product
  These routines take two vectors and calculates their dot product.
  Two versions of routines are available: regular versions (vec_dot64x32,
  vec_dot64x64, vec_dot64x64i, vec_dot32x16, vec_dot32x32,vec_dot16x16, 
  vec_dotf) work with arbitrary arguments, faster versions (vec_dot64x32_fast, 
  vec_dot64x64_fast, vec_dot64x64i_fast, vec_dot32x16_fast, vec_dot32x32_fast,
  vec_dot16x16_fast) apply some restrictions.  
  NOTE:
  vec_dot16x16_fast utilizes 32-bit saturating accumulator, so input data 
  should be scaled properly to avoid erroneous results.

  Precision: 
  64x32  64x32-bit data, 64-bit output (fractional multiply Q63xQ31->Q63)
  64x64  64x64-bit data, 64-bit output (fractional multiply Q63xQ63->Q63)
  64x64i 64x64-bit data, 64-bit output (low 64 bit of integer multiply)
  32x32  32x32-bit data, 64-bit output
  32x16  32x16-bit data, 64-bit output
  16x16  16x16-bit data, 64-bit output for regular version and 32-bit for 
                        fast version
  f      single precision floating point

  Input:
  x[N]  input data, Q15, Q31, Q63 or floating point
  y[N]  input data, Q15, Q31, Q63 or floating point
  N	    length of vectors
  Returns:
  dot product of all data pairs, Q31, Q63 or floating point

  Restrictions:
  Regular versions:
    none
  Faster versions:
    x,y - aligned on 16-byte boundary
    N   - multiple of 4
-------------------------------------------------------------------------*/
int64_t vec_dot32x32_fast (const int32_t * restrict x,const int32_t * restrict y,int N)
{
    const ae_int32x4* restrict pX;
    const ae_int32x4* restrict pY;
    ae_int32x2 x0,x1,x2,x3;
    ae_int32x2 y0,y1,y2,y3;
    ae_int64 a0,a1,a2,a3;
    int n;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(N%4==0);
    if(N<=0) return 0;
    a0=a1=a2=a3=0;
    pX=(const ae_int32x4 *)(x);
    pY=(const ae_int32x4 *)(y);
    for (n=0; n<(N>>3); n++) 
    {
        AE_L32X2X2_I (x2,x3,pX,1*sizeof(ae_int32x4));
        AE_L32X2X2_IP(x0,x1,pX,2*sizeof(ae_int32x4));
        AE_L32X2X2_I (y2,y3,pY,1*sizeof(ae_int32x4));
        AE_L32X2X2_IP(y0,y1,pY,2*sizeof(ae_int32x4));
        AE_MULAAF2D32RA_HH_LL(a0,a1,x0,x1,y0,y1);
        AE_MULAAF2D32RA_HH_LL(a2,a3,x2,x3,y2,y3);
    }
    if (N&7)
    {
        AE_L32X2X2_IP(x0,x1,pX,sizeof(ae_int32x4));
        AE_L32X2X2_IP(y0,y1,pY,sizeof(ae_int32x4));
        AE_MULAAF2D32RA_HH_LL(a0,a1,x0,x1,y0,y1);
    }
    a0=AE_ADD64(a0,a1);
    a2=AE_ADD64(a2,a3);
    a0=AE_ADD64(a0,a2);
    return (int64_t)AE_SRAI64(a0, 16);
}
