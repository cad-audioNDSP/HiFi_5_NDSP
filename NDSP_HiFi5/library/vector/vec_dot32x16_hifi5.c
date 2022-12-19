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
int64_t vec_dot32x16 (const int32_t * restrict x, const int16_t * restrict y, int N)
{
    static const int16_t ALIGN(16) seq[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    xtbool4 b0123,b4567,b89ab,bcdef;
    ae_int16x4 s0123,s4567,s89ab,scdef;
    ae_int64    z0,z1,z2,z3;
    ae_int32x2  x0,x1,x2,x3,x4,x5,x6,x7;
    ae_int16x4  y0,y1,y2,y3;
    const ae_int32x4 * restrict pX0;
    const ae_int16x8 * restrict pY0;
    const ae_int32x4 * restrict pX1;
    const ae_int16x8 * restrict pY1;
    ae_valignx2 ax0,ay0;
    ae_valignx2 ax1,ay1;

    int n,N0;
    if(N<=0) return 0;
    if (N<=15)
    {
        z0=AE_ZERO64();
        __Pragma("no_unroll")
        for (n=0; n<N; n++) 
        {
            AE_L32_IP(x0,castxcc(ae_int32,x),sizeof(int32_t));
            AE_L16_IP(y0,castxcc(ae_int16,y),sizeof(int16_t));
            AE_MULAAAAFQ32X16(z0,x0,x0,y0);
        }
        return (int64_t)AE_SRAI64(z0,18);
    }
    AE_L16X4X2_I(s0123,s4567,(const ae_int16x8*)seq,0*sizeof(ae_int16x8));
    AE_L16X4X2_I(s89ab,scdef,(const ae_int16x8*)seq,1*sizeof(ae_int16x8));
    N0=((N-1)&15)+1;
    b0123=AE_LT16(s0123,N0);    // mask unnessesary elements on the first iteration
    b4567=AE_LT16(s4567,N0);
    b89ab=AE_LT16(s89ab,N0);
    bcdef=AE_LT16(scdef,N0);
    pX0=(const ae_int32x4 *)x;
    pY0=(const ae_int16x8 *)y;
    ax0=AE_LA128_PP(pX0);
    ay0=AE_LA128_PP(pY0);
    AE_LA32X2X2_IP(x0,x1,ax0,pX0);
    AE_LA32X2X2_IP(x2,x3,ax0,pX0);
    AE_LA32X2X2_IP(x4,x5,ax0,pX0);
    AE_LA32X2X2_IP(x6,x7,ax0,pX0);
    AE_LA16X4X2_IP(y0,y1,ay0,pY0);
    AE_LA16X4X2_IP(y2,y3,ay0,pY0);
    AE_MOVF16X4(y0,0,b0123);
    AE_MOVF16X4(y1,0,b4567);
    AE_MOVF16X4(y2,0,b89ab);
    AE_MOVF16X4(y3,0,bcdef);
    z0=AE_MULZAAAAFQ32X16(x0,x1,y0);
    z1=AE_MULZAAAAFQ32X16(x2,x3,y1);
    z2=AE_MULZAAAAFQ32X16(x4,x5,y2);
    z3=AE_MULZAAAAFQ32X16(x6,x7,y3);
    x+=N0;
    y+=N0;
    N-=N0;
    pX0=(const ae_int32x4 *)(x);
    pY0=(const ae_int16x8 *)(y);
    pX1=(const ae_int32x4 *)(x+(N>>1));
    pY1=(const ae_int16x8 *)(y+(N>>1));
    ax0=AE_LA128_PP(pX0);
    ay0=AE_LA128_PP(pY0);
    ax1=AE_LA128_PP(pX1);
    ay1=AE_LA128_PP(pY1);
    for (n=0; n<(N>>4); n++) 
    {
        AE_LA32X2X2_IP(x0,x1,ax0,pX0);
        AE_LA32X2X2_IP(x2,x3,ax0,pX0);
        AE_LA32X2X2_IP(x4,x5,ax1,pX1);
        AE_LA32X2X2_IP(x6,x7,ax1,pX1);
        AE_LA16X4X2_IP(y0,y1,ay0,pY0);
        AE_LA16X4X2_IP(y2,y3,ay1,pY1);
        AE_MULAAAAFQ32X16(z0,x0,x1,y0);
        AE_MULAAAAFQ32X16(z1,x2,x3,y1);
        AE_MULAAAAFQ32X16(z2,x4,x5,y2);
        AE_MULAAAAFQ32X16(z3,x6,x7,y3);
    }
    z0=AE_ADD64(z0,z2);
    z1=AE_ADD64(z1,z3);
    return (int64_t)AE_SRAI64(AE_ADD64(z0,z1), 16);
}
