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
int64_t vec_dot32x32 (const int32_t * restrict x,const int32_t * restrict y,int N)
{
    static const int32_t ALIGN(16) seq[]={0,1,2,3,4,5,6,7};
    xtbool2 b01,b23,b45,b67;
    ae_int32x2 s01,s23,s45,s67;
    const ae_int32x4* restrict pX;
    const ae_int32x4* restrict pY;
    ae_int32x2 x0,x1,x2,x3;
    ae_int32x2 y0,y1,y2,y3;
    ae_int64 a0,a1,a2,a3;
    ae_valignx2 ax,ay;
    int n,N0;
    if(N<=0) return 0;
    if (N<=7)
    {
        a0=AE_ZERO64();
        __Pragma("no_unroll")
        for (n=0; n<N; n++) 
        {
            AE_L32_IP(x0,castxcc(ae_int32,x),sizeof(int32_t));
            AE_L32_IP(y0,castxcc(ae_int32,y),sizeof(int32_t));
            AE_MULAF32R_HH(a0,x0,y0);
        }
        return (int64_t)AE_SRAI64(a0, 16);
    }
    AE_L32X2X2_I(s01,s23,(const ae_int32x4*)seq,0*sizeof(ae_int32x4));
    AE_L32X2X2_I(s45,s67,(const ae_int32x4*)seq,1*sizeof(ae_int32x4));
    N0=((N-1)&7)+1;
    b01=AE_LT32(s01,N0);    // mask unnessesary elements on the first iteration
    b23=AE_LT32(s23,N0);
    b45=AE_LT32(s45,N0);
    b67=AE_LT32(s67,N0);
    pX=(const ae_int32x4 *)x;
    ax=AE_LA128_PP(pX);
    pY=(const ae_int32x4 *)y;
    ay=AE_LA128_PP(pY);
    AE_LA32X2X2_IP(x0,x1,ax,pX);
    AE_LA32X2X2_IP(x2,x3,ax,pX);
    AE_LA32X2X2_IP(y0,y1,ay,pY);
    AE_LA32X2X2_IP(y2,y3,ay,pY);
    AE_MOVF32X2(x0,0,b01);
    AE_MOVF32X2(x1,0,b23);
    AE_MOVF32X2(x2,0,b45);
    AE_MOVF32X2(x3,0,b67);
    AE_MULZAAF2D32RA_HH_LL(a0,a1,x0,x1,y0,y1);
    AE_MULZAAF2D32RA_HH_LL(a2,a3,x2,x3,y2,y3);
    N-=N0;
    pX=(const ae_int32x4 *)(x+N0);
    ax=AE_LA128_PP(pX);
    pY=(const ae_int32x4 *)(y+N0);
    ay=AE_LA128_PP(pY);
    for (n=0; n<(N>>3); n++) 
    {
        AE_LA32X2X2_IP(x0,x1,ax,pX);
        AE_LA32X2X2_IP(x2,x3,ax,pX);
        AE_LA32X2X2_IP(y0,y1,ay,pY);
        AE_LA32X2X2_IP(y2,y3,ay,pY);
        AE_MULAAF2D32RA_HH_LL(a0,a1,x0,x1,y0,y1);
        AE_MULAAF2D32RA_HH_LL(a2,a3,x2,x3,y2,y3);
    }
    a0=AE_ADD64(a0,a1);
    a2=AE_ADD64(a2,a3);
    a0=AE_ADD64(a0,a2);
    return (int64_t)AE_SRAI64(a0, 16);
}
