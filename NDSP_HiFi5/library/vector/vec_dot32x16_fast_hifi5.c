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

#define SMALLER_CODESIZE 1

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
int64_t vec_dot32x16_fast (const int32_t * restrict x, const int16_t * restrict y, int N)
{
    int n;
    ae_int32x2 x0,x1,x2,x3;
    ae_int16x4 y0,y1,y2,y3;
    ae_f64 a0,a1,a2,a3;
    const ae_int32x4 *    restrict  px = (const ae_int32x4 *)   x;
    const ae_int16x8 *  restrict  py = (const ae_int16x8 *) y;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT((N&3)==0);
    if(N<=0) return 0;
    a0=a1=a2=a3= AE_ZERO64();
    for (n=0; n<(N>>4); n++)
    {
        AE_L32X2X2_I (x2,x3, px,1*sizeof(ae_int32x4));
        AE_L32X2X2_IP(x0,x1, px,2*sizeof(ae_int32x4));
        AE_L16X4X2_I (y2,y3, py,1*sizeof(ae_int16x8));
        AE_L16X4X2_IP(y0,y1, py,2*sizeof(ae_int16x8));
        AE_MULAAAAFQ32X16(a0,x0,x1,y0);
        AE_MULAAAAFQ32X16(a1,x2,x3,y1);
        AE_L32X2X2_I (x2,x3, px,1*sizeof(ae_int32x4));
        AE_L32X2X2_IP(x0,x1, px,2*sizeof(ae_int32x4));
        AE_MULAAAAFQ32X16(a2,x0,x1,y2);
        AE_MULAAAAFQ32X16(a3,x2,x3,y3);
    }
    a0=AE_ADD64(a0,a2);
    a1=AE_ADD64(a1,a3);
    N&=15;
#if SMALLER_CODESIZE
    __Pragma("concurrent")
    __Pragma("no_unroll")
    __Pragma("loop_count max=3")
    for (n=0; n<(N>>2); n++)
    {
        AE_L32X2X2_IP(x0,x1, px,sizeof(ae_int32x4));
        AE_L16X4_IP(y0, castxcc(ae_int16x4,py),sizeof(ae_int16x4));
        AE_MULAAAAFQ32X16(a0,x0,x1,y0);
    }
    return (int64_t)AE_SRAI64(AE_ADD64(a0,a1), 16);
#else
    if (N==0) return (int64_t)AE_SRAI64(AE_ADD64(a0,a1), 16);
    {
        ae_int32x2 x4,x5;
        xtbool4 b1,b2;
        int off1, off2;
        off1=N>4 ?  8:0;
        off2=N>8 ? 16:0;
        b1=AE_EQ16(0,off1);
        b2=AE_EQ16(0,off2);
        y0=AE_L16X4_I ( (const ae_int16x4*)py,    0);
        y1=AE_L16X4_X ( (const ae_int16x4*)py, off1);
        y2=AE_L16X4_X ( (const ae_int16x4*)py, off2);
        AE_MOVT16X4(y1,0,b1);
        AE_MOVT16X4(y2,0,b2);
        AE_L32X2X2_I (x0,x1, px,    0<<1);
        AE_L32X2X2_X (x2,x3, px, off1<<1);
        AE_L32X2X2_X (x4,x5, px, off2<<1);
        a0=AE_ADD64(a0,a1);
        AE_MULAAAAFQ32X16(a0,x0,x1,y0);
        AE_MULAAAAFQ32X16(a0,x2,x3,y1);
        AE_MULAAAAFQ32X16(a0,x4,x5,y2);
        return (int64_t)AE_SRAI64(a0, 16);
    }
#endif
}
