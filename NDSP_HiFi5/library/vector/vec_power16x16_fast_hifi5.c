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
#include "NatureDSP_types.h"
#include "common.h"
// code optimized for HiFi5 core

/*===========================================================================
  Vector matematics:
  vec_power            Power of a Vector
===========================================================================*/

/*-------------------------------------------------------------------------
  Power of a Vector
  These routines compute power of vector with scaling output result by rsh 
  bits. Fixed point rountines make accumulation in the 64-bit wide 
  accumulator and output may scaled down with saturation by rsh bits. 
  So, if representation of x input is Qx, result will be represented in 
  Q(2x-rsh) format.
  Two versions of routines are available: regular versions (vec_power32x32, 
  vec_power16x16, vec_powerf) work with arbitrary arguments, faster versions 
  (vec_power32x32_fast, vec_power16x16_fast) apply some restrictions.

  Precision: 
  32x32 32x32-bit data, 64-bit output
  16x16 16x16-bit data, 64-bit output
  f     single precision floating point

  Input:
  x[N]  input data, Q31, Q15 or floating point
  rsh   right shift of result
  N     length of vector
  Returns: 
  Sum of squares of a vector, Q(2x-rsh)

  Restrictions:
  for vec_power32x32(): rsh in range 31...62
  for vec_power16x16(): rsh in range 0...31
  For regular versions (vec_power32x32, vec_power16x16, vec_powerf):
  none
  For faster versions (vec_power32x32_fast, vec_power16x16_fast ):
  x - aligned on 16-byte boundary
  N - multiple of 4
-------------------------------------------------------------------------*/
int64_t vec_power16x16_fast (const int16_t * restrict x, int rsh, int N)
{
    int n;
    static const int16_t ALIGN(16) seq[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

    ae_int64    vai,z0,z1,z2,z3;
    ae_int16x4  x0,x1,x2,x3;
    const ae_int16x8 * restrict px = (const ae_int16x8*)x;

    if(N<=0) return 0;
    NASSERT_ALIGN16(x);
    NASSERT((N&3)==0);
    NASSERT(rsh>=0 && rsh<=31);
    z0=z1=z2=z3=AE_ZERO64();
    for (n=0; n<(N>>4); n++)
    {
        AE_L16X4X2_I (x2,x3, px,  1*sizeof(ae_int16x8));
        AE_L16X4X2_IP(x0,x1, px,  2*sizeof(ae_int16x8));
        AE_MULAAAA2Q16(z0,z1,x0,x1,x0,x1);
        AE_MULAAAA2Q16(z2,z3,x2,x3,x2,x3);
    }
    N&=15;
    if (N==0)
    {
        return (int64_t)AE_SRAA64(AE_ADD64(AE_ADD64(z0,z1),AE_ADD64(z2,z3)), rsh);
    }
    {
        int off1,off2;
        xtbool4 b0123,b4567,b89ab,bcdef;
        ae_int16x4 s0123,s4567,s89ab,scdef;
        AE_L16X4X2_I(s0123,s4567,(const ae_int16x8*)seq,0*sizeof(ae_int16x8));
        AE_L16X4X2_I(s89ab,scdef,(const ae_int16x8*)seq,1*sizeof(ae_int16x8));
        b0123=AE_LT16(s0123,N);    // mask unnessesary elements on the first iteration
        b4567=AE_LT16(s4567,N);
        b89ab=AE_LT16(s89ab,N);
        bcdef=AE_LT16(scdef,N);
        off1=N>4?4*sizeof(int16_t):0;
        off2=N>8?8*sizeof(int16_t):0;
        x0=AE_L16X4_I((const ae_int16x4*)px,0);
        x1=AE_L16X4_X((const ae_int16x4*)px,off1);
        x2=AE_L16X4_X((const ae_int16x4*)px,off2);
        AE_MOVF16X4(x0,0,b0123);
        AE_MOVF16X4(x1,0,b4567);
        AE_MOVF16X4(x2,0,b89ab);
        AE_MULAAAA2Q16(z0,z1,x0,x1,x0,x1);
        AE_MULAAAAQ16(z2,x2,x2);
    }
    vai = AE_ADD64(AE_ADD64(z0,z1),AE_ADD64(z2,z3));
    vai = AE_SRAA64(vai, rsh);
    return_int64(vai);
}
