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

#include "common.h"
#include "NatureDSP_types.h"
/* Library API */
#include "NatureDSP_Signal_vector.h"
// code optimized for HiFi5 core

/*-------------------------------------------------------------------------
  Vector Min/Max
  These routines find maximum/minimum value in a vector.
  Two versions of functions available: regular version (vec_min32x32, 
  vec_max32x32, vec_min16x16, vec_max16x16, vec_minf, vec_maxf) with 
  arbitrary arguments and faster version (vec_min32x32_fast, vec_max32x32_fast, 
  vec_min16x16_fast, vec_max16x16_fast) that apply some restrictions
  NOTE: functions return zero if N is less or equal to zero

  Precision: 
  32x32 32-bit data, 32-bit output
  16x16 16-bit data, 16-bit output
  f     single precision floating point
  
  Input:
  x[N]  input data
  N     length of vector
  Returned value:
  Minimum or maximum value correspondingly

  Restriction:
  For regular routines:
  none
  For faster routines:
  x aligned on 16-byte boundary
  N   - multiple of 4
-------------------------------------------------------------------------*/
int16_t vec_min16x16_fast (const int16_t* restrict x, int N)
#if 0
{
    int         n;
    const ae_int16x4 * restrict px = (const ae_int16x4 *)x;
    ae_int16x4  vxh, h0,h1;
    xtbool4     cprt;

    h0=h1 = AE_MOVDA16(0x7fff);
    NASSERT_ALIGN8(x);
    NASSERT((N&3)==0);
  if (N <= 0) return 0;
    for (n=0;n<(N>>3); n++)
    {
        AE_L16X4_IP(vxh, px, +8);
        cprt = AE_LT16(h0, vxh);
        AE_MOVF16X4(h0, vxh, cprt);

        AE_L16X4_IP(vxh, px, +8);
        cprt = AE_LT16(h1, vxh);
        AE_MOVF16X4(h1, vxh, cprt);
    }
    if(N&4)
    {
        AE_L16X4_IP(vxh, px, +8);
        cprt = AE_LT16(h0, vxh);
        AE_MOVF16X4(h0, vxh, cprt);
    }
    cprt = AE_LT16(h0, h1);
    AE_MOVF16X4(h0, h1, cprt);

    int m01,m23;
    m01=XT_MIN(AE_MOVAD16_0(h0),AE_MOVAD16_1(h0));
    m23=XT_MIN(AE_MOVAD16_2(h0),AE_MOVAD16_3(h0));
    return XT_MIN(m01,m23);
}
#else
{
    const ae_int16x8* restrict pX;
    ae_int16x4 x0,x1,x2,x3,min0,min1;
    int n;
    NASSERT_ALIGN16(x);
    NASSERT((N&3)==0);
    if (N <= 0) return 0;
    pX=(const ae_int16x8 *)x;
    min0=MAX_INT16; min1=MAX_INT16;
    for (n=0; n<(N>>4); n++) 
    {
        AE_L16X4X2_IP(x0,x1,pX,sizeof(ae_int16x8));
        AE_L16X4X2_IP(x2,x3,pX,sizeof(ae_int16x8));
        min0=AE_MIN16(min0,AE_MIN16(x0,x1));
        min1=AE_MIN16(min1,AE_MIN16(x2,x3));
    }
    if (N&15)
    {
        N&=15;
        int off1,off2;
        off1=N>4?4*sizeof(int16_t):0;
        off2=N>8?8*sizeof(int16_t):0;
        x0=AE_L16X4_I((const ae_int16x4*)pX,0);
        x1=AE_L16X4_X((const ae_int16x4*)pX,off1);
        x2=AE_L16X4_X((const ae_int16x4*)pX,off2);
        min0=AE_MIN16(min0,AE_MIN16(x0,x1));
        min1=AE_MIN16(min1,x2);
    }
    min0=AE_MIN16(min0,min1);
    return AE_RMIN16X4(min0);
}
#endif
