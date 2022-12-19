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
int16_t vec_max16x16_fast (const int16_t* restrict x, int N)
{
    const ae_int16x8* restrict pX;
    ae_int16x4 x0,x1,x2,x3,max0,max1;
    int n;
    NASSERT_ALIGN16(x);
    NASSERT((N&3)==0);
    if (N <= 0) return 0;
    pX=(const ae_int16x8 *)x;
    max0=MIN_INT16; max1=MIN_INT16;
    for (n=0; n<(N>>4); n++) 
    {
        AE_L16X4X2_IP(x0,x1,pX,sizeof(ae_int16x8));
        AE_L16X4X2_IP(x2,x3,pX,sizeof(ae_int16x8));
        max0=AE_MAX16(max0,AE_MAX16(x0,x1));
        max1=AE_MAX16(max1,AE_MAX16(x2,x3));
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
        max0=AE_MAX16(max0,AE_MAX16(x0,x1));
        max1=AE_MAX16(max1,x2);
    }
    max0=AE_MAX16(max0,max1);
    return AE_RMAX16X4(max0);
}
