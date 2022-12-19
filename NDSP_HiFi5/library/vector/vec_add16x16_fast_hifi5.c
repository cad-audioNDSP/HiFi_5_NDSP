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
  vec_add              Vector Sum
===========================================================================*/

/*-------------------------------------------------------------------------
  Vector Sum
  This routine makes pair wise saturated summation of vectors.
  Two versions of routines are available: regular versions (vec_add32x32, 
  vec_add16x16, vec_addf) work with arbitrary arguments, faster versions 
  (vec_add32x32_fast, vec_add16x16_fast) apply some restrictions.

  Precision: 
  32x32 32-bit inputs, 32-bit output
  16x16 16-bit inputs, 16-bit output
  f     single precision floating point

  Input:
  x[N]   input data
  y[N]   input data
  N      length of vectors
  Output:
  z[N]   output data

  Restriction:
  Regular versions (vec_add32x32, vec_add16x16, vec_addf):
  x,y,z - should not be overlapped
  Faster versions (vec_add32x32_fast, vec_add16x16_fast):
  z,x,y - aligned on 16-byte boundary
  N     - multiple of 4
-------------------------------------------------------------------------*/
void vec_add16x16_fast ( int16_t * restrict z,
              const int16_t * restrict x,
              const int16_t * restrict y,
              int N)
#if 0
{
  int n;
  ae_int16x4  x0, y0, z0;

  const ae_int16x4 * restrict px = (const ae_int16x4 *)x;
  const ae_int16x4 * restrict py = (const ae_int16x4 *)y;
        ae_int16x4 * restrict pz = (      ae_int16x4 *)z;
 
  NASSERT_ALIGN8(x);
  NASSERT_ALIGN8(y);
  NASSERT((N&3)==0);
  if(N<=0) return ;
  __Pragma("loop_count min=1")
  for (n=0; n<(N>>2); n++)
  {
    AE_L16X4_IP(x0, px,8);
    AE_L16X4_IP(y0, py,8);
    z0 = AE_ADD16S(x0, y0);
    AE_S16X4_IP(z0, pz,8);
  }
}
#else
{
    int n;

    ae_int16x4 x0,x1,x2,x3;
    ae_int16x4 y0,y1,y2,y3;
    const ae_int16x8 * restrict pX = (const ae_int16x8 *)x;
    const ae_int16x8 * restrict pY = (const ae_int16x8 *)y;
          ae_int16x8 * restrict pZ = (      ae_int16x8 *)z;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT((N&3)==0);
    if(N<=0) return ;
    pX=(const ae_int16x8 *)x;
    pY=(const ae_int16x8 *)y;
    pZ=(      ae_int16x8 *)z;

    for (n=0; n<(N>>4); n++)
    {
        AE_L16X4X2_IP(x0,x1,pX,sizeof(ae_int16x8));
        AE_L16X4X2_IP(x2,x3,pX,sizeof(ae_int16x8));
        AE_L16X4X2_IP(y0,y1,pY,sizeof(ae_int16x8));
        AE_L16X4X2_IP(y2,y3,pY,sizeof(ae_int16x8));
        y0=AE_ADD16S(x0,y0);
        y1=AE_ADD16S(x1,y1);
        y2=AE_ADD16S(x2,y2);
        y3=AE_ADD16S(x3,y3);
        AE_S16X4X2_IP(y0,y1,pZ,sizeof(ae_int16x8));
        AE_S16X4X2_IP(y2,y3,pZ,sizeof(ae_int16x8));
    }
    N&=15;
    if (N)
    {
        __Pragma("no_unroll")
        for (n=0; n<(N>>2); n++)
        {
            AE_L16X4_IP(x0,castxcc(ae_int16x4,pX),sizeof(ae_int16x4));
            AE_L16X4_IP(y0,castxcc(ae_int16x4,pY),sizeof(ae_int16x4));
            y0=AE_ADD16S(x0,y0);
            AE_S16X4_IP(y0,castxcc(ae_int16x4,pZ),sizeof(ae_int16x4));
        }
    }
}
#endif
