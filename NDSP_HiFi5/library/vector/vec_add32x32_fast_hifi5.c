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
void vec_add32x32_fast ( int32_t * restrict z,
              const int32_t * restrict x,
              const int32_t * restrict y,
              int N)
{
    int n;

    ae_int32x2 x0,x1,x2,x3;
    ae_int32x2 y0,y1,y2,y3;
    const ae_int32x4 * restrict pX = (const ae_int32x4 *)x;
    const ae_int32x4 * restrict pY = (const ae_int32x4 *)y;
          ae_int32x4 * restrict pZ = (      ae_int32x4 *)z;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT((N&3)==0);
    if(N<=0) return ;
    pX=(const ae_int32x4 *)x;
    pY=(const ae_int32x4 *)y;
    pZ=(      ae_int32x4 *)z;

    for (n=0; n<(N>>3); n++)
    {
        AE_L32X2X2_IP(x0,x1,pX,sizeof(ae_int32x4));
        AE_L32X2X2_IP(x2,x3,pX,sizeof(ae_int32x4));
        AE_L32X2X2_IP(y0,y1,pY,sizeof(ae_int32x4));
        AE_L32X2X2_IP(y2,y3,pY,sizeof(ae_int32x4));
        y0=AE_ADD32S(x0,y0);
        y1=AE_ADD32S(x1,y1);
        y2=AE_ADD32S(x2,y2);
        y3=AE_ADD32S(x3,y3);
        AE_S32X2X2_IP(y0,y1,pZ,sizeof(ae_int32x4));
        AE_S32X2X2_IP(y2,y3,pZ,sizeof(ae_int32x4));
    }
    if (N&7)
    {
        AE_L32X2X2_IP(x0,x1,pX,sizeof(ae_int32x4));
        AE_L32X2X2_IP(y0,y1,pY,sizeof(ae_int32x4));
        y0=AE_ADD32S(x0,y0);
        y1=AE_ADD32S(x1,y1);
        AE_S32X2X2_IP(y0,y1,pZ,sizeof(ae_int32x4));
    }
}
