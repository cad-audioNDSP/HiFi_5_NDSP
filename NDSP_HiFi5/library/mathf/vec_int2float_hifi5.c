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

/* DSP Library API */
/*    Code optimized for HiFi5 core */
#include "NatureDSP_Signal_math.h"
/* Common helper macros. */
#include "common_fpu.h"

#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(void, vec_int2float,( float32_t * restrict y, const int32_t * restrict x, int t, int N ))
#elif HAVE_VFPU
/*
  NatureDSP Signal Processing Library. Vector Mathematics
   Integer to float conversion
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/
/*-------------------------------------------------------------------------
  Integer to float conversion
  Routines convert integer to float and scale result up by 2^t.

  Precision: 
  f     32-bit input, floating point output

  Input:
  x[N]  input data, integer
  t     scale factor
  N     length of vector
  Output:
  y[N]  conversion results, floating point

  Restriction:
  t should be in range -126...126
-------------------------------------------------------------------------*/
#define sz_f32    (int)sizeof(float32_t)
#define sz_i32    (int)sizeof(int32_t)
void   vec_int2float ( float32_t  * restrict y, const int32_t  * restrict x, int t, int N)
{
    const ae_int32x4 * restrict pX=(const ae_int32x4 *)x;
          xtfloatx4  * restrict pY=(      xtfloatx4  *)y;
    int n,N0;
    ae_int32x2 x0,x1,x2,x3;
    xtfloatx2 y0,y1,y2,y3,s;
    ae_valignx2 aX,aY;
    if (N<=0) return;
    if (N<=8)
    {
        xtfloat   s=FLOATEXP_S(AE_MOVINT8_FROMINT8X8(AE_MOVDA8(t)));
        __Pragma("loop_count max=8")
        for (n=0; n<N; n++)
        {
            y[n]=XT_MUL_S(XT_FLOAT_S(x[n],0),s);
        }
        return;
    }
    s=FLOATEXP_SX2(AE_MOVDA8(t));
    aX=AE_LA128_PP(pX);
    aY=AE_ZALIGN128();
    AE_LA32X2X2_IP(x0,x1,aX,pX);
    AE_LA32X2X2_IP(x2,x3,aX,pX);
    y0=XT_FLOAT_SX2(x0,0);
    y1=XT_FLOAT_SX2(x1,0);
    y2=XT_FLOAT_SX2(x2,0);
    y3=XT_FLOAT_SX2(x3,0);
    MULQ_S(y0,y1,y0,y1,s);
    MULQ_S(y2,y3,y2,y3,s);
    AE_SASX2X2_IP(y0,y1,aY,pY);
    AE_SASX2X2_IP(y2,y3,aY,pY);
    AE_SA128POS_FP(aY,pY);
    N0=((N-1)&7)+1;
    N-=N0;
    if (N<=0) return;
    pX=(const ae_int32x4 *)(x+N0);
    pY=(      xtfloatx4  *)(y+N0);
    aX=AE_LA128_PP(pX);
    for (n=0; n<(N>>3); n++)
    {
        AE_LA32X2X2_IP(x0,x1,aX,pX);
        AE_LA32X2X2_IP(x2,x3,aX,pX);
        y0=XT_FLOAT_SX2(x0,0);
        y1=XT_FLOAT_SX2(x1,0);
        y2=XT_FLOAT_SX2(x2,0);
        y3=XT_FLOAT_SX2(x3,0);
        MULQ_S(y0,y1,y0,y1,s);
        MULQ_S(y2,y3,y2,y3,s);
        AE_SASX2X2_IP(y0,y1,aY,pY);
        AE_SASX2X2_IP(y2,y3,aY,pY);
    }
    AE_SA128POS_FP(aY,pY);
} /* vec_int2float() */
#elif HAVE_FPU
#define sz_f32    (int)sizeof(float32_t) 

void vec_int2float( float32_t * restrict y, const int32_t * restrict x, int t, int N )
{
	int n;
	xtfloat f0, y0;
	int32_t x0;

	NASSERT(x);
	NASSERT(y);
	NASSERT(t >= -126 && t <= 126);
	if (N <= 0) return;
	y0 = XT_WFR((uint32_t)(t + 127) << 23);
	for (n=0; n<N; n++)
	{
		x0 = XT_L32I((const int *)x, 0);
		x++;
		f0 = XT_FLOAT_S(x0, 0);
		f0 = XT_MUL_S(f0, y0);
		XT_SSIP(f0, castxcc(xtfloat, y), sz_f32);
	}
} /* vec_int2float() */
#endif
