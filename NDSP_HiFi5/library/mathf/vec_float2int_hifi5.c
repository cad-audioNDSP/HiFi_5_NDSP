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
#include "common_fpu.h"

#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(void,vec_float2int,( int32_t * restrict y, const float32_t * restrict x, int t, int N ))
#elif HAVE_VFPU
#define sz_f32    (int)sizeof(float32_t) 
#define sz_i32    (int)sizeof(int32_t)
/*
  NatureDSP Signal Processing Library. Vector Mathematics
  Float to integer conversion
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/
/*-------------------------------------------------------------------------
  Float to integer conversion
  routines scale floating point input down by 2^t and convert it to integer 
  with saturation

  Precision: 
  f     single precision floating point

  Input:
  x[N]  input data, floating point
  t     scale factor
  N     length of vector
  Output:
  y[N]  conversion results, integer

  Restriction:
  t should be in range -126...126
-------------------------------------------------------------------------*/
void   vec_float2int (  int32_t * restrict y, const float32_t * restrict x, int t, int N)
{
    const xtfloatx4  * restrict pX=(const xtfloatx4 *)x;
          ae_int32x4 * restrict pY=(      ae_int32x4*)y;
    int n,N0;
    ae_int32x2  y0,y1,y2,y3;
    xtfloatx2 x0,x1,x2,x3,s;
    ae_valignx2 aX,aY;
    if (N<=0) return;
    if (N<=8)
    {
        xtfloat   s=FLOATEXP_S(AE_MOVINT8_FROMINT8X8(AE_MOVDA8(-t)));
        __Pragma("loop_count max=8")
        for (n=0; n<N; n++)
        {
            y[n]=XT_TRUNC_S(XT_MUL_S(x[n],s),0);
        }
        return;
    }
    s=FLOATEXP_SX2(AE_MOVDA8(-t));
    aX=AE_LA128_PP(pX);
    aY=AE_ZALIGN128();
    AE_LASX2X2_IP(x0,x1,aX,pX);
    AE_LASX2X2_IP(x2,x3,aX,pX);
    MULQ_S(x0,x1,x0,x1,s);
    MULQ_S(x2,x3,x2,x3,s);
    y0=XT_TRUNC_SX2(x0,0);
    y1=XT_TRUNC_SX2(x1,0);
    y2=XT_TRUNC_SX2(x2,0);
    y3=XT_TRUNC_SX2(x3,0);
    AE_SA32X2X2_IP(y0,y1,aY,pY);
    AE_SA32X2X2_IP(y2,y3,aY,pY);
    AE_SA128POS_FP(aY,pY);
    N0=((N-1)&7)+1;
    N-=N0;
    if (N<=0) return;
    pX=(const xtfloatx4 * )(x+N0);
    pY=(      ae_int32x4* )(y+N0);
    aX=AE_LA128_PP(pX);
    for (n=0; n<(N>>3); n++)
    {
        AE_LASX2X2_IP(x0,x1,aX,pX);
        AE_LASX2X2_IP(x2,x3,aX,pX);
        MULQ_S(x0,x1,x0,x1,s);
        MULQ_S(x2,x3,x2,x3,s);
        y0=XT_TRUNC_SX2(x0,0);
        y1=XT_TRUNC_SX2(x1,0);
        y2=XT_TRUNC_SX2(x2,0);
        y3=XT_TRUNC_SX2(x3,0);
        AE_SA32X2X2_IP(y0,y1,aY,pY);
        AE_SA32X2X2_IP(y2,y3,aY,pY);
    }
    AE_SA128POS_FP(aY,pY);
} /* vec_float2int() */
#elif HAVE_FPU
#define sz_f32    (int)sizeof(float32_t) 

void vec_float2int( int32_t * restrict y, const float32_t * restrict x, int t, int N )
{
    int n;
	xtfloat f0, y0;
	int32_t x0;

	NASSERT(x);
	NASSERT(y);
	NASSERT(t >= -126 && t <= 126);
	if (N <= 0) return;
	y0 = XT_WFR((uint32_t)(-t + 127) << 23);
	for (n=0; n<N; n++)
	{
        XT_LSIP(f0, castxcc(xtfloat, x), sz_f32);
        f0 = XT_MUL_S(f0, y0);
        x0 = XT_TRUNC_S(f0, 0);
		XT_S32I(x0, (int *)y, 0);
		y++;
	}
} /* vec_float2int() */
#endif
