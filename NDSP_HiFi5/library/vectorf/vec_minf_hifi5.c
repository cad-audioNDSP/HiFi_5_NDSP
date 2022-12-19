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
/*
  NatureDSP Signal Processing Library. Vector Operations
    Vector Min/Max
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/
#include "common.h"
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_vector.h"
#include "common_fpu.h"
#include "inff_tbl.h"

#if (HAVE_VFPU==0 && HAVE_FPU==0)
DISCARD_FUN(float32_t,vec_minf,(const float32_t* restrict x, int N))
#elif (HAVE_VFPU)

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
float32_t vec_minf     (const float32_t* restrict x, int N)
{
    const xtfloatx4* restrict pX;
    xtfloatx2 x0,x1,x2,x3,min0,min1;
    xtfloat minval;
    ae_valignx2 ax;
    int n,N0;
    NASSERT(x);
    if (N <= 0) return 0;
    if (N<=7)
    {
        minval=plusInff.f;
        __Pragma("no_unroll")
        for (n=0; n<N; n++) minval=MINNUM_S(x[n],minval);
        return minval;
    }
    N0=((N-1)&7)+1;
    pX=(const xtfloatx4 *)x;
    ax=AE_LA128_PP(pX);
    AE_LASX2X2_IP(x0,x1,ax,pX);
    AE_LASX2X2_IP(x2,x3,ax,pX);
    pX=(const xtfloatx4 *)(x+N0);
    ax=AE_LA128_PP(pX);
    N-=N0;
    min0=MINNUM_SX2(x0,x1);
    min1=MINNUM_SX2(x2,x3);
    for (n=0; n<(N>>3); n++) 
    {
        AE_LASX2X2_IP(x0,x1,ax,pX);
        AE_LASX2X2_IP(x2,x3,ax,pX);
        min0=MINNUM_SX2(min0,MINNUM_SX2(x0,x1));
        min1=MINNUM_SX2(min1,MINNUM_SX2(x2,x3));
    }
    min0=MINNUM_SX2(min0,min1);
    min1=AE_SEL32_LH_SX2(min0,min0);
    min0=MINNUM_SX2(min0,min1);
    return XT_LOW_S(min0);
}

#else
float32_t vec_minf     (const float32_t* restrict x, int N)
{
    const xtfloat * pX=(const xtfloat * )x;
    static const union ufloat32uint32 plusInff ={0x7f800000}; /* +Inf */
    xtfloat min =plusInff.f;
    int n;
    NASSERT( x );
    for ( n=0; n<N; n++ )
    {
        xtfloat x0;
        xtbool lt;
        XT_LSIP(x0,pX,sizeof(xtfloat));
        lt=XT_OLT_S(x0,min);
        XT_MOVT_S(min,x0,lt);
    }
    XT_MOVF_S(min,XT_CONST_S(0) ,(xtbool)(N>0));
    return min;
}
#endif
