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
int64_t vec_power32x32 (const int32_t * restrict x, int rsh, int N)
{
    static const int32_t ALIGN(16) seq[]={0,1,2,3,4,5,6,7};
    xtbool2 b01,b23,b45,b67;
    ae_int32x2 s01,s23,s45,s67;
    const ae_int32x4* restrict pX;
    const ae_int32x4* restrict pY;
    ae_int32x2 x0,x1,x2,x3;
    ae_int64 a0,a1,a2,a3;
    ae_valignx2 ax,ay;
    int n,N0;
    NASSERT(rsh>=31 && rsh<=62);
    if(N<=0) return 0;
    if (N<=7)
    {
        a0=AE_ZERO64();
        __Pragma("no_unroll")
        for (n=0; n<N; n++) 
        {
            AE_L32_IP(x0,castxcc(ae_int32,x),sizeof(int32_t));
            AE_MULAF32R_HH(a0,x0,x0);
        }
        return (int64_t)AE_SLAA64S(a0, 15-rsh);;
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
    AE_LA32X2X2_IP(x0,x1,ax,pX);
    AE_LA32X2X2_IP(x2,x3,ax,pX);
    AE_MOVF32X2(x0,0,b01);
    AE_MOVF32X2(x1,0,b23);
    AE_MOVF32X2(x2,0,b45);
    AE_MOVF32X2(x3,0,b67);
    AE_MULZAAF2D32RA_HH_LL(a0,a1,x0,x1,x0,x1);
    AE_MULZAAF2D32RA_HH_LL(a2,a3,x2,x3,x2,x3);
    N-=N0;
    pX=(const ae_int32x4 *)(x+N0);
    ax=AE_LA128_PP(pX);
#if 1
    pY=(const ae_int32x4 *)((x+N0)+(N>>1));
    ay=AE_LA128_PP(pY);
    for (n=0; n<(N>>3); n++) 
    {
        AE_LA32X2X2_IP(x0,x1,ax,pX);
        AE_LA32X2X2_IP(x2,x3,ay,pY);
        AE_MULAAF2D32RA_HH_LL(a0,a1,x0,x1,x0,x1);
        AE_MULAAF2D32RA_HH_LL(a2,a3,x2,x3,x2,x3);
    }
#else
    pY=(const ae_int32x4 *)((x+N0)+(N>>1));
    ay=AE_LA128_PP(pY);
    if ((N>>4)>0)
    {
        AE_LA32X2X2_IP(x0,x1,ax,pX);
        AE_LA32X2X2_IP(x2,x3,ay,pY);
        AE_MULAAF2D32RA_HH_LL(a0,a1,x0,x1,x0,x1);
        AE_MULAAF2D32RA_HH_LL(a2,a3,x2,x3,x2,x3);
        AE_LA32X2X2_IP(x0,x1,ax,pX);
        for (n=1; n<(N>>4); n++) 
        {
            AE_LA32X2X2_IP(x2,x3,ay,pY);
            AE_MULAAF2D32RA_HH_LL(a0,a1,x0,x1,x0,x1);
            AE_MULAAF2D32RA_HH_LL(a2,a3,x2,x3,x2,x3);
            AE_LA32X2X2_IP(x0,x1,ax,pX);
            AE_LA32X2X2_IP(x2,x3,ay,pY);
            AE_MULAAF2D32RA_HH_LL(a0,a1,x0,x1,x0,x1);
            AE_MULAAF2D32RA_HH_LL(a2,a3,x2,x3,x2,x3);
            AE_LA32X2X2_IP(x0,x1,ax,pX);
        }
        AE_LA32X2X2_IP(x2,x3,ay,pY);
        AE_MULAAF2D32RA_HH_LL(a0,a1,x0,x1,x0,x1);
        AE_MULAAF2D32RA_HH_LL(a2,a3,x2,x3,x2,x3);
    }
    if (N&8)
    {
        AE_LA32X2X2_IP(x0,x1,ax,pX);
        AE_LA32X2X2_IP(x2,x3,ay,pY);
        AE_MULAAF2D32RA_HH_LL(a0,a1,x0,x1,x0,x1);
        AE_MULAAF2D32RA_HH_LL(a2,a3,x2,x3,x2,x3);
    }
#endif
    a0=AE_ADD64(a0,a1);
    a2=AE_ADD64(a2,a3);
    return (int64_t)AE_SLAA64S(AE_ADD64(a0,a2), 15-rsh);;
}
