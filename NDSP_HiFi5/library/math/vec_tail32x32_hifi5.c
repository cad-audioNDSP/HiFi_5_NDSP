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
#include "vec_tail32x32.h"

/*
standard tail processing for 32x32 math functions (N<8)
*/

/*-------------------------------------------------------------------------
  Input:
  x[N]  input data
  N     length of vectors (<7)
  fxn   function to be called (with N==8)
  Output:
  y[N]  output data
-------------------------------------------------------------------------*/
void vec_tail32x32 (int32_t * restrict y, const int32_t * restrict x, void (*fxn)(int32_t*,const int32_t*,int), int N)
{
   int32_t ALIGN(16) buff[8];
   int n;
   AE_S32X2X2_I(0,0,(ae_int32x4*)buff,0);
   AE_S32X2X2_I(0,0,(ae_int32x4*)buff,sizeof(ae_int32x4));
   __Pragma("loop_count min=1,max=7")
   __Pragma("no_unroll")
   __Pragma("no_simd")
   for (n=0; n<N; n++) buff[n]=x[n];
   fxn(buff,buff,8);
   __Pragma("loop_count min=1,max=7")
   __Pragma("no_unroll")
   __Pragma("no_simd")
   for (n=0; n<N; n++) y[n]=buff[n];
} 
