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
#ifndef MTX_MPY8X16_HELPER__
#define MTX_MPY8X16_HELPER__
/* helper functions for matrix 8x16 multiplies  */
#include "NatureDSP_types.h"
#include "common.h"

/* Define conditions for using 'Neural Network' instruction extension */
#if !defined(AE_MULA8Q8X16) || !defined(AE_MULA4O8X16) || !defined(AE_SAV16X4X2_XP) || !defined(AE_LAV16X4X2_XP)
#define USE_NN_EXTENSION_8X16 0
#else
#define USE_NN_EXTENSION_8X16 1
#endif

/* code optimized for HiFi5 core */
void mtx_mpyt8x16_copyz(  
                     int16_t* restrict z,
               const int16_t* restrict v,
               int M, int P );

/*-------------------------------------------------------------------------
  Matrix Multiply

  multiplies x by y' where y is written with paddings

  input:
  x[M][N]
  y[P][N0]

  output:
  z[M][N0] -- zero padded output

  N0 - closest bigger multiple of 8 of N

  z,y - aligned
-------------------------------------------------------------------------*/
void mtx_mpyt8x16_partiallyaligned (  
                     int16_t* restrict z,
               const int8_t * restrict x,
               const int16_t* restrict y,
               int M, int N, int P, int lsh );

/*-------------------------------------------------------------------------
copy small x to the r with padding

input:
x[M][N]
Output:
r[M][8]

r - aligned
N<8
-------------------------------------------------------------------------*/
void mtx_mpy8x16_copysmallx(int8_t * restrict r,const int8_t * restrict x,int M,int N);

#endif
