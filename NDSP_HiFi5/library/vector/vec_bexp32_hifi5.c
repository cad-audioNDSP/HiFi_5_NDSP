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
  Common Exponent 
  These functions determine the number of redundant sign bits for each value 
  (as if it was loaded in a 32-bit register) and returns the minimum number 
  over the whole vector. This may be useful for a FFT implementation to 
  normalize data.  
  NOTES:
  Faster version of functions make the same task but in a different manner - 
  they compute exponent of maximum absolute value in the array. It allows 
  faster computations but not bitexact results - if minimum value in the 
  array will be -2^n , fast function returns max(0,30-n) while non-fast 
  function returns (31-n).
  Floating point function returns 0-floor(log2(max(abs(x)))). Returned 
  result will be always in range [-129...146]. 
  Special cases
  x       | result
  --------+-------
  0       |    0
  +/-Inf  | -129
  NaN     |    0

  If dimension N<=0 functions return 0

  Precision: 
  32 32-bit inputs 
  16 16-bit inputs 
  f  single precision floating point

  Input:
  x[N]    input data
  N       length of vector
  Returned value:
  minimum exponent

  Restriction:
  Regular versions (vec_bexp16, vec_bexp32, vec_bexpf):
  none
  Faster versions (vec_bexp16_fast, vec_bexp32_fast):
  x   - aligned on 16-byte boundary
  N   - a multiple of 4
-------------------------------------------------------------------------*/
int vec_bexp32 (const int32_t * restrict x, int N)
{
    const ae_int32x4 * restrict pX;
    ae_valignx2 ax;
    ae_int32x2 x0,x1,x2,x3;
    ae_int16x4 nsa0,nsa1;
    int nsa,n,N0;
    if (N<=0) return 0;
    if (N<=7) 
    {
        __Pragma("no_unroll");
        for (n=0,nsa=31; n<N;n++) nsa=XT_MIN(nsa,NSA(x[n]));
        return nsa;
    }
    N0=((N-1)&7)+1;
    pX=(const ae_int32x4 *)x;
    ax=AE_LA128_PP(pX);
    AE_LA32X2X2_IP(x0,x1,ax,pX);
    AE_LA32X2X2_IP(x2,x3,ax,pX);
    pX=(const ae_int32x4 *)(x+N0);
    ax=AE_LA128_PP(pX);
    N-=N0;
    nsa0=AE_NSA32X4(x0,x1);
    nsa1=AE_NSA32X4(x2,x3);
    for (n=0; n<(N>>3); n++) 
    {
        AE_LA32X2X2_IP(x0,x1,ax,pX);
        AE_LA32X2X2_IP(x2,x3,ax,pX);
        nsa0=AE_MIN16(nsa0,AE_NSA32X4(x0,x1));
        nsa1=AE_MIN16(nsa1,AE_NSA32X4(x2,x3));
    }
    nsa=AE_RMIN16X4(AE_MIN16(nsa0,nsa1));
    return nsa;
}
