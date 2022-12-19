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
  NatureDSP Signal Processing Library. FIR part
    complex data circular cross-correlation, complex floating point
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "common_fpu.h"

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void,cxfir_xcorrf,(   complex_float * restrict r,
               const complex_float * restrict x,
               const complex_float * restrict y,
               int N, int M ))
#elif (HAVE_VFPU)

/*-------------------------------------------------------------------------
  Circular Correlation
  Estimates the circular cross-correlation between vectors x (of length N) 
  and y (of length M)  resulting in vector r of length N. It is a similar 
  to correlation but x is read in opposite direction.

  Precision: 
  16x16     16x16-bit data, 16-bit outputs
  32x16     32x16-bit data, 32-bit outputs
  32x32     32x32-bit data, 32-bit outputs
  32x32ep   the same as above but using 72-bit accumulator for intermediate 
            computations
  f         floating point 


  Input:
  x[N]      input data Q15, Q31 or floating point
  y[M]      input data Q15, Q31 or floating point
  N         length of x
  M         length of y
  Output:
  r[N]      output data, Q15, Q31 or floating point

  Restriction:
  x,y,r     should not overlap
  x,y,r     aligned on an 16-bytes boundary
  N,M       multiples of 4 and >0
-------------------------------------------------------------------------*/
void cxfir_xcorrf(   complex_float * restrict r,
               const complex_float * restrict x,
               const complex_float * restrict y,
               int N, int M )
{
    //
    // Circular cross-correlation algorithm:
    //
    //   r[n] = sum( x[mod(n+m,N)]*y[m] )
    //        m=0..M-1
    //
    //   where n = 0..N-1
    //
    const xtfloatx4 * restrict pX;
    const xtfloatx4 * restrict S0;
    const xtfloatx4 * restrict S1;
    const xtfloatx2 * restrict pY;
          xtfloatx4 * restrict pR;
    xtfloatx2 r00, r10, r20, r30;
    xtfloatx2 r01, r11, r21, r31;
    xtfloatx2 r02, r12, r22, r32;
    xtfloatx2 r03, r13, r23, r33;
    xtfloatx2 x0, x1, x2, x3, x4;
    xtfloatx2 y0, y1;
    int n, m;

    NASSERT(r);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(r, 16);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT(M > 0 && M % 4 == 0);
    NASSERT(N > 0 && N % 4 == 0);

    pX = (const xtfloatx4 *)x;
    pR = (      xtfloatx4 *)r;
    /* set circular buffer boundaries */
    WUR_AE_CBEGIN0((uintptr_t)(x + 0));
    WUR_AE_CEND0  ((uintptr_t)(x + N));

    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 2); n++)
    {
        pY = (const xtfloatx2 *)y;
        /* preload data from x */
        AE_LSX2X2_XC(x0, x1, pX, 2 * sizeof(complex_float));
        S0 = pX; S1 = pX;
        AE_ADDCIRC_XC(castxcc(ae_int64, pX), 2 * sizeof(complex_float));

        CONST_SX2X2(r00, r01, 0); CONST_SX2X2(r02, r03, 0);
        CONST_SX2X2(r10, r11, 0); CONST_SX2X2(r12, r13, 0);
        CONST_SX2X2(r20, r21, 0); CONST_SX2X2(r22, r23, 0);
        CONST_SX2X2(r30, r31, 0); CONST_SX2X2(r32, r33, 0);

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 1); m++)
        {
            /* load data from x */
            AE_LSX2X2_XC(x2, x3, S1, 2 * sizeof(complex_float));
            x4 = AE_LSX2I((xtfloatx2*)S1, 0);
            /* load data from y */
            XT_LSX2IP(y0, pY, sizeof(complex_float));
            XT_LSX2IP(y1, pY, sizeof(complex_float));
            /* compute correlation of 4 values */
            MADDMUX_SX2X2(r00, r10, y0, y0, x0, x1, 4);
            MADDMUX_SX2X2(r01, r11, y0, y0, x0, x1, 5);
            MADDMUX_SX2X2(r20, r30, y0, y0, x2, x3, 4);
            MADDMUX_SX2X2(r21, r31, y0, y0, x2, x3, 5);
            MADDMUX_SX2X2(r02, r12, y1, y1, x1, x2, 4);
            MADDMUX_SX2X2(r03, r13, y1, y1, x1, x2, 5);
            MADDMUX_SX2X2(r22, r32, y1, y1, x3, x4, 4);
            MADDMUX_SX2X2(r23, r33, y1, y1, x3, x4, 5);
            /* reload data for the next iteration */
            AE_LSX2X2_XC(x0, x1, S0, 2 * sizeof(complex_float));
        }
        /* reduction add */
        ADD_SX2X2(r00, r10, r00, r10, r01, r11);
        ADD_SX2X2(r20, r30, r20, r30, r21, r31);
        ADD_SX2X2(r02, r12, r02, r12, r03, r13);
        ADD_SX2X2(r22, r32, r22, r32, r23, r33);
        ADD_SX2X2(r00, r10, r00, r10, r02, r12);
        ADD_SX2X2(r20, r30, r20, r30, r22, r32);
        /* save computed samples */
        AE_SSX2X2_IP(r00, r10, pR, 2 * sizeof(complex_float));
        AE_SSX2X2_IP(r20, r30, pR, 2 * sizeof(complex_float));
    }
} // cxfir_xcorrf()
#else
// for scalar FPU
void cxfir_xcorrf(   complex_float * restrict r,
               const complex_float * restrict x,
               const complex_float * restrict y,
               int N, int M )
{
  /*
  * Circular cross-correlation algorithm:
  *
  *   r[n] = sum( x[mod(n+m,N)]*y[m] )
  *        m=0..M-1
  *
  *   where n = 0..N-1
  */
  xtfloat A0_re, A1_re, A0_im, A1_im;
  xtfloat X0_re, X1_re, X0_im, X1_im,Y0_re, Y0_im;
  const xtfloat * restrict pX;
  const xtfloat* restrict pY;
  xtfloat * restrict pR;
  int n, m;

  NASSERT(r);
  NASSERT(x);
  NASSERT(y);
  NASSERT_ALIGN(r,8);
  NASSERT_ALIGN(x,8);
  NASSERT_ALIGN(y,8);
  NASSERT(M>0 && M%4==0);
  NASSERT(N>0 && N%4==0);
  if (N <= 0 || M <= 0) return;
  pR=(      xtfloat*)r;
  /* set circular buffer boundaries */
  WUR_AE_CBEGIN0( (uintptr_t)( x + 0 ) );
  WUR_AE_CEND0  ( (uintptr_t)( x + N ) );
  for ( n=0; n<N; n+=2,x+=2 )
  {
    A0_re=A1_re=A0_im=A1_im=XT_CONST_S(0);
    pX = (const xtfloat*)(x);
    pY = (const xtfloat*)y;
    XT_LSXC(X0_re, pX, 4);
    XT_LSXC(X0_im, pX, 4);
    XT_LSXC(X1_re, pX, 4);
    XT_LSXC(X1_im, pX, 4);
    __Pragma("loop_count min=1");
    for ( m=0; m<(M); m++ )
    {
      XT_LSIP(Y0_re, pY, 4);
      XT_LSIP(Y0_im, pY, 4);
      XT_MADD_S(A0_re, X0_re, Y0_re);
      XT_MADD_S(A0_re, X0_im, Y0_im);
      XT_MADD_S(A0_im, X0_re, Y0_im);
      XT_MSUB_S(A0_im, X0_im, Y0_re);
      XT_MADD_S(A1_re, X1_re, Y0_re);
      XT_MADD_S(A1_re, X1_im, Y0_im);
      XT_MADD_S(A1_im, X1_re, Y0_im);
      XT_MSUB_S(A1_im, X1_im, Y0_re);
      X0_re = X1_re;
      X0_im = X1_im;
      XT_LSXC(X1_re, pX, 4);
      XT_LSXC(X1_im, pX, 4);
    }
    XT_SSIP(A0_re,pR,4);
    XT_SSIP(A0_im,pR,4);
    XT_SSIP(A1_re,pR,4);
    XT_SSIP(A1_im,pR,4);
  }
} /* cxfir_xcorrf() */

#endif
