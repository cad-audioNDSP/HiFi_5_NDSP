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
    Real data circular cross-correlation, 32x16-bit, no requirements on vectors
    length and alignment.
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "raw_corr32x16.h"

/*-------------------------------------------------------------------------
  Circular Correlation
  Estimates the circular cross-correlation between vectors x (of length N) 
  and y (of length M)  resulting in vector r of length N. It is a similar 
  to correlation but x is read in opposite direction.
  These functions implement the circular correlation algorithm with no 
  limitations on x and y vectors length and alignment at the cost of 
  increased processing complexity. In addition, this implementation variant
  requires scratch memory area.

  Precision: 
  16x16     16x16-bit data, 16-bit outputs
  32x16     32x16-bit data, 32-bit outputs
  32x32     32x32-bit data, 32-bit outputs
  32x32ep   the same as above but using 72-bit accumulator for intermediate 
            computations
  f         floating point

  Input:
  s[]       scratch area, 
              FIR_XCORRA16X16_SCRATCH_SIZE( N, M )
              FIR_XCORRA32X16_SCRATCH_SIZE( N, M ) or
              FIR_XCORRA32X32_SCRATCH_SIZE( N, M ) or
              FIR_XCORRA32X32EP_SCRATCH_SIZE( N, M ) or
              FIR_XCORRAF_SCRATCH_SIZE( N, M ) or
              CXFIR_XCORRA32X32_SCRATCH_SIZE( N, M ) or
              CXFIR_XCORRAF_SCRATCH_SIZE( N, M ) bytes

  x[N]      input data Q15, Q31 or floating point
  y[M]      input data Q15, Q31 or floating point
  N         length of x
  M         length of y
  Output:
  r[N]      output data, Q15, Q31 or floating point

  Restrictions:
  x,y,r,s   should not overlap
  s         must be aligned on an 16-bytes boundary
  N,M       must be >0
  N >= M-1  minimum allowed length of vector x is the length of y minus one

-------------------------------------------------------------------------*/

void fir_xcorra32x16( void    * restrict s,
                      int32_t * restrict r,
                const int32_t * restrict x,
                const int16_t * restrict y,
                int N, int M )
{
          ae_int32x4 * restrict S;
    const ae_int32x4 * restrict X;
    const ae_int16x8 * restrict Y;
    ae_valignx2 as, ax, ay;
    ae_int32x2 t0, t1;
    ae_int16x4 p0, p1;

    //
    // Circular cross-correlation algorithm:
    //
    //   r[n] = sum( x[mod(n+m,N)]*y[m] )
    //        m=0..M-1
    //
    //   where n = 0..N-1
    //
    int32_t * x_buf;
    int16_t * y_buf;

    int n, m;

    NASSERT(s && r && x && y && N > 0 && M > 0 && N >= M - 1);
    NASSERT_ALIGN(s, 16);

    // Partition the scratch memory area.
    x_buf = (int32_t*)s;
    y_buf = (int16_t*)(x_buf + (((N + M - 1) + 1)&~1)); // closest bigger aligned

    // Copy x[N] data into the scratch memory in a way that simplifies the
    // correlation calculation:
    // x[0]..x[N-1] x[0]..x[M-2]
    S = (      ae_int32x4*)x_buf;
    X = (const ae_int32x4*)x;
    ax = AE_LA128_PP(X);
    for (n = 0; n < (N >> 2); n++)
    {
        AE_LA32X2X2_IP(t0, t1, ax, X);
        AE_S32X2X2_IP(t0, t1, S, 4 * sizeof(int32_t));
    }
    for (n = 0; n < (N & 3); n++)
    {
        AE_L32_IP(t0, castxcc(ae_int32, X), sizeof(int32_t));
        AE_S32_L_IP(t0, castxcc(ae_int32, S), sizeof(int32_t));
    }
    X = (const ae_int32x4*)x;
    for (n = 0; n < ((M - 1) & 3); n++)
    {
        AE_L32_IP(t0, castxcc(ae_int32, X), sizeof(int32_t));
        AE_S32_L_IP(t0, castxcc(ae_int32, S), sizeof(int32_t));
    }
    ax = AE_LA128_PP(X);
    as = AE_ZALIGN128();
    for (n = 0; n < ((M - 1) >> 2); n++)
    {
        AE_LA32X2X2_IP(t0, t1, ax, X);
        AE_SA32X2X2_IP(t0, t1, as, S);
    }
    AE_SA128POS_FP(as, S);
    AE_S32_L_IP(AE_ZERO32(), castxcc(ae_int32, S), sizeof(int32_t));


    S = (ae_int32x4*)(y_buf);
    NASSERT_ALIGN(S, 8);
    Y = (const ae_int16x8*)y;
    ay = AE_LA128_PP(Y);
    as = AE_ZALIGN128();
    for (m = 0; m < (M >> 3); m++)
    {
        AE_LA16X4X2_IP(p0, p1, ay, Y);
        AE_SA16X4X2_IP(p0, p1, as, castxcc(ae_int16x8, S));
    }
    AE_SA128POS_FP(as, S);
    if (M & 4)
    {
        ae_valign ay1;
        ay1 = AE_LA64_PP(Y);
        AE_LA16X4_IP(p0, ay1, castxcc(ae_int16x4, Y));
        AE_S16X4_IP(p0, castxcc(ae_int16x4, S), 4 * sizeof(int16_t));
    }
    if (M & 3)
    {
        AE_S16X4_I(AE_ZERO16(), (ae_int16x4 *)S, 0);
        __Pragma("loop_count min=1");
        for (m = 0; m < (M & 3); m++)
        {
            AE_L16_IP(p0, castxcc(ae_int16, Y), sizeof(int16_t));
            AE_S16_0_IP(p0, castxcc(ae_int16, S), sizeof(int16_t));
        }
    }

    raw_corr32x16(r, x_buf, y_buf, N, M);
} // fir_xcorra32x16()
