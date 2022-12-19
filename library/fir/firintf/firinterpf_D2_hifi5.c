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
    Interpolating block real FIR filter, floating point
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
#include "common_fpu.h"
#include "firinterpf_common.h"

#if (HAVE_VFPU)

#define SMALLER_CODESIZE 1
/*-----------------------------------------------------------------------------
 * Data processing function of a particular interpolating filter. Stores a
 * block of N input samples to the circular delay line buffer and computes
 * N*D samples of interpolating FIR filter's response.
 * Input:
 *   delayLine - circular delay line buffer start address
 *   delayLen  - Delay line buffer length
 *   wrIx    - next position in the buffer to be filled with an input sample
 *   x[N]    - input samples
 *   h[]     - decimating FIR filter coefficients, array layout varies
 * Output:
 *   y[N*D]  - output samples
 *   retval  - updated index of the oldest sample
 * Notes and restrictions:
 *   1. Most of data processing functions feature a single, hard-coded
 *      interpolation factor, so they expect a determined value for parameter D.
 *   2. All pointers with the exception of y[N] must be aligned on an 16-bytes
 *      boundary.
 *   3. N - must be a multiple of 8.
 *   4. M - must be a multiple of 4.
 -----------------------------------------------------------------------------*/

/* Data processing function for a factor 2 interpolating FIR filter. */
int fir_interpf_2x( float32_t * restrict y,
                    float32_t * delayLine, int delayLen,
              const float32_t * restrict x,
              const float32_t * restrict h,
              int wrIx, int D, int N, int M )
{
    const xtfloatx4 * restrict pX;
          xtfloatx4 * restrict pDw;
    const xtfloatx4 * restrict pDr;
    const xtfloatx4 * restrict S0;
    const xtfloatx4 * restrict S1;
    const xtfloatx4 * restrict S2;
    const xtfloatx4 * restrict pH0;
    const xtfloatx4 * restrict pH1;
    const xtfloatx4 * restrict pH2;
    const xtfloatx4 * restrict pH3;
          xtfloatx4  * restrict pY;

    ae_valignx2 aY;
    ae_valignx2 aH1, aH3;

    xtfloatx2 q0, q1, q2, q3, q4, q5, q6, q7;
    xtfloatx2 q8, q9, qa, qb, qc, qd, qe, qf;
    xtfloatx2 d0, d1, d2, d3, d4, d5;
    xtfloatx2 h0, h1, h2, h3, h4, h5, h6, h7;

    int m, n;

    static const ALIGN(16) int16_t Sel[4] = { 0x0705, 0x0604, 0x0301, 0x0200 };
    ae_int16x4 sel;
    sel = AE_L16X4_I((ae_int16x4*)&Sel, 0);

    NASSERT(y && delayLine && x && h);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(h, 16);
    NASSERT_ALIGN(delayLine, 16);
    NASSERT(D == 2);
    NASSERT(M > 0 && M % 4 == 0);
    NASSERT(N > 0 && N % 8 == 0);

    //
    // Setup pointers and circular delay line buffer.
    //
    pX  = (const xtfloatx4 *)x;
    pY  = (      xtfloatx4  *)y;
    pDw = (      xtfloatx4 *)(delayLine + wrIx);
    WUR_AE_CBEGIN0((uintptr_t)(delayLine));
    WUR_AE_CEND0  ((uintptr_t)(delayLine + delayLen));
    aY = AE_ZALIGN128();

    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 3); n++)
    {
        pDr = pDw;
        AE_LSX2X2_IP(d0, d1, pX, 4 * sizeof(float32_t));
        AE_LSX2X2_IP(d2, d3, pX, 4 * sizeof(float32_t));
        AE_SSX2X2_XC(d0, d1, pDw, 4 * sizeof(float32_t));
        AE_SSX2X2_XC(d2, d3, pDw, 4 * sizeof(float32_t));

        pH0 = (const xtfloatx4 *)h;
        pH1 = (const xtfloatx4 *)(h + 1);
        pH2 = (const xtfloatx4 *)(h + (M + 4));
        pH3 = (const xtfloatx4 *)(h + (M + 4) + 1);
        aH1 = AE_LA128_PP(pH1);
        aH3 = AE_LA128_PP(pH3);

        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 8 * sizeof(float32_t)); S0 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(float32_t)); S1 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(float32_t)); S2 = pDr;

        CONST_SX2X2(q0, q1, 0); CONST_SX2X2(q2, q3, 0);
        CONST_SX2X2(q4, q5, 0); CONST_SX2X2(q6, q7, 0);
        CONST_SX2X2(q8, q9, 0); CONST_SX2X2(qa, qb, 0);
        CONST_SX2X2(qc, qd, 0); CONST_SX2X2(qe, qf, 0);
#if SMALLER_CODESIZE
        __Pragma("no_unroll");
#endif
        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            AE_LSX2X2_IP(h0, h1, pH0, 4 * sizeof(float32_t));
            AE_LASX2X2_IP(h2, h3, aH1, pH1);
            AE_LSX2X2_IP(h4, h5, pH2, 4 * sizeof(float32_t));
            AE_LASX2X2_IP(h6, h7, aH3, pH3);

            AE_LSX2X2_XC(d0, d1, S0, 4 * sizeof(float32_t));
            AE_LSX2X2_XC(d2, d3, S1, 4 * sizeof(float32_t));
            AE_LSX2X2_XC(d4, d5, S2, 4 * sizeof(float32_t));

            MADD_SX2X2(q0, q1, d0, d0, h0, h4);
            MADD_SX2X2(q2, q3, d1, d1, h2, h6);
            MADD_SX2X2(q4, q5, d1, d1, h0, h4);
            MADD_SX2X2(q6, q7, d2, d2, h2, h6);
            MADD_SX2X2(q8, q9, d2, d2, h0, h4);
            MADD_SX2X2(qa, qb, d3, d3, h2, h6);
            MADD_SX2X2(qc, qd, d3, d3, h0, h4);
            MADD_SX2X2(qe, qf, d4, d4, h2, h6);

            MADD_SX2X2(q0, q1, d1, d1, h1, h5);
            MADD_SX2X2(q2, q3, d2, d2, h3, h7);
            MADD_SX2X2(q4, q5, d2, d2, h1, h5);
            MADD_SX2X2(q6, q7, d3, d3, h3, h7);
            MADD_SX2X2(q8, q9, d3, d3, h1, h5);
            MADD_SX2X2(qa, qb, d4, d4, h3, h7);
            MADD_SX2X2(qc, qd, d4, d4, h1, h5);
            MADD_SX2X2(qe, qf, d5, d5, h3, h7);
        }

        AE_DSELSX2(q0, q1, q0, q1, sel);
        AE_DSELSX2(q2, q3, q2, q3, sel);
        ADD_SX2X2(d0, d1, q0, q2, q1, q3);
        AE_SASX2X2_IP(d0, d1, aY, pY);

        AE_DSELSX2(q4, q5, q4, q5, sel);
        AE_DSELSX2(q6, q7, q6, q7, sel);
        ADD_SX2X2(d0, d1, q4, q6, q5, q7);
        AE_SASX2X2_IP(d0, d1, aY, pY);

        AE_DSELSX2(q8, q9, q8, q9, sel);
        AE_DSELSX2(qa, qb, qa, qb, sel);
        ADD_SX2X2(d0, d1, q8, qa, q9, qb);
        AE_SASX2X2_IP(d0, d1, aY, pY);

        AE_DSELSX2(qc, qd, qc, qd, sel);
        AE_DSELSX2(qe, qf, qe, qf, sel);
        ADD_SX2X2(d0, d1, qc, qe, qd, qf);
        AE_SASX2X2_IP(d0, d1, aY, pY);
    }
    AE_SA128POS_FP(aY, pY);
    return (int)((float32_t *)pDw - delayLine);
}

#elif HAVE_FPU
// for scalar FPU
int fir_interpf_2x( float32_t * restrict z,
                    float32_t * delay, int delayLen,
              const float32_t * restrict x,
              const float32_t * restrict h,
              int wrIx, int D, int N, int M )
{
  int n, m;
  const xtfloat*  restrict pX = (const xtfloat*)x;
  const xtfloat*  restrict px = (const xtfloat*)x;
  const xtfloat* restrict pD  = (const xtfloat*)(delay + wrIx);
  const xtfloat*   restrict pH = (const xtfloat*)h;
  xtfloat*          pZ = (xtfloat*)z;
  NASSERT(x);
  NASSERT(z);
  NASSERT(N>0);
  NASSERT(M>0);
  NASSERT(M % 4 == 0);
  NASSERT(N % 8 == 0);
  WUR_AE_CBEGIN0((uintptr_t)(delay));
  WUR_AE_CEND0((uintptr_t)(delay + M));
  for (n = 0; n<N; n +=2)
  {
    xtfloat x0, x1;
    xtfloat H0, H1;
    xtfloat A0, A1, A2, A3;
    xtfloat s0, s1;
    pH = (const xtfloat*)h;

    {
      xtfloat temp;
      XT_LSIP(x0, castxcc(xtfloat, pX), 4);
      XT_LSIP(x1, castxcc(xtfloat, pX), 4);
      A0 = A1 = XT_CONST_S(0);
      A2 = A3 = XT_CONST_S(0);
      XT_LSXC(temp, castxcc(xtfloat, pD), -4);
      __Pragma("loop_count min=1")
        for (m = 0; m<M; m ++)
        {
          H0 = pH[m + 0*M];
          H1 = pH[m + 1*M];
          XT_MADD_S(A0, H0, x0);
          XT_MADD_S(A1, H0, x1);
          XT_MADD_S(A2, H1, x0);
          XT_MADD_S(A3, H1, x1);
          x1 = x0;
          XT_LSXC(x0, castxcc(xtfloat, pD), -4);
        }

      XT_LSXC(temp, castxcc(xtfloat, pD), 4);
      XT_SSXP(A0, pZ, 4 );
      XT_SSXP(A2, pZ, 4);
      XT_SSXP(A1, pZ, 4);
      XT_SSXP(A3, pZ, 4);
    }
    XT_LSIP(s0, castxcc(xtfloat, px), 4);
    XT_LSIP(s1, castxcc(xtfloat, px), 4);
    XT_SSXC(s0, castxcc(xtfloat, pD), 4);
    XT_SSXC(s1, castxcc(xtfloat, pD), 4);
  }
  return (int)((float32_t *)pD - delay);
}
#endif
