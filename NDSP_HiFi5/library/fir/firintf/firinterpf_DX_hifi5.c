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

/* Generic data processing function for a interpolating FIR filter. */
int fir_interpf_Dx( float32_t * restrict y,
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
          ae_int32  * restrict pY;

    ae_valignx2 aH1;

    xtfloatx2 q0, q1, q2, q3, q4, q5, q6, q7;
    xtfloatx2 q8, q9, qa, qb, qc, qd, qe, qf;
    xtfloatx2 d0, d1, d2, d3, d4, d5;
    xtfloatx2 h0, h1, h2, h3;

    int m, n, d;

    static const ALIGN(16) int16_t Sel[4] = { 0x0705, 0x0604, 0x0301, 0x0200 };
    ae_int16x4 sel;
    sel = AE_L16X4_I((ae_int16x4*)&Sel, 0);

    NASSERT(y && delayLine && x && h);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(h, 16);
    NASSERT_ALIGN(delayLine, 16);
    NASSERT(D > 1);
    NASSERT(M > 0 && M % 4 == 0);
    NASSERT(N > 0 && N % 8 == 0);

    //
    // Setup pointers and circular delay line buffer.
    //
    pX  = (const xtfloatx4 *)x;
    pY  = (      ae_int32  *)y;
    pDw = (      xtfloatx4 *)(delayLine + wrIx);
    WUR_AE_CBEGIN0((uintptr_t)(delayLine));
    WUR_AE_CEND0  ((uintptr_t)(delayLine + delayLen));

    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 3); n++)
    {
        AE_LSX2X2_IP(d0, d1, pX, 4 * sizeof(float32_t));
        AE_LSX2X2_IP(d2, d3, pX, 4 * sizeof(float32_t));
        AE_SSX2X2_XC(d0, d1, pDw, 4 * sizeof(float32_t));
        AE_SSX2X2_XC(d2, d3, pDw, 4 * sizeof(float32_t));

        pH0 = (const xtfloatx4 *)h;
        pH1 = (const xtfloatx4 *)(h + 1);
        aH1 = AE_LA128_PP(pH1);

        __Pragma("loop_count min=1");
        for (d = 0; d < D; d++)
        {
            pDr = pDw;
            S0 = pDr;
            AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(float32_t)); S1 = pDr;
            AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(float32_t)); S2 = pDr;

            CONST_SX2X2(q0, q1, 0); CONST_SX2X2(q2, q3, 0);
            CONST_SX2X2(q4, q5, 0); CONST_SX2X2(q6, q7, 0);
            CONST_SX2X2(q8, q9, 0); CONST_SX2X2(qa, qb, 0);
            CONST_SX2X2(qc, qd, 0); CONST_SX2X2(qe, qf, 0);

            __Pragma("loop_count min=2");
            for (m = 0; m < (M >> 2) + 1; m++)
            {
                AE_LSX2X2_IP(h0, h1, pH0, 4 * sizeof(float32_t));
                AE_LASX2X2_IP(h2, h3, aH1, pH1);
                AE_LSX2X2_XC(d0, d1, S0, 4 * sizeof(float32_t));
                AE_LSX2X2_XC(d2, d3, S1, 4 * sizeof(float32_t));
                AE_LSX2X2_XC(d4, d5, S2, 4 * sizeof(float32_t));
                MADD_SX2X2(q0, q8, d0, d1, h0, h1);
                MADD_SX2X2(q2, qa, d1, d2, h0, h1);
                MADD_SX2X2(q4, qc, d2, d3, h0, h1);
                MADD_SX2X2(q6, qe, d3, d4, h0, h1);
                MADD_SX2X2(q1, q9, d1, d2, h2, h3);
                MADD_SX2X2(q3, qb, d2, d3, h2, h3);
                MADD_SX2X2(q5, qd, d3, d4, h2, h3);
                MADD_SX2X2(q7, qf, d4, d5, h2, h3);
            }

            ADD_SX2X2(q0, q1, q0, q1, q8, q9);
            ADD_SX2X2(q2, q3, q2, q3, qa, qb);
            ADD_SX2X2(q4, q5, q4, q5, qc, qd);
            ADD_SX2X2(q6, q7, q6, q7, qe, qf);
            AE_DSELSX2(q0, q1, q0, q1, sel);
            AE_DSELSX2(q2, q3, q2, q3, sel);
            ADD_SX2X2(d0, d1, q0, q2, q1, q3);
            AE_DSELSX2(q4, q5, q4, q5, sel);
            AE_DSELSX2(q6, q7, q6, q7, sel);
            ADD_SX2X2(d2, d3, q4, q6, q5, q7);
            AE_S32_H_XP(AE_MOVINT32X2_FROMXTFLOATX2(d0), pY, D*sizeof(float32_t));
            AE_S32_L_XP(AE_MOVINT32X2_FROMXTFLOATX2(d0), pY, D*sizeof(float32_t));
            AE_S32_H_XP(AE_MOVINT32X2_FROMXTFLOATX2(d1), pY, D*sizeof(float32_t));
            AE_S32_L_XP(AE_MOVINT32X2_FROMXTFLOATX2(d1), pY, D*sizeof(float32_t));
            AE_S32_H_XP(AE_MOVINT32X2_FROMXTFLOATX2(d2), pY, D*sizeof(float32_t));
            AE_S32_L_XP(AE_MOVINT32X2_FROMXTFLOATX2(d2), pY, D*sizeof(float32_t));
            AE_S32_H_XP(AE_MOVINT32X2_FROMXTFLOATX2(d3), pY, D*sizeof(float32_t));
            AE_S32_L_XP(AE_MOVINT32X2_FROMXTFLOATX2(d3), pY, (1 - 7 * D)*sizeof(float32_t));
        }
        pY = (ae_int32 *)XT_ADDX4(7 * D, (uintptr_t)pY);
    }
    return (int)((float32_t *)pDw - delayLine);
}

#elif HAVE_FPU
// for scalar FPU
int fir_interpf_Dx( float32_t * restrict z,
                    float32_t * delay, int delayLen,
              const float32_t * restrict x,
              const float32_t * restrict h,
              int wrIx, int D, int N, int M )
{
#define P 4
  int n, m, j, _N;
  const xtfloat*  restrict pX = (const xtfloat*)x;
  const xtfloat*  restrict px = (const xtfloat*)x;
  const xtfloat* restrict pD  = (const xtfloat*)(delay + wrIx);
  const xtfloat*   restrict pH = (const xtfloat*)h;
  xtfloat*          pZ = (xtfloat*)z;
  ae_valign ax, az, ad;
  NASSERT(x);
  NASSERT(z);
  NASSERT(N>0);
  NASSERT(M>0);
  NASSERT(D>1);
  NASSERT(M % 4 == 0);
  NASSERT(N % 8 == 0);
  _N = N&(~3);
  ax = AE_LA64_PP(pX);
  az = AE_ZALIGN64();
  WUR_AE_CBEGIN0((uintptr_t)(delay));
  WUR_AE_CEND0((uintptr_t)(delay + M));
  for (n = 0; n<_N; n +=P)
  {
    xtfloat x0, x1, x2, x3, H0;
    xtfloat A0, A1, A2, A3;
    xtfloat s0, s1, s2, s3;
    pH = (const xtfloat*)h;
    __Pragma("loop_count min=5")
      for (j = 0; j<D; j++)
      {
        xtfloat temp;
        pX = (xtfloat*)((uintptr_t)(x + n));
        ax = AE_LA64_PP(pX);
        XT_LSIP(x0, castxcc(xtfloat, pX), 4);
        XT_LSIP(x1, castxcc(xtfloat, pX), 4);
        XT_LSIP(x2, castxcc(xtfloat, pX), 4);
        XT_LSIP(x3, castxcc(xtfloat, pX), 4);
        A0 =
        A1 =
        A2 =
        A3 = XT_CONST_S(0);
        ad = AE_LA64_PP(pD);
        XT_LSXC(temp, castxcc(xtfloat, pD), -4);

        __Pragma("loop_count min=1")
        for (m = 0; m<M; m ++)
        {
          XT_LSIP(H0, pH, 4);
          XT_MADD_S(A0, H0, x0);
          XT_MADD_S(A1, H0, x1);
          XT_MADD_S(A2, H0, x2);
          XT_MADD_S(A3, H0, x3);
          x3=x2;x2=x1;x1=x0;
          XT_LSXC(x0, castxcc(xtfloat, pD), -4);

        }
        XT_LSXC(temp, castxcc(xtfloat, pD), 4);
        pZ = (xtfloat*)((uintptr_t*)z + n*D + j);
        XT_SSXP(A0, pZ, 4 * D);
        XT_SSXP(A1, pZ, 4 * D);
        XT_SSXP(A2, pZ, 4 * D);
        XT_SSXP(A3, pZ, -3 * 4 * D + 4);
      }

      XT_LSIP(s0, castxcc(xtfloat, px), 4);
      XT_LSIP(s1, castxcc(xtfloat, px), 4);
      XT_LSIP(s2, castxcc(xtfloat, px), 4);
      XT_LSIP(s3, castxcc(xtfloat, px), 4);

      XT_SSXC(s0, castxcc(xtfloat, pD), 4);
      XT_SSXC(s1, castxcc(xtfloat, pD), 4);
      XT_SSXC(s2, castxcc(xtfloat, pD), 4);
      XT_SSXC(s3, castxcc(xtfloat, pD), 4);
    }
  return (int)((float32_t *)pD - delay);
}
#endif
