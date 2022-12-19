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
    Decimating block real FIR filter, floating point
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
#include "common_fpu.h"
#include "firdecf_common.h"

#if (HAVE_VFPU)

#define SMALLER_CODESIZE 0

/*-----------------------------------------------------------------------------
 * Data processing function of a particular decimating filter. Stores a
 * block of input samples to the circular delay line buffer and computes
 * decimating FIR filter's response.
 * Input:
 *   delayLine - circular delay line buffer start address
 *   delayLen  - Delay line buffer length
 *   wrIx    - next position in the buffer to be filled with an input sample
 *   x[N*D]  - input samples
 *   h[]     - decimating FIR filter coefficients, array layout varies
 * Output:
 *   y[N]    - output samples
 *   retval  - updated index of the oldest sample
 * Notes and restrictions:
 *   1. Most of data processing functions feature a single, hard-coded
 *      decimation factor, so they expect a determined value for parameter D.
 *   2. All pointers with the exception of y[N] must be aligned on an 16-bytes
 *      boundary.
 *   3. N - must be a multiple of 8.
 *   4. M - must be a multiple of 4.
 -----------------------------------------------------------------------------*/

int fir_decimaf_3x(float32_t * restrict y,
                   float32_t * delayLine, int delayLen,
             const float32_t * restrict x,
             const float32_t * restrict h,
             int wrIx, int D, int N, int M)
{
    const xtfloatx4  *          pX;
          xtfloatx4  * restrict pDw;
    const ae_int64   *          pDr;
    const xtfloatx4  *          pD0;
    const xtfloatx4  *          pD1;
    const xtfloatx2  *          pD2;
    const xtfloatx4  *          pH0;
    const xtfloatx4  *          pH1;
          xtfloatx4  * restrict pY;

    ae_valignx2 aY, aH1;

    xtfloatx2 q0, q1, q2, q3, q4, q5, q6, q7;
    xtfloatx2 q8, q9, qa, qb, qc, qd, qe, qf;
    xtfloatx2 d0, d1, d2, d3, d4, d5, d6, d7, d8;
    xtfloatx2 h0, h1, h2, h3;

    int m, n;

    static const ALIGN(16) int16_t Sel[4] = { 0x0705, 0x0604, 0x0301, 0x0200 };
    ae_int16x4 sel;
    sel = AE_L16X4_I((ae_int16x4*)&Sel, 0);

    NASSERT(y && delayLine && x && h);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(h, 16);
    NASSERT_ALIGN(delayLine, 16);
    NASSERT(D == 3);
    NASSERT(M > 0 && M % 4 == 0);
    NASSERT(N > 0 && N % 8 == 0);

    //
    // Setup pointers and circular delay line buffer.
    //
    pX  = (const xtfloatx4 *)x;
    pY  = (      xtfloatx4 *)y;
    pDw = (      xtfloatx4 *)(delayLine + wrIx);
    WUR_AE_CBEGIN0((uintptr_t)(delayLine));
    WUR_AE_CEND0  ((uintptr_t)(delayLine + delayLen));
    aY = AE_ZALIGN128();

    //
    // Break the input signal into 8*D-samples blocks. For each block, store
    // 8*D samples to the delay line buffer, and compute 8 samples of decimated
    // response signal.
    //
    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 2); n++)
    {
        pDr = (const ae_int64 *)pDw;
        AE_LSX2X2_IP(d0, d1, pX, 4 * sizeof(float32_t));
        AE_LSX2X2_IP(d2, d3, pX, 4 * sizeof(float32_t));
        AE_LSX2X2_IP(d4, d5, pX, 4 * sizeof(float32_t));
        AE_SSX2X2_XC(d0, d1, pDw, 4 * sizeof(float32_t));
        AE_SSX2X2_XC(d2, d3, pDw, 4 * sizeof(float32_t));
        AE_SSX2X2_XC(d4, d5, pDw, 4 * sizeof(float32_t));

        pH0 = (const xtfloatx4 *)h;
        pH1 = (const xtfloatx4 *)(h + 1);
        aH1 = AE_LA128_PP(pH1);
        AE_ADDCIRC_XC(pDr, 24 * sizeof(float32_t)); pD0 = (const xtfloatx4 *)pDr;
        AE_ADDCIRC_XC(pDr,  8 * sizeof(float32_t)); pD1 = (const xtfloatx4 *)pDr;
        AE_ADDCIRC_XC(pDr,  8 * sizeof(float32_t)); pD2 = (const xtfloatx2 *)pDr;

        CONST_SX2X2(q0, q1, 0); CONST_SX2X2(q2, q3, 0);
        CONST_SX2X2(q4, q5, 0); CONST_SX2X2(q6, q7, 0);
        CONST_SX2X2(q8, q9, 0); CONST_SX2X2(qa, qb, 0);
        CONST_SX2X2(qc, qd, 0); CONST_SX2X2(qe, qf, 0);

        AE_LSX2X2_XC(d0, d1, pD0, 4 * sizeof(float32_t));
        AE_LSX2X2_XC(d2, d3, pD0, 4 * sizeof(float32_t));
        AE_LSX2X2_XC(d4, d5, pD1, 4 * sizeof(float32_t));
        AE_LSX2X2_XC(d6, d7, pD1, 4 * sizeof(float32_t));
        #if SMALLER_CODESIZE
        __Pragma("no_unroll")
        #endif
        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 3) + 1; m++)
        {
            AE_LSX2XC(d8, pD2, 8 * sizeof(float32_t));

            AE_LSX2X2_IP(h0, h1, pH0, 4 * sizeof(float32_t));
            AE_LASX2X2_IP(h2, h3, aH1, pH1);
            MADD_SX2X2(q0, q1, d0, d2, h0, h2);
            MADD_SX2X2(q2, q3, d3, d5, h0, h2);
            MADD_SX2X2(q4, q5, d1, d3, h1, h3);
            MADD_SX2X2(q6, q7, d4, d6, h1, h3);

            AE_LSX2X2_IP(h0, h1, pH0, 4 * sizeof(float32_t));
            AE_LASX2X2_IP(h2, h3, aH1, pH1);
            MADD_SX2X2(q8, q9, d2, d4, h0, h2);
            MADD_SX2X2(qa, qb, d5, d7, h0, h2);
            MADD_SX2X2(qc, qd, d3, d5, h1, h3);
            MADD_SX2X2(qe, qf, d6, d8, h1, h3);

            d0 = d4; d1 = d5; d2 = d6; d3 = d7;
            AE_LSX2X2_XC(d4, d5, pD1, 4 * sizeof(float32_t));
            AE_LSX2X2_XC(d6, d7, pD1, 4 * sizeof(float32_t));
        }
        ADD_SX2X2(q0, q1, q0, q1, q4, q5);
        ADD_SX2X2(q2, q3, q2, q3, q6, q7);
        ADD_SX2X2(q8, q9, q8, q9, qc, qd);
        ADD_SX2X2(qa, qb, qa, qb, qe, qf);
        ADD_SX2X2(q0, q1, q0, q1, q8, q9);
        ADD_SX2X2(q2, q3, q2, q3, qa, qb);
        AE_DSELSX2(q0, q1, q0, q1, sel);
        AE_DSELSX2(q2, q3, q2, q3, sel);
        ADD_SX2X2(d0, d1, q0, q2, q1, q3);
        AE_SASX2X2_IP(d0, d1, aY, pY);
    }
    AE_SA128POS_FP(aY, pY);

    return (int)((float32_t *)pDw - delayLine);
}

#elif (HAVE_FPU)
// for scalar FPU
int fir_decimaf_3x(float32_t * restrict z,
                   float32_t * delay, int delayLen,
             const float32_t * restrict x,
             const float32_t * restrict h,
             int wrIx, int D, int N, int M)
{
  xtfloat* restrict pZ;
  const xtfloat *restrict pH;
  const xtfloat *restrict pX;
  xtfloat * pp;
  xtfloat* restrict pD;
  xtfloat x0, x1, x2, x3, x4, x5;
  xtfloat s0, s1, s2, s3, s4, s5;
  xtfloat h0, h1;
  int n, m;
  NASSERT(x);
  NASSERT(z);
  NASSERT(h);
  NASSERT(delay);
  NASSERT_ALIGN(delay, 8);
  NASSERT(N>0);
  NASSERT(M>0);
  NASSERT(M % 2 == 0);
  NASSERT(N % 8 == 0);
  /* set circular buffer boundaries */
  WUR_AE_CBEGIN0((uintptr_t)(delay + 0));
  WUR_AE_CEND0((uintptr_t)(delay + M));

  pp = (xtfloat*)(delay + wrIx);
  pZ = (xtfloat*)z;
  pX = (const xtfloat*)x;

  /* process by 4 input samples */
  for (n = 0; n<(N >> 1); n++)
  {
    xtfloat A0, A1, A2, A3, xx;
    pH = (const xtfloat*)h;
    pD = (xtfloat*)pp;
    A0 = A1 = A2 = A3 = XT_CONST_S(0);
    XT_LSIP(x0, pX, 4);
    XT_LSIP(x1, pX, 4);
    XT_LSIP(x2, pX, 4);
    XT_LSIP(x3, pX, 4);
    XT_LSIP(x4, pX, 4);
    XT_LSIP(x5, pX, 4);
    s0 = x0;
    s1 = x1;
    s2 = x2;
    s3 = x3;
    s4 = x4;
    s5 = x5;
    { xtfloat dummy;  XT_LSXC(dummy, pD, -4); }
    { xtfloat dummy;  XT_LSXC(dummy, pD, -4); xx = dummy; }
    for (m = 0; m<M; m += 2)
    {
      h0 = pH[M-m-0];
      h1 = pH[M-m-1];

      XT_MADD_S(A0, x0, h0);
      XT_MADD_S(A1, x3, h0);
      XT_MADD_S(A2, xx, h1);
      XT_MADD_S(A3, x2, h1);

      x3 = x1;
      x2 = x0;
      x1 = xx;
      XT_LSXC(x0, pD, -4);
      XT_LSXC(xx, pD, -4);
    }
    A0 = XT_ADD_S(A0, A2);
    A1 = XT_ADD_S(A1, A3);
    { xtfloat dummy;  XT_LSXC(dummy, pD, 4); }
    XT_SSXC(s0, pp, 4);
    XT_SSXC(s1, pp, 4);
    XT_SSXC(s2, pp, 4);
    XT_SSXC(s3, pp, 4);
    XT_SSXC(s4, pp, 4);
    XT_SSXC(s5, pp, 4);

    XT_SSIP(A0, pZ, 4);
    XT_SSIP(A1, pZ, 4);
  }
  return (int)((float32_t *)pp - delay);
}
#endif
