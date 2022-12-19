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
    Complex block FIR filter, floating point
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
#include "cxfirf.h"

/*-------------------------------------------------------------------------
  Complex Block FIR Filter
  Computes a complex FIR filter (direct-form) using complex IR stored in 
  vector h. The complex data input is stored in vector x. The filter output
  result is stored in vector r. The filter calculates N output samples using 
  M coefficients and requires last M-1 samples in the delay line. Real and
  imaginary parts are interleaved and real parts go first (at even indexes).
  NOTE: 
  1. User application is not responsible for management of delay lines
  2. User has an option to set IR externally or copy from original location 
     (i.e. from the slower constant memory). In the first case, user is 
     responsible for right alignment, ordering and zero padding of filter 
     coefficients - usually array is composed from zeroes (left padding), 
     reverted IR and right zero padding.

  Precision: 
  16x16     16-bit data, 16-bit coefficients, 16-bit outputs
  32x16     32-bit data, 16-bit coefficients, 32-bit outputs
  32x32     32-bit data, 32-bit coefficients, 32-bit outputs
  32x32ep   the same as above but using 72-bit accumulator for intermediate 
            computations
  f         floating point

  Input:
  h[M]      complex filter coefficients; h[0] is to be multiplied with the 
            newest sample, Q15, Q31, floating point
  x[N]      input samples, Q15, Q31, floating point
  N         length of sample block (in complex samples) 
  M         length of filter 
  extIR     if zero, IR is copied from original location, otherwise not
            but user should keep alignment, order of coefficients 
            and zero padding requirements shown below
  Output:			
  y[N]      output samples, Q15, Q31, floating point

  Alignment, ordering and zero padding for external IR  (extIR!=0)
  -----------------+----------+--------------+--------------+----------------
  Function	       |Alignment,|Left zero     |   Coefficient| Right zero 
                   | bytes    |padding, bytes|   order      | padding, bytes
  -----------------+----------+--------------+--------------+----------------
  cxfir16x16_init, |    16    |  2 before    |  *           |  6 after
                   |          |  each copy   |              |  each copy
  cxfir32x16_init  |    16    |  2 before    |  *           |  6 after
                   |          |  each copy   |              |  each copy
  cxfir32x32_init  |    16    |    8         |  inverted    |  0
  cxfir32x32ep_init|    16    |    0         |  inv,conj    |  0
  cxfirf_init      |    16    |    0         |  direct      |  0
  -----------------+----------+--------------+--------------+----------------
  * inverted: conjugated copy and (imaginary; real) copy at 4*(M+4) bytes offset

  Restriction:
  x,y       should not overlap
  x,h       aligned on a 16-bytes boundary
  N,M       multiples of 4
-------------------------------------------------------------------------*/

/* Circular load with using CBEGIN1/CEND1 */
#define XT_LSX2XC1(reg, addr, offs)\
{\
    ae_int32x2 t;\
    AE_L32X2_XC1(t, addr, offs);\
    reg = XT_AE_MOVXTFLOATX2_FROMINT32X2(t);\
}

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void,cxfirf_process,( cxfirf_handle_t handle,
                         complex_float * restrict  y,
                   const complex_float * restrict  x, int N ))
#elif HAVE_VFPU

/* process block of samples */
void cxfirf_process( cxfirf_handle_t handle,
                         complex_float * restrict  y,
                   const complex_float * restrict  x, int N )
{
#if 1
  cxfirf_ptr_t cxfir = (cxfirf_ptr_t)handle;

  const xtfloatx4 * restrict pX;
  xtfloatx4 * restrict pDw;
  const xtfloatx2 * restrict S0;
  const xtfloatx4 * restrict S1;
  const xtfloatx4 * restrict S2;
  const xtfloatx4 * restrict pH;
  xtfloatx4 * restrict pY;
  ae_valignx2 aY;
  xtfloatx2 r00, r10, r20, r30;
  xtfloatx2 r01, r11, r21, r31;
  xtfloatx2 r02, r12, r22, r32;
  xtfloatx2 r03, r13, r23, r33;
  xtfloatx2 x0, x1, x2, x3, x4;
  xtfloatx2 x5, x6, x7, x8, x9, t0;
  xtfloatx2 y0, y1;

  int M;
  int n, m;

  M = cxfir->M;
  NASSERT(cxfir && cxfir->magic == CXFIRF_MAGIC && y && x);
  NASSERT_ALIGN(x, 16);
  NASSERT_ALIGN((cxfir->coef), 16);
  NASSERT_ALIGN((cxfir->delayLine), 16);
  NASSERT_ALIGN((cxfir->delayLine + cxfir->delayLen), 16);
  NASSERT_ALIGN((cxfir->delayPos), 16);
  NASSERT(N % 4 == 0 && M % 4 == 0);
  if (N <= 0) return;

  // Setup pointers and circular delay line buffer.
  pX = (const xtfloatx4 *)x;
  pY = (      xtfloatx4 *)y;
  pDw= (      xtfloatx4 *)cxfir->delayPos;
  WUR_AE_CBEGIN0((uintptr_t)(cxfir->delayLine));
  WUR_AE_CEND0  ((uintptr_t)(cxfir->delayLine + cxfir->delayLen));
  /* set bounds of the IR coeffs */
  WUR_AE_CBEGIN1((uintptr_t)(cxfir->coef));
  WUR_AE_CEND1((uintptr_t)(cxfir->coef + M));
  aY = AE_ZALIGN128();
  pH = (const xtfloatx4 *)(cxfir->coef + M - 2);

  for (n = 0; n < (N >> 3); n++)
  {
    AE_LSX2X2_IP(x0, x1, pX, 2 * sizeof(complex_float));
    AE_LSX2X2_IP(x2, x3, pX, 2 * sizeof(complex_float));
    AE_LSX2X2_IP(x4, x5, pX, 2 * sizeof(complex_float));
    AE_LSX2X2_IP(x6, x7, pX, 2 * sizeof(complex_float));
    AE_SSX2X2_XC(x0, x1, pDw, 2 * sizeof(complex_float));
    AE_SSX2X2_XC(x2, x3, pDw, 2 * sizeof(complex_float));
    AE_SSX2X2_XC(x4, x5, pDw, 2 * sizeof(complex_float));
    AE_SSX2X2_XC(x6, x7, pDw, 2 * sizeof(complex_float));

    S1 =  S2 = pDw;
    AE_ADDCIRC32X2_XC(castxcc(ae_int32x2, S2), 5 * 2 *(int)sizeof(complex_float));
    /* preload data from x */
    AE_LSX2X2_XC(x7, x0, S1, 2 * sizeof(complex_float));
    AE_LSX2X2_XC(x1, x2, S1, 2 * sizeof(complex_float));
    AE_LSX2X2_XC(x3, x4, S1, 2 * sizeof(complex_float));
    AE_LSX2X2_XC(x5, x6, S1, 2 * sizeof(complex_float));
    AE_LSX2X2_XC(x7, x8, S1, 2 * sizeof(complex_float));

    CONST_SX2X2(r00, r01, 0); CONST_SX2X2(r02, r03, 0);
    CONST_SX2X2(r10, r11, 0); CONST_SX2X2(r12, r13, 0);
    CONST_SX2X2(r20, r21, 0); CONST_SX2X2(r22, r23, 0);
    CONST_SX2X2(r30, r31, 0); CONST_SX2X2(r32, r33, 0);

    __Pragma("loop_count min=2 factor=2")
    for (m = 0; m < (M >> 1); m++) //(8 cycles, unroll 1, MAC bound)
    {
      /* load data from y */
      AE_LSX2X2_XC1(y1, y0, castxcc(xtfloatx4, pH), -2 * (int)sizeof(complex_float));

      /* compute correlation of 8 values */
      MADDMUXQ_S(r00, r10, x0, x1, y0, 0);
      MADDMUXQ_S(r01, r11, x0, x1, y0, 1); 
      MADDMUXQ_S(r20, r30, x2, x3, y0, 0);
      MADDMUXQ_S(r21, r31, x2, x3, y0, 1);

      MADDMUXQ_S(r02, r12, x4, x5, y0, 0);
      MADDMUXQ_S(r03, r13, x4, x5, y0, 1);
      MADDMUXQ_S(r22, r32, x6, x7, y0, 0);
      MADDMUXQ_S(r23, r33, x6, x7, y0, 1);

      /* reload data for the next iteration */
      AE_LSX2X2_XC(x9, t0, S2, 2 * sizeof(complex_float));

      /* compute correlation of 8 values */
      MADDMUXQ_S(r00, r10, x1, x2, y1, 0);
      MADDMUXQ_S(r01, r11, x1, x2, y1, 1);
      MADDMUXQ_S(r20, r30, x3, x4, y1, 0);
      MADDMUXQ_S(r21, r31, x3, x4, y1, 1);

      MADDMUXQ_S(r02, r12, x5, x6, y1, 0);
      MADDMUXQ_S(r03, r13, x5, x6, y1, 1);
      MADDMUXQ_S(r22, r32, x7, x8, y1, 0);
      MADDMUXQ_S(r23, r33, x7, x8, y1, 1);

      x0 = x2; x1 = x3; x2 = x4; x3 = x5;
      x4 = x6; x5 = x7; x6 = x8; x7 = x9;
      x8 = t0;
    }
    /* reduction add */
    ADD_SX2X2(r00, r10, r00, r10, r01, r11);
    ADD_SX2X2(r20, r30, r20, r30, r21, r31);

    ADD_SX2X2(r02, r12, r02, r12, r03, r13);
    ADD_SX2X2(r22, r32, r22, r32, r23, r33);
      
    /* save computed samples */
    AE_SASX2X2_IP(r00, r10, aY, pY);
    AE_SASX2X2_IP(r20, r30, aY, pY);
    AE_SASX2X2_IP(r02, r12, aY, pY);
    AE_SASX2X2_IP(r22, r32, aY, pY);
  }
  AE_SA128POS_FP(aY, pY);
  if (N & 4)
  {
    pH = (const xtfloatx4 *)(cxfir->coef + M - 1);
    AE_LSX2X2_IP(x0, x1, pX, 2 * sizeof(complex_float));
    AE_LSX2X2_IP(x2, x3, pX, 2 * sizeof(complex_float));
    AE_SSX2X2_XC(x0, x1, pDw, 2 * sizeof(complex_float));
    AE_SSX2X2_XC(x2, x3, pDw, 2 * sizeof(complex_float));

    S1 = pDw;
    /* preload data from x */
    AE_ADDCIRC32X2_XC(castxcc(ae_int32x2, S1), 2 * 2 * (int)sizeof(complex_float));
    AE_LSX2X2_XC(x7, x0, S1, 2 * sizeof(complex_float));
    AE_LSX2X2_XC(x1, x2, S1, 2 * sizeof(complex_float));
    S0 = (const xtfloatx2 *)S1;
    AE_LSX2XC(x3, S0, sizeof(complex_float));

    CONST_SX2X2(r00, r01, 0); CONST_SX2X2(r02, r03, 0);
    CONST_SX2X2(r10, r11, 0); CONST_SX2X2(r12, r13, 0);
    CONST_SX2X2(r20, r21, 0); CONST_SX2X2(r22, r23, 0);
    CONST_SX2X2(r30, r31, 0); CONST_SX2X2(r32, r33, 0);

    __Pragma("loop_count min=4 factor=4")
    for (m = 0; m < (M ); m++)
    {
      /* load data from y */
      XT_LSX2XC1(y0, castxcc(ae_int32x2, pH), -1*(int)sizeof(complex_float));

      /* compute correlation of 4 values */
      MADDMUXQ_S(r00, r10, x0, x1, y0, 0);
      MADDMUXQ_S(r01, r11, x0, x1, y0, 1);
      MADDMUXQ_S(r20, r30, x2, x3, y0, 0);
      MADDMUXQ_S(r21, r31, x2, x3, y0, 1);

      x0 = x1; x1 = x2; x2 = x3; x3 = x4;
      /* reload data for the next iteration */
      AE_LSX2XC(x3, S0, sizeof(complex_float));
    }
    /* reduction add */
    ADD_SX2X2(r00, r10, r00, r10, r01, r11);
    ADD_SX2X2(r20, r30, r20, r30, r21, r31);

    /* save computed samples */
    AE_SASX2X2_IP(r00, r10, aY, pY);
    AE_SASX2X2_IP(r20, r30, aY, pY);
  }
  AE_SA128POS_FP(aY, pY);
  cxfir->delayPos = (complex_float*)pDw;

#else
    cxfirf_ptr_t cxfir = (cxfirf_ptr_t)handle;

    const xtfloatx4 * restrict pX;
          xtfloatx4 * restrict pDw;
    const xtfloatx2 * restrict S0;
    const xtfloatx4 * restrict S1;
    const xtfloatx4 * restrict S2;
    const xtfloatx4 * restrict pH;
          xtfloatx4 * restrict pY;
    ae_valignx2 aY;
    xtfloatx2 r00, r10, r20, r30;
    xtfloatx2 r01, r11, r21, r31;
    xtfloatx2 r02, r12, r22, r32;
    xtfloatx2 r03, r13, r23, r33;
    xtfloatx2 x0, x1, x2, x3, x4;
    xtfloatx2 y0, y1;

    int M;
    int n, m;

    M = cxfir->M;
    NASSERT(cxfir && cxfir->magic == CXFIRF_MAGIC && y && x);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN((cxfir->coef), 16);
    NASSERT_ALIGN((cxfir->delayLine), 16);
    NASSERT_ALIGN((cxfir->delayLine + cxfir->delayLen), 16);
    NASSERT_ALIGN((cxfir->delayPos), 16);
    NASSERT(N % 4 == 0 && M % 4 == 0);
    if (N <= 0) return;

    // Setup pointers and circular delay line buffer.
    pX = (const xtfloatx4 *)x;
    pY = (      xtfloatx4 *)y;
    pDw= (      xtfloatx4 *)cxfir->delayPos;
    WUR_AE_CBEGIN0((uintptr_t)(cxfir->delayLine));
    //!!!!WUR_AE_CEND0  ((uintptr_t)(cxfir->delayLine + cxfir->delayLen));
    WUR_AE_CEND0  ((uintptr_t)(cxfir->delayLine + cxfir->delayLen - 4));
    aY = AE_ZALIGN128();

    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 2); n++)
    {
        AE_LSX2X2_IP(x0, x1, pX, 2 * sizeof(complex_float));
        AE_LSX2X2_IP(x2, x3, pX, 2 * sizeof(complex_float));
        AE_SSX2X2_XC(x0, x1, pDw, 2 * sizeof(complex_float));
        AE_SSX2X2_XC(x2, x3, pDw, 2 * sizeof(complex_float));

        pH = (const xtfloatx4 *)(cxfir->coef + M - 2);
        S0 = (const xtfloatx2 *)XT_ADDX8(1, (uintptr_t)pDw);
        S1 = (const xtfloatx4 *)XT_ADDX8(2, (uintptr_t)pDw);
        S2 = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, S2), 4 * sizeof(complex_float));
        /* preload data from x */
        AE_LSX2XC(x0, S0, 2 * sizeof(complex_float));

        CONST_SX2X2(r00, r01, 0); CONST_SX2X2(r02, r03, 0);
        CONST_SX2X2(r10, r11, 0); CONST_SX2X2(r12, r13, 0);
        CONST_SX2X2(r20, r21, 0); CONST_SX2X2(r22, r23, 0);
        CONST_SX2X2(r30, r31, 0); CONST_SX2X2(r32, r33, 0);

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 1); m++)
        {
            /* load data from x */
            AE_LSX2X2_XC(x1, x2, S1, 2 * sizeof(complex_float));
            AE_LSX2X2_XC(x3, x4, S2, 2 * sizeof(complex_float));
            /* load data from y */
            AE_LSX2X2_IP(y1, y0, pH, -2 * (int)sizeof(complex_float));
            /* compute correlation of 4 values */
            MADDMUX_SX2X2(r00, r10, y0, y0, x0, x1, 0);
            MADDMUX_SX2X2(r01, r11, y0, y0, x0, x1, 1);
            MADDMUX_SX2X2(r20, r30, y0, y0, x2, x3, 0);
            MADDMUX_SX2X2(r21, r31, y0, y0, x2, x3, 1);
            MADDMUX_SX2X2(r02, r12, y1, y1, x1, x2, 0);
            MADDMUX_SX2X2(r03, r13, y1, y1, x1, x2, 1);
            MADDMUX_SX2X2(r22, r32, y1, y1, x3, x4, 0);
            MADDMUX_SX2X2(r23, r33, y1, y1, x3, x4, 1);
            /* reload data for the next iteration */
            AE_LSX2XC(x0, S0, 2 * sizeof(complex_float));
        }
        /* reduction add */
        ADD_SX2X2(r00, r10, r00, r10, r01, r11);
        ADD_SX2X2(r20, r30, r20, r30, r21, r31);
        ADD_SX2X2(r02, r12, r02, r12, r03, r13);
        ADD_SX2X2(r22, r32, r22, r32, r23, r33);
        ADD_SX2X2(r00, r10, r00, r10, r02, r12);
        ADD_SX2X2(r20, r30, r20, r30, r22, r32);
        /* save computed samples */
        AE_SASX2X2_IP(r00, r10, aY, pY);
        AE_SASX2X2_IP(r20, r30, aY, pY);
    }
    AE_SA128POS_FP(aY, pY);

    cxfir->delayPos = (complex_float*)pDw;
#endif
} // cxfirf_process()
#else

/* process block of samples */
void cxfirf_process( cxfirf_handle_t _cxfir,
                         complex_float * restrict  y,
                   const complex_float * restrict  x, int N )
{
  int n, m, M;
  xtfloat x0, x1, x2, x3;
  xtfloat s0, s1, s2, s3;
  xtfloat acc0, acc1, acc2, acc3;
  const xtfloat* restrict pX = (const xtfloat*)x;
        xtfloat* restrict pD;
  const xtfloat* restrict pH;
        xtfloat*          pZ = (xtfloat*)y;
  cxfirf_t* state;
  NASSERT(_cxfir);
  state=(cxfirf_t*)_cxfir;
  NASSERT(state->coef);
  NASSERT(state->delayLine);
  NASSERT(state->delayPos);
  NASSERT_ALIGN(state->coef,8);
  NASSERT_ALIGN(state->delayLine,8);
  NASSERT_ALIGN(state->delayPos,8);
  NASSERT(N%4==0);
  NASSERT_ALIGN(x,8);
  NASSERT((state->M%4)==0);
  NASSERT(x);
  NASSERT(y);
  if(N<=0) return;
  M=state->M;
  pD = ( xtfloat*)state->delayPos;
  pH = (const xtfloat*)state->coef;
  NASSERT(N>0);
  NASSERT(M>0);
  WUR_AE_CBEGIN0((uintptr_t)(state->delayLine));
  WUR_AE_CEND0((uintptr_t)(state->delayLine + M));
  for (n = 0; n<(N>>1); n ++)
  {
    XT_LSIP(x0, pX, 4);
    XT_LSIP(x1, pX, 4);
    XT_LSIP(x2, pX, 4);
    XT_LSIP(x3, pX, 4);
    acc0 = XT_CONST_S(0);
    acc1 = XT_CONST_S(0);
    acc2 = XT_CONST_S(0);
    acc3 = XT_CONST_S(0);
    s0 = x0;
    s1 = x1;
    s2 = x2;
    s3 = x3;
    { xtfloat dummy;  XT_LSXC(dummy, pD, -4); }
    for (m = 0; m<M; m ++)
    {
      xtfloat h0, h1;
      h0 = pH[2*m + 0];
      h1 = pH[2*m + 1];
      XT_MADD_S(acc0, x0, h0); /*re*re*/
      XT_MADD_S(acc1, x0, h1); /*re*im*/

      XT_MADD_S(acc2, x2, h0);
      XT_MADD_S(acc3, x2, h1);

      XT_MSUB_S(acc0, x1, h1); /*-im*im*/
      XT_MADD_S(acc1, x1, h0); /*im*re*/

      XT_MSUB_S(acc2, x3, h1);
      XT_MADD_S(acc3, x3, h0);
      x3 = x1;
      x2 = x0;
      XT_LSXC(x1, pD, -4);
      XT_LSXC(x0, pD, -4);

    }

    { xtfloat dummy;  XT_LSXC(dummy, pD, 4); }
    XT_SSXC(s0, pD, 4);
    XT_SSXC(s1, pD, 4);
    XT_SSXC(s2, pD, 4);
    XT_SSXC(s3, pD, 4);

    XT_SSIP(acc0, pZ, 4);
    XT_SSIP(acc1, pZ, 4);
    XT_SSIP(acc2, pZ, 4);
    XT_SSIP(acc3, pZ, 4);

  }
  state->delayPos = (complex_float*)pD;
} // cxfirf_process()
#endif
