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
    Real block FIR filter, 16x16-bit, unaligned data and arbitrary M/N allowed
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

/*-------------------------------------------------------------------------
  Real FIR filter.
  Computes a real FIR filter (direct-form) using IR stored in vector h. The 
  real data input is stored in vector x. The filter output result is stored 
  in vector y. The filter calculates N output samples using M coefficients 
  and requires last M-1 samples in the delay line.
  These functions implement FIR filter with no limitation on size of data
  block, alignment and length of impulse response at the cost of increased
  processing complexity.
  NOTE: 
  User application is not responsible for management of delay lines.

  Precision: 
  16x16    16-bit data, 16-bit coefficients, 16-bit outputs
  32x16    32-bit data, 16-bit coefficients, 32-bit outputs
  32x32    32-bit data, 32-bit coefficients, 32-bit outputs
  32x32ep  the same as above but using 72-bit accumulator for intermediate 
           computations
  f        floating point
  Input:
  x[N]     input samples, Q15, Q31, floating point
  h[M]     filter coefficients in normal order, Q15, Q31, floating point
  N        length of sample block
  M        length of filter
  Output:
  y[N]     input samples, Q15, Q31, floating point 

  Restrictions:
  x,y      should not be overlapping
-------------------------------------------------------------------------*/

/* Instance pointer validation number. */
#define MAGIC     0x29EDBF10

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
    ((size_t)(size)+(align)-1)

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

#define sz_i16   sizeof(int16_t)

/* Filter instance structure. */
typedef struct tag_bkfira16x16_t
{
    int32_t           magic;     // Instance pointer validation number
    int               M;         // Number of filter coefficients
    const int16_t   * coef;      // Filter coefficients
    int16_t         * delayLine; // Delay line for samples, aligned
    int               delayLen;  // Delay line length, in samples
    int               wrIx;      // Index of the oldest sample
} bkfira16x16_t, *bkfira16x16_ptr_t;

/* Calculate the memory block size for an FIR filter with given attributes. */
size_t bkfira16x16_alloc( int M )
{
    int _M;
    _M = M + (-M & 3);
    NASSERT(M > 0);

    return (ALIGNED_SIZE(sizeof(bkfira16x16_t), 4)
        + // Delay line
        ALIGNED_SIZE((_M + 16)*sz_i16, 16)
        + // Filter coefficients
        ALIGNED_SIZE((_M + 4)*sz_i16, 16));
} /* bkfira16x16_alloc() */

/* Initialize the filter structure. The delay line is zeroed. */
bkfira16x16_handle_t bkfira16x16_init(void * objmem, int M, const int16_t * h)
{
    bkfira16x16_t * bkfir;
    void          * ptr;
    int16_t       * delLine;
    int16_t       * coef;

    int m, _M, n;

    NASSERT(objmem &&  M > 0 && h);
    _M = M + (-M & 3);

    //
    // Partition the memory block
    //
    ptr     = objmem;
    bkfir   = (bkfira16x16_ptr_t)ALIGNED_ADDR(ptr, 4);
    ptr     = bkfir + 1;
    delLine = (int16_t *)ALIGNED_ADDR(ptr, 16);
    ptr     = delLine + _M + 16;
    coef    = (int16_t *)ALIGNED_ADDR(ptr, 16);
    ptr     = coef + _M + 4;

    ASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)bkfira16x16_alloc(M));

    //
    // Copy the filter coefficients in reverted order and zero the delay line.
    //

    for (n = 0; n<(-M & 3) + 1; n++)
    {
        coef[n] = 0;
    }
    for (m = 0; m<M; m++, n++)
    {
        coef[n] = h[M - m - 1];
    }

    for (m = 0; m<3; m++, n++)
    {
        coef[n] = 0;
    }

    for (m = 0; m < _M + 16; m++)
    {
        delLine[m] = 0;
    }

    //
    // Initialize the filter instance.
    //
    bkfir->magic     = MAGIC;
    bkfir->M         = _M;
    bkfir->coef      = coef;
    bkfir->delayLine = delLine;
    bkfir->delayLen  = _M + 16;
    bkfir->wrIx      = 0;

    return (bkfir);
} /* bkfira16x16_init() */

/* process block of samples */
void bkfira16x16_process( bkfira16x16_handle_t handle,
                         int16_t * restrict  y,
                   const int16_t * restrict  x, int N )
{
    bkfira16x16_ptr_t bkfir = (bkfira16x16_ptr_t)handle;

    const ae_int16x4 * restrict pX;
    const ae_int16x4 * restrict pDr;
          ae_int16x4 * restrict pDw;
    const ae_int16x4 * restrict S0;
    const ae_int16x4 * restrict S1;
    const ae_int16x4 * restrict S2;
    const ae_int16x4 * restrict pH;
          ae_int16x4 * restrict pY;
    ae_valign aX, aD, aY;
    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
    ae_int16x4 d0, d1, d2, d3, d4;
    ae_int16x4 h0;
    ae_int32x2 t0, t1;

    int M;
    int m, n;
    int wrIx, ua_cnt;

    M    = bkfir->M;
    wrIx = bkfir->wrIx;

    NASSERT(bkfir && bkfir->magic == MAGIC && y && x);
    NASSERT_ALIGN((bkfir->coef), 16);
    NASSERT_ALIGN((bkfir->delayLine), 8);
    NASSERT_ALIGN((bkfir->delayLine + bkfir->delayLen), 8);
    if (N <= 0) return;

    //
    // Setup pointers and circular delay line buffer.
    //
    pX = (const ae_int16x4 *)x;
    pY = (      ae_int16x4 *)y;
    pDw= (      ae_int16x4 *)(bkfir->delayLine + wrIx);
    WUR_AE_CBEGIN0((uintptr_t)(bkfir->delayLine));
    WUR_AE_CEND0  ((uintptr_t)(bkfir->delayLine + bkfir->delayLen));

    ua_cnt = (4 - wrIx) & 3;
    ua_cnt = XT_MIN(ua_cnt, N);
    //process first 0..3 samples (until pDw is aligned to 8 bytes)
    if (ua_cnt)
    {
        // Insert 0..3 input samples into the delay line one-by-one.
        __Pragma("loop_count min=1, max=3");
        for (n = 0; n < ua_cnt; n++)
        {
            AE_L16_IP(d0, castxcc(ae_int16, pX), sizeof(int16_t));
            AE_S16_0_XC(d0, castxcc(ae_int16, pDw), sizeof(int16_t));
        }

        // Reset the coefficients pointer. Now it looks at the tap corresponding
        // to the oldest sample in the delay line.
        pH = (const ae_int16x4 *)bkfir->coef;

        // Circular buffer write pointer looks at the oldest sample: M+4 samples
        // back from the newest one.
        pDr = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), (16 - ua_cnt) * sizeof(int16_t));
        AE_LA16X4POS_PC(aD, pDr);
        AE_LA16X4_IC(d0, aD, pDr);
        q0 = q1 = q2 = q3 = AE_ZERO64();

        // Inner loop: process 4 taps for 4 accumulators on each trip. Actually we
        // perform M+4 MACs for each accumulator, 4 of which fall on zero taps
        // inserted into the impulse response during initialization.
        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            AE_L16X4_IP(h0, pH, 4 * sizeof(int16_t));
            AE_LA16X4_IC(d1, aD, pDr);
            AE_MULAFQ16X2_FIR_3(q0, q1, d0, d1, h0);
            AE_MULAFQ16X2_FIR_1(q2, q3, d0, d1, h0);
            d0 = d1;
        }
        t0 = AE_TRUNCA32X2F64S(q1, q0, 32);
        t1 = AE_SAT32X2(q2, q2);
        d0 = AE_ROUND16X4F32SASYM(t1, t0);
        AE_S16_0_IP(d0, castxcc(ae_int16, pY), sizeof(int16_t));
        if (ua_cnt > 1) { d0 = AE_SEL16_4321(d0, d0); AE_S16_0_IP(d0, castxcc(ae_int16, pY), sizeof(int16_t)); }
        if (ua_cnt > 2) { d0 = AE_SEL16_4321(d0, d0); AE_S16_0_IP(d0, castxcc(ae_int16, pY), sizeof(int16_t)); }
        N -= ua_cnt;
    }

    aX = AE_LA64_PP(pX);
    aY = AE_ZALIGN64();
    for (n = 0; n < (N >> 4); n++)
    {
        AE_LA16X4_IP(d0, aX, pX);
        AE_LA16X4_IP(d1, aX, pX);
        AE_LA16X4_IP(d2, aX, pX);
        AE_LA16X4_IP(d3, aX, pX);
        AE_S16X4_XC(d0, pDw, 4 * sizeof(int16_t));
        AE_S16X4_XC(d1, pDw, 4 * sizeof(int16_t));
        AE_S16X4_XC(d2, pDw, 4 * sizeof(int16_t));
        AE_S16X4_XC(d3, pDw, 4 * sizeof(int16_t));

        pH = (const ae_int16x4 *)bkfir->coef;
        AE_L16X4_IP(h0, pH, 4 * sizeof(int16_t));

        pDr = pDw;
        AE_L16X4_XC(d0, pDr, 4 * sizeof(int16_t));
        S0 = pDr;
        AE_L16X4_XC(d1, pDr, 4 * sizeof(int16_t));
        AE_L16X4_XC(d2, pDr, 4 * sizeof(int16_t));
        S1 = pDr;
        AE_L16X4_XC(d3, pDr, 4 * sizeof(int16_t));
        AE_L16X4_XC(d4, pDr, 4 * sizeof(int16_t));
        S2 = pDr;

        AE_MULFQ16X2_FIR_3(q0, q1, d0, d1, h0);
        AE_MULFQ16X2_FIR_1(q2, q3, d0, d1, h0);
        AE_MULFQ16X2_FIR_3(q4, q5, d1, d2, h0);
        AE_MULFQ16X2_FIR_1(q6, q7, d1, d2, h0);
        AE_MULFQ16X2_FIR_3(q8, q9, d2, d3, h0);
        AE_MULFQ16X2_FIR_1(qa, qb, d2, d3, h0);
        AE_MULFQ16X2_FIR_3(qc, qd, d3, d4, h0);
        AE_MULFQ16X2_FIR_1(qe, qf, d3, d4, h0);

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            AE_L16X4_IP(h0, pH, 4 * sizeof(int16_t));

            AE_L16X4_XC(d0, S0, 4 * sizeof(int16_t));
            d1 = AE_L16X4_I(S0, 0);
            AE_L16X4_XC(d2, S1, 4 * sizeof(int16_t));
            d3 = AE_L16X4_I(S1, 0);
            AE_L16X4_XC(d4, S2, 4 * sizeof(int16_t));

            AE_MULAFQ16X2_FIR_3(q0, q1, d0, d1, h0);
            AE_MULAFQ16X2_FIR_1(q2, q3, d0, d1, h0);
            AE_MULAFQ16X2_FIR_3(q4, q5, d1, d2, h0);
            AE_MULAFQ16X2_FIR_1(q6, q7, d1, d2, h0);
            AE_MULAFQ16X2_FIR_3(q8, q9, d2, d3, h0);
            AE_MULAFQ16X2_FIR_1(qa, qb, d2, d3, h0);
            AE_MULAFQ16X2_FIR_3(qc, qd, d3, d4, h0);
            AE_MULAFQ16X2_FIR_1(qe, qf, d3, d4, h0);
        }

        t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
        t1 = AE_SAT32X2(q2, q3);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
        t0 = AE_TRUNCA32X2F64S(q4, q5, 32);
        t1 = AE_SAT32X2(q6, q7);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
        t0 = AE_TRUNCA32X2F64S(q8, q9, 32);
        t1 = AE_SAT32X2(qa, qb);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
        t0 = AE_TRUNCA32X2F64S(qc, qd, 32);
        t1 = AE_SAT32X2(qe, qf);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
    }
    if (N & 8)
    {
        AE_LA16X4_IP(d0, aX, pX);
        AE_LA16X4_IP(d1, aX, pX);
        AE_S16X4_XC(d0, pDw, 4 * sizeof(int16_t));
        AE_S16X4_XC(d1, pDw, 4 * sizeof(int16_t));

        pH = (const ae_int16x4 *)bkfir->coef;

        pDr = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 8 * sizeof(int16_t));
        S0 = pDr;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 4 * sizeof(int16_t));
        S1 = pDr;
        q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            AE_L16X4_IP(h0, pH, 4 * sizeof(int16_t));

            AE_L16X4_XC(d0, S0, 4 * sizeof(int16_t));
            AE_L16X4_XC(d1, S1, 4 * sizeof(int16_t));
            d2 = AE_L16X4_I(S1, 0);

            AE_MULAFQ16X2_FIR_3(q0, q1, d0, d1, h0);
            AE_MULAFQ16X2_FIR_1(q2, q3, d0, d1, h0);
            AE_MULAFQ16X2_FIR_3(q4, q5, d1, d2, h0);
            AE_MULAFQ16X2_FIR_1(q6, q7, d1, d2, h0);
        }

        t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
        t1 = AE_SAT32X2(q2, q3);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
        t0 = AE_TRUNCA32X2F64S(q4, q5, 32);
        t1 = AE_SAT32X2(q6, q7);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
    }
    if (N & 4)
    {
        // Load 4 input samples.
        AE_LA16X4_IP(d0, aX, pX);
        // Store 4 samples to the delay line buffer with circular address update.
        AE_S16X4_XC(d0, pDw, 4 * sizeof(int16_t));

        // Reset the coefficients pointer. Now it looks at the tap corresponding
        // to the oldest sample in the delay line.
        pH = (const ae_int16x4 *)bkfir->coef;

        // Circular buffer write pointer looks at the oldest sample: M+4 samples
        // back from the newest one.
        pDr = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), 12 * sizeof(int16_t));
        AE_L16X4_XC(d0, pDr, 4 * sizeof(int16_t));
        q0 = q1 = q2 = q3 = AE_ZERO64();

        // Inner loop: process 4 taps for 4 accumulators on each trip. Actually we
        // perform M+4 MACs for each accumulator, 4 of which fall on zero taps
        // inserted into the impulse response during initialization.
        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            AE_L16X4_IP(h0, pH, 4 * sizeof(int16_t));
            AE_L16X4_XC(d1, pDr, 4 * sizeof(int16_t));
            AE_MULAFQ16X2_FIR_3(q0, q1, d0, d1, h0);
            AE_MULAFQ16X2_FIR_1(q2, q3, d0, d1, h0);
            d0 = d1;
        }
        t0 = AE_TRUNCA32X2F64S(q0, q1, 32);
        t1 = AE_SAT32X2(q2, q3);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
    }
    AE_SA64POS_FP(aY, pY);

    N &= 3;
    if (N)
    {
        // Insert 0..3 input samples into the delay line one-by-one.
        __Pragma("loop_count min=1, max=3");
        for (n = 0; n < N; n++)
        {
            AE_L16_IP(d0, castxcc(ae_int16, pX), sizeof(int16_t));
            AE_S16_0_XC(d0, castxcc(ae_int16, pDw), sizeof(int16_t));
        }

        // Reset the coefficients pointer. Now it looks at the tap corresponding
        // to the oldest sample in the delay line.
        pH = (const ae_int16x4 *)bkfir->coef;

        // Circular buffer write pointer looks at the oldest sample: M+4 samples
        // back from the newest one.
        pDr = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pDr), (16 - N) * sizeof(int16_t));
        AE_L16X4_XC(d0, pDr, 4 * sizeof(int16_t));
        q0 = q1 = q2 = q3 = AE_ZERO64();

        // Inner loop: process 4 taps for 4 accumulators on each trip. Actually we
        // perform M+4 MACs for each accumulator, 4 of which fall on zero taps
        // inserted into the impulse response during initialization.
        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            AE_L16X4_IP(h0, pH, 4 * sizeof(int16_t));
            AE_L16X4_XC(d1, pDr, 4 * sizeof(int16_t));
            AE_MULAFQ16X2_FIR_3(q0, q1, d0, d1, h0);
            AE_MULAFQ16X2_FIR_1(q2, q3, d0, d1, h0);
            d0 = d1;
        }
        t0 = AE_TRUNCA32X2F64S(q1, q0, 32);
        t1 = AE_SAT32X2(q2, q2);
        d0 = AE_ROUND16X4F32SASYM(t1, t0);
        AE_S16_0_IP(d0, castxcc(ae_int16, pY), sizeof(int16_t));
        if (N > 1) { d0 = AE_SEL16_4321(d0, d0); AE_S16_0_IP(d0, castxcc(ae_int16, pY), sizeof(int16_t)); }
        if (N > 2) { d0 = AE_SEL16_4321(d0, d0); AE_S16_0_IP(d0, castxcc(ae_int16, pY), sizeof(int16_t)); }
    }

    bkfir->wrIx = (int)((int16_t *)pDw - bkfir->delayLine);
} /* bkfira16x16_process() */
