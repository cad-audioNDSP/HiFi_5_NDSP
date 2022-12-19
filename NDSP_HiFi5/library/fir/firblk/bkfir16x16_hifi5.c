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
    Real block FIR filter, 16x16-bit
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
  NOTE: 
  1. User application is not responsible for management of delay lines
  2. User has an option to set IR externally or copy from original location 
     (i.e. from the slower constant memory). In the first case, user is 
     responsible for right alignment, ordering and zero padding of filter 
     coefficients - usually array is composed from zeroes (left padding), 
     reverted IR and right zero padding.


  Precision: 
  16x16    16-bit data, 16-bit coefficients, 16-bit outputs. Ordinary variant 
           and stereo
  32x16    32-bit data, 16-bit coefficients, 32-bit outputs
  32x32    32-bit data, 32-bit coefficients, 32-bit outputs. Ordinary variant 
           and stereo
  32x32ep  the same as above but using 72-bit accumulator for intermediate 
           computations
  f        floating point. Ordinary variant and stereo

  Input:
  x[N*S]   input samples, Q31, Q15, floating point
  h[M]     filter coefficients in normal order, Q31, Q15, floating point
  hl[M]    for stereo filters: filter coefficients for left channel
  hr[M]    for stereo filters: filter coefficients for right channel
  N        length of sample block, should be a multiple of 4
  M        length of filter, should be a multiple of 4
  extIR    if zero, IR is copied from original location, otherwise not
           but user should keep alignment, order of coefficients 
           and zero padding requirements shown below
  S        1 for ordinary (single channel) filters, 2 - for stereo variant
  
  Output:
  y[N*S]   output samples, Q31, Q15, floating point

  Alignment, ordering and zero padding for external IR  (extIR!=0)
  ------------------------+----------+--------------+--------------+----------------
  Function                |Alignment,|Left zero     |   Coefficient| Right zero 
                          | bytes    |padding, bytes|   order      | padding, bytes
  ------------------------+----------+--------------+--------------+----------------
  bkfir16x16_init         |    16    |      2       |  inverted    |  6
  bkfir32x16_init         |    16    |      2       |  inverted    |  6
  bkfir32x32_init         |    16    |      4       |  inverted    |  12
  bkfir32x32ep_init       |    16    |      4       |  inverted    |  12
  bkfirf_init             |    16    |      4       |  inverted    |  12
  stereo_bkfir16x16_init  |    16    |      2       |  inverted    |  6
  stereo_bkfir32x32_init  |    16    |      4       |  inverted    |  12
  stereo_bkfirf_init      |    16    |      4       |  inverted    |  12
  ------------------------+----------+--------------+--------------+----------------

  Restrictions:
  x, y     should not be overlapping
  x, h     aligned on a 16-bytes boundary
  N, M     multiples of 4 
-------------------------------------------------------------------------*/

#define MAGIC     0x5CDED6E2

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
    ((size_t)(size)+(align)-1)

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

#define sz_i16   sizeof(int16_t)

/* Filter instance structure. */
typedef struct tag_bkfir16x16_t
{
    int32_t           magic;     // Instance pointer validation number
    int               M;         // Number of filter coefficients
    const int16_t   * coef;      // Filter coefficients
    int16_t         * delayLine; // Delay line for samples, aligned
    int               delayLen;  // Delay line length, in samples
    int16_t         * delayPos;  // Delay line slot to be filled next
} bkfir16x16_t, *bkfir16x16_ptr_t;

/* Calculate the memory block size for an FIR filter with given attributes. */
size_t bkfir16x16_alloc(int M, int extIR)
{
    NASSERT(M > 0 && M % 4 == 0);
    return (ALIGNED_SIZE(sizeof(bkfir16x16_t), 4)
        + // Delay line
        ALIGNED_SIZE((M + 16)*sz_i16, 16)
        + // Filter coefficients
        (extIR?0:ALIGNED_SIZE((M + 4)*sz_i16, 16)));
} /* bkfir16x16_alloc() */

/* Initialize the filter structure. The delay line is zeroed. */
bkfir16x16_handle_t bkfir16x16_init(void * objmem, int M, int extIR, const int16_t * h)
{
    bkfir16x16_t * bkfir;
    void         * ptr;
    int16_t      * delLine;
    int16_t      * coef;

    int m;

    NASSERT(objmem && h);
    NASSERT_ALIGN(h, 16);
    NASSERT(M > 0 && M % 4 == 0);

    //
    // Partition the memory block
    //

    ptr     = objmem;
    bkfir   = (bkfir16x16_ptr_t)ALIGNED_ADDR(ptr, 4);
    ptr     = bkfir + 1;
    delLine = (int16_t *)ALIGNED_ADDR(ptr, 16);
    ptr     = delLine + M + 16;
    if (extIR)
    {
        coef = (int16_t *)h;
    }
    else
    {
        coef = (int16_t *)ALIGNED_ADDR(ptr, 16);
        ptr  = coef + M + 4;
    }
    ASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)bkfir16x16_alloc(M, extIR));

    //
    // Copy the filter coefficients and zero the delay line. Original impulse
    // response is padded with zeros: three zeros go before the first tap
    // (corresponds to the newest sample), one zero follows the last tap,
    // which matches the oldest sample. After that the order of filter
    // coefficients is reverted.
    //
    if (extIR == 0)
    {
        coef[0] = 0;

        for (m = 1; m<M + 1; m++)
        {
            coef[m] = h[M - m];
        }

        for (; m<M + 4; m++)
        {
            coef[m] = 0;
        }
    }

    for (m = 0; m < M + 16; m++)
    {
        delLine[m] = 0;
    }

    //
    // Initialize the filter instance.
    //

    bkfir->magic     = MAGIC;
    bkfir->M         = M;
    bkfir->coef      = coef;
    bkfir->delayLine = delLine;
    bkfir->delayLen  = M + 16;
    bkfir->delayPos  = delLine;

    return (bkfir);
} /* bkfir16x16_init() */

void bkfir16x16_process( bkfir16x16_handle_t handle,
                         int16_t * restrict  y,
                   const int16_t * restrict  x, int N)
{
    bkfir16x16_ptr_t bkfir = (bkfir16x16_ptr_t)handle;

    const ae_int16x4 * restrict pX;
    const ae_int16x4 * restrict pDr;
          ae_int16x4 * restrict pDw;
    const ae_int16x4 * restrict S0;
    const ae_int16x4 * restrict S1;
    const ae_int16x4 * restrict S2;
    const ae_int16x4 * restrict pH;
          ae_int16x4 * restrict pY;
    ae_valign aY;
    ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
    ae_f64 q8, q9, qa, qb, qc, qd, qe, qf;
    ae_int16x4 d0, d1, d2, d3, d4;
    ae_int16x4 h0;
    ae_int32x2 t0, t1;

    int M;
    int m, n;

    M = bkfir->M;
    NASSERT(bkfir && bkfir->magic == MAGIC && y && x);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN((bkfir->coef), 16);
    NASSERT_ALIGN((bkfir->delayLine), 8);
    NASSERT_ALIGN((bkfir->delayLine + bkfir->delayLen), 8);
    NASSERT_ALIGN((bkfir->delayPos), 8);
    NASSERT(N % 4 == 0 && M % 4 == 0);
    if (N <= 0) return;

    // Setup pointers and circular delay line buffer.
    pX = (const ae_int16x4 *)x;
    pY = (      ae_int16x4 *)y;
    pDw= (      ae_int16x4 *)bkfir->delayPos;
    WUR_AE_CBEGIN0((uintptr_t)(bkfir->delayLine));
    WUR_AE_CEND0  ((uintptr_t)(bkfir->delayLine + bkfir->delayLen));
    aY = AE_ZALIGN64();

    for (n = 0; n < (N >> 4); n++)
    {
        AE_L16X4_IP(d0, pX, 4 * sizeof(int16_t));
        AE_L16X4_IP(d1, pX, 4 * sizeof(int16_t));
        AE_L16X4_IP(d2, pX, 4 * sizeof(int16_t));
        AE_L16X4_IP(d3, pX, 4 * sizeof(int16_t));
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
        AE_L16X4_IP(d0, pX, 4 * sizeof(int16_t));
        AE_L16X4_IP(d1, pX, 4 * sizeof(int16_t));
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
        AE_L16X4_IP(d0, pX, 4 * sizeof(int16_t));
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

    bkfir->delayPos = (int16_t*)pDw;
} /* bkfir16x16_process() */
