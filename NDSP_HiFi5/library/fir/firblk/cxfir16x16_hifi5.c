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
    Complex block FIR filter, 16x16-bit
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

/* Instance pointer validation number. */
#define MAGIC     0x6B2D04F2

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
    ((size_t)(size)+(align)-1)

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

#define sz_i16c     sizeof(complex_fract16)

/* Filter instance structure. */
typedef struct tag_cxfir16x16_t
{
    uint32_t                magic;     // Instance pointer validation number
    int                     M;         // Number of complex filter coefficients
    const complex_fract16 * coef;      // M complex coefficients
    complex_fract16       * delayLine; // Delay line for complex samples
    int                     delayLen;  // Delay line length, in complex samples
    complex_fract16       * delayPos;  // Delay line slot to be filled next
} cxfir16x16_t, *cxfir16x16_ptr_t;

/* Calculate the memory block size for an FIR filter with given attributes. */
size_t cxfir16x16_alloc( int M, int extIR )
{
    NASSERT(M > 0 && M % 4 == 0);
    return (ALIGNED_SIZE(sizeof(cxfir16x16_t), 4)
            + // Delay line
            ALIGNED_SIZE((M + 4)*sz_i16c, 16)
            + // Filter coefficients
            (extIR ? 0 : ALIGNED_SIZE(2 * (M + 4)*sz_i16c, 16)));
} /* cxfir16x16_alloc() */

/* Initialize the filter structure. The delay line is zeroed. */
cxfir16x16_handle_t cxfir16x16_init( void *             objmem,
                                     int                M,
                                     int extIR,
                                const complex_fract16 * restrict h)
{
    cxfir16x16_ptr_t  cxfir;
    void            * ptr;
    complex_fract16 * coef;
    complex_fract16 * coefB;
    complex_fract16 * delLine;

    int m;

    NASSERT(objmem && h);
    NASSERT_ALIGN(h, 16);
    NASSERT(M > 0 && M % 4 == 0);

    //
    // Partition the memory block.
    //
    ptr     = objmem;
    cxfir   = (cxfir16x16_ptr_t)ALIGNED_ADDR(ptr, 4);
    ptr     = cxfir + 1;
    delLine = (complex_fract16 *)ALIGNED_ADDR(ptr, 16);
    ptr     = delLine + (M + 4);
    if (extIR)
    {
        coef = (complex_fract16 *)h;
    }
    else
    {
        coef = (complex_fract16 *)ALIGNED_ADDR(ptr, 16);
        ptr = coef + 2 * (M + 4);
    }
    ASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)cxfir16x16_alloc(M, extIR));

    //
    // Copy the filter coefficients in reverted order and zero the delay line.
    //
    if (extIR == 0)
    {
        coefB = coef + (M + 4);

        coef[0].a = coefB[0].a = 0;
        for (m = 1; m < M + 1; m++)
        {
            coef[m].s.re = h[(M - m)].s.re;
            coef[m].s.im = AE_NEG16S(h[(M - m)].s.im);

            coefB[m].s.re = h[(M - m)].s.im;
            coefB[m].s.im = h[(M - m)].s.re;
        }
        for (; m < M + 4; m++)
        {
            coef[m].a = 0;
            coefB[m].a = 0;
        }
    }

    for (m = 0; m < M + 4; m++)
    {
        delLine[m].a = 0;
    }

    //
    // Initialize the filter instance.
    //
    cxfir->magic     = MAGIC;
    cxfir->M         = M;
    cxfir->coef      = coef;
    cxfir->delayLine = delLine;
    cxfir->delayLen  = M + 4;
    cxfir->delayPos  = delLine;

    return (cxfir);
} /* cxfir16x16_init() */

/* Put a chunk of input signal into the delay line and compute the filter
 * response. */
void cxfir16x16_process( cxfir16x16_handle_t handle,
                         complex_fract16 * restrict  y,
                   const complex_fract16 * restrict  x, int N)
{
    cxfir16x16_ptr_t cxfir = (cxfir16x16_ptr_t)handle;

    const ae_int16x8 *          pX;
          ae_int16x8 * restrict pDw;
    const ae_int16x8 *          S0;
    const ae_int16x4 *          S1;
    const ae_int16x8 *          S2;
    const ae_int16x8 *          pH;
          ae_int16x4 * restrict pY;
    ae_valign aS1, aY;
    ae_valignx2 aS2;
    ae_f64 q0r, q1r, q2r, q3r;
    ae_f64 q0i, q1i, q2i, q3i;
    ae_int16x4 d0, d1, d2, d3, d4, d5;
    ae_int16x4 c0, c0i, c1, c1i;
    ae_int32x2 t0, t1;

    int M;
    int n, m;

    M = cxfir->M;
    NASSERT(cxfir && cxfir->magic == MAGIC && y && x);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN((cxfir->coef), 16);
    NASSERT_ALIGN((cxfir->delayLine), 16);
    NASSERT_ALIGN((cxfir->delayLine + cxfir->delayLen), 16);
    NASSERT_ALIGN((cxfir->delayPos), 16);
    NASSERT(N % 4 == 0 && M % 4 == 0);
    if (N <= 0) return;

    // Setup pointers and circular delay line buffer.
    pX = (const ae_int16x8 *)x;
    pY = (      ae_int16x4 *)y;
    pDw= (      ae_int16x8 *)(cxfir->delayPos);
    WUR_AE_CBEGIN0((uintptr_t)(cxfir->delayLine));
    WUR_AE_CEND0  ((uintptr_t)(cxfir->delayLine + cxfir->delayLen));
    aY = AE_ZALIGN64();

    //
    // Break the input signal into 4-samples blocks. For each block, store 4
    // samples to the delay line and compute the filter response.
    //
    __Pragma("loop_count min=1");
    for (n = 0; n < (N >> 2); n++)
    {
        AE_L16X4X2_IP(d0, d1, pX, 4 * sz_i16c);
        AE_S16X4X2_XC(d0, d1, pDw, 4 * sz_i16c);

        pH = (const ae_int16x8 *)cxfir->coef;
        S0 = pDw;
        S1 = (const ae_int16x4 *)((complex_fract16 *)pDw + 1);
        S2 = (const ae_int16x8 *)((complex_fract16 *)pDw + 3);

        AE_LA16X4POS_PC(aS1, S1);
        AE_LA16X4_IC(d1, aS1, S1);
        AE_LA16X4X2POS_PC(aS2, S2);

        q0r = q1r = q2r = q3r = AE_ZERO64();
        q0i = q1i = q2i = q3i = AE_ZERO64();

        //__Pragma("loop_count min=2");
        __Pragma("no_unroll");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            AE_L16X4X2_XC(d0, d2, S0, 4 * sz_i16c);
            d4 = AE_L16X4_I((ae_int16x4*)S0, 0);
            AE_LA16X4X2_IC(d3, d5, aS2, S2);

            AE_L16X4X2_X(c0i, c1i, pH, (M + 4)*sz_i16c);
            AE_L16X4X2_IP(c0, c1, pH, 4 * sz_i16c);

            AE_MULAAAA2Q16(q0r, q0i, d0, d0, c0, c0i);
            AE_MULAAAA2Q16(q1r, q1i, d1, d1, c0, c0i);
            AE_MULAAAA2Q16(q2r, q2i, d2, d2, c0, c0i);
            AE_MULAAAA2Q16(q3r, q3i, d3, d3, c0, c0i);

            AE_MULAAAA2Q16(q0r, q0i, d2, d2, c1, c1i);
            AE_MULAAAA2Q16(q1r, q1i, d3, d3, c1, c1i);
            AE_MULAAAA2Q16(q2r, q2i, d4, d4, c1, c1i);
            AE_MULAAAA2Q16(q3r, q3i, d5, d5, c1, c1i);

            d1 = d5;
        }

        t0 = AE_TRUNCA32X2F64S(q0r, q0i, 33);
        t1 = AE_TRUNCA32X2F64S(q1r, q1i, 33);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
        t0 = AE_TRUNCA32X2F64S(q2r, q2i, 33);
        t1 = AE_TRUNCA32X2F64S(q3r, q3i, 33);
        AE_SA16X4_IP(AE_ROUND16X4F32SASYM(t0, t1), aY, pY);
    }
    AE_SA64POS_FP(aY, pY);

    cxfir->delayPos = (complex_fract16*)pDw;
} /* cxfir16x16_process */
