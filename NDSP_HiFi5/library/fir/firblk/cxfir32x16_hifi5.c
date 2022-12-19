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
    Complex block FIR filter, 32x16-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "cxfir32x16_common.h"

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

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
      ( (size_t)(size) + (align) - 1 )

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
      (void*)( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

#define sz_i16c     sizeof(complex_fract16)
#define sz_i32c     sizeof(complex_fract32)

/* Allocation routine for filters. Returns: size of memory in bytes to be allocated */
size_t cxfir32x16_alloc( int M, int extIR )
{
    NASSERT(M > 0 && M % 4 == 0);
    return (ALIGNED_SIZE(sizeof(cxfir32x16_t), 4)
            + // Delay line
            ALIGNED_SIZE((M + 4)*sz_i32c, 16)
            + // Filter coefficients
            (extIR ? 0 : ALIGNED_SIZE(2 * (M + 4)*sz_i16c, 16)));
} // cxfir32x16_alloc()

/* Initialization for filters. Returns: handle to the object */
cxfir32x16_handle_t cxfir32x16_init( void *             objmem,
                                     int                M,
                                     int                extIR,
                               const complex_fract16 * restrict h )
{
    cxfir32x16_ptr_t  cxfir;
    void            * ptr;
    complex_fract16 * coef;
    complex_fract16 * coefB;
    complex_fract32 * delLine;

    int m;

    NASSERT(objmem && h);
    NASSERT_ALIGN(h, 16);
    NASSERT(M > 0 && M % 4 == 0);

    //
    // Partition the memory block.
    //
    ptr     = objmem;
    cxfir   = (cxfir32x16_ptr_t)ALIGNED_ADDR(ptr, 4);
    ptr     = cxfir + 1;
    delLine = (complex_fract32*)ALIGNED_ADDR(ptr, 16);
    ptr     = delLine + (M + 4);
    if (extIR)
    {
        coef = (complex_fract16*)h;
    }
    else
    {
        coef = (complex_fract16*)ALIGNED_ADDR(ptr, 16);
        ptr = coef + 2 * (M + 4);
    }
    ASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)cxfir32x16_alloc(M, extIR));

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
    cxfir->magic     = CXFIR32X16_MAGIC;
    cxfir->M         = M;
    cxfir->coef      = coef;
    cxfir->delayLine = delLine;
    cxfir->delayLen  = M + 4;
    cxfir->delayPos  = delLine;

    return (cxfir);
} // cxfir32x16_init()

