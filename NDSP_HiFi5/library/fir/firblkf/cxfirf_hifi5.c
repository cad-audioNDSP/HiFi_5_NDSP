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

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(size_t,cxfirf_alloc,( int M , int extIR))
DISCARD_FUN(cxfirf_handle_t,cxfirf_init,( void *             objmem,
                                     int                M, int extIR,
                               const complex_float * restrict h ))
#else

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
      ( (size_t)(size) + (align) - 1 )

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
      (void*)( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

#define sz_f32c     sizeof(complex_float)

/* Calculate the memory block size for an FIR filter with given attributes. */
size_t cxfirf_alloc( int M , int extIR)
{
    NASSERT(M > 0 && M % 4 == 0);
    return (ALIGNED_SIZE(sizeof(cxfirf_t), 4)
            + // Delay line
           // ALIGNED_SIZE((M + 4)*sz_f32c, 16)
            ALIGNED_SIZE((M + 8)*sz_f32c, 16)
            + // Filter coefficients
            (extIR ? 0 : ALIGNED_SIZE(M*sz_f32c, 16)));
} // cxfirf_alloc()

/* Initialize the filter structure. The delay line is zeroed. */
cxfirf_handle_t cxfirf_init( void *             objmem,
                                     int                M, int extIR,
                               const complex_float * restrict h )
{
    cxfirf_t      * cxfir;
    void          * ptr;
    complex_float * delLine;
    complex_float * coef;
    union { struct { float re, im; }s; complex_float c; } zero = { { 0.f, 0.f } };
    int m;

    NASSERT(objmem && h);
    NASSERT_ALIGN(h, 16);
    NASSERT(M > 0 && M % 4 == 0);

    //
    // Partition the memory block.
    //
    ptr     = objmem;
    cxfir   = (cxfirf_ptr_t)ALIGNED_ADDR(ptr, 4);
    ptr     = cxfir + 1;
    delLine = (complex_float *)ALIGNED_ADDR(ptr, 16);
    //ptr     = delLine + (M + 4);
    ptr = delLine + (M + 8);
    if (extIR)
    {
        coef = (complex_float *)h;
    }
    else
    {
        coef = (complex_float *)ALIGNED_ADDR(ptr, 16);
        ptr  = coef + M;
    }
    ASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)cxfirf_alloc(M, extIR));

    //
    // Copy the filter coefficients and zero the delay line. Original impulse
    // response is padded with zeros: one zero follows the last tap,
    // which matches the oldest sample. After that the order of filter
    // coefficients is reverted.
    //
    if (extIR == 0)
    {
        for (m = 0; m < M; m++)
        {
            coef[m] = h[m];
        }
    }
    //for (m = 0; m < M + 4; m++)
    for (m = 0; m < M + 8; m++)
    {
        delLine[m] = zero.c;
    }

    //
    // Initialize the filter instance.
    //
    cxfir->magic     = CXFIRF_MAGIC;
    cxfir->M         = M;
    cxfir->coef      = coef;
    cxfir->delayLine = delLine;
   // cxfir->delayLen  = M + 4;
    cxfir->delayLen  = M + 8;
    cxfir->delayPos  = delLine;

    return (cxfir);
} // cxfirf_init()
#endif
