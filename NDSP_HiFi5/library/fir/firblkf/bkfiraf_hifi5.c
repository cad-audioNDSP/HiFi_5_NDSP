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
    Real block FIR filter, floating point
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
#include "bkfiraf.h"

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

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(size_t,bkfiraf_alloc,( int M ))
DISCARD_FUN(bkfirf_handle_t,bkfiraf_init,( void * objmem, int M, const float32_t * h ))
#else

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
      ( (size_t)(size) + (align) - 1 )

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
      (void*)( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

#define sz_f32   sizeof(float32_t)

/* Calculate the memory block size for an FIR filter with given attributes. */
size_t bkfiraf_alloc( int M )
{
    int _M;
    _M = M + (-M & 3);
    NASSERT(M > 0);
    return (ALIGNED_SIZE(sizeof(bkfiraf_t), 4)
            + // Delay line
            ALIGNED_SIZE((_M + 8)*sz_f32, 16)
            + // Filter coefficients
            ALIGNED_SIZE((_M + 4)*sz_f32, 16));
} // bkfiraf_alloc()

/* Initialize the filter structure. The delay line is zeroed. */
bkfirf_handle_t bkfiraf_init( void * objmem, int M, const float32_t * h )
{
    bkfiraf_t * bkfir;
    void      * ptr;
    float32_t * delLine;
    float32_t * coef;

    int m, _M, n;

    NASSERT(objmem &&  M > 0 && h);
    _M = M + (-M & 3);

    //
    // Partition the memory block
    //
    ptr     = objmem;
    bkfir   = (bkfiraf_ptr_t)ALIGNED_ADDR(ptr, 4);
    ptr     = bkfir + 1;
    delLine = (float32_t *)ALIGNED_ADDR(ptr, 16);
    ptr     = delLine + _M + 8;
    coef    = (float32_t *)ALIGNED_ADDR(ptr, 16);
    ptr     = coef + _M + 4;
    ASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)bkfiraf_alloc(M));

    //
    // Copy the filter coefficients in reverted order and zero the delay line.
    //

    for (n = 0; n < (-M & 3) + 1; n++)
    {
        coef[n] = 0;
    }
    for (m = 0; m < M; m++, n++)
    {
        coef[n] = h[M - m - 1];
    }

    for (m = 0; m < 3; m++, n++)
    {
        coef[n] = 0;
    }

    for (m = 0; m < _M + 8; m++)
    {
        delLine[m] = 0;
    }

    //
    // Initialize the filter instance.
    //
    bkfir->magic     = BKFIRAF_MAGIC;
    bkfir->M         = _M;
    bkfir->coef      = coef;
    bkfir->delayLine = delLine;
    bkfir->delayLen  = _M + 8;
    bkfir->delayPos  = delLine;

    return (bkfir);
} // bkfiraf_init()
#endif
