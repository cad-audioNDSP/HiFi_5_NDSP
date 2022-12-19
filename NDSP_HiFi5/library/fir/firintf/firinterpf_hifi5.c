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

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "common_fpu.h"
#include "firinterpf_common.h"

/*-------------------------------------------------------------------------
  Interpolating Block Real FIR Filter
  Computes a real FIR filter (direct-form) with interpolation using IR stored 
  in vector h. The real data input is stored in vector x. The filter output 
  result is stored in vector y. The filter calculates N*D output samples 
  using M*D coefficients and requires last N+M*D-1 samples on the delay line.
  NOTE:
  user application is not responsible for management of delay lines

  Precision: 
  16x16     16-bit data, 16-bit coefficients, 16-bit outputs
  32x16     32-bit data, 16-bit coefficients, 32-bit outputs
  32x32     32-bit data, 32-bit coefficients, 32-bit outputs
  32x32ep   the same as above but using 72-bit accumulator for intermediate 
            computations
  f         floating point

  Input:
  h[M*D]    filter coefficients; h[0] is to be multiplied with the 
            newest sample, Q15, Q31, floating point
  D         interpolation ratio
  N         length of input sample block
  M         length of subfilter. Total length of filter is M*D
  x[N]      input samples, Q15, Q31, floating point
  Output:
  y[N*D]    output samples, Q15, Q31, floating point

  Restrictions:
  x,h,y     should not overlap
  x,h       aligned on an 16-bytes boundary
  N         multiple of 8
  M         multiple of 4
  D         should be >1

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  D - 2, 3 or 4
-------------------------------------------------------------------------*/

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(size_t,firinterpf_alloc,( int D, int M ))
DISCARD_FUN(firinterpf_handle_t,firinterpf_init,( void * objmem,
                                             int D, int M,
                                             const float32_t * restrict h ))
DISCARD_FUN(void,firinterpf_process,( firinterpf_handle_t _firinterp,
                             float32_t * restrict      y,
                       const float32_t * restrict      x, int N ))
#else

/* Instance pointer validation number. */
#define MAGIC     0x283db462

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
    ((size_t)(size)+(align)-1)

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

#define sz_f32    sizeof(float32_t)

/* Interpolator instance structure. */
typedef struct tag_firinterpf_t
{
    uint32_t          magic;     // Instance pointer validation number
    int               D;         // Interpolation factor
    int               M;         // Number of filter coefficients
    proc_fxn_t      * procFxn;   // Filter data processing function
    const float32_t * coef;      // Filter coefficients
    float32_t       * delayLine; // Delay line
    int               delayLen;  // Delay line length, in samples
    int               wrIx;      // Index of the oldest sample
} firinterpf_t, *firinterpf_ptr_t;

/* Calculate the memory block size for an interpolator with given attributes. */
size_t firinterpf_alloc( int D, int M )
{
    int delLen, coefNum;
    NASSERT(D > 1);
    NASSERT(M > 0 && M % 4 == 0);

    if (/*D == 3 ||*/ D == 4) delLen = M + 4;
    else delLen = M + 8;
    coefNum = (M + 4)*D + 1;

    return (ALIGNED_SIZE(sizeof(firinterpf_t), 4)
            + // Delay line
            ALIGNED_SIZE(delLen*sz_f32, 16)
            + // Coefficients
            ALIGNED_SIZE(coefNum*sz_f32, 16));
} /* firinterpf_alloc() */

/* Initialize the interpolator structure. The delay line is zeroed. */
firinterpf_handle_t firinterpf_init( void * objmem, int D, int M,
                                     const float32_t * restrict h )
{
    firinterpf_ptr_t firinterp;
    void           * ptr;
    float32_t      * delLine;
    int              delLen;
    int              coefNum;
    float32_t      * coef;
    float32_t      * coefB;
    proc_fxn_t     * procFxn;
	
    int d, m;

    NASSERT(objmem && h);
    NASSERT_ALIGN(h, 16);
    NASSERT(D > 1);
    NASSERT(M > 0 && M % 4 == 0);

    if (/*D == 3 ||*/ D == 4) delLen = M + 4;
    else delLen = M + 8;
    coefNum = (M + 4)*D + 1;

    procFxn = ( D == 2 ? &fir_interpf_2x :
                D == 3 ? &fir_interpf_3x :
                D == 4 ? &fir_interpf_4x :
                         &fir_interpf_Dx);

    //
    // Partition the memory block.
    //

    ptr         = objmem;
    firinterp   = (firinterpf_ptr_t)ALIGNED_ADDR(ptr, 4);
    ptr         = firinterp + 1;
    delLine     = (float32_t*)ALIGNED_ADDR(ptr, 16);
    ptr         = delLine + delLen;
    coef        = (float32_t*)ALIGNED_ADDR(ptr, 16);
    ptr         = coef + coefNum;

    NASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)firinterpf_alloc(D, M));

    //
    // Break the impulse response into D coefficients banks and copy them in
    // reverted order.
    //

#if HAVE_VFPU
    for (d = 0; d < D; d++)
    {
        coefB = coef + d*(M + 4);

        // To avoid a 1-sample delay, insert a zero coefficient that will match the
        // oldest sample in the delay line. To keep the filter length a multiple of
        // 4, append 3 zeros after the last coefficient.
        coefB[0] = coefB[M + 1] = coefB[M + 2] = coefB[M + 3] = 0.0f;

        // Copy bank's coefficients in reverted order.
        for (m = 0; m < M; m++)
        {
            coefB[m + 1] = h[D*(M - m - 1) + d] * D;
        }
    }
    coef[coefNum - 1] = 0.0f;
#else
    coefB = coef;
    for (m = 0; m < M; m++)
        for (d = 0; d < D; d++)
            coefB[d*M + m] = h[m*D + d] * D;
#endif

    //
    // Zero the delay line.
    //

    for (m = 0; m < delLen; m++)
    {
        delLine[m] = 0;
    }

    //
    // Initialize the interpolator instance.
    //

    firinterp->magic      = MAGIC;
    firinterp->D          = D;
    firinterp->M          = M;
    firinterp->procFxn    = procFxn;
    firinterp->coef       = coef;
    firinterp->delayLine  = delLine;
    firinterp->delayLen   = delLen;
    firinterp->wrIx       = 0;

    return (firinterp);
} /* firinterpf_init() */

/* Put a chunk of input signal into the delay line and compute the filter response. */
void firinterpf_process( firinterpf_handle_t  handle,
                         float32_t * restrict y,
                   const float32_t * restrict x, int N )
{
    firinterpf_ptr_t firinterp = (firinterpf_ptr_t)handle;

    NASSERT(firinterp && firinterp->magic == MAGIC && y && x);
    NASSERT(firinterp->procFxn);
    NASSERT(N % 8 == 0);
    if (N <= 0) return;

    //
    // Call filter's data processing function. It will store the block of input
    // samples to the delay line, and compute the filter response. Returns the
    // updated next position pointer into the delay line buffer.
    //

    firinterp->wrIx = (*firinterp->procFxn)(
        y,
        firinterp->delayLine,
        firinterp->delayLen,
        x,
        firinterp->coef,
        firinterp->wrIx,
        firinterp->D,
        N,
        firinterp->M);
} /* firinterpf_process() */

#endif
