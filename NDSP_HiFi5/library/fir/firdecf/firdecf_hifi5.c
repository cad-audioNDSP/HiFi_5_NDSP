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

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "common_fpu.h"
#include "firdecf_common.h"

/*-------------------------------------------------------------------------
  Decimating Block Real FIR Filter
  Computes a real FIR filter (direct-form) with decimation using IR stored 
  in vector h. The real data input is stored in vector x. The filter output 
  result is stored in vector r. The filter calculates N output samples using
  M coefficients and requires last D*N+M-1 samples on the delay line.
  NOTE:
  - To avoid aliasing IR should be synthesized in such a way to be narrower 
    than input sample rate divided to 2D.
  - user application is not responsible for management of delay lines

  Precision: 
  16x16     16-bit data, 16-bit coefficients, 16-bit outputs
  32x16     32-bit data, 16-bit coefficients, 32-bit outputs
  32x32     32-bit data, 32-bit coefficients, 32-bit outputs
  32x32ep   the same as above but using 72-bit accumulator for intermediate 
            computations
  f         floating point

  Input:
  h[M]      filter coefficients; h[0] is to be multiplied with the newest 
            sample, Q15, Q31, floating point
  D         decimation factor 
  N         length of output sample block
  M         length of filter
  x[D*N]    input samples, Q15, Q31, floating point
  Output:
  y[N]      output samples, Q15, Q31, floating point

  Restriction:
  x,h,r     should not overlap
  x,h       aligned on an 16-bytes boundary
  N         multiple of 8
  D         should exceed 1

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  D - 2, 3 or 4
-------------------------------------------------------------------------*/

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(size_t,firdecf_alloc,( int D, int M ))
DISCARD_FUN(firdecf_handle_t,firdecf_init,( void * objmem, int D,
                               int M, const float32_t * restrict h ))
DISCARD_FUN(void,firdecf_process,( firdecf_handle_t _firdec,
                          float32_t * restrict       y,
                    const float32_t *                x, int N ))
#else

/* Instance pointer validation number. */
#define MAGIC     0x38fa74d2

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
    ((size_t)(size)+(align)-1)

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
    (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

#define sz_f32    sizeof(float32_t)

/* Decimator instance structure. */
typedef struct tag_firdecf_t
{
    uint32_t          magic;     // Instance pointer validation number
    int               D;         // Decimation factor
    int               M;         // Number of filter coefficients
    const float32_t * coef;      // Filter coefficients
    proc_fxn_t      * procFxn;   // Filter data processing function
    float32_t       * delayLine; // Delay line
    int               delayLen;  // Delay line length, in samples
    int               wrIx;      // Index of the oldest sample
} firdecf_t, *firdecf_ptr_t;

/* Calculate the memory block size for a decimator with given attributes. */
size_t firdecf_alloc( int D, int M )
{
    int bM;
    int delayLen, coefNum;
    NASSERT(D > 1 && M > 0);

    bM = (M + 3)&~3;
    coefNum = bM + 9;
    delayLen = bM + 8 * D;

    return (ALIGNED_SIZE(sizeof(firdecf_t), 4)
            + // Delay line
            ALIGNED_SIZE(delayLen*sz_f32, 16)
            + // Coefficients
            ALIGNED_SIZE(coefNum*sz_f32, 16));
} /* firdecf_alloc() */

/* Initialize the decimator structure. The delay line is zeroed. */
firdecf_handle_t firdecf_init( void * objmem, int D, int M,
                               const float32_t * restrict h)
{
    firdecf_ptr_t firdec;
    void *        ptr;
    float32_t *   delLine;
    int           delLen;
    int           coefNum;
    float32_t *   coef;
    float32_t *   coefB;
    proc_fxn_t *  procFxn;

    int bM;
    int m;

    NASSERT(objmem && h);
    NASSERT_ALIGN(h, 16);
    NASSERT(D > 1 && M > 0);

    //
    // Select the processing function, delay line length and coefficients
    // block layout.
    //

    bM = (M + 3)&~3;
    coefNum = bM + 9;
    delLen = bM + 8 * D;
    procFxn = (D == 2 ? &fir_decimaf_2x :
               D == 3 ? &fir_decimaf_3x :
               D == 4 ? &fir_decimaf_4x :
                        &fir_decimaf_Dx);

    //
    // Partition the memory block.
    //

    ptr     = objmem;
    firdec  = (firdecf_ptr_t)(ALIGNED_ADDR(ptr, 4));
    ptr     = firdec + 1;
    delLine = (float32_t*)(ALIGNED_ADDR(ptr, 16));
    ptr     = delLine + delLen;
    coef    = (float32_t*)(ALIGNED_ADDR(ptr, 16));
    ptr     = coef + coefNum;

    NASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)firdecf_alloc(D, M));

    //
    // Copy the filter coefficients in reverted order and zero the delay line.
    //

    coefB = coef;

    // Begin by a few zeros to make the number of coefficients a multiple of 4.
    for (m = M; m < bM; m++)
    {
        *coefB++ = 0;
    }

    // Pad the coefficients with 9 zeros (one at the beginning, 8 following
    // the end) to avoid a 1-sample delay of filter response.
    coefB[0] = 0;
    for (m = 1; m < 9; m++)
    {
        coefB[M + m] = 0;
    }

    // Copy coefficients in reverted order.
    for (m = 1; m <= M; m++)
    {
        coefB[m] = h[M - m];
    }

    //
    // Zero the delay line.
    //

    for (m = 0; m < delLen; m++)
    {
        delLine[m] = 0;
    }

    //
    // Initialize the decimator instance.
    //

    firdec->magic     = MAGIC;
    firdec->D         = D;
    firdec->M         = bM;
    firdec->procFxn   = procFxn;
    firdec->coef      = coef;
    firdec->delayLine = delLine;
    firdec->delayLen  = delLen;
    firdec->wrIx      = 0;

    return (firdec);
} /* firdecf_init() */

/* process block of samples */
void firdecf_process( firdecf_handle_t     handle,
                      float32_t * restrict y,
                const float32_t *          x, int N )
{
    firdecf_ptr_t firdec = (firdecf_ptr_t)handle;

    NASSERT(firdec->magic == MAGIC);
    NASSERT(firdec->procFxn);
    NASSERT(N % 8 == 0);
    if (N <= 0) return;

    //
    // Call filter's data processing function. It will store the block of input
    // samples to the delay line, and compute the filter response. Returns the
    // updated next position pointer into the delay line buffer.
    //

    firdec->wrIx = (*firdec->procFxn)(
        y,
        firdec->delayLine,
        firdec->delayLen,
        x,
        firdec->coef,
        firdec->wrIx,
        firdec->D,
        N,
        firdec->M);
} /* firdecf_process() */
#endif