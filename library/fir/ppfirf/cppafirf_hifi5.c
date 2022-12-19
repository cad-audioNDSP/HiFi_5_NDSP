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
    Complex polyphase analysis filter, floating point
    C code optimized for HiFi5 with FPU/VFPU
  IntegrIT, 2006-2019
*/

#include <string.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "common_fpu.h"
/* Internal declarations for polyphase FIR filters. */
#include "ppfir_internal.h"

#define ALIGN_SIZE     (HIFI_SIMD_WIDTH)
#define ALIGN_PAD      (ALIGN_SIZE-1)
#define ALIGN_PTR(p)   (void*)(((uintptr_t)(p)+ALIGN_PAD)&~ALIGN_PAD)

#define sz_f32         sizeof(float32_t)
#define sz_cf32        sizeof(complex_float)

/*-------------------------------------------------------------------------
  Polyphase analysis filter 
  Polyphase analysis filters perform filtering of input sample block x[M].
  Each input sample is passed through a dedicated subfilter (altogether
  M subfilters, each of length N). Output samples are organized into 
  type-1 polyphase decomposition frame y[M], that is the output sample
  of 0th subfilter is stored to y[0], output sample of 1st subfilter is
  stored to y[1] and so forth. 
  NOTE: 
  1. This is formal description of algorithm, in reality processing
     in done using circular buffers, so user application is not responsible 
     for management of delay lines

  Precision: 
  32x32    32-bit data, 32-bit coefficients, 32-bit outputs
  f        floating point

  Input:
  x[M]     input samples, Q31, floating point
  h[M*N]   filter coefficients in normal order, Q31, floating point
  N        length of subfilter
  M        length of sample block
  
  Output:
  y[M]     output samples, Q31, floating point

  Restrictions:
  x,y   should not overlap
  x,h   aligned on a 16-bytes boundary
  N     4..24, multiplies of 2
  M     32...640, multiplies of 32
-------------------------------------------------------------------------*/

#if (!HAVE_FPU && !HAVE_VFPU)
DISCARD_FUN(size_t           , cppafirf_alloc  , (int M, int N));
DISCARD_FUN(cppafirf_handle_t, cppafirf_init   , (void * objmem, int M, int N, const float32_t * h));
DISCARD_FUN(void             , cppafirf_process, (cppafirf_handle_t handle, complex_float * y, const complex_float * x));
#else
/* Polyphase filter instance structure */
typedef struct tag_cppafirf_t
{
    int                   M; /* Sample block length          */
    int                   N; /* Number of taps per subfilter */
    const float32_t     * h; /* Filter coefficients (M*N)    */
          complex_float * d; /* Delay line (M*N)             */
          complex_float * p; /* Pointer into the delay line  */
} cppafirf_t;

/* Returns: size of memory in bytes to be allocated */
size_t cppafirf_alloc(int M, int N)
{
    NASSERT(0==(M%32) && 32<=M && M<=640);
    NASSERT(0==(N%2) && 4<=N && N<=24);
    return ALIGN_PAD + M*N*sz_f32 + M*N*sz_cf32 + sizeof(cppafirf_t);
} /* cppafirf_alloc() */

/* Returns: handle to the object */
cppafirf_handle_t cppafirf_init(void * objmem, int M, int N, const float32_t * h)
{
    cppafirf_t *ppfir;
    float32_t *ph;
    complex_float *pd;
    void *p = objmem;
    int m, n;
    NASSERT(0==(M%32) && 32<=M && M<=640);
    NASSERT(0==(N%2) && 4<=N && N<=24);
    NASSERT_ALIGN(h, ALIGN_SIZE);
    /* Partition the object memory */
    ph = (float32_t*)ALIGN_PTR(p); p = ph + M*N;
    pd = (complex_float*)p; p = pd + M*N;
    ppfir = (cppafirf_t*)p; p = ppfir + 1;
#ifdef _DEBUG
    NASSERT((int8_t*)p - (int8_t*)objmem <= (int)cppafirf_alloc(M, N));
#endif
    /* Copy filter coefficients in polyphase order (analysis). */
    for ( m=0; m<M; m++ ) {
        for ( n=0; n<N; n++ ) {
            ph[(N-1-n)*M+m] = h[n*M+m];
        } /* n */
    } /* m */
    /* Zero the delay line and initialize the filter instance. */
    memset(pd, 0, M*N*sz_cf32);
    memset(ppfir, 0, sizeof(*ppfir));
    ppfir->M = M;
    ppfir->N = N;
    ppfir->h = ph;
    ppfir->d = pd;
    ppfir->p = pd;
    return ppfir;
} /* cppafirf_init() */

/* Update the delay line and compute filter output */
void cppafirf_process(cppafirf_handle_t handle, complex_float * y, const complex_float * x)
{
    cppafirf_t * ppfir = (cppafirf_t*)handle;
    NASSERT_ALIGN(x, ALIGN_SIZE);
    ppfir->p = cppfirf(y, ppfir->d, ppfir->p, ppfir->h, x, ppfir->M, ppfir->N);
} /* cppafirf_process() */
#endif /* HAVE_FPU || HAVE_VFPU */
