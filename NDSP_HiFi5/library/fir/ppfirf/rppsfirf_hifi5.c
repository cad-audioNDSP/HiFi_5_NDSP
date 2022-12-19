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
    Real polyphase synthesis filter, floating point
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

/*-------------------------------------------------------------------------
  Polyphase synthesis filter 
  Polyphase synthesis filter reconstructs a block of M samples y[M] from 
  type-1 polyphase decomposition frame x[M]. First, each sample of the input
  frame x[M] is passed through a dedicated subfilter (altogether M subfilters,
  each of length N). After that, output samples of every subfilter are 
  interleaved such that the first-in-time sample of the reconstructed block
  is stored to y[0], the next-in-time sample – to y[1], and so forth.
  Normally, polyphase decomposition frames x[M] come from IFFT routines which
  may have scaled data to prevent overflow in the case of fixed point data.
  Input argument lsh  accommodates the need to scale the reconstructed signal
  back to the desired representation.

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
  lsh      scale factor (positive for left shift, negative for right 
           shift) for scaling resulted output. For fixed point API only.
  
  Output:
  y[M]     output samples, Q31, floating point

  Restrictions:
  x,y   should not overlap
  x,h   aligned on a 16-bytes boundary
  N     4..24, multiplies of 2
  M     32...640, multiplies of 32
  -------------------------------------------------------------------------*/

#if (!HAVE_FPU && !HAVE_VFPU)
DISCARD_FUN(size_t           , rppsfirf_alloc  , (int M, int N));
DISCARD_FUN(rppsfirf_handle_t, rppsfirf_init   , (void * objmem, int M, int N, const float32_t* h));
DISCARD_FUN(void             , rppsfirf_process, (rppsfirf_handle_t handle, float32_t * y, const float32_t * x));
#else
/* Polyphase filter instance structure */
typedef struct tag_rppsfirf_t
{
    int               M; /* Sample block length          */
    int               N; /* Number of taps per subfilter */
    const float32_t * h; /* Filter coefficients (M*N)    */
          float32_t * d; /* Delay line (M*N)             */
          float32_t * p; /* Pointer into the delay line  */
} rppsfirf_t;

/* Returns: size of memory in bytes to be allocated */
size_t rppsfirf_alloc(int M, int N)
{
    NASSERT(0==(M%32) && 32<=M && M<=640);
    NASSERT(0==(N%2) && 4<=N && N<=24);
    return ALIGN_PAD + 2*M*N*sz_f32 + sizeof(rppsfirf_t) + ADD_PAD_F32*sizeof(float32_t);
} /* rppsfirf_alloc() */

/* Returns: handle to the object */
rppsfirf_handle_t rppsfirf_init(void * objmem, int M, int N, const float32_t* h)
{
    rppsfirf_t *ppfir;
    float32_t *ph, *pd;
    void *p = objmem;
    int m, n;
    NASSERT(0==(M%32) && 32<=M && M<=640);
    NASSERT(0==(N%2) && 4<=N && N<=24);
    NASSERT_ALIGN(h, ALIGN_SIZE);
    /* Partition the object memory */
    ph = (float32_t*)ALIGN_PTR(p); p = ph + M*N;
    pd = (float32_t*)p + ADD_PAD_F32; p = pd + M*N;
    ppfir = (rppsfirf_t*)p; p = ppfir + 1;
#ifdef _DEBUG
    NASSERT((int8_t*)p - (int8_t*)objmem <= (int)rppsfirf_alloc(M, N));
#endif
    /* Copy filter coefficients in polyphase order (synthesis). */
    for ( m=0; m<M; m++ ) {
        for ( n=0; n<N; n++ ) {
            ph[(N-1-n)*M+(M-1-m)] = h[n*M+m];
        } /* n */
    } /* m */
    /* Zero the delay line and initialize the filter instance. */
    memset(pd, 0, M*N*sz_f32);
    memset(ppfir, 0, sizeof(*ppfir));
    ppfir->M = M;
    ppfir->N = N;
    ppfir->h = ph;
    ppfir->d = pd;
    ppfir->p = pd;
    return ppfir;
} /* rppsfirf_init() */

/* Update the delay line and compute filter output */
void rppsfirf_process(rppsfirf_handle_t handle, float32_t * y, const float32_t * x)
{
    rppsfirf_t * ppfir = (rppsfirf_t*)handle;
    NASSERT_ALIGN(x, ALIGN_SIZE);
    ppfir->p = rppfirf(y, ppfir->d, ppfir->p, ppfir->h, x, ppfir->M, ppfir->N);
} /* rppsfirf_process() */
#endif /* HAVE_FPU || HAVE_VFPU */
