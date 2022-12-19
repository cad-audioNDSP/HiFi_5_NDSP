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
    Complex polyphase analysis filter, 32x32-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

#include <string.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Internal declarations for polyphase FIR filters. */
#include "ppfir_internal.h"

#define ALIGN_SIZE     (HIFI_SIMD_WIDTH)
#define ALIGN_PAD      (ALIGN_SIZE-1)
#define ALIGN_PTR(p)   (void*)(((uintptr_t)(p)+ALIGN_PAD)&~ALIGN_PAD)

#define sz_i32         sizeof(int32_t)
#define sz_ci32        sizeof(complex_fract32)

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

/* Polyphase filter instance structure */
typedef struct tag_cppafir32x32_t
{
    int                     M; /* Sample block length          */
    int                     N; /* Number of taps per subfilter */
    const int32_t         * h; /* Filter coefficients (M*N)    */
          complex_fract32 * d; /* Delay line (M*N)             */
          complex_fract32 * p; /* Pointer into the delay line  */
} cppafir32x32_t;

/* Returns: size of memory in bytes to be allocated */
size_t cppafir32x32_alloc(int M, int N)
{
    NASSERT(0==(M%32) && 32<=M && M<=640);
    NASSERT(0==(N%2) && 4<=N && N<=24);
    return ALIGN_PAD + M*N*sz_i32 + M*N*sz_ci32 + sizeof(cppafir32x32_t) + ADD_PAD_CI32*sz_ci32;
} /* cppafir32x32_alloc() */

/* Returns: handle to the object */
cppafir32x32_handle_t cppafir32x32_init(void * objmem, int M, int N, const int32_t * h)
{
    cppafir32x32_t *ppfir;
    int32_t *ph;
    complex_fract32 *pd;
    void *p = objmem;
    int m, n;
    NASSERT(0==(M%32) && 32<=M && M<=640);
    NASSERT(0==(N%2) && 4<=N && N<=24);
    NASSERT_ALIGN(h, ALIGN_SIZE);
    /* Partition the object memory */
    ph = (int32_t*)ALIGN_PTR(p); p = ph + M*N;
    pd = (complex_fract32*)p + ADD_PAD_CI32; p = pd + M*N;
    ppfir = (cppafir32x32_t*)p; p = ppfir + 1;
#ifdef _DEBUG
    NASSERT((int8_t*)p - (int8_t*)objmem <= (int)cppafir32x32_alloc(M, N));
#endif
    /* Copy filter coefficients in polyphase order (analysis). */
    for ( n=0; n<N; n+=2 ) {
        for ( m=0; m<M; m++ ) {
            ph[(N-2-n)*M+2*m+1] = h[n*M  +m];
            ph[(N-2-n)*M+2*m  ] = h[n*M+M+m];
        } /* n */
    } /* m */
    /* Zero the delay line and initialize the filter instance. */
    memset(pd, 0, M*N*sz_ci32);
    memset(ppfir, 0, sizeof(*ppfir));
    ppfir->M = M;
    ppfir->N = N;
    ppfir->h = ph;
    ppfir->d = pd;
    ppfir->p = pd;
    return ppfir;
} /* cppafir32x32_init() */

/* Update the delay line and compute filter output */
void cppafir32x32_process(cppafir32x32_handle_t handle, complex_fract32 * y, const complex_fract32 * x)
{
    cppafir32x32_t * ppfir = (cppafir32x32_t*)handle;
    NASSERT_ALIGN(x, ALIGN_SIZE);
    ppfir->p = cppfir32x32(y, ppfir->d, ppfir->p, ppfir->h, x, ppfir->M, ppfir->N, 0);
} /* cppafir32x32_process() */
