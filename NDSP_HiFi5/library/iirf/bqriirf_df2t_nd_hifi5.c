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
  NatureDSP Signal Processing Library. IIR part
    Bi-quad Real Block IIR, floating point, Direct Form II transposed
    Reference C code
  IntegrIT, 2006-2018
*/

/*-------------------------------------------------------------------------
  Bi-quad Real Block IIR
  Computes a real IIR filter (cascaded IIR direct form I or II using 5 
  coefficients per bi-quad + gain term). Real input data are stored
  in vector x. The filter output result is stored in vector r. The filter 
  calculates N output samples using SOS and G matrices.
  NOTE:
  1. Bi-quad coefficients may be derived from standard SOS and G matrices 
  generated by MATLAB. However, typically biquad stages have big peaks 
  in their step response which may cause undesirable overflows at the 
  intermediate outputs. To avoid that the additional scale factors 
  coef_g[M] may be applied. These per-section scale factors may require 
  some tuning to find a compromise between quantization noise and possible 
  overflows. Output of the last section is directed to an additional 
  multiplier, with the gain factor being a power of two, either negative 
  or non-negative. It is specified through the total gain shift amount 
  parameter gain of each filter initialization function.
  2. 16x16 filters may suffer more from accumulation of the roundoff errors,
  so filters should be properly designed to match noise requirements
  3. Due to the performance reasons, IIR biquad filters may introduce 
  additional algorithmic delay of several sampless. Amount of that delay
  might be requested by the  xxx_groupDelay API. For sensitive applications
  all the filters have delayless implementations (with  _nd  suffix in the name).
  Formally, the xxx_groupDelay APIs is also implemented for that kind of filters,
  but return zero.
  
  Precision: 
  16x16         16-bit data, 16-bit coefficients, 16-bit intermediate 
                stage outputs (DF1, DF1 stereo, DF II form)
  32x16         32-bit data, 16-bit coefficients, 32-bit intermediate 
                (DF1, DF1 stereo, DF II form) stage outputs
  32x32         32-bit data, 32-bit coefficients, 32-bit intermediate 
                (DF I, DF1 stereo,  DF II form) stage outputs 
  f             floating point (DF I, DF1 stereo, DF II and DF IIt)

  Input:
  N             length of input sample block
  M             number of bi-quad sections
  S             1 for mono, 2 for stereo API
  s[]           scratch memory area (for fixed-point functions only). 
                Minimum number of bytes depends on selected filter structure 
                and precision:
                  BQRIIR16X16_DF1_SCRATCH_SIZE( N,  M ), or
                  BQRIIR16X16_DF2_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X16_DF1_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X16_DF2_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X32_DF1_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X32_DF2_SCRATCH_SIZE( N,  M ),
                  STEREO_BQRIIR16X16_DF1_SCRATCH_SIZE( N, M ) , or
                  STEREO_BQRIIR32X32_DF1_SCRATCH_SIZE( N, M ) , or
                  STEREO_BQRIIRF_DF1_SCRATCH_SIZE    ( N, M )
                 If a particular macro returns zero, then the corresponding
                 IIR doesn't require a scratch area and parameter s may 
                 hold zero.
  coef_sos[M*5]  filter coefficients stored in blocks of 5 numbers: 
                 b0 b1 b2 a1 a2. 
                 For fixed-point funcions, fixed point format of filter 
                 coefficients is Q1.14 for 16x16 and 32x16, or Q1.30 for 32x32 
  coef_sosl[M*5] filter coefficients for the left channel (stereo filters only)
  coef_sosr[M*5] filter coefficients for the left channel (stereo filters only)
  coef_g[M]      scale factor for each section, Q15 (for fixed-point 
                 functions only)
  coef_gl[M]     scale factor for the left channel (stereo filters only)
  coef_gr[M]     scale factor for the right channel (stereo filters only)
  gain           total gain shift amount applied to output signal of the
                 last section, -48..15
  gainl          total gain shift amount  for the left channel (stereo filters 
                 only)
  gainr          total gain shift amount for the right channel (stereo filters 
                 only)

  x[N*S]         input samples, Q15, Q31 or floating point. Stereo samples 
                 go in interleaved order (left, right)
  Output:
  r[N*S]         output data, Q15, Q31 or floating point. Stereo samples go 
                 in interleaved order (left, right) 

  Restriction:
  x,r,s,coef_g,coef_sos should not overlap
  N   - must be a multiple of 2
  s[] - whenever supplied must be aligned on an 16-bytes boundary
-------------------------------------------------------------------------*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_iir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "bqriirf_df2t_nd.h"
#include "common_fpu.h"

#if (HAVE_VFPU==0 && HAVE_FPU==0)
DISCARD_FUN(size_t,bqriirf_df2t_nd_alloc,( int M ))
DISCARD_FUN(size_t, bqriirf_df2t_nd_groupDelay, (bqriirf_df2t_nd_handle_t  _bqriir))
DISCARD_FUN(bqriirf_df2t_nd_handle_t,bqriirf_df2t_nd_init,( void * objmem, int M,
                                               const float32_t     * coef_sos,
                                               int16_t         gain ))
#elif (HAVE_VFPU)
// Determine the memory area size for a filter instance.
size_t bqriirf_df2t_nd_alloc( int M )
{
  NASSERT( M >= 0 );
  return sizeof(bqriirf_df2t_nd_t)+7*M*sizeof(float32_t)+14;
} // bqriirf_df2t_nd_alloc()

// Given a memory area of sufficient size, initialize a filter instance and 
// return a handle to it.
bqriirf_df2t_nd_handle_t bqriirf_df2t_nd_init( void * objmem, int M,
                                               const float32_t     * coef_sos,
                                               int16_t         gain )
{
#define VLEN 4  /* vector length */
    bqriirf_df2t_nd_t* iir;
    float32_t      * cf; /* 5*M filter coefficients   */
    float32_t      * st; /* 2*M   delay elements      */
    int m;
    st = (float32_t*)((((uintptr_t)objmem)+7)&~7);
    cf = (float32_t*)(((uintptr_t)(st+2*M)+7)&~7);
    iir = (bqriirf_df2t_nd_t*)(cf+5*M);
    iir->M=M;
    iir->st=st;
    iir->cf=cf;
    iir->gain=gain;
    /* clean up and copy coefficients */
    for(m=0; m<2*M; m++) st[m]=0.f;     
    for (m = 0; m<M; m += VLEN)
    {
        int p, P = XT_MIN(VLEN, M - m);
        if (P == VLEN)
        {
            for (p = 0; p < P; p++)
            {
                cf[0 * P + p + m * 5] = coef_sos[0 + (p + m) * 5]; //b0
                cf[1 * P + p + m * 5] = coef_sos[1 + (p + m) * 5]; //b1
                cf[2 * P + p + m * 5] = coef_sos[2 + (p + m) * 5]; //b2
                cf[3 * P + p + m * 5] = coef_sos[3 + (p + m) * 5]; //a1
                cf[4 * P + p + m * 5] = coef_sos[4 + (p + m) * 5]; //a2
            }
        }
        else
        {
            int t = 0;
            if (P & 2)
            {
                t = 2;
                for (p = 0; p < 2; p++)
                {
                    cf[0 * 2 + p + m * 5] = coef_sos[0 + (p + m) * 5]; //b0
                    cf[1 * 2 + p + m * 5] = coef_sos[1 + (p + m) * 5]; //b1
                    cf[2 * 2 + p + m * 5] = coef_sos[2 + (p + m) * 5]; //b2
                    cf[3 * 2 + p + m * 5] = coef_sos[3 + (p + m) * 5]; //a1
                    cf[4 * 2 + p + m * 5] = coef_sos[4 + (p + m) * 5]; //a2
                }
            }
            for (p = 0; p < P%2; p++)
            {
                cf[0 + (p+t + m) * 5] = coef_sos[0 + (p+t + m) * 5]; //b0
                cf[1 + (p+t + m) * 5] = coef_sos[1 + (p+t + m) * 5]; //b1
                cf[2 + (p+t + m) * 5] = coef_sos[2 + (p+t + m) * 5]; //b2
                cf[3 + (p+t + m) * 5] = coef_sos[3 + (p+t + m) * 5]; //a1
                cf[4 + (p+t + m) * 5] = coef_sos[4 + (p+t + m) * 5]; //a2
            }
        }
    }
    return iir;
#undef VLEN
} // bqriirf_df2t_nd_init()

size_t bqriirf_df2t_nd_groupDelay(bqriirf_df2t_nd_handle_t _bqriir)
{
    return 0;
} // bqriirf_df2t_nd_groupDelay()

#else
// for scalar FPU use simpler coefficient layout
size_t bqriirf_df2t_nd_alloc( int M )
{
  NASSERT( M >= 0 );
  return sizeof(bqriirf_df2t_nd_t)+7*M*sizeof(float32_t)+14;
} // bqriirf_df2t_nd_alloc()

// Given a memory area of sufficient size, initialize a filter instance and 
// return a handle to it.
bqriirf_df2t_nd_handle_t bqriirf_df2t_nd_init( void * objmem, int M,
                                               const float32_t     * coef_sos,
                                               int16_t         gain )
{
    bqriirf_df2t_nd_t* iir;
    float32_t      * cf; /* 5*M filter coefficients   */
    float32_t      * st; /* 2*M   delay elements      */
    int m;
    st = (float32_t*)((((uintptr_t)objmem)+7)&~7);
    cf = (float32_t*)(((uintptr_t)(st+2*M)+7)&~7);
    iir = (bqriirf_df2t_nd_t*)(cf+5*M);
    iir->M=M;
    iir->st=st;
    iir->cf=cf;
    iir->gain=gain;
    /* clean up and copy coefficients */
    for(m=0; m<2*M; m++) st[m]=0.f;
    for(m=0; m<5*M; m++) cf[m]=coef_sos[m];
    return iir;
} // bqriirf_df2t_init()

size_t bqriirf_df2t_nd_groupDelay(bqriirf_df2t_nd_handle_t _bqriir)
{
    return 0;
} // bqriirf_df2t_nd_groupDelay()

#endif
