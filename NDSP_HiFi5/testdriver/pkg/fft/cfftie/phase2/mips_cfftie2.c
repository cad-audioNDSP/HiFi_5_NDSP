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
* Test module for testing cycle performance (Complex FFT with Optimized Memory)
*/

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* Filters and transformations API. */
#include LIBRARY_HEADER(fft)
/* Measurement utilities */
#include "mips.h"

#define PROFILE_CFFT(cond,verb,fun,N,suffix,...)                 PROFILE_INVERTED(cond,verb,fun,(mips.out1.suffix,mips.inp0.suffix,mips.inp1.suffix,1,N), \
                                                                 fout,"N=" #N,prf_ptscycle3,(N));

#define PROFILE_BATCH( _PROFILER, fxn, suffix, ... ) {   \
              _PROFILER( isFull, isVerbose,fxn,     8,suffix, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,    16,suffix, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,    32,suffix, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,    64,suffix, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,   128,suffix, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,   256,suffix, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,   512,suffix, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,  1024,suffix, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,  2048,suffix, __VA_ARGS__ ); \
              _PROFILER( 1, isVerbose,fxn,  4096,suffix, __VA_ARGS__ ); \
             }

void mips_cfftie2(int isFull, int isVerbose, FILE* fout)
{
    /* Stage 2 */
    PROFILE_BATCH( PROFILE_CFFT, fft_cplxf_ie        , cf32 );
    PROFILE_BATCH( PROFILE_CFFT, ifft_cplxf_ie       , cf32 );

    PROFILE_BATCH( PROFILE_CFFT, stereo_fft_cplxf_ie , cf32 );
    PROFILE_BATCH( PROFILE_CFFT, stereo_ifft_cplxf_ie, cf32 );
}
