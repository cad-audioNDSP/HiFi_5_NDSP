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
* Test module for testing cycle performance (Real FFT with Optimized Memory)
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

#define PROFILE_RFFT(cond,verb,fun,N,suffix1,suffix2,suffix3,...)                 PROFILE_INVERTED(cond,verb,fun,(mips.out1.suffix2,mips.inp0.suffix1,mips.inp1.suffix3,1,N), \
                                                                 fout,"N=" #N,prf_ptscycle3,(N));

#define PROFILE_BATCH( _PROFILER, fxn,suffix1,suffix2,suffix3, ... ) {   \
              _PROFILER( isFull, isVerbose,fxn,     8,suffix1,suffix2,suffix3, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,    16,suffix1,suffix2,suffix3, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,    32,suffix1,suffix2,suffix3, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,    64,suffix1,suffix2,suffix3, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,   128,suffix1,suffix2,suffix3, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,   256,suffix1,suffix2,suffix3, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,   512,suffix1,suffix2,suffix3, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,  1024,suffix1,suffix2,suffix3, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose,fxn,  2048,suffix1,suffix2,suffix3, __VA_ARGS__ ); \
              _PROFILER( 1, isVerbose,fxn,  4096,suffix1,suffix2,suffix3, __VA_ARGS__ ); \
}

void mips_rfftie2(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_BATCH( PROFILE_RFFT         , fft_realf_ie ,f32,cf32,cf32);
    PROFILE_BATCH( PROFILE_RFFT         ,ifft_realf_ie ,cf32,f32,cf32);
}
