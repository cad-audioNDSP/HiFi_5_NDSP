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
 * Test module for testing cycle performance (FFT spectrum functions)
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

#define PROFILE_FFT_MAG(isFull,isVerbose,fun,N,mode,suffix,...) \
                            PROFILE_INVERTED(isFull,isVerbose,fun,(mips.out1.suffix,mips.inp0.c##suffix,N,mode), \
                                             fout,"N=" #N "[mode=" #mode "]", \
                                             prf_ptscycle2,(N));

#define PROFILE_FFT_MAG_INPL(isFull,isVerbose,fun,N,mode,suffix,...) \
                            PROFILE_INVERTED(isFull,isVerbose,fun,(mips.out0.suffix,mips.out0.c##suffix,N,mode), \
                                             fout,"N=" #N "[mode=" #mode " inplace]", \
                                             prf_ptscycle2,(N));

#define PROFILE_BATCH( _PROFILER, isFull, isVerbose, fxn, mode, suffix, suffix1,... ) {   \
              _PROFILER( isFull, isVerbose, fxn,     2, mode, suffix,suffix1, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose, fxn,     4, mode, suffix,suffix1, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose, fxn,     8, mode, suffix,suffix1, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose, fxn,    16, mode, suffix,suffix1, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose, fxn,    32, mode, suffix,suffix1, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose, fxn,    64, mode, suffix,suffix1, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose, fxn,   128, mode, suffix,suffix1, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose, fxn,   256, mode, suffix,suffix1, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose, fxn,   512, mode, suffix,suffix1, __VA_ARGS__ ); \
/*              _PROFILER( isFull, isVerbose, fxn,  1024, mode, suffix, __VA_ARGS__ ); */\
              _PROFILER( isFull, isVerbose, fxn,  2048, mode, suffix,suffix1, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose, fxn,  4096, mode, suffix,suffix1, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose, fxn,  8192, mode, suffix,suffix1, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose, fxn, 16384, mode, suffix,suffix1, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose, fxn, 32768, mode, suffix,suffix1, __VA_ARGS__ ); \
              _PROFILER( isFull, isVerbose, fxn, 65536, mode, suffix,suffix1, __VA_ARGS__ ); }

void mips_spectrum2(int isFull, int isVerbose, FILE * fout)
{
  PROFILE_BATCH( PROFILE_FFT_MAG          , isFull, isVerbose, fft_spectrumf     , 0, f32,f32     );
  PROFILE_FFT_MAG     (1     , isVerbose, fft_spectrumf, 1024, 0, f32);
  PROFILE_BATCH( PROFILE_FFT_MAG          , isFull, isVerbose, fft_spectrumf     , 1, f32,f32     );
  PROFILE_FFT_MAG     (isFull, isVerbose, fft_spectrumf, 1024, 1, f32);
  PROFILE_BATCH( PROFILE_FFT_MAG_INPL     , isFull, isVerbose, fft_spectrumf     , 0, f32,f32     );
  PROFILE_FFT_MAG_INPL(isFull, isVerbose, fft_spectrumf, 1024, 0, f32);
  PROFILE_BATCH( PROFILE_FFT_MAG_INPL     , isFull, isVerbose, fft_spectrumf     , 1, f32,f32     );
  PROFILE_FFT_MAG_INPL(isFull, isVerbose, fft_spectrumf, 1024, 1, f32);
}
