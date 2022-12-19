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
* Test module for testing cycle performance (Complex FFT)
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
 
#define PROFILE_INVERTED_FFT_cfft16x16(cond,verb,N,i) PROFILE_INVERTED_FFT_SC   ( cond,verb,fft_cplx16x16,N,i,( mips.out0.i16, mips.inp1.i16, cfft16_##N, i),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_cfft24x24(cond,verb,N,i)                                                                       \
{                                                                                                                           \
    int k;                                                                                                                  \
    if (i == 2) /* For 32-bit scaling on the first stage, need initialize a input array by small value for correct test  */ \
    for (k = 0; k < 2 * N; k++) mips.inp1.i32[k] = k;                                                                            \
    PROFILE_INVERTED_FFT_SC( cond,verb,fft_cplx24x24,N,i,( mips.out0.i32, mips.inp1.i32, cfft24_##N, i),fout,prf_ptscycle3,N );   \
};
#define PROFILE_INVERTED_FFT_cfft32x16(cond,verb,N,i)   PROFILE_INVERTED_FFT_SC( cond,verb,fft_cplx32x16,N,i,( mips.out0.i32, mips.inp1.i32, cfft16_##N, i),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_cfft32x32(cond,verb,N,i)   PROFILE_INVERTED_FFT_SC( cond,verb,fft_cplx32x32,N,i,( mips.out0.i32, mips.inp1.i32, cfft32_##N, i),fout,prf_ptscycle3,N );

#define PROFILE_INVERTED_FFT_cifft16x16(cond,verb,N,i) PROFILE_INVERTED_FFT_SC   ( cond,verb,ifft_cplx16x16,N,i,( mips.out0.i16, mips.inp1.i16, cifft16_##N, i),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_cifft24x24(cond,verb,N,i)                                                                              \
{                                                                                                                                   \
    int k;                                                                                                                          \
    if (i == 2) /* For 32-bit scaling on the first stage, need initialize a input array by small value for correct test  */         \
    for (k = 0; k < 2 * N; k++) mips.inp1.i32[k] = k;                                                                                    \
    PROFILE_INVERTED_FFT_SC(cond, verb, ifft_cplx24x24, N, i, (mips.out0.i32, mips.inp1.i32, cifft24_##N, i), fout, prf_ptscycle3, N)     \
}; 
#define PROFILE_INVERTED_FFT_cifft32x16(cond,verb,N,i)   PROFILE_INVERTED_FFT_SC( cond,verb,ifft_cplx32x16,N,i,( mips.out0.i32, mips.inp1.i32, cifft16_##N, i),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_cifft32x32(cond,verb,N,i)   PROFILE_INVERTED_FFT_SC( cond,verb,ifft_cplx32x32,N,i,( mips.out0.i32, mips.inp1.i32, cifft32_##N, i),fout,prf_ptscycle3,N );

void mips_cfft1(int isFull, int isVerbose, FILE* fout)
{
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 16, 3);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 16, 2);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 32, 3);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 32, 2);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 64, 3);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 64, 2);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 128, 3);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 128, 2);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 256, 3);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 256, 2);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 512, 3);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 512, 2);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 1024, 3);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 1024, 2);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 2048, 3);
    PROFILE_INVERTED_FFT_cfft16x16(isFull, isVerbose, 2048, 2);
    PROFILE_INVERTED_FFT_cfft16x16(1     , isVerbose, 4096, 3);
    PROFILE_INVERTED_FFT_cfft16x16(1     , isVerbose, 4096, 2);
#if 0 //HiFi3/3z API
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 16, 0);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 16, 1);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 16, 2);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 16, 3);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 32, 0);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 32, 1);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 32, 2);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 32, 3);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 64, 0);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 64, 1);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 64, 2);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 64, 3);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 128, 0);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 128, 1);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 128, 2);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 128, 3);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 256, 0);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 256, 1);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 256, 2);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 256, 3);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 512, 0);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 512, 1);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 512, 2);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 512, 3);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 1024, 0);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 1024, 1);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 1024, 2);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 1024, 3);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 2048, 0);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 2048, 1);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 2048, 2);
    PROFILE_INVERTED_FFT_cfft24x24(isFull, isVerbose, 2048, 3);
    PROFILE_INVERTED_FFT_cfft24x24(1, isVerbose, 4096, 0);
    PROFILE_INVERTED_FFT_cfft24x24(1, isVerbose, 4096, 1);
    PROFILE_INVERTED_FFT_cfft24x24(1, isVerbose, 4096, 2);
    PROFILE_INVERTED_FFT_cfft24x24(1, isVerbose, 4096, 3);
#endif

    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 16, 3);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 16, 2);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 32, 3);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 32, 2);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 64, 3);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 64, 2);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 128, 3);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 128, 2);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 256, 3);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 256, 2);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 512, 3);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 512, 2);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 1024, 3);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 1024, 2);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 2048, 3);
    PROFILE_INVERTED_FFT_cfft32x16(isFull, isVerbose, 2048, 2);
    PROFILE_INVERTED_FFT_cfft32x16(1     , isVerbose, 4096, 3);
    PROFILE_INVERTED_FFT_cfft32x16(1     , isVerbose, 4096, 2);

    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 16, 3);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 16, 2);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 32, 3);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 32, 2);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 64, 3);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 64, 2);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 128, 3);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 128, 2);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 256, 3);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 256, 2);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 512, 3);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 512, 2);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 1024, 3);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 1024, 2);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 2048, 3);
    PROFILE_INVERTED_FFT_cfft32x32(isFull, isVerbose, 2048, 2);
    PROFILE_INVERTED_FFT_cfft32x32(1     , isVerbose, 4096, 3);
    PROFILE_INVERTED_FFT_cfft32x32(1     , isVerbose, 4096, 2);

    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 16, 3);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 16, 2);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 32, 3);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 32, 2);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 64, 3);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 64, 2);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 128, 3);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 128, 2);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 256, 3);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 256, 2);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 512, 3);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 512, 2);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 1024, 3);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 1024, 2);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 2048, 3);
    PROFILE_INVERTED_FFT_cifft16x16(isFull, isVerbose, 2048, 2);
    PROFILE_INVERTED_FFT_cifft16x16(1     , isVerbose, 4096, 3);
    PROFILE_INVERTED_FFT_cifft16x16(1     , isVerbose, 4096, 2);
#if 0 //HiFi3/3z API
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 16, 0);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 16, 1);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 16, 2);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 16, 3);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 32, 0);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 32, 1);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 32, 2);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 32, 3);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 64, 0);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 64, 1);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 64, 2);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 64, 3);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 128, 0);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 128, 1);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 128, 2);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 128, 3);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 256, 0);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 256, 1);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 256, 2);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 256, 3);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 512, 0);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 512, 1);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 512, 2);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 512, 3);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 1024, 0);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 1024, 1);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 1024, 2);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 1024, 3);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 2048, 0);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 2048, 1);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 2048, 2);
    PROFILE_INVERTED_FFT_cifft24x24(isFull, isVerbose, 2048, 3);
    PROFILE_INVERTED_FFT_cifft24x24(1, isVerbose, 4096, 0);
    PROFILE_INVERTED_FFT_cifft24x24(1, isVerbose, 4096, 1);
    PROFILE_INVERTED_FFT_cifft24x24(1, isVerbose, 4096, 2);
    PROFILE_INVERTED_FFT_cifft24x24(1, isVerbose, 4096, 3);
#endif
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 16, 3);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 16, 2);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 32, 3);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 32, 2);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 64, 3);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 64, 2);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 128, 3);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 128, 2);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 256, 3);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 256, 2);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 512, 3);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 512, 2);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 1024, 3);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 1024, 2);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 2048, 3);
    PROFILE_INVERTED_FFT_cifft32x16(isFull, isVerbose, 2048, 2);
    PROFILE_INVERTED_FFT_cifft32x16(1     , isVerbose, 4096, 3);
    PROFILE_INVERTED_FFT_cifft32x16(1     , isVerbose, 4096, 2);

    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 16, 3);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 16, 2);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 32, 3);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 32, 2);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 64, 3);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 64, 2);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 128, 3);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 128, 2);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 256, 3);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 256, 2);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 512, 3);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 512, 2);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 1024, 3);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 1024, 2);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 2048, 3);
    PROFILE_INVERTED_FFT_cifft32x32(isFull, isVerbose, 2048, 2);
    PROFILE_INVERTED_FFT_cifft32x32(1     , isVerbose, 4096, 3);
    PROFILE_INVERTED_FFT_cifft32x32(1     , isVerbose, 4096, 2);
}
