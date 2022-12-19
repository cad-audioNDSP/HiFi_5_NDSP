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

#define PROFILE_INVERTED_FFT_rfft16x16_ie(cond,verb,N, twstep)     PROFILE_INVERTED_FFT( cond,verb,fft_real16x16_ie,N,( mips.out0.ci16, mips.inp1.i16, mips.inp0.ci16, twstep, N, 2),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_rfft32x16_ie(cond,verb,N, twstep, s)  PROFILE_INVERTED_FFT_SC( cond,verb,fft_real32x16_ie,N,s,( mips.out0.ci32, mips.inp1.i32, mips.inp0.ci16, twstep, N, s),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_rfft24x24_ie(cond,verb,N, twstep)     PROFILE_INVERTED_FFT( cond,verb,fft_real24x24_ie,N,( mips.out0.ci32, mips.inp1.i32, mips.inp0.ci32, twstep, N, 3),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_rfft32x16_ie_24p(cond,verb,N, twstep) PROFILE_INVERTED_FFT( cond,verb,fft_real32x16_ie_24p,N,( mips.out0.u8, mips.inp1.u8, mips.inp0.ci16, twstep, N, 3),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_rfft24x24_ie_24p(cond,verb,N, twstep) PROFILE_INVERTED_FFT( cond,verb,fft_real24x24_ie_24p,N,( mips.out0.u8, mips.inp1.u8, mips.inp0.ci32, twstep, N, 1),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_rfft32x32_ie(cond,verb,N, twstep, s)  PROFILE_INVERTED_FFT_SC( cond,verb,fft_real32x32_ie,N,s,( mips.out0.ci32, mips.inp1.i32, mips.inp0.ci32, twstep, N, s),fout,prf_ptscycle3,N );

#define PROFILE_INVERTED_FFT_rifft16x16_ie(cond,verb,N, twstep)     PROFILE_INVERTED_FFT( cond,verb,ifft_real16x16_ie,N,( mips.out0.i16, mips.inp1.ci16, mips.inp0.ci16, twstep, N, 2),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_rifft32x16_ie(cond,verb,N, twstep, s)  PROFILE_INVERTED_FFT_SC( cond,verb,ifft_real32x16_ie,N,s,( mips.out0.i32, mips.inp1.ci32, mips.inp0.ci16, twstep, N, s),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_rifft24x24_ie(cond,verb,N, twstep)     PROFILE_INVERTED_FFT( cond,verb,ifft_real24x24_ie,N,( mips.out0.i32, mips.inp1.ci32, mips.inp0.ci32, twstep, N, 3),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_rifft32x16_ie_24p(cond,verb,N, twstep) PROFILE_INVERTED_FFT( cond,verb,ifft_real32x16_ie_24p,N,( mips.out0.u8, mips.inp1.u8, mips.inp0.ci16, twstep, N, 3),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_rifft24x24_ie_24p(cond,verb,N, twstep) PROFILE_INVERTED_FFT( cond,verb,ifft_real24x24_ie_24p,N,( mips.out0.u8, mips.inp1.u8, mips.inp0.ci32, twstep, N, 1),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_rifft32x32_ie(cond,verb,N, twstep, s)  PROFILE_INVERTED_FFT_SC( cond,verb,ifft_real32x32_ie,N,s,( mips.out0.i32, mips.inp1.ci32, mips.inp0.ci32, twstep, N, s),fout,prf_ptscycle3,N );

#define PROFILE_RFFT(cond,verb,fun,N,suffix1,suffix2,suffix3,...)                 PROFILE_INVERTED(cond,verb,fun,(mips.out1.suffix2,mips.inp0.suffix1,mips.inp1.suffix3,1,N), \
                                                                 fout,"N=" #N,prf_ptscycle3,(N));

void mips_rfftie1(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_INVERTED_FFT_rfft16x16_ie(isFull, isVerbose, 256, 1);
    PROFILE_INVERTED_FFT_rfft16x16_ie(isFull, isVerbose, 512, 1);
    PROFILE_INVERTED_FFT_rfft16x16_ie(1, isVerbose, 1024, 1);

    PROFILE_INVERTED_FFT_rfft32x16_ie(isFull, isVerbose, 256 , 1, 3);
    PROFILE_INVERTED_FFT_rfft32x16_ie(isFull, isVerbose, 256 , 1, 2);
    PROFILE_INVERTED_FFT_rfft32x16_ie(isFull, isVerbose, 512 , 1, 3);
    PROFILE_INVERTED_FFT_rfft32x16_ie(isFull, isVerbose, 512 , 1, 2);
    PROFILE_INVERTED_FFT_rfft32x16_ie(     1, isVerbose, 1024, 1, 3);
    PROFILE_INVERTED_FFT_rfft32x16_ie(     1, isVerbose, 1024, 1, 2);
#if 0 //HiFi3/3z API
    PROFILE_INVERTED_FFT_rfft32x16_ie_24p(isFull, isVerbose, 256, 1);
    PROFILE_INVERTED_FFT_rfft32x16_ie_24p(isFull, isVerbose, 512, 1);
    PROFILE_INVERTED_FFT_rfft32x16_ie_24p(1, isVerbose, 1024,1);
    PROFILE_INVERTED_FFT_rfft24x24_ie(isFull, isVerbose, 256, 1);
    PROFILE_INVERTED_FFT_rfft24x24_ie(isFull, isVerbose, 512, 1);
    PROFILE_INVERTED_FFT_rfft24x24_ie(1, isVerbose, 1024,1);

    PROFILE_INVERTED_FFT_rfft24x24_ie_24p(isFull, isVerbose, 256, 1);
    PROFILE_INVERTED_FFT_rfft24x24_ie_24p(isFull, isVerbose, 512, 1);
    PROFILE_INVERTED_FFT_rfft24x24_ie_24p(1, isVerbose, 1024,1);
#endif
    PROFILE_INVERTED_FFT_rfft32x32_ie(isFull, isVerbose, 256, 1, 3);
    PROFILE_INVERTED_FFT_rfft32x32_ie(isFull, isVerbose, 256, 1, 2);
    PROFILE_INVERTED_FFT_rfft32x32_ie(isFull, isVerbose, 512, 1, 3);
    PROFILE_INVERTED_FFT_rfft32x32_ie(isFull, isVerbose, 512, 1, 2);
    PROFILE_INVERTED_FFT_rfft32x32_ie(1, isVerbose, 1024, 1, 3);
    PROFILE_INVERTED_FFT_rfft32x32_ie(1, isVerbose, 1024, 1, 2);

    PROFILE_INVERTED_FFT_rifft16x16_ie(isFull, isVerbose, 256, 1);
    PROFILE_INVERTED_FFT_rifft16x16_ie(isFull, isVerbose, 512, 1);
    PROFILE_INVERTED_FFT_rifft16x16_ie(1, isVerbose, 1024, 1);

    PROFILE_INVERTED_FFT_rifft32x16_ie(isFull, isVerbose, 256 , 1, 3);
    PROFILE_INVERTED_FFT_rifft32x16_ie(isFull, isVerbose, 256 , 1, 2);
    PROFILE_INVERTED_FFT_rifft32x16_ie(isFull, isVerbose, 512 , 1, 3);
    PROFILE_INVERTED_FFT_rifft32x16_ie(isFull, isVerbose, 512 , 1, 2);
    PROFILE_INVERTED_FFT_rifft32x16_ie(     1, isVerbose, 1024, 1, 3);
    PROFILE_INVERTED_FFT_rifft32x16_ie(     1, isVerbose, 1024, 1, 2);
#if 0 //HiFi3/3z API
    PROFILE_INVERTED_FFT_rifft32x16_ie_24p(isFull, isVerbose, 256, 1);
    PROFILE_INVERTED_FFT_rifft32x16_ie_24p(isFull, isVerbose, 512, 1);
    PROFILE_INVERTED_FFT_rifft32x16_ie_24p(1, isVerbose, 1024,1);

    PROFILE_INVERTED_FFT_rifft24x24_ie(isFull, isVerbose, 256, 1);
    PROFILE_INVERTED_FFT_rifft24x24_ie(isFull, isVerbose, 512, 1);
    PROFILE_INVERTED_FFT_rifft24x24_ie(1, isVerbose, 1024,1);

    PROFILE_INVERTED_FFT_rifft24x24_ie_24p(isFull, isVerbose, 256, 1);
    PROFILE_INVERTED_FFT_rifft24x24_ie_24p(isFull, isVerbose, 512, 1);
    PROFILE_INVERTED_FFT_rifft24x24_ie_24p(1, isVerbose, 1024,1);
#endif
    PROFILE_INVERTED_FFT_rifft32x32_ie(isFull, isVerbose, 256, 1, 3);
    PROFILE_INVERTED_FFT_rifft32x32_ie(isFull, isVerbose, 256, 1, 2);
    PROFILE_INVERTED_FFT_rifft32x32_ie(isFull, isVerbose, 512, 1, 3);
    PROFILE_INVERTED_FFT_rifft32x32_ie(isFull, isVerbose, 512, 1, 2);
    PROFILE_INVERTED_FFT_rifft32x32_ie(1, isVerbose, 1024, 1, 3);
    PROFILE_INVERTED_FFT_rifft32x32_ie(1, isVerbose, 1024, 1, 2);
}
