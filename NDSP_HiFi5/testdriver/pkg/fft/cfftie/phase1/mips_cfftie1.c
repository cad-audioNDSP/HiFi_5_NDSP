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
#define PROFILE_CFFT_SCL(cond,verb,fun,N,suffix,scl_mtd)       { int bexp; \
                                                PROFILE_INVERTED(cond,verb,fun,(mips.out1.suffix,mips.inp0.suffix,mips.inp1.suffix,1,N,&bexp,scl_mtd), \
                                                                 fout,"N=" #N " [scl=" #scl_mtd "]", prf_ptscycle3,(N)) };

#define PROFILE_INVERTED_FFT_cfft16x16_ie(cond,verb,N, twstep)    PROFILE_INVERTED_FFT( cond,verb,fft_cplx16x16_ie,N,( mips.out0.ci16, mips.inp1.ci16, mips.inp0.ci16, twstep, N, 2),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_cfft32x16_ie(cond,verb,N, twstep, s) PROFILE_INVERTED_FFT_SC( cond,verb,fft_cplx32x16_ie,N,s,( mips.out0.ci32, mips.inp1.ci32, mips.inp0.ci16, twstep, N, s),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_cfft24x24_ie(cond,verb,N, twstep)    PROFILE_INVERTED_FFT( cond,verb,fft_cplx24x24_ie,N,( mips.out0.ci32, mips.inp1.ci32, mips.inp0.ci32, twstep, N, 3),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_cfft32x32_ie(cond,verb,N, twstep, s) PROFILE_INVERTED_FFT_SC( cond,verb,fft_cplx32x32_ie,N,s,( mips.out0.ci32, mips.inp1.ci32, mips.inp0.ci32, twstep, N, s),fout,prf_ptscycle3,N );

#define PROFILE_INVERTED_FFT_cifft16x16_ie(cond,verb,N, twstep)    PROFILE_INVERTED_FFT( cond,verb,ifft_cplx16x16_ie,N,( mips.out0.ci16, mips.inp1.ci16, mips.inp0.ci16, twstep, N, 2),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_cifft32x16_ie(cond,verb,N, twstep, s) PROFILE_INVERTED_FFT_SC( cond,verb,ifft_cplx32x16_ie,N,s,( mips.out0.ci32, mips.inp1.ci32, mips.inp0.ci16, twstep, N, s),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_cifft24x24_ie(cond,verb,N, twstep)    PROFILE_INVERTED_FFT( cond,verb,ifft_cplx24x24_ie,N,( mips.out0.ci32, mips.inp1.ci32, mips.inp0.ci32, twstep, N, 3),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_cifft32x32_ie(cond,verb,N, twstep, s) PROFILE_INVERTED_FFT_SC( cond,verb,ifft_cplx32x32_ie,N,s,( mips.out0.ci32, mips.inp1.ci32, mips.inp0.ci32, twstep, N, s),fout,prf_ptscycle3,N );

#define PROFILE_INVERTED_FFT_stereo_cfft16x16_ie(cond,verb,N, twstep)    PROFILE_INVERTED_FFT( cond,verb,stereo_fft_cplx16x16_ie,N,( mips.out0.ci16, mips.inp1.ci16, mips.inp0.ci16, twstep, N, 2),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_stereo_cfft32x16_ie(cond,verb,N, twstep, s) PROFILE_INVERTED_FFT_SC( cond,verb,stereo_fft_cplx32x16_ie,N,s,( mips.out0.ci32, mips.inp1.ci32, mips.inp0.ci16, twstep, N, s),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_stereo_cfft32x32_ie(cond,verb,N, twstep, s) PROFILE_INVERTED_FFT_SC( cond,verb,stereo_fft_cplx32x32_ie,N,s,( mips.out0.ci32, mips.inp1.ci32, mips.inp0.ci32, twstep, N, s),fout,prf_ptscycle3,N );

#define PROFILE_INVERTED_FFT_stereo_cifft16x16_ie(cond,verb,N, twstep)    PROFILE_INVERTED_FFT( cond,verb,stereo_ifft_cplx16x16_ie,N,( mips.out0.ci16, mips.inp1.ci16, mips.inp0.ci16, twstep, N, 2),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_stereo_cifft32x16_ie(cond,verb,N, twstep, s) PROFILE_INVERTED_FFT_SC( cond,verb,stereo_ifft_cplx32x16_ie,N,s,( mips.out0.ci32, mips.inp1.ci32, mips.inp0.ci16, twstep, N, s),fout,prf_ptscycle3,N );
#define PROFILE_INVERTED_FFT_stereo_cifft32x32_ie(cond,verb,N, twstep, s) PROFILE_INVERTED_FFT_SC( cond,verb,stereo_ifft_cplx32x32_ie,N,s,( mips.out0.ci32, mips.inp1.ci32, mips.inp0.ci32, twstep, N, s),fout,prf_ptscycle3,N );


void mips_cfftie1(int isFull, int isVerbose, FILE* fout)
{
    /*  Stage 1*/
    PROFILE_INVERTED_FFT_cfft16x16_ie(isFull, isVerbose, 128, 1);
    PROFILE_INVERTED_FFT_cfft16x16_ie(isFull, isVerbose, 256, 1);
    PROFILE_INVERTED_FFT_cfft16x16_ie(isFull, isVerbose, 512, 1);
    PROFILE_INVERTED_FFT_cfft16x16_ie(     1, isVerbose, 1024, 1);

    PROFILE_INVERTED_FFT_cfft32x16_ie(isFull, isVerbose,256 , 1, 3);
    PROFILE_INVERTED_FFT_cfft32x16_ie(isFull, isVerbose,256 , 1, 2);
    PROFILE_INVERTED_FFT_cfft32x16_ie(isFull, isVerbose,512 , 1, 3);
    PROFILE_INVERTED_FFT_cfft32x16_ie(isFull, isVerbose,512 , 1, 2);
    PROFILE_INVERTED_FFT_cfft32x16_ie(     1, isVerbose,1024, 1, 3);
    PROFILE_INVERTED_FFT_cfft32x16_ie(     1, isVerbose,1024, 1, 2);
#if 0 //HiFi3/3z API
    PROFILE_INVERTED_FFT_cfft24x24_ie(isFull, isVerbose,256, 1);
    PROFILE_INVERTED_FFT_cfft24x24_ie(isFull, isVerbose,512, 1);
    PROFILE_INVERTED_FFT_cfft24x24_ie(     1, isVerbose,1024,1);
#endif
    PROFILE_INVERTED_FFT_cfft32x32_ie(isFull, isVerbose, 128,  1, 3);
    PROFILE_INVERTED_FFT_cfft32x32_ie(isFull, isVerbose, 128,  1, 2);
    PROFILE_INVERTED_FFT_cfft32x32_ie(isFull, isVerbose, 256,  1, 3);
    PROFILE_INVERTED_FFT_cfft32x32_ie(isFull, isVerbose, 256,  1, 2);
    PROFILE_INVERTED_FFT_cfft32x32_ie(isFull, isVerbose, 512,  1, 3);
    PROFILE_INVERTED_FFT_cfft32x32_ie(isFull, isVerbose, 512,  1, 2);
    PROFILE_INVERTED_FFT_cfft32x32_ie(1,      isVerbose, 1024, 1, 3);
    PROFILE_INVERTED_FFT_cfft32x32_ie(1,      isVerbose, 1024, 1, 2);

    PROFILE_INVERTED_FFT_cifft16x16_ie(isFull, isVerbose, 128, 1);
    PROFILE_INVERTED_FFT_cifft16x16_ie(isFull, isVerbose, 256, 1);
    PROFILE_INVERTED_FFT_cifft16x16_ie(isFull, isVerbose, 512, 1);
    PROFILE_INVERTED_FFT_cifft16x16_ie(     1, isVerbose, 1024, 1);

    PROFILE_INVERTED_FFT_cifft32x16_ie(isFull, isVerbose,256 , 1, 3);
    PROFILE_INVERTED_FFT_cifft32x16_ie(isFull, isVerbose,256 , 1, 2);
    PROFILE_INVERTED_FFT_cifft32x16_ie(isFull, isVerbose,512 , 1, 3);
    PROFILE_INVERTED_FFT_cifft32x16_ie(isFull, isVerbose,512 , 1, 2);
    PROFILE_INVERTED_FFT_cifft32x16_ie(     1, isVerbose,1024, 1, 3);
    PROFILE_INVERTED_FFT_cifft32x16_ie(     1, isVerbose,1024, 1, 2);
#if 0 //HiFi3/3z API
    PROFILE_INVERTED_FFT_cifft24x24_ie(isFull, isVerbose,256 , 1);
    PROFILE_INVERTED_FFT_cifft24x24_ie(isFull, isVerbose,512 , 1);
    PROFILE_INVERTED_FFT_cifft24x24_ie(     1, isVerbose,1024, 1);
#endif
    PROFILE_INVERTED_FFT_cifft32x32_ie(isFull, isVerbose, 128, 1, 3);
    PROFILE_INVERTED_FFT_cifft32x32_ie(isFull, isVerbose, 128, 1, 2);
    PROFILE_INVERTED_FFT_cifft32x32_ie(isFull, isVerbose, 256, 1, 3);
    PROFILE_INVERTED_FFT_cifft32x32_ie(isFull, isVerbose, 256, 1, 2);
    PROFILE_INVERTED_FFT_cifft32x32_ie(isFull, isVerbose, 512, 1, 3);
    PROFILE_INVERTED_FFT_cifft32x32_ie(isFull, isVerbose, 512, 1, 2);
    PROFILE_INVERTED_FFT_cifft32x32_ie(1, isVerbose, 1024, 1, 3);
    PROFILE_INVERTED_FFT_cifft32x32_ie(1, isVerbose, 1024, 1, 2);

    PROFILE_INVERTED_FFT_stereo_cfft16x16_ie(isFull, isVerbose, 128, 1);
    PROFILE_INVERTED_FFT_stereo_cfft16x16_ie(isFull, isVerbose, 256, 1);
    PROFILE_INVERTED_FFT_stereo_cfft16x16_ie(isFull, isVerbose, 512, 1);
    PROFILE_INVERTED_FFT_stereo_cfft16x16_ie(     1, isVerbose, 1024, 1);

    PROFILE_INVERTED_FFT_stereo_cfft32x16_ie(isFull, isVerbose,256 , 1, 3);
    PROFILE_INVERTED_FFT_stereo_cfft32x16_ie(isFull, isVerbose,256 , 1, 2);
    PROFILE_INVERTED_FFT_stereo_cfft32x16_ie(isFull, isVerbose,512 , 1, 3);
    PROFILE_INVERTED_FFT_stereo_cfft32x16_ie(isFull, isVerbose,512 , 1, 2);
    PROFILE_INVERTED_FFT_stereo_cfft32x16_ie(     1, isVerbose,1024, 1, 3);
    PROFILE_INVERTED_FFT_stereo_cfft32x16_ie(     1, isVerbose,1024, 1, 2);

    PROFILE_INVERTED_FFT_stereo_cfft32x32_ie(isFull, isVerbose, 128,  1, 3);
    PROFILE_INVERTED_FFT_stereo_cfft32x32_ie(isFull, isVerbose, 128,  1, 2);
    PROFILE_INVERTED_FFT_stereo_cfft32x32_ie(isFull, isVerbose, 256,  1, 3);
    PROFILE_INVERTED_FFT_stereo_cfft32x32_ie(isFull, isVerbose, 256,  1, 2);
    PROFILE_INVERTED_FFT_stereo_cfft32x32_ie(isFull, isVerbose, 512,  1, 3);
    PROFILE_INVERTED_FFT_stereo_cfft32x32_ie(isFull, isVerbose, 512,  1, 2);
    PROFILE_INVERTED_FFT_stereo_cfft32x32_ie(1,      isVerbose, 1024, 1, 3);
    PROFILE_INVERTED_FFT_stereo_cfft32x32_ie(1,      isVerbose, 1024, 1, 2);

    PROFILE_INVERTED_FFT_stereo_cifft16x16_ie(isFull, isVerbose, 128, 1);
    PROFILE_INVERTED_FFT_stereo_cifft16x16_ie(isFull, isVerbose, 256, 1);
    PROFILE_INVERTED_FFT_stereo_cifft16x16_ie(isFull, isVerbose, 512, 1);
    PROFILE_INVERTED_FFT_stereo_cifft16x16_ie(     1, isVerbose, 1024, 1);

    PROFILE_INVERTED_FFT_stereo_cifft32x16_ie(isFull, isVerbose,256 , 1, 3);
    PROFILE_INVERTED_FFT_stereo_cifft32x16_ie(isFull, isVerbose,256 , 1, 2);
    PROFILE_INVERTED_FFT_stereo_cifft32x16_ie(isFull, isVerbose,512 , 1, 3);
    PROFILE_INVERTED_FFT_stereo_cifft32x16_ie(isFull, isVerbose,512 , 1, 2);
    PROFILE_INVERTED_FFT_stereo_cifft32x16_ie(     1, isVerbose,1024, 1, 3);
    PROFILE_INVERTED_FFT_stereo_cifft32x16_ie(     1, isVerbose,1024, 1, 2);

    PROFILE_INVERTED_FFT_stereo_cifft32x32_ie(isFull, isVerbose, 128, 1, 3);
    PROFILE_INVERTED_FFT_stereo_cifft32x32_ie(isFull, isVerbose, 128, 1, 2);
    PROFILE_INVERTED_FFT_stereo_cifft32x32_ie(isFull, isVerbose, 256, 1, 3);
    PROFILE_INVERTED_FFT_stereo_cifft32x32_ie(isFull, isVerbose, 256, 1, 2);
    PROFILE_INVERTED_FFT_stereo_cifft32x32_ie(isFull, isVerbose, 512, 1, 3);
    PROFILE_INVERTED_FFT_stereo_cifft32x32_ie(isFull, isVerbose, 512, 1, 2);
    PROFILE_INVERTED_FFT_stereo_cifft32x32_ie(1, isVerbose, 1024, 1, 3);
    PROFILE_INVERTED_FFT_stereo_cifft32x32_ie(1, isVerbose, 1024, 1, 2);
}
