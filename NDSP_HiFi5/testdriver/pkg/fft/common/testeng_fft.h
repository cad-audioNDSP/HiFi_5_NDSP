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
 * Test engine extension for FFT tests
 */

#ifndef __TESTENG_FFT_H
#define __TESTENG_FFT_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API. */
#include LIBRARY_HEADER(fft)
/* Fixed point arithmetics. */
#include "NatureDSP_Math.h"
/* Test engine API. */ 
#include "testeng.h" 
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"
/* Test environment utils. */
#include "utils.h"

#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )
#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )

const void* GetFFT_handle(int N, int FFT_Id);
 
/* Suppress Visual C warnings on +/-INFINITY macros. */
#ifdef COMPILER_MSVC
#pragma warning(disable:4056)
#pragma warning(disable:4756)
#endif

#define sz_fp64c    sizeof(complex_double)


#define NUM_RFFT32X32 43
#define NUM_CFFT32X32 41
#define NUM_RFFT16X16 15
#define NUM_CFFT16X16 15
#define NUM_RFFT32X16 15
#define NUM_CFFT32X16 15
#define NUM_RFFT24X24 9
#define NUM_CFFT24X24 9
#define H_TAB_SZ 18

typedef struct
{
    int N; 
    const fft_handle_t h;
} N_handle_pair_t;

typedef struct
{
    int FFT_Id;
    int numN; 
    const N_handle_pair_t *pNh; 
} FFT_handle_tab_t;

typedef struct
{
    int N;
    void* inc;
    void* twd;
}
tFftDescr_tmp;


#define M_DIM_ dim[2]
#define N_DIM_ dim[0]
#define P_DIM_ dim[3]
#define L_DIM_ dim[1] 

/* Test case types related to FFT. */
#define TE_FFT_TESTCASE_FORWARD    6
#define TE_FFT_TESTCASE_INVERSE    7
#define TE_FFT_TESTCASE_RECONST    8

/* tTestEngDesc::extraParam options for FFT tests. */
#define TE_FFT_REAL             0x0001 /* Forward/inverse real-valued FFT               */
#define TE_FFT_CPLX             0x0002 /* Forward/inverse complex-valued FFT            */
#define TE_FFT_BLOCKWISE        0x0004 /* Blockwise FFT                                 */
#define TE_FFT_FULL_SPECTRUM    0x0008 /* Real FFT forming full symmetric spectrum      */
#define TE_FFT_OPT_SCALE_METH   0x0010 /* Scale method support (fixed point only)       */
#define TE_FFT_OPT_INPLACE      0x0020 /* Input/outbut buffers may be aliased           */
#define TE_FFT_OPT_REUSE_INBUF  0x0040 /* FFT reuses input buffer for intermediate data */
#define TE_DCT                  0x0080 /* DCT type II                                   */
#define TE_FFT_TWD16            0x0100 /* use 16-bit twiddle tables                     */
#define TE_FFT_PACKED24         0x0200 /* packed 24-bit inputs/outputs                  */
#define TE_FFT_UNPACKED24       0x0400 /* unpacked 24-bit inputs/outputs                */
#define TE_FFT_32X16            0x0800 /* 32-bit inputs/outputs, 16-bit twiddles        */
#define TE_MDCT                 0x1000 /* MDCT                                          */
#define TE_DCT2D                0x2000 /* forward DCT-2D                                */
#define TE_IDCT2D               0x4000 /* inverse DCT-2D                                */
#define TE_FFT_STEREO           0x8000 /* FFT with stereo input/output                  */


/* Number of dimensional arguments of a test case (M,N,etc.) */
#define TE_ARGNUM_1    1
#define TE_ARGNUM_2    2
#define TE_ARGNUM_3    3
#define TE_ARGNUM_4    4

/* Make a unique identifier of the FFT function */
#define MAKE_FFT_ID(_testcase_fwd_inv, _extraParams, _dataFormat)\
    (((_testcase_fwd_inv) << 24) | ((_dataFormat) << 16) | ((_extraParams)&(TE_FFT_REAL | TE_FFT_CPLX | TE_FFT_PACKED24 | TE_FFT_UNPACKED24 | TE_FFT_32X16)))




/*
 * Test data frame processing for a standalone test. Apply an FFT routine to 16-bit
 * fixed point data read from the input file and convert the result to 64-bit
 * floating point samples.
 */

#if 1
typedef int tTestEngFrameProcFxn_fft_stnd( tTestEngContext * context,
                                     const int16_t         * in,  /* 16-bit fixed point complex    */
                                           float64_t       * out, /* 64-bit floating point complex */
                                     tVec * xVec, tVec * yVec);	/* In/out vectors for target FFT   */

tTestEngFrameProcFxn_fft_stnd te_frameProc_stnd_fft;     /* <r|c>[i]fft[d], <r|c>[i]fftf[d|_fr16|_fr32] */
tTestEngFrameProcFxn_fft_stnd te_frameProc_stnd_scl_fft; /* <r|c>[i]fft_fr<16|32>                       */
tTestEngFrameProcFxn_fft_stnd te_frameProc_stnd_blkfft;  /* blk<r|c>[i]fft[_fr16|_fr32]                 */

#endif

/* FFT test target */
typedef void (*tFrwTransFxn)( const void * x, void * y, const void * twdTbl,
                     int twdStep, int N, int * bexp, int scl_mtd );
typedef void (*tInvTransFxn)( const void * x, void * y, const void * twdTbl,
                     int twdStep, int N, int * bexp, int scl_mtd );

typedef struct tagTestEngDesc_fft
{
    tTestEngDesc desc;
    tFrwTransFxn frwTransFxn;  /* Forward transform function */
    tInvTransFxn invTransFxn;  /* Inverse transform function */

} tTestEngDesc_fft;

/* FFT test context; is accessible through tTestEngContext::target::handle */
typedef struct tagTestEngContext_fft
{
  char fInName [64]; /* Input data filename                                   */
  char fRefName[64]; /* Reference data filename                               */ 
  int  scale_method; /* Scale method (for feature rich fixed point FFTs only) */
  int  frameCnt;     /* Frame counter                                         */
  size_t dataSize;   /* data size in bytes allocated for input/output         */

} tTestEngContext_fft;

/* FFT test context. */
typedef struct tagTestEngContext_fft_int
{
  tTestEngContext_fft ext; /* Externally visible part of FFT context. */

  struct {                 /* Twiddle factor tables.                  */
    int    baseSize;       /* First twiddle table FFT size.           */
    int    tblNum;         /* Number of twiddle tables.               */
    tVec * tblSet;         /* Twiddle tables, with FFT size doubled   */
  }                        /* in succession.                          */
  twd;

} tTestEngContext_fft_int;

/*
 * Test engine methods for FFT tests.
 */
 
/* Create a target algorithm instance and set tTestEngContext::target fields.
 * In particular, load twiddle factor tables. Return zero if failed. */
int te_create_fft( tTestEngContext * context );

/* Destroy the target algorithm instance and free memory block(s) allocated
 * for the target object. Return zero whenever any violation of memory area
 * bounds is detected. */
int te_destroy_fft( tTestEngContext * context );

/* Allocate in/out vectors for the next test case, and load the data set
 * from the SEQ-file. Return zero if failed. */
int te_load_fft( tTestEngContext * context );

/* Return a pointer to twiddle factor table. If step parameter
 * is non-zero, then the table is selected from the set of available
 * tables in dependence of the test frame counter, with appropriate
 * stride amount returned through step. If step is zero, return the
 * twiddle table such that stride is 1. Return zero if found no table
 * for the requested FFT size. */
void * te_get_twd_tbl( tTestEngContext * context, int fftSize, int * step );

/*
 * Test engine methods for 2D FFT tests.
 */
 
/* Allocate vectors for the next in turn test case, and load the data
 * set from a SEQ-file. Return zero if failed. */
int te_load_fft2d( tTestEngContext * context );

/* Apply the target function to test case data. */
void te_process_fft2d( tTestEngContext * context );

#endif /* __TESTENG_FFT_H */

