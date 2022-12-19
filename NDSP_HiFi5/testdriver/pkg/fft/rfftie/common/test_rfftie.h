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
* Test procedures for complex FFT functions.
*/
#ifndef __TEST_RFFTIE_H__
#define __TEST_RFFTIE_H__
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
/* Test engine extension for FFT. */
#include "../../common/testeng_fft.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

/* Suppress Visual C warnings on +/-HUGE_VALF macros. */
#ifdef COMPILER_MSVC
#pragma warning(disable:4056)
#pragma warning(disable:4756)
#endif

#define sz_fr16c    sizeof(complex_int16_t)
#define sz_fl64c    sizeof(complex_double)

/* Period, in frames, between reallocation of in/out buffers consumed by the
 * target FFT routine. */
#define XY_VEC_REALLOC_PERIOD     16

/* FFT attributes for feature rich and fast real FFTs. */
#define ATTR_FEATURE_RICH      (TE_FFT_REAL|TE_FFT_OPT_INPLACE|TE_FFT_FULL_SPECTRUM)
#define ATTR_FEATURE_RICH_SCL  (TE_FFT_REAL|TE_FFT_OPT_INPLACE|TE_FFT_FULL_SPECTRUM|TE_FFT_OPT_SCALE_METH)
#define ATTR_FAST              (TE_FFT_REAL|TE_FFT_OPT_REUSE_INBUF)

typedef struct tagTestDef
{
  int                phaseNum;
  const char *       seqFilename;
  tTestEngTarget     target;
  tTestEngDesc_fft   desc;

}TestDef_t;

/* Perform all tests for fft_realMxN_ie, ifft_realMxN_ie API functions. */
int main_rfftie(int phaseNum, int isFull, int isVerbose, int breakOnError, const TestDef_t *tbl);

#endif  /* __TEST_RFFTIE_H__ */
