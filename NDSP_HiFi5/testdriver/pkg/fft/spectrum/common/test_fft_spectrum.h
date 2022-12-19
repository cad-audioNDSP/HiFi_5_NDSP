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
 * Test procedures for FFT magnitude functions
 */
#ifndef __TEST_FFT_SPECTRUM_H__
#define __TEST_FFT_SPECTRUM_H__
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
/* Test engine extension for FFT. */
#include "../../common/testeng_fft.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"

/* Test Engine method: allocate vectors and load the data set: */
int loadFxn_fft_spectrum( tTestEngContext * context );

/* Test Engine method: apply the function under test to test case data set. */
void processFxn_fft_spectrum( tTestEngContext * context );

/* Test Engine context extension for auxilliary fft_magnitude() arguments. */
#define fft_mag_aux_bexp   1
#define fft_mag_aux_mode   2

/* Initializer for a test description. */
#define TEST_DESC( fmt )    { (fmt)|FMT_CPLX, 0, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, \
                              NULL, NULL, &loadFxn_fft_spectrum, &processFxn_fft_spectrum }

/* Test definition table. Each entry ties together a function under test (a.k.a.
 * the test target), a SEQ-file containing test data, and a set of attributes
 * (i.e. test description) needed by the Test Engine to perform the test. */
typedef struct tagTestDef 
{
  int            phaseNum;
  const char *   seqFilename;
  tTestEngTarget target;
  tTestEngDesc   desc;

} TestDef_t;

/* Perform all tests for FFT magnitude API functions. */
int main_spectrum( int phaseNum, int isFull, int isVerbose, int breakOnError, const TestDef_t *tbl );

#endif  /* __TEST_FFT_SPECTRUM_H__ */
