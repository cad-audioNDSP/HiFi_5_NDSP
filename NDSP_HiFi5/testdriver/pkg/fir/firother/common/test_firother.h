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
 * Test procedures for FIR
 */
#ifndef __TEST_FIROTHER_H__
#define __TEST_FIROTHER_H__
#include <string.h>
#include <stdlib.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API */
#include LIBRARY_HEADER(fir)
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"
#include "../../common/testeng_fir.h"
#include "../../common/testeng_fir_old.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

#define MAX_FUNC_NUM  4
/* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }

/* Initializer for a test description structure. */
#define TEST_DESC( fmt, extraParam, dimNum, align, loadFxn, procFxn ) { (fmt),extraParam,NULL, (dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }
#define TEST_DESC_CONVOLVE( fmt, extraParam, align )  { (fmt), (extraParam)|TE_FIR_CONVOLVE_API , NULL, TE_DIM_NUM_2, (align), NULL,NULL, &te_loadFxn_crosscorr, &te_processFxn_crosscorr }
#define TEST_DESC_CROSSCORR( fmt, extraParam, align ) { (fmt), (extraParam)|TE_FIR_CROSSCORR_API, NULL, TE_DIM_NUM_2, (align), NULL,NULL, &te_loadFxn_crosscorr, &te_processFxn_crosscorr }
#define TEST_DESC_AUTOCORR( fmt, extraParam, align )  { (fmt), (extraParam)|TE_FIR_AUTOCORR_API , NULL, TE_DIM_NUM_1, (align), NULL,NULL, &te_loadFxn_autocorr,  &te_processFxn_autocorr  }

#define TEST_DESC_GCCPHAT( fmt, extraParam, api, align )  { (fmt), (extraParam), api, TE_DIM_NUM_1, (align), NULL, NULL, &te_loadFxn_gccphat,  &te_processFxn_gccphat  }

/*-----------------------------------------------
wrapper for LMS to test convergence:
input:
x[M+N*P-1]  far end
r[N*P]      near end
norm[P]     norm factor for P blocks
mu          adapt rate
M           IR length
N           block size
P           number of blocks
output:
h[M]        estimated impulse response
temporary
e[N]        error
-----------------------------------------------*/
void fir_blmsf_convergence( 
                      float32_t * e, float32_t * h, const float32_t * r,
                const float32_t * x, const float32_t* norm, float32_t mu, int N, int M, int P );

void cxfir_blmsf_convergence( complex_float * e, complex_float * h, 
                const complex_float * r,
                const complex_float * x, 
                const float32_t* norm, float32_t mu, 
                int          N, int       M , int P);

void fir_blms16x16_convergence (  int16_t * e, int16_t *  h,
                const int16_t * r,
                const int16_t * x,
                const int16_t * norm,int16_t   mu,
                int       N,   int       M, int P);

void fir_blms16x32_convergence (  int32_t * e, int32_t *  h,
                const int16_t * r,
                const int16_t * x,
                const int32_t * norm,int16_t   mu,
                int       N,   int       M, int P);

void fir_blms32x32_convergence(  int32_t * e, int32_t *  h,
                const int32_t * r,
                const int32_t * x,
                const int32_t * norm, int32_t mu,
                int       N,   int       M, int P);

void fir_blms32x32ep_convergence(  int32_t * e, int32_t *  h,
                const int32_t * r,
                const int32_t * x,
                const int32_t * norm, int32_t mu,
                int       N,   int       M, int P);

void cxfir_blms32x32_convergence (complex_fract32 *  e, complex_fract32 *  h,
                const complex_fract32 *  r,
                const complex_fract32 *  x,
                const int32_t *  norm, int32_t mu,
                int       N,   int       M, int P);

/* check if target function is not available in the target configuration already */
// static int te_create_convergence(tTestEngContext * context);
int te_create_convergence(tTestEngContext * context);

typedef struct
{
  int                 phaseNum;
  const tTestEngDesc *pFirDescr;
  tTestEngTarget      fxns;
  int                 runAlways;   /* 1 - brief & full & sanity, 0 - full only */
  const char*         seqFile;
}
tTbl;

int test_firother( int phaseNum, int isFull, int isVerbose, int breakOnError, const tTbl * tbl, int szTbl );

#endif /* __TEST_FIROTHER_H__ */
