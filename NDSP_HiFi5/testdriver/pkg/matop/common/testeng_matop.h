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
 * Test-engine add-on for matrix categories 
 */

#ifndef TESTENG_MTX_H__
#define TESTENG_MTX_H__

/* Cross-platform data type definitions. */
#include "types.h"
/* Test engine API. */
#include "testeng.h"

typedef struct {
    size_t (*getScratchSize)(int M, int N, int P);
} tMatmulApi_MNP;

typedef struct {
    size_t (*getScratchSize)(int N);
} tMatmulApi_N;

/* extraParam flags */
#define MTX_PLAIN          0    /* non-streaming data                            */
#define MTX_PLAIN_LSH      1    /* non-streaming data with aditional lsh         */
#define MTX_PLAIN_RSH_LSH  2    /* non-streaming data with aditional rsh and lsh */
#define MTX_8X16_LSH       3    /* non-streaming 8X16 data with aditional lsh    */

int te_loadFxn_mmlt( tTestEngContext * context ); /* Allocate vectors and load the data set for [c]matmmlt: * X[M][N], Y[N][P], Z[M][P] */
void te_processFxn_mmlt( tTestEngContext * context ); /* Apply the function under test to test case data set: [c]matmmlt */
void te_processFxn_vecmmlt( tTestEngContext * context ); /* Apply the function under test to test case data set: vecmpy */

int  te_loadFxn_vecmpyt   ( tTestEngContext * context ); /* Allocate vectors and load the data set for [c]mtx_vecmpyt: X[N], Z[N][N] */
void te_processFxn_vecmpyt( tTestEngContext * context ); /* Apply the function under test to test case data set: [c]mtx_vecmpyt */

int  te_loadFxn_lrmpy   ( tTestEngContext * context ); /* Allocate vectors and load the data set for [c]mtx_lrmpy: X[N][N], Y[N][N], Z[N][N] */
void te_processFxn_lrmpy( tTestEngContext * context ); /* Apply the function under test to test case data set: [c]mtx_lrmpy */

void te_processFxn_mtx_transpose( tTestEngContext * context ); /* Apply the function under test to test case data set: matrix transpose */

#endif

