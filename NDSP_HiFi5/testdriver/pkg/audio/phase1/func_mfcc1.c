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
 * Test procesdures for MFCC features extraction APIs.
 */

#include <stdio.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */ 
#include "config.h"
#include "packages.h"
/* DSP Library API. */
#include LIBRARY_HEADER(audio)
/* Test engine API. */
#include "testeng.h"
/* Test engine add-on for log mel filterbank and MFCC features extractor tests. */
#include "../common/testeng_logmel.h"
#include "../common/testeng_mfcc.h"

#define MAX(a,b)    ((a)>(b) ? (a) : (b))
#if 1
static const te_logmel_api_t logmel32x32_api = {
    (te_logmel_alloc_fxn_t          *)&logmel32x32_alloc,
    (te_logmel_init_fxn_t           *)&logmel32x32_init,
    (te_logmel_process_fxn_t        *)&logmel32x32_process,
    (te_logmel_getScratchSize_fxn_t *)&logmel32x32_getScratchSize,
};

static const te_mfcc_api_t mfcc32x32_api = {
    (te_mfcc_getDefaultParams_fxn_t *)&mfcc_getDefaultParams,
    (te_mfcc_alloc_fxn_t            *)&mfcc32x32_alloc,
    (te_mfcc_init_fxn_t             *)&mfcc32x32_init,
    (te_mfcc_process_fxn_t          *)&mfcc32x32_process,
    (te_mfcc_getScratchSize_fxn_t   *)&mfcc32x32_getScratchSize
};
#endif
static int  te_loadFxn_htkdelta   ( tTestEngContext * context );
static void te_processFxn_htkdelta( tTestEngContext * context );

#define LOGMEL_TEST_DESC(fmt, align)   { (fmt), 0, NULL, TE_DIM_NUM_1, (align), &te_createFxn_logmel, &te_destroyFxn_logmel, \
                                                                                &te_loadFxn_logmel, &te_processFxn_logmel }
#define MFCC_TEST_DESC(fmt, align)     { (fmt), 0, NULL, TE_DIM_NUM_2, (align), &te_createFxn_mfcc, &te_destroyFxn_mfcc, \
                                                                                &te_loadFxn_mfcc, &te_processFxn_mfcc }
#define HTKDELTA_TEST_DESC(fmt, align) { (fmt), 0, NULL, TE_DIM_NUM_2, (align), NULL, NULL, \
                                                                                &te_loadFxn_htkdelta, &te_processFxn_htkdelta }

static const tTestEngDesc logmel32x32_desc   = LOGMEL_TEST_DESC(FMT_REAL|FMT_FRACT32, TE_ALIGN_YES);
static const tTestEngDesc mfcc32x32_desc     = MFCC_TEST_DESC(FMT_REAL|FMT_FRACT32, TE_ALIGN_YES);
static const tTestEngDesc htkdelta32x32_desc = HTKDELTA_TEST_DESC(FMT_REAL|FMT_FRACT32, TE_ALIGN_YES);

#define DO_TEST( api, desc, seqFile )                                                                   \
    if ( res || !breakOnError ) res &= ( 0 != TestEngRun((tTestEngTarget)&api, &desc, "mfcc1/" seqFile,  \
                                                         isFull, isVerbose, breakOnError, 0) )

/* Perform all tests for MFCC features extraction APIs. */
int func_mfcc1(int isFull, int isVerbose, int breakOnError)
{
    int res = 1;
    DO_TEST(logmel32x32_api, logmel32x32_desc  , "logmel32x32.seq"  );
    DO_TEST(mfcc32x32_api  , mfcc32x32_desc    , "mfcc32x32.seq"    );
    DO_TEST(htkdelta32x32  , htkdelta32x32_desc, "htkdelta32x32.seq");
    return (res);
}

int te_loadFxn_htkdelta( tTestEngContext * context )
{
    int M, N;
    int szX, szZ;
    int res=0;
    /* Retrieve test data dimensions. */
    M = context->args.dim[0];
    N = context->args.dim[1];
    szX = (2*MAX(0, N)+1)*MAX(0, M);
    szZ = MAX(0, M);
    /* Allocate vectors for input data, results and lower/upper result thresholds. */
    if (4!=vecsAlloc(context->desc->isAligned, context->desc->fmt,
                     &context->dataSet.X, szX,
                     &context->dataSet.Z, szZ,
                     &context->dataSet.Zlo, szZ,
                     &context->dataSet.Zhi, szZ, NULL)) {
        printf("te_loadFxn_htkdelta: failed to allocate vectors\n");
    /* Load input data and lower/upper result thresholds from the SEQ-file. */
    } else if (!seqFileReadVecs(context->seqFile, &context->dataSet.X,
                                &context->dataSet.Zlo, &context->dataSet.Zhi, NULL)) {
        printf("te_loadFxn_htkdelta: failed to load vectors from the SEQ-file\n");
    } else {
        res = 1;
    }
    return res;
} /* te_loadFxn_htkdelta() */

void te_processFxn_htkdelta( tTestEngContext * context )
{
    typedef void tFxn( void * d, const void * c, int M, int N );
    const void *pc = vecGetElem(&context->dataSet.X, 0);
    void *pd = vecGetElem(&context->dataSet.Z, 0);
    int M = context->args.dim[0];
    int N = context->args.dim[1];
    te_vReportStd(context);
    ((tFxn*)context->target.fut)(pd, pc, M, N);
} /* te_processFxn_htkdelta() */

