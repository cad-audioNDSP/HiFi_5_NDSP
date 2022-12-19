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
 * Test-engine add-on for log mel filterbank APIs.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* DSP Library API. */
#include LIBRARY_HEADER(audio)
/* Test engine API. */
#include "testeng.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"
/* Test engine add-on for log mel filterbank tests. */
#include "testeng_logmel.h"
/* Aligning wrapper over standard malloc/free. */
#include "malloc16.h"


/* Load a logmel params structure from a SEQ-file. Return zero if failed. */
static int loadParams( logmel_params_t * logmelParams, tSeqFile seqFile );

/* Test host encapsulates various resources needed to run a whole SEQ-file. */
typedef struct te_logmel_host_tag {
    logmel_params_t logmelParams;
    int scaleExp;
} te_logmel_host_t;

/* Prepare to run test cases. */
int te_createFxn_logmel( tTestEngContext * context )
{
    te_logmel_api_t * logmelApi = (te_logmel_api_t*)context->target.fut;
    te_logmel_host_t * logmelHost;
    /* Check if the logmel API is implemented for the current build configuration. */
    if (!IS_PRESENT(logmelApi->alloc) || !IS_PRESENT(logmelApi->init) ||
        !IS_PRESENT(logmelApi->process) || !IS_PRESENT(logmelApi->getScratchSize)) return -1;
    /* Setup the host structure. */
    if (NULL==(logmelHost = (te_logmel_host_t*)mallocAlign(sizeof(*logmelHost), 0))) {
        printf("te_createFxn_logmel: failed to allocate the host structure\n"); return 0;
    }
    memset(logmelHost, 0, sizeof(*logmelHost));
    context->target.handle = logmelHost;
    return 1;
} /* te_createFxn_logmel() */

/* Release resources after the last test case, */
int te_destroyFxn_logmel( tTestEngContext * context )
{
    te_logmel_host_t * logmelHost = (te_logmel_host_t*)context->target.handle;
    if (NULL!=logmelHost) freeAlign(logmelHost);
    return 1;
} /* te_destroyFxn_logmel() */

/* Allocate vectors and load the data set for a test case. */
int te_loadFxn_logmel( tTestEngContext * context )
{
    te_logmel_host_t * logmelHost = (te_logmel_host_t*)context->target.handle;
    logmel_params_t * logmelParams;
    int res=0;
    NASSERT(NULL!=logmelHost);
    logmelParams = &logmelHost->logmelParams;
    /* Load the logmel configuration and scale exponent argument. */
    if (0==loadParams(logmelParams, context->seqFile)) {
        printf("te_loadFxn_logmel: failed to load the params structure\n");
    } else if (1!=seqFileScanf(context->seqFile, "%d", &logmelHost->scaleExp)) {
        printf("te_loadFxn_logmel: failed to load scale exponent argument\n");
    /* Allocate the test spectra vector */
    } else if (0==vecAlloc(&context->dataSet.X, logmelParams->fftSize/2+1,  TE_ALIGN_YES, 
                           (context->desc->fmt&FMT_DTYPE_MASK)|FMT_CPLX, NULL)) {
        printf("te_loadFxn_logmel: failed to allocate the input vector\n");
    /* Allocate vectors for logFbe results, and lower/upper threshold vectors */
    } else if (3!=vecsAlloc(TE_ALIGN_YES, (context->desc->fmt&FMT_DTYPE_MASK)|FMT_REAL, 
                            &context->dataSet.Z, logmelParams->mfbBandNum,
                            &context->dataSet.Zlo, logmelParams->mfbBandNum,
                            &context->dataSet.Zhi, logmelParams->mfbBandNum, NULL)) {
        printf("te_loadFxn_logmel: failed to allocate output/threshold vectors\n");
    /* Load test spectra vector and logFbe threshold vectors */
    } else if (0==seqFileReadVecs(context->seqFile, &context->dataSet.X, 
                                  &context->dataSet.Zlo, &context->dataSet.Zhi, NULL)) {
        printf("te_loadFxn_logmel: failed to load vectors from the SEQ-file\n");
    } else {
        res = 1;
    }
    return res;
} /* te_loadFxn_logmel() */

/* Apply the target function to the test case data set. */
void te_processFxn_logmel( tTestEngContext * context )
{
    te_logmel_host_t * logmelHost = (te_logmel_host_t*)context->target.handle;
    te_logmel_api_t * logmelApi = (te_logmel_api_t*)context->target.fut;
    logmel_params_t * logmelParams;
    te_logmel_handle_t logmel;
    size_t szObj, szScr;
    tVec vecObj, vecScr, vecLogFbe, vecTmp;
    void *pObj, *pScr, *spectra, *logFbe;
    NASSERT(NULL!=logmelHost);
    NASSERT(NULL!=logmelApi);
    logmelParams = &logmelHost->logmelParams;
    /* Allocate object and scratch memory areas, then instantiate the logmel filterbank. */
    szObj = logmelApi->alloc(logmelParams);
    szScr = logmelApi->getScratchSize(logmelParams);
    if (2!=vecsAlloc(TE_ALIGN_NO, FMT_UINT8, &vecObj, szObj, &vecScr, szScr, NULL)) {
        printf("te_processFxn_logmel: failed to allocate object and/or scratch memory\n"); return;
    }
    pObj = vecGetElem(&vecObj, 0);
    pScr = vecGetElem(&vecScr, 0);
    if (NULL==(logmel = logmelApi->init(pObj, logmelParams))) {
        printf("te_processFxn_logmel: failed to initialize the log mel filterbank\n"); return;
    }
    /* Allocate the logFbe results vector. */
    if (0==vecAlloc(&vecLogFbe, logmelParams->mfbBandNum, TE_ALIGN_YES, 
                    (context->desc->fmt&FMT_DTYPE_MASK)|FMT_REAL, NULL)) {
        printf("te_processFxn_logmel: failed to allocate the results vector\n"); return;
    }
    /* Add to the log. */
    {
        char str[200];
        tReportFUT fut[3];
        fut[0]=(tReportFUT)logmelApi->init;
        fut[1]=(tReportFUT)logmelApi->alloc;
        fut[2]=(tReportFUT)logmelApi->process;
        sprintf(str,"Bands: %d; opt: 0x%02x", logmelParams->mfbBandNum, (unsigned)logmelParams->opt);
        vReportAdd(fut,3,str,context->seqFile->filename ,context->args.caseType,vecGetSize(&context->dataSet.X));
    }
    /* Process the test data. */
    spectra = vecGetElem(&context->dataSet.X, 0);
    logFbe = vecGetElem(&vecLogFbe, 0);
    logmelApi->process(logmel, pScr, logFbe, spectra, logmelHost->scaleExp);
    /* Move results to a temporal storage. */
    if (0==vecClone(&vecTmp, &vecLogFbe)) {
        printf("te_processFxn_logmel: failed to clone the results vector\n"); return;
    }
    /* Iff there were no errors encountered during cleanup, copy logFbe results from the temporal
     * storage to the output vector. Such an approach ensures that a test case is invalidated 
     * whenever the target function corrupts the memory, even if has managed to produce valid
     * results. */
    if (0!=vecsFree(&vecObj, &vecScr, &vecLogFbe, NULL)) {
        size_t szElem = vecElemSize((context->desc->fmt&FMT_DTYPE_MASK)|FMT_REAL);
        void *pZ = vecGetElem(&context->dataSet.Z, 0);
        const void * pTmp = vecGetElem(&vecTmp, 0);
        memcpy(pZ, pTmp, logmelParams->mfbBandNum*szElem);
        vecFree(&vecTmp);
    }
} /* te_processFxn_logmel() */

/* Load a logmel params structure from a SEQ-file. Return zero if failed. */
int loadParams( logmel_params_t * logmelParams, tSeqFile seqFile )
{
    int mfbLowFreqQ8Int=0, mfbUppFreqQ8Int=0;
    int melScaleOpt=0, fbeLogOpt=0, fbNormOpt=0, res=0;
    memset(logmelParams, 0, sizeof(*logmelParams));
    if (1==seqFileScanf(seqFile, "%d", &logmelParams->Fs        ) &&
        1==seqFileScanf(seqFile, "%d", &logmelParams->fftSize   ) &&
        1==seqFileScanf(seqFile, "%d", &mfbLowFreqQ8Int         ) &&
        1==seqFileScanf(seqFile, "%d", &mfbUppFreqQ8Int         ) &&
        1==seqFileScanf(seqFile, "%d", &logmelParams->mfbBandNum) &&
        1==seqFileScanf(seqFile, "%d", &melScaleOpt             ) &&
        1==seqFileScanf(seqFile, "%d", &fbeLogOpt               ) &&
        1==seqFileScanf(seqFile, "%d", &fbNormOpt               ))
    {
        logmelParams->mfbLowFreqQ8 = (fract32)mfbLowFreqQ8Int;
        logmelParams->mfbUppFreqQ8 = (fract32)mfbUppFreqQ8Int;
        logmelParams->opt = 
          (melScaleOpt ? LOGMEL_OPT_MELSCALE_AUDITORY : LOGMEL_OPT_MELSCALE_HTK  ) |
          (fbeLogOpt   ? LOGMEL_OPT_FBELOG_BASE10     : LOGMEL_OPT_FBELOG_NATURAL) |
          (fbNormOpt   ? LOGMEL_OPT_FBNORM_AREA       : LOGMEL_OPT_FBNORM_NONE   );
        res = 1;
    }
    return res;
} /* loadParams() */

