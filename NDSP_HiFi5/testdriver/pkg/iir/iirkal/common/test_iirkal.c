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
 * Test procedures for Kalman Filter functions
 */

#include "test_iirkal.h"

int te_loadFxn_kalmanupd1( tTestEngContext * context )
{
    int N;
    int szXYZ, szU;
    int res=0;
    /* Retrieve test data dimensions. */
    N = context->args.dim[0];
    szXYZ = MAX(0, N); szU = 1;
    /* Allocate vectors for input data, results and lower/upper result thresholds. */
    if (6!=vecsAlloc(context->desc->isAligned, context->desc->fmt,
                     &context->dataSet.X, szXYZ,
                     &context->dataSet.Y, szXYZ,
                     &context->dataSet.U, szU,
                     &context->dataSet.Z, szXYZ,
                     &context->dataSet.Zlo, szXYZ,
                     &context->dataSet.Zhi, szXYZ, NULL)) {
        printf("te_loadFxn_kalmanupd1: failed to allocate vectors\n");
    /* Load input data and lower/upper result thresholds from the SEQ-file. */
    } else if (!seqFileReadVecs(context->seqFile, 
                                &context->dataSet.X, &context->dataSet.Y, &context->dataSet.U,
                                &context->dataSet.Zlo, &context->dataSet.Zhi, NULL)) {
        printf("te_loadFxn_kalmanupd1: failed to load vectors from the SEQ-file\n");
    } else {
        res = 1;
    }
    return res;
} /* te_loadFxn_kalmanupd1() */

void te_processFxn_kalmanupd1( tTestEngContext * context )
{
    typedef void tFxn_fxp(void * pScr, void * K, const void * U, const void * H, const void * R, int N, int qK, int qU, int qH, int qR);
    typedef void tFxn_flp(void * pScr, void * K, const void * U, const void * H, const void * R, int N);
    void *pK, *pU, *pH, *pR, *pScr;
    size_t szScr;
    int N, qK, qU, qH, qR;
    ASSERT( context && context->target.fut );
    /* Retrieven pointers to in/out arrays together with respective fixed point positions. */
    pH = vecGetElem( &context->dataSet.X, 0 ); qH = context->args.dim[1];
    pU = vecGetElem( &context->dataSet.Y, 0 ); qU = context->args.dim[2];
    pR = vecGetElem( &context->dataSet.U, 0 ); qR = context->args.dim[3];
    pK = vecGetElem( &context->dataSet.Z, 0 ); qK = context->args.dim[4];
    /* Test data dimension */
    N = context->args.dim[0];
    /* Allocate the scratch memory area, if required. */
    szScr = 0;
    if (context->desc->extraPtr) {
        const t_kalmanupd1_api * api = (t_kalmanupd1_api*)context->desc->extraPtr;
        szScr = api->getScratchSize(N);
    }
    pScr = szScr==0 ? NULL : mallocAlign(szScr, 1);
    /* Update the Test Coverage report. */
    te_vReportStd(context);
    /* Invoke the test target function */
    switch (context->desc->fmt & FMT_DTYPE_MASK) {
    case FMT_FRACT32: ((tFxn_fxp*)context->target.fut)(pScr, pK, pU, pH, pR, N, qK, qU, qH, qR); break;
    case FMT_FLOAT32: ((tFxn_flp*)context->target.fut)(pScr, pK, pU, pH, pR, N                ); break;
    }
    /* Cleanup */
    if (pScr) freeAlign(pScr);
} /* te_processFxn_kalmanupd1() */

