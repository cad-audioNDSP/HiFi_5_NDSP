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
 * Test-engine add-on for old IIR categories 
 */

#include <string.h>
#include <stdlib.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test engine extension for IIR filters. */
#include "testeng_iir_old.h"
/* Memory allocator w/ alignment support. */
#include "malloc16.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

typedef struct
{
    int  M;           /* number of sections                */
    void* objMem;     /* allocated memory for filter       */
    void *handle;     /* filter handle                     */
    tVec delay;       /* vector containing dealy elements  */
    tVec scratch;     /* temporal data, same size as delay */
}
tIirContext;

/* function reads impulse response from file and creates IIR structure. returns 0 if failed */
int te_create_iir(tTestEngContext * context)
{
    int res;
    size_t szObj;
    tVec coef, gain, scale;
    tIirDescr *api;
    int M, L, coefLen = 0;
    tIirContext *iirContext;
    if (seqFileScanf(context->seqFile, "%d %d", &M, &L) != 2)
    {
        printf("bad SEQ-file format\n");
        return 0;
    }
    ASSERT(L == 1);
    iirContext = (tIirContext *)malloc(sizeof(tIirContext));
    context->target.handle = (void*)iirContext;
    if (iirContext == NULL) return 0;
    memset(iirContext, 0, sizeof(*iirContext));
    iirContext->M = M;
    switch (context->desc->extraParam & 3)
    {
    case TE_IIR_DF1:
    case TE_IIR_DF2:
    case TE_IIR_DF2T:
        coefLen = M * 5; break;
    default: ASSERT(0);
    }
    api = (tIirDescr *)context->target.fut;
    if (!IS_PRESENT(api->init) ||
        !IS_PRESENT(api->alloc) ||
        !IS_PRESENT(api->process) ||
        !IS_PRESENT(api->delay))
    {
        return -1;  // FUT is not defined
    }
    if (context->desc->extraParam & TE_IIR_FLOAT)
    {
        int fmt = context->desc->fmt & FMT_DTYPE_MASK;
        if (!vecAlloc(&coef, coefLen, 0, fmt | FMT_REAL, NULL)) return 0;
        if (!vecAlloc(&gain, 0, 0, FMT_FRACT16 | FMT_REAL, NULL)) return 0;
        if (!vecAlloc(&scale, 1, 0, FMT_FLOAT32 | FMT_REAL, NULL)) return 0;
        res = seqFileReadVec(context->seqFile, &coef) &
            seqFileReadVec(context->seqFile, &scale);
    }
    else
    {
        int fmt = context->desc->fmt & FMT_DTYPE_MASK;
        if (!vecAlloc(&coef, coefLen, 0, fmt | FMT_REAL, NULL)) return 0;
        if (!vecAlloc(&gain, M, 0, FMT_FRACT16 | FMT_REAL, NULL)) return 0;
        if (!vecAlloc(&scale, 1, 0, FMT_FRACT16 | FMT_REAL, NULL)) return 0;
        res = seqFileReadVec(context->seqFile, &coef) &
            seqFileReadVec(context->seqFile, &gain) &
            seqFileReadVec(context->seqFile, &scale);
    }

    switch (context->desc->extraParam & 3)
    {
    case TE_IIR_DF1:
    case TE_IIR_DF2:
    case TE_IIR_DF2T:
        szObj = api->alloc(M);
        iirContext->objMem = malloc(szObj);
        if (context->desc->extraParam & TE_IIR_FLOAT)
        {
            iirContext->handle = ((tIirFxnInitFloat*)api->init)(iirContext->objMem, M, vecGetElem(&coef, 0), (int16_t)*vecGetElem_fl32(&scale, 0));
        }
        else
        {
            iirContext->handle = api->init(iirContext->objMem, M, vecGetElem(&coef, 0), vecGetElem_fr16(&gain, 0), *vecGetElem_fr16(&scale, 0));
        }
        break;
    }
    //initialise delay here, for now specifically for fract16
    int delaySize = api->delay(iirContext->handle);
    int fmt = context->desc->fmt;
    vecAlloc(&(iirContext->delay), 2*delaySize, 1, fmt, NULL);
    vecAlloc(&(iirContext->scratch), delaySize, 1, fmt, NULL);
    int i = 0;
    

    if (fmt == (FMT_REAL | FMT_FRACT16))
    {
        int16_t * delay = (int16_t*)(iirContext->delay.ptr) + iirContext->delay.offset;
        for (i = 0; i < delaySize; i++)
        {
            delay[i] = -10;
        }
        for (i = delaySize; i < 2 * delaySize; i++)
        {
            delay[i] = 10;
        }
    }
    else if (fmt == (FMT_REAL | FMT_FRACT32))
    {
        int32_t * delay = (int32_t*)(iirContext->delay.ptr) + iirContext->delay.offset;
        for (i = 0; i < delaySize; i++)
        {
            delay[i] = -10;
        }
        for (i = delaySize; i < 2 * delaySize; i++)
        {
            delay[i] = 10;
        }
    }
    else if (fmt == (FMT_REAL | FMT_FLOAT32))
    {
        float32_t * delay = (float32_t*)(iirContext->delay.ptr) + iirContext->delay.offset;
        for (i = 0; i < delaySize; i++)
        {
            delay[i] = -3e-5;
        }
        for (i = delaySize; i < 2 * delaySize; i++)
        {
            delay[i] = 3e-5;
        }
    }
    else if (fmt == (FMT_CPLX | FMT_FLOAT32))
    {
#if defined(COMPILER_MSVC)
        complex_float * delay = (complex_float*)(iirContext->delay.ptr) + iirContext->delay.offset;
        for (i = 0; i < delaySize; i++)
        {
            delay[i].s.re = -3e-5;
            delay[i].s.im = -3e-5;
        }
        for (i = delaySize; i < 2 * delaySize; i++)
        {
            delay[i].s.re = 3e-5;
            delay[i].s.im = 3e-5;
        }
#else
        float32_t * delay = (float32_t*)(iirContext->delay.ptr) + iirContext->delay.offset*2;
        for (i = 0; i < delaySize; i++)
        {
            delay[2*i]   = -3e-5;
            delay[2*i+1] = -3e-5;
        }
        for (i = delaySize; i < 2 * delaySize; i++)
        {
            delay[2*i]   = 3e-5;
            delay[2*i+1] = 3e-5;
        }
#endif
    }
    vecsFree(&coef, &gain, &scale, NULL);
    return res;
}

/* function destroys IIR structure, returns 0 if failed */
int te_destroy_iir(tTestEngContext * context)
{
    tIirContext *iirContext;
    iirContext = (tIirContext *)context->target.handle;
    if (iirContext)
    {
        free(iirContext->objMem);
        if (iirContext->delay.nElem != 0)
        {
            vecFree(&iirContext->delay);
        }
        if (iirContext->scratch.nElem != 0)
        {
            vecFree(&iirContext->scratch);
        }
        free(iirContext);
    }
    return 1;
}

/* 
   Allocate vectors and load the data set for IIR:
*  vector X (in), vector Z (out) */
int te_loadFxn_iir(tTestEngContext * context)
{
    int type;
    int M, N, L;
    int nElemIn, nElemOut, res;
    tIirContext *iirContext;


    ASSERT(context && context->seqFile);

    iirContext = (tIirContext *)context->target.handle;
    M = context->args.dim[1];
    N = context->args.dim[0];
    L = 1;
    type = (context->desc->extraParam & TE_IIR_FLOAT) ? FMT_FLOAT32 : ((context->desc->extraParam & TE_IIR_16X16) ? FMT_FRACT16 : FMT_FRACT32);
    if (context->desc->fmt & FMT_CPLX) type |= FMT_CPLX;


    nElemIn = MAX(0, M*N*L);
    nElemIn = nElemOut = MAX(0, M*N*L);

    memset(&context->dataSet, 0, sizeof(context->dataSet));

    /* Allocate data vectors memory. */
    res = (4 == vecsAlloc(context->desc->isAligned, type,
        &context->dataSet.X, nElemIn,
        &context->dataSet.Z, nElemOut,
        &context->dataSet.Zlo, nElemOut,
        &context->dataSet.Zhi, nElemOut, 0));
    if (res)
    {
        /* Load vectors data from the SEQ-file. */
        if (!(res = seqFileReadVecs(context->seqFile,
            &context->dataSet.X,
            &context->dataSet.Zlo,
            &context->dataSet.Zhi, 0)))
        {
            printf("te_loadFxn_iir_old(): failed to read vectors data; "
                "fmt = 0x%02x, nElemIn = %d, nElemOut = %d\n",
                (unsigned)context->desc->fmt, nElemIn, nElemOut);
        }
    }
    else
    {
        printf("te_loadFxn_iir_old(): failed to allocate vectors; "
            "fmt = 0x%02x, nElemIn = %d, nElemOut = %d\n",
            (unsigned)context->desc->fmt, nElemIn, nElemOut);
    }

    //apply delay
    int sz = (iirContext->delay.nElem)/2;
    if (sz > 0)
    {
        int fmt = iirContext->delay.fmt;
        int n_ = context->dataSet.Zlo.nElem;
        void* delay = (void*)((uintptr_t)(iirContext->delay.ptr) + vecElemSize(fmt) * iirContext->delay.offset);
        void* zlo = (void*)((uintptr_t)(context->dataSet.Zlo.ptr) + vecElemSize(fmt) * context->dataSet.Zlo.offset);
        void* zhi = (void*)((uintptr_t)(context->dataSet.Zhi.ptr) + vecElemSize(fmt) * context->dataSet.Zhi.offset);
        void* tmp = (void*)((uintptr_t)(iirContext->scratch.ptr) + vecElemSize(fmt) * iirContext->scratch.offset);
        if (sz <= n_)
        {
            memmove(tmp, (void*)((uintptr_t)zlo + (n_ - sz)*vecElemSize(fmt)), sz*vecElemSize(fmt));
            memmove((void*)((uintptr_t)zlo + sz*vecElemSize(fmt)), zlo, (n_ - sz)*vecElemSize(fmt));
            memmove(zlo, delay, sz*vecElemSize(fmt));
            memmove(delay, tmp, sz*vecElemSize(fmt));

            memmove(tmp, (void*)((uintptr_t)zhi + (n_ - sz)*vecElemSize(fmt)), sz*vecElemSize(fmt));
            memmove((void*)((uintptr_t)zhi + sz*vecElemSize(fmt)), zhi, (n_ - sz)*vecElemSize(fmt));
            memmove(zhi, (void*)((uintptr_t)delay + sz*vecElemSize(fmt)), sz*vecElemSize(fmt));
            memmove((void*)((uintptr_t)delay + sz*vecElemSize(fmt)), tmp, sz*vecElemSize(fmt));

        }
        else
        {
            memmove(tmp, zlo, n_*vecElemSize(fmt));
            memmove(zlo, delay, n_*vecElemSize(fmt));
            memmove(delay, (void*)((uintptr_t)delay + n_*vecElemSize(fmt)), (sz - n_)*vecElemSize(fmt));
            memmove((void*)((uintptr_t)delay + (sz - n_)*vecElemSize(fmt)), tmp, n_*vecElemSize(fmt));

            memmove(tmp, zhi, n_*vecElemSize(fmt));
            memmove(zhi, (void*)((uintptr_t)delay + sz*vecElemSize(fmt)), n_*vecElemSize(fmt));
            memmove((void*)((uintptr_t)delay + sz*vecElemSize(fmt)), (void*)((uintptr_t)delay + (sz + n_)*vecElemSize(fmt)), (sz - n_)*vecElemSize(fmt));
            memmove((void*)((uintptr_t)delay + (sz + (sz - n_))*vecElemSize(fmt)), tmp, n_*vecElemSize(fmt));

        }
    }

    /* Free vectors data if failed. */
    if (!res) te_freeVectors(context);
    return (res);

} /* te_loadFxn_iir() */


/* Apply IIR function to the test case data set.
*  vector X (in), vector Z (out) */
void te_processFxn_iir(tTestEngContext * context)
{
    void *X, *Z;
    tIirContext* iirContext;
    tIirDescr *api;
    int N;

    ASSERT(context && context->target.fut);

    iirContext = (tIirContext *)context->target.handle;
    api = (tIirDescr *)context->target.fut;
    X = vecGetElem(&context->dataSet.X, 0);
    Z = vecGetElem(&context->dataSet.Z, 0);

    N = context->args.dim[0];
    /* call IIR function */
    if (context->desc->extraParam & TE_IIR_FLOAT)
    {
        { /* logging */
            tReportFUT fut[3];
            fut[0] = (tReportFUT)api->process;
            fut[1] = (tReportFUT)api->alloc;
            fut[2] = (tReportFUT)api->init;
            vReportAdd(fut, 3, "", context->seqFile->filename, context->args.caseType, te_vGetDataSize(context));
        }
        ((tIirFxnProcessFloat*)api->process)(iirContext->handle, Z, X, N);
    }
    else
    {
        void* pScr;
        { /* logging */
            tReportFUT fut[4];
            fut[0] = (tReportFUT)api->process;
            fut[1] = (tReportFUT)api->getScratch;
            fut[2] = (tReportFUT)api->alloc;
            fut[3] = (tReportFUT)api->init;
            vReportAdd(fut, 4, "", context->seqFile->filename, context->args.caseType, te_vGetDataSize(context));
        }
        /* allocate scratch */
        pScr = mallocAlign(api->getScratch(N, iirContext->M), 0);
        api->process(iirContext->handle, pScr, Z, X, N);
        /* free scratch */
        freeAlign(pScr);
    }
} /* te_processFxn_iir() */

/* function reads impulse response from file and creates IIR structure. returns 0 if failed */
int te_create_iir_stereo(tTestEngContext * context)
{
    int res;
    size_t szObj;
    tVec coefl, gainl, scalel, coefr, gainr, scaler;
    tIirStereoDescr *api;
    int M, L, coefLen = 0;
    tIirContext *iirContext;
    if (seqFileScanf(context->seqFile, "%d %d", &M, &L) != 2)
    {
        printf("bad SEQ-file format\n");
        return 0;
    }
    ASSERT(L == 1);
    iirContext = (tIirContext *)malloc(sizeof(tIirContext));
    context->target.handle = (void*)iirContext;
    if (iirContext == NULL) return 0;
    memset(iirContext, 0, sizeof(*iirContext));
    iirContext->M = M;
    switch (context->desc->extraParam & 3)
    {
    case TE_IIR_DF1:
    case TE_IIR_DF2:
    case TE_IIR_DF2T:
        coefLen = M * 5; break;
    default: ASSERT(0);
    }
    api = (tIirStereoDescr *)context->target.fut;
    if (!NatureDSP_Signal_isPresent(api->init) ||
        !NatureDSP_Signal_isPresent(api->alloc) ||
        !NatureDSP_Signal_isPresent(api->process))
    {
        return -1;  // FUT is not defined
    }
    if (context->desc->extraParam & TE_IIR_FLOAT)
    {
        if (!vecAlloc(&coefl, coefLen, 0, context->desc->fmt & ~(FMT_CPLX), NULL)) return 0;
        if (!vecAlloc(&coefr, coefLen, 0, context->desc->fmt & ~(FMT_CPLX), NULL)) return 0;
        if (!vecAlloc(&gainl, 0, 0, FMT_FRACT16, NULL)) return 0;
        if (!vecAlloc(&gainr, 0, 0, FMT_FRACT16, NULL)) return 0;
        if (!vecAlloc(&scalel, 1, 0, FMT_FLOAT32, NULL)) return 0;
        if (!vecAlloc(&scaler, 1, 0, FMT_FLOAT32, NULL)) return 0;
        res = seqFileReadVec(context->seqFile, &coefl) &
            seqFileReadVec(context->seqFile, &coefr) &
            seqFileReadVec(context->seqFile, &scalel) &
            seqFileReadVec(context->seqFile, &scaler);
    }
    else
    {
        if (!vecAlloc(&coefl, coefLen, 0, context->desc->fmt & ~(FMT_CPLX), NULL)) return 0;
        if (!vecAlloc(&coefr, coefLen, 0, context->desc->fmt & ~(FMT_CPLX), NULL)) return 0;
        if (!vecAlloc(&gainl, M, 0, FMT_FRACT16, NULL)) return 0;
        if (!vecAlloc(&gainr, M, 0, FMT_FRACT16, NULL)) return 0;
        if (!vecAlloc(&scalel, 1, 0, FMT_FRACT16, NULL)) return 0;
        if (!vecAlloc(&scaler, 1, 0, FMT_FRACT16, NULL)) return 0;
        res = seqFileReadVec(context->seqFile, &coefl) &
            seqFileReadVec(context->seqFile, &coefr) &
            seqFileReadVec(context->seqFile, &gainl) &
            seqFileReadVec(context->seqFile, &gainr) &
            seqFileReadVec(context->seqFile, &scalel) &
            seqFileReadVec(context->seqFile, &scaler);
    }

    switch (context->desc->extraParam & 3)
    {
    case TE_IIR_DF1:
    case TE_IIR_DF2:
    case TE_IIR_DF2T:
        szObj = api->alloc(M);
        iirContext->objMem = malloc(szObj);
        if (context->desc->extraParam & TE_IIR_FLOAT)
        {
            iirContext->handle = ((tIirStereoFxnInitFloat*)api->init)(iirContext->objMem, M, vecGetElem(&coefl, 0), (int16_t)*vecGetElem_fl32(&scalel, 0), vecGetElem(&coefr, 0), (int16_t)*vecGetElem_fl32(&scaler, 0));
        }
        else
        {
            iirContext->handle = api->init(iirContext->objMem, M, vecGetElem(&coefl, 0), vecGetElem_fr16(&gainl, 0), *vecGetElem_fr16(&scalel, 0), vecGetElem(&coefr, 0), vecGetElem_fr16(&gainr, 0), *vecGetElem_fr16(&scaler, 0));
        }
        break;
    }

    //initialise delay here, for now specifically for fract16
    int delaySize = 2*api->delay(iirContext->handle);
    int fmt = context->desc->fmt;
    vecAlloc(&(iirContext->delay), 2 * delaySize, 1, fmt, NULL);
    vecAlloc(&(iirContext->scratch), delaySize, 1, fmt, NULL);
    int i = 0;


    if (fmt == (FMT_REAL | FMT_FRACT16))
    {
        int16_t * delay = (int16_t*)(iirContext->delay.ptr) + iirContext->delay.offset;
        for (i = 0; i < delaySize; i++)
        {
            delay[i] = -10;
        }
        for (i = delaySize; i < 2 * delaySize; i++)
        {
            delay[i] = 10;
        }
    }
    else if (fmt == (FMT_REAL | FMT_FRACT32))
    {
        int32_t * delay = (int32_t*)(iirContext->delay.ptr) + iirContext->delay.offset;
        for (i = 0; i < delaySize; i++)
        {
            delay[i] = -10;
        }
        for (i = delaySize; i < 2 * delaySize; i++)
        {
            delay[i] = 10;
        }
    }
    else if (fmt == (FMT_REAL | FMT_FLOAT32))
    {
        float32_t * delay = (float32_t*)(iirContext->delay.ptr) + iirContext->delay.offset;
        for (i = 0; i < delaySize; i++)
        {
            delay[i] = -3e-5;
        }
        for (i = delaySize; i < 2 * delaySize; i++)
        {
            delay[i] = 3e-5;
        }
    }
    else if (fmt == (FMT_CPLX | FMT_FLOAT32))
    {
#if defined(COMPILER_MSVC)
        complex_float * delay = (complex_float*)(iirContext->delay.ptr) + iirContext->delay.offset;
        for (i = 0; i < delaySize; i++)
        {
            delay[i].s.re = -3e-5;
            delay[i].s.im = -3e-5;
        }
        for (i = delaySize; i < 2 * delaySize; i++)
        {
            delay[i].s.re = 3e-5;
            delay[i].s.im = 3e-5;
        }
#else
        float32_t * delay = (float32_t*)(iirContext->delay.ptr) + iirContext->delay.offset * 2;
        for (i = 0; i < delaySize; i++)
        {
            delay[2 * i] = -3e-5;
            delay[2 * i + 1] = -3e-5;
        }
        for (i = delaySize; i < 2 * delaySize; i++)
        {
            delay[2 * i] = 3e-5;
            delay[2 * i + 1] = 3e-5;
        }
#endif
    }

    vecsFree(&coefl, &gainl, &scalel, &coefr, &gainr, &scaler, NULL);
    return res;
} /* te_create_iir_stereo() */

/* function destroys IIR structure, returns 0 if failed */
int te_destroy_iir_stereo(tTestEngContext * context)
{
    tIirContext *iirContext;
    iirContext = (tIirContext *)context->target.handle;
    if (iirContext)
    {
        free(iirContext->objMem);
        if (iirContext->delay.nElem != 0)
        {
            vecFree(&iirContext->delay);
        }
        if (iirContext->scratch.nElem != 0)
        {
            vecFree(&iirContext->scratch);
        }
        free(iirContext);
    }
    return 1;
} /* te_destroy_iir_stereo() */

/* 
   Allocate vectors and load the data set for IIR:
*  vector X (in), vector Z (out) */
int te_loadFxn_iir_stereo(tTestEngContext * context)
{
    int type;
    int M, N, L;
    int nElemIn, nElemOut, res;
    tIirContext *iirContext;    

    ASSERT(context && context->seqFile);

    iirContext = (tIirContext *)context->target.handle;

    M = context->args.dim[1];
    N = context->args.dim[0];
    L = 1;
    type = (context->desc->extraParam & TE_IIR_FLOAT) ? FMT_FLOAT32 : ((context->desc->extraParam & TE_IIR_16X16) ? FMT_FRACT16 : FMT_FRACT32);
    if (context->desc->fmt & FMT_CPLX) type |= FMT_CPLX;

    nElemIn = MAX(0, 2 * M*N*L);
    nElemIn = nElemOut = MAX(0, 2 * M*N*L);

    memset(&context->dataSet, 0, sizeof(context->dataSet));

    /* Allocate data vectors memory. */
    res = (4 == vecsAlloc(context->desc->isAligned, type,
        &context->dataSet.X, nElemIn,
        &context->dataSet.Z, nElemOut,
        &context->dataSet.Zlo, nElemOut,
        &context->dataSet.Zhi, nElemOut, 0));
    if (res)
    {
        /* Load vectors data from the SEQ-file. */
        if (!(res = seqFileReadVecs(context->seqFile,
            &context->dataSet.X,
            &context->dataSet.Zlo,
            &context->dataSet.Zhi, 0)))
        {
            printf("te_loadFxn_iir_old(): failed to read vectors data; "
                "fmt = 0x%02x, nElemIn = %d, nElemOut = %d\n",
                (unsigned)context->desc->fmt, nElemIn, nElemOut);
        }
    }
    else
    {
        printf("te_loadFxn_iir_old(): failed to allocate vectors; "
            "fmt = 0x%02x, nElemIn = %d, nElemOut = %d\n",
            (unsigned)context->desc->fmt, nElemIn, nElemOut);
    }

    //apply delay
    int sz = (iirContext->delay.nElem) / 2;
    if (sz > 0)
    {
        int fmt = iirContext->delay.fmt;
        int n_ = context->dataSet.Zlo.nElem;
        void* delay = (void*)((uintptr_t)(iirContext->delay.ptr) + vecElemSize(fmt) * iirContext->delay.offset);
        void* zlo = (void*)((uintptr_t)(context->dataSet.Zlo.ptr) + vecElemSize(fmt) * context->dataSet.Zlo.offset);
        void* zhi = (void*)((uintptr_t)(context->dataSet.Zhi.ptr) + vecElemSize(fmt) * context->dataSet.Zhi.offset);
        void* tmp = (void*)((uintptr_t)(iirContext->scratch.ptr) + vecElemSize(fmt) * iirContext->scratch.offset);
        if (sz <= n_)
        {
            memmove(tmp, (void*)((uintptr_t)zlo + (n_ - sz)*vecElemSize(fmt)), sz*vecElemSize(fmt));
            memmove((void*)((uintptr_t)zlo + sz*vecElemSize(fmt)), zlo, (n_ - sz)*vecElemSize(fmt));
            memmove(zlo, delay, sz*vecElemSize(fmt));
            memmove(delay, tmp, sz*vecElemSize(fmt));

            memmove(tmp, (void*)((uintptr_t)zhi + (n_ - sz)*vecElemSize(fmt)), sz*vecElemSize(fmt));
            memmove((void*)((uintptr_t)zhi + sz*vecElemSize(fmt)), zhi, (n_ - sz)*vecElemSize(fmt));
            memmove(zhi, (void*)((uintptr_t)delay + sz*vecElemSize(fmt)), sz*vecElemSize(fmt));
            memmove((void*)((uintptr_t)delay + sz*vecElemSize(fmt)), tmp, sz*vecElemSize(fmt));

        }
        else
        {
            memmove(tmp, zlo, n_*vecElemSize(fmt));
            memmove(zlo, delay, n_*vecElemSize(fmt));
            memmove(delay, (void*)((uintptr_t)delay + n_*vecElemSize(fmt)), (sz - n_)*vecElemSize(fmt));
            memmove((void*)((uintptr_t)delay + (sz - n_)*vecElemSize(fmt)), tmp, n_*vecElemSize(fmt));

            memmove(tmp, zhi, n_*vecElemSize(fmt));
            memmove(zhi, (void*)((uintptr_t)delay + sz*vecElemSize(fmt)), n_*vecElemSize(fmt));
            memmove((void*)((uintptr_t)delay + sz*vecElemSize(fmt)), (void*)((uintptr_t)delay + (sz + n_)*vecElemSize(fmt)), (sz - n_)*vecElemSize(fmt));
            memmove((void*)((uintptr_t)delay + (sz + (sz - n_))*vecElemSize(fmt)), tmp, n_*vecElemSize(fmt));

        }
    }

    /* Free vectors data if failed. */
    if (!res) te_freeVectors(context);
    return (res);
} /* te_loadFxn_iir_stereo() */

/* Apply IIR function to the test case data set.
*  vector X (in), vector Z (out) */
void te_processFxn_iir_stereo(tTestEngContext * context)
{
    void *X, *Z;
    tIirContext* iirContext;
    tIirStereoDescr *api;
    int N;

    ASSERT(context && context->target.fut);

    iirContext = (tIirContext *)context->target.handle;
    api = (tIirStereoDescr *)context->target.fut;
    X = vecGetElem(&context->dataSet.X, 0);
    Z = vecGetElem(&context->dataSet.Z, 0);

    N = context->args.dim[0];
    /* call IIR function */
    if (context->desc->extraParam & TE_IIR_FLOAT)
    {
        { /* logging */
            tReportFUT fut[3];
            fut[0] = (tReportFUT)api->process;
            fut[1] = (tReportFUT)api->alloc;
            fut[2] = (tReportFUT)api->init;
            vReportAdd(fut, 3, "", context->seqFile->filename, context->args.caseType, te_vGetDataSize(context));
        }
        ((tIirStereoFxnProcessFloat*)api->process)(iirContext->handle, Z, X, N);
    }
    else
    {
        void* pScr;
        { /* logging */
            tReportFUT fut[4];
            fut[0] = (tReportFUT)api->process;
            fut[1] = (tReportFUT)api->getScratch;
            fut[2] = (tReportFUT)api->alloc;
            fut[3] = (tReportFUT)api->init;
            vReportAdd(fut, 4, "", context->seqFile->filename, context->args.caseType, te_vGetDataSize(context));
        }
        /* allocate scratch */
        pScr = mallocAlign(api->getScratch(N, iirContext->M), 0);
        api->process(iirContext->handle, pScr, Z, X, N);
        /* free scratch */
        freeAlign(pScr);
    }
} /* te_processFxn_iir_stereo() */

