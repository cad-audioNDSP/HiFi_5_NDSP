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
 * Test-engine add-on for MFCC features extraction APIs.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* DSP Library API. */
#include LIBRARY_HEADER(audio)
#include LIBRARY_HEADER(fft) 
/* Test engine API. */
#include "testeng.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"
/* Test engine add-on for MFCC tests. */
#include "testeng_mfcc.h"
/* Aligning wrapper over standard malloc/free. */
#include "malloc16.h"
/* Fixed point arithmetics. */
#include "NatureDSP_Math.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Min/Max FFT size supported by real-to-complex FFT routines. */
#define TE_MFCC_MIN_FFT_SIZE           256
#define TE_MFCC_MAX_FFT_SIZE          2048
/* If non-zero, save cepstra results to a log file. To visualize test results, invoke the
 * following MATLAB script: NatureDSP_HiFi4_vectors/Matlab/testmfcc/analyse_cepsta_log.m */
#define TE_MFCC_CEPSTRA_LOG_ENABLE    0
#define TE_MFCC_CEPSTRA_LOG_FILENAME  "te_mfcc_cepstra.log"
/* Select either periodic (1) or symmetric (0) flavor of Hamming weighting window. */
#define TE_MFCC_USE_PERIODIC_HAMMING  1

static complex_float _makecomplexf( float32_t re, float32_t im )
{
  union { float32_t r[2]; complex_float c; } u;
  u.r[0]=re;
  u.r[1]=im;
  return (u.c);
}

/* Load the MFCC params structure from a SEQ-file. Return zero if failed. */
static int loadParams( mfcc_params_t * mfccParams, const te_mfcc_api_t * mfccApi, tSeqFile seqFile );

/* Compare computed MFCC features with reference results, and return the maximum absolute
 * difference over the frame. */
typedef float64_t comp_cepstra_mismatch_fxn_t(const void * cepstra, const float64_t * ref, int cepstraNum);
static float64_t comp_cepstra_mismatch32x32(const fract32 * cepstra, const float64_t * ref, int cepstraNum);
static float64_t comp_cepstra_mismatchf(const float32_t * cepstra, const float64_t * ref, int cepstraNum);

/* Print cepstra features to a log file. */
static void save_cepstra(FILE * fLog, const tVec * vecCepstra);

/* Generate Hamming weighting window, 32-bit fixed-point. */
static void hamming32x32_cbfxn( void * host, fract32 * w, int N );
/* Generate Hamming weighting window, single precision floating-point. */
static void hammingf_cbfxn( void * host, float32_t * w, int N );

/* Real-to-complex FFT, 32-bit fixed-point. */
static int rfft32x32_cbfxn( void * host, complex_fract32 * restrict y, fract32 * restrict x, int fftSize );
/* Real-to-complex FFT, single precision floating-point. */
static void rfftf_cbfxn( void * host, complex_float * restrict y, float32_t * restrict x, int fftSize );

/* Test host encapsulates various resources needed to run a whole SEQ-file. */
typedef struct te_mfcc_host_tag {
    mfcc_params_t mfccParams;
    union {
        mfcc32x32_callback_t _32x32;
        mfccf_callback_t     _f;
    } callback[2];
    comp_cepstra_mismatch_fxn_t * compCepstraMismatch;
    complex_float * twdTbl;     /* Precomputed twiddle factor table for the floating-point FFT. */
    int             twdSize;    /* Max FFT size allowed by the pre-computed table.              */
    FILE *          cepstraLog; /* File handle for logging cepstra results.                     */
} te_mfcc_host_t;

/* Prepare to run test cases. */
int te_createFxn_mfcc( tTestEngContext * context )
{
    const te_mfcc_api_t * mfccApi = (te_mfcc_api_t*)context->target.fut;
    te_mfcc_host_t * mfccHost;
    complex_float * twdTblf = NULL;
    /* Check if the MFCC API is implemented for the current build configuration. */
    if (!IS_PRESENT(mfccApi->getDefaultParams) || !IS_PRESENT(mfccApi->alloc) || !IS_PRESENT(mfccApi->init) || 
        !IS_PRESENT(mfccApi->process) || !IS_PRESENT(mfccApi->getScratchSize)) return -1;
    /* For a floating-point MFCC test, precompute the twiddle factor table to be used by the FFT callback. */
    if ((context->desc->fmt & FMT_DTYPE_MASK)==FMT_FLOAT32) {
        const int N = TE_MFCC_MAX_FFT_SIZE;
        int m, n;
        twdTblf = (complex_float*)mallocAlign(3*N/4*sizeof(complex_float), 0);
        if (NULL==twdTblf) {
            printf("te_createFxn_mfcc: failed to allocate the twiddle factor table\n"); return 0;
        }
        for ( m=0; m<3; m++ ) {
            for ( n=0; n<N/4; n++ ) {
                float64_t phi = -2*M_PI*(m+1)*n/N, cs = cos(phi), sn = sin(phi);
                twdTblf[n*3+m] = _makecomplexf((float32_t)cs, (float64_t)sn);
            } /* m */
        } /* n */
    }
    /* Setup the host structure. */
    if (NULL==(mfccHost = (te_mfcc_host_t*)mallocAlign(sizeof(*mfccHost), 0))) {
        printf("te_createFxn_mfcc: failed to allocate the host structure\n"); return 0;
    }
    memset(mfccHost, 0, sizeof(*mfccHost));
    switch (context->desc->fmt) {
    case FMT_REAL|FMT_FRACT32:
        /* useHamming == 0 */
        mfccHost->callback[0]._32x32.host = mfccHost;
        mfccHost->callback[0]._32x32.genWindow = NULL;
        mfccHost->callback[0]._32x32.rfft = &rfft32x32_cbfxn;
        /* useHamming == 1 */
        mfccHost->callback[1]._32x32.host = mfccHost;
        mfccHost->callback[1]._32x32.genWindow = &hamming32x32_cbfxn;
        mfccHost->callback[1]._32x32.rfft = &rfft32x32_cbfxn;
        mfccHost->compCepstraMismatch = (comp_cepstra_mismatch_fxn_t*)&comp_cepstra_mismatch32x32;
        break;
    case FMT_REAL|FMT_FLOAT32:
        /* useHamming == 0 */
        mfccHost->callback[0]._f.host = mfccHost;
        mfccHost->callback[0]._f.genWindow = NULL;
        mfccHost->callback[0]._f.rfft = &rfftf_cbfxn;
        /* useHamming == 1 */
        mfccHost->callback[1]._f.host = mfccHost;
        mfccHost->callback[1]._f.genWindow = &hammingf_cbfxn;
        mfccHost->callback[1]._f.rfft = &rfftf_cbfxn;
        mfccHost->twdTbl = twdTblf;
        mfccHost->twdSize = TE_MFCC_MAX_FFT_SIZE;
        mfccHost->compCepstraMismatch = (comp_cepstra_mismatch_fxn_t*)&comp_cepstra_mismatchf;
        break;
    default:
        NASSERT(!"Unknown data format");
        return 0;
    }
#if TE_MFCC_CEPSTRA_LOG_ENABLE
    mfccHost->cepstraLog = fopen(TE_MFCC_CEPSTRA_LOG_FILENAME, "wt");
    if (NULL==mfccHost->cepstraLog) {
        printf("te_createFxn_mfcc: failed to open %s for output\n", TE_MFCC_CEPSTRA_LOG_FILENAME); return 0;
    }
    fprintf(mfccHost->cepstraLog, "%s\n", context->seqFile->filename);
#endif
    context->target.handle = mfccHost;
    return 1;
} /* te_createFxn_mfcc() */

/* Release resources after the last test case, */
int te_destroyFxn_mfcc( tTestEngContext * context )
{
    te_mfcc_host_t * mfccHost = (te_mfcc_host_t*)context->target.handle;
    if (NULL==mfccHost) return 1;
    if (mfccHost->cepstraLog) fclose(mfccHost->cepstraLog);
    if (NULL!=mfccHost->twdTbl) freeAlign(mfccHost->twdTbl);
    freeAlign(mfccHost);
    return 1;
} /* te_destroyFxn_mfcc() */

/* Allocate vectors and load the data set for a test case. */
int te_loadFxn_mfcc( tTestEngContext * context )
{
    te_mfcc_host_t * mfccHost = (te_mfcc_host_t*)context->target.handle;
    const te_mfcc_api_t * mfccApi = (te_mfcc_api_t*)context->target.fut;
    mfcc_params_t * mfccParams;
    int testLen, res=0;
    NASSERT(NULL!=mfccHost);
    NASSERT(NULL!=mfccApi);
    testLen = context->args.dim[1];
    mfccParams = &mfccHost->mfccParams;
    /* Load the MFCC configuration. */
    if (0==loadParams(mfccParams, mfccApi, context->seqFile)) {
        printf("te_loadFxn_mfcc: failed to load the params structure\n");
    /* Allocate the test speech vector, reference MFCC features vector, mismatch and threshold vectors. */
    } else if (0==vecAlloc(&context->dataSet.X, testLen, 
                           TE_ALIGN_YES, context->desc->fmt, NULL) ||
               4!=vecsAlloc(TE_ALIGN_YES, FMT_REAL|FMT_FLOAT64, 
                            &context->dataSet.Y, testLen/mfccParams->stftHopLen*mfccParams->cepstraNum,
                            &context->dataSet.Z, testLen/mfccParams->stftHopLen,
                            &context->dataSet.Zlo, testLen/mfccParams->stftHopLen,
                            &context->dataSet.Zhi, testLen/mfccParams->stftHopLen, NULL)) {
        printf("te_loadFxn_mfcc: failed to allocate vectors\n");
    /* Load the test speech vector, reference MFCC features vector, and mismatch threshold vectors. */
    } else if (0==seqFileReadVecs(context->seqFile, &context->dataSet.X, &context->dataSet.Y, 
                                  &context->dataSet.Zlo, &context->dataSet.Zhi, NULL)) {
        printf("te_loadFxn_mfcc: failed to load vectors from the SEQ-file\n");
    } else {
        res = 1;
    }
    return res;
} /* te_loadFxn_mfcc() */

/* Apply the target function to the test case data set. */
void te_processFxn_mfcc( tTestEngContext * context )
{
    te_mfcc_host_t * mfccHost = (te_mfcc_host_t*)context->target.handle;
    const te_mfcc_api_t *mfccApi = (te_mfcc_api_t*)context->target.fut;
    mfcc_params_t * mfccParams;
    te_mfcc_handle_t mfcc;
    size_t szObj, szScr, szElem;
    tVec vecObj, vecScr, vecSpeech, vecCepstra;
    void *pObj, *pScr, *speech, *cepstra, *pX;
    float64_t *pY, *pZ;
    int useHamming, testLen, hopLen, hopNum, cepstraNum;
    int k, res=1;
    NASSERT(NULL!=mfccHost);
    NASSERT(NULL!=mfccApi);
    mfccParams = &mfccHost->mfccParams;
    szElem = vecElemSize(context->desc->fmt);
    /* Test data dimensions. */
    useHamming = (0!=context->args.dim[0]); 
    testLen = context->args.dim[1]; hopLen = mfccParams->stftHopLen; 
    hopNum = testLen/hopLen; NASSERT(hopNum*hopLen==testLen);
    cepstraNum = mfccParams->cepstraNum;
    /* If cepstra results logging is enabled, print a test case header. */
    if (mfccHost->cepstraLog) {
        fprintf(mfccHost->cepstraLog, "%d %d %d; test case # cepstraNum hopNum\n", 
                context->args.caseNum, cepstraNum, hopNum);
    }
    /* Allocate input/output vectors. */
    if (2!=vecsAlloc(context->desc->isAligned, context->desc->fmt, 
                     &vecSpeech, hopLen, &vecCepstra, cepstraNum, NULL)) {
        printf("te_processFxn_mfcc: failed to allocate vectors for in/out data"); return;
    }
    speech = vecGetElem(&vecSpeech, 0);
    cepstra = vecGetElem(&vecCepstra, 0);
    /* Allocate object and scratch memory areas, then instantiate the MFCC features extractor. */
    szObj = mfccApi->alloc(mfccParams);
    szScr = mfccApi->getScratchSize(mfccParams);
    if (!vecsAlloc(TE_ALIGN_NO, FMT_UINT8, &vecObj, szObj, &vecScr, szScr, NULL)) {
        printf("te_processFxn_mfcc: failed to allocate object and/or scratch memory\n"); return;
    }
    pObj = vecGetElem(&vecObj, 0); 
    pScr = vecGetElem(&vecScr, 0);
    if (NULL==(mfcc = mfccApi->init(pObj, mfccParams, &mfccHost->callback[useHamming]))) {
        printf("te_processFxn_mfcc: failed to initialize the MFCC features extractor\n"); return;
    }
    /* Add to the log. */
    {
        char str[200];
        tReportFUT fut[3];
        fut[0]=(tReportFUT)mfccApi->init;
        fut[1]=(tReportFUT)mfccApi->alloc;
        fut[2]=(tReportFUT)mfccApi->process;
        sprintf(str, "WeightingWindow: %s; Pre-emph: %s; Bands: %d; Ceps: %d; Lifter: %s; opt: 0x%02x", 
                useHamming ? "Yes" : "No",
                mfccParams->preemph != 0 ? "Yes" : "No",
                mfccParams->mfbBandNum, mfccParams->cepstraNum,
                mfccParams->lifter != 0 ? "Yes" : "No",
                (unsigned)mfccParams->opt);
        vReportAdd(fut,3,str,context->seqFile->filename ,context->args.caseType,vecGetSize(&context->dataSet.X));
    }
    /* Process the test data. For each chunk of stftHopLen input samples, compute cepstraNum MFCC
     * features. Each MFCC vector is compared with reference result, and peak value of the absolute
     * difference (i.e. the L-inf norm of the eror vector) is appended to the output vector Z. */
    pX = vecGetElem(&context->dataSet.X, 0); /* (In) Test speech */
    pY = vecGetElem_fl64(&context->dataSet.Y, 0); /* (Out) Reference MFCC features */
    pZ = vecGetElem_fl64(&context->dataSet.Z, 0); /* (Out) L1 norm of error vector */
    for ( k=0; k<hopNum; k++ ) {
        memcpy(speech, pX, hopLen*szElem);
        mfccApi->process(mfcc, pScr, cepstra, speech);
        *pZ++ = mfccHost->compCepstraMismatch(cepstra, pY, cepstraNum);
        if (mfccHost->cepstraLog) save_cepstra(mfccHost->cepstraLog, &vecCepstra);
        pX = (uint8_t*)pX + hopLen*szElem; pY += cepstraNum;
    } /* n */
    /* When in verbose mode, print the maximum difference between MFCC results and reference cepstra. */
    if (context->isVerbose) {
        float64_t peakErr = 0;
        pZ = vecGetElem_fl64(&context->dataSet.Z, 0);
        for ( k=0; k<hopNum; k++, pZ++ ) {
            if (!(isnan(peakErr) || peakErr>=*pZ)) peakErr = *pZ;
        }
        printf("Peak error: %.2e ", peakErr);
    }
    /* Cleanup */
    res &= (0!=vecsFree(&vecObj, &vecScr, &vecSpeech, &vecCepstra, NULL)); 
    /* If encountered any error during cleanup, deliberately distort the output vector so
     * that the test case is invalidated. */
    if (!res) {
        pZ = vecGetElem_fl64(&context->dataSet.Z, 0);
        for ( k=0; k<hopNum; k++ ) *pZ++ = FP_NAN;
    }
} /* te_processFxn_mfcc() */

/* Load the MFCC params structure from a SEQ-file. Return zero if failed. */
int loadParams( mfcc_params_t * mfccParams, const te_mfcc_api_t * mfccApi, tSeqFile seqFile )
{
    int preemphInt=0, mfbLowFreqQ8Int=0, mfbUppFreqQ8Int=0;
    int melScaleOpt=0, fbeLogOpt=0, fbNormOpt=0, dcOpt=0, preemphOpt=0, dctNormOpt=0, res=0;
    memset(mfccParams, 0, sizeof(*mfccParams));
    mfccApi->getDefaultParams(mfccParams);
    if (1==seqFileScanf(seqFile, "%d", &mfccParams->Fs        ) &&
        1==seqFileScanf(seqFile, "%d", &mfccParams->scaleExp  ) &&
        1==seqFileScanf(seqFile, "%d", &preemphInt            ) &&
        1==seqFileScanf(seqFile, "%d", &mfccParams->stftWinLen) &&
        1==seqFileScanf(seqFile, "%d", &mfccParams->stftHopLen) &&
        1==seqFileScanf(seqFile, "%d", &mfccParams->fftSize   ) &&
        1==seqFileScanf(seqFile, "%d", &mfbLowFreqQ8Int       ) &&
        1==seqFileScanf(seqFile, "%d", &mfbUppFreqQ8Int       ) &&
        1==seqFileScanf(seqFile, "%d", &mfccParams->mfbBandNum) &&
        1==seqFileScanf(seqFile, "%d", &mfccParams->cepstraNum) &&
        1==seqFileScanf(seqFile, "%d", &mfccParams->lifter    ) &&
        1==seqFileScanf(seqFile, "%d", &melScaleOpt           ) &&
        1==seqFileScanf(seqFile, "%d", &fbeLogOpt             ) &&
        1==seqFileScanf(seqFile, "%d", &fbNormOpt             ) &&
        1==seqFileScanf(seqFile, "%d", &dcOpt                 ) &&
        1==seqFileScanf(seqFile, "%d", &preemphOpt            ) &&
        1==seqFileScanf(seqFile, "%d", &dctNormOpt            ))
    {
        mfccParams->preemph = (fract16)preemphInt;
        mfccParams->mfbLowFreqQ8 = (fract32)mfbLowFreqQ8Int;
        mfccParams->mfbUppFreqQ8 = (fract32)mfbUppFreqQ8Int;
        mfccParams->opt = 
          (melScaleOpt ? LOGMEL_OPT_MELSCALE_AUDITORY : LOGMEL_OPT_MELSCALE_HTK      ) |
          (fbeLogOpt   ? LOGMEL_OPT_FBELOG_BASE10     : LOGMEL_OPT_FBELOG_NATURAL    ) |
          (fbNormOpt   ? LOGMEL_OPT_FBNORM_AREA       : LOGMEL_OPT_FBNORM_NONE       ) |
          (dcOpt       ? MFCC_OPT_DC_MEAN_DONT_REMOVE : MFCC_OPT_DC_MEAN_REMOVE      ) |
          (preemphOpt  ? MFCC_OPT_PREEMPH_CONTINUOUS  : MFCC_OPT_PREEMPH_FRAMEBYFRAME) |
          (dctNormOpt  ? MFCC_OPT_DCT_ORTHOGONAL      : MFCC_OPT_DCT_NORMALIZED      );
        res = 1;
    }
    return res;
} /* loadParams() */

/* Compare computed MFCC features (32-bit fixed-point) with reference results, and return the
 * maximum absolute difference over the frame. */
 float64_t comp_cepstra_mismatch32x32(const fract32 * cepstra, const float64_t * ref, int cepstraNum)
{
    float64_t diff, diff_Linf=0;
    int n;
    for ( n=0; n<cepstraNum; n++ ) {
        diff = fabs(ldexp((float64_t)cepstra[n], -MFCC_CEPSTRA_FRACT_BITS) - ref[n]);
        if (diff_Linf<diff) diff_Linf = diff;
    }
    return diff_Linf;
} /* comp_cepstra_mismatch32x32() */

/* Compare computed MFCC features (single precision floating-point) with reference results, and 
 * return the maximum absolute difference over the frame. */
float64_t comp_cepstra_mismatchf(const float32_t * cepstra, const float64_t * ref, int cepstraNum)
{
    float64_t diff, diff_Linf=0;
    int n;
    for ( n=0; n<cepstraNum; n++ ) {
        diff = fabs((float64_t)cepstra[n]-ref[n]);
        /* This condition clause should allow for NaN propagation. */
        if (!(isnan(diff_Linf) || diff_Linf>=diff)) diff_Linf = diff;
    }
    return diff_Linf;
} /* comp_cepstra_mismatchf() */ 

/* Print cepstra features to a log file. */
static void save_cepstra(FILE * fLog, const tVec * vecCepstra)
{
    float64_t * cepstra;
    int n, cepstraNum;
    NASSERT(fLog && vecCepstra);
    cepstraNum = vecCepstra->nElem;
    cepstra = (float64_t*)malloc(cepstraNum*sizeof(float64_t));
    if (NULL==cepstra) return;
    vecToFp64(cepstra, vecCepstra, MFCC_CEPSTRA_FRACT_BITS-31);
    for ( n=0; n<cepstraNum; n++ ) {
        fprintf(fLog, "%+.17e ", cepstra[n]);
    }
    fprintf(fLog, "\n");
    free(cepstra);
} /* save_cepstra() */

/* Generate Hamming weighting window, 32-bit fixed-point. */
void hamming32x32_cbfxn( void * host, fract32 * w, int N )
{
    float64_t f;
    int n;
    for ( n=0; n<N; n++ ) {
#if TE_MFCC_USE_PERIODIC_HAMMING
        /* MATLAB: hamming(N,'periodic') */
        f = 0.54-0.46*cos(2*M_PI*n/N);
#else
        /* MATLAB: hamming(N) or hamming(N,'symmetric') */
        f = 0.54-0.46*cos(2*M_PI*n/(N-1));
#endif
        f = round(ldexp(f, 31));
        if (f>2147483647.) f = 2147483647.;
        if (f<-2147483648.) f = -2147483648.;
        w[n] = (fract32)f;
    }
} /* hamming32x32_cbfxn() */ 

/* Generate Hamming weighting window, single precision floating-point. */
void hammingf_cbfxn( void * host, float32_t * w, int N )
{
    float64_t f;
    int n;
    for ( n=0; n<N; n++ ) {
#if TE_MFCC_USE_PERIODIC_HAMMING
        /* MATLAB: hamming(N,'periodic') */
        f = 0.54-0.46*cos(2*M_PI*n/N);
#else
        /* MATLAB: hamming(N) or hamming(N,'symmetric') */
        f = 0.54-0.46*cos(2*M_PI*n/(N-1));
#endif
        w[n] = (float32_t)f;
    }
} /* hammingf_cbfxn() */

/* Real-to-complex FFT, 32-bit fixed-point. */
int rfft32x32_cbfxn( void * host, complex_fract32 * restrict y, fract32 * restrict x, int fftSize )
{
    int hidx;
    fft_handle_t hfft[] = {rfft32_128, rfft32_256, rfft32_512, rfft32_1024, rfft32_2048};
    (void)host;
    NASSERT(0==(fftSize & (fftSize-1)));
    hidx = 30-S_exp0_l(fftSize)-7;
    NASSERT(0<=hidx && hidx<(int)sizeof(hfft)/sizeof(hfft[0]));
    return fft_real32x32((int32_t*)y, x, hfft[hidx], 2); /* scalingOpt == 2: dynamic scaling */
} /* rfft32x32_cbfxn() */

/* Real-to-complex FFT, single precision floating-point. */
void rfftf_cbfxn( void * host, complex_float * restrict y, float32_t * restrict x, int fftSize )
{
    te_mfcc_host_t * mfccHost = (te_mfcc_host_t*)host;
    int twdStep;
    NASSERT(NULL!=mfccHost && NULL!=mfccHost->twdTbl);
    NASSERT(0==(fftSize & (fftSize-1)));
    NASSERT(0<fftSize && fftSize<=mfccHost->twdSize);
#if 0 /* Enable to store FFT input to a text file, then run the following MATLAB script:
       * NatureDSP_HiFi4_vectors/Matlab/testmfcc/external_rfftf.m */
    {
        static FILE * f = NULL;
        int n;
        if (NULL==f) {
            f = fopen("rfft_x.txt","wt");
        }
        NASSERT(f);
        for ( n=0; n<fftSize; n++ ) {
            fprintf(f, "%+.17e ", x[n]);
        }
        fprintf(f, "\n");
    }
#endif
    twdStep = mfccHost->twdSize >> (30-S_exp0_l(fftSize));
    fft_realf_ie(y, x, mfccHost->twdTbl, twdStep, fftSize);
#if 0 /* Enable to load FFT results from the file generated by MATLAB. */
    {
        static FILE * f = NULL;
        int n, res;
        if (NULL==f) {
            f = fopen("rfft_y.txt","rt");
        }
        if (NULL!=f) {
            for ( n=0; n<fftSize/2+1; n++ ) {
                float32_t re, im;
                res = fscanf(f, "(%f,%f) ", &re, &im);
                NASSERT(2==res);
                y[n] = _makecomplexf(re, im);
            }
        }
    }
#endif
} /* rfftf_cbfxn() */

