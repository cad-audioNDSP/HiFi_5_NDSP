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
 * Test procedures for image processing APIs.
 */

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API. */
#include LIBRARY_HEADER(img)
/* Test engine API. */
#include "testeng.h"
#include "../../../fft/common/testeng_fft.h"
#include <stdlib.h>
#include <string.h>

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

/***************************** imgfft tests*******************************************/
typedef int tTestEngFrameProcFxn_rfft2d_stnd(tTestEngContext * context,
    const int16_t         * in,     /* 16-bit fixed point complex     */
    float64_t       * out,          /* 64-bit floating point complex  */
    tVec * xVec, tVec * yVec,	      /* In/out vectors for target FFT  */
    tVec * sVec);                   /* Scratch vector for target FFT  */

tTestEngFrameProcFxn_rfft2d_stnd te_frameProc_imgfft;

/* FFT test context; is accessible through tTestEngContext::target::handle */
typedef struct tagTestEngContext_imgfft
{
    char fInName[64]; /* Input data filename                                    */
    char fRefName[64]; /* Reference data filename                               */
    int  scale_method; /* Scale method (for feature rich fixed point FFTs only) */
    int  frameCnt;     /* Frame counter                                         */
    size_t dataSize;   /* data size in bytes allocated for input/output         */

} tTestEngContext_imgfft;

/* FFT test context. */
typedef struct tagTestEngContext_fft_int_imgfft
{
    tTestEngContext_imgfft ext; /* Externally visible part of FFT context. */

    struct {                    /* Twiddle factor tables.                  */
        int    baseSize;        /* First twiddle table FFT size.           */
        int    tblNum;          /* Number of twiddle tables.               */
        tVec * tblSet;          /* Twiddle tables, with FFT size doubled   */
    }                           /* in succession.                          */
    twd;

} tTestEngContext_fft_int_imgfft;

#include <math.h> // HUGE_VALF

/* tTestEngDesc::extraParam options for FFT tests. */
#define TE_FFT_REAL             0x0001 /* Forward/inverse real-valued FFT                                               */
#define TE_FFT_CPLX             0x0002 /* Forward/inverse complex-valued FFT                                            */
#define TE_FFT_BLOCKWISE        0x0004 /* Blockwise FFT                                                                 */
#define TE_FFT_FULL_SPECTRUM    0x0008 /* Real FFT forming full symmetric spectrum                                      */
#define TE_FFT_OPT_SCALE_METH   0x0010 /* Scale method support (fixed point only)                                       */
#define TE_FFT_OPT_INPLACE      0x0020 /* Input/outbut buffers may be aliased                                           */
#define TE_FFT_OPT_REUSE_INBUF  0x0040 /* FFT reuses input buffer for intermediate data                                 */
#define TE_FFT_NORM             0x0080 /* FFT input data to be normalized                                               */
#define TE_FFT_ISMIXED          0x0100 /* Mixed radix FFT                                                               */  
#define TE_FFT_STREAMWISE       0x0200 /* Streamwise FFT                                                                */
#define TE_FFT_AS               0x0400 /* Block/stream FFT with autoscaling                                             */
#define TE_FFT_NOSCALING        0x0800 /* Block/stream FFT without scaling                                              */
#define TE_FFT_IE_EXT_TWD       0x1000 /* Improved memory Efficiency FFT with externally supplied twiddle factors table */

/* Test case types related to FFT. */
#define TE_FFT_TESTCASE_FORWARD    6
#define TE_FFT_TESTCASE_INVERSE    7
/* Number of dimensional arguments of a test case (M,N,etc.) */
#define TE_ARGNUM_1    1
#define TE_ARGNUM_2    2
#define TE_ARGNUM_3    3
#define TE_ARGNUM_4    4

#define TEST_DESC_BRFFT_NOSCALING( fmt,fwdfft,invfft )   {{ (fmt),TE_FFT_REAL|TE_FFT_BLOCKWISE|TE_FFT_NOSCALING, NULL, TE_ARGNUM_2,TE_ALIGN_YES,   \
    &te_create_fft, &te_destroy_fft, &te_load_fft, &fileProc_imgfft }, (tFrwTransFxn)fwdfft, (tInvTransFxn)invfft}

/* Convert 64-bit floating point values to vector data format.
Return scaling power.
*/
static int vecPartFract16_FromFp64(tVec * vec, float64_t * x, int begin, int lenPart)
{
    int n, N;
    fract16 * z = (fract16 *)vecGetElem(vec, 0);
    int i, Exp;
    double maxX, minX;

    begin *= 1 + ((vec->fmt & FMT_CPLX) != 0); 

    N = ((vec->fmt & FMT_CPLX) ? 2 * lenPart : lenPart);

    ASSERT(vec && x);
    NASSERT((vec->fmt & FMT_DTYPE_MASK) == FMT_FRACT16);

    minX = maxX = x[0];

    for (i = 0; i < N; i++)
    {
        maxX = MAX(maxX, x[i]);
        minX = MIN(minX, x[i]);
    }
    maxX = MAX(abs(maxX), abs(minX));
    Exp = ceil(log(maxX) / log(2));
    if (Exp < 0)
    {
        printf("aa"); 
    }
    ASSERT(Exp >= 0);

    for (n = 0; n<N; n++)
    {
        float64_t t = ldexp(x[n], 15 - Exp);
        z[n + begin] = (fract16)MAX(MIN_INT16, MIN(MAX_INT16, round(t)));
    }
    return Exp;
}
#if 0
static void vecUint8ToFp64(float64_t * z, const tVec * vec, int shift)
{
    int n, N;

    ASSERT(z && vec);
    N = ((vec->fmt & FMT_CPLX) ? 2 * vec->nElem : vec->nElem);

    {
        uint8_t * x = (uint8_t *)vecGetElem(vec, 0);
        for (n = 0; n<N; n++) z[n] = ldexp(x[n], -shift - 15);
    }
}
#endif
/*
    Convert part of vector to double 
*/
static void vecPartToFp64(float64_t * z, 
                          const tVec * vec, 
                          int shift, 
                          int begin, /* index of the first element */
                          int Num    /* number elements to be converted */)
{
    int n, N;

    ASSERT(z && vec);
    N = ((vec->fmt & FMT_CPLX) ? 2 * Num : Num);
    begin = ((vec->fmt & FMT_CPLX) ? 2 * begin : begin);

    if ((vec->fmt & FMT_DTYPE_MASK)==FMT_UINT8)
    {
        uint8_t * x = (uint8_t *)vecGetElem(vec, 0);
        for (n = 0; n<N; n++) z[n] = ldexp(x[n+begin], -shift - 15);
    }
    else if ((vec->fmt & FMT_DTYPE_MASK) == FMT_INT8)
    {
      int8_t * x = (int8_t *)vecGetElem(vec, 0);
      for (n = 0; n<N; n++)
      {
        uint8_t tmp;
        tmp = x[n + begin] ^ 128;
        z[n] = ldexp((tmp), -shift - 15);
      }
    }
    else if ((vec->fmt & FMT_DTYPE_MASK) == FMT_FRACT16)
    {
        fract16 * x = (fract16 *)vecGetElem(vec, 0);
        for (n = 0; n<N; n++) z[n] = ldexp(x[n+begin], -shift - 15);
    }
    else
    {
        NASSERT(!"vecPartToFp64 unsupported vector format"); 
    }
}

/* Apply the FFT function to a single frame of test data, fft2d_real, ifft2d_real routines. */
int te_frameProc_imgfft(tTestEngContext * context,
    const int16_t         * in,
    float64_t       * out,
    tVec * xVec, tVec * yVec, tVec *sVec)
{
    tTestEngContext_imgfft* context_fft;
    typedef int(*tFxn)(void *scr, void * y, void * x, int mean, const imgsize_t* sz, ...);

    tTestEngTarget fxn = NULL;
    void *px, *py, *pscr;
    int Exp = 0;
    int  shift2;
    int32_t meanQ15_16=0, meanQ31_0=0;
    int M, N;
    int isRealForward, isRealInverse;
    int logMN;
    uint32_t paramFFT;

    int32_t histValues[1]; 
    int32_t totalScaling; 
    imghist_t h = {0, 0, 1, histValues};
    imgsize_t sz; 

    NASSERT(context && context->desc && context->target.handle);
    NASSERT(in && out && xVec && yVec);
    NASSERT(0 != (context->desc->extraParam & (TE_FFT_BLOCKWISE | TE_FFT_STREAMWISE)));

    context_fft = (tTestEngContext_imgfft *)context->target.handle;
    paramFFT = context->desc->extraParam;
    N = context->args.dim[0];
    M = context->args.dim[1];

    sz.width = N; 
    sz.height = M; 
    sz.stride = N;

    isRealForward = (paramFFT & TE_FFT_REAL) &&
        (context->args.caseType == TE_FFT_TESTCASE_FORWARD);
    isRealInverse = (paramFFT & TE_FFT_REAL) &&
        (context->args.caseType == TE_FFT_TESTCASE_INVERSE);
    if (isRealInverse)
    {
        logMN = (int)(0.5 + log(M*N) / log(2));
    }
    else
    {
        logMN = 0;
    }
    /* Convert 16-bit PCM input data to target FFT format. */
    if (context->desc->fmt == FMT_FRACT16 || isRealInverse != 0)
    {
        int i;
        int16_t *px = (int16_t *)vecGetElem(xVec, 0);
        if (isRealForward)
        {
            //mean = 1 << 14;
            for (i = 0; i < M*N; i++)
            {
                px[i] = in[i];
            }
            imghist_gs16(&h, (const void*)px, &sz, 1); 
            meanQ15_16 = h.mean;
            meanQ31_0 = (h.mean + 0x8000L) >> 16;  // Q15.16 -> Q31.0
        }
        else if (isRealInverse)
        {
            double *pref = (double*)in;
            double X0 = pref[0];
            //mean = (context->desc->fmt == FMT_UINT8) ? 1 << 7 : 1 << 14;

            meanQ31_0 = (int)(0.5 + X0 / (M*N));
            meanQ15_16 = (int)(0.5 + (1<<16) * X0 / (M*N));

            if ((context->desc->fmt == FMT_UINT8) || (context->desc->fmt == FMT_INT8))
            {
                meanQ31_0 >>= 7;
            }
#if 0
            pref[0] = 0;
            pref[1] = 0;
            Exp = -15 + vecPartFract16_FromFp64(xVec, pref, 0, xVec->nElem);

            pref[0] = X0;
#else
            Exp = (int32_t) pref[1]; 
#endif
        }
        else
        {
            ASSERT("ERROR:: te_frameProc_stnd_blkfft2");
        }
    }
    else if (context->desc->fmt == FMT_UINT8 && isRealForward != 0)
    {
       // meanQ31_0 = 1 << 7;

        int i;
        uint8_t *px = (uint8_t *)vecGetElem(xVec, 0);
        for (i = 0; i < M*N*(1 + isRealInverse); i++)
        {
            px[i] = (uint8_t)(in[i] >> 7);
        }
        imghist_gu8(&h, (const void*)px, &sz, 1);
        meanQ15_16 = h.mean;
        meanQ31_0 = (h.mean + 0x8000L) >> 16;  // Q15.16 -> Q31.0

    }
    else if (context->desc->fmt == FMT_INT8 && isRealForward != 0)
    {
      // meanQ31_0 = 1 << 7;
      int i;
      int8_t *px = (int8_t *)vecGetElem(xVec, 0);
      for (i = 0; i < M*N*(1 + isRealInverse); i++)
      {
        px[i] = (int8_t)((in[i] >> 7)^128);
      }
      imghist_gs8(&h, (const void*)px, &sz, 1);
      meanQ15_16 = h.mean;
      meanQ31_0 = (h.mean + 0x8000L) >> 16;  // Q15.16 -> Q31.0
    }

    /* Select in/out buffers for the FFT, and wipe the output buffer. */
    px = vecGetElem(xVec, 0);
    py = vecGetElem(yVec, 0);
    pscr = vecGetElem(sVec, 0);
    memset(py, 0, vecGetSize(yVec));

    /* Select the target FFT routine (either forward or inverse). */
    if (context->args.caseType == TE_FFT_TESTCASE_FORWARD)
    {
        fxn = (tTestEngTarget)((const tTestEngDesc_fft *)context->desc)->frwTransFxn;
    }
    else if (context->args.caseType == TE_FFT_TESTCASE_INVERSE)
    {
        fxn = (tTestEngTarget)((const tTestEngDesc_fft *)context->desc)->invTransFxn;
    }
    else
    {
        NASSERT(!"Bad test case type!");
    }

    if (((tTestEngContext_fft_int_imgfft *)context->target.handle)->twd.baseSize == 0)
    {
        /* Some FFTs functions not use external twiddle tables, in this case baseSize
        (first twiddle table fft size) set to 0 in the *.seq file  */
    }
    else
    {
        NASSERT(!"NO external twiddles for imgfft!"); 
    }

    { /* logging */ //shift2 = ((tFxn)fxn)(NULL, py, px, M, N);
        char info[20];
        tReportFUT fut[1];
        fut[0] = (tReportFUT)fxn;
        sprintf(info, "N=%d", N);
        vReportAdd(fut, 1, info, context_fft->fInName, context->args.caseType, vecGetSize(xVec) + vecGetSize(yVec));
    }

    shift2 = ((tFxn)fxn)(pscr, py, px, meanQ15_16, &sz, Exp);

    if (((context->desc->fmt == FMT_UINT8) || (context->desc->fmt == FMT_INT8) )&& context->args.caseType == TE_FFT_TESTCASE_INVERSE)
    {
        // vecUint8ToFp64(out, yVec, -15-7/*-15 - shift2 - Exp + logMN*/);
        totalScaling = -15 - 7; 
    }
    else
    {
        if (context->args.caseType == TE_FFT_TESTCASE_FORWARD)
        {
            double dcBias = (context->desc->fmt == FMT_UINT8) ? (128.0 * meanQ31_0*M * N) : 
              (context->desc->fmt == FMT_INT8) ? (128.0 * (meanQ31_0 + 128) *M * N) : (double)meanQ31_0*M*N;
            // vecToFp64((float64_t*)out, yVec, -15 - shift2 - Exp + logMN);
            out[0] += dcBias;
            totalScaling = -15 - shift2 - Exp + logMN;
        }
        else
        {
            // vecToFp64((float64_t*)out, yVec, -15);
            totalScaling = -15;
        }
    }

    * vecGetElem_i32(&context->dataSet.Y, 0) = totalScaling;
    * vecGetElem_i32(&context->dataSet.Y, 1) = meanQ31_0;

    return (1);
} /* te_frameProc_stnd_blkfft() */

/* Test data file processing for a standalone test. Applies the target FFT function
* to test data, compares the FFT output with reference data and estimates the SINAD. */
static void fileProc_imgfft(tTestEngContext * context)
{
    tTestEngContext_imgfft           * context_fft;
    tTestEngFrameProcFxn_rfft2d_stnd * frameProcFxn;

    char fInNameBuf[256], fRefNameBuf[256];
#define __MAX_DIM 640
    tVec inVec, outVec, refVec; /* In/out/reference vectors related to test data files. */
    tVec xVec, yVec;            /* In/out vectors for the target FFT routine. */
    tVec sVec;                  /* Scratch vector */
    FILE *fIn = 0, *fRef = 0;
    int lenIn = 0, lenOut = 0;
    int lenInBuffer = 0;       /* Size of buffer for input samples */
    int fmtX = 0, fmtY = 0;
    int res = 0;
    int n, N, M;
    int isForward;
    int blockExp[__MAX_DIM];
    int maxExp = -1000;

    NASSERT(context && context->desc);
    NASSERT(context->target.fut && context->target.handle);

    context_fft = (tTestEngContext_imgfft*)context->target.handle;

    /* FFT size */
    N = context->args.dim[0];
    M = context->args.dim[1];

    NASSERT(__MAX_DIM >= MAX(N, M)); 

    if (!vecAlloc(&context->dataSet.Y, 2, TE_ALIGN_YES, FMT_INT32, 0))
    {
        printf("fileProc_stnd(): failed to allocate dataSet.Y\n");
    }
    else
        res = 1;

    /* Is forward transform under test? Or inverse? */
    isForward = (context->args.caseType == TE_FFT_TESTCASE_FORWARD);

    /* Select length and format (real/complex) for all vectors. */
    lenIn = lenOut = M*N;
    lenInBuffer = (isForward) ? lenIn : N; 

    fmtX = (isForward ? FMT_REAL : FMT_CPLX);
    fmtY = (isForward ? FMT_CPLX : FMT_REAL);

    /* Preset an invalid SINAD for the test result. */
    *vecGetElem_fl32(&context->dataSet.Z, 0) = -HUGE_VALF;

    if (context->isVerbose)
    {
        printf("%-40s  ", context_fft->fRefName);
    }

    NASSERT(strlen(context->seqDir) + 8 + 1 + strlen(context_fft->fInName) < sizeof(fInNameBuf));
    sprintf(fInNameBuf, "%s/imgfft1/%s", context->seqDir, context_fft->fInName);
    NASSERT(strlen(context->seqDir) + 8 + 1 + strlen(context_fft->fRefName) < sizeof(fRefNameBuf));
    sprintf(fRefNameBuf, "%s/imgfft1/%s", context->seqDir, context_fft->fRefName);

    memset(&inVec, 0, sizeof(inVec));
    memset(&outVec, 0, sizeof(outVec));
    memset(&refVec, 0, sizeof(refVec));
    memset(&xVec, 0, sizeof(xVec));
    memset(&yVec, 0, sizeof(yVec));
    memset(&sVec, 0, sizeof(sVec));

    /* Allocate vectors for in/out/reference data. */
    if (!vecAlloc(&inVec, lenInBuffer, TE_ALIGN_YES /*| ALLOC_DRAM0*/, ((isForward) ? FMT_FRACT16 : FMT_FLOAT64) | fmtX, 0) ||
        !vecAlloc(&outVec, N, TE_ALIGN_YES,  FMT_FLOAT64 | fmtY, 0) ||
        !vecAlloc(&refVec, N, TE_ALIGN_YES,  FMT_FLOAT64 | fmtY, 0))          
    {
        printf("fileProc_stnd(): failed to allocate in/out/ref vectors, NxM=%d\n", N*M);
    }
    /* Open input data file. */
    else if (!(fIn = fopen(fInNameBuf, "rb")))
    {
        printf("fileProc_stnd(): failed to open %s for reading\n", context_fft->fInName);
    }
    /* Open reference data file. */
    else if (!(fRef = fopen(fRefNameBuf, "rb")))
    {
        printf("fileProc_stnd(): failed to open %s for reading\n", context_fft->fRefName);
    }
    else res = 1;

    /*
    * Process input data file with increasing blocks number.
    */
    if (res)
    {
        int isDone = 1;
        int baseFmtX = (isForward) ? (context->desc->fmt & FMT_DTYPE_MASK) : FMT_FRACT16;
        int baseFmtY = (isForward) ? FMT_FRACT16 : (context->desc->fmt & FMT_DTYPE_MASK);

        int isAligned = context->desc->isAligned;

        int16_t   * pin = (int16_t   *)vecGetElem(&inVec, 0);
        float64_t * pout = (float64_t *)vecGetElem(&outVec, 0);
        float64_t * pref = (float64_t *)vecGetElem(&refVec, 0);

        float32_t sinadAvg, sinadMin = HUGE_VALF;
        float64_t errSum = 0, refSum = 0;

        int efbMin = 32;
        context_fft->frameCnt = 0;
        frameProcFxn = ((tTestEngFrameProcFxn_rfft2d_stnd*)context->target.fut);

        /* Read blkNum*N 16-bit real/complex samples from the input file. */
        while (isDone)
        {
            /* Allocate in/out vectors. */
            if (!vecAlloc(&xVec, lenIn, isAligned | ALLOC_DRAM0, baseFmtX | fmtX, 0) ||             // ALLOC_DRAM0
                !vecAlloc(&yVec, lenOut, isAligned | ALLOC_DRAM0, baseFmtY | fmtY, 0))              // ALLOC_DRAM0
            {
                printf("fileProc_stnd(): failed to allocate xVec/yVec, frameCnt=%d, "
                    "fmtX=0x%x, fmtY=0x%x, MxN=%dx%d\n", context_fft->frameCnt, baseFmtX,
                    baseFmtY, M, N);
                res = 0; break;
            }

            if (context->args.caseType == TE_FFT_TESTCASE_FORWARD)
            {
                n = fread(pin, inVec.szElem, lenIn, fIn); 
                if (n == 0) break;

                if (n < lenIn) isDone = 0;
                /* Augment the (last) incomplete blocks frame with zeros. */
                memset((uint8_t*)pin + n*inVec.szElem, 0, (lenIn - n)*inVec.szElem);
            }
            else
            {
                int i;
                float64_t X0 = 0; 
                float64_t *pBlock; 
                int16_t *px = (int16_t*)vecGetElem_fr16c(&xVec, 0);
                n = 0;
                for (i = 0; i < M; i++)
                {
                    pBlock = (float64_t*)pin;// +i * N * 2;
                    n += fread(pBlock, inVec.szElem, N, fIn);
                    if (i == 0)
                    {
                        // Store DC component
                        X0 = pBlock[0];
                        pBlock[0] = 0;
                        if (n == 0) break;

                    }

                    blockExp[i] = -15 + vecPartFract16_FromFp64(&xVec, pBlock, i*N, N);

                    if (i == 0)
                    {
                        pBlock[0] = X0;
                    }
                    maxExp = MAX(maxExp, blockExp[i]);
                }

                for (i = 0; i < M; i++)
                {
                    int j;
                    int s = maxExp - blockExp[i];
                    int16_t rnd = (1 << s) >> 1;
                    for (j = 0; j < 2 * N; j++)
                    {
                        px[i * 2 * N + j] = (px[i * 2 * N + j] + rnd) >> s;
                    }
                }

                ((float64_t*)pin)[0] = X0;
                ((float64_t*)pin)[1] = maxExp;

                if (n == 0) break;
                if (n < lenIn) isDone = 0;
                /* Augment the (last) incomplete blocks frame with zeros. */
                memset((uint8_t*)pin + n*inVec.szElem, 0, (lenIn - n)*inVec.szElem);
                pin[0] = X0;
            }

            /* Zero the output buffer. */
            memset(pout, 0, N*outVec.szElem);

            int scrSize = 0;
            imgsize_t imSize = { N, M, N }; 
            if (isForward)
            {
                scrSize = (baseFmtX == FMT_FRACT16) ?   imgifft_gs16_getScratchSize(&imSize) :
                          (baseFmtX == FMT_INT8) ?      imgifft_gs8_getScratchSize(&imSize) :
                                                        imgifft_gu8_getScratchSize(&imSize);  
            }
            else
            {
                scrSize = (baseFmtX == FMT_FRACT16) ?   imgifft_gs16_getScratchSize(&imSize) :
                          (baseFmtX == FMT_INT8) ?      imgifft_gs8_getScratchSize(&imSize) :
                                                        imgifft_gu8_getScratchSize(&imSize);  
            }

            if (!vecAlloc(&sVec, scrSize, isAligned, FMT_CPLX | FMT_FRACT16, 0))
            {
                printf("fileProc_stnd(): failed to allocate sVec, frameCnt=%d, "
                    "fmt=0x%x, MxN=%dx%d\n", context_fft->frameCnt, FMT_CPLX | FMT_FRACT16, M, N);
                res = 0; break;
            }

            context_fft->dataSize = (size_t)xVec.nElem*xVec.nElem + (size_t)yVec.nElem*yVec.nElem + (size_t)sVec.nElem*sVec.nElem;

            /* Use a proprietary frame processing function to perform the target FFT. */
            if (!(res = frameProcFxn(context, (int16_t*)pin, (float64_t*)pout, &xVec, &yVec, &sVec))) break;

            /* Free in vectors. */
            if (!(res = vecsFree(&xVec, 0)))
            {
                printf("fileProc_stnd(): vecsFree() failed for xVec, frameCnt=%d, "
                    "fmtX=0x%x, fmtY=0x%x,  MxN=%dx%d\n", context_fft->frameCnt, baseFmtX, baseFmtY, M, N);
                break;
            }

            if (!(res = vecsFree(&sVec, 0)))
            {
                printf("fileProc_stnd(): vecsFree() failed for sVec, frameCnt=%d, "
                    "MxN=%dx%d\n", context_fft->frameCnt, M, N);
                break;
            }

            /* Estimate the SINAD ratio for the current blocks sequence, and update error/reference
            * power sums for the whole file SINAD. */
            //            for (blkIx = 0; blkIx<blkNum; blkIx++)
            {
                int m, n;
                float64_t refMax = 0, errMax = 0;
                float64_t p, q;
                int k;// , LenOut2;
                int lenRow = (fmtY & FMT_CPLX) ? 2 * N : N;
                int32_t totalScaling = * vecGetElem_i32(&context->dataSet.Y, 0); 
                int32_t meanQ31_0    = *vecGetElem_i32(&context->dataSet.Y, 1);

                p = q = 0;
                k = 0;

                for (m = 0; m < M; m++)
                {
                    /* Read double precision reference samples. Reference block always contains
                    * exactly N samples, real or complex! */
                    int nRead = fread(pref, refVec.szElem, N, fRef);
                    /* Augment the (last) incomplete block frame with zeros. */
                    memset((uint8_t*)pref + nRead*refVec.szElem, 0, (N - nRead)*refVec.szElem);

                    vecPartToFp64(pout, &yVec, totalScaling, m*N, N);

                    if (context->args.caseType == TE_FFT_TESTCASE_FORWARD && m==0)
                    {
                      double dcBias = (context->desc->fmt == FMT_UINT8) ? (128.0 * meanQ31_0*M * N) :
                        (context->desc->fmt == FMT_INT8) ? (128.0 * (meanQ31_0 + 128) *M * N) : (double)meanQ31_0*M*N;
                        pout[0] += dcBias;
                    }

                    for (n = 0; n < lenRow; n++)
                    {
                        float64_t err = pout[n] - pref[n];

                        /* |signal|^2 */
                        p += pref[n] * pref[n];
                        /* |noise+distortion|^2 */
                        q += err*err;

                        refMax = MAX(refMax, fabs(pref[n]));
                        errMax = MAX(errMax, fabs(err));
                        k++;
                    }
                }

                refSum += p;
                errSum += q;

                if (p>0)
                {
                    sinadMin = MIN(sinadMin, (float32_t)(p / q));
                    efbMin = MAX(0, MIN(efbMin, ((int)logb(refMax) - (int)logb(errMax))));
                }
            }
            /* Free out vectors. */
            if (!(res = vecsFree(&yVec, 0)))
            {
                printf("fileProc_stnd(): vecsFree() failed for yVec, frameCnt=%d, "
                    "fmtX=0x%x, fmtY=0x%x,  MxN=%dx%d\n", context_fft->frameCnt, baseFmtX, baseFmtY, M, N);
                break;
            }
            context_fft->frameCnt++;
        }

        /*
        * Finalize the min SINAD estimation and print a summary.
        */
        if (res)
        {
            sinadMin = 10.f*log10f(sinadMin);

            /* Set the test result for test verification. */
            *vecGetElem_fl32(&context->dataSet.Z, 0) = sinadMin;

            if (context->isVerbose)
            {
                sinadAvg = (refSum > 0 ? 10.f*log10f((float32_t)(refSum / errSum)) : -HUGE_VALF);

                printf("Error-Free Bits %2d  SINAD min %5.1f dB  avg %5.1f dB  ",
                        efbMin, sinadMin, sinadAvg);
            }
        }
    }

    /*
    * Close files and free vectors.
    */
    if (fIn) fclose(fIn);
    if (fRef) fclose(fRef);

    if (context->dataSet.Y.szBulk > 0)
    {
        vecFree(&context->dataSet.Y); 
        memset(&context->dataSet.Y, 0, sizeof(context->dataSet.Y));
    }

    if (inVec.szBulk > 0) vecFree(&inVec);
    if (outVec.szBulk > 0) vecFree(&outVec);
    if (refVec.szBulk > 0) vecFree(&refVec);
    if (xVec.szBulk > 0) vecFree(&xVec);
    if (yVec.szBulk > 0) vecFree(&yVec);
    if (sVec.szBulk > 0) vecFree(&sVec);
} /* fileProc_stnd() */

#include "NatureDSP_Signal_fft.h"  

int func_imgfft1(int isFull, int isVerbose, int breakOnError)
{
    int res = 1;
    int runAlways;
    vecInitRegion(ALLOC_DRAM0);
    #define DO_TESTU8(runAlways, fxn, seqFile)                                                      \
    {                                                                                               \
      if (runAlways) {                                                                              \
        tTestEngTarget     target = (tTestEngTarget)&te_frameProc_imgfft;                           \
        tTestEngDesc_fft desc = TEST_DESC_BRFFT_NOSCALING(FMT_UINT8, fxn, fxn);                     \
        res &= (0 != TestEngRun(target, &desc.desc, "imgfft1/" seqFile, isFull, isVerbose, breakOnError, 0));  \
        if (!(res || !breakOnError))                                                                \
            return res;                                                                             \
      }                                                                                             \
    }
    #define DO_TESTS8(runAlways, fxn, seqFile)                                                      \
    {                                                                                               \
      if (runAlways) {                                                                              \
        tTestEngTarget     target = (tTestEngTarget)&te_frameProc_imgfft;                           \
        tTestEngDesc_fft desc = TEST_DESC_BRFFT_NOSCALING(FMT_INT8, fxn, fxn);                      \
        res &= (0 != TestEngRun(target, &desc.desc, "imgfft1/" seqFile, isFull, isVerbose, breakOnError, 0));  \
        if (!(res || !breakOnError))                                                                \
            return res;                                                                             \
      }                                                                                             \
    }
    #define DO_TEST16(runAlways, fxn, seqFile)                                                      \
    {                                                                                               \
      if (runAlways) {                                                                              \
        tTestEngTarget     target = (tTestEngTarget)&te_frameProc_imgfft;                           \
        tTestEngDesc_fft desc = TEST_DESC_BRFFT_NOSCALING(FMT_FRACT16, fxn, fxn);                   \
        res &= (0 != TestEngRun(target, &desc.desc, "imgfft1/" seqFile, isFull, isVerbose, breakOnError, 0));  \
        if (!(res || !breakOnError))                                                                \
            return res;                                                                             \
      }                                                                                             \
    }

    runAlways = (isFull == 2) ? 0 : 1;/* 0 - for sanity test; 1 - for sanity & brief & full */
    DO_TESTU8 (         1, &imgfft_gu8,"imgfft_gu8_64x64.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_128x128.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_256x256.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_512x512.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_64x640.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_96x576.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_128x512.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_144x480.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_176x384.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_240x352.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_256x320.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_288x288.seq");       
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_640x64.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_576x96.seq");        
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_512x128.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_480x144.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_384x176.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_352x240.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_320x256.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_288x288.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_SQCIF.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_QCIF.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_QVGA.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_CIF.seq");
    DO_TESTU8 ( runAlways, &imgfft_gu8,"imgfft_gu8_VGA.seq");

    DO_TESTS8 (         1, &imgfft_gs8,"imgfft_gs8_64x64.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_128x128.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_256x256.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_512x512.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_64x640.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_96x576.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_128x512.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_144x480.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_176x384.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_240x352.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_256x320.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_288x288.seq");       
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_640x64.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_576x96.seq");        
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_512x128.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_480x144.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_384x176.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_352x240.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_320x256.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_288x288.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_SQCIF.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_QCIF.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_QVGA.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_CIF.seq");
    DO_TESTS8 ( runAlways, &imgfft_gs8,"imgfft_gs8_VGA.seq");

    DO_TEST16 (         1, &imgfft_gs16,"imgfft_gs16_64x64.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_128x128.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_256x256.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_512x512.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_64x640.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_96x576.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_128x512.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_144x480.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_176x384.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_240x352.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_256x320.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_288x288.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_640x64.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_576x96.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_512x128.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_480x144.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_384x176.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_352x240.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_320x256.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_288x288.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_SQCIF.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_QCIF.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_QVGA.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_CIF.seq");
    DO_TEST16 ( runAlways, &imgfft_gs16,"imgfft_gs16_VGA.seq");

    DO_TESTU8 (         1, &imgifft_gu8,"imgifft_gu8_64x64.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_128x128.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_256x256.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_512x512.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_64x640.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_96x576.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_128x512.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_144x480.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_176x384.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_240x352.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_256x320.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_288x288.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_640x64.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_576x96.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_512x128.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_480x144.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_384x176.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_352x240.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_320x256.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_288x288.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_SQCIF.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_QCIF.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_QVGA.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_CIF.seq");
    DO_TESTU8 ( runAlways, &imgifft_gu8,"imgifft_gu8_VGA.seq");

    DO_TESTS8 (         1, &imgifft_gs8,"imgifft_gs8_64x64.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_128x128.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_256x256.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_512x512.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_64x640.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_96x576.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_128x512.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_144x480.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_176x384.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_240x352.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_256x320.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_288x288.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_640x64.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_576x96.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_512x128.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_480x144.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_384x176.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_352x240.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_320x256.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_288x288.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_SQCIF.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_QCIF.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_QVGA.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_CIF.seq");
    DO_TESTS8 ( runAlways, &imgifft_gs8,"imgifft_gs8_VGA.seq");

    DO_TEST16 (         1, &imgifft_gs16,"imgifft_gs16_64x64.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_128x128.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_256x256.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_512x512.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_64x640.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_96x576.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_128x512.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_144x480.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_176x384.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_240x352.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_256x320.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_288x288.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_640x64.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_576x96.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_512x128.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_480x144.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_384x176.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_352x240.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_320x256.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_288x288.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_SQCIF.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_QCIF.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_QVGA.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_CIF.seq");
    DO_TEST16 ( runAlways, &imgifft_gs16,"imgifft_gs16_VGA.seq");

#undef DO_TESTU8
#undef DO_TESTS8
#undef DO_TEST16

    return (res);
}
