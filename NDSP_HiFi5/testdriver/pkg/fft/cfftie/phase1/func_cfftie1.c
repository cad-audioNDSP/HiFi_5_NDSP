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
#include "../common/test_cfftie.h"

/* Test data file processing for a standalone test. Applies the target FFT function
* to test data, compares the FFT output with reference data and estimates the SINAD. */
static void fileProc_stnd_cfftie(tTestEngContext * context)
{
    tTestEngContext_fft           * context_fft;
    tTestEngFrameProcFxn_fft_stnd * frameProcFxn;

    tVec inVec, outVec, refVec; /* In/out/reference vectors related to test data files. */
    tVec xVec, yVec;            /* In/out vectors for the target FFT routine. */
    FILE *fIn = 0, *fRef = 0;
    int res = 0;
    int n, N;

    NASSERT(context && context->desc);
    NASSERT(context->target.fut && context->target.handle);

    context_fft = (tTestEngContext_fft*)context->target.handle;

    /* FFT size */
    N = context->args.N_DIM_;
    if (context->desc->extraParam & TE_FFT_STEREO) N*=2;

    /* Preset an invalid SINAD for the test result. */
    *vecGetElem_fl32(&context->dataSet.Z, 0) = -HUGE_VALF;

    if (context->isVerbose)
    {
        if (context->desc->extraParam & TE_FFT_OPT_SCALE_METH)
        {
            printf("scale_mtd %d ", context_fft->scale_method);
        }

        printf("%-40s  ", context_fft->fRefName);
    }

    memset(&inVec, 0, sizeof(inVec));
    memset(&outVec, 0, sizeof(outVec));
    memset(&refVec, 0, sizeof(refVec));
    memset(&xVec, 0, sizeof(xVec));
    memset(&yVec, 0, sizeof(yVec));
    {
        char * fname;
        const char *dir;        
		dir = getVectorsDir(context->isFull);
        fname=(char *)malloc(strlen(context_fft->fInName)+strlen(dir)+2+6);
        sprintf(fname,"%s/%s",dir,context_fft->fInName);
        fIn = fopen(fname, "rb");
        free(fname);
        fname=(char *)malloc(strlen(context_fft->fRefName)+strlen(dir)+2+6);
        sprintf(fname,"%s/%s",dir,context_fft->fRefName);
        fRef = fopen(fname, "rb");
        free(fname);
    }

    /* Allocate vectors for in/out/reference data. */
    if (!vecAlloc(&inVec, N, TE_ALIGN_YES, FMT_FRACT16 | FMT_CPLX, 0) ||
        !vecAlloc(&outVec, N, TE_ALIGN_YES, FMT_FLOAT64 | FMT_CPLX, 0) ||
        !vecAlloc(&refVec, N, TE_ALIGN_YES, FMT_FLOAT64 | FMT_CPLX, 0))
    {
        printf("fileProc_stnd_cfftie(): failed to allocate in/out/ref vectors, N=%d\n", N);
    }
    /* Open input data file. */
    else if (fIn==NULL)
    {
        printf("fileProc_stnd_cfftie(): failed to open %s for reading\n", context_fft->fInName);
    }
    /* Open reference data file. */
    else if (fRef==NULL)
    {
        printf("fileProc_stnd_cfftie(): failed to open %s for reading\n", context_fft->fRefName);
    }
    else res = 1;

    /*
    * Process the input file frame-by-frame.
    */

    if (res)
    {
        int fmt = context->desc->fmt | FMT_CPLX;
        int isAligned = context->desc->isAligned;

        complex_fract16 * pin = vecGetElem_fr16c(&inVec, 0);
        complex_double  * pout = vecGetElem_fl64c(&outVec, 0);
        complex_double  * pref = vecGetElem_fl64c(&refVec, 0);

        float32_t sinadAvg, sinadMin = HUGE_VALF;
        float64_t errSum = 0, refSum = 0;

        int efbMin = 32;

        memset(&xVec, 0, sizeof(xVec));
        memset(&yVec, 0, sizeof(yVec));

        context_fft->frameCnt = 0;

        frameProcFxn = ((tTestEngFrameProcFxn_fft_stnd*)context->target.fut);

        /* Read N 16-bit complex samples from the input file. */
        while ((n = fread(pin, sz_fr16c, N, fIn)) > 0)
        {
            /* Augment the (last) incomplete frame with zeros. */
            memset((uint8_t*)pin + n*sz_fr16c, 0, (N - n)*sz_fr16c);
            /* Zero the output frame. */
            memset(pout, 0, N*sz_fl64c);

            /* Allocate in/out buffers for the target FFT routine. */
            if (!xVec.szBulk && (2 != vecsAlloc(isAligned, fmt, &xVec, N, &yVec, N, 0)))
            {
                printf("fileProc_stnd_cfftie(): failed to allocate xVec/yVec, "
                    "frameCnt=%d, fmt=0x%x, N=%d\n", context_fft->frameCnt, fmt, N);
                res = 0; break;
            }
            context_fft->dataSize = (size_t)xVec.szElem*xVec.nElem + (size_t)yVec.szElem*yVec.nElem;

            /* Use a proprietary frame processing function to perform the target FFT. */
            if (!(res = frameProcFxn(context, (int16_t*)pin, (float64_t*)pout, &xVec, &yVec))) break;

            /* When in unaligned mode, in/out vectors should be periodically reallocated to try various
            * address offsets. */
            if (!(++context_fft->frameCnt % XY_VEC_REALLOC_PERIOD) && !isAligned)
            {
                if (!(res = vecsFree(&xVec, &yVec, 0)))
                {
                    printf("fileProc_stnd_cfftie(): vecsFree() failed for xVec/yVec, "
                        "frameCnt=%d, fmt=0x%x, N=%d\n", context_fft->frameCnt, fmt, N);
                    break;
                }

                memset(&xVec, 0, sizeof(xVec));
                memset(&yVec, 0, sizeof(yVec));
            }

            /* Read N double precision complex reference samples. */
            if ((int)fread(pref, sz_fl64c, N, fRef) < N)
            {
                printf("fileProc_stnd_cfftie(): failed to read reference data, "
                    "frameCnt=%d, fmt=0x%x, N=%d\n", context_fft->frameCnt, fmt, N);
                res = 0; break;
            }

            /* Estimate the SINAD ratio for the current frame, and update error/reference
            * power sums for the whole file SINAD. */
            {
                float64_t refMax = 0, errMax = 0;
                float64_t p, q;

                for (p = q = 0, n = 0; n<N; n++)
                {
                    double err_re, err_im;
                    err_re = creal(pout[n]) - creal(pref[n]);
                    err_im = cimag(pout[n]) - cimag(pref[n]);

                    /* |signal|^2 */
                    p += creal(pref[n])*creal(pref[n]) + cimag(pref[n])*cimag(pref[n]);
                    /* |noise+distortion|^2 */
                    q += err_re*err_re + err_im*err_im;

                    refMax = MAX(refMax, MAX(fabs(creal(pref[n])), fabs(cimag(pref[n]))));
                    errMax = MAX(errMax, MAX(fabs(err_re), fabs(err_im)));
                }

                refSum += p;
                errSum += q;

                if (p>0)
                {
                    sinadMin = MIN(sinadMin, (float32_t)(p / q));
                    efbMin = MAX(0, MIN(efbMin, L_sub_ll((int)logb(refMax), (int)logb(errMax))));
                }
            }
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

                if ((fmt & FMT_DTYPE_MASK) == FMT_FRACT16 || (fmt & FMT_DTYPE_MASK) == FMT_FRACT32)
                {
                    printf("Error-Free Bits %2d  SINAD min %5.1f dB  avg %5.1f dB  ",
                        efbMin, sinadMin, sinadAvg);
                }
                else
                {
                    printf("SINAD min %5.1f dB  avg %5.1f dB  ", sinadMin, sinadAvg);
                }
            }
        }
    }

    /*
    * Close files and free vectors.
    */

    if (fIn) fclose(fIn);
    if (fRef) fclose(fRef);

    if (inVec.szBulk > 0) vecFree(&inVec);
    if (outVec.szBulk > 0) vecFree(&outVec);
    if (refVec.szBulk > 0) vecFree(&refVec);
    if (xVec.szBulk > 0) vecFree(&xVec); //!!! xVec not initialized !!!
    if (yVec.szBulk > 0) vecFree(&yVec);

} /* fileProc_stnd_cfftie() */

/* Initializer for a test description structure. */
#define TEST_DESC_STND( fmt, extraParam, align,fwdfft,invfft )  { { (fmt),(extraParam), NULL, TE_ARGNUM_1,(align), \
    &te_create_fft, &te_destroy_fft, \
    &te_load_fft, &fileProc_stnd_cfftie }, (tFrwTransFxn)fwdfft, (tInvTransFxn)invfft}

TestDef_t testTbl_cfftie1[] =
{
#if 1
    { 1, "cfftie1/fft128_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx16x16_ie, NULL) },
    { 1, "cfftie1/fft256_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx16x16_ie, NULL) },
    { 1, "cfftie1/fft512_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx16x16_ie, NULL) },
    { 1, "cfftie1/fft1024_cplx16x16_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx16x16_ie, NULL) },

    { 1, "cfftie1/ifft128_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx16x16_ie) },
    { 1, "cfftie1/ifft256_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx16x16_ie) },
    { 1, "cfftie1/ifft512_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx16x16_ie) },
    { 1, "cfftie1/ifft1024_cplx16x16_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx16x16_ie) },
#endif
#if 0 //HiFi3/3z API
    { 1, "cfftie1/fft_cplx24x24_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx24x24_ie, NULL) },
    { 1, "cfftie1/ifft_cplx24x24_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx24x24_ie) },
#endif
#if 1
    { 1, "cfftie1/fft256_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx32x16_ie, NULL) },
    { 1, "cfftie1/fft512_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx32x16_ie, NULL) },
    { 1, "cfftie1/fft1024_cplx32x16_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx32x16_ie, NULL) },

    { 1, "cfftie1/fft256_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx32x16_ie, NULL) },
    { 1, "cfftie1/fft512_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx32x16_ie, NULL) },
    { 1, "cfftie1/fft1024_cplx32x16_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, &fft_cplx32x16_ie, NULL) },

    { 1, "cfftie1/ifft256_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx32x16_ie) },
    { 1, "cfftie1/ifft512_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx32x16_ie) },
    { 1, "cfftie1/ifft1024_cplx32x16_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx32x16_ie) },

    { 1, "cfftie1/ifft256_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx32x16_ie) },
    { 1, "cfftie1/ifft512_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx32x16_ie) },
    { 1, "cfftie1/ifft1024_cplx32x16_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &ifft_cplx32x16_ie) },
#endif
#if 1
    { 1, "cfftie1/fft128_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },
    { 1, "cfftie1/fft256_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },
    { 1, "cfftie1/fft512_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },
    { 1, "cfftie1/fft1024_cplx32x32_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },

    { 1, "cfftie1/fft128_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },
    { 1, "cfftie1/fft256_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },
    { 1, "cfftie1/fft512_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },
    { 1, "cfftie1/fft1024_cplx32x32_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32_ie, NULL) },

    { 1, "cfftie1/ifft128_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },
    { 1, "cfftie1/ifft256_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },
    { 1, "cfftie1/ifft512_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },
    { 1, "cfftie1/ifft1024_cplx32x32_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },

    { 1, "cfftie1/ifft128_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },
    { 1, "cfftie1/ifft256_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },
    { 1, "cfftie1/ifft512_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },
    { 1, "cfftie1/ifft1024_cplx32x32_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32_ie) },
#endif
#if 1
    { 1, "cfftie1/stereo_fft128_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx16x16_ie, NULL) },
    { 1, "cfftie1/stereo_fft256_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx16x16_ie, NULL) },
    { 1, "cfftie1/stereo_fft512_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx16x16_ie, NULL) },
    { 1, "cfftie1/stereo_fft1024_cplx16x16_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx16x16_ie, NULL) },

    { 1, "cfftie1/stereo_ifft128_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx16x16_ie) },
    { 1, "cfftie1/stereo_ifft256_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx16x16_ie) },
    { 1, "cfftie1/stereo_ifft512_cplx16x16_ie.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx16x16_ie) },
    { 1, "cfftie1/stereo_ifft1024_cplx16x16_ie.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx16x16_ie) },
#endif
#if 1
    { 1, "cfftie1/stereo_fft256_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx32x16_ie, NULL) },
    { 1, "cfftie1/stereo_fft512_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx32x16_ie, NULL) },
    { 1, "cfftie1/stereo_fft1024_cplx32x16_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx32x16_ie, NULL) },

    { 1, "cfftie1/stereo_fft256_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx32x16_ie, NULL) },
    { 1, "cfftie1/stereo_fft512_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx32x16_ie, NULL) },
    { 1, "cfftie1/stereo_fft1024_cplx32x16_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, &stereo_fft_cplx32x16_ie, NULL) },

    { 1, "cfftie1/stereo_ifft256_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x16_ie) },
    { 1, "cfftie1/stereo_ifft512_cplx32x16_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x16_ie) },
    { 1, "cfftie1/stereo_ifft1024_cplx32x16_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x16_ie) },
              
    { 1, "cfftie1/stereo_ifft256_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x16_ie) },
    { 1, "cfftie1/stereo_ifft512_cplx32x16_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x16_ie) },
    { 1, "cfftie1/stereo_ifft1024_cplx32x16_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16 | TE_FFT_STEREO | TE_FFT_TWD16), TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x16_ie) },
#endif
#if 1
    { 1, "cfftie1/stereo_fft128_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },
    { 1, "cfftie1/stereo_fft256_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },
    { 1, "cfftie1/stereo_fft512_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },
    { 1, "cfftie1/stereo_fft1024_cplx32x32_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },
              
    { 1, "cfftie1/stereo_fft128_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },
    { 1, "cfftie1/stereo_fft256_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },
    { 1, "cfftie1/stereo_fft512_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },
    { 1, "cfftie1/stereo_fft1024_cplx32x32_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, &stereo_fft_cplx32x32_ie, NULL) },
              
    { 1, "cfftie1/stereo_ifft128_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
    { 1, "cfftie1/stereo_ifft256_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
    { 1, "cfftie1/stereo_ifft512_cplx32x32_ies3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
    { 1, "cfftie1/stereo_ifft1024_cplx32x32_ies3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
              
    { 1, "cfftie1/stereo_ifft128_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
    { 1, "cfftie1/stereo_ifft256_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
    { 1, "cfftie1/stereo_ifft512_cplx32x32_ies2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
    { 1, "cfftie1/stereo_ifft1024_cplx32x32_ies2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_STEREO, TE_ALIGN_YES, NULL, &stereo_ifft_cplx32x32_ie) },
#endif
    { 0 } /* End of table */
}; // testTbl_cfftie

int func_cfftie1(int isFull, int isVerbose, int breakOnError)
{
    return main_cfftie(1, isFull, isVerbose, breakOnError, testTbl_cfftie1);
}
