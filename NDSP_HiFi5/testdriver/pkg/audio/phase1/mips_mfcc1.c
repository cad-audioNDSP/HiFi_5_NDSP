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
* Test module for testing cycle performance (MFCC features extraction APIs)
*/

#include <string.h>
#include <math.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API. */
#include LIBRARY_HEADER(audio)
#include LIBRARY_HEADER(fft)
#include LIBRARY_HEADER(vector)
/* MIPS measurement means. */
#include "mips.h"
/* Utility functions and macros. */
#include "utils.h"

#define REAL_f      f32
#define REAL_32x32  i32
#define CPLX_f      cf32
#define CPLX_32x32  ci32

#define IS_API_PRESENT(apiName, prec)  (IS_PRESENT(apiName##prec##_alloc) && \
                                        IS_PRESENT(apiName##prec##_init) && \
                                        IS_PRESENT(apiName##prec##_getScratchSize))

#define OPT_HTK       (LOGMEL_OPT_MELSCALE_HTK|LOGMEL_OPT_FBELOG_NATURAL|LOGMEL_OPT_FBNORM_NONE| \
                       MFCC_OPT_DC_MEAN_REMOVE|MFCC_OPT_PREEMPH_FRAMEBYFRAME|MFCC_OPT_DCT_NORMALIZED)
#define OPT_AUDITORY  (LOGMEL_OPT_MELSCALE_AUDITORY|LOGMEL_OPT_FBELOG_BASE10 |LOGMEL_OPT_FBNORM_AREA| \
                       MFCC_OPT_DC_MEAN_DONT_REMOVE|MFCC_OPT_PREEMPH_CONTINUOUS|MFCC_OPT_DCT_ORTHOGONAL)

#define PROFILE_LOGMEL(prec, isFull, isVerbose, _Fs, _fftSize, lowFreq, uppFreq, bandNum, flavor) { \
        logmel##prec##_handle_t logmel = NULL; \
        size_t szObj, szScr; \
        logmel_params_t params; \
        if (IS_API_PRESENT(logmel, prec)) { \
            memset(&params, 0, sizeof(params)); \
            params.Fs = _Fs; params.fftSize = _fftSize; params.mfbBandNum = bandNum; \
            params.mfbLowFreqQ8 = (fract32)round(ldexp(lowFreq, 8)); \
            params.mfbUppFreqQ8 = (fract32)round(ldexp(uppFreq, 8)); \
            params.opt = OPT_##flavor; \
            szObj = logmel##prec##_alloc(&params); szScr = logmel##prec##_getScratchSize(&params); \
            (void)szObj; NASSERT(szObj>0 && szObj<=sizeof(objinstance_memory)); /* Just to validate the parameters. */ \
            (void)szScr; NASSERT(szScr<=sizeof(mips.scratch0)); \
            logmel = logmel##prec##_init(objinstance_memory, &params); NASSERT(NULL!=logmel); \
        } \
        PROFILE_SIMPLE(isFull, isVerbose, logmel##prec##_process, \
                       (logmel, &mips.scratch0.u8[0], &mips.out0.REAL_##prec[0], &mips.inp0.CPLX_##prec[0], 15), fout, \
                       "Fs: "#_Fs"; fftSize: "#_fftSize"; Range: "#lowFreq"-"#uppFreq" Hz; Bands: "#bandNum"; Flavor: "#flavor, \
                       "%d (cycles per STFT hop)"); \
    }

#define PROFILE_MFCC(prec, isFull, isVerbose, _Fs, _fftSize, winLenMs, hopLenMs, lowFreq, uppFreq, bandNum, _cepstraNum, flavor) { \
        mfcc##prec##_handle_t mfcc = NULL; \
        size_t szObj, szScr; \
        mfcc_params_t params; \
        if (IS_API_PRESENT(mfcc, prec) && IS_PRESENT(mfcc_getDefaultParams)) { \
            mfcc_getDefaultParams(&params); \
            params.Fs = _Fs; params.scaleExp = 15; params.preemph = 31785; params.fftSize = _fftSize; \
            params.stftWinLen = (fract32)round((float64_t)winLenMs*_Fs/1000); \
            params.stftHopLen = (fract32)round((float64_t)hopLenMs*_Fs/1000); \
            params.mfbLowFreqQ8 = (fract32)round(ldexp(lowFreq, 8)); \
            params.mfbUppFreqQ8 = (fract32)round(ldexp(uppFreq, 8)); \
            params.mfbBandNum = bandNum; params.cepstraNum = _cepstraNum; \
            params.lifter = (OPT_##flavor == OPT_HTK ? 22 : 0); \
            params.opt = OPT_##flavor; \
            szObj = mfcc##prec##_alloc(&params); szScr = mfcc##prec##_getScratchSize(&params); \
            (void)szObj; NASSERT(szObj>0 && szObj<=sizeof(objinstance_memory)); /* Just to validate the parameters. */ \
            (void)szScr; NASSERT(szScr<=sizeof(mips.scratch0)); \
            mfcc = mfcc##prec##_init(objinstance_memory, &params, &mfcc##prec##_cbk); NASSERT(NULL!=mfcc); \
        } \
        PROFILE_SIMPLE(isFull, isVerbose, mfcc##prec##_process, \
                       (mfcc, &mips.scratch0.u8[0], &mips.out0.REAL_##prec[0], &mips.inp0.REAL_##prec[0]), fout, \
                       "Fs: "#_Fs"; fftSize: "#_fftSize"; Win: "#winLenMs" ms; Hop: "#hopLenMs" ms; " \
                       "Range: "#lowFreq"-"#uppFreq" Hz; Bands: "#bandNum"; Ceps: "#_cepstraNum"; Flavor: "#flavor, \
                       "%d (cycles per STFT hop)"); \
    }

#define PROFILE_htkdelta(prec,isFull,isVerbose,fun,M,N)     PROFILE_INVERTED(isFull,isVerbose,fun,(mips.out0.REAL_##prec, mips.inp1.REAL_##prec, M, N),fout,"M="#M" N="#N,prf_maccycle, (2*N*M));

/* STFT weighting window generator callback function. */
static void mfcc32x32_genWindow_cbk( void * host, fract32  * window, int len )
{
    (void)host; (void)window; (void)len;
}

/* Real-to-complex FFT callback function. */
static int mfcc32x32_rfft_cbk( void * host, complex_fract32 * restrict y, fract32 * restrict x, int fftSize )
{
    int hidx;
    fft_handle_t hfft[] = {rfft32_256, rfft32_512, rfft32_1024, rfft32_2048};
    (void)host;
    hidx = (30-8)-scl_bexp32(fftSize);
    NASSERT(0<=hidx && hidx<(int)sizeof(hfft)/sizeof(hfft[0]));
    return fft_real32x32((int32_t*)y, x, hfft[hidx], 2); /* scalingOpt == 2: dynamic scaling */
} /* mfcc32x32_rfft() */

/* User-defined callbacks for an MFCC extractor. */
static const mfcc32x32_callback_t mfcc32x32_cbk = {
  NULL, &mfcc32x32_genWindow_cbk, &mfcc32x32_rfft_cbk
};

void mips_mfcc1(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_LOGMEL(32x32, isFull, isVerbose,  8000,  256,         133.3, 3700.0, 20,          HTK);
    PROFILE_LOGMEL(32x32,      1, isVerbose, 16000,  512,         133.3, 6853.8, 20,          HTK);
    PROFILE_LOGMEL(32x32, isFull, isVerbose, 24000, 1024,         133.3, 6853.8, 20,          HTK);
    PROFILE_LOGMEL(32x32, isFull, isVerbose, 32000, 2048,         133.3, 6853.8, 20,          HTK);
    PROFILE_LOGMEL(32x32, isFull, isVerbose,  8000,  256,         133.3, 3700.0, 40,     AUDITORY);
    PROFILE_LOGMEL(32x32,      1, isVerbose, 16000,  512,         133.3, 6853.8, 40,     AUDITORY);
    PROFILE_LOGMEL(32x32, isFull, isVerbose, 24000, 1024,         133.3, 6853.8, 40,     AUDITORY);
    PROFILE_LOGMEL(32x32, isFull, isVerbose, 32000, 2048,         133.3, 6853.8, 40,     AUDITORY);

    PROFILE_MFCC  (32x32, isFull, isVerbose,  8000,  256, 25, 10, 133.3, 3700.0, 20, 13,      HTK);
    PROFILE_MFCC  (32x32,      1, isVerbose, 16000,  512, 25, 10, 133.3, 6853.8, 20, 13,      HTK);
    PROFILE_MFCC  (32x32, isFull, isVerbose, 24000, 1024, 25, 10, 133.3, 6853.8, 20, 13,      HTK);
    PROFILE_MFCC  (32x32, isFull, isVerbose, 32000, 2048, 30, 10, 133.3, 6853.8, 20, 13,      HTK);
    PROFILE_MFCC  (32x32, isFull, isVerbose,  8000,  256, 16, 10, 133.3, 3700.0, 40, 13, AUDITORY);
    PROFILE_MFCC  (32x32,      1, isVerbose, 16000,  512, 16, 10, 133.3, 6853.8, 40, 13, AUDITORY);
    PROFILE_MFCC  (32x32, isFull, isVerbose, 24000, 1024, 16, 10, 133.3, 6853.8, 40, 13, AUDITORY);
    PROFILE_MFCC  (32x32, isFull, isVerbose, 32000, 2048, 30, 10, 133.3, 6853.8, 40, 13, AUDITORY);

    PROFILE_htkdelta(32x32,      1, isVerbose, htkdelta32x32, 13, 2);
    PROFILE_htkdelta(32x32, isFull, isVerbose, htkdelta32x32, 15, 4);
    PROFILE_htkdelta(32x32, isFull, isVerbose, htkdelta32x32, 20, 6);
    PROFILE_htkdelta(32x32, isFull, isVerbose, htkdelta32x32, 24, 8);
}
