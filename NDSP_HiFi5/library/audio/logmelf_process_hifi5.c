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
  NatureDSP Signal Processing Library. Audio processing part
    Compute log mel filterbank energies
    Single precision floating-point variant
    C code optimized for HiFi5 core
  IntegrIT, 2006-2019
*/
/*
 * References:
 * [1] S. Young, G. Evermann, M. Gales, T. Hain, D. Kershaw, X. Liu, G. Moore,
 *     J. Odell, D. Ollason, D. Povey, V. Valtchev, P. Woodland,
 *     The HTK Book (for HTK version 3.4), 
 *     Cambridge University Engineering Department, 2009. 
 *     http://htk.eng.cam.ac.uk/docs/docs.shtml
 * [2] Auditory Toolbox for MATLAB by Malcolm Slaney
 *     https://uk.mathworks.com/matlabcentral/linkexchange/links/38-auditory-toolbox
 * [3] Reproducing the feature outputs of common programs using Matlab and 
 *     melfcc.m, by Dan Ellis
 *     http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/mfccs.html
 */

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "common_fpu.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_audio.h"
#include "NatureDSP_Signal_math.h"
/*  Log mel filterbank internal definitions. */
#include "logmel_internal.h"
/* +/-infinity floating-point constants */
#include "inff_tbl.h"

#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(size_t, logmelf_getScratchSize, (const logmel_params_t * params));
DISCARD_FUN(void  , logmelf_process       , (logmelf_handle_t handle, void * restrict pScr, float32_t * restrict logFbe,
                                             const complex_float * restrict spectra, int scaleExp));
#else

#define PROFILE_ENABLE  0 /* If non-zero, measure cycles and print a report to stdout. */
#define ALIGN_SIZE      (HIFI_SIMD_WIDTH)
#define ALIGN_PAD       (ALIGN_SIZE-1)
#define ALIGN_PTR(p)    (void*)(((uintptr_t)(p)+ALIGN_PAD)&~ALIGN_PAD)
#define sz_i16          sizeof(int16_t)
#define sz_f32          sizeof(float32_t)
#define sz_f32c         sizeof(complex_float)

/* Integrated profiler. */
#include "profile.h"

PROFILE_CREATE(logmelf_vec_complex2mag);
PROFILE_CREATE(logmelf_computeFbe);
PROFILE_CREATE(logmelf_vec_logXf);

/*-------------------------------------------------------------------------
  Complex magnitude
  Routine computes complex magnitude

  Input:
  N            length of vector
  x[N]         input complex data
  Output:
  y[N]         output data
  Temporary:
  pScr[N]      scratch buffer of MAX(8,N) elements

  Restriction:
  x,y,pScr     Must be 16-bytes aligned
-------------------------------------------------------------------------*/
static void vec_complex2mag_logmel(float32_t * pScr, float32_t * restrict y, const complex_float * restrict x, int N);

/*-------------------------------------------------------------------------
  Compute log mel filterbank energies
  In the initialization stage, split the frequency range [lowFreq..uppFreq]
  into 1/2 overlapping bands of equal mel-frequency width, and compute triangular
  weighting function for each band. 
  Data processing functions are applied to Fourier image of real signal, specified
  through the input argument spectra[fftSize/2+1]. Fourier image is converted to
  magnitude spectrum. For every mel-frequency band, magnitude samples are multiplied
  by the corresponding triangular weighting function, and summed together to form 
  filterbank energies (FBEs). Finally, log-scaled FBEs are stored to the output
  argument logFbe[mfbBandNum].
  Log mel filterbank routines follow the algorithm used in the Hidden Markov Models
  Toolkit (HTK), as descibed in:
  [1] S. Young, G. Evermann, M. Gales, T. Hain, D. Kershaw, X. Liu, G. Moore,
      J. Odell, D. Ollason, D. Povey, V. Valtchev, P. Woodland,
      The HTK Book (for HTK version 3.4), 
      Cambridge University Engineering Department, 2009. 
      http://htk.eng.cam.ac.uk/docs/docs.shtml
  In a few aspects, such as the linear to mel-scale frequency mapping, the implementation
  may optionally mimic another popular package for speech analysis:
  [2] The Auditory Toolbox for MATLAB by Malcolm Slaney, Version 2
      Interval Research Corporation
      https://engineering.purdue.edu/~malcolm/interval/1998-010/
  Precision: 
  32x32                       32-bit fixed-point input/output data
  f                           Single precision floating-point input/output data
  Input
  objmem                      Memory block allocated for the instance object:
  params                      Parameters of log mel filterbank operation
  spectra[fftSize/2+1]        Fourier image of real signal, positive frequencies only; 
                              Q31 for 32x32
  scaleExp                    Exponent value to scale the Fourier image by a factor 
                              of 2^scaleExp:
                                32x32  For full-scale Q31 real signal the scale exponent
                                       should be set to 15 plus the sum of bit shifts 
                                       applied to data throughout the real-to-complex FFT
                                       transform, as indicated by the respective FFT routine
                                f      For real signal varying in the range [-1,1] the 
                                       scale exponent should be set to 15
  Temporary:
  pScr                        Scratch memory area for the processing function. To 
                              determine the scratch area size, use the respective
                              helper function: logmel<32x32|f>_getScratchSize()
  Output:
  logFbe[mfbBandNum]          Log-scaled filterbank energies; Q6.25 for 32x32
  Restrictions:
  logFbe[],spectra[]          Must not overlap, and must be aligned by 16-bytes
  Fs                          8000 <= Fs <= 48000
  fftSize                     256, 512, 1024 or 2048
  mfbLowFreqQ8, mfbUppFreqQ8  0 <= mfbLowFreqQ8 < mfbUppFreqQ8 <= 16000*256
  mfbBandNum                  0 < mfbBandNum <= 40
-------------------------------------------------------------------------*/

/* Compute log mel filterbank energies */ 
void logmelf_process( logmelf_handle_t handle, void * restrict pScr, float32_t * restrict logFbe, const complex_float * restrict spectra, int scaleExp )
{
    logmelf_t * logmel = (logmelf_t*)handle;
    void * p = pScr;
    float32_t *cplx2mag_scr;
    float32_t *magspec, *fbe;
    int binNum, bandNum;
    NASSERT(logmel && logmel->magic==LOGMELF_MAGIC);
    NASSERT_ALIGN(logFbe, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(spectra, HIFI_SIMD_WIDTH);
    /* Profiler scores should be explicitly reset because this module does not invoke 
     * PROFILE_REPORT macro for them. */
    PROFILE_RESET(logmelf_vec_complex2mag);
    PROFILE_RESET(logmelf_computeFbe);
    PROFILE_RESET(logmelf_vec_logXf);
    binNum = logmel->binUpp-logmel->binLow;
    bandNum = logmel->params.mfbBandNum;
    /* Partition the scratch memory area. */
    magspec      = (float32_t*)ALIGN_PTR(p); p = magspec + binNum;
    fbe          = (float32_t*)ALIGN_PTR(p); p = fbe + bandNum;
    cplx2mag_scr = (float32_t*)ALIGN_PTR(p); p = cplx2mag_scr + binNum; (void)p;
#ifdef _DEBUG
    /* Check that the scratch size is enough to fit all temporary arrays. 
     * This step is skipped for non-debug builds, because _getScratchSize()
     * function may involve time-consuming computations. */
    NASSERT((uint8_t*)p - (uint8_t*)pScr <= (int)logmelf_getScratchSize(&logmel->params));
#endif
    /* Marginal initialization for the memory debugger. */
#if HAVE_VFPU
    XT_SSX2I(XT_CONST_S(0), (xtfloatx2*)magspec, 0);
#else
    XT_SSI(XT_CONST_S(0), (xtfloat*)magspec, 0*sz_f32);
    XT_SSI(XT_CONST_S(0), (xtfloat*)magspec, 1*sz_f32);
#endif
    /* Compute the magnitude spectrum, for the frequencies of interest only. */
    PROFILE_START(logmelf_vec_complex2mag);
    vec_complex2mag_logmel(cplx2mag_scr, magspec, spectra+logmel->binLow, binNum);
    PROFILE_STOP(logmelf_vec_complex2mag);
    PROFILE_START(logmelf_computeFbe);
    logmelf_calcFbe(fbe, magspec, 
                    logmel->weights, logmel->segments, logmel->fbeScales, 
                    binNum, bandNum, scaleExp);
    PROFILE_STOP(logmelf_computeFbe);
    /* Log-scale filterbank energies. */
    PROFILE_START(logmelf_vec_logXf);
    if (LOGMEL_OPT_EQ(logmel->params.opt, FBELOG, NATURAL)) {
        vec_lognf(logFbe, fbe, bandNum);
    } else {
        vec_log10f(logFbe, fbe, bandNum);
    }
    PROFILE_STOP(logmelf_vec_logXf);
} /* logmelf_process() */

/* Returns: size of scratch memory area, in bytes. */
size_t logmelf_getScratchSize( const logmel_params_t * params )
{
    int binLow, binUpp, binNum; /* FFT bin indices; number of FFT bins */
    /* Determine the number of FFT bins in the filterbank frequency range. */
    if (LOGMEL_OPT_EQ(params->opt, MELSCALE, HTK)) {
        logmel_binRange_htk(&binLow, &binUpp, params->mfbLowFreqQ8, params->mfbUppFreqQ8, params->Fs, params->fftSize);
    } else {
        logmel_binRange_auditory(&binLow, &binUpp, params->mfbLowFreqQ8, params->mfbUppFreqQ8, params->Fs, params->fftSize);
    }
    binNum = binUpp-binLow;
    /* Compute the allocation size. */
    return ALIGN_PAD + binNum            *sz_f32 + /* magspec             */
           ALIGN_PAD + XT_MAX(8, binNum) *sz_f32 + /* complex2mag scratch */
           ALIGN_PAD + params->mfbBandNum*sz_f32;  /* fbe                 */
} /* logmelf_getScratchSize() */

void vec_complex2mag_logmel(float32_t * pScr, float32_t * restrict y, const complex_float * restrict x, int N)
{
#if HAVE_VFPU
    const xtfloatx4 * restrict X_rd;
    const xtfloatx4 * restrict Y_rd;
          xtfloatx4 * restrict Y_wr;
    const xtfloatx4 * restrict S_rd;
          xtfloatx2 * restrict S_wr;
    ae_valignx2 X_va;
    int n, tailLen = N&7;
    xtfloatx2 x0, x1, y0, y1, z0, z1, xre, xim;
    ae_int32x2 u0, e0, t0, t1, expx, expy;
    NASSERT(pScr);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(pScr, ALIGN_SIZE);
    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, ALIGN_SIZE);
    S_wr = (xtfloatx2*)pScr;
    Y_wr = (xtfloatx4*)y;
    X_rd = (xtfloatx4*)x;
    X_va = AE_LA128_PP(X_rd);
    __Pragma("loop_count factor=4");
    /* 10 cycles per pipeline stage in steady state with unroll=2 */
    for ( n=0; n<((N>>3)<<2); n++ ) {
        AE_LASX2X2_IP(x0, x1, X_va, X_rd);
        ABS_SX2X2(x0, x1, x0, x1);
        xre = XT_SEL32_HH_SX2(x0, x1);
        xim = XT_SEL32_LL_SX2(x0, x1);
        t0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xre);
        t1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xim);
        expx = AE_MAX32(t0, t1);
        expx = AE_AND32(expx, ~((1<<23)-1));
        expx = AE_MIN32(expx, 254<<23);
        e0 = AE_SUB32(254<<23, expx);
        expy = AE_MAX32(e0, 1<<22);
        y0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(expy);
        MUL_SX2X2(xre, xim, xre, xim, y0, y0);
        MUL_SX2X2(x0, x1, xre, xim, xre, xim);
        x0 = XT_ADD_SX2(x0, x1);
        u0 = AE_MAX32(expx, 1<<22);
        x1 = XT_AE_MOVXTFLOATX2_FROMINT32X2(u0);
        XT_SSX2IP(x0, castxcc(xtfloatx2, S_wr), 2 * sz_f32);
        XT_SSX2IP(x1, castxcc(xtfloatx2, Y_wr), 2*sz_f32);
    } /* n */
    __Pragma("no_reorder");
    S_rd = (xtfloatx4*)pScr;
    Y_rd = (xtfloatx4*)y;
    Y_wr = (xtfloatx4*)y;
    __Pragma("loop_count factor=2");
    /* 28 cycles per pipeline stage in steady state with unroll=2 */
    for ( n=0; n<((N>>3)<<1); n++ ) {
        AE_LSX2X2_IP(x0, x1, S_rd, 4*sz_f32);
        AE_LSX2X2_IP(y0, y1, Y_rd, 4*sz_f32);
        XT_MOVT_SX2(x0, XT_CONST_S(1), XT_OEQ_SX2(x0, plusInff.f));
        XT_MOVT_SX2(x1, XT_CONST_S(1), XT_OEQ_SX2(x1, plusInff.f));
        z0 = XT_FSQRT_SX2(x0); z1 = XT_FSQRT_SX2(x1);
        MUL_SX2X2(y0, y1, y0, y1, z0, z1);
        AE_SSX2X2_IP(y0, y1, Y_wr, 4*sz_f32);
    } /* n */
    if (tailLen>0) {
        xtfloatx2 ALIGN(16) a_x[8];
        float32_t ALIGN(16) a_y[8];
        __Pragma("loop_count min=1");
        for ( n=0; n<tailLen; n++ ) {
            a_x[n] = *(xtfloatx2*)&x[N-tailLen+n];
        }
        __Pragma("loop_count min=1");
        for ( ; n<8; n++ ) {
            a_x[n] = XT_CONST_S(0);
        }
        vec_complex2mag_logmel(pScr, a_y, (complex_float*)a_x, 8);
        __Pragma("loop_count min=1");
        __Pragma("no_simd");
        for ( n=0; n<tailLen; n++ ) {
            y[N-tailLen+n] = a_y[n];
        }
    } /* (tailLen>0) */
#else /* HAVE_VFPU */
    int n;
    const xtfloat* restrict pX = (const xtfloat*)x;
    const xtfloat* restrict pZrd;
          xtfloat* restrict pZwr;
    const xtfloat* restrict pYrd = (const xtfloat*)y;
          xtfloat* restrict pYwr  =(      xtfloat*)y;
    NASSERT(pScr);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(pScr, ALIGN_SIZE);
    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, ALIGN_SIZE);
    pZwr = (xtfloat*)pScr;
    pYwr = (xtfloat*)y;
    /* 41 cycles per pipeline stage in steady state with unroll=2 */
    for (n=0; n<N; n++) {
        xtfloat x0, x1, y0;
        xtfloat xre, xim;
        int t0, t1, nsa;
        int e0;
        int nsa0;

        XT_LSIP(xre,pX,1*sizeof(xtfloat));
        XT_LSIP(xim,pX,1*sizeof(xtfloat));
        xre = XT_ABS_S(xre);
        xim = XT_ABS_S(xim);
        t0 = XT_RFR(xre);
        t1 = XT_RFR(xim);
        nsa = XT_MAX(t0, t1);
        nsa = ((uint32_t)nsa)>> 23;
        nsa = (nsa-127);
        nsa = XT_MIN(nsa, 127);
        e0 = (127-nsa);
        nsa0 = (e0<<23);
        XT_MOVEQZ(nsa0,0x00400000,e0);
        y0 = XT_WFR(nsa0);

        xre = XT_MUL_S(xre, y0);
        xim = XT_MUL_S(xim, y0);

        x0 = XT_MUL_S(xre, xre);
        x1 = XT_MUL_S(xim, xim);

        x0 = XT_ADD_S(x0, x1);
        XT_SSIP(x0,pYwr,sizeof(xtfloat));

        e0 = (127+nsa);
        nsa0 = (e0<<23);
        XT_MOVEQZ(nsa0, 0x00400000, e0);
        x0 = XT_WFR(nsa0);
        XT_SSIP(x0,pZwr,sizeof(xtfloat));
    }
    __Pragma("no_reorder")
    pZrd = (xtfloat*)pScr;
    pYrd = (xtfloat*)y;
    pYwr = (xtfloat*)y;
    /* 43 cycles per pipeline stage in steady state with unroll=2 */
    for (n=0; n<N; n++) {
        xtfloat z0,y0,x0;
        XT_LSIP(y0,pYrd,sizeof(xtfloat));
        z0 = XT_FSQRT_S(y0);
        XT_MOVT_S(z0, plusInff.f, XT_OEQ_S(y0, plusInff.f));
        XT_LSIP(x0,pZrd,sizeof(xtfloat));
        z0 = XT_MUL_S(z0, x0);
        XT_SSIP(z0,pYwr,sizeof(xtfloat));
    }
#endif /* !HAVE_VFPU */
} /* vec_complex2mag_logmel() */

#endif /* HAVE_ */
