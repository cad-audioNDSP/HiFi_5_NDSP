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
    32-bit fixed-point variant
    C code optimized for HiFi5 core
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Log mel filterbank internal definitions. */
#include "logmel_internal.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_audio.h"
#include "NatureDSP_Signal_math.h"
#include "NatureDSP_Signal_vector.h"

#define PROFILE_ENABLE  0 /* If non-zero, measure cycles and print a report to stdout. */
#define ALIGN_SIZE      (HIFI_SIMD_WIDTH)
#define ALIGN_PAD       (ALIGN_SIZE-1)
#define ALIGN_PTR(p)    (void*)(((uintptr_t)(p)+ALIGN_PAD)&~ALIGN_PAD)
#define sz_i16          sizeof(int16_t)
#define sz_i32          sizeof(int32_t)
#define sz_i32c         sizeof(complex_fract32)

#define MAX(a,b)        ((a)>(b) ? (a) : (b))
#define MIN(a,b)        ((a)<(b) ? (a) : (b))

/* Integrated profiler. */
#include "profile.h"

PROFILE_CREATE(logmel32x32_normalizeSpectra);
PROFILE_CREATE(logmel32x32_vec_complex2mag);
PROFILE_CREATE(logmel32x32_computeFbe);
PROFILE_CREATE(logmel32x32_logScaleFbe);

/*-------------------------------------------------------------------------
  Complex magnitude
  Routines compute complex magnitude or its reciprocal

  Precision: 
  f     single precision floating point

  Input:
  x[N]  input complex data
  N     length of vector
  Output:
  y[N]  output data

  Restriction:
  none
-------------------------------------------------------------------------*/
static void vec_complex2mag( fract32 * restrict y, const complex_fract32 * restrict x, int N );

/* Compute base-2 logarithm of each value from the input argument x[N] and store
 * results to the output argument y[N]. Input values should be either positive normalized
 * Q31 numbers from the range [0.5,1), or zero. Results belong to the range [-1,0), Q31,
 * and the left boundary also corresponds to zero input values. */
static void vec_log2( fract32 * restrict y, const fract32 * restrict x, int N );

/* Log-scale filterbank energies.
 * Logarithm base is selected via the input argument base10:
 *   base10==0: logFbe(n) = log(max(1,fbe(n)*2^-fbeExp(n)));
 *   base10!=0: logFbe(n) = log10(max(1,fbe(n)*2^-fbeExp(n)));
 * Input energy estimates must be either zero, or positive normal numbers
 * (i.e. with no redundant sign bits). */
static void logScaleFbe( 
                     fract32 * restrict logFbe, /* [bandNum], Out, Q6.25 */
               const fract32 * restrict fbe,    /* [bandNum], In, Q(31+fbeExp[bandNum]) */
               const int16_t * restrict fbeExp, /* [bandNum], In */
               int bandNum, int base10 );

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
void logmel32x32_process( logmel32x32_handle_t handle, void * restrict pScr, fract32 * restrict logFbe, const complex_fract32 * restrict spectra, int scaleExp )
{
    logmel32x32_t * logmel = (logmel32x32_t*)handle;
    void * p = pScr;
    complex_fract32 *normspec; /* Normalized spectra */
    fract32 *magspec, *fbe; /* Magnitude spectrum, filterbank energies */
    int16_t *fbeExp; /* Filterbank energy exponents */
    int specExp, base10, nsa;
    int binNum, bandNum;
    NASSERT(logmel && logmel->magic==LOGMEL32X32_MAGIC);
    NASSERT_ALIGN(logFbe, ALIGN_SIZE);
    NASSERT_ALIGN(spectra, ALIGN_SIZE);
    /* Profiler scores should be explicitly reset because this module does not invoke 
     * PROFILE_REPORT macro for them. */
    PROFILE_RESET(logmel32x32_normalizeSpectra);
    PROFILE_RESET(logmel32x32_vec_complex2mag);
    PROFILE_RESET(logmel32x32_computeFbe);
    PROFILE_RESET(logmel32x32_logScaleFbe);
    binNum = logmel->binUpp-logmel->binLow;
    bandNum = logmel->params.mfbBandNum;
    /* Partition the scratch memory area. */
    normspec = (complex_fract32*)ALIGN_PTR(p); p = normspec + binNum;
    magspec  = (fract32        *)ALIGN_PTR(p); p = magspec  + binNum;
    fbe      = (fract32        *)ALIGN_PTR(p); p = fbe      + bandNum;
    fbeExp   = (int16_t        *)ALIGN_PTR(p); p = fbeExp   + bandNum; (void)p;
#ifdef _DEBUG
    /* Check that the scratch size is enough to fit all temporary arrays. 
     * This step is skipped for non-debug builds, because _getScratchSize()
     * function may involve time-consuming computations. */
    NASSERT((uint8_t*)p - (uint8_t*)pScr <= (int)logmel32x32_getScratchSize(&logmel->params));
#endif
    /* Marginal initialization for the memory debugger. */
    AE_S32X2_I(AE_ZERO32(), (ae_int32x2*)magspec, 0);
    /* Normally, the input spectrum is either normal or close to normal (block exponent is
     * snall), although frequency bins that belong to the filterbank range should be 
     * re-normalized. */
    PROFILE_START(logmel32x32_normalizeSpectra);
    nsa = vec_bexp32((int32_t*)(spectra+logmel->binLow), 2*binNum);
    vec_shift32x32((int32_t*)normspec, (int32_t*)(spectra+logmel->binLow), nsa, 2*binNum);
    PROFILE_STOP(logmel32x32_normalizeSpectra);
    /* Spectra scale exponent */
    specExp = scaleExp-nsa;
    /* Compute the magnitude spectrum; CQ31 in, Q30 out. */
    PROFILE_START(logmel32x32_vec_complex2mag);
    vec_complex2mag(magspec, normspec, binNum);
    PROFILE_STOP(logmel32x32_vec_complex2mag);
    /* Compute filterbank energies and optionally normalize the results. */
    PROFILE_START(logmel32x32_computeFbe);
    logmel32x32_calcFbe(fbe, fbeExp, magspec, logmel->weights, logmel->segments, 
                        logmel->fbeScales_fract, logmel->fbeScales_exp,
                        binNum, bandNum, specExp);
    PROFILE_STOP(logmel32x32_computeFbe);
    /* Base-10 or natural logarithm scale for filterbank energies? */
    base10 = LOGMEL_OPT_EQ(logmel->params.opt, FBELOG, BASE10);
    /* Log-scale filterbank energies. */
    PROFILE_START(logmel32x32_logScaleFbe);
    logScaleFbe(logFbe, fbe, fbeExp, bandNum, base10);
    PROFILE_STOP(logmel32x32_logScaleFbe);
} /* logmel32x32_process() */

/* Returns: size of scratch memory area, in bytes. */
size_t logmel32x32_getScratchSize( const logmel_params_t * params )
{
    int binLow, binUpp;
    NASSERT(NULL!=params);
    /* Determine the number of FFT bins in the filterbank frequency range. */
    if (LOGMEL_OPT_EQ(params->opt, MELSCALE, HTK)) {
        logmel_binRange_htk(&binLow, &binUpp, params->mfbLowFreqQ8, params->mfbUppFreqQ8, params->Fs, params->fftSize);
    } else {
        logmel_binRange_auditory(&binLow, &binUpp, params->mfbLowFreqQ8, params->mfbUppFreqQ8, params->Fs, params->fftSize);
    }
    /* Compute the allocation size. */
    return ALIGN_PAD + (binUpp-binLow)   *sz_i32c + /* normspec */
           ALIGN_PAD + (binUpp-binLow)   *sz_i32  + /* magspec  */
           ALIGN_PAD + params->mfbBandNum*sz_i32  + /* fbe      */
           ALIGN_PAD + params->mfbBandNum*sz_i16;   /* fbeExp   */
} /* logmel32x32_getScratchSize() */

/*
 * For each complex number in the input argument x[N], compute its magnitude and
 * store the result to respective position in the output argument y[N].
 * Input:
 *   N          Vectors size
 *   x[N]       Input vector, CQ31
 * Output:
 *   y[N]       Output vector, Q30
 * Restrictions:
 *   x[N],y[N]  Must not overlap, and must be 16-bytes aligned
 */

void vec_complex2mag( fract32 * restrict y, const complex_fract32 * restrict x, int N )
{
#if 1
    const ae_int32x4 * restrict pX = (ae_int32x4*)x;
          ae_int32x4 * restrict pY = (ae_int32x4*)y;
          ae_int32x2 * restrict pScr;
    int n;
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT_ALIGN(y, ALIGN_SIZE);
    /* 29 cycles per pipeline stage in steady state with unroll=2 */
    for ( n=0; n<(N>>2); n++ ) {
        ae_int32x2 x0, x1, x2, x3;
        ae_int32x2 f0, f1, s0, s1;
        ae_int16x4 f, g, r, s;
        ae_int64 w0, w1, w2, w3;
        int nsa0, nsa1, nsa2, nsa3;
        /* Q31 */
        AE_L32X2X2_IP(x0, x1, pX, sizeof(*pX));
        AE_L32X2X2_IP(x2, x3, pX, sizeof(*pX));
        /* Q62 = Q31*Q31 */
        w0 = AE_MULZAAD32S_HH_LL(x0, x0);
        w1 = AE_MULZAAD32S_HH_LL(x1, x1);
        w2 = AE_MULZAAD32S_HH_LL(x2, x2);
        w3 = AE_MULZAAD32S_HH_LL(x3, x3);
        /* Data are implicitly shifted right by 1 bit to make the fixed point
         * position odd; then normalized such that the normalization shift
         * amount is even. */
        nsa0 = (AE_NSA64(w0)+1) & ~1;
        nsa1 = (AE_NSA64(w1)+1) & ~1;
        nsa2 = (AE_NSA64(w2)+1) & ~1;
        nsa3 = (AE_NSA64(w3)+1) & ~1;
        /* Q31*2^(2-nsa) = Q(29+nsa) = Q62 + nsa-1 - 32
         * f: signed Q31, 0.25<=f<1 */
        f0 = AE_TRUNCA32F64S(w0, nsa0-1);
        f0 = AE_TRUNCA32F64S_L(f0, w1, nsa1-1);
        f1 = AE_TRUNCA32F64S(w2, nsa2-1);
        f1 = AE_TRUNCA32F64S_L(f1, w3, nsa3-1);
        /*
            A note on sqrt evaluation.
            This variant of algorihm computes almost all in 16-bit domain
            We compute 2 values: r ~ -0.5/sqrt(f) and s ~ sqrt(f), via iterations:
                r = r + r*(1+2*r*s)
                s = s + r*(s*s-f);
            Initial approximation:
                r = 0.5*f-1
                s = 0.658*f+0.352
            We perform the first 2 iterations in 16-bit arithmetic, then the final
            refinement of s is computed in 32-bit precision to achieve the relative
            accuracy of ~4e-8.
        */
        /* Q15 = Q31 - 16 w/ rounding */
        f = AE_ROUND16X4F32SASYM(f0, f1);
        /*** Initial approximations ***/
        /* r <- 0.5*f-1; Q15 */
        r = AE_MULFD16X16X4RAS(AE_MOVDA16(-32768), AE_MOVDA16(16384), AE_MOVDA16(32767), f);
        /* s <- 0.658*f + 0.352; Q15 */
        s = AE_MULFD16X16X4RAS(AE_MOVDA16(21561), AE_MOVDA16(11534), f, AE_MOVDA16(32767));
        /*** First refining iteration ***/
        /* g <- 2*(0.5+r*s); Q15 */
        g = AE_MULFD16X16X4RAS(AE_MOVDA16(-16384), r, AE_MOVDA16(-32768), s);
        g = AE_ADD16(g, g);
        /* r <- r + r*g; Q15 */
        r = AE_MULFD16X16X4RAS(r, r, AE_MOVDA16(32767), g);
        /* g <- s*s-f; Q15 */
        g = AE_MULFD16X16X4RAS(s, AE_MOVDA16(-32768), s, f);
        /* s <- s + r*g; Q15 */
        s = AE_MULFD16X16X4RAS(s, r, AE_MOVDA16(32767), g);
        /*** Second refinement iteration ***/
        g = AE_MULFD16X16X4RAS(AE_MOVDA16(-16384), r, AE_MOVDA16(-32768), s);
        g = AE_ADD16(g, g);
        r = AE_MULFD16X16X4RAS(r, r, AE_MOVDA16(32767), g);
        g = AE_MULFD16X16X4RAS(s, AE_MOVDA16(-32768), s, f);
        s = AE_MULFD16X16X4RAS(s, r, AE_MOVDA16(32767), g);
        /*** Final refinement iteration ***/
        /* f <- f - s*s; Q31 = Q31 - [Q15*Q15 + 1] */
        AE_MULSF16X4SS(f0, f1, s, s);
        /* s <- s - f*r; Q31 <- [Q15 + 16] - [Q31*Q15 - 15] */
        AE_CVTI32X4F16(s0, s1, s, 16);
        AE_MULSF2P32X16X4S(s0, s1, f0, f1, r);
        /*** Format and save the results. ***/
        /* Q30 <- Q31*2^(1-nsa/2) - 1 w/ rounding */
        s0 = AE_SRAV32RS(s0, AE_SRAI32(AE_MOVDA32X2(nsa0, nsa1), 1));
        s1 = AE_SRAV32RS(s1, AE_SRAI32(AE_MOVDA32X2(nsa2, nsa3), 1));
        AE_S32X2X2_IP(s0, s1, pY, sizeof(*pY));
    } /* n */
    if (N&3) {
        ae_int32x2 x0, x1, x2;
        ae_int32x2 f0, f1, s0, s1;
        ae_int16x4 f, g, r, s;
        ae_int64 w0, w1, w2;
        int nsa0, nsa1, nsa2;
        int32_t ALIGN(16) aScr[3*2];
        pScr = (ae_int32x2*)aScr;
        /* Initialize scratch for the memory debugger. */
        AE_S32X2X2_I(AE_ZERO32(), AE_ZERO32(), (ae_int32x4*)pScr, 0);
        AE_S32X2_I(AE_ZERO32(), pScr, 2*sz_i32c);
        __Pragma("loop_count min=1");
        for ( n=0; n<(N&3); n++ ) {
            AE_L32X2_IP(x0, castxcc(ae_int32x2, pX), sz_i32c);
            AE_S32X2_IP(x0, pScr, sz_i32c);
        }
        pScr = (ae_int32x2*)XT_ADDX8(-(N&3), (uintptr_t)pScr);
        /* Q31 */
        AE_L32X2X2_I(x0, x1, (ae_int32x4*)pScr, 0);
        x2 = AE_L32X2_I(pScr, 2*sz_i32c);
        /* Q62 = Q31*Q31 */
        w0 = AE_MULZAAD32S_HH_LL(x0, x0);
        w1 = AE_MULZAAD32S_HH_LL(x1, x1);
        w2 = AE_MULZAAD32S_HH_LL(x2, x2);
        /* Data are implicitly shifted right by 1 bit to make the fixed point
         * position odd; then normalized such that the normalization shift
         * amount is even. */
        nsa0 = (AE_NSA64(w0)+1) & ~1;
        nsa1 = (AE_NSA64(w1)+1) & ~1;
        nsa2 = (AE_NSA64(w2)+1) & ~1;
        /* Q31*2^(2-nsa) = Q(29+nsa) = Q62 + nsa-1 - 32
         * f: signed Q31, 0.25<=f<1 */
        f0 = AE_TRUNCA32F64S(w0, nsa0-1);
        f0 = AE_TRUNCA32F64S_L(f0, w1, nsa1-1);
        f1 = AE_TRUNCA32F64S(w2, nsa2-1);
        /* Q15 = Q31 - 16 w/ rounding */
        f = AE_ROUND16X4F32SASYM(f0, f1);
        /*** Initial approximations ***/
        /* r <- 0.5*f-1; Q15 */
        r = AE_MULFD16X16X4RAS(AE_MOVDA16(-32768), AE_MOVDA16(16384), AE_MOVDA16(32767), f);
        /* s <- 0.658*f + 0.352; Q15 */
        s = AE_MULFD16X16X4RAS(AE_MOVDA16(21561), AE_MOVDA16(11534), f, AE_MOVDA16(32767));
        /*** First refining iteration ***/
        /* g <- 2*(0.5+r*s); Q15 */
        g = AE_MULFD16X16X4RAS(AE_MOVDA16(-16384), r, AE_MOVDA16(-32768), s);
        g = AE_ADD16(g, g);
        /* r <- r + r*g; Q15 */
        r = AE_MULFD16X16X4RAS(r, r, AE_MOVDA16(32767), g);
        /* g <- s*s-f; Q15 */
        g = AE_MULFD16X16X4RAS(s, AE_MOVDA16(-32768), s, f);
        /* s <- s + r*g; Q15 */
        s = AE_MULFD16X16X4RAS(s, r, AE_MOVDA16(32767), g);
        /*** Second refinement iteration ***/
        g = AE_MULFD16X16X4RAS(AE_MOVDA16(-16384), r, AE_MOVDA16(-32768), s);
        g = AE_ADD16(g, g);
        r = AE_MULFD16X16X4RAS(r, r, AE_MOVDA16(32767), g);
        g = AE_MULFD16X16X4RAS(s, AE_MOVDA16(-32768), s, f);
        s = AE_MULFD16X16X4RAS(s, r, AE_MOVDA16(32767), g);
        /*** Final refinement iteration ***/
        /* f <- f - s*s; Q31 = Q31 - [Q15*Q15 + 1] */
        AE_MULSF16X4SS(f0, f1, s, s);
        /* s <- s - f*r; Q31 <- [Q15 + 16] - [Q31*Q15 - 15] */
        AE_CVTI32X4F16(s0, s1, s, 16);
        AE_MULSF2P32X16X4S(s0, s1, f0, f1, r);
        /*** Format and save the results. ***/
        /* Q30 <- Q31*2^(1-nsa/2) - 1 w/ rounding */
        s0 = AE_SRAV32RS(s0, AE_SRAI32(AE_MOVDA32X2(nsa0, nsa1), 1));
        s1 = AE_SRAV32RS(s1, AE_SRAI32(AE_MOVDA32(nsa2), 1));
        AE_S32X2X2_I(s0, s1, (ae_int32x4*)pScr, 0);
        __Pragma("loop_count min=1");
        for ( n=0; n<(N&3); n++ ) {
            AE_L32_IP(s0, castxcc(ae_int32, pScr), sz_i32);
            AE_S32_L_IP(s0, castxcc(ae_int32, pY), sz_i32);
        } /* n */
    } /* (N&3) */
#else /* Functional stub */
    int n;
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT_ALIGN(y, ALIGN_SIZE);
    for ( n=0; n<N; n++ ) {
        float64_t re, im, mag;
        re = ldexp((float64_t)x[n].s.re, -31);
        im = ldexp((float64_t)x[n].s.im, -31);
        mag = sqrt(re*re+im*im);
        y[n] = (fract32)MAX(MIN_INT32, MIN(MAX_INT32, round(ldexp(mag, 30))));
    } /* n */
#endif
} /* vec_complex2mag() */ 

/*
 * Compute base-2 logarithm of each value from the input argument x[N] and store
 * results to the output argument y[N]. Input values should be either positive normalized
 * Q31 numbers from the range [0.5,1), or zero. Results belong to the range [-1,0), Q31,
 * and the left boundary also corresponds to zero input values.
 * 
 * Input:
 *   N          Vectors size
 *   x[N]       Input vector, Q31
 * Output:
 *   y[N]       Output vector, Q31
 * Restrictions:
 *   N          Must be positive
 *   x[N],y[N]  Must be 16-bytes aligned
 */

void vec_log2( fract32 * restrict y, const fract32 * restrict x, int N )
{
#if 1
    /*
     Matlab code for computing the polynomial:
     x=(sqrt(0.5):pow2(1,-16):sqrt(2));
     z=1-x;
     y=log(x)./z;
     p=polyfit(z,y,8);
    */

    static const int32_t ALIGN(32) polylog_tbl[9] = {
        -161926367,-273781379,-283444439,-304997874,-356837178,-429521585,-536898174,-715827933,-1073741641
    };
    const ae_int32x4 * restrict pX0  = (ae_int32x4*)x;
    const ae_int32x4 * restrict pX1  = (ae_int32x4*)x;
          ae_int32x4 * restrict pY   = (ae_int32x4*)y;
    const ae_int32   * restrict pTbl = (ae_int32*)polylog_tbl;
    const ae_int32   * restrict p_rd;
          ae_int32   * restrict p_wr;
    int n, tailLen = N&7;
    NASSERT(N>0);
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    __Pragma("loop_count factor=2");
    /* 26 cycles per pipeline stage in steady state with unroll=2 */
    for ( n=0; n<(N>>3)<<1; n++ ) {
        ae_int32x2 c0, c1, c2, c3, c4, c5, c6, c7, c8; /* Polynomial coeffs */
        ae_int32x2 x0, x1, y0, y1; /* Input data; approximation results */
        ae_int32x2 f0, f1, g0, g1; /* Polynomial arguments; polynomial results */
        ae_int32x2 h00, h01, h10, h11; /* Polynomial accumulators */
        ae_int32x2 h20, h21, h30, h31;
        ae_int32x2 h40, h41, h50, h51;
        ae_int32x2 h60, h61, h70, h71;
        ae_int32x2 h80, h81;
        ae_int32x2 t0, t1;
        AE_L32X2X2_IP(x0, x1, pX0, 4*sz_i32);
        /* To fit 0.707...1.414 */
        f0 = x0; AE_MOVT32X2(f0, AE_MULP32X2(x0, 2), AE_LT32(x0, (int32_t)1518500250L));
        f1 = x1; AE_MOVT32X2(f1, AE_MULP32X2(x1, 2), AE_LT32(x1, (int32_t)1518500250L));
        f0 = AE_SUB32(MIN_INT32, f0); f1 = AE_SUB32(MIN_INT32, f1);
        /* Here _XP provides a better schedule than _IP. */
        AE_L32_XP(c0, pTbl, sz_i32); AE_L32_XP(c1, pTbl, sz_i32);
        AE_L32_XP(c2, pTbl, sz_i32); AE_L32_XP(c3, pTbl, sz_i32);
        AE_L32_XP(c4, pTbl, sz_i32); AE_L32_XP(c5, pTbl, sz_i32);
        AE_L32_XP(c6, pTbl, sz_i32); AE_L32_XP(c7, pTbl, sz_i32);
        AE_L32_XP(c8, pTbl, -8*(int)sz_i32); h80 = c8;
        /* We use simple assignment for those accumulators that do not
         * actually require register moves. */
        h00 = c0; h01 = c0;
        h10 = c1; h20 = c2; AE_MOVD32X4(h11, h21, c1, c2);
        h30 = c3; h40 = c4; AE_MOVD32X4(h31, h41, c3, c4);
        h50 = c5; h60 = c6; AE_MOVD32X4(h51, h61, c5, c6);
        h70 = c7; h80 = c8; AE_MOVD32X4(h71, h81, c7, c8);
        AE_MULAF2P32X4RAS(h10, h11, h00, h01, f0, f1);
        AE_MULAF2P32X4RAS(h20, h21, h10, h11, f0, f1);
        AE_MULAF2P32X4RAS(h30, h31, h20, h21, f0, f1);
        AE_MULAF2P32X4RAS(h40, h41, h30, h31, f0, f1);
        AE_MULAF2P32X4RAS(h50, h51, h40, h41, f0, f1);
        AE_MULAF2P32X4RAS(h60, h61, h50, h51, f0, f1);
        AE_MULAF2P32X4RAS(h70, h71, h60, h61, f0, f1);
        AE_MULAF2P32X4RAS(h80, h81, h70, h71, f0, f1);
        AE_MULF2P32X4RAS(t0, t1, f0, f1, f0, f1);  
        AE_MULF2P32X4RAS(g0, g1, t0, t1, h80, h81);
        AE_MULAF2P32X4RAS(g0, g1, f0, f1, MIN_INT32, MIN_INT32);
        AE_MUL2P32X4S(g0, g1, g0, g1, 2, 2);
        AE_MULF2P32X4RAS(g0, g1, g0, g1, (int32_t)1549082005L, (int32_t)1549082005L);
        /* Reload data to break the dependency path. */
        AE_L32X2X2_IP(x0, x1, pX1, sizeof(ae_int32x4));
        /* To fit 0.707...1.414 */
        y0 = g0; AE_MOVT32X2(y0, AE_XOR32(y0, MIN_INT32), AE_LT32(x0, (int32_t)1518500250L));
        y1 = g1; AE_MOVT32X2(y1, AE_XOR32(y1, MIN_INT32), AE_LT32(x1, (int32_t)1518500250L));
        AE_S32X2X2_IP(y0, y1, pY, sizeof(ae_int32x4));
    } /* n */
    if (tailLen>0) {
        int32_t ALIGN(16) a_x[8];
        int32_t ALIGN(16) a_y[8];
        /* Initialize the local array for the memory debugger. */
        AE_S32X2X2_I(MAX_INT32, MAX_INT32, (ae_int32x4*)a_x, 0*4*sz_i32);
        AE_S32X2X2_I(MAX_INT32, MAX_INT32, (ae_int32x4*)a_x, 1*4*sz_i32);
        p_rd = (ae_int32*)XT_ADDX4(N-tailLen, (uintptr_t)x); 
        p_wr = (ae_int32*)a_x;
        __Pragma("no_simd");
        __Pragma("loop_count min=1");
        for ( n=0; n<tailLen; n++ ) *p_wr++ = *p_rd++;
        vec_log2(a_y, a_x, 8);
        p_rd = (ae_int32*)a_y; 
        p_wr = (ae_int32*)XT_ADDX4(N-tailLen, (uintptr_t)y); 
        __Pragma("no_simd");
        __Pragma("loop_count min=1");
        for ( n=0; n<tailLen; n++ ) *p_wr++ = *p_rd++;
    } /* (tailLen>0) */
#else /* Functional stub */
    int n;
    NASSERT(N>0);
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    for ( n=0; n<N; n++ ) {
        if (x[n]>0) {
            float64_t f = ldexp((float64_t)x[n], -31);
            NASSERT(S_exp0_l(x[n])==0);
            y[n] = (fract32)MAX(MIN_INT32, round(ldexp(log2(f), 31)));
        } else {
            y[n] = MIN_INT32;
        }
    } /* n */
#endif
} /* vec_log2() */

/*
 * Log-scale filterbank energies.
 * Logarithm base is selected via the input argument base10:
 *   base10==0: logFbe(n) = log(max(1,fbe(n)*2^-fbeExp(n)));
 *   base10!=0: logFbe(n) = log10(max(1,fbe(n)*2^-fbeExp(n)));
 * Input energy estimates must be either zero, or positive normal numbers
 * (i.e. with no redundant sign bits). 
 * Input:
 *   bandNum          Number of filterbank bands
 *   base10           If non-zero, use base-10 logarithm, otherwise use the natural logarithm
 *   fbe[bandNum]     Filterbank energies, Q(31+fbeExp[bandNum])
 *   fbeExp[bandNum]  Energy exponents 
 * Output:
 *   logFbe[bandNum]  Log-scaled filterbank energies, Q6.25
 * Restrictions:
 *   bandNum                  Must be positive
 *   logFbe[],fbe[],fbeExp[]  Must be 16-bytes aligned, must not overlap
 *   fbe[]                    Either positive normalized, or zero
 */

void logScaleFbe( fract32 * restrict logFbe, 
            const fract32 * restrict fbe, 
            const int16_t * restrict fbeExp, 
            int bandNum, int base10 )
{
    const ae_int32x4 * restrict p_logFbe_rd;
          ae_int32x4 * restrict p_logFbe_wr;
    const ae_int16x4 * restrict p_fbeExp;
    const ae_int32   * restrict p32_r;
          ae_int32   * restrict p32_w;
    const ae_int16   * restrict p16_r;
          ae_int16   * restrict p16_w;
    ae_int32x2 vlogFbe0, vlogFbe1, vfbeExp0, vfbeExp1, vg0, vg1, vlog;
    ae_int16x4 vfbeExp_fr16;
    fract32 f;
    int n, tailLen = bandNum&3;

    NASSERT(bandNum>0);
    NASSERT_ALIGN(logFbe, ALIGN_SIZE);
    NASSERT_ALIGN(fbe, ALIGN_SIZE);
    NASSERT_ALIGN(fbeExp, ALIGN_SIZE);

    /* logFbe[n] <- log2(fbe[n]); Q31 in, Q31 out. */
    vec_log2(logFbe, fbe, bandNum);

    f = 646456993; /* log10(2), Q31 */
    XT_MOVEQZ(f, 1488522236, base10); /* log(2), Q31 */
    vlog = f;

    p_logFbe_rd = (ae_int32x4*)logFbe;
    p_logFbe_wr = (ae_int32x4*)logFbe;
    p_fbeExp    = (ae_int16x4*)fbeExp;

    /* 5 cycles per pipeline stage in steady state with unroll=2 */
    for ( n=0; n<(bandNum>>2); n++ ) {
        /* g <- log2(fbe(n)) - fbeExp; Q6.25 <- Q31 - 6 */
        AE_L32X2X2_IP(vlogFbe0, vlogFbe1, p_logFbe_rd, 4*sz_i32);
        AE_L16X4_IP(vfbeExp_fr16, p_fbeExp, 4*sz_i16);
        AE_CVTI32X4F16S(vfbeExp0, vfbeExp1, vfbeExp_fr16, 25);
        /* logFbe <- logFbe*2^-6 w/ rounding */
        AE_MULF2P32X4RAS(vlogFbe0, vlogFbe1, vlogFbe0, vlogFbe1, 1<<(31-6), 1<<(31-6));
        /* g <- MAX(logFbe-fbeExp, 0) */
        AE_MULAF2P32X4RAS(vlogFbe0, vlogFbe1, vfbeExp0, vfbeExp1, MIN_INT32, MIN_INT32);
        vg0 = AE_MAX32(vlogFbe0, 0);
        vg1 = AE_MAX32(vlogFbe1, 0);
        /* Q25 <- Q31*Q25 - 31 w/ rounding */
        AE_MULF2P32X4RAS(vg0, vg1, vg0, vg1, vlog, vlog);
        AE_S32X2X2_IP(vg0, vg1, p_logFbe_wr, 4*sz_i32);
    } /* n */

    if (tailLen>0) {
        int32_t ALIGN(16) a_logFbe[4];
        int16_t ALIGN(8) a_fbeExp[4];
        /* Initialize local arrays for the memory debugger. */
        AE_S32X2X2_I(0, 0, (ae_int32x4*)a_logFbe, 0);
        AE_S16X4_I(0, (ae_int16x4*)a_fbeExp, 0);
        p32_r = (ae_int32*)XT_ADDX4(bandNum-tailLen, (uintptr_t)logFbe);
        p32_w = (ae_int32*)a_logFbe;
        p16_r = (ae_int16*)XT_ADDX2(bandNum-tailLen, (uintptr_t)fbeExp);
        p16_w = (ae_int16*)a_fbeExp;
        __Pragma("no_simd");
        __Pragma("loop_count min=1");
        for ( n=0; n<tailLen; n++ ) {*p32_w++ = *p32_r++; *p16_w++ = *p16_r++;}
        __Pragma("no_reorder");
        /* g <- log2(fbe(n)) - fbeExp; Q6.25 <- Q31 - 6 */
        AE_L32X2X2_I(vlogFbe0, vlogFbe1, (ae_int32x4*)a_logFbe, 0);
        vfbeExp_fr16 = AE_L16X4_I((ae_int16x4*)a_fbeExp, 0);
        AE_CVTI32X4F16S(vfbeExp0, vfbeExp1, vfbeExp_fr16, 25);
        /* logFbe <- logFbe*2^-6 w/ rounding */
        AE_MULF2P32X4RAS(vlogFbe0, vlogFbe1, vlogFbe0, vlogFbe1, 1<<(31-6), 1<<(31-6));
        /* g <- MAX(logFbe-fbeExp, 0) */
        AE_MULAF2P32X4RAS(vlogFbe0, vlogFbe1, vfbeExp0, vfbeExp1, MIN_INT32, MIN_INT32);
        vg0 = AE_MAX32(vlogFbe0, 0);
        vg1 = AE_MAX32(vlogFbe1, 0);
        /* Q25 <- Q31*Q25 - 31 w/ rounding */
        AE_MULF2P32X4RAS(vg0, vg1, vg0, vg1, vlog, vlog);
        AE_S32X2X2_I(vg0, vg1, (ae_int32x4*)a_logFbe, 0);
        __Pragma("no_reorder");
        p32_r = (ae_int32*)a_logFbe;
        p32_w = (ae_int32*)XT_ADDX4(bandNum-tailLen, (uintptr_t)logFbe);
        __Pragma("no_simd");
        __Pragma("loop_count min=1");
        for ( n=0; n<tailLen; n++ ) *p32_w++ = *p32_r++;
    } /* (tailLen>0) */
} /* logScaleFbe() */
