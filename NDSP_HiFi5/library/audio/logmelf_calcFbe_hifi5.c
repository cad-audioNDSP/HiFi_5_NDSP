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
#include "common_fpu.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_math.h"
/* Log mel filterbank internal definitions. */
#include "logmel_internal.h"

#if HAVE_FPU || HAVE_VFPU

#define USE_REFERENCE_CODE  0
#define ALIGN_SIZE          (HIFI_SIMD_WIDTH)
#define sz_i16              sizeof(int16_t)
#define sz_f32              sizeof(float32_t)

/* 
 * MATLAB reference:
 *   fbe = zeros(bandNum,1);
 *   n = 1;
 *   for p=0:logMel.segments(1)-1
 *     w = logMel.weights(n+p); 
 *     fbe(1) = fbe(1)+w*magspec(n+p);
 *   end
 *   n = n+logMel.segments(1);
 *   for m=1:bandNum-1
 *     for p=0:logMel.segments(m+1)-1
 *       w = logMel.weights(n+p);
 *       fbe(m) = fbe(m)+(1-w)*magspec(n+p);
 *       fbe(m+1) = fbe(m+1)+w*magspec(n+p);
 *     end
 *     n = n+logMel.segments(m+1);
 *   end
 *   for p=0:logMel.segments(bandNum+1)-1
 *     w = logMel.weights(n+p);
 *     fbe(bandNum) = fbe(bandNum)+w*magspec(n+p);
 *   end
 *   % FBE normalization
 *   fbe = pow2(fbe,scaleExp).*logMel.fbeScales(:);
 */

#if USE_REFERENCE_CODE
#include <math.h>
#define MAX(a,b)  ((a)>(b) ? (a) : (b))
#endif

/* Immediate options for inlined kernel routine */
#define CFK_SKIP_NEG_NO    0  /* Do not update the negative-slope accumulator */
#define CFK_SKIP_NEG_YES   1  /* Update the negative-slope accumulator */

static void ATTRIBUTE_ALWAYS_INLINE calcFbe_Kernel( float32_t *  p_negAcc,
                                                    float32_t *  p_posAcc,
                                              const float32_t ** pp_magspec,
                                              const float32_t ** pp_weights,
                                              int segmentLen, int imm_skip_neg_acc );

static void calcFbe_Scale(
                      float32_t * restrict fbe,       /* In/Out */
                const float32_t * restrict fbeScales, /* In, optional */
                int bandNum, int specExp );

/* Compute filterbank energies from the magnitude spectrum, and optionally normalize the results.
 * Input:
 *   binNum                    Number of magnitude spectrum bins
 *   bandNum                   Number of filterbank bands
 *   specExp                   Exponent of the magnitude spectrum
 *   magspec[binNum]           Magnitude spectrum, Q30*2^specExp (32x32)
 *   weights[binNum]           Filterbank weights, Q31 (32x32)
 *   segments[bandNum+1]       Segmentation of spectrum bins w.r.t. filter bands
 *   fbeScales_fract[bandNum]  (32x32) FBE normalization factors, or NULL if normalization is disabled
 *   fbeScales_exp[bandNum]    (32x32) Exponent of normalization factors:
 *                                fbeScales[k] = fbeScales_fract[k]*2^-fbeScales_exp[k]
 *   fbeScales[bandNum]        (f) FBE normalization factors
 * Output:
 *   fbe[bandNum]              Filterbank energies, Q(31+fbeExp[bandNum]) (32x32)
 *   fbeExp[bandNum]           (32x32) Energy exponents 
 * Restrictions:
 *   bandNum                   Must be positive
 *   fbe[],fbeExp[],magspec[],weights[],segments[], fbeScales_fract[],fbeScales_exp[], fbeScales[]
 *                             Must be 16-bytes aligned, must not overlap. */

/* Single precision floating-point variant */
void logmelf_calcFbe( 
                float32_t * restrict fbe,
          const float32_t * restrict magspec,
          const float32_t * restrict weights,
          const int16_t   * restrict segments,
          const float32_t * restrict fbeScales,
          int binNum, int bandNum, int specExp )
{
    xtfloat * restrict p_fbe = (xtfloat*)fbe;
    const float32_t * p_magspec = magspec;
    const float32_t * p_weights = weights;
    float32_t negAcc, posAcc; 
    int m, segmentLen;
    NASSERT(bandNum>0);
    NASSERT_ALIGN(fbe, ALIGN_SIZE);
    NASSERT_ALIGN(magspec, ALIGN_SIZE);
    NASSERT_ALIGN(weights, ALIGN_SIZE);
    NASSERT_ALIGN(segments, ALIGN_SIZE);
    NASSERT_ALIGN(fbeScales, ALIGN_SIZE);
    /* Positive-slope side of the left-most triangle. */
    segmentLen = *segments++; posAcc = 0;
    calcFbe_Kernel(NULL, (float32_t*)&posAcc, &p_magspec, &p_weights, segmentLen, CFK_SKIP_NEG_YES);
    /* Because of 1/2 bands overlapping, each inner segment simultaneously updates
     * a pair of bands: negative-slope side for the left triangle, positive-slope 
     * side for the right triangle. */
    for ( m=0; m<bandNum-1; m++ ) {
        segmentLen = *segments++; negAcc = posAcc; posAcc = 0;
        calcFbe_Kernel((float32_t*)&negAcc, (float32_t*)&posAcc, &p_magspec, &p_weights, segmentLen, CFK_SKIP_NEG_NO);
        *p_fbe++ = negAcc;
    } /* m */
    /* Negative-slope side of the right-most triangle. */
    segmentLen = *segments++;
    calcFbe_Kernel(NULL, (float32_t*)&posAcc, &p_magspec, &p_weights, segmentLen, CFK_SKIP_NEG_YES);
    *p_fbe = posAcc;
    /* Apply spectra scale factor, energy scale factors and the lower limit. */
    calcFbe_Scale(fbe, fbeScales, bandNum, specExp);
} /* logmelf_calcFbe() */

void calcFbe_Kernel( float32_t *  p_negAcc,
                     float32_t *  p_posAcc,
               const float32_t ** pp_magspec,
               const float32_t ** pp_weights,
               int segmentLen, int imm_skip_neg_acc )
{
#if USE_REFERENCE_CODE
    const float32_t * p_magspec = *pp_magspec;
    const float32_t * p_weights = *pp_weights;
    float32_t wp, wn, s; /* Positive-slope weight; negative-slope weight; spectrum sample */
    float32_t negAcc, posAcc;
    int n;
    NASSERT(imm_skip_neg_acc || p_negAcc);
    if (!imm_skip_neg_acc) {
        negAcc = *p_negAcc; posAcc = *p_posAcc;
        for ( n=0; n<segmentLen; n++ ) {
            s = *p_magspec++; wp = *p_weights++; wn = 1.f-wp;
            negAcc += wn*s; posAcc += wp*s;
        } /* n */
    } else { /* (imm_skip_neg_acc) */
        posAcc = *p_posAcc;
        for ( n=0; n<segmentLen; n++ ) {
            s = *p_magspec++; wp = *p_weights++;
            posAcc += wp*s;
        } /* n */
    } /* (imm_skip_neg_acc) */
    if (!imm_skip_neg_acc) {
        *p_negAcc = negAcc;
    }
    *p_posAcc = posAcc;
    *pp_magspec = p_magspec;
    *pp_weights = p_weights;
#elif HAVE_VFPU
    const xtfloatx4 * p_magspec;
    const xtfloatx4 * p_weights;
    ae_valignx2 magspec_va2, weights_va2;
    float32_t negAcc=0.f, posAcc;
    int n, headLen;
    NASSERT(imm_skip_neg_acc || p_negAcc);
    if (!imm_skip_neg_acc) {
        if (segmentLen>=8) {
            static const int32_t ALIGN(16) seq[] = {0,1,2,3,4,5,6,7};
            xtbool2 b0, b1, b2, b3;
            ae_int32x2 seq0, seq1, seq2, seq3;
            xtfloatx2 s0, s1, s2, s3;
            xtfloatx2 wp0, wp1, wp2, wp3;
            xtfloatx2 wn0, wn1, wn2, wn3;
            xtfloatx2 negAcc0, negAcc1, negAcc2, negAcc3, t;
            xtfloatx2 posAcc0, posAcc1, posAcc2, posAcc3;
            AE_L32X2X2_I(seq0, seq1, (const ae_int32x4*)seq, 0*sizeof(ae_int32x4));
            AE_L32X2X2_I(seq2, seq3, (const ae_int32x4*)seq, 1*sizeof(ae_int32x4));
            headLen = ((segmentLen-1)&7)+1;
            b0 = AE_LT32(seq0, headLen); /* Mask useful elements for the prologue. */
            b1 = AE_LT32(seq1, headLen);
            b2 = AE_LT32(seq2, headLen);
            b3 = AE_LT32(seq3, headLen);
            negAcc0 = *(xtfloat*)p_negAcc; negAcc1 = negAcc2 = negAcc3 = 0;
            posAcc0 = *(xtfloat*)p_posAcc; posAcc1 = posAcc2 = posAcc3 = 0;
            negAcc0 = AE_SEL32_HH_SX2(negAcc0, negAcc1); /* Zero the LSW! */
            posAcc0 = AE_SEL32_HH_SX2(posAcc0, posAcc1);
            p_magspec = (xtfloatx4*)*pp_magspec;
            p_weights = (xtfloatx4*)*pp_weights;
            magspec_va2 = AE_LA128_PP(p_magspec);
            weights_va2 = AE_LA128_PP(p_weights);
            AE_LASX2X2_IP(wp0, wp1, weights_va2, p_weights);
            AE_LASX2X2_IP(wp2, wp3, weights_va2, p_weights);
            AE_LASX2X2_IP(s0, s1, magspec_va2, p_magspec);
            AE_LASX2X2_IP(s2, s3, magspec_va2, p_magspec);
            XT_MOVF_SX2(s0, 0.f, b0); XT_MOVF_SX2(s1, 0.f, b1);
            XT_MOVF_SX2(s2, 0.f, b2); XT_MOVF_SX2(s3, 0.f, b3);
            SUB_SX2X2(wn0, wn1, 1.f, 1.f, wp0, wp1);
            SUB_SX2X2(wn2, wn3, 1.f, 1.f, wp2, wp3);
            MADD_SX2X2(negAcc0, negAcc1, wn0, wn1, s0, s1);
            MADD_SX2X2(negAcc2, negAcc3, wn2, wn3, s2, s3);
            MADD_SX2X2(posAcc0, posAcc1, wp0, wp1, s0, s1);
            MADD_SX2X2(posAcc2, posAcc3, wp2, wp3, s2, s3);
            segmentLen -= headLen;
            *pp_magspec = (float32_t*)XT_ADDX4(headLen, (uintptr_t)*pp_magspec);
            *pp_weights = (float32_t*)XT_ADDX4(headLen, (uintptr_t)*pp_weights);
            magspec_va2 = AE_LA128_PP(*pp_magspec);
            weights_va2 = AE_LA128_PP(*pp_weights);
            /* 4 cycles per pipeline stage in steady state with unroll=1 */
            for ( n=0; n<(segmentLen>>3); n++ ) {
                AE_LASX2X2_IP(wp0, wp1, weights_va2, castxcc(xtfloatx4, *pp_weights));
                AE_LASX2X2_IP(wp2, wp3, weights_va2, castxcc(xtfloatx4, *pp_weights));
                AE_LASX2X2_IP(s0, s1, magspec_va2, castxcc(xtfloatx4, *pp_magspec));
                AE_LASX2X2_IP(s2, s3, magspec_va2, castxcc(xtfloatx4, *pp_magspec));
                SUB_SX2X2(wn0, wn1, 1.f, 1.f, wp0, wp1);
                SUB_SX2X2(wn2, wn3, 1.f, 1.f, wp2, wp3);
                MADD_SX2X2(negAcc0, negAcc1, wn0, wn1, s0, s1);
                MADD_SX2X2(negAcc2, negAcc3, wn2, wn3, s2, s3);
                MADD_SX2X2(posAcc0, posAcc1, wp0, wp1, s0, s1);
                MADD_SX2X2(posAcc2, posAcc3, wp2, wp3, s2, s3);
            } /* n */
            ADD_SX2X2(negAcc0, negAcc1, negAcc0, negAcc1, negAcc2, negAcc3);
            ADD_SX2X2(posAcc0, posAcc1, posAcc0, posAcc1, posAcc2, posAcc3);
            ADD_SX2X2(negAcc0, posAcc0, negAcc0, posAcc0, negAcc1, posAcc1);
            t = XT_ADD_SX2(AE_SEL32_HH_SX2(negAcc0, posAcc0), AE_SEL32_LL_SX2(negAcc0, posAcc0));
            negAcc = XT_HIGH_S(t); posAcc = XT_LOW_S(t);
        } else { /* (segmentLen<8) */
            negAcc = *p_negAcc; posAcc = *p_posAcc;
            /* 4 cycles per pipeline stage in steady state with unroll=1 */
            for ( n=0; n<segmentLen; n++ ) {
                xtfloat wn, wp, s;
                AE_LSIP(wp, castxcc(xtfloat, *pp_weights), sz_f32);
                AE_LSIP(s, castxcc(xtfloat, *pp_magspec), sz_f32);
                wn = XT_SUB_S(1.f, wp);
                XT_MADD_S(negAcc, wn, s);
                XT_MADD_S(posAcc, wp, s);
            } /* n */
        } /* (segmentLen<8) */
    } else { /* (imm_skip_neg_acc) */
        if (segmentLen>=8) {
            static const int32_t ALIGN(16) seq[] = {0,1,2,3,4,5,6,7};
            xtbool2 b0, b1, b2, b3;
            ae_int32x2 seq0, seq1, seq2, seq3;
            xtfloatx2 s0, s1, s2, s3;
            xtfloatx2 wp0, wp1, wp2, wp3;
            xtfloatx2 posAcc0, posAcc1, posAcc2, posAcc3;
            AE_L32X2X2_I(seq0, seq1, (const ae_int32x4*)seq, 0*sizeof(ae_int32x4));
            AE_L32X2X2_I(seq2, seq3, (const ae_int32x4*)seq, 1*sizeof(ae_int32x4));
            headLen = ((segmentLen-1)&7)+1;
            b0 = AE_LT32(seq0, headLen); /* Mask useful elements for the prologue. */
            b1 = AE_LT32(seq1, headLen);
            b2 = AE_LT32(seq2, headLen);
            b3 = AE_LT32(seq3, headLen);
            posAcc0 = *(xtfloat*)p_posAcc; posAcc1 = posAcc2 = posAcc3 = 0;
            posAcc0 = AE_SEL32_HH_SX2(posAcc0, posAcc1);
            p_magspec = (xtfloatx4*)*pp_magspec;
            p_weights = (xtfloatx4*)*pp_weights;
            magspec_va2 = AE_LA128_PP(p_magspec);
            weights_va2 = AE_LA128_PP(p_weights);
            AE_LASX2X2_IP(wp0, wp1, weights_va2, p_weights);
            AE_LASX2X2_IP(wp2, wp3, weights_va2, p_weights);
            AE_LASX2X2_IP(s0, s1, magspec_va2, p_magspec);
            AE_LASX2X2_IP(s2, s3, magspec_va2, p_magspec);
            XT_MOVF_SX2(s0, 0.f, b0); XT_MOVF_SX2(s1, 0.f, b1);
            XT_MOVF_SX2(s2, 0.f, b2); XT_MOVF_SX2(s3, 0.f, b3);
            MADD_SX2X2(posAcc0, posAcc1, wp0, wp1, s0, s1);
            MADD_SX2X2(posAcc2, posAcc3, wp2, wp3, s2, s3);
            segmentLen -= headLen;
            *pp_magspec = (float32_t*)XT_ADDX4(headLen, (uintptr_t)*pp_magspec);
            *pp_weights = (float32_t*)XT_ADDX4(headLen, (uintptr_t)*pp_weights);
            magspec_va2 = AE_LA128_PP(*pp_magspec);
            weights_va2 = AE_LA128_PP(*pp_weights);
            /* 4 cycles per pipeline stage in steady state with unroll=1 */
            for ( n=0; n<(segmentLen>>3); n++ ) {
                AE_LASX2X2_IP(wp0, wp1, weights_va2, castxcc(xtfloatx4, *pp_weights));
                AE_LASX2X2_IP(wp2, wp3, weights_va2, castxcc(xtfloatx4, *pp_weights));
                AE_LASX2X2_IP(s0, s1, magspec_va2, castxcc(xtfloatx4, *pp_magspec));
                AE_LASX2X2_IP(s2, s3, magspec_va2, castxcc(xtfloatx4, *pp_magspec));
                MADD_SX2X2(posAcc0, posAcc1, wp0, wp1, s0, s1);
                MADD_SX2X2(posAcc2, posAcc3, wp2, wp3, s2, s3);
            } /* n */
            ADD_SX2X2(posAcc0, posAcc1, posAcc0, posAcc1, posAcc2, posAcc3);
            posAcc0 = XT_ADD_SX2(posAcc0, posAcc1);
            posAcc = XT_RADD_SX2(posAcc0);
        } else { /* (segmentLen<8) */
            posAcc = *p_posAcc;
            /* 4 cycles per pipeline stage in steady state with unroll=1 */
            for ( n=0; n<segmentLen; n++ ) {
                xtfloat wp, s;
                AE_LSIP(wp, castxcc(xtfloat, *pp_weights), sz_f32);
                AE_LSIP(s, castxcc(xtfloat, *pp_magspec), sz_f32);
                XT_MADD_S(posAcc, wp, s);
            } /* n */
        } /* (segmentLen<8) */
    } /* (imm_skip_neg_acc) */
    if (!imm_skip_neg_acc) {
        *p_negAcc = negAcc;
    }
    *p_posAcc = posAcc;
#else /* HAVE_FPU */
    float32_t wp, wn, s; /* Positive-slope weight; negative-slope weight; spectrum sample */
    float32_t negAcc, posAcc;
    int n;
    NASSERT(imm_skip_neg_acc || p_negAcc);
    if (!imm_skip_neg_acc) {
        negAcc = *p_negAcc; posAcc = *p_posAcc;
        /* 5 cycles per pipeline stage in steady state with unroll=1 */
        for ( n=0; n<segmentLen; n++ ) {
            XT_LSIP(s, *pp_magspec, sz_f32);
            XT_LSIP(wp, *pp_weights, sz_f32);
            wn = XT_SUB_S(XT_CONST_S(1), wp);
            XT_MADD_S(negAcc, wn, s);
            XT_MADD_S(posAcc, wp, s);
        } /* n */
    } else { /* (imm_skip_neg_acc) */
        posAcc = *p_posAcc;
        /* 4 cycles per pipeline stage in steady state with unroll=1 */
        for ( n=0; n<segmentLen; n++ ) {
            XT_LSIP(s, *pp_magspec, sz_f32);
            XT_LSIP(wp, *pp_weights, sz_f32);
            XT_MADD_S(posAcc, wp, s);
        } /* n */
    } /* (imm_skip_neg_acc) */
    if (!imm_skip_neg_acc) {
        *p_negAcc = negAcc;
    }
    *p_posAcc = posAcc;
#endif /* HAVE_FPU */
} /* calcFbe_Kernel() */

void calcFbe_Scale( float32_t * restrict fbe,       /* In/Out */
              const float32_t * restrict fbeScales, /* In, optional */
              int bandNum, int specExp )
{
    /*
     * % FBE normalization
     * fbe = pow2(fbe,scaleExp).*logMel.fbeScales(:);
     * % The lower limit is necessary to avoid huge negatives at spectral zeros.
     * fbe = max(1,fbe);
     */
#if USE_REFERENCE_CODE
    int m;
    NASSERT_ALIGN(fbe, ALIGN_SIZE);
    NASSERT_ALIGN(fbeScales, ALIGN_SIZE);
    if (NULL!=fbeScales) {
        for ( m=0; m<bandNum; m++ ) {
            fbe[m] = MAX(1.f, ldexpf(fbe[m]*fbeScales[m], specExp));
        }
    } else {
        for ( m=0; m<bandNum; m++ ) {
            fbe[m] = MAX(1.f, ldexpf(fbe[m], specExp));
        }
    }
#elif HAVE_VFPU
    const xtfloatx4 * restrict p_fbe_r  = (xtfloatx4*)fbe;
          xtfloatx4 * restrict p_fbe_w  = (xtfloatx4*)fbe;
    const xtfloatx4 * restrict p_scales = (xtfloatx4*)fbeScales;
    xtfloat specScale = XT_WFR((uint32_t)(specExp + 127) << 23);
    int m;
    NASSERT_ALIGN(fbe, ALIGN_SIZE);
    NASSERT_ALIGN(fbeScales, ALIGN_SIZE);
    if (NULL!=fbeScales) {
        __Pragma("loop_count factor=4");
        /* 8 cycles per pipeline stage in steady state with unroll=4 */
        for ( m=0; m<((bandNum>>4)<<2); m++ ) {
            xtfloatx2 f0, f1, s0, s1;
            AE_LSX2X2_IP(f0, f1, p_fbe_r, 4*sz_f32);
            AE_LSX2X2_IP(s0, s1, p_scales, 4*sz_f32);
            MUL_SX2X2(f0, f1, f0, f1, s0, s1);
            MULQ_S(f0, f1, f0, f1, specScale);
            f0 = XT_MAX_SX2(XT_CONST_S(1), f0);
            f1 = XT_MAX_SX2(XT_CONST_S(1), f1);
            AE_SSX2X2_IP(f0, f1, p_fbe_w, 4*sz_f32);
        } /* m */
        __Pragma("no_unroll");
        /* 5 cycles per pipeline stage in steady state with unroll=1 */
        for ( m=0; m<((bandNum>>1)&7); m++ ) {
            xtfloatx2 f, s;
            XT_LSX2IP(f, castxcc(xtfloatx2, p_fbe_r), 2*sz_f32);
            XT_LSX2IP(s, castxcc(xtfloatx2, p_scales), 2*sz_f32);
            f = XT_MUL_SX2(specScale, XT_MUL_SX2(f, s));
            f = XT_MAX_SX2(XT_CONST_S(1), f);
            XT_SSX2IP(f, castxcc(xtfloatx2, p_fbe_w), 2*sz_f32);
        } /* m */
        if (bandNum&1) {
            xtfloat f, s;
            f = XT_LSI((xtfloat*)p_fbe_r, 0);
            s = XT_LSI((xtfloat*)p_scales, 0);
            f = XT_MUL_S(specScale, XT_MUL_S(f, s));
            f = XT_MAX_S(XT_CONST_S(1), f);
            XT_SSI(f, (xtfloat*)p_fbe_w, 0);
        } /* (bandNum&1) */
    } else { /* (NULL==fbeScales) */
        __Pragma("loop_count factor=4");
        /* 6 cycles per pipeline stage in steady state with unroll=4 */
        for ( m=0; m<((bandNum>>4)<<2); m++ ) {
            xtfloatx2 f0, f1;
            AE_LSX2X2_IP(f0, f1, p_fbe_r, 4*sz_f32);
            MULQ_S(f0, f1, f0, f1, specScale);
            f0 = XT_MAX_SX2(XT_CONST_S(1), f0);
            f1 = XT_MAX_SX2(XT_CONST_S(1), f1);
            AE_SSX2X2_IP(f0, f1, p_fbe_w, 4*sz_f32);
        } /* m */
        __Pragma("no_unroll");
        /* 4 cycles per pipeline stage in steady state with unroll=1 */
        for ( m=0; m<((bandNum>>1)&7); m++ ) {
            xtfloatx2 f;
            XT_LSX2IP(f, castxcc(xtfloatx2, p_fbe_r), 2*sz_f32);
            f = XT_MUL_SX2(specScale, f);
            f = XT_MAX_SX2(XT_CONST_S(1), f);
            XT_SSX2IP(f, castxcc(xtfloatx2, p_fbe_w), 2*sz_f32);
        } /* m */
        if (bandNum&1) {
            xtfloat f;
            f = XT_LSI((xtfloat*)p_fbe_r, 0);
            f = XT_MUL_S(specScale, f);
            f = XT_MAX_S(XT_CONST_S(1), f);
            XT_SSI(f, (xtfloat*)p_fbe_w, 0);
        } /* (bandNum&1) */
    } /* (NULL==fbeScales) */
#else /* HAVE_FPU */
    const xtfloat * restrict p_fbe_r  = (xtfloat*)fbe;
          xtfloat * restrict p_fbe_w  = (xtfloat*)fbe;
    const xtfloat * restrict p_scales = (xtfloat*)fbeScales;
    xtfloat specScale = XT_WFR((uint32_t)(specExp + 127) << 23);
    xtfloat c1f = XT_CONST_S(1);
    int m;
    NASSERT_ALIGN(fbe, ALIGN_SIZE);
    NASSERT_ALIGN(fbeScales, ALIGN_SIZE);
    if (NULL!=fbeScales) {
        /* 7 cycles per pipeline stage in steady state with unroll=1 */
        for ( m=0; m<bandNum; m++ ) {
            xtfloat f, s;
            XT_LSIP(f, p_fbe_r, sz_f32);
            XT_LSIP(s, p_scales, sz_f32);
            f = XT_MUL_S(specScale, XT_MUL_S(f, s));
            XT_MOVT_S(f, c1f, XT_OLT_S(f, c1f));
            XT_SSIP(f, p_fbe_w, sz_f32);
        }
    } else {
        /* 20 cycles per pipeline stage in steady state with unroll=4 */
        for ( m=0; m<bandNum; m++ ) {
            xtfloat f;
            XT_LSIP(f, p_fbe_r, sz_f32);
            f = XT_MUL_S(specScale, f);
            XT_MOVT_S(f, c1f, XT_OLT_S(f, c1f));
            XT_SSIP(f, p_fbe_w, sz_f32);
        }
    }
#endif /* HAVE_FPU */
} /* calcFbe_Scale() */

#endif /* HAVE_FPU || HAVE_VFPU */
