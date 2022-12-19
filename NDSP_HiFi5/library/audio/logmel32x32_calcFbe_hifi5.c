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

#define USE_REFERENCE_CODE  0
#define ALIGN_SIZE          (HIFI_SIMD_WIDTH)
#define sz_i16              sizeof(int16_t)
#define sz_i32              sizeof(int32_t)

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
#include "baseop.h"

/* Compute the Normalization Shift Amount for a 64-bit signed input argument. 
 * Return 63 if the input argument is zero. */
static int16_t nsa64(int64_t x)
{
    int16_t z = 0;
    if (x == 0)  return 0x3F;
    while ((int64_t)(x ^ (x << 1))>0) //MSB != MSB-1
    {
        x <<= 1;
        z++;
    }
    return z;
} /* nsa64() */

/* Q(x+y-15) <- Qx*Qy - 15 w/ asym. rounding */
static int32_t mulf32x16ras(int32_t x, int16_t y)
{
    int64_t z;
    z=(int64_t)x*y+(1L<<14);
    z>>=15;
    if (z>MAX_INT32)z=MAX_INT32;
    if (z<MIN_INT32)z=MIN_INT32;
    return (int32_t)z;
} /* mulf32x16ras() */

/* Q(x+y-31) <- Qx*Qy - 31 w/ asym. rounding */
static int32_t mulf32x32ras(int32_t x, int32_t y)
{
    int64_t z;
    z=(int64_t)x*y+(1L<<30);
    z>>=31;
    if (z>MAX_INT32)z=MAX_INT32;
    if (z<MIN_INT32)z=MIN_INT32;
    return (int32_t)z;
} /* mulf32x32ras() */

/* simplified fractional multiply Q63*Q31->Q63 */
static int64_t mulf64x32(int64_t x,int32_t y)
{
    int64_t z;
    uint32_t xlo;
    int32_t xhi;
    xlo=(uint32_t)x;
    xhi=(int32_t)(x>>32);
    z=((int64_t)xlo*y)>>32;
    z+=(int64_t)xhi*y;
    z<<=1;
    return z;
} /* mulf64x32() */
#endif /* USE_REFERENCE_CODE */

/* Immediate options for inlined kernel routine */
#define CFK_SKIP_NEG_NO    0  /* Do not update the negative-slope accumulator */
#define CFK_SKIP_NEG_YES   1  /* Update the negative-slope accumulator */
#define CFK_SCALE_NEG_NO   0  /* Do not scale the negative-slope result */
#define CFK_SCALE_NEG_YES  1  /* Normalize the negative-slope result */
#define CFK_SCALE_POS_NO   0  /* Do not scale the positive-slope result */
#define CFK_SCALE_POS_YES  1  /* Normalize the positive-slope result */

static int ATTRIBUTE_ALWAYS_INLINE calcFbe_Kernel( int64_t *  p_negAcc,
                                                   int64_t *  p_posAcc,
                                             const fract32 ** pp_magspec,
                                             const fract32 ** pp_weights,
                                             const fract32 ** pp_fbeScale_fract,
                                             const int16_t ** pp_fbeScale_exp,
                                             int segmentLen, int imm_skip_neg_acc, 
                                             int imm_scale_neg_acc, int imm_scale_pos_acc );

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

/* 32-bit fixed-point variant */
void logmel32x32_calcFbe( fract32 * restrict fbe,
                          int16_t * restrict fbeExp,
                    const fract32 * restrict magspec,
                    const fract32 * restrict weights,
                    const int16_t * restrict segments,
                    const fract32 * restrict fbeScales_fract,
                    const int16_t * restrict fbeScales_exp,
                    int binNum, int bandNum, int specExp )
{
#if !USE_REFERENCE_CODE
    const fract32 * p_magspec        = magspec;
    const fract32 * p_weights        = weights;
    const fract32 * p_fbeScale_fract = fbeScales_fract;
    const int16_t * p_fbeScale_exp   = fbeScales_exp;
    ae_int64 negAcc, posAcc; 
    int fbNormOpt, normExp=0, nsa; 
    int m, segmentLen;
    NASSERT(bandNum>0);
    NASSERT_ALIGN(fbe, ALIGN_SIZE);
    NASSERT_ALIGN(fbeExp, ALIGN_SIZE);
    NASSERT_ALIGN(magspec, ALIGN_SIZE);
    NASSERT_ALIGN(weights, ALIGN_SIZE);
    NASSERT_ALIGN(segments, ALIGN_SIZE);
    NASSERT_ALIGN(fbeScales_fract, ALIGN_SIZE);
    NASSERT_ALIGN(fbeScales_exp, ALIGN_SIZE);
    /* FBE normalization option. */
    fbNormOpt = (NULL!=fbeScales_fract);
    /* Positive-slope side of the left-most triangle. */
    segmentLen = *segments++; posAcc = 0;
    calcFbe_Kernel(NULL, (int64_t*)&posAcc, 
                   &p_magspec, &p_weights, NULL, NULL, segmentLen, 
                   CFK_SKIP_NEG_YES, CFK_SCALE_NEG_NO, CFK_SCALE_POS_NO);
    /* Because of 1/2 bands overlapping, each inner segment simultaneously updates
     * a pair of bands: negative-slope side for the left triangle, positive-slope 
     * side for the right triangle. */
    if (fbNormOpt) {
        for ( m=0; m<bandNum-1; m++ ) {
            segmentLen = *segments++; negAcc = posAcc; posAcc = 0;
            normExp = calcFbe_Kernel((int64_t*)&negAcc, (int64_t*)&posAcc, 
                                     &p_magspec, &p_weights, &p_fbeScale_fract, &p_fbeScale_exp, segmentLen,
                                     CFK_SKIP_NEG_NO, CFK_SCALE_NEG_YES, CFK_SCALE_POS_NO);
            nsa = AE_NSA64(negAcc);
            /* Pack the accumulated energy to a 32-bit value:
             * Q(31+(nsa-normExp-specExp-17)) <- Q(46-normExp)*2^specExp + nsa - specExp - 32 */
            AE_S32_L_IP(AE_ROUND32F64SASYM(AE_SLAA64S(negAcc, nsa)), castxcc(ae_int32, fbe), sz_i32);
            *fbeExp++ = nsa-normExp-specExp+(46-32)-31;
        } /* m */
    } else {
        for ( m=0; m<bandNum-1; m++ ) {
            segmentLen = *segments++; negAcc = posAcc; posAcc = 0;
            normExp = calcFbe_Kernel((int64_t*)&negAcc, (int64_t*)&posAcc, 
                                     &p_magspec, &p_weights, &p_fbeScale_fract, &p_fbeScale_exp, segmentLen,
                                     CFK_SKIP_NEG_NO, CFK_SCALE_NEG_NO, CFK_SCALE_POS_NO);
            nsa = AE_NSA64(negAcc);
            /* Pack the accumulated energy to a 32-bit value:
             * Q(31+(nsa-normExp-specExp-17)) <- Q(46-normExp)*2^specExp + nsa - specExp - 32 */
            AE_S32_L_IP(AE_ROUND32F64SASYM(AE_SLAA64S(negAcc, nsa)), castxcc(ae_int32, fbe), sz_i32);
            *fbeExp++ = nsa-normExp-specExp+(46-32)-31;
        } /* m */
    }
    /* Negative-slope side of the right-most triangle. */
    segmentLen = *segments++;
    normExp = calcFbe_Kernel(NULL, (int64_t*)&posAcc, 
                             &p_magspec, &p_weights, &p_fbeScale_fract, &p_fbeScale_exp, segmentLen,
                             CFK_SKIP_NEG_YES, CFK_SCALE_NEG_NO, fbNormOpt ? CFK_SCALE_POS_YES : CFK_SCALE_POS_NO);
    nsa = AE_NSA64(posAcc); 
    /* Q(31+(nsa-normExp-specExp-17)) <- Q(46-normExp)*2^specExp + nsa - specExp - 32 */
    AE_S32_L_I(AE_ROUND32F64SASYM(AE_SLAA64S(posAcc, nsa)), (ae_int32*)fbe, 0);
    *fbeExp++ = nsa-normExp-specExp+(46-32)-31;
#else /* USE_REFERENCE_CODE */
    const fract32 *ps, *pw;
    fract32 wp, wn, s; /* Positive-slope weight; negative-slope weight; sample */
    int64_t accL, accR; 
    int fbNormOpt, normExp=0, nsa; 
    int m, n;
    NASSERT(bandNum>0);
    NASSERT_ALIGN(fbe, ALIGN_SIZE);
    NASSERT_ALIGN(fbeExp, ALIGN_SIZE);
    NASSERT_ALIGN(magspec, ALIGN_SIZE);
    NASSERT_ALIGN(weights, ALIGN_SIZE);
    NASSERT_ALIGN(segments, ALIGN_SIZE);
    NASSERT_ALIGN(fbeScales_fract, ALIGN_SIZE);
    NASSERT_ALIGN(fbeScales_exp, ALIGN_SIZE);
    /* FBE normalization option. */
    fbNormOpt = (NULL!=fbeScales_fract);
    /* Positive-slope side of the left-most triangle. */
    ps = magspec; pw = weights; accL = 0;
    for ( n=0; n<segments[0]; n++ ) {
        wp = *pw++; s = *ps++;
        /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
        accL += ((int64_t)wp*s+(1<<14))>>15;
    } /* n */
    /* Because of 1/2 bands overlapping, each inner segment simultaneously updates
     * a pair of bands: negative-slope side for the left triangle, positive-slope 
     * side for the right triangle. */
    for ( m=0; m<bandNum-1; m++ ) {
        accR = 0;
        for ( n=0; n<segments[m+1]; n++ ) {
            wp = *pw++; s = *ps++;
            wn = 0x80000000-wp;
            /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
            accL += ((int64_t)wn*s+(1<<14))>>15;
            accR += ((int64_t)wp*s+(1<<14))>>15;
        } /* n */
        /* Optionally normalize the filterbank energu  */
        if (fbNormOpt) {
            /* Q(46-(31-scalesExp[m])) <- Q46*Q(scalesExp[m]) - 31 */
            accL = mulf64x32(accL, fbeScales_fract[m]);
            normExp = 31-fbeScales_exp[m];
        }
        nsa = nsa64(accL); 
        /* Pack the accumulated energy to a 32-bit value:
         * Q(31+(nsa-normExp-specExp-17)) <- Q(46-normExp)*2^specExp + nsa - specExp - 32 */
        fbe[m] = roundQ63_Q31(accL<<nsa); fbeExp[m] = nsa-normExp-specExp+(46-32)-31;
        accL = accR;
    } /* m */
    /* Negative-slope side of the right-most triangle. */
    for ( n=0; n<segments[bandNum]; n++ ) {
        wn = *pw++; s = *ps++;
        /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
        accR += ((int64_t)wn*s+(1<<14))>>15;
    } /* n */
    if (fbNormOpt) {
        /* Q(46-(31-scalesExp[m])) <- Q46*Q(scalesExp[m]) - 31 */
        accR = mulf64x32(accR, fbeScales_fract[bandNum-1]);
        normExp = 31-fbeScales_exp[bandNum-1];
    }
    nsa = nsa64(accR); 
    /* Q(31+(nsa-normExp-specExp-17)) <- Q(46-normExp)*2^specExp + nsa - specExp - 32 */
    fbe[bandNum-1] = roundQ63_Q31(accR<<nsa); fbeExp[m] =  nsa-normExp-specExp+(46-32)-31;
#endif /* USE_REFERENCE_CODE */
} /* logmel32x32_calcFbe() */

int calcFbe_Kernel( int64_t *  p_negAcc,
                    int64_t *  p_posAcc,
              const fract32 ** pp_magspec,
              const fract32 ** pp_weights,
              const fract32 ** pp_fbeScale_fract,
              const int16_t ** pp_fbeScale_exp,
              int segmentLen, int imm_skip_neg_acc, 
              int imm_scale_neg_acc, int imm_scale_pos_acc )
{
#if 0
    const fract32 * p_magspec = *pp_magspec;
    const fract32 * p_weights = *pp_weights;
    int64_t negAcc, posAcc;
    fract32 fbeScale_fract = 0;
    int16_t fbeScale_exp = 0;
    int n, normExp=0;
    NASSERT(imm_skip_neg_acc || p_negAcc);
    NASSERT((!imm_scale_neg_acc && !imm_scale_pos_acc) || (pp_fbeScale_fract && pp_fbeScale_exp));
    if (!imm_skip_neg_acc) {
        negAcc = *p_negAcc;
    }
    posAcc = *p_posAcc;
    for ( n=0; n<segmentLen; n++ ) {
        fract32 wp, wn, s; /* Positive-slope weight; negative-slope weight; spectrum sample */
        s = *p_magspec++; /* Q30 */
        wp = *p_weights++; /* Q31 */
        if (!imm_skip_neg_acc) {
            /* wn <- 1-wp; Q31 */
            wn = (1UL<<31)-wp;
            /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
            negAcc += ((int64_t)wn*s+(1<<14))>>15;
        }
        /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
        posAcc += ((int64_t)wp*s+(1<<14))>>15;
    } /* n */
    if (imm_scale_neg_acc || imm_scale_pos_acc) {
        fbeScale_fract = *(*pp_fbeScale_fract)++;
        fbeScale_exp = *(*pp_fbeScale_exp)++;
        normExp = 31-fbeScale_exp;
    }
    if (!imm_skip_neg_acc && imm_scale_neg_acc) {
        /* Q(46-(31-scalesExp[m])) <- Q46*Q(scalesExp[m]) - 31 */
        negAcc = mulf64x32(negAcc, fbeScale_fract);
    }
    if (imm_scale_pos_acc) {
        /* Q(46-(31-scalesExp[m])) <- Q46*Q(scalesExp[m]) - 31 */
        posAcc = mulf64x32(posAcc, fbeScale_fract);
    }
    if (!imm_skip_neg_acc) {
        *p_negAcc = negAcc;
    }
    *p_posAcc = posAcc;
    *pp_magspec = p_magspec;
    *pp_weights = p_weights;
    return normExp;
#else
    const ae_int32x4 * p_magspec;
    const ae_int32x4 * p_weights;
    ae_valignx2 magspec_va2, weights_va2;
    /* Positive-slope weights; negative-slope weights; spectrum samples. */
    ae_int32x2 wp0, wp1, wn0, wn1, s0, s1;
    ae_int64 negAcc0, negAcc1, posAcc0, posAcc1;
    ae_int64 wlo, whi;
    ae_ep ep;
    ae_int32x2 fbeScale_fract;
    int n, headLen, normExp=0;
    NASSERT(imm_skip_neg_acc || p_negAcc);
    NASSERT((!imm_scale_neg_acc && !imm_scale_pos_acc) || (pp_fbeScale_fract && pp_fbeScale_exp));
    if (!imm_skip_neg_acc) {
        negAcc0 = *(ae_int64*)p_negAcc; negAcc1 = 0;
    }
    posAcc0 = *(ae_int64*)p_posAcc; posAcc1 = 0;
    if (!imm_skip_neg_acc) {
        if (segmentLen<8) {
            __Pragma("no_unroll");
            /* 2 cycles per pipeline stage in steady state with unroll=1 */
            for ( n=0; n<segmentLen; n++ ) {
                AE_L32_IP(s0, castxcc(ae_int32, *pp_magspec), sz_i32); /* Q30 */
                AE_L32_IP(wp0, castxcc(ae_int32, *pp_weights), sz_i32); /* Q31 */
                /* wn <- 1-wp; Q31 */
                wn0 = AE_SUB32(MIN_INT32, wp0);
                /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
                AE_MULAF32RA_HH(negAcc0, s0, wn0);
                AE_MULAF32RA_HH(posAcc0, s0, wp0);
            }
        } else { /* (segmentLen<8) */
            static const int32_t ALIGN(16) seq[] = {0,1,2,3,4,5,6,7};
            xtbool2 b01, b23, b45, b67;
            ae_int32x2 s01, s23, s45, s67;
            AE_L32X2X2_I(s01, s23, (const ae_int32x4*)seq, 0*sizeof(ae_int32x4));
            AE_L32X2X2_I(s45, s67, (const ae_int32x4*)seq, 1*sizeof(ae_int32x4));
            headLen = ((segmentLen-1)&7)+1;
            b01 = AE_LT32(s01, headLen); /* Mask useful elements for the prologue. */
            b23 = AE_LT32(s23, headLen);
            b45 = AE_LT32(s45, headLen);
            b67 = AE_LT32(s67, headLen);
            p_magspec = (ae_int32x4*)*pp_magspec;
            p_weights = (ae_int32x4*)*pp_weights;
            magspec_va2 = AE_LA128_PP(p_magspec);
            weights_va2 = AE_LA128_PP(p_weights);
            AE_LA32X2X2_IP(s0, s1, magspec_va2, p_magspec); /* Q30 */
            AE_LA32X2X2_IP(wp0, wp1, weights_va2, p_weights); /* Q31 */
            AE_MOVF32X2(s0, 0, b01);
            AE_MOVF32X2(s1, 0, b23);
            /* wn <- 1-wp; Q31 */
            wn0 = AE_SUB32(1U<<31, wp0);
            wn1 = AE_SUB32(1U<<31, wp1);
            /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
            AE_MULAAF2D32RA_HH_LL(negAcc0, posAcc0, s0, s0, wn0, wp0);
            AE_MULAAF2D32RA_HH_LL(negAcc1, posAcc1, s1, s1, wn1, wp1);
            AE_LA32X2X2_IP(s0, s1, magspec_va2, p_magspec); /* Q30 */
            AE_LA32X2X2_IP(wp0, wp1, weights_va2, p_weights); /* Q31 */
            AE_MOVF32X2(s0, 0, b45);
            AE_MOVF32X2(s1, 0, b67);
            /* wn <- 1-wp; Q31 */
            wn0 = AE_SUB32(1U<<31, wp0);
            wn1 = AE_SUB32(1U<<31, wp1);
            /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
            AE_MULAAF2D32RA_HH_LL(negAcc0, posAcc0, s0, s0, wn0, wp0);
            AE_MULAAF2D32RA_HH_LL(negAcc1, posAcc1, s1, s1, wn1, wp1);
            segmentLen -= headLen;
            *pp_magspec = (fract32*)XT_ADDX4(headLen, (uintptr_t)*pp_magspec);
            *pp_weights = (fract32*)XT_ADDX4(headLen, (uintptr_t)*pp_weights);
            magspec_va2 = AE_LA128_PP(*pp_magspec);
            weights_va2 = AE_LA128_PP(*pp_weights);
            __Pragma("no_unroll");
            /* 4 cycles per pipeline stage in steady state with unroll=1 */
            for ( n=0; n<(segmentLen>>3); n++ ) {
                AE_LA32X2X2_IP(s0, s1, magspec_va2, castxcc(ae_int32x4, *pp_magspec)); /* Q30 */
                AE_LA32X2X2_IP(wp0, wp1, weights_va2, castxcc(ae_int32x4, *pp_weights)); /* Q31 */
                /* wn <- 1-wp; Q31 */
                wn0 = AE_SUB32(1U<<31, wp0);
                wn1 = AE_SUB32(1U<<31, wp1);
                /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
                AE_MULAAF2D32RA_HH_LL(negAcc0, posAcc0, s0, s0, wn0, wp0);
                AE_MULAAF2D32RA_HH_LL(negAcc1, posAcc1, s1, s1, wn1, wp1);

                AE_LA32X2X2_IP(s0, s1, magspec_va2, castxcc(ae_int32x4, *pp_magspec)); /* Q30 */
                AE_LA32X2X2_IP(wp0, wp1, weights_va2, castxcc(ae_int32x4, *pp_weights)); /* Q31 */
                /* wn <- 1-wp; Q31 */
                wn0 = AE_SUB32(1U<<31, wp0);
                wn1 = AE_SUB32(1U<<31, wp1);
                /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
                AE_MULAAF2D32RA_HH_LL(negAcc0, posAcc0, s0, s0, wn0, wp0);
                AE_MULAAF2D32RA_HH_LL(negAcc1, posAcc1, s1, s1, wn1, wp1);
            } /* n */
            negAcc0 = AE_ADD64S(negAcc0, negAcc1);
            posAcc0 = AE_ADD64S(posAcc0, posAcc1);
        } /* (segmentLen>=8) */
    } else { /* (imm_skip_neg_acc) */
        if (segmentLen<8) {
            __Pragma("no_unroll");
            /* 1 cycles per pipeline stage in steady state with unroll=1 */
            for ( n=0; n<segmentLen; n++ ) {
                AE_L32_IP(s0, castxcc(ae_int32, *pp_magspec), sz_i32); /* Q30 */
                AE_L32_IP(wp0, castxcc(ae_int32, *pp_weights), sz_i32); /* Q31 */
                /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
                AE_MULAF32RA_HH(posAcc0, s0, wp0);
            }
        } else { /* (segmentLen<8) */
            static const int32_t ALIGN(16) seq[] = {0,1,2,3,4,5,6,7};
            xtbool2 b01, b23, b45, b67;
            ae_int32x2 s01, s23, s45, s67;
            AE_L32X2X2_I(s01, s23, (const ae_int32x4*)seq, 0*sizeof(ae_int32x4));
            AE_L32X2X2_I(s45, s67, (const ae_int32x4*)seq, 1*sizeof(ae_int32x4));
            headLen = ((segmentLen-1)&7)+1;
            b01 = AE_LT32(s01, headLen); /* Mask useful elements for the prologue. */
            b23 = AE_LT32(s23, headLen);
            b45 = AE_LT32(s45, headLen);
            b67 = AE_LT32(s67, headLen);
            p_magspec = (ae_int32x4*)*pp_magspec;
            p_weights = (ae_int32x4*)*pp_weights;
            magspec_va2 = AE_LA128_PP(p_magspec);
            weights_va2 = AE_LA128_PP(p_weights);
            AE_LA32X2X2_IP(s0, s1, magspec_va2, p_magspec); /* Q30 */
            AE_LA32X2X2_IP(wp0, wp1, weights_va2, p_weights); /* Q31 */
            AE_MOVF32X2(s0, 0, b01);
            AE_MOVF32X2(s1, 0, b23);
            /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
            AE_MULAAF2D32RA_HH_LL(posAcc0, posAcc1, s0, s1, wp0, wp1);
            AE_LA32X2X2_IP(s0, s1, magspec_va2, p_magspec); /* Q30 */
            AE_LA32X2X2_IP(wp0, wp1, weights_va2, p_weights); /* Q31 */
            AE_MOVF32X2(s0, 0, b45);
            AE_MOVF32X2(s1, 0, b67);
            /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
            AE_MULAAF2D32RA_HH_LL(posAcc0, posAcc1, s0, s1, wp0, wp1);
            segmentLen -= headLen;
            *pp_magspec = (fract32*)XT_ADDX4(headLen, (uintptr_t)*pp_magspec);
            *pp_weights = (fract32*)XT_ADDX4(headLen, (uintptr_t)*pp_weights);
            magspec_va2 = AE_LA128_PP(*pp_magspec);
            weights_va2 = AE_LA128_PP(*pp_weights);
            __Pragma("no_unroll");
            /* 2 cycles per pipeline stage in steady state with unroll=1 */
            for ( n=0; n<(segmentLen>>3); n++ ) {
                AE_LA32X2X2_IP(s0, s1, magspec_va2, castxcc(ae_int32x4, *pp_magspec)); /* Q30 */
                AE_LA32X2X2_IP(wp0, wp1, weights_va2, castxcc(ae_int32x4, *pp_weights)); /* Q31 */
                /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
                AE_MULAAF2D32RA_HH_LL(posAcc0, posAcc1, s0, s1, wp0, wp1);

                AE_LA32X2X2_IP(s0, s1, magspec_va2, castxcc(ae_int32x4, *pp_magspec)); /* Q30 */
                AE_LA32X2X2_IP(wp0, wp1, weights_va2, castxcc(ae_int32x4, *pp_weights)); /* Q31 */
                /* Q46 <- Q31*Q30 - 15 w/ asym. rounding */
                AE_MULAAF2D32RA_HH_LL(posAcc0, posAcc1, s0, s1, wp0, wp1);
            } /* n */
            posAcc0 = AE_ADD64S(posAcc0, posAcc1);
        } /* (segmentLen>=8) */
    } /* (imm_skip_neg_acc) */
    if (imm_scale_neg_acc || imm_scale_pos_acc) {
        /* fbeScale_fract[m]: Q(scalesExp[m]) */
        AE_L32_IP(fbeScale_fract, castxcc(ae_int32, *pp_fbeScale_fract), sz_i32);
        normExp = 31 - *(*pp_fbeScale_exp)++;
    }
    if (!imm_skip_neg_acc && imm_scale_neg_acc) {
        /* Q(46-(31-scalesExp[m])) <- Q46*Q(scalesExp[m]) - 31 */
        whi = AE_MULF32S_HH(AE_MOVINT32X2_FROMF64(negAcc0), fbeScale_fract);
        AE_MUL32USEP_LH(ep, wlo, AE_MOVINT32X2_FROMF64(negAcc0), fbeScale_fract);
        wlo = AE_SRAI72(ep, wlo, 31);
        negAcc0 = AE_ADD64S(wlo, whi);
    }
    if (imm_scale_pos_acc) {
        /* Q(46-(31-scalesExp[m])) <- Q46*Q(scalesExp[m]) - 31 */
        whi = AE_MULF32S_HH(AE_MOVINT32X2_FROMF64(posAcc0), fbeScale_fract);
        AE_MUL32USEP_LH(ep, wlo, AE_MOVINT32X2_FROMF64(posAcc0), fbeScale_fract);
        wlo = AE_SRAI72(ep, wlo, 31);
        posAcc0 = AE_ADD64S(wlo, whi);
    }
    if (!imm_skip_neg_acc) {
        *(ae_int64*)p_negAcc = negAcc0;
    }
    *(ae_int64*)p_posAcc = posAcc0;
    return normExp;
#endif
} /* calcFbe_Kernel() */
