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
    Compute Mel-Frequency Cepstrum Coefficients (MFCC) from speech signal
    32-bit fixed-point variant
    C code optimized for HiFi5 core
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
/* MFCC extractor internal declarations. */
#include "mfcc_internal.h"

#define USE_REFERENCE_CODE    0
#define ALIGN_SIZE            (HIFI_SIMD_WIDTH)

#if USE_REFERENCE_CODE
/* NON OPTIMIZED REFERENCE CODE: to be use it for educational purposes only */
#include "baseop.h"

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
#endif /* USE_REFERENCE_CODE */

/*
 * Perform pre-emphasis filtering (1st order FIR), similarly to the 
 * following MATLAB code:
 *   y = filter([1,-alpha],1,[st;x]); st = x(end);
 * Input and output arguments x[N] and y[N] may refer to the same array.
 * Input:
 *   alpha      Pre-emphasis coefficient; Q15 for 32x32
 *   st         Initial filter state
 *   N          Number of samples to be processed
 *   x[N]       Input samples
 * Output:
 *   y[N]       Output samples
 *   return     Updated filter state
 * Restrictions:
 *   N          Must be a multiple of 2
 *   x[N],y[N]  Must be 16-bytes aligned
 */

fract32 mfcc32x32_preemph( fract32 * y, const fract32 * x, fract16 alpha, fract32 st, int N )
{
#if !USE_REFERENCE_CODE
    const ae_int32x4 * restrict pX = (ae_int32x4*)&x[0];
    ae_int32x4 * restrict pY = (ae_int32x4*)&y[0];
    ae_int32x2 f, h;
    ae_int32x2 x0, x1, x2, x3, x4, x5, x6, x7;
    ae_int32x2 y0, y1, y2, y3;
    ae_int64 w0, w1, w2, w3, w4, w5, w6, w7;
    int n;
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(0==(N%2));
    f = AE_MOVDA32(st);
    /* Q31 = Q15 + 16 */
    h = AE_NEG32S(AE_MOVDA32X2((int32_t)alpha<<16, MIN_INT32));
    /* 7 cycles per pipeline stage in steady state with unroll=1 */
    for ( n=0; n<(N>>4); n++ ) {
        AE_L32X2X2_IP(x0, x1, pX, 4*sizeof(ae_int32x4)); /* Q31 */
        AE_L32X2X2_I (x2, x3, pX, -3*(int)sizeof(ae_int32x4));
        AE_L32X2X2_I (x4, x5, pX, -2*(int)sizeof(ae_int32x4));
        AE_L32X2X2_I (x6, x7, pX, -1*(int)sizeof(ae_int32x4));
        /* Q63 = Q31*Q31 + 1 */
        AE_MULFD32X2S_FIR_L(w0, w1, f, x0, h);
        AE_MULFD32X2S_FIR_L(w2, w3, x0, x1, h);
        AE_MULFD32X2S_FIR_L(w4, w5, x1, x2, h);
        AE_MULFD32X2S_FIR_L(w6, w7, x2, x3, h);
        /* Q31 = Q63 - 32 w/ asym. rounding */
        y0 = AE_ROUND32X2F64SASYM(w0, w1);
        y1 = AE_ROUND32X2F64SASYM(w2, w3);
        y2 = AE_ROUND32X2F64SASYM(w4, w5);
        y3 = AE_ROUND32X2F64SASYM(w6, w7);
        AE_S32X2X2_IP(y0, y1, pY, sizeof(ae_int32x4)); /* Q31 */
        AE_S32X2X2_IP(y2, y3, pY, sizeof(ae_int32x4));
        /* Q63 = Q31*Q31 + 1 */
        AE_MULFD32X2S_FIR_L(w0, w1, x3, x4, h);
        AE_MULFD32X2S_FIR_L(w2, w3, x4, x5, h);
        AE_MULFD32X2S_FIR_L(w4, w5, x5, x6, h);
        AE_MULFD32X2S_FIR_L(w6, w7, x6, x7, h); f = x7;
        /* Q31 = Q63 - 32 w/ asym. rounding */
        y0 = AE_ROUND32X2F64SASYM(w0, w1);
        y1 = AE_ROUND32X2F64SASYM(w2, w3);
        y2 = AE_ROUND32X2F64SASYM(w4, w5);
        y3 = AE_ROUND32X2F64SASYM(w6, w7);
        AE_S32X2X2_IP(y0, y1, pY, sizeof(ae_int32x4));
        AE_S32X2X2_IP(y2, y3, pY, sizeof(ae_int32x4));
    } /* n */
    /* 2 cycles per pipeline stage in steady state with unroll=2 */
    for ( n=0; n<((N>>1)&7); n++ ) {
        AE_L32X2_IP(x0, castxcc(ae_int32x2, pX), sizeof(ae_int32x2)); /* Q31 */
        /* Q63 = Q31*Q31 + 1 */
        AE_MULFD32X2S_FIR_L(w0, w1, f, x0, h); f = x0;
        /* Q31 = Q63 - 32 w/ asym. rounding */
        y0 = AE_ROUND32X2F64SASYM(w0, w1);
        AE_S32X2_IP(y0, castxcc(ae_int32x2, pY), sizeof(ae_int32x2)); /* Q31 */
    } /* n */
    return AE_MOVAD32_L(f);
#else
    int n;
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(0==(N%2));
    for ( n=0; n<N; n++ ) {
        /* Note that x[] and y[] may refer to the same array! */
        fract32 s = x[n];
        /* Qx <- Qx*Q15 - 15 w/ rounding */
        y[n] = L_sub_ll(s, mulf32x16ras(st, alpha)); 
        st = s;
    }
    return st;
#endif /* USE_REFERENCE_CODE */
} /* mfcc32x32_preemph() */
