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
    Single precision floating-point variant
    C code optimized for HiFi5 core
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common_fpu.h"
/* MFCC extractor internal declarations. */
#include "mfcc_internal.h"

#if HAVE_VFPU || HAVE_FPU

#define USE_REFERENCE_CODE    0
#define ALIGN_SIZE            (HIFI_SIMD_WIDTH)
#define sz_f32                sizeof(float32_t)

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

float32_t mfccf_preemph( float32_t * y, const float32_t * x, float32_t alpha, float32_t st, int N )
{
#if !USE_REFERENCE_CODE
#if HAVE_VFPU
    const xtfloatx4 * restrict pX0 = (xtfloatx4*)&x[0];
    const xtfloatx4 * restrict pX1 = (xtfloatx4*)&x[1];
          xtfloatx4 * restrict pY = (xtfloatx4*)&y[0];
    ae_valign va;
    ae_valignx2 va_x2;
    xtfloatx2 h = st;
    int n;
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(0==(N%2));
    if ((N>>4)>0) {
        xtfloatx2 f0, f1, f2, f3, f4, f5, f6, f7;
        xtfloatx2 g0, g1, g2, g3, g4, g5, g6, g7;
        g0 = XT_SEL32_LH_SX2(h, XT_LSI((xtfloat*)pX1, -(int)sz_f32));
        va_x2 = AE_LA128_PP(pX1);
        __Pragma("ymemory(pX0)");
        /* 6 cycles per pipeline stage in steady state with unroll=1 */
        for ( n=0; n<((N>>4)-1); n++ ) {
            AE_LSX2X2_I (f2, f3, pX0, 1*4*sz_f32);
            AE_LSX2X2_I (f4, f5, pX0, 2*4*sz_f32);
            AE_LSX2X2_I (f6, f7, pX0, 3*4*sz_f32);
            AE_LSX2X2_IP(f0, f1, pX0, 4*4*sz_f32);
            AE_LASX2X2_IP(g1, g2, va_x2, pX1);
            AE_LASX2X2_IP(g3, g4, va_x2, pX1);
            AE_LASX2X2_IP(g5, g6, va_x2, pX1);
            AE_LASX2X2_IP(g7, h , va_x2, pX1);
            MSUBQ_S(f0, f1, g0, g1, alpha); g0 = h;
            MSUBQ_S(f2, f3, g2, g3, alpha);
            MSUBQ_S(f4, f5, g4, g5, alpha);
            MSUBQ_S(f6, f7, g6, g7, alpha);
            AE_SSX2X2_IP(f0, f1, pY, 4*sz_f32);
            AE_SSX2X2_IP(f2, f3, pY, 4*sz_f32);
            AE_SSX2X2_IP(f4, f5, pY, 4*sz_f32);
            AE_SSX2X2_IP(f6, f7, pY, 4*sz_f32);
        } /* n */
        AE_LSX2X2_I (f2, f3, pX0, 1*4*sz_f32);
        AE_LSX2X2_I (f4, f5, pX0, 2*4*sz_f32);
        AE_LSX2X2_I (f6, f7, pX0, 3*4*sz_f32);
        AE_LSX2X2_IP(f0, f1, pX0, 4*4*sz_f32);
        AE_LASX2X2_IP(g1, g2, va_x2, pX1);
        AE_LASX2X2_IP(g3, g4, va_x2, pX1);
        AE_LASX2X2_IP(g5, g6, va_x2, pX1);
        va = AE_LA64_PP(pX1);
        AE_LASX2IP(g7, va, castxcc(xtfloatx2, pX1));
        MSUBQ_S(f0, f1, g0, g1, alpha); h = f7;
        MSUBQ_S(f2, f3, g2, g3, alpha);
        MSUBQ_S(f4, f5, g4, g5, alpha);
        MSUBQ_S(f6, f7, g6, g7, alpha);
        AE_SSX2X2_IP(f0, f1, pY, 4*sz_f32);
        AE_SSX2X2_IP(f2, f3, pY, 4*sz_f32);
        AE_SSX2X2_IP(f4, f5, pY, 4*sz_f32);
        AE_SSX2X2_IP(f6, f7, pY, 4*sz_f32);
    } /* ((N>>4)>0) */
    __Pragma("no_unroll");
    /* 6 cycles per pipeline stage in steady state with unroll=1 */
    for ( n=0; n<((N>>1)&7); n++ ) {
        xtfloatx2 f, g;
        XT_LSX2IP(f, castxcc(xtfloatx2, pX0), 2*sz_f32);
        g = XT_SEL32_LH_SX2(h, f); h = f;
        XT_MSUBN_SX2(f, g, alpha);
        XT_SSX2IP(f, castxcc(xtfloatx2, pY), 2*sz_f32);
    } /* n */
    return XT_LOW_S(h);
#elif HAVE_FPU
    const xtfloat * restrict X;
    xtfloat * restrict Y;
    float32_t f, g, h;
    int n;
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(0==(N%2));
    X = (xtfloat*)&x[N-1];
    Y = (xtfloat*)&y[N-1];
    XT_LSXP(f, X, -(int)sz_f32); h = f;
    /* 24 cycles per pipeline stage in steady state with unroll=8 */
    for ( n=N-1; n>0; n-- ) {
        XT_LSXP(g, X, -(int)sz_f32);
        XT_MSUB_S(f, g, alpha);
        XT_SSXP(f, Y, -(int)sz_f32); f = g;
    }
    XT_MSUB_S(f, st, alpha);
    XT_SSI(f, Y, 0);
    return h;
#endif /* HAVE_FPU */
#else /* USE_REFERENCE_CODE */
    int n;
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(0==(N%2));
    for ( n=0; n<N; n++ ) {
        /* Note that x[] and y[] may refer to the same array! */
        float32_t s = x[n];
        y[n] = s-alpha*st; st = s;
    }
    return st;
#endif
} /* mfccf_preemph() */

#endif /* HAVE_VFPU || HAVE_FPU */
