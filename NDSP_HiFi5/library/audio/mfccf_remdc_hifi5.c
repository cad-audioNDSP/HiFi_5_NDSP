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
    C code optimized for HiF5 core
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

#if USE_REFERENCE_CODE || (!HAVE_VFPU && HAVE_FPU)
#include <math.h>
#endif

/*
 * Compute the mean value over the input argument x[N], then subtract it from 
 * each element of x[N] and store the result to the output argument y[N]. 
 * Input and output arguments may refer to the same array.
 * Note:
 *   For the fixed-point variant, input data should be properly scaled, 
 *   otherwise results may saturate. It is sufficient to ensure that the
 *   minimum number of redunsant sign bits over input data is non-zero to
 *   prevent overflow.
 * Input:
 *   N          Vectors size
 *   x[N]       Input vector
 * Output:
 *   y[N]       Output vector
 * Restrictions:
 *   N          Must be a multiple of 2
 *   x[N],y[N]  Must be 16-bytes aligned
 */

void mfccf_remdc( float32_t * y, const float32_t * x, int N )
{
#if !USE_REFERENCE_CODE
    /*
     * Reference code:
     *   float32_t s, t, acch, accl, mean;
     *   int n;
     *   NASSERT_ALIGN(y, ALIGN_SIZE);
     *   NASSERT_ALIGN(x, ALIGN_SIZE);
     *   NASSERT(0==(N%2));
     *   acch = accl = 0;
     *   for ( n=0; n<N; n++ ) {
     *       s = x[n];
     *       t = s+accl+acch;
     *       // It is essential to accumulate the sum in extended precision.
     *       accl = -(t-acch-accl-s);
     *       acch = t;
     *   }
     *   mean = acch/N;
     *   for ( n=0; n<N; n++ ) {
     *       y[n] = x[n]-mean;
     *   }
     */
#if HAVE_VFPU
    const xtfloatx4 * restrict pX;
    xtfloatx4 * restrict pY;
    ae_valignx2 X_u;
    xtfloat suml, sumh, mean;
    xtfloatx2 accl, acch;
    int n, headLen;
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(0==(N%2));
    /*
     * Compute the sum over input vector: acc <- sum(x)
     */
    if (N>=16) {
        static const int32_t ALIGN(16) seq[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
        xtbool2 b0, b1, b2, b3, b4, b5, b6, b7;
        ae_int32x2 seq0, seq1, seq2, seq3, seq4, seq5, seq6, seq7;
        xtfloatx2 acch0, acch1, acch2, acch3, acch4, acch5, acch6, acch7;
        xtfloatx2 accl0, accl1, accl2, accl3, accl4, accl5, accl6, accl7;
        xtfloatx2 r0, r1, r2, r3, r4, r5, r6, r7;
        xtfloatx2 s0, s1, s2, s3, s4, s5, s6, s7;
        xtfloatx2 t0, t1, t2, t3, t4, t5, t6, t7;
        AE_L32X2X2_I(seq0, seq1, (const ae_int32x4*)seq, 0*sizeof(ae_int32x4));
        AE_L32X2X2_I(seq2, seq3, (const ae_int32x4*)seq, 1*sizeof(ae_int32x4));
        AE_L32X2X2_I(seq4, seq5, (const ae_int32x4*)seq, 2*sizeof(ae_int32x4));
        AE_L32X2X2_I(seq6, seq7, (const ae_int32x4*)seq, 3*sizeof(ae_int32x4));
        /* Mask useful elements for the prologue. */
        headLen = ((N-1)&15)+1;
        b0 = AE_LT32(seq0, headLen); b1 = AE_LT32(seq1, headLen);
        b2 = AE_LT32(seq2, headLen); b3 = AE_LT32(seq3, headLen);
        b4 = AE_LT32(seq4, headLen); b5 = AE_LT32(seq5, headLen);
        b6 = AE_LT32(seq6, headLen); b7 = AE_LT32(seq7, headLen);
        pX = (xtfloatx4*)x;
        AE_LSX2X2_IP(s0, s1, pX, 4*sz_f32); XT_MOVF_SX2(s0, 0, b0); XT_MOVF_SX2(s1, 0, b1);
        AE_LSX2X2_IP(s2, s3, pX, 4*sz_f32); XT_MOVF_SX2(s2, 0, b2); XT_MOVF_SX2(s3, 0, b3);
        AE_LSX2X2_IP(s4, s5, pX, 4*sz_f32); XT_MOVF_SX2(s4, 0, b4); XT_MOVF_SX2(s5, 0, b5);
        AE_LSX2X2_IP(s6, s7, pX, 4*sz_f32); XT_MOVF_SX2(s6, 0, b6); XT_MOVF_SX2(s7, 0, b7);
        acch0 = s0; accl0 = 0; acch1 = s1; accl1 = 0;
        acch2 = s2; accl2 = 0; acch3 = s3; accl3 = 0;
        acch4 = s4; accl4 = 0; acch5 = s5; accl5 = 0;
        acch6 = s6; accl6 = 0; acch7 = s7; accl7 = 0;
        pX = (xtfloatx4*)XT_ADDX4(headLen, (uintptr_t)x);
        X_u = AE_LA128_PP(pX);
        /* 20 cycles per pipeline stage in steady state with unroll=1 */
        for ( n=0; n<((N-headLen)>>4); n++ ) {
            AE_LASX2X2_IP(s0, s1, X_u, pX);
            AE_LASX2X2_IP(s2, s3, X_u, pX);
            AE_LASX2X2_IP(s4, s5, X_u, pX);
            AE_LASX2X2_IP(s6, s7, X_u, pX);
            /* t <- s+accl+acch */
            ADD_SX2X2(t0, t1, s0, s1, accl0, accl1); ADD_SX2X2(t0, t1, t0, t1, acch0, acch1);
            ADD_SX2X2(t2, t3, s2, s3, accl2, accl3); ADD_SX2X2(t2, t3, t2, t3, acch2, acch3);
            ADD_SX2X2(t4, t5, s4, s5, accl4, accl5); ADD_SX2X2(t4, t5, t4, t5, acch4, acch5);
            ADD_SX2X2(t6, t7, s6, s7, accl6, accl7); ADD_SX2X2(t6, t7, t6, t7, acch6, acch7);
            /* accl <- -(t-acch-accl-s); acch <- t */
            SUB_SX2X2(r0, r1, acch0, acch1, t0, t1); ADD_SX2X2(r0, r1, accl0, accl1, r0, r1);
            SUB_SX2X2(r2, r3, acch2, acch3, t2, t3); ADD_SX2X2(r2, r3, accl2, accl3, r2, r3);
            SUB_SX2X2(r4, r5, acch4, acch5, t4, t5); ADD_SX2X2(r4, r5, accl4, accl5, r4, r5);
            SUB_SX2X2(r6, r7, acch6, acch7, t6, t7); ADD_SX2X2(r6, r7, accl6, accl7, r6, r7);
            ADD_SX2X2(accl0, accl1, s0, s1, r0, r1); acch0 = t0; acch1 = t1;
            ADD_SX2X2(accl2, accl3, s2, s3, r2, r3); acch2 = t2; acch3 = t3;
            ADD_SX2X2(accl4, accl5, s4, s5, r4, r5); acch4 = t4; acch5 = t5;
            ADD_SX2X2(accl6, accl7, s6, s7, r6, r7); acch6 = t6; acch7 = t7;
        } /* n */
        {
            xtfloatx2 f, g;
            xtbool2 bm2;
            BMAXNUMABS_SX2(bm2, f, acch0, acch1);
            MOV_SX2X2(f, g, acch1, acch0);
            XT_MOVT_SX2(f, acch0, bm2); /* Max */
            XT_MOVT_SX2(g, acch1, bm2); /* Min */
            acch0 = XT_ADD_SX2(f, XT_ADD_SX2(g, XT_ADD_SX2(accl0, accl1)));
            accl0 = XT_ADD_SX2(accl0, XT_ADD_SX2(accl1, XT_ADD_SX2(g, XT_SUB_SX2(f, acch0))));

            BMAXNUMABS_SX2(bm2, f, acch2, acch3);
            MOV_SX2X2(f, g, acch3, acch2);
            XT_MOVT_SX2(f, acch2, bm2); /* Max */
            XT_MOVT_SX2(g, acch3, bm2); /* Min */
            acch1 = XT_ADD_SX2(f, XT_ADD_SX2(g, XT_ADD_SX2(accl2, accl3)));
            accl1 = XT_ADD_SX2(accl2, XT_ADD_SX2(accl3, XT_ADD_SX2(g, XT_SUB_SX2(f, acch1))));

            BMAXNUMABS_SX2(bm2, f, acch4, acch5);
            MOV_SX2X2(f, g, acch5, acch4);
            XT_MOVT_SX2(f, acch4, bm2); /* Max */
            XT_MOVT_SX2(g, acch5, bm2); /* Min */
            acch2 = XT_ADD_SX2(f, XT_ADD_SX2(g, XT_ADD_SX2(accl4, accl5)));
            accl2 = XT_ADD_SX2(accl4, XT_ADD_SX2(accl5, XT_ADD_SX2(g, XT_SUB_SX2(f, acch2))));

            BMAXNUMABS_SX2(bm2, f, acch6, acch7);
            MOV_SX2X2(f, g, acch7, acch6);
            XT_MOVT_SX2(f, acch6, bm2); /* Max */
            XT_MOVT_SX2(g, acch7, bm2); /* Min */
            acch3 = XT_ADD_SX2(f, XT_ADD_SX2(g, XT_ADD_SX2(accl6, accl7)));
            accl3 = XT_ADD_SX2(accl6, XT_ADD_SX2(accl7, XT_ADD_SX2(g, XT_SUB_SX2(f, acch3))));

            BMAXNUMABS_SX2(bm2, f, acch0, acch1);
            MOV_SX2X2(f, g, acch1, acch0);
            XT_MOVT_SX2(f, acch0, bm2); /* Max */
            XT_MOVT_SX2(g, acch1, bm2); /* Min */
            acch0 = XT_ADD_SX2(f, XT_ADD_SX2(g, XT_ADD_SX2(accl0, accl1)));
            accl0 = XT_ADD_SX2(accl0, XT_ADD_SX2(accl1, XT_ADD_SX2(g, XT_SUB_SX2(f, acch0))));

            BMAXNUMABS_SX2(bm2, f, acch2, acch3);
            MOV_SX2X2(f, g, acch3, acch2);
            XT_MOVT_SX2(f, acch2, bm2); /* Max */
            XT_MOVT_SX2(g, acch3, bm2); /* Min */
            acch1 = XT_ADD_SX2(f, XT_ADD_SX2(g, XT_ADD_SX2(accl2, accl3)));
            accl1 = XT_ADD_SX2(accl2, XT_ADD_SX2(accl3, XT_ADD_SX2(g, XT_SUB_SX2(f, acch1))));

            BMAXNUMABS_SX2(bm2, f, acch0, acch1);
            MOV_SX2X2(f, g, acch1, acch0);
            XT_MOVT_SX2(f, acch0, bm2); /* Max */
            XT_MOVT_SX2(g, acch1, bm2); /* Min */
            acch = XT_ADD_SX2(f, XT_ADD_SX2(g, XT_ADD_SX2(accl0, accl1)));
            accl = XT_ADD_SX2(accl0, XT_ADD_SX2(accl1, XT_ADD_SX2(g, XT_SUB_SX2(f, acch))));
        }
    } else { /* (N<16) */
        xtfloatx2 s, t;
        acch = accl = 0;
        pX = (xtfloatx4*)x;
        /* 15 cycles per pipeline stage in steady state with unroll=1 */
        for ( n=0; n<(N>>1); n++ ) {
            XT_LSX2IP(s, castxcc(xtfloatx2, pX), 2*sz_f32);
            /* t <- s+accl+acch */
            t = XT_ADD_SX2(XT_ADD_SX2(s, accl), acch);
            /* accl <- -(t-acch-accl-s); acch <- t */
            accl = XT_ADD_SX2(XT_ADD_SX2(XT_SUB_SX2(acch, t), accl), s);
            acch = t;
        } /* n */
    } /* (N<16) */
    /*
     * Reduce 2-element accumulators to 1-element sum.
     */
    {
        xtfloat f, g, h, r;
        xtbool2 bm2;
        xtbool bm;
        h = XT_LOW_S(acch); r = XT_HIGH_S(acch);
        BMAXNUMABS_S(bm2, f, h, r); bm = xtbool2_extract_0(bm2);
        f = r; XT_MOVT_S(f, h, bm); /* Max */
        g = r; XT_MOVF_S(g, h, bm); /* Min */
        h = XT_LOW_S(accl); r = XT_HIGH_S(accl);
        sumh = XT_ADD_S(f, XT_ADD_S(g, XT_ADD_S(h, r)));
        suml = XT_ADD_S(r, XT_ADD_S(h, XT_ADD_S(g, XT_SUB_S(f, sumh))));
    }
    /*
     * Compute mean <- sum/N
     */
    {
        xtfloat f, h, r;
        /* r <- 1.f/N */
        r = XT_RECIP_S(XT_FLOAT_S(N, 0));
        /* f <- (suml*r)+sumh*r */
        f = XT_MUL_S(suml, r); XT_MADD_S(f, sumh, r);
        /* h <- (sumh-f*N)+suml */
        h = sumh; XT_MSUB_S(h, f, XT_FLOAT_S(N, 0)); h = XT_ADD_S(h, suml);
        /* f <- f+r*h */
        XT_MADD_S(f, r, h); mean = f;
    }
    /* If commented out, the division code is duplicated and merged into prologues of
     * the two loops below. */
    __Pragma("no_reorder"); 
    /*
     * Subract the mean: y <- x-mean
     */
    pX = (xtfloatx4*)x;
    pY = (xtfloatx4*)y;
    __Pragma("loop_count factor=4");
    /* 4 cycles per pipeline stage in steady state with unroll=4 */
    for ( n=0; n<(N>>4)<<2; n++ ) {
        xtfloatx2 s0, s1;
        AE_LSX2X2_IP(s0, s1, pX, 4*sz_f32);
        SUB_SX2X2(s0, s1, s0, s1, mean, mean);
        AE_SSX2X2_IP(s0, s1, pY, 4*sz_f32);
    } /* n */
    __Pragma("no_unroll");
    /* 3 cycles per pipeline stage in steady state with unroll=1 */
    for ( n=0; n<((N&15)>>1); n++ ) {
        xtfloatx2 s;
        XT_LSX2IP(s, castxcc(xtfloatx2, pX), 2*sz_f32);
        s = XT_SUB_SX2(s, mean);
        XT_SSX2IP(s, castxcc(xtfloatx2, pY), 2*sz_f32);
    } /* n */
#elif HAVE_FPU
    /*
     * NOTE: optimized mfccf_remdc is not available for scalar FPU, use the reference code instead
     */
    float32_t r, s, t, acch, accl;
    int n;
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(0==(N%2));
    acch = accl = 0;
    for ( n=0; n<N; n++ ) {
        s = x[n];
        t = s+accl+acch;
        /* It is essential to accumulate the sum in extended precision. */
        accl = -(t-acch-accl-s);
        acch = t;
    }
    r = 1.f/N;
    t = fmaf(acch, r, accl*r);
    s = fmaf(-t, N, acch)+accl;
    t = fmaf(r, s, t);
    for ( n=0; n<N; n++ ) {
        y[n] = x[n]-t;
    }
#endif /* HAVE_FPU */
#else /* USE_REFERENCE_CODE */
    float32_t r, s, t, acch, accl;
    int n;
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(0==(N%2));
    acch = accl = 0;
    for ( n=0; n<N; n++ ) {
        s = x[n];
        t = s+accl+acch;
        /* It is essential to accumulate the sum in extended precision. */
        accl = -(t-acch-accl-s);
        acch = t;
    }
    r = 1.f/N;
    t = fmaf(acch, r, accl*r);
    s = fmaf(-t, N, acch)+accl;
    t = fmaf(r, s, t);
    for ( n=0; n<N; n++ ) {
        y[n] = x[n]-t;
    }
#endif
} /* mfccf_remdc() */

#endif /* HAVE_VFPU || HAVE_FPU */
