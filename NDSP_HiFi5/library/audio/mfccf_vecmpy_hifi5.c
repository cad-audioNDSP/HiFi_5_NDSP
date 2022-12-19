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

#define ALIGN_SIZE     (HIFI_SIMD_WIDTH)
#define sz_f32         sizeof(float32_t)

/*
 * Pairwise multiplication of input vector arguments x[n] and y[n], with 
 * resulting values stored to the output argument z[N]. Input arguments 
 * x[N] nad y[N] must be distinct, but the output argument z[N] may refer
 * to any of input arguments.
 * Input:
 *   N          Vectors size
 *   x[N]       First multiplicand; Qx for 32x32
 *   y[N]       Second multiplicand; Qy for 32x32
 * Output:
 *   z[N]       Pairwise multiplication results; Q(x+y-31) for 32x32
 * Restrictions:
 *   x[N],y[N]  Must not overlap, and must be aligned by 16-bytes
 *   z[N]       Must be aligned by 16-bytes
 */

void mfccf_vecmpy( float32_t * z, const float32_t * x, const float32_t * y, int N )
{
#if HAVE_VFPU
    const xtfloatx4 * restrict pX = (xtfloatx4*)x;
    const xtfloatx4 * restrict pY = (xtfloatx4*)y;
          xtfloatx4 * restrict pZ = (xtfloatx4*)z;
    int n;
    NASSERT_ALIGN(z, ALIGN_SIZE);
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(x!=y);
    __Pragma("loop_count factor=4");
    /* 6 cycles per pipeline stage in steady state with unroll=4 */
    for ( n=0; n<((N>>4)<<2); n++ ) {
        xtfloatx2 a0, a1, b0, b1, c0, c1;
        AE_LSX2X2_IP(a0, a1, pX, 4*sz_f32);
        AE_LSX2X2_IP(b0, b1, pY, 4*sz_f32);
        MUL_SX2X2(c0, c1, a0, a1, b0, b1);
        AE_SSX2X2_IP(c0, c1, pZ, 4*sz_f32);
    } /* n */
    __Pragma("no_unroll");
    /* 4 cycles per pipeline stage in steady state with unroll=1 */
    for ( n=0; n<((N>>1)&7); n++ ) {
        xtfloatx2 a, b, c;
        XT_LSX2IP(a, castxcc(xtfloatx2, pX), 2*sz_f32);
        XT_LSX2IP(b, castxcc(xtfloatx2, pY), 2*sz_f32);
        c = XT_MUL_SX2(a, b);
        XT_SSX2IP(c, castxcc(xtfloatx2, pZ), 2*sz_f32);
    } /* n */
    if (N&1) {
        xtfloat a, b, c;
        a = XT_LSI((xtfloat*)pX, 0);
        b = XT_LSI((xtfloat*)pY, 0);
        c = XT_MUL_S(a, b);
        XT_SSI(c, (xtfloat*)pZ, 0);
    } /* (N&1) */
#else /* !HAVE_VFPU */
    int n;
    NASSERT_ALIGN(z, ALIGN_SIZE);
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(x!=y);
    __Pragma("aligned(x, 16)");
    __Pragma("aligned(y, 16)");
    __Pragma("aligned(z, 16)");
    __Pragma("concurrent");
    /* SFPU: 16 cycles per pipeline stage in steady state with unroll=4 (w/o vectorization) */
    for ( n=0; n<N; n++ ) {
        z[n] = x[n]*y[n];
    }
#endif /* !HAVE_VFPU */
} /* mfccf_vecmpy() */

#endif /* HAVE_VFPU || HAVE_FPU */
