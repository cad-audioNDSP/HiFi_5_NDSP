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

#define USE_REFERENCE_CODE  0
#define ALIGN_SIZE          (HIFI_SIMD_WIDTH)

#if USE_REFERENCE_CODE
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
#endif

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

void mfcc32x32_vecmpy( fract32 * z, const fract32 * x, const fract32 * y, int N )
{
#if !USE_REFERENCE_CODE
    const ae_int32x4 * restrict pX = (ae_int32x4*)x;
    const ae_int32x4 * restrict pY = (ae_int32x4*)y;
          ae_int32x4 * restrict pZ = (ae_int32x4*)z;
    int n;
    NASSERT_ALIGN(z, ALIGN_SIZE);
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(x!=y);
    /* 3 cycles per pipeline stage in steady state with unroll=2 */
    for ( n=0; n<(N>>2); n++ ) {
        ae_int32x2 a0, a1, b0, b1, c0, c1;
        AE_L32X2X2_IP(a0, a1, pX, sizeof(ae_int32x4));
        AE_L32X2X2_IP(b0, b1, pY, sizeof(ae_int32x4));
        /* Q(x+y-31) <- Qx*Qy - 31 w/ rounding */
        AE_MULF2P32X4RAS(c0, c1, a0, a1, b0, b1);
        AE_S32X2X2_IP(c0, c1, pZ, sizeof(ae_int32x4));
    } /* n */
    __Pragma("no_unroll");
    /* 2 cycles per pipeline stage in steady state with unroll=1 */
    for ( n=0; n<(N&3); n++ ) {
        ae_int32x2 a, b, c;
        AE_L32_IP(a, castxcc(ae_int32, pX), sizeof(ae_int32));
        AE_L32_IP(b, castxcc(ae_int32, pY), sizeof(ae_int32));
        /* Q(x+y-31) <- Qx*Qy - 31 w/ rounding */
        c = AE_MULFP32X2RAS(a, b);
        AE_S32_H_IP(c, castxcc(ae_int32, pZ), sizeof(ae_int32));
    }
#else
    int n;
    NASSERT_ALIGN(z, ALIGN_SIZE);
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT(x!=y);
    for ( n=0; n<N; n++ ) {
        /* Q(x+y-31) <- Qx*Qy - 31 w/ rounding */
        z[n] = mulf32x32ras(x[n], y[n]);
    }
#endif /* USE_REFERENCE_CODE */
} /* mfcc32x32_vecmpy() */
