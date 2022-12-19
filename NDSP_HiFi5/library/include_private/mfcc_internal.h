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
    Internal declarations
  IntegrIT, 2006-2018
*/

#ifndef __MFCC_INTERNAL_H
#define __MFCC_INTERNAL_H

#include <string.h>

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_audio.h"

/* Fast copying of 32-bit data. */
inline_ void ATTRIBUTE_ALWAYS_INLINE mfcc_copy_32b( void * dst, const void * src, int len ) 
{
    int n, headLen;
    __Pragma("no_reorder");
    if (len<4) {
        const ae_int32 * _src = (ae_int32*)src;
        ae_int32 * _dst = (ae_int32*)dst;
        __Pragma("concurrent");
        for ( n=0; n<len; n++ ) {
            ae_int32x2 t;
            AE_L32_IP(t, _src, sizeof(*_src));
            AE_S32_H_IP(t, _dst, sizeof(*_dst));
        }
    } else {
        const ae_int32x4 * _src = (ae_int32x4*)src;
        ae_int32x4 * _dst = (ae_int32x4*)dst;
        ae_valignx2 src_va = AE_LA128_PP(_src);
        ae_valignx2 dst_va = AE_ZALIGN128();
        ae_int32x2 t0, t1;
        AE_LA32X2X2_IP(t0, t1, src_va, _src);
        AE_SA32X2X2_IP(t0, t1, dst_va, _dst);
        AE_SA128POS_FP(dst_va, _dst);
        headLen = ((len-1)&3)+1; len -= headLen;
        _src = (ae_int32x4*)XT_ADDX4(headLen, (uintptr_t)src);
        _dst = (ae_int32x4*)XT_ADDX4(headLen, (uintptr_t)dst);
        src_va = AE_LA128_PP(_src);
        __Pragma("concurrent");
        for ( n=0; n<(len>>2); n++ ) {
            AE_LA32X2X2_IP(t0, t1, src_va, _src);
            AE_SA32X2X2_IP(t0, t1, dst_va, _dst);
        }
        AE_SA128POS_FP(dst_va, _dst);
    }
    __Pragma("no_reorder");
} /* mfcc_copy_32b_aligned() */

#define MFCC_COPY_32b(dst, src, len)      mfcc_copy_32b(dst, src, len)

inline_ void ATTRIBUTE_ALWAYS_INLINE mfcc_zero_32b( void * dst, int len )
{
    int n, headLen;
    __Pragma("no_reorder");
    if (len<4) {
        ae_int32 * _dst = (ae_int32*)dst;
        for ( n=0; n<len; n++ ) AE_S32_L_IP(0, _dst, sizeof(*_dst));
    } else {
        ae_int32x4 * _dst = (ae_int32x4*)dst;
        ae_valignx2 va = AE_ZALIGN128();
        AE_SA32X2X2_IP(0, 0, va, _dst);
        AE_SA128POS_FP(va, _dst);
        headLen = ((len-1)&3)+1; len -= headLen;
        _dst = (ae_int32x4*)XT_ADDX4(headLen, (uintptr_t)dst);
        for ( n=0; n<(len>>2); n++ ) AE_SA32X2X2_IP(0, 0, va, _dst);
        AE_SA128POS_FP(va, _dst);
    }
    __Pragma("no_reorder");
} /* mfcc_zero_32b() */

#define MFCC_ZERO_32b(dst, len)      mfcc_zero_32b(dst, len)

#define LOGMEL_PARAMS_FROM_MFCC_PARAMS(p_logmelParams, p_mfccParams) { \
    memset(p_logmelParams, 0, sizeof(*(p_logmelParams)));              \
    (p_logmelParams)->Fs           = (p_mfccParams)->Fs;               \
    (p_logmelParams)->fftSize      = (p_mfccParams)->fftSize;          \
    (p_logmelParams)->mfbLowFreqQ8 = (p_mfccParams)->mfbLowFreqQ8;     \
    (p_logmelParams)->mfbUppFreqQ8 = (p_mfccParams)->mfbUppFreqQ8;     \
    (p_logmelParams)->mfbBandNum   = (p_mfccParams)->mfbBandNum;       \
    (p_logmelParams)->opt          = (p_mfccParams)->opt; }

/* Helper to retrieve individual flags from options bitmask. */
#define MFCC_OPT_DC_MEAN_MASK            (MFCC_OPT_DC_MEAN_REMOVE|MFCC_OPT_DC_MEAN_DONT_REMOVE)
#define MFCC_OPT_PREEMPH_MASK            (MFCC_OPT_PREEMPH_FRAMEBYFRAME|MFCC_OPT_PREEMPH_CONTINUOUS)
#define MFCC_OPT_DCT_MASK                (MFCC_OPT_DCT_NORMALIZED|MFCC_OPT_DCT_ORTHOGONAL)
#define MFCC_OPT_EQ(opt,feature,value)   (((opt)&(MFCC_OPT_##feature##_MASK))==MFCC_OPT_##feature##_##value)

/* Instance pointer validation numbers. */
#define MFCC32X32_MAGIC    0x5dcb387a
#define MFCCF_MAGIC        0x112f2bc5

/* MFCC extractor performs the DCT via matrix-vector multiplication. The DSP Library
 * provides regular and fast variants of matrix-vector mupltiply routines. If the flag
 * defined below is non-zero, then internal data arrays are zero-padded so that the
 * fast routine is selected irrespective of MFCC configuration, otherwise the regular
 * variant is used whenever restrictions of the fast routine are not met. */
#define MFCC32X32_ALWAYS_USE_MTXVECMPY_FAST    1 
#define MFCCF_ALWAYS_USE_MTXVECMPY_FAST        1

/* MFCC features extractor instance structure, 32-bit fixed-point variant. */
typedef struct mfcc32x32_tag
{
    uint32_t             magic;        /* Instance pointer validation number                                */
    mfcc_params_t        params;       /* MFCC extractor parameters                                         */
    mfcc32x32_callback_t callback;     /* User-provided callback functions                                  */
    logmel32x32_handle_t logmel;       /* Log mel filterbank handle                                         */
    const fract32      * stftWeights;  /* STFT window weights (Q31), or NULL if no window function          */
    fract32            * stftFrame;    /* STFT signal frame, zero-padded up to the FFT size                 */
    fract32              preemphState; /* Pre-emphasis filter state from the last STFT frame                */
    int                  dctRows;      /* DCT matrix vertical dimension                                     */
    int                  dctCols;      /* DCT matrix horizontal dimension                                   */
    const fract32      * dctMtx;       /* Precomputed DCT matrix, Q(31-dctExp)                              */
    int                  dctExp;       /* Precomputed DCT matrix exponent                                   */
    const fract32      * lifterCoefs;  /* Lifter coeffs, Q(31-lifterExp), or NULL if lifter is not selected */
    int                  lifterExp;    /* Lifter coeffs exponent                                            */
    /* Matrix-vector multiply routine for the DCT. */
    void (*mtxvecmpy)(fract32 * z, const fract32 * x, const fract32 * y, int M, int N, int lsh);
} mfcc32x32_t;

/* MFCC features extractor instance structure, single precision floating-point variant. */
typedef struct mfccf_tag
{
    uint32_t             magic;        /* Instance pointer validation number                 */
    mfcc_params_t        params;       /* MFCC extractor parameters                          */
    mfccf_callback_t     callback;     /* User-provided callback functions                   */
    logmelf_handle_t     logmel;       /* Log mel filterbank handle                          */
    const float32_t    * stftWeights;  /* STFT window weights, or NULL if no window function */
    float32_t          * stftFrame;    /* STFT signal frame, zero-padded up to the FFT size  */
    float32_t            preemphState; /* Pre-emphasis filter state from the last STFT frame */
    int                  dctRows;      /* DCT matrix vertical dimension                      */
    int                  dctCols;      /* DCT matrix horizontal dimension                    */
    const float32_t    * dctMtx;       /* Precomputed DCT matrix                             */
    const float32_t    * lifterCoefs;  /* Lifter coeffs, or NULL if lifter is not selected   */
    /* Matrix-vector multiply routine for the DCT. */
    void (*mtxvecmpy)(float32_t * z, const float32_t * x, const float32_t * y, int M, int N);
} mfccf_t;

/* Polynomial coeffs for 32-bit sin/cos evaluation. */
extern const fract32 mfcc32x32_polysin_tbl[6*8];

/* 
 * Calculate the DCT matrix, as defined in The HTK Book, see [1] p.77 Eq.(5.14).
 * MATLAB reference code:
 *   dctMtx = (2/N)^0.5*cos(pi/N*(0:M-1)'*((1:N)-0.5));
 *   if orthNorm, dctMtx(1,:) = dctMtx(1,:)/2^0.5; end;
 * For efficiency reasons, the DCT matrix may be stored as a submatrix of a larger
 * matrix, as specified by input arguemnts M,N (DCT dimensions) and rows,cols (matrix
 * dimensions): 0<M<=rows, 0<N<=cols. Padding elements of the larger matrix are zeroed.
 * Input:
 * M                  Vertical dimension of the DCT transform
 * N                  Horizontal dimension of the DCT transform
 * rows               Vertical dimension of the storage matrix
 * cols               Horizontal dimension of the storage matrix
 * orthNorm           If non-zero, then in addition to regular normalization of the whole matrix
 *                    by a factor of (2/N)^0.5, its first row is further multiplied by 1/2^0.5, 
 *                    so that the resulting DCT matrix becomes orthogonally normalized.
 * Output:
 * dctMtx[rows*cols]  DCT matrix of size MxN, padded with zeros to form a storage matrix
 *                    of size rows x cols; Q(31-dctExp) for 32x32
 * dctExp             Returned value (32x32 only), base-2 exponent of the DCT matrix
 * Restrictions:
 * dctMtx[]           Aligned by 8-bytes
 * 0<M<=rows
 * 0<N<=cols 
 */

int  mfcc32x32_compDctMatrix( fract32   * restrict dctMtx, int M, int N, int rows, int cols, int orthNorm );
void mfccf_compDctMatrix    ( float32_t * restrict dctMtx, int M, int N, int rows, int cols, int orthNorm );

/*
 * Compute the lifter coefficients, as defined in The HTK Book, see [1] p.75 Eq.(5.12).
 * MATLAB reference code for the liftering operation:
 *   lifter = @(c,L)((1+L/2*sin(pi/L*(0:length(c)-1)')).*c);
 * Input:
 * L               Lifter parameter
 * N               Number of coeffs to compute (equals the number of cepstral features)
 * Output:
 * lifterCoefs[N]  Lifter coeffs, Q(31-lifterExp) for 32x32
 * lifterExp       Returned value (32x32 only), base-2 block exponent of lifter coeffs.
 * Restrictions:
 * lifterCoefs[]   Aligned by 8-bytes
 * N>0
 */

int  mfcc32x32_compLifterCoefs( fract32   * restrict lifterCoefs, int L, int N );
void mfccf_compLifterCoefs    ( float32_t * restrict lifterCoefs, int L, int N );

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
 *   x[N],y[N]  Must not overlap, and must be aligned by 8-bytes
 *   z[N]       Must be aligned by 8-bytes
 */

void mfcc32x32_vecmpy( fract32   * z, const fract32   * x, const fract32   * y, int N );
void mfccf_vecmpy    ( float32_t * z, const float32_t * x, const float32_t * y, int N );

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
 *   x[N],y[N]  Must be 8-bytes aligned
 */

void mfcc32x32_remdc( fract32   * y, const fract32   * x, int N );
void mfccf_remdc    ( float32_t * y, const float32_t * x, int N );

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
 *   x[N],y[N]  Must be 8-bytes aligned
 */

fract32   mfcc32x32_preemph( fract32   * y, const fract32   * x, fract16   alpha, fract32   st, int N );
float32_t mfccf_preemph    ( float32_t * y, const float32_t * x, float32_t alpha, float32_t st, int N );

#endif /* __MFCC_INTERNAL_H */
