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
  NatureDSP Signal Processing Library. FFT part
    Complex-valued FFT stages with butterflies radix-4, radix-8
    with static data scaling: 32-bit data, 16-bit twiddle factors
    C code optimized for HiFi4
  IntegrIT, 2006-2019
*/
#include "NatureDSP_types.h"
#include "common.h"
#include "fft_32x16_stages.h"
#include "fft_x16_common.h"

/* radix-4 butterfly with normalization */
#define DFT4X1RNG(x0, x1, x2, x3)\
{   \
    ae_int32x2 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG32_L(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32_L(s1, d1, x1, x3); \
    AE_ADDANDSUB32S(x0, x2, s0, s1);   \
    AE_ADDANDSUB32JS(x3, x1, d0, d1);   \
}

ALIGN(32) static const int32_t __fft16_tw1[] =
{
    (int32_t)0x7FFFFFFF, (int32_t)0x00000000, (int32_t)0x7FFFFFFF, (int32_t)0x00000000, (int32_t)0x7FFFFFFF, (int32_t)0x00000000,
    (int32_t)0x7641AF3D, (int32_t)0xCF043AB3, //t5,
    (int32_t)0x5A82799A, (int32_t)0xA57D8666, //t9,
    (int32_t)0x30FBC54D, (int32_t)0x89BE50C3, //t13,
    (int32_t)0x5A82799A, (int32_t)0xA57D8666, //t6,
    (int32_t)0x00000000, (int32_t)0x80000000, //t10,
    (int32_t)0xA57D8666, (int32_t)0xA57D8666, //t14,
    (int32_t)0x30FBC54D, (int32_t)0x89BE50C3, //t7,
    (int32_t)0xA57D8666, (int32_t)0xA57D8666, //t11,
    (int32_t)0x89BE50C3, (int32_t)0x30FBC54D, //t15,
};

int ifft_32x16_stage_last_scl3_DFT16(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    const ae_int32x2 * restrict px0;
    ae_int64 * restrict py0;
    ae_int32x2 * restrict ptw = (ae_int32x2*)(__fft16_tw1 + 6);

    const int stride = (N >> 4);
    int j;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);

    px0 = (const ae_int32x2 *)x;
    py0 = (ae_int64 *)y;

    ae_int32x2 x0, x1, x2, x3;
    ae_int32x2 x4, x5, x6, x7;
    ae_int32x2 x8, x9, x10, x11;
    ae_int32x2 x12, x13, x14, x15;
    ae_int32x2 t5, t9, t13;
    ae_int32x2 t6, t10, t14;
    ae_int32x2 t7, t11, t15;

    AE_L32X2_IP(t5, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t9, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t13, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t6, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t10, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t14, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t7, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t11, ptw, sizeof(complex_fract32));
    AE_L32X2_IP(t15, ptw, sizeof(complex_fract32));

    WUR_AE_SAR(2);

    __Pragma("loop_count min=1");
    for (j = 0; j < stride; j++)
    {
        AE_L32X2_XP(x0, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x1, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x2, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x3, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x4, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x5, px0, stride*sizeof(complex_fract32)); 
        AE_L32X2_XP(x6, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x7, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x8, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x9, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x10, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x11, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x12, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x13, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x14, px0, stride*sizeof(complex_fract32));
        AE_L32X2_XP(x15, px0, (1 - 15 * stride)*sizeof(complex_fract32));

        DFT4X1RNG(x0, x4, x8, x12);
        DFT4X1RNG(x1, x5, x9, x13);
        DFT4X1RNG(x2, x6, x10, x14);
        DFT4X1RNG(x3, x7, x11, x15);

        x5 = AE_MULFC32RAS(x5, t5);
        x6 = AE_MULFC32RAS(x6, t6);
        x7 = AE_MULFC32RAS(x7, t7);
        x9 = AE_MULFC32RAS(x9, t6);
        x10 = AE_MULFC32RAS(x10, t10);
        x11 = AE_MULFC32RAS(x11, t11);
        x13 = AE_MULFC32RAS(x13, t7);
        x14 = AE_MULFC32RAS(x14, t11);
        x15 = AE_MULFC32RAS(x15, t15);

        DFT4X1RNG(x0, x1, x2, x3);
        DFT4X1RNG(x4, x5, x6, x7);
        DFT4X1RNG(x8, x9, x10, x11);
        DFT4X1RNG(x12, x13, x14, x15);

        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x0), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x4), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x8), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x12), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x1), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x5), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x9), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x13), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x2), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x6), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x10), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x14), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x3), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x7), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x11), py0, stride*sizeof(complex_fract32));
        AE_S64_XP(AE_MOVINT64_FROMINT32X2(x15), py0, (1 - 15 * stride)*sizeof(complex_fract32));
    }
    return 4;
} /* ifft_32x16_stage_last_scl3_DFT16() */
