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
    with dynamic data scaling: 32-bit data, 16-bit twiddle factors
    C code optimized for HiFi5
*/
#include "NatureDSP_types.h"
#include "common.h"
#include "fft_32x16_stages.h"
#include "fft_x16_common.h"

/* radix-4 butterfly with normalization */
#define DFT4X1RNG(x0, x1, x2, x3)\
{   \
    ae_int32x2 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG32(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32(s1, d1, x1, x3); \
    d1 = AE_MUL32JS(d1);               \
    AE_ADDANDSUB32S(x0, x2, s0, s1);   \
    AE_ADDANDSUB32S(x3, x1, d0, d1);   \
}

/* Table of twiddle factors for radix-8 butterfly */
ALIGN(32) static const int16_t __fft8_tw1[] =
{
  (int16_t)0x4000, (int16_t)0x4000,/* first 2 numbers are used for scaling */
  (int16_t)0x2D41, (int16_t)0xD2BF,
  (int16_t)0x0000, (int16_t)0xC000,
  (int16_t)0xD2BF, (int16_t)0xD2BF,
};

/*
 *  32x16 FFT/IFFT last stage Radix 8, scalingOption=2
 */
int fft_32x16_stage_last_scl2_DFT8(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    int shift;
    int j;
    const ae_int32x2 * restrict px0;
    const ae_int32x2 * restrict px1;
          ae_int32x2 * restrict py0;
          ae_int32x2 * restrict py1;
    const ae_int16x4 * restrict ptwd;
    ae_int16x4 tw1, tw2;

    const int stride = (N >> 3);
    const int min_shift = 3;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(*bexp >= 0);
    NASSERT(v[0] == stride);

    shift = XT_MAX(0, min_shift - *bexp);
    WUR_AE_SAR(shift);
    px0 = (ae_int32x2 *)x;
    px1 = px0 + 4*stride;
    py0 = (ae_int32x2 *)y;
    py1 = py0 + 4*stride;
    ptwd = (const ae_int16x4 *)__fft8_tw1;
    tw1 = AE_L16X4_I(ptwd, 0*sizeof(ae_int16x4));
    tw2 = AE_L16X4_I(ptwd, 1*sizeof(ae_int16x4));

    __Pragma("loop_count min=3");
    for (j = 0; j < (stride>>1); j++)
    {
        ae_int32x2 x00, x10, x20, x30;
        ae_int32x2 x40, x50, x60, x70;
        ae_int32x2 x01, x11, x21, x31;
        ae_int32x2 x41, x51, x61, x71;

        AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), stride*sizeof(ae_int32x2));
        AE_L32X2X2_XP(x10, x11, castxcc(ae_int32x4, px0), stride*sizeof(ae_int32x2));
        AE_L32X2X2_XP(x20, x21, castxcc(ae_int32x4, px0), stride*sizeof(ae_int32x2));
        AE_L32X2X2_XP(x30, x31, castxcc(ae_int32x4, px0), (2 - 3 * stride)*sizeof(ae_int32x2));
        AE_L32X2X2_XP(x40, x41, castxcc(ae_int32x4, px1), stride*sizeof(ae_int32x2));
        AE_L32X2X2_XP(x50, x51, castxcc(ae_int32x4, px1), stride*sizeof(ae_int32x2));
        AE_L32X2X2_XP(x60, x61, castxcc(ae_int32x4, px1), stride*sizeof(ae_int32x2));
        AE_L32X2X2_XP(x70, x71, castxcc(ae_int32x4, px1), (2 - 3 * stride)*sizeof(ae_int32x2));

        DFT4X1RNG(x00, x20, x40, x60);
        DFT4X1RNG(x10, x30, x50, x70);
        x30 = AE_MULFC32X16RAS_L(x30, tw1);
        x50 = AE_MULFC32X16RAS_H(x50, tw2);
        x70 = AE_MULFC32X16RAS_L(x70, tw2);

        {
            ae_int32x2 s0, s1, s2, s3;
            ae_int32x2 d0, d1, d2, d3;

            x00 = AE_SRAI32(x00, 1);
            x10 = AE_SRAI32(x10, 1);
            x20 = AE_MULFP32X16X2S_H(x20, tw1);
            x40 = AE_MULFP32X16X2S_H(x40, tw1);
            x60 = AE_MULFP32X16X2S_H(x60, tw1);

            AE_ADDANDSUB32S(s0, d0, x00, x10);
            AE_ADDANDSUB32S(s1, d1, x20, x30);
            AE_ADDANDSUB32S(s2, d2, x40, x50);
            AE_ADDANDSUB32S(s3, d3, x60, x70);

            x00 = s0;        x40 = d0;
            x10 = s1;        x50 = d1;
            x20 = s2;        x60 = d2;
            x30 = s3;        x70 = d3;
        }

        DFT4X1RNG(x01, x21, x41, x61);
        DFT4X1RNG(x11, x31, x51, x71);
        x31 = AE_MULFC32X16RAS_L(x31, tw1);
        x51 = AE_MULFC32X16RAS_H(x51, tw2);
        x71 = AE_MULFC32X16RAS_L(x71, tw2);

        {
            ae_int32x2 s0, s1, s2, s3;
            ae_int32x2 d0, d1, d2, d3;

            x01 = AE_SRAI32(x01, 1);
            x11 = AE_SRAI32(x11, 1);
            x21 = AE_MULFP32X16X2S_H(x21, tw1);
            x41 = AE_MULFP32X16X2S_H(x41, tw1);
            x61 = AE_MULFP32X16X2S_H(x61, tw1);

            AE_ADDANDSUB32S(s0, d0, x01, x11);
            AE_ADDANDSUB32S(s1, d1, x21, x31);
            AE_ADDANDSUB32S(s2, d2, x41, x51);
            AE_ADDANDSUB32S(s3, d3, x61, x71);

            x01 = s0;        x41 = d0;
            x11 = s1;        x51 = d1;
            x21 = s2;        x61 = d2;
            x31 = s3;        x71 = d3;
        }

        AE_S32X2X2_XP(x00, x01, castxcc(ae_int32x4, py0), stride*sizeof(ae_int32x2));
        AE_S32X2X2_XP(x10, x11, castxcc(ae_int32x4, py0), stride*sizeof(ae_int32x2));
        AE_S32X2X2_XP(x20, x21, castxcc(ae_int32x4, py0), stride*sizeof(ae_int32x2));
        AE_S32X2X2_XP(x30, x31, castxcc(ae_int32x4, py0), (2 - 3 * stride)*sizeof(ae_int32x2));
        AE_S32X2X2_XP(x40, x41, castxcc(ae_int32x4, py1), stride*sizeof(ae_int32x2));
        AE_S32X2X2_XP(x50, x51, castxcc(ae_int32x4, py1), stride*sizeof(ae_int32x2));
        AE_S32X2X2_XP(x60, x61, castxcc(ae_int32x4, py1), stride*sizeof(ae_int32x2));
        AE_S32X2X2_XP(x70, x71, castxcc(ae_int32x4, py1), (2 - 3 * stride)*sizeof(ae_int32x2));
    }

    return shift+1;
} /* fft_32x16_stage_last_scl2_DFT8() */
