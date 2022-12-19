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
 

/* radix-4 butterfly with normalization applied twice */
#define DFT4X1RNG2(x0, x1, x2, x3)\
{   \
    ae_int32x2 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG32(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32(s1, d1, x1, x3); \
    d1 = AE_MUL32JS(d1);               \
    AE_ADDANDSUBRNG32(x0, x2, s0, s1); \
    AE_ADDANDSUBRNG32(x3, x1, d0, d1); \
}

/*
 *  32x16 IFFT last stage Radix 4, scalingOption = 3
 */
int ifft_32x16_stage_last_scl3_DFT4(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    const ae_int32x2 * restrict px0;
    const ae_int32x2 * restrict px1;
    const ae_int32x2 * restrict px2;
    const ae_int32x2 * restrict px3;
    ae_int32x2 * restrict py0;
    ae_int32x2 * restrict py1;
    ae_int32x2 * restrict py2;
    ae_int32x2 * restrict py3;
    const int stride = (N >> 2);
    int j, shift;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);

    shift = 2;
    WUR_AE_SAR(shift >> 1);

    px0 = (const ae_int32x2 *)x;
    px1 = (const ae_int32x2 *)px0 + stride;
    px2 = (const ae_int32x2 *)px1 + stride;
    px3 = (const ae_int32x2 *)px2 + stride;
    py0 = (ae_int32x2 *)y;
    py1 = (ae_int32x2 *)py0 + stride;
    py2 = (ae_int32x2 *)py1 + stride;
    py3 = (ae_int32x2 *)py2 + stride;

    NASSERT((stride & 1) == 0);
    __Pragma("loop_count min=2");
    for (j = 0; j < (stride >> 1); j++)
    {
        /* 4 cycles per pipeline stage in steady state with unroll=1 */
        ae_int32x2 x0, x1, x2, x3;
        ae_int32x2 x4, x5, x6, x7;

        AE_L32X2X2_IP(x0, x4, castxcc(ae_int32x4, px0), 2 * sizeof(ae_int32x2));
        AE_L32X2X2_IP(x1, x5, castxcc(ae_int32x4, px1), 2 * sizeof(ae_int32x2));
        AE_L32X2X2_IP(x2, x6, castxcc(ae_int32x4, px2), 2 * sizeof(ae_int32x2));
        AE_L32X2X2_IP(x3, x7, castxcc(ae_int32x4, px3), 2 * sizeof(ae_int32x2));

        DFT4X1RNG2(x0, x1, x2, x3);
        DFT4X1RNG2(x4, x5, x6, x7);

        AE_S64X2_IP(AE_MOVINT64_FROMINT32X2(x0), AE_MOVINT64_FROMINT32X2(x4), castxcc(ae_int64x2, py0), 2 * sizeof(ae_int32x2));
        AE_S64X2_IP(AE_MOVINT64_FROMINT32X2(x1), AE_MOVINT64_FROMINT32X2(x5), castxcc(ae_int64x2, py1), 2 * sizeof(ae_int32x2));
        AE_S64X2_IP(AE_MOVINT64_FROMINT32X2(x2), AE_MOVINT64_FROMINT32X2(x6), castxcc(ae_int64x2, py2), 2 * sizeof(ae_int32x2));
        AE_S64X2_IP(AE_MOVINT64_FROMINT32X2(x3), AE_MOVINT64_FROMINT32X2(x7), castxcc(ae_int64x2, py3), 2 * sizeof(ae_int32x2));
    }
    return shift;
}