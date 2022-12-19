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
    with static data scaling: 16-bit data, 16-bit twiddle factors
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/
#include "NatureDSP_types.h"
#include "common.h"
#include "fft_16x16_stages.h"

/*
    Set scaling for DFT4
    Range of the 'scale' is 0...3 
*/
#define SetDFT4_Scaling(scale)          \
{                                       \
    int sar = 0;                        \
    if (scale == 3)        sar = 0x285; \
    else if (scale == 2)   sar = 0x183; \
    else if (scale == 1)   sar = 0x102; \
    else sar = 0;                       \
    WUR_AE_SAR(sar);                    \
}

/*  16-bit radix-4 butterfly with scaling. 
    Call SetDFT4_Scaling() before */
#define DFT4XI2(_x0, _x1, _x2, _x3)               \
{                                                 \
    ae_int16x4 s0, s1, d0, d1;                    \
    s0 = _x0;    s1 = _x1;                        \
    d0 = _x2;    d1 = _x3;                        \
    AE_ADDANDSUBRNG16RAS_S1(s0, d0);              \
    AE_ADDANDSUBRNG16RAS_S1(s1, d1);              \
    d1 = AE_MUL16JS(d1);                          \
    AE_ADDANDSUBRNG16RAS_S2(s0, s1);              \
    AE_ADDANDSUBRNG16RAS_S2(d0, d1);              \
    _x0 = s0;    _x2 = s1;                        \
    _x3 = d0;    _x1 = d1;                        \
}

/*
 *  Intermediate stage of FFT/IFFT 16x16, radix-4, static scaling
 */
int fft_16x16_stage_inner_scl3_DFT4x2_v4(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    #define shift 2
    int i;
    const ae_int16x4 * restrict px0;
    const ae_int16x4 * restrict px1; 
    const ae_int16x4 * restrict px2; 
    const ae_int16x4 * restrict px3; 
        ae_int16x4 * restrict py0;
    const ae_int32 * restrict ptwd32;

    const int R = 4; /* stage radix */
    const int stride = (N >> 2);
    ae_int16x4 x00, x10, x20, x30;
    ae_int16x4 x01, x11, x21, x31;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(v[0] == 4);
    NASSERT(N > 16);

    SetDFT4_Scaling(shift); 

    ptwd32 = (const ae_int32*)tw;
    px0 = (ae_int16x4 *)x ;
    px1 = px0 + stride / 2;
    px2 = px1 + stride / 2;
    px3 = px2 + stride / 2;
    py0 = (ae_int16x4 *)y ;
    __Pragma("loop_count min=1");
    for (i = 0; i < (N >> 4); i++)
    {
        /* 8 cycles per pipeline stage in steady state with unroll=1 */
        ae_int32x2 t1, t2, t3;
        ae_int16x4 tw1, tw2, tw3;

        AE_L16X4X2_XP(x00, x01, castxcc(ae_int16x8, px0), 4 * 2 * sizeof(int16_t));
        AE_L16X4X2_XP(x10, x11, castxcc(ae_int16x8, px1), 4 * 2 * sizeof(int16_t));
        AE_L16X4X2_XP(x20, x21, castxcc(ae_int16x8, px2), 4 * 2 * sizeof(int16_t));
        AE_L16X4X2_XP(x30, x31, castxcc(ae_int16x8, px3), 4 * 2 * sizeof(int16_t));

        AE_L32_XP(t1, ptwd32, sizeof(complex_fract16));
        AE_L32_XP(t2, ptwd32, sizeof(complex_fract16));
        AE_L32_XP(t3, ptwd32, sizeof(complex_fract16));

        tw1 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t1));
        tw2 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t2));
        tw3 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t3));

        DFT4XI2(x00, x10, x20, x30);
        DFT4XI2(x01, x11, x21, x31);

        x10 = AE_MULFC16RAS(x10, tw1);
        x20 = AE_MULFC16RAS(x20, tw2);
        x30 = AE_MULFC16RAS(x30, tw3);
        x11 = AE_MULFC16RAS(x11, tw1);
        x21 = AE_MULFC16RAS(x21, tw2);
        x31 = AE_MULFC16RAS(x31, tw3);

        AE_S16X4X2_XP(x00, x01, castxcc(ae_int16x8, py0), 4 * 2 * sizeof(int16_t));
        AE_S16X4X2_XP(x10, x11, castxcc(ae_int16x8, py0), 4 * 2 * sizeof(int16_t));
        AE_S16X4X2_XP(x20, x21, castxcc(ae_int16x8, py0), 4 * 2 * sizeof(int16_t));
        AE_S16X4X2_XP(x30, x31, castxcc(ae_int16x8, py0), 4 * 2 * sizeof(int16_t));
    } 

    *v *= R;
    return shift;
} /* fft_16x16_stage_inner_scl3_DFT4x2_v4 */
#undef shift
