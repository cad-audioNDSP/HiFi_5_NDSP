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

ALIGN(32) static const int16_t __fft16_tw2_[] =
{
    (int16_t)0x7642, (int16_t)0xcf04, (int16_t)0x7642, (int16_t)0xcf04,
    (int16_t)0x5a82, (int16_t)0xa57e, (int16_t)0x5a82, (int16_t)0xa57e,
    (int16_t)0x30fc, (int16_t)0x89be, (int16_t)0x30fc, (int16_t)0x89be,
    (int16_t)0x5a82, (int16_t)0xa57e, (int16_t)0x5a82, (int16_t)0xa57e,
    (int16_t)0x0000, (int16_t)0x8000, (int16_t)0x0000, (int16_t)0x8000,
    (int16_t)0xa57e, (int16_t)0xa57e, (int16_t)0xa57e, (int16_t)0xa57e,
    (int16_t)0x30fc, (int16_t)0x89be, (int16_t)0x30fc, (int16_t)0x89be,
    (int16_t)0xa57e, (int16_t)0xa57e, (int16_t)0xa57e, (int16_t)0xa57e,
    (int16_t)0x89be, (int16_t)0x30fc, (int16_t)0x89be, (int16_t)0x30fc
};

int fft_16x16_stage_last_scl3_DFT16(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    ae_int16x4 x0, x1, x2, x3;
    ae_int16x4 x4, x5, x6, x7;
    ae_int16x4 x8, x9, x10, x11;
    ae_int16x4 x12, x13, x14, x15;
    ae_int16x4 t5, t9, t13;
    ae_int16x4 t6, t10, t14;
    ae_int16x4 t7, t11, t15;

    const ae_int16x4 * restrict px0;
    const ae_int16x4 * restrict px1;
    const ae_int16x4 * restrict px2;
    const ae_int16x4 * restrict px3;
    ae_int16x4 * restrict py0;

    const int stride = (N >> 4);
    int j, shift;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);

    shift = 2;

    px0 = (const ae_int16x4 *)x;
    px1 = (const ae_int16x4 *)(stride*sizeof(complex_fract16) + (uintptr_t)px0);
    px2 = (const ae_int16x4 *)(stride*sizeof(complex_fract16) + (uintptr_t)px1);
    px3 = (const ae_int16x4 *)(stride*sizeof(complex_fract16) + (uintptr_t)px2);
    py0 = (ae_int16x4 *)y;

    AE_L16X4X2_I(t5 ,t9 ,(const ae_int16x8*)__fft16_tw2_,0*sizeof(ae_int16x8));
    AE_L16X4X2_I(t13,t6 ,(const ae_int16x8*)__fft16_tw2_,1*sizeof(ae_int16x8));
    AE_L16X4X2_I(t10,t14,(const ae_int16x8*)__fft16_tw2_,2*sizeof(ae_int16x8));
    AE_L16X4X2_I(t7 ,t11,(const ae_int16x8*)__fft16_tw2_,3*sizeof(ae_int16x8));
    t15=AE_L16X4_X((const ae_int16x4*)__fft16_tw2_,4*sizeof(ae_int16x8));

    SetDFT4_Scaling(shift);

    __Pragma("loop_count min=1");
    for (j = 0; j < (stride>>1); j++)
    {
        /* 22 cycles with unroll=1 */
        AE_L16X4_XP(x0,  px0, 4*stride*sizeof(complex_fract16));
        AE_L16X4_XP(x1,  px1, 4*stride*sizeof(complex_fract16));
        AE_L16X4_XP(x2,  px2, 4*stride*sizeof(complex_fract16));
        AE_L16X4_XP(x3,  px3, 4*stride*sizeof(complex_fract16));
        AE_L16X4_XP(x4,  px0, 4*stride*sizeof(complex_fract16));
        AE_L16X4_XP(x5,  px1, 4*stride*sizeof(complex_fract16));
        AE_L16X4_XP(x6,  px2, 4*stride*sizeof(complex_fract16));
        AE_L16X4_XP(x7,  px3, 4*stride*sizeof(complex_fract16));
        AE_L16X4_XP(x8,  px0, 4*stride*sizeof(complex_fract16));
        AE_L16X4_XP(x9,  px1, 4*stride*sizeof(complex_fract16));
        AE_L16X4_XP(x10, px2, 4*stride*sizeof(complex_fract16));
        AE_L16X4_XP(x11, px3, 4*stride*sizeof(complex_fract16));
        AE_L16X4_XP(x12, px0, (2 - 12 * stride)*sizeof(complex_fract16));
        AE_L16X4_XP(x13, px1, (2 - 12 * stride)*sizeof(complex_fract16));
        AE_L16X4_XP(x14, px2, (2 - 12 * stride)*sizeof(complex_fract16));
        AE_L16X4_XP(x15, px3, (2 - 12 * stride)*sizeof(complex_fract16));

        DFT4XI2(x0, x4, x8, x12);
        DFT4XI2(x1, x5, x9, x13);
        DFT4XI2(x2, x6, x10, x14);
        DFT4XI2(x3, x7, x11, x15);

        x5 = AE_MULFC16RAS(x5, t5);
        x6 = AE_MULFC16RAS(x6, t6);
        x7 = AE_MULFC16RAS(x7, t7);
        x9 = AE_MULFC16RAS(x9, t6);
        x10 = AE_MULFC16RAS(x10, t10);
        x11 = AE_MULFC16RAS(x11, t11);
        x13 = AE_MULFC16RAS(x13, t7);
        x14 = AE_MULFC16RAS(x14, t11);
        x15 = AE_MULFC16RAS(x15, t15);

        DFT4XI2(x0, x1, x2, x3);
        DFT4XI2(x4, x5, x6, x7);
        DFT4XI2(x8, x9, x10, x11);
        DFT4XI2(x12, x13, x14, x15);

        AE_S16X4_XP(x0, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x4, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x8, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x12, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x1, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x5, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x9, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x13, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x2, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x6, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x10, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x14, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x3, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x7, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x11, py0, stride*sizeof(complex_fract16));
        AE_S16X4_XP(x15, py0, (2 - 15 * stride)*sizeof(complex_fract16));
    }
    return 2*shift;
} /* fft_16x16_stage_last_scl3_DFT16() */
