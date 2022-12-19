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
    Complex-valued FFT stages with butterflies radix-2, radix-3, radix-5
    with dynamic data scaling: 16-bit data, 16-bit twiddle factors
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

#include "NatureDSP_types.h"
#include "common.h"
#include "fft_16x16_stages.h"


ALIGN(32) static const int16_t __dft5_tw[] =
{
    (int16_t)0x278E, (int16_t)0x278E, (int16_t)0x278E, (int16_t)0x278E,
    (int16_t)0x8644, (int16_t)0x79BC, (int16_t)0x8644, (int16_t)0x79BC,
    (int16_t)0x9872, (int16_t)0x9872, (int16_t)0x9872, (int16_t)0x9872,
    (int16_t)0xB4C3, (int16_t)0x4B3D, (int16_t)0xB4C3, (int16_t)0x4B3D
};

/* twiddles should be loaded from the table above
   requires setting bits 7,0 of SAR register to zero
*/
#define DFT5X2_NOSCALING(x0, x1, x2, x3, x4, w1, w2, w3, w4)\
{ \
    ae_int16x4 s1, s2, d1, d2;             \
    ae_int16x4 t0, t2;                     \
    ae_int16x4 y0, y1, y2, y3;             \
    /*  AE_SEL16_6745(); AE_SEL16_2301();    */ \
    ae_int16x4 dsel = AE_MOVINT16X4_FROMINT64(0x0602070304000501);\
    AE_ADDANDSUBRNG16RAS_S1(x1,x4);        \
    AE_ADDANDSUBRNG16RAS_S1(x2,x3);        \
    s1=x1;d1=x4; s2=x2; d2=x3;             \
    t0 = AE_MULFD16X16X4RAS(s1,s2, w1,w3); \
    t2 = AE_MULFD16X16X4RAS(s1,s2, w3,w1); \
    y0 = AE_ADD16S(x0, t0);                \
    y1 = AE_ADD16S(x0, t2);                \
    y2 = AE_MULFD16X16X4RAS(d1,d2, w2,w4); \
    y3 = AE_MULFD16X16X4RAS(d1,d2, w4,AE_NEG16S(w2)); \
    AE_DSEL16X4(y2, y3, y2, y3, dsel);      \
    x0 = AE_ADD16S(x0, AE_ADD16S(s1, s2));\
    AE_ADDANDSUBRNG16RAS_S1(y0, y2);      \
    AE_ADDANDSUBRNG16RAS_S1(y1, y3);      \
    x1=y0; x4=y2; x2=y1; x3=y3;           \
}



/*
*  Last stage of FFT 16x16, radix-5, dynamic scaling
*/
int fft_16x16_stage_last_scl2_DFT5(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    int i;
    const int shift = XT_MAX(0, 3 - *bexp);
    const int N5=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),429496730));//N/5
    const int stride = N5;
    const ae_int16x4 * restrict px0;
    const ae_int16x4 * restrict px1;
    const ae_int16x4 * restrict px2;
    const ae_int16x4 * restrict px3;
    const ae_int16x4 * restrict px4;
    const ae_int16x8 * restrict ptwd_dft;
    ae_int16x4 * restrict py0;
    ae_int16x4 x0, x1, x2, x3, x4;
    ae_int16x4 x01, x11, x21, x31, x41;
    ae_int16x4 w1, w2, w3, w4;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT((stride & 1) == 0);

    ae_int16x4 scale=AE_SLAA16S(AE_MOVDA16(1),15-shift);

    /* Load twiddles for DFT5 */
    ptwd_dft = (const ae_int16x8 *)__dft5_tw;
    AE_L16X4X2_IP(w1, w2,ptwd_dft, sizeof(ae_int16x8));
    AE_L16X4X2_IP(w3, w4,ptwd_dft, sizeof(ae_int16x8));

    px0 = (ae_int16x4 *)x;
    px1 = px0 + stride / 2;
    px2 = px1 + stride / 2;
    px3 = px2 + stride / 2;
    px4 = px3 + stride / 2;
    py0 = (ae_int16x4 *)y;
    WUR_AE_SAR(0);
    __Pragma("loop_count min=1");
    for (i = 0; i <(stride >> 2); i++)
    { /* 14 cycles per pipeline stage in steady state with unroll=1 */

        AE_L16X4X2_IP(x0, x01, castxcc(ae_int16x8, px0), 2 * sizeof(ae_int16x4));
        AE_L16X4X2_IP(x1, x11, castxcc(ae_int16x8, px1), 2 * sizeof(ae_int16x4));
        AE_L16X4X2_IP(x2, x21, castxcc(ae_int16x8, px2), 2 * sizeof(ae_int16x4));
        AE_L16X4X2_IP(x3, x31, castxcc(ae_int16x8, px3), 2 * sizeof(ae_int16x4));
        AE_L16X4X2_IP(x4, x41, castxcc(ae_int16x8, px4), 2 * sizeof(ae_int16x4));
        /******************* DFT5 *********************/
        x0 = AE_SRAA16RS(x0, shift);
        x1=AE_MULFP16X4RS(x1,scale);
        x2=AE_MULFP16X4RS(x2,scale);
        x3=AE_MULFP16X4RS(x3,scale);
        x4=AE_MULFP16X4RS(x4,scale);
        DFT5X2_NOSCALING(x0, x1, x2, x3, x4, w1, w2, w3, w4);

        /******************* DFT5 *********************/
        x01 = AE_SRAA16RS(x01, shift);
        x11=AE_MULFP16X4RS(x11,scale);
        x21=AE_MULFP16X4RS(x21,scale);
        x31=AE_MULFP16X4RS(x31,scale);
        x41=AE_MULFP16X4RS(x41,scale);
        DFT5X2_NOSCALING(x01, x11, x21, x31, x41, w1, w2, w3, w4);

        AE_S16X4X2_XP(x0, x01, castxcc(ae_int16x8, py0), stride*sizeof(complex_fract16));
        AE_S16X4X2_XP(x1, x11, castxcc(ae_int16x8, py0), stride*sizeof(complex_fract16));
        AE_S16X4X2_XP(x2, x21, castxcc(ae_int16x8, py0), stride*sizeof(complex_fract16));
        AE_S16X4X2_XP(x3, x31, castxcc(ae_int16x8, py0), stride*sizeof(complex_fract16));
        AE_S16X4X2_XP(x4, x41, castxcc(ae_int16x8, py0), (4 - 4 * stride)*sizeof(complex_fract16));
    }
    *v *= 5;
    return shift;
} /* fft_16x16_stage_last_scl2_DFT5() */
