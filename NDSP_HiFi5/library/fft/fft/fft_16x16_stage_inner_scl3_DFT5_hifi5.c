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
    with static data scaling: 16-bit data, 16-bit twiddle factors
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
*  Intermediate stage of FFT/IFFT 16x16, radix-5, static scaling
*/
int fft_16x16_stage_inner_scl3_DFT5(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    /*  AE_SEL16_7632();
        AE_SEL16_5410();    */
    ALIGN(32) static const int16_t sel_tab[4] = {0x705,0x604,0x301,0x200};
    ae_int16x4 dsel=AE_L16X4_I((const ae_int16x4*)sel_tab,0);
    #define shift 3
    int i, j, _v;
    const int N5=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),429496730));//N/5
    const int stride = N5;
    const ae_int16x4 * restrict px0;
    ae_int16x4 * restrict py0;

    const ae_int16x4 * restrict ptwd;
    ae_int16x4 x00, x10, x20, x30, x40;
    ae_int16x4 x01, x11, x21, x31, x41;
    const ae_int16x8 * restrict ptwd_dft;
    ae_int16x4 w1, w2, w3, w4;
    ae_int16x4 tw1, tw2, tw3, tw4;
    
    _v = *v;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT((_v & 1) == 0);

    /* Load twiddles for DFT5 */
    ptwd_dft = (const ae_int16x8 *)__dft5_tw;
    AE_L16X4X2_IP(w1, w2,ptwd_dft, sizeof(ae_int16x8));
    AE_L16X4X2_IP(w3, w4,ptwd_dft, sizeof(ae_int16x8));

    WUR_AE_SAR(0);
    __Pragma("loop_count min=1");
    for (j = 0; j < (_v >> 1); j+=2)
    {
        ptwd = (const ae_int16x4 *)tw;
        px0 = (const ae_int16x4 *)x + j;
        py0 = (ae_int16x4 *)y + j;

        __Pragma("loop_count min=1");
        for (i = 0; i < (stride / _v); i++)
        {
            /* 18 cycles per pipeline stage in steady state with unroll=1 */

            AE_L16X4X2_XP(x00, x01, castxcc(ae_int16x8, px0), stride*sizeof(complex_fract16));
            AE_L16X4X2_XP(x10, x11, castxcc(ae_int16x8, px0), stride*sizeof(complex_fract16));
            AE_L16X4X2_XP(x20, x21, castxcc(ae_int16x8, px0), stride*sizeof(complex_fract16));
            AE_L16X4X2_XP(x30, x31, castxcc(ae_int16x8, px0), stride*sizeof(complex_fract16));
            AE_L16X4X2_XP(x40, x41, castxcc(ae_int16x8, px0), (4 - 4 * stride)*sizeof(complex_fract16));

            /******************* DFT5 *********************/
            x00 = AE_SRAI16R(x00, shift);
            x10 = AE_SRAI16R(x10, shift);
            x20 = AE_SRAI16R(x20, shift);
            x30 = AE_SRAI16R(x30, shift);
            x40 = AE_SRAI16R(x40, shift);

            DFT5X2_NOSCALING(x00, x10, x20, x30, x40, w1, w2, w3, w4);

            AE_L16X4X2_IP(tw1, tw3, castxcc(ae_int16x8, ptwd), sizeof(ae_int16x8));
            AE_DSEL16X4(tw1,tw2,tw1,tw1,dsel);
            AE_DSEL16X4(tw3,tw4,tw3,tw3,dsel);

            x10 = AE_MULFC16RAS(x10, tw1);
            x20 = AE_MULFC16RAS(x20, tw2);
            x30 = AE_MULFC16RAS(x30, tw3);
            x40 = AE_MULFC16RAS(x40, tw4);

            /******************* DFT5 *********************/
            x01 = AE_SRAI16R(x01, shift);
            x11 = AE_SRAI16R(x11, shift);
            x21 = AE_SRAI16R(x21, shift);
            x31 = AE_SRAI16R(x31, shift);
            x41 = AE_SRAI16R(x41, shift);

            DFT5X2_NOSCALING(x01, x11, x21, x31, x41, w1, w2, w3, w4);

            x11 = AE_MULFC16RAS(x11, tw1);
            x21 = AE_MULFC16RAS(x21, tw2);
            x31 = AE_MULFC16RAS(x31, tw3);
            x41 = AE_MULFC16RAS(x41, tw4);

            AE_S16X4X2_XP(x00, x01, castxcc(ae_int16x8,py0), _v*sizeof(complex_fract16));
            AE_S16X4X2_XP(x10, x11, castxcc(ae_int16x8,py0), _v*sizeof(complex_fract16));
            AE_S16X4X2_XP(x20, x21, castxcc(ae_int16x8,py0), _v*sizeof(complex_fract16));
            AE_S16X4X2_XP(x30, x31, castxcc(ae_int16x8,py0), _v*sizeof(complex_fract16));
            AE_S16X4X2_XP(x40, x41, castxcc(ae_int16x8,py0), _v*sizeof(complex_fract16));
        }
    }

    *v *= 5;
    return shift;
#undef shift
} /* fft_16x16_stage_inner_scl3_DFT5() */
