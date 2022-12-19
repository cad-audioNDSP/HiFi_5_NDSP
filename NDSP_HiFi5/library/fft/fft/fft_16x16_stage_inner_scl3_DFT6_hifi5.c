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

inline_  void __DFT3xI2_s2(ae_int16x4 *x0, ae_int16x4 *x1, ae_int16x4 *x2, const ae_int16x4 *c, int shift)
{
    ae_int16x4 s0, s1, d0;

    s0 = *x1;
    d0 = *x2;
    *x0 = AE_SRAA16RS(*x0, shift);
    AE_ADDANDSUBRNG16RAS_S2(s0, d0);

    s1 = AE_ADD16S(*x0, s0);
    s0 = AE_SRAI16(s0, 1);
    s0 = AE_SUB16S(*x0, s0);
    d0 = AE_MULFC16RAS(*c, d0);
    *x0 = s1;

    *x1 = AE_SUB16S(s0, d0);
    *x2 = AE_ADD16S(s0, d0);
}
#define DFT3XI2_S2(__x0, __x1, __x2, __c, __s) __DFT3xI2_s2(&__x0, &__x1, &__x2, (const ae_int16x4*)&__c, __s)

/*
*  Intermediate stage of FFT/IFFT 16x16, radix-6, static scaling
*/
int fft_16x16_stage_inner_scl3_DFT6(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    int i, j, _v;
    const int N6=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),357913941));//N/6
    int shift;
    const int stride = N6;
    const ae_int16x4 * restrict px0;
    ae_int16x4 * restrict py0;
    const ae_int16x4 * restrict ptwd;
    ae_int16x4 x0, x1, x2, x3, x4, x5;
    ae_int16x4 tw1, tw2, tw3, tw4, tw5, c = AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32(0x00006EDA));
    ae_int32x2 t1, t2, t3, t4, t5;
    _v = *v;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT((_v & 1) == 0);
    int s2 = 2,
        s1 = 1;
    shift = s1 + s2;

    WUR_AE_SAR(s2 * 0x102 + s1 * 0x81);
    
    if (_v&2)
    {
    __Pragma("loop_count min=1");
    for (j = 0; j < (_v >> 1); j++)
    {
        ptwd = (const ae_int16x4 *)tw;
        px0 = (const ae_int16x4 *)x + j;
        py0 = (ae_int16x4 *)y + j;
        __Pragma("loop_count min=1");
        for (i = 0; i < (stride / _v); i++)
        {
            /* 12 cycles per pipeline stage in steady state with unroll=1 */
            x1 = AE_L16X4_X(px0, 1 * stride*sizeof(complex_fract16));
            x2 = AE_L16X4_X(px0, 2 * stride*sizeof(complex_fract16));
            x3 = AE_L16X4_X(px0, 3 * stride*sizeof(complex_fract16));
            x4 = AE_L16X4_X(px0, 4 * stride*sizeof(complex_fract16));
            x5 = AE_L16X4_X(px0, 5 * stride*sizeof(complex_fract16));
            AE_L16X4_XP(x0, px0, _v*sizeof(complex_fract16));

            {
                /* DFT6. Prime-factor FFT algorithm */
                ae_int16x4 _s0, _s1, _s2, _d0, _d1, _d2;

                _s0 = x0;    _d0 = x3;    _s1 = x4;    _d1 = x1;
                _s2 = x2;    _d2 = x5;
                /* Stage1: DFT2 */
                AE_ADDANDSUBRNG16RAS_S1(_s0, _d0);
                AE_ADDANDSUBRNG16RAS_S1(_s1, _d1);
                AE_ADDANDSUBRNG16RAS_S1(_s2, _d2);
                /* Stage2: DFT3 */
                DFT3XI2_S2(_s0, _s1, _s2, c, s2);
                DFT3XI2_S2(_d0, _d1, _d2, c, s2);

                x0 = _s0;    x1 = _d2;    x2 = _s1;    x3 = _d0;
                x4 = _s2;    x5 = _d1;
            }

            AE_L32_XP(t1, castxcc(ae_int32, ptwd), sizeof(complex_fract16));
            AE_L32_XP(t2, castxcc(ae_int32, ptwd), sizeof(complex_fract16));
            AE_L32_XP(t3, castxcc(ae_int32, ptwd), sizeof(complex_fract16));
            AE_L32_XP(t4, castxcc(ae_int32, ptwd), sizeof(complex_fract16));
            AE_L32_XP(t5, castxcc(ae_int32, ptwd), sizeof(complex_fract16));

            tw1 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t1));
            tw2 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t2));
            tw3 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t3));
            tw4 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t4));
            tw5 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t5));

            x1 = AE_MULFC16RAS(x1, tw1);
            x2 = AE_MULFC16RAS(x2, tw2);
            x3 = AE_MULFC16RAS(x3, tw3);
            x4 = AE_MULFC16RAS(x4, tw4);
            x5 = AE_MULFC16RAS(x5, tw5);

            AE_S16X4_XP(x0, py0, _v * sizeof(complex_fract16));
            AE_S16X4_XP(x1, py0, _v * sizeof(complex_fract16));
            AE_S16X4_XP(x2, py0, _v * sizeof(complex_fract16));
            AE_S16X4_XP(x3, py0, _v * sizeof(complex_fract16));
            AE_S16X4_XP(x4, py0, _v * sizeof(complex_fract16));
            AE_S16X4_XP(x5, py0, _v * sizeof(complex_fract16));
        }
    }
    }
    else
    {
        NASSERT((_v&3) == 0);
        __Pragma("loop_count min=1");
        for (j = 0; j < (_v >> 1); j+=2)
        {
            ae_int16x4 x0_, x1_, x2_, x3_, x4_, x5_;
            ptwd = (const ae_int16x4 *)tw;
            px0 = (const ae_int16x4 *)x + j;
            py0 = (ae_int16x4 *)y + j;
            __Pragma("loop_count min=1");
            for (i = 0; i < (stride / _v); i++)
            {
                /* 21 cycles per pipeline stage in steady state with unroll=1 */
                AE_L16X4X2_X(x1, x1_, (const ae_int16x8*)px0, 1 * stride*sizeof(complex_fract16));
                AE_L16X4X2_X(x2, x2_, (const ae_int16x8*)px0, 2 * stride*sizeof(complex_fract16));
                AE_L16X4X2_X(x3, x3_, (const ae_int16x8*)px0, 3 * stride*sizeof(complex_fract16));
                AE_L16X4X2_X(x4, x4_, (const ae_int16x8*)px0, 4 * stride*sizeof(complex_fract16));
                AE_L16X4X2_X(x5, x5_, (const ae_int16x8*)px0, 5 * stride*sizeof(complex_fract16));
                AE_L16X4X2_XP(x0, x0_, castxcc(ae_int16x8, px0), _v*sizeof(complex_fract16));

                {
                    /* DFT6. Prime-factor FFT algorithm */
                    ae_int16x4 _s0, _s1, _s2, _d0, _d1, _d2;

                    _s0 = x0;    _d0 = x3;    _s1 = x4;    _d1 = x1;
                    _s2 = x2;    _d2 = x5;
                    /* Stage1: DFT2 */
                    AE_ADDANDSUBRNG16RAS_S1(_s0, _d0);
                    AE_ADDANDSUBRNG16RAS_S1(_s1, _d1);
                    AE_ADDANDSUBRNG16RAS_S1(_s2, _d2);
                    /* Stage2: DFT3 */
                    DFT3XI2_S2(_s0, _s1, _s2, c, s2);
                    DFT3XI2_S2(_d0, _d1, _d2, c, s2);

                    x0 = _s0;    x1 = _d2;    x2 = _s1;    x3 = _d0;
                    x4 = _s2;    x5 = _d1;
                }

                {
                    /* DFT6. Prime-factor FFT algorithm */
                    ae_int16x4 _s0, _s1, _s2, _d0, _d1, _d2;

                    _s0 = x0_;    _d0 = x3_;    _s1 = x4_;    _d1 = x1_;
                    _s2 = x2_;    _d2 = x5_;
                    /* Stage1: DFT2 */
                    AE_ADDANDSUBRNG16RAS_S1(_s0, _d0);
                    AE_ADDANDSUBRNG16RAS_S1(_s1, _d1);
                    AE_ADDANDSUBRNG16RAS_S1(_s2, _d2);
                    /* Stage2: DFT3 */
                    DFT3XI2_S2(_s0, _s1, _s2, c, s2);
                    DFT3XI2_S2(_d0, _d1, _d2, c, s2);

                    x0_ = _s0;    x1_ = _d2;    x2_ = _s1;    x3_ = _d0;
                    x4_ = _s2;    x5_ = _d1;
                }

                AE_L32_XP(t1, castxcc(ae_int32, ptwd), sizeof(complex_fract16));
                AE_L32_XP(t2, castxcc(ae_int32, ptwd), sizeof(complex_fract16));
                AE_L32_XP(t3, castxcc(ae_int32, ptwd), sizeof(complex_fract16));
                AE_L32_XP(t4, castxcc(ae_int32, ptwd), sizeof(complex_fract16));
                AE_L32_XP(t5, castxcc(ae_int32, ptwd), sizeof(complex_fract16));

                tw1 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t1));
                tw2 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t2));
                tw3 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t3));
                tw4 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t4));
                tw5 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t5));

                x1 = AE_MULFC16RAS(x1, tw1);
                x2 = AE_MULFC16RAS(x2, tw2);
                x3 = AE_MULFC16RAS(x3, tw3);
                x4 = AE_MULFC16RAS(x4, tw4);
                x5 = AE_MULFC16RAS(x5, tw5);

                x1_ = AE_MULFC16RAS(x1_, tw1);
                x2_ = AE_MULFC16RAS(x2_, tw2);
                x3_ = AE_MULFC16RAS(x3_, tw3);
                x4_ = AE_MULFC16RAS(x4_, tw4);
                x5_ = AE_MULFC16RAS(x5_, tw5);

                AE_S16X4X2RNG_XP(x0, x0_, castxcc(ae_int16x8, py0), _v * sizeof(complex_fract16));
                AE_S16X4X2RNG_XP(x1, x1_, castxcc(ae_int16x8, py0), _v * sizeof(complex_fract16));
                AE_S16X4X2RNG_XP(x2, x2_, castxcc(ae_int16x8, py0), _v * sizeof(complex_fract16));
                AE_S16X4X2RNG_XP(x3, x3_, castxcc(ae_int16x8, py0), _v * sizeof(complex_fract16));
                AE_S16X4X2RNG_XP(x4, x4_, castxcc(ae_int16x8, py0), _v * sizeof(complex_fract16));
                AE_S16X4X2RNG_XP(x5, x5_, castxcc(ae_int16x8, py0), _v * sizeof(complex_fract16));
            }
        }
        }

    *v *= 6;
    return shift;
} /* fft_16x16_stage_inner_scl3_DFT6() */
