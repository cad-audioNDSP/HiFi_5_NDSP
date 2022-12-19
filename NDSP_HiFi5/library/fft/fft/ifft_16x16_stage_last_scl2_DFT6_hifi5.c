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
*  Last stage of FFT 16x16, radix-6, dynamic scaling
*/
int ifft_16x16_stage_last_scl2_DFT6(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{

    int i;
    const int N6=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),357913941));//N/6
    const int stride = N6;
    const ae_int16x4 * restrict px0 = (const ae_int16x4 *)x;
    ae_int16x4 * restrict py0 = (ae_int16x4 *)y;
    ae_int16x4 x00, x10, x20, x30, x40, x50,
        c = AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32(0x00006EDA));

    ae_int16x4 x01, x11, x21, x31, x41, x51;
    // const int shift = XT_MAX(0, 3 - *bexp);

    int s2 = 2,
        s1 = 1;
    int shift = s1 + s2;
    WUR_AE_SAR(s2 * 0x102 + s1 * 0x81);

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT((stride & 1) == 0);


    __Pragma("loop_count min=1");
    for (i = 0; i <(stride >> 2); i++)
    { /* 23 cycles per pipeline stage in steady state with unroll=1 */

        AE_L16X4X2_XP(x00, x01, castxcc(ae_int16x8, px0), stride * sizeof(complex_fract16));
        AE_L16X4X2_XP(x10, x11, castxcc(ae_int16x8, px0), stride * sizeof(complex_fract16));
        AE_L16X4X2_XP(x20, x21, castxcc(ae_int16x8, px0), stride * sizeof(complex_fract16));
        AE_L16X4X2_XP(x30, x31, castxcc(ae_int16x8, px0), stride * sizeof(complex_fract16));
        AE_L16X4X2_XP(x40, x41, castxcc(ae_int16x8, px0), stride * sizeof(complex_fract16));
        AE_L16X4X2_XP(x50, x51, castxcc(ae_int16x8, px0), (4 - 5 * stride) * sizeof(complex_fract16));

        {
            /* DFT6. Prime-factor FFT algorithm */
            ae_int16x4 _s0, _s1, _s2, _d0, _d1, _d2;

            _s0 = x00;    _d0 = x30;    _s1 = x40;
            _s2 = x20;    _d2 = x50;    _d1 = x10;
            /* Stage1: DFT2 */
            AE_ADDANDSUBRNG16RAS_S1(_s0, _d0);
            AE_ADDANDSUBRNG16RAS_S1(_s1, _d1);
            AE_ADDANDSUBRNG16RAS_S1(_s2, _d2);
            /* Stage2: DFT3 */
            DFT3XI2_S2(_s0, _s1, _s2, c, s2);
            DFT3XI2_S2(_d0, _d1, _d2, c, s2);

            x00 = _s0;    x10 = _d2;    x20 = _s1;
            x40 = _s2;    x50 = _d1;    x30 = _d0;
        }
        {
            /* DFT6. Prime-factor FFT algorithm */
            ae_int16x4 _s0, _s1, _s2, _d0, _d1, _d2;

            _s0 = x01;    _d0 = x31;    _s1 = x41;
            _s2 = x21;    _d2 = x51;    _d1 = x11;
            /* Stage1: DFT2 */
            AE_ADDANDSUBRNG16RAS_S1(_s0, _d0);
            AE_ADDANDSUBRNG16RAS_S1(_s1, _d1);
            AE_ADDANDSUBRNG16RAS_S1(_s2, _d2);
            /* Stage2: DFT3 */
            DFT3XI2_S2(_s0, _s1, _s2, c, s2);
            DFT3XI2_S2(_d0, _d1, _d2, c, s2);

            x01 = _s0;    x11 = _d2;    x21 = _s1;
            x41 = _s2;    x51 = _d1;    x31 = _d0;
        }

        AE_S32X2X2_XP(AE_MOVINT32X2_FROMINT16X4(x00), AE_MOVINT32X2_FROMINT16X4(x01), castxcc(ae_int32x4, py0), stride * sizeof(complex_fract16));
        AE_S32X2X2_XP(AE_MOVINT32X2_FROMINT16X4(x10), AE_MOVINT32X2_FROMINT16X4(x11), castxcc(ae_int32x4, py0), stride * sizeof(complex_fract16));
        AE_S32X2X2_XP(AE_MOVINT32X2_FROMINT16X4(x20), AE_MOVINT32X2_FROMINT16X4(x21), castxcc(ae_int32x4, py0), stride * sizeof(complex_fract16));
        AE_S32X2X2_XP(AE_MOVINT32X2_FROMINT16X4(x30), AE_MOVINT32X2_FROMINT16X4(x31), castxcc(ae_int32x4, py0), stride * sizeof(complex_fract16));
        AE_S32X2X2_XP(AE_MOVINT32X2_FROMINT16X4(x40), AE_MOVINT32X2_FROMINT16X4(x41), castxcc(ae_int32x4, py0), stride * sizeof(complex_fract16));
        AE_S32X2X2_XP(AE_MOVINT32X2_FROMINT16X4(x50), AE_MOVINT32X2_FROMINT16X4(x51), castxcc(ae_int32x4, py0), (4 - 5 * stride) * sizeof(complex_fract16));
    }
    *v *= 6;
    return shift;
} /* ifft_16x16_stage_last_scl2_DFT6() */
