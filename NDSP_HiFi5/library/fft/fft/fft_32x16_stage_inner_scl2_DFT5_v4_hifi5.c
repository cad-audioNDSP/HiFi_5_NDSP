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
    with dynamic data scaling: 32-bit data, 16-bit twiddle factors
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/
#include "NatureDSP_types.h"
#include "common.h"
#include "fft_32x16_stages.h"
 
inline_ void _cmult32x32(ae_int32x2 *result, ae_int32x2 *x, ae_int32x2 *y)
{
    ae_f32x2 z;
    z = AE_MULFC32RAS(AE_MOVF32X2_FROMINT32X2(*x), AE_MOVF32X2_FROMINT32X2(*y));
    *result = AE_MOVINT32X2_FROMF32X2(z);
}
/*
    DFT5 algorithm with scaling:
    x - input complex vector
    y - output complex vector
    y = fft(x)
    w1 =  exp(-1j*2*pi/5);
    w2 =  exp(-1j*2*pi*2/5);

    y = zeros(5,1);
    s1 = (x1+x4);
    s2 = (x2+x3);
    d1 = (x1-x4);
    d2 = (x2-x3);

    y(1) = x0 + s1 + s2;
    y(2) = x0 + (s1*real(w1) + s2*real(w2)) + 1j*(d1*imag(w1) + d2*imag(w2));
    y(5) = x0 + (s1*real(w1) + s2*real(w2)) - 1j*(d1*imag(w1) + d2*imag(w2));
    y(3) = x0 + (s1*real(w2) + s2*real(w1)) + 1j*(d1*imag(w2)  - d2*imag(w1));
    y(4) = x0 + (s1*real(w2) + s2*real(w1)) - 1j*(d1*imag(w2)  - d2*imag(w1));
*/
#define DFT5X1S(x0, x1, x2, x3, x4, shift)                       \
{                                                                \
    ae_int32x2 s1, s2, d1, d2;                                   \
    ae_int32x2 y0, y1, y2;                                       \
    ae_int32x2 t0, t1, t2, t3;                                   \
    ae_int32x2 real_w1, jimag_w1, real_w2, jimag_w2;             \
    real_w1 = AE_MOVDA32X2(0x278DDE6E, 0x0);                     \
    jimag_w1 = AE_MOVDA32X2(0x0, 0x8643C7B3);                    \
    real_w2 = AE_MOVDA32X2(0x98722192, 0x0);                     \
    jimag_w2 = AE_MOVDA32X2(0x0, 0xB4C373EE);                    \
    x0 = AE_SRAA32RS(x0, shift);                                  \
    AE_ADDANDSUBRNG32(s1, d1, x1, x4);                           \
    AE_ADDANDSUBRNG32(s2, d2, x2, x3);                           \
    _cmult32x32(&t0, &s1, &real_w1);                             \
    _cmult32x32(&t1, &s2, &real_w2);                             \
    _cmult32x32(&t2, &s1, &real_w2);                             \
    _cmult32x32(&t3, &s2, &real_w1);                             \
    y0 = AE_ADD32S(x0, AE_ADD32S(s1, s2));                       \
    y1 = AE_ADD32S(x0, AE_ADD32S(t0, t1));                       \
    y2 = AE_ADD32S(x0, AE_ADD32S(t2, t3));                       \
    _cmult32x32(&t0, &d1, &jimag_w1);                            \
    _cmult32x32(&t1, &d2, &jimag_w2);                            \
    _cmult32x32(&t2, &d2, &jimag_w1);                            \
    _cmult32x32(&t3, &d1, &jimag_w2);                            \
    t0 = AE_ADD32S(t0, t1);                                      \
    t1 = AE_SUB32S(t3, t2);                                      \
    x0 = y0;                                                     \
    AE_ADDANDSUB32S(x1, x4, y1, t0);                             \
    AE_ADDANDSUB32S(x2, x3, y2, t1);                             \
}

/*
 *  32x16 FFT/IFFT intermediate stage Radix 5, scalingOption=2, v=4
 */
int fft_32x16_stage_inner_scl2_DFT5_v4(const int16_t *tw, const int32_t *x, int32_t *y, int N, int *v, int tw_step, int *bexp)
{
    const int shift = XT_MAX(0, 3 - *bexp);
    int i;
    const int R = 5; // stage radix
    const int stride = N / R;
    ae_int32x2 * restrict px0 = (ae_int32x2 *)x + 0;
    ae_int32x2 * restrict py0 = (ae_int32x2 *)y + 0;
    ae_int32x2 * restrict px1 = (ae_int32x2 *)x + 2;
    ae_int32x2 * restrict py1 = (ae_int32x2 *)y + 2;
    const ae_int16x4 * restrict ptwd = (const ae_int16x4 *)tw;

    ae_int32x2 x00, x10, x20, x30, x40;
    ae_int32x2 x01, x11, x21, x31, x41;
    ae_int16x4 tw12, tw34;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(v[0] == 4);
    NASSERT(tw_step == 1);

    /* Set scaling for AE_ADDANDSUBRNG32 */
    WUR_AE_SAR(shift);

    __Pragma("loop_count min=1");
    for (i = 0; i < (stride >> 2); i++)
    {
        AE_L16X4X2_IP(tw12, tw34, castxcc(ae_int16x8, ptwd), 2 * sizeof(ae_int16x4));

        AE_L32X2X2_X(x10, x11, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride);
        AE_L32X2X2_X(x20, x21, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride * 2);
        AE_L32X2X2_X(x30, x31, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride * 3);
        AE_L32X2X2_X(x40, x41, (ae_int32x4*)px0, sizeof(ae_int32x2)* stride * 4);
        AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px0), 4*sizeof(ae_int32x2));

        DFT5X1S(x00, x10, x20, x30, x40, shift);
        DFT5X1S(x01, x11, x21, x31, x41, shift);

        x10 = AE_MULFC32X16RAS_H(x10, tw12);
        x20 = AE_MULFC32X16RAS_L(x20, tw12);
        x30 = AE_MULFC32X16RAS_H(x30, tw34);
        x40 = AE_MULFC32X16RAS_L(x40, tw34);
        x11 = AE_MULFC32X16RAS_H(x11, tw12);
        x21 = AE_MULFC32X16RAS_L(x21, tw12);
        x31 = AE_MULFC32X16RAS_H(x31, tw34);
        x41 = AE_MULFC32X16RAS_L(x41, tw34);

        AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x40, x41, castxcc(ae_int32x4, py0), 4 * sizeof(ae_int32x2));

        /////////

        //AE_L16X4X2_IP(tw12, tw34, castxcc(ae_int16x8, ptwd), 2 * sizeof(ae_int16x4));

        AE_L32X2X2_X(x10, x11, (ae_int32x4*)px1, sizeof(ae_int32x2)* stride);
        AE_L32X2X2_X(x20, x21, (ae_int32x4*)px1, sizeof(ae_int32x2)* stride * 2);
        AE_L32X2X2_X(x30, x31, (ae_int32x4*)px1, sizeof(ae_int32x2)* stride * 3);
        AE_L32X2X2_X(x40, x41, (ae_int32x4*)px1, sizeof(ae_int32x2)* stride * 4);
        AE_L32X2X2_XP(x00, x01, castxcc(ae_int32x4, px1), 4*sizeof(ae_int32x2));

        DFT5X1S(x00, x10, x20, x30, x40, shift);
        DFT5X1S(x01, x11, x21, x31, x41, shift);

        x10 = AE_MULFC32X16RAS_H(x10, tw12);
        x20 = AE_MULFC32X16RAS_L(x20, tw12);
        x30 = AE_MULFC32X16RAS_H(x30, tw34);
        x40 = AE_MULFC32X16RAS_L(x40, tw34);
        x11 = AE_MULFC32X16RAS_H(x11, tw12);
        x21 = AE_MULFC32X16RAS_L(x21, tw12);
        x31 = AE_MULFC32X16RAS_H(x31, tw34);
        x41 = AE_MULFC32X16RAS_L(x41, tw34);

        AE_S32X2X2RNG_XP(x00, x01, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x10, x11, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x20, x21, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x30, x31, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
        AE_S32X2X2RNG_XP(x40, x41, castxcc(ae_int32x4, py1), 4 * sizeof(ae_int32x2));
    }


  AE_CALCRNG3();
  *bexp = 3 - RUR_AE_SAR();
  *v = 20;
  return shift;
} /* fft_32x16_stage_inner_scl2_DFT5() */
