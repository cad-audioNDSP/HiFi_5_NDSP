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
#if 0 /* not used in library for now */

#include "NatureDSP_types.h"
#include "common.h"
#include "fft_16x16_stages.h" 

/*
DFT3 algorithm:
x - input complex vector
y - output complex vector
y = fft(x)
y = [ x(1) + x(2)  + x(3);
x(1) + (x(2) + x(3))*cos(2*pi/3) - 1j*(x(2) - x(3))*sin(2*pi/3);
x(1) + (x(2) + x(3))*cos(2*pi/3) + 1j*(x(2) - x(3))*sin(2*pi/3) ]

*/
#define DFT3X2(x0, x1, x2)\
{\
    ae_int16x4 s0, s1, d0, c;\
    ae_int32x2 c32;          \
    c32 = AE_MOVDA32(0x6EDA9126);      \
    c = AE_MOVINT16X4_FROMINT32X2(c32);\
    s0 = AE_ADD16S(x1, x2);            \
    d0 = AE_SUB16S(x1, x2);            \
    s1 = AE_ADD16S(x0, s0);            \
    s0 = AE_SRAI16(s0, 1);             \
    s0 = AE_SUB16S(x0, s0);            \
    d0 = AE_MULFP16X4RAS(d0, c);       \
    d0 = AE_SEL16_2301(d0, d0);        \
    x0 = s1;                           \
    x1 = AE_SUB16S(s0, d0);            \
    x2 = AE_ADD16S(s0, d0);            \
}
/*
*  Last stage of IFFT 16x16, radix-3, static scaling
*/
int ifft_16x16_stage_last_scl3_DFT3(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
#define shift 3
    int i;
    const ae_int16x4 * restrict px0;
    const ae_int16x4 * restrict px1;
    const ae_int16x4 * restrict px2;
    ae_int16x4 * restrict py0;
    const int N3=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),715827883));//N/3
    const int stride = N3;
    ae_int16x4 x0, x1, x2;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT((stride & 1) == 0);

    px0 = (ae_int16x4 *)x;
    px1 = px0 + stride / 2;
    px2 = px1 + stride / 2;
    py0 = (ae_int16x4 *)y;

    __Pragma("loop_count min=1");
    for (i = 0; i < (stride >> 1); i++)
    {
        AE_L16X4_IP(x0, px0, sizeof(ae_int16x4));
        AE_L16X4_IP(x1, px1, sizeof(ae_int16x4));
        AE_L16X4_IP(x2, px2, sizeof(ae_int16x4));

        x0 = AE_SRAI16R(x0, shift);
        x1 = AE_SRAI16R(x1, shift);
        x2 = AE_SRAI16R(x2, shift);

        DFT3X2(x0, x1, x2);

        AE_S32X2_XP(AE_MOVINT32X2_FROMINT16X4(x0), castxcc(ae_int32x2, py0), stride*sizeof(complex_fract16));
        AE_S32X2_XP(AE_MOVINT32X2_FROMINT16X4(x1), castxcc(ae_int32x2, py0), stride*sizeof(complex_fract16));
        AE_S32X2_XP(AE_MOVINT32X2_FROMINT16X4(x2), castxcc(ae_int32x2, py0), (2-2*stride)*sizeof(complex_fract16));
    }

    *v = *v * 3;
    return shift;
#undef shift
} /* ifft_16x16_stage_last_scl3_DFT3() */

#endif
