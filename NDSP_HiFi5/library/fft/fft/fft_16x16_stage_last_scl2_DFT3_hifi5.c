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

#if 0 /* not used in library for now */

#include "NatureDSP_types.h"
#include "common.h"
#include "fft_16x16_stages.h"



/*
*  Last stage of FFT 16x16, radix-3, dynamic scaling
*/
int fft_16x16_stage_last_scl2_DFT3(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{

    int i;
    const ae_int16x4 * restrict px0;
    const ae_int16x4 * restrict px1;
    const ae_int16x4 * restrict px2;
    ae_int16x4 * restrict py0;
    const int N3=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),715827883));//N/3
    const int stride = N3;
    ae_int16x4 x0, x1, x2, c = AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32(0x00006EDA));
    ae_int16x4 x01, x11, x21; 
    const int shift = XT_MAX(0, 2 - *bexp);

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT((stride & 1) == 0);

    WUR_AE_SAR(shift * 0x102);
    px0 = (ae_int16x4 *)x;
    px1 = px0 + stride / 2;
    px2 = px1 + stride / 2;
    py0 = (ae_int16x4 *)y;


    __Pragma("loop_count min=1");
    for (i = 0; i < (stride >> 2); i++)
    {
        AE_L16X4X2_IP(x0, x01, castxcc(ae_int16x8, px0), 2 * sizeof(ae_int16x4));
        AE_L16X4X2_IP(x1, x11, castxcc(ae_int16x8, px1), 2 * sizeof(ae_int16x4));
        AE_L16X4X2_IP(x2, x21, castxcc(ae_int16x8, px2), 2 * sizeof(ae_int16x4));
        /* 6 cycles per pipeline stage in steady state with unroll=1 */


        x0 = AE_SRAA16RS(x0, shift);
        {
            ae_int16x4 s0, s1, d0;

            s0 = x1;
            d0 = x2;
            AE_ADDANDSUBRNG16RAS_S2(s0, d0);

            s1 = AE_ADD16S(x0, s0);
            s0 = AE_SRAI16(s0, 1);
            s0 = AE_SUB16S(x0, s0);

            d0 = AE_MULFC16RAS(c, d0);
            x0 = s1;

            AE_ADDANDSUBRNG16RAS_S1(s0, d0);
            x1 = d0;
            x2 = s0;
        }

        x01 = AE_SRAA16RS(x01, shift);
        {
            ae_int16x4 s0, s1, d0;

            s0 = x11;
            d0 = x21;
            AE_ADDANDSUBRNG16RAS_S2(s0, d0);

            s1 = AE_ADD16S(x01, s0);
            s0 = AE_SRAI16(s0, 1);
            s0 = AE_SUB16S(x01, s0);

            d0 = AE_MULFC16RAS(c, d0);
            x01 = s1;

            AE_ADDANDSUBRNG16RAS_S1(s0, d0);
            x11 = d0;
            x21 = s0;
        }

        AE_S16X4X2_XP(x0, x01, castxcc(ae_int16x8, py0), stride * sizeof(complex_fract16));
        AE_S16X4X2_XP(x1, x11, castxcc(ae_int16x8, py0), stride * sizeof(complex_fract16));
        AE_S16X4X2_XP(x2, x21, castxcc(ae_int16x8, py0), (4 - 2 * stride) * sizeof(complex_fract16));
    }

    v[0] *= 3;
    return shift;
} /* fft_16x16_stage_last_scl2_DFT3() */

#endif

