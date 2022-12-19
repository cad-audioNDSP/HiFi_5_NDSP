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
    Reference C code
	Integrit, 2006-2019
*/
#include "NatureDSP_Signal_fft.h"
#include "NatureDSP_Signal_vector.h"
#include "common.h"
/* Twiddle factor tables and FFT descriptor structure. */
#include "fft_x16_common.h"
#include "fft_16x16_stages.h"

/*
 *  Intermediate stage of FFT/IFFT 16x16, radix-3, dynamic scaling
 */
int fft_16x16_stage_inner_scl2_DFT3(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{ 
    int i, j, _v;
    const int N3=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),715827883));//N/3
    const int shift = XT_MAX(0, 2 - *bexp);
    const int stride = N3;
    const ae_int16x4 * restrict px0;
        ae_int16x4 * restrict py0;
    const ae_int16x4 * restrict ptwd;
    ae_int16x4 x00, x10, x20;
    ae_int16x4 x01, x11, x21;
    ae_int16x4  c = AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32(0x00006EDA));
    int Ninner;
    _v = *v;
    Ninner=stride / _v;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT((_v & 1) == 0);
    {
        /* Reset RANGE register */
        int a, b;
        AE_CALCRNG16(a, b, 0, 3);
        (void)b;
        (void)a;

    }

    WUR_AE_SAR( shift * 0x102); 
    if (_v&2)
    {   // the code if _v not a multiple of 4: in that case load/stores can not be done with 16-byte alignement
        __Pragma("loop_count min=1");
        for (j = 0; j < ((_v>>1)&~1); j+=2)
        {
            ptwd = (const ae_int16x4 *)tw;
            px0 = (const ae_int16x4 *)x + j;
            py0 = (      ae_int16x4 *)y + j;
            __Pragma("loop_count min=1");
            for (i = 0; i < (Ninner); i++)
            {
                ae_int16x4 tw1, tw2, tmp;
                ae_int16x4 s0, s1, d0;
                AE_L16X4_IP(tmp, ptwd, sizeof(tmp));
                tw2 = AE_SEL16_5410(tmp, tmp);
                tw1 = AE_SEL16_7632(tmp, tmp);

                x01=AE_L16X4_I(px0,sizeof(ae_int16x4)); AE_L16X4_XP(x00, px0, stride*sizeof(complex_fract16));
                x11=AE_L16X4_I(px0,sizeof(ae_int16x4)); AE_L16X4_XP(x10, px0, stride*sizeof(complex_fract16));
                x21=AE_L16X4_I(px0,sizeof(ae_int16x4)); AE_L16X4_XP(x20, px0, (_v - 2 * stride)*sizeof(complex_fract16));

                /************************* DFT3 ******************/
                x00 = AE_SRAA16RS(x00, shift);
                s0 = x10;
                d0 = x20;
                AE_ADDANDSUBRNG16RAS_S2(s0, d0);
                s1 = AE_ADD16S(x00, s0);
                s0 = AE_SRAI16(s0, 1);
                s0 = AE_SUB16S(x00, s0);
                d0 = AE_MULFC16RAS(c, d0);
                AE_ADDANDSUBRNG16RAS_S1(s0, d0);
                x00 = s1;
                x10 = d0;
                x20 = s0;
                /************************* DFT3 ******************/
                x01 = AE_SRAA16RS(x01, shift);
                s0 = x11;
                d0 = x21;
                AE_ADDANDSUBRNG16RAS_S2(s0, d0);
                s1 = AE_ADD16S(x01, s0);
                s0 = AE_SRAI16(s0, 1);
                s0 = AE_SUB16S(x01, s0);
                d0 = AE_MULFC16RAS(c, d0);
                AE_ADDANDSUBRNG16RAS_S1(s0, d0);
                x01 = s1;
                x11 = d0;
                x21 = s0;

                x10 = AE_MULFC16RAS(x10, tw1);
                x20 = AE_MULFC16RAS(x20, tw2);
                x11 = AE_MULFC16RAS(x11, tw1);
                x21 = AE_MULFC16RAS(x21, tw2);

                AE_S16X4RNG_I(x01, py0, sizeof(ae_int16x4)); AE_S16X4RNG_XP(x00, py0, _v * sizeof(complex_fract16));
                AE_S16X4RNG_I(x11, py0, sizeof(ae_int16x4)); AE_S16X4RNG_XP(x10, py0, _v * sizeof(complex_fract16));
                AE_S16X4RNG_I(x21, py0, sizeof(ae_int16x4)); AE_S16X4RNG_XP(x20, py0, _v * sizeof(complex_fract16));
            }
        }
        // tail
        ptwd = (const ae_int16x4 *)tw;
        px0 = (const ae_int16x4 *)x + j;
        py0 = (      ae_int16x4 *)y + j;
        __Pragma("loop_count min=1");
        for (i = 0; i < (Ninner); i++)
        {
            ae_int16x4 tw1, tw2, tmp;
            ae_int16x4 s0, s1, d0;
            AE_L16X4_IP(tmp, ptwd, sizeof(tmp));
            tw2 = AE_SEL16_5410(tmp, tmp);
            tw1 = AE_SEL16_7632(tmp, tmp);

            AE_L16X4_XP(x00, px0, stride*sizeof(complex_fract16));
            AE_L16X4_XP(x10, px0, stride*sizeof(complex_fract16));
            AE_L16X4_XP(x20, px0, (_v - 2 * stride)*sizeof(complex_fract16));

            /************************* DFT3 ******************/
            x00 = AE_SRAA16RS(x00, shift);
            s0 = x10;
            d0 = x20;
            AE_ADDANDSUBRNG16RAS_S2(s0, d0);
            s1 = AE_ADD16S(x00, s0);
            s0 = AE_SRAI16(s0, 1);
            s0 = AE_SUB16S(x00, s0);
            d0 = AE_MULFC16RAS(c, d0);
            AE_ADDANDSUBRNG16RAS_S1(s0, d0);
            x00 = s1;
            x10 = d0;
            x20 = s0;
            x10 = AE_MULFC16RAS(x10, tw1);
            x20 = AE_MULFC16RAS(x20, tw2);
            AE_S16X4RNG_XP(x00, py0, _v * sizeof(complex_fract16));
            AE_S16X4RNG_XP(x10, py0, _v * sizeof(complex_fract16));
            AE_S16X4RNG_XP(x20, py0, _v * sizeof(complex_fract16));
        }
    }
    else
    {   // faster variant for _v multiple of 4
        NASSERT(_v%4==0);
        __Pragma("loop_count min=1");
        for (j = 0; j < (_v>>1); j+=2)
        {
            ptwd = (const ae_int16x4 *)tw;
            px0 = (const ae_int16x4 *)x + j;
            py0 = (      ae_int16x4 *)y + j;
            __Pragma("loop_count min=1");
            for (i = 0; i < (Ninner); i++)
            {
                ae_int16x4 tw1, tw2, tmp;
                AE_L16X4_IP(tmp, ptwd, sizeof(tmp));
                tw2 = AE_SEL16_5410(tmp, tmp);
                tw1 = AE_SEL16_7632(tmp, tmp);

                AE_L16X4X2_XP(x00, x01, castxcc(ae_int16x8, px0), stride*sizeof(complex_fract16));
                AE_L16X4X2_XP(x10, x11, castxcc(ae_int16x8, px0), stride*sizeof(complex_fract16));
                AE_L16X4X2_XP(x20, x21, castxcc(ae_int16x8, px0), (_v - 2 * stride)*sizeof(complex_fract16));

                /************************* DFT3 ******************/
                x00 = AE_SRAA16RS(x00, shift);
                {
                    ae_int16x4 s0, s1, d0;

                    s0 = x10;
                    d0 = x20;
                    AE_ADDANDSUBRNG16RAS_S2(s0, d0);
                    s1 = AE_ADD16S(x00, s0);
                    s0 = AE_SRAI16(s0, 1);
                    s0 = AE_SUB16S(x00, s0);
                    d0 = AE_MULFC16RAS(c, d0);
                    AE_ADDANDSUBRNG16RAS_S1(s0, d0);
                    x00 = s1;
                    x10 = d0;
                    x20 = s0;
                }
                /************************* DFT3 ******************/
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
                    AE_ADDANDSUBRNG16RAS_S1(s0, d0);
                    x01 = s1;
                    x11 = d0;
                    x21 = s0;
                }

                x10 = AE_MULFC16RAS(x10, tw1);
                x20 = AE_MULFC16RAS(x20, tw2);
                x11 = AE_MULFC16RAS(x11, tw1);
                x21 = AE_MULFC16RAS(x21, tw2);

                AE_S16X4X2RNG_XP(x00, x01, castxcc(ae_int16x8, py0), _v * sizeof(complex_fract16));
                AE_S16X4X2RNG_XP(x10, x11, castxcc(ae_int16x8, py0), _v * sizeof(complex_fract16));
                AE_S16X4X2RNG_XP(x20, x21, castxcc(ae_int16x8, py0), _v * sizeof(complex_fract16));
            }
        }
    }

    {
        int a, b;
        AE_CALCRNG16(a, b, 0, 3);
        *bexp = 3 - a;
        (void)b;
    }
    *v *= 3;
    return shift;
} /* fft_16x16_stage_inner_scl2_DFT3() */
