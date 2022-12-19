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
    C code optimized for HiFi4
    Integrit, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fft.h"
/* Common utility and macros declarations. */
#include "common.h"

/*
Set scaling for DFT4
Range of the 'scale' is 0...3
*/
#define SetDFT4_Scaling(scale)                \
{                                             \
    int sar;                                  \
    NASSERT(scale>=0 && scale<=3);            \
    /*(!"DFT4XI2: scale is out of range"); */ \
    if (scale == 3)        sar = 0x285;       \
    else if (scale == 2)   sar = 0x183;       \
    else if (scale == 1)   sar = 0x102;       \
    else sar = 0;                             \
    WUR_AE_SAR(sar);                          \
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
    Internal stages of fft_cplx16x16_ie, ifft_cplx16x16_ie.
*/
int stage_inner_DFT4_16x16_ie(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp, int isStereo)
{
    int i, j, shift;
    const int R = 4; // stage radix
    const int min_shift = 3;
    const int stride = N / R;
    const int _v = v[0];
    const int tw_inc0 = N / _v / 4 * tw_step* sizeof(complex_fract16);
    int i_count = N / _v / R; 
    int j_count = (_v*(isStereo + 1))>>1;
    ae_int16x8 * restrict _py;
    ae_int16x8 * restrict _px;
    ae_int32  * restrict ptw;

    shift = XT_MAX(0, min_shift - *bexp);
    XT_MOVLTZ(shift, 0, shift);
    ASSERT(shift >= 0 && shift <= 3);

    _py = (ae_int16x8 *)y;
    _px = (ae_int16x8 *)x;

    ptw = (ae_int32*)tw;
    {
        /* Reset RANGE register */
        int a, b;
        AE_CALCRNG16(a, b, 0, 3); (void)b; (void)a;
    }
    SetDFT4_Scaling(shift);
    ASSERT(_v >= 2);
    if (j_count*2 >= i_count)
    {
        __Pragma("loop_count min=2");
        for (i = 0; i < i_count; i++)
        {
            ae_int16x4 tw1, tw2, tw3;
            _py = (ae_int16x8*)(4 * _v * (isStereo + 1) * i  * sizeof(complex_fract16)+(uintptr_t)y);

            tw1 = AE_MOVINT16X4_FROMF32X2(AE_L32_X(ptw, 0));
            tw2 = AE_MOVINT16X4_FROMF32X2(AE_L32_X(ptw, tw_inc0));
            tw3 = AE_MOVINT16X4_FROMF32X2(AE_L32_X(ptw, 2 * tw_inc0));
            ptw += tw_step;

            __Pragma("loop_count min=1");
            for (j = 0; j < (j_count>>1); j ++)
            {
                ae_int16x4 x01, x11, x21, x31;
                ae_int16x4 x00, x10, x20, x30;

                AE_L16X4X2_XP(x00, x01, _px, stride*(isStereo + 1)*sizeof(complex_fract16));
                AE_L16X4X2_XP(x10, x11, _px, stride*(isStereo + 1)*sizeof(complex_fract16));
                AE_L16X4X2_XP(x20, x21, _px, stride*(isStereo + 1)*sizeof(complex_fract16));
                AE_L16X4X2_XP(x30, x31, _px, (4 - 3 * stride*(isStereo + 1))*sizeof(complex_fract16));

                DFT4XI2(x00, x10, x20, x30);
                DFT4XI2(x01, x11, x21, x31);

                x10 = AE_MULFC16RAS(x10, tw1);
                x20 = AE_MULFC16RAS(x20, tw2);
                x30 = AE_MULFC16RAS(x30, tw3);

                x11 = AE_MULFC16RAS(x11, tw1);
                x21 = AE_MULFC16RAS(x21, tw2);
                x31 = AE_MULFC16RAS(x31, tw3);

                AE_S16X4X2RNG_XP(x00, x01, _py, _v*(isStereo + 1)*sizeof(complex_fract16));
                AE_S16X4X2RNG_XP(x10, x11, _py, _v*(isStereo + 1)*sizeof(complex_fract16));
                AE_S16X4X2RNG_XP(x20, x21, _py, _v*(isStereo + 1)*sizeof(complex_fract16));
                AE_S16X4X2RNG_XP(x30, x31, _py, (4  -3 * _v*(isStereo + 1))*sizeof(complex_fract16));
            }
        }
    }
    else /*if (j_count >= i_count) */
    {
        ae_int16x4 tw1, tw2, tw3;

        __Pragma("loop_count min=1");
        for (j = 0; j < j_count; j+=2)
        {
            ptw = (ae_int32*)tw;
            _px = (ae_int16x8*)((0 * j_count + j) * 2 * sizeof(complex_fract16)+(uintptr_t)x);
            _py = (ae_int16x8*)((4 * _v * (isStereo + 1)*0 + 2 * j) * sizeof(complex_fract16)+(uintptr_t)y);
            __Pragma("loop_count min=1");
            for (i = 0; i < i_count; i++)
            {
                ae_int32x2 t32; 
                ae_int16x4 x01, x11, x21, x31;
                ae_int16x4 x00, x10, x20, x30;

                /*
                _py = (ae_int16x4*)((4 *_v * (isStereo + 1)*i + 2*j) * sizeof(complex_fract16)+(uintptr_t)y);
                _px = (ae_int16x4*)((i * j_count + j) * 2 * sizeof(complex_fract16)+(uintptr_t)x);
                ptw = (ae_int32*)((complex_fract16*)tw + tw_step*i);
                */
                tw2 = AE_MOVINT16X4_FROMF32X2(AE_L32_X(ptw, tw_inc0));
                tw3 = AE_MOVINT16X4_FROMF32X2(AE_L32_X(ptw, 2 * tw_inc0));
                AE_L32_XP(t32, ptw, tw_step*sizeof(complex_fract16));
                tw1 = AE_MOVINT16X4_FROMF32X2(t32);

                AE_L16X4X2_XP(x00, x01, _px, stride*(isStereo + 1)*sizeof(complex_fract16));
                AE_L16X4X2_XP(x10, x11, _px, stride*(isStereo + 1)*sizeof(complex_fract16));
                AE_L16X4X2_XP(x20, x21, _px, stride*(isStereo + 1)*sizeof(complex_fract16));
                AE_L16X4X2_XP(x30, x31, _px, 2 * j_count *sizeof(complex_fract16)-3 * stride*(isStereo + 1)*sizeof(complex_fract16));

                DFT4XI2(x00, x10, x20, x30);
                DFT4XI2(x01, x11, x21, x31);

                x10 = AE_MULFC16RAS(x10, tw1);
                x20 = AE_MULFC16RAS(x20, tw2);
                x30 = AE_MULFC16RAS(x30, tw3);

                x11 = AE_MULFC16RAS(x11, tw1);
                x21 = AE_MULFC16RAS(x21, tw2);
                x31 = AE_MULFC16RAS(x31, tw3);

                AE_S16X4X2RNG_XP(x00, x01, _py, _v*(isStereo + 1)*sizeof(complex_fract16));
                AE_S16X4X2RNG_XP(x10, x11, _py, _v*(isStereo + 1)*sizeof(complex_fract16));
                AE_S16X4X2RNG_XP(x20, x21, _py, _v*(isStereo + 1)*sizeof(complex_fract16));
                AE_S16X4X2RNG_XP(x30, x31, _py, _v*(isStereo + 1)*sizeof(complex_fract16));
            }
        }
    } /* if (j_count >= i_count) ... else... */

    {
        int a, b;
        AE_CALCRNG16(a, b, 0, 3); (void)b;
        *bexp = 3 - a;
    }

    *v *= R;
    return shift;
} //stage_inner_DFT4_16x16_ie
