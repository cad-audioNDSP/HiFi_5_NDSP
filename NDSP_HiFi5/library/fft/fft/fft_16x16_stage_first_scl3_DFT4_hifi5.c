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

/*
 *  First stage of FFT 16x16, radix-4, static scaling
 */
int fft_16x16_stage_first_scl3_DFT4(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    /*  AE_SEL16_7632();
        AE_SEL16_5410();    */
    ALIGN(32) static const int16_t sel_tab[4] = {0x705,0x604,0x301,0x200};
    ae_int16x4 dsel=AE_L16X4_I((const ae_int16x4*)sel_tab,0);
    #define shift 3
    int i;
    const ae_int16x4 * restrict px0;
    const ae_int16x4 * restrict px1;
    const ae_int16x4 * restrict px2;
    const ae_int16x4 * restrict px3;
        ae_int16x4 * restrict py0;
    const ae_int16x4 * restrict ptwd;
    const int stride = (N >> 2);
    ae_int16x4 x0, x1, x2, x3;
    ae_int16x4 z0, z1, z2, z3;
    ae_int16x4 tw0102, tw0311, tw1213;
    ae_valignx2 a1,a3;
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT(tw_step == 1);

    SetDFT4_Scaling(shift); 
    px0 = (const ae_int16x4 *)x;
    px1 = (const ae_int16x4 *)((complex_fract16*)px0 + stride);
    px2 = (const ae_int16x4 *)((complex_fract16*)px1 + stride);
    px3 = (const ae_int16x4 *)((complex_fract16*)px2 + stride);
    py0 = (      ae_int16x4 *)y;
    ptwd = (const ae_int16x4 *)tw;
    a1=AE_LA128_PP(px1);
    a3=AE_LA128_PP(px3);

    NASSERT(stride >= 4 ); 
    __Pragma("loop_count min=1");
    for (i = 0; i < (stride >> 2); i++)
    { 
        ae_int16x4 x01, x11, x21, x31;
        ae_int16x4 z01, z11, z21, z31;
        ae_int16x4 tw0102_1, tw0311_1, tw1213_1;

        AE_L16X4X2_IP (x0, x01, castxcc(ae_int16x8, px0), 2 * sizeof(ae_int16x4));
        AE_LA16X4X2_IP(x1, x11, a1,castxcc(ae_int16x8, px1));
        AE_L16X4X2_IP (x2, x21, castxcc(ae_int16x8, px2), 2 * sizeof(ae_int16x4));
        AE_LA16X4X2_IP(x3, x31, a3,castxcc(ae_int16x8, px3));

        AE_L16X4X2_IP(tw0102,     tw0311, castxcc(ae_int16x8, ptwd), 2 * sizeof(ae_int16x4)); 
        AE_L16X4X2_IP(tw1213,   tw0102_1, castxcc(ae_int16x8, ptwd), 2 * sizeof(ae_int16x4));
        AE_L16X4X2_IP(tw0311_1, tw1213_1, castxcc(ae_int16x8, ptwd), 2 * sizeof(ae_int16x4));

        //  Btf 0                         
        DFT4XI2(x0, x1, x2, x3);

        x1 = AE_MULFC16RAS(x1, tw0102);
        x2 = AE_MULFC16RAS(x2, tw0311);
        x3 = AE_MULFC16RAS(x3, tw1213);

        //  Btf 1
        DFT4XI2(x01, x11, x21, x31);

        x11 = AE_MULFC16RAS(x11, tw0102_1);
        x21 = AE_MULFC16RAS(x21, tw0311_1);
        x31 = AE_MULFC16RAS(x31, tw1213_1);

        AE_DSEL16X4(z0,z2,x0,x1,dsel);
        AE_DSEL16X4(z1,z3,x2,x3,dsel);
        AE_DSEL16X4(z01,z21,x01,x11,dsel);
        AE_DSEL16X4(z11,z31,x21,x31,dsel);

        AE_S16X4X2_IP(z0,   z1, castxcc(ae_int16x8, py0), 2 * sizeof(ae_int16x4));
        AE_S16X4X2_IP(z2,   z3, castxcc(ae_int16x8, py0), 2 * sizeof(ae_int16x4));
        AE_S16X4X2_IP(z01, z11, castxcc(ae_int16x8, py0), 2 * sizeof(ae_int16x4));
        AE_S16X4X2_IP(z21, z31, castxcc(ae_int16x8, py0), 2 * sizeof(ae_int16x4));
    }
        // tail
    if (stride & 3)
    { 
        AE_L16X4_IP(x0, px0, sizeof(ae_int16x4));
        AE_L16X4_IP(x1, px1, sizeof(ae_int16x4));
        AE_L16X4_IP(x2, px2, sizeof(ae_int16x4));
        AE_L16X4_IP(x3, px3, sizeof(ae_int16x4));

        AE_L16X4_IP(tw0102, ptwd, sizeof(ae_int16x4));
        AE_L16X4_IP(tw0311, ptwd, sizeof(ae_int16x4));
        AE_L16X4_IP(tw1213, ptwd, sizeof(ae_int16x4));

        //  Btf 0                         
        DFT4XI2(x0, x1, x2, x3);

        x1 = AE_MULFC16RAS(x1, tw0102);
        x2 = AE_MULFC16RAS(x2, tw0311);
        x3 = AE_MULFC16RAS(x3, tw1213);

        AE_DSEL16X4(z0,z2,x0,x1,dsel);
        AE_DSEL16X4(z1,z3,x2,x3,dsel);

        AE_S16X4X2_IP(z0, z1, castxcc(ae_int16x8, py0), 2 * sizeof(ae_int16x4));
        AE_S16X4X2_IP(z2, z3, castxcc(ae_int16x8, py0), 2 * sizeof(ae_int16x4));
    }

  *v = v[0] * 4;
  return shift;
  #undef shift
} /* fft_16x16_stage_first_scl3_DFT4() */
