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
with dynamic data scaling: 16-bit data, 16-bit twiddle factors
C code optimized for HiFi5
IntegrIT, 2006-2019
*/
#include "NatureDSP_types.h"
#include "common.h"
#include "fft_16x16_stages.h"

/*  Table of twiddle factors for radix-8 butterfly
    kron(exp( -1j*2*pi/8*(1:3)'), [1;1] ) */
ALIGN(32) static const int16_t __fft8_tw1_v2_[] =
{
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x5A82, (int16_t)0xA57E,
    (int16_t)0x0000, (int16_t)0x8000, (int16_t)0x0000, (int16_t)0x8000,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0xA57E,
};

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
*  Last stage of IFFT 16x16, radix-8, dynamic scaling
*/
int ifft_16x16_stage_last_scl2_DFT8(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    int j;
    const ae_int16x8 * restrict px0;
    const ae_int16x8 * restrict px1;
    const ae_int16x8 * restrict px2;
    const ae_int16x8 * restrict px3;
    ae_int32x4 * restrict py0;
    ae_int32x4 * restrict py1;
    ae_int32x4 * restrict py2;
    ae_int32x4 * restrict py3;
    ae_int16x4 t1, t2, t3;
    int stride = N>>3;
    int shift;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);

    if (v[0] == stride)
    {
        /* Last stage of the ordinary fft_cplx_16x16 */
    }
    else
    {
        /* Last stage of the stero_fft_cplx_16x16 */
        stride = v[0];
    }

    {
        int sar;
        int S1 = 0, S2 = 0;
        if (*bexp == 0)
        {
            S1 = 1; S2 = 2;
        }
        else if (*bexp == 1)
        {
            S1 = 1; S2 = 1;
        }
        else if (*bexp == 2)
        {
            S1 = 0; S2 = 1;
        }

        sar = (S1 << 7) + S1 + (S2 << 1) + (S2 << 8);
        WUR_AE_SAR(sar);
        shift = 2 * S1 + S2;
    }
    px0 = (ae_int16x8 *)x;
    px1 = px0 + (stride>>1);
    px2 = px1 + (stride>>1);
    px3 = px2 + (stride>>1);
    py0 = (ae_int32x4 *)y;
    py1 = py0 + (stride>>1);
    py2 = py1 + (stride>>1);
    py3 = py2 + (stride>>1);

    t1 = AE_L16X4_I((ae_int16x4*)__fft8_tw1_v2_, 0);
    t2 = AE_L16X4_I((ae_int16x4*)__fft8_tw1_v2_, sizeof(t1));
    t3 = AE_L16X4_I((ae_int16x4*)__fft8_tw1_v2_, 2 * sizeof(t1));

    __Pragma("loop_count min=1");
    for (j = 0; j < (stride>>2); j++)
    {
        ae_int16x4 x0[2], x1[2], x2[2], x3[2], x4[2], x5[2], x6[2], x7[2];
        ae_int16x4 s0[2], s1[2], s2[2], s3[2];
        ae_int16x4 d0[2], d1[2], d2[2], d3[2];
        AE_L16X4X2_X (x1[0],x1[1], px0, stride*2*sizeof(int16_t));
        AE_L16X4X2_IP(x0[0],x0[1], px0, sizeof(ae_int16x8));
        AE_L16X4X2_X (x3[0],x3[1], px1, stride*2*sizeof(int16_t));
        AE_L16X4X2_IP(x2[0],x2[1], px1, sizeof(ae_int16x8));
        AE_L16X4X2_X (x5[0],x5[1], px2, stride*2*sizeof(int16_t));
        AE_L16X4X2_IP(x4[0],x4[1], px2, sizeof(ae_int16x8));
        AE_L16X4X2_X (x7[0],x7[1], px3, stride*2*sizeof(int16_t));
        AE_L16X4X2_IP(x6[0],x6[1], px3, sizeof(ae_int16x8));

        DFT4XI2(x0[0], x2[0], x4[0], x6[0]);  
        DFT4XI2(x1[0], x3[0], x5[0], x7[0]);
        DFT4XI2(x0[1], x2[1], x4[1], x6[1]); 
        DFT4XI2(x1[1], x3[1], x5[1], x7[1]);

        x3[0] = AE_MULFC16RAS(x3[0], t1); 
        x5[0] = AE_MULFC16RAS(x5[0], t2);
        x7[0] = AE_MULFC16RAS(x7[0], t3);
        x3[1] = AE_MULFC16RAS(x3[1], t1); 
        x5[1] = AE_MULFC16RAS(x5[1], t2);
        x7[1] = AE_MULFC16RAS(x7[1], t3);

        s0[0] = x0[0];    s2[0] = x4[0];    s1[0] = x2[0];    s3[0] = x6[0];
        d0[0] = x1[0];    d2[0] = x5[0];    d1[0] = x3[0];    d3[0] = x7[0];
        s0[1] = x0[1];    s2[1] = x4[1];    s1[1] = x2[1];    s3[1] = x6[1];
        d0[1] = x1[1];    d2[1] = x5[1];    d1[1] = x3[1];    d3[1] = x7[1];

        AE_ADDANDSUBRNG16RAS_S1(s0[0], d0[0]); 
        AE_ADDANDSUBRNG16RAS_S1(s1[0], d1[0]);
        AE_ADDANDSUBRNG16RAS_S1(s2[0], d2[0]);
        AE_ADDANDSUBRNG16RAS_S1(s3[0], d3[0]);
        AE_ADDANDSUBRNG16RAS_S1(s0[1], d0[1]); 
        AE_ADDANDSUBRNG16RAS_S1(s1[1], d1[1]);
        AE_ADDANDSUBRNG16RAS_S1(s2[1], d2[1]);
        AE_ADDANDSUBRNG16RAS_S1(s3[1], d3[1]);

        x0[0] = s0[0];    x4[0] = d0[0];    x2[0] = s2[0];    x6[0] = d2[0];
        x1[0] = s1[0];    x5[0] = d1[0];    x3[0] = s3[0];    x7[0] = d3[0];
        x0[1] = s0[1];    x4[1] = d0[1];    x2[1] = s2[1];    x6[1] = d2[1];
        x1[1] = s1[1];    x5[1] = d1[1];    x3[1] = s3[1];    x7[1] = d3[1];

        AE_S32X2X2_X (AE_MOVINT32X2_FROMINT16X4(x1[0]),AE_MOVINT32X2_FROMINT16X4(x1[1]), py0, stride*2*sizeof(int16_t));
        AE_S32X2X2_IP(AE_MOVINT32X2_FROMINT16X4(x0[0]),AE_MOVINT32X2_FROMINT16X4(x0[1]), py0, sizeof(ae_int16x8));
        AE_S32X2X2_X (AE_MOVINT32X2_FROMINT16X4(x3[0]),AE_MOVINT32X2_FROMINT16X4(x3[1]), py1, stride*2*sizeof(int16_t));
        AE_S32X2X2_IP(AE_MOVINT32X2_FROMINT16X4(x2[0]),AE_MOVINT32X2_FROMINT16X4(x2[1]), py1, sizeof(ae_int16x8));
        AE_S32X2X2_X (AE_MOVINT32X2_FROMINT16X4(x5[0]),AE_MOVINT32X2_FROMINT16X4(x5[1]), py2, stride*2*sizeof(int16_t));
        AE_S32X2X2_IP(AE_MOVINT32X2_FROMINT16X4(x4[0]),AE_MOVINT32X2_FROMINT16X4(x4[1]), py2, sizeof(ae_int16x8));
        AE_S32X2X2_X (AE_MOVINT32X2_FROMINT16X4(x7[0]),AE_MOVINT32X2_FROMINT16X4(x7[1]), py3, stride*2*sizeof(int16_t));
        AE_S32X2X2_IP(AE_MOVINT32X2_FROMINT16X4(x6[0]),AE_MOVINT32X2_FROMINT16X4(x6[1]), py3, sizeof(ae_int16x8));
    }

    v[0]<<=3;
    return shift;
} /* ifft_16x16_stage_last_scl2_DFT8() */
