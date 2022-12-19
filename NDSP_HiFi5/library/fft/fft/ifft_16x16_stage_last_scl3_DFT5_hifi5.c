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


/*  16-bit radix-4 butterfly with scaling.
assumes SAR=0x81  */
#define DFT4XI2_scale2(_x0, _x1, _x2, _x3)                 \
{                                       \
    ae_int16x4 s0, s1, d0, d1;          \
    s0 = _x0;    s1 = _x1;              \
    d0 = _x2;    d1 = _x3;              \
    AE_ADDANDSUBRNG16RAS_S1(s0, d0);    \
    AE_ADDANDSUBRNG16RAS_S1(s1, d1);    \
    d1 = AE_MUL16JS(d1);                \
    AE_ADDANDSUBRNG16RAS_S1(s0, s1);    \
    AE_ADDANDSUBRNG16RAS_S1(d0, d1);    \
    _x0 = s0;    _x2 = s1;              \
    _x3 = d0;    _x1 = d1;              \
}

#define DFT4XI2_scale1(_x0, _x1, _x2, _x3)                 \
{                                      \
    ae_int16x4 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG16RAS_S1(_x0, _x2); \
    AE_ADDANDSUBRNG16RAS_S1(_x1, _x3); \
    s0=_x0;d0=_x2;s1=_x1;d1=_x3;       \
    d1 = AE_MUL16JS(d1);               \
    AE_ADDANDSUBRNG16RAS_S2(s0, s1);   \
    AE_ADDANDSUBRNG16RAS_S2(d0, d1);   \
    _x0 = s0;    _x2 = s1;             \
    _x3 = d0;    _x1 = d1;             \
}
#if 0
/*
    Set scaling for DFT4
    Range of the 'scale' is 0...3
*/
#define SetDFT4_Scaling(scale)         \
{                                      \
    int sar = 0;                       \
    NASSERT(scale>=0 && scale<=3);     \
    if (scale == 3)        sar = 0x285;\
    else if (scale == 2)   sar = 0x183;\
    else if (scale == 1)   sar = 0x102;\
    else sar = 0;                      \
    WUR_AE_SAR(sar);                   \
}

/*  16-bit radix-4 butterfly with scaling.
    Call SetDFT4_Scaling() before */
#define DFT4XI2(_x0, _x1, _x2, _x3)                 \
{                                                   \
    ae_int16x4 s0, s1, d0, d1;                      \
    s0 = _x0;    s1 = _x1;                          \
    d0 = _x2;    d1 = _x3;                          \
    AE_ADDANDSUBRNG16RAS_S1(s0, d0);                \
    AE_ADDANDSUBRNG16RAS_S1(s1, d1);                \
    d1 = AE_MUL16JS(d1);                            \
    AE_ADDANDSUBRNG16RAS_S2(s0, s1);                \
    AE_ADDANDSUBRNG16RAS_S2(d0, d1);                \
    _x0 = s0;    _x2 = s1;                          \
    _x3 = d0;    _x1 = d1;                          \
}

/*  16-bit radix-4 butterfly without scaling.     */
#define DFT4XI2_NOSCALING(_x0, _x1, _x2, _x3)       \
{                                                   \
    ae_int16x4 s0, s1, d0, d1;                      \
    ae_int16x4 s10, s11, d10, d11;                  \
    s0 = AE_ADD16S(_x0, _x2);                       \
    s1 = AE_ADD16S(_x1, _x3);                       \
    d0 = AE_SUB16S(_x0, _x2);                       \
    d1 = AE_SUB16S(_x1, _x3);                       \
    d1 = AE_MUL16JS(d1);                            \
    s10 = AE_ADD16S(s0, s1);                        \
    s11 = AE_SUB16S(s0, s1);                        \
    d10 = AE_ADD16S(d0, d1);                        \
    d11 = AE_SUB16S(d0, d1);                        \
    _x0 = s10;    _x2 = s11;                        \
    _x3 = d10;    _x1 = d11;                        \
}                                   
#endif
ALIGN(32) static const int16_t DFT5_twd[16] = 
{
     -8192 ,           0,    -8192 ,           0,
      15582,        9630,     15582,        9630,
    - 18318,           0,   - 18318,           0,
    - 15582,        9630,   - 15582,        9630,
};
int ifft_16x16_stage_last_scl3_DFT5(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
#define shift 3
	int i;
	const int N5 = AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N), 429496730));//N/5
	const int stride = N5;
	const ae_int16x4 * restrict px0;
	const ae_int16x4 * restrict px1;
	const ae_int16x4 * restrict px2;
	const ae_int16x4 * restrict px3;
	const ae_int16x4 * restrict px4;
	ae_int16x4 * restrict py0;
	ae_int16x4 x0, x1, x2, x3, x4;
	ae_int16x4 x01, x11, x21, x31, x41;
	ae_int16x4 a0, a1, a2, a3;

	NASSERT_ALIGN16(x);
	NASSERT_ALIGN16(y);
	NASSERT((stride & 3) == 0);

	/* Load twiddles for DFT5 */
	AE_L16X4X2_I(a0, a1, (ae_int16x8*)DFT5_twd, 0 * sizeof(a0));
	AE_L16X4X2_I(a2, a3, (ae_int16x8*)DFT5_twd, 1 * sizeof(ae_int16x8));

	WUR_AE_SAR(0x81);
	px0 = (ae_int16x4 *)x;
	px1 = px0 + stride / 2;
	px2 = px1 + stride / 2;
	px3 = px2 + stride / 2;
	px4 = px3 + stride / 2;
	py0 = (ae_int16x4 *)y;

	__Pragma("loop_count min=1");
	for (i = 0; i <(stride >> 2); i++)
	{
		ae_int16x4 z00, z01, z02, z03, z0, z1, z2, z3;

		AE_L16X4X2_IP(x0, x01, castxcc(ae_int16x8, px0), 2 * sizeof(ae_int16x4));
		AE_L16X4X2_IP(x1, x11, castxcc(ae_int16x8, px1), 2 * sizeof(ae_int16x4));
		AE_L16X4X2_IP(x2, x21, castxcc(ae_int16x8, px2), 2 * sizeof(ae_int16x4));
		AE_L16X4X2_IP(x3, x31, castxcc(ae_int16x8, px3), 2 * sizeof(ae_int16x4));
		AE_L16X4X2_IP(x4, x41, castxcc(ae_int16x8, px4), 2 * sizeof(ae_int16x4));
		/******************* DFT5 *********************/
		x0 = AE_SRAI16R(x0, shift);
		z00 = x1;    z01 = x2;    z02 = x4;    z03 = x3;
		/* DFT4 with scaling 3*/
		DFT4XI2_scale2(z00, z01, z02, z03);

		/* z1, z3 are swaped! */
		z0 = AE_MULFC16RAS(z00, a0);
		z3 = AE_MULFC16RAS(z01, a1);
		z2 = AE_MULFC16RAS(z02, a2);
		z1 = AE_MULFC16RAS(z03, a3);
		DFT4XI2_scale1(z0, z1, z2, z3);

		x1 = AE_ADD16S(x0, z3);
		x2 = AE_ADD16S(x0, z2);
		x3 = AE_ADD16S(x0, z0);
		x4 = AE_ADD16S(x0, z1);
		x0 = AE_ADD16S(x0, AE_SRAI16R(z00, 1));

		/******************* DFT5 *********************/
		x01 = AE_SRAI16R(x01, shift);
		z00 = x11;    z01 = x21;    z02 = x41;    z03 = x31;
		/* DFT4 with scaling 3*/
		DFT4XI2_scale2(z00, z01, z02, z03);

		/* z1, z3 are swaped! */
		z0 = AE_MULFC16RAS(z00, a0);
		z3 = AE_MULFC16RAS(z01, a1);
		z2 = AE_MULFC16RAS(z02, a2);
		z1 = AE_MULFC16RAS(z03, a3);
		DFT4XI2_scale1(z0, z1, z2, z3);

		x11 = AE_ADD16S(x01, z3);
		x21 = AE_ADD16S(x01, z2);
		x31 = AE_ADD16S(x01, z0);
		x41 = AE_ADD16S(x01, z1);
		x01 = AE_ADD16S(x01, AE_SRAI16R(z00, 1));

		AE_S32X2X2_XP(AE_MOVINT32X2_FROMINT16X4(x0), AE_MOVINT32X2_FROMINT16X4(x01), castxcc(ae_int32x4, py0), stride*sizeof(complex_fract16));
		AE_S32X2X2_XP(AE_MOVINT32X2_FROMINT16X4(x1), AE_MOVINT32X2_FROMINT16X4(x11), castxcc(ae_int32x4, py0), stride*sizeof(complex_fract16));
		AE_S32X2X2_XP(AE_MOVINT32X2_FROMINT16X4(x2), AE_MOVINT32X2_FROMINT16X4(x21), castxcc(ae_int32x4, py0), stride*sizeof(complex_fract16));
		AE_S32X2X2_XP(AE_MOVINT32X2_FROMINT16X4(x3), AE_MOVINT32X2_FROMINT16X4(x31), castxcc(ae_int32x4, py0), stride*sizeof(complex_fract16));
		AE_S32X2X2_XP(AE_MOVINT32X2_FROMINT16X4(x4), AE_MOVINT32X2_FROMINT16X4(x41), castxcc(ae_int32x4, py0), (4 - 4 * stride)*sizeof(complex_fract16));
	}
	*v *= 5;
    return shift;
#undef shift
} /* ifft_16x16_stage_last_scl3_DFT5() */
