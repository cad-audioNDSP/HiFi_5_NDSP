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
 *  First stage of IFFT 16x16, radix-4, dynamic scaling
 */
int ifft_16x16_stage_first_scl2_DFT4(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
	int i;
	ALIGN(32) static const int16_t sel_tab[4] = { 0x705, 0x604, 0x301, 0x200 };
	ae_int16x4 dsel=AE_L16X4_I((const ae_int16x4*)sel_tab,0);
	const ae_int16x4 * restrict px0;
	const ae_int16x4 * restrict px1;
	const ae_int16x4 * restrict px2;
	const ae_int16x4 * restrict px3;
	ae_int16x8 * restrict py0;
	const ae_int16x8 * restrict ptwd;
	const int stride = (N >> 2);
	const int shift = 3 - *bexp;
	ae_valignx2 a1,a3;

	ae_int16x4 x0, x1, x2, x3, scale=AE_SLAA16S(AE_MOVDA16(1),bexp[0]);
	ae_int16x4 tw0102, tw0311, tw1213;
	NASSERT_ALIGN16(x);
	NASSERT_ALIGN16(y);
	NASSERT(tw_step == 1);

	px0 = (const ae_int16x4 *)x;
	px1 = (const ae_int16x4 *)((complex_fract16*)px0 + stride);
	px2 = (const ae_int16x4 *)((complex_fract16*)px1 + stride);
	px3 = (const ae_int16x4 *)((complex_fract16*)px2 + stride);
	py0 = (      ae_int16x8 *)y;
	ptwd = (const ae_int16x8 *)tw;
	/* Reset RANGE register */
	SetDFT4_Scaling(3); 
	a1=AE_LA128_PP(px1);a3=AE_LA128_PP(px3);
	for (i = 0; i < (stride>>2); i++)
	{   
		ae_int16x4 x01, x11, x21, x31;
		ae_int16x4 z0, z1, z2, z3;
		ae_int16x4 z01, z11, z21, z31;
		ae_int16x4 tw0102_1, tw0311_1, tw1213_1;
		ae_int32x2 t0, t1;
		/* load by 32 to swap re <-> im */
		AE_L32X2X2_IP (t0, t1, castxcc(ae_int32x4, px0), 2 * sizeof(ae_int16x4));
		x0 = AE_MOVINT16X4_FROMINT32X2(t0); x01 = AE_MOVINT16X4_FROMINT32X2(t1);
		AE_LA32X2X2_IP(t0, t1, a1, castxcc(ae_int32x4, px1));
		x1 = AE_MOVINT16X4_FROMINT32X2(t0); x11 = AE_MOVINT16X4_FROMINT32X2(t1);
		AE_L32X2X2_IP(t0, t1, castxcc(ae_int32x4, px2), 2 * sizeof(ae_int16x4));
		x2 = AE_MOVINT16X4_FROMINT32X2(t0); x21 = AE_MOVINT16X4_FROMINT32X2(t1);
		AE_LA32X2X2_IP(t0, t1, a3, castxcc(ae_int32x4, px3));
		x3 = AE_MOVINT16X4_FROMINT32X2(t0); x31 = AE_MOVINT16X4_FROMINT32X2(t1);

		AE_L16X4X2_IP(tw0102,     tw0311, ptwd, 2 * sizeof(ae_int16x4)); 
		AE_L16X4X2_IP(tw1213,   tw0102_1, ptwd, 2 * sizeof(ae_int16x4));
		AE_L16X4X2_IP(tw0311_1, tw1213_1, ptwd, 2 * sizeof(ae_int16x4));

		x21 = AE_SLAA16(x21, bexp[0]);
		x31 = AE_SLAA16(x31, bexp[0]);
		//x31 =AE_MULP16X16X4S(x31 ,scale);
		x0 =AE_MULP16X16X4S(x0 ,scale);
		x1 =AE_MULP16X16X4S(x1 ,scale);
		x2 =AE_MULP16X16X4S(x2 ,scale);
		x3 =AE_MULP16X16X4S(x3 ,scale);
		x01=AE_MULP16X16X4S(x01,scale);
		x11=AE_MULP16X16X4S(x11,scale);

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

		AE_S16X4X2RNG_IP(z0, z1,   py0, 2 * sizeof(ae_int16x4));
		AE_S16X4X2RNG_IP(z2, z3,   py0, 2 * sizeof(ae_int16x4));
		AE_S16X4X2RNG_IP(z01, z11, py0, 2 * sizeof(ae_int16x4));
		AE_S16X4X2RNG_IP(z21, z31, py0, 2 * sizeof(ae_int16x4));
	}
	if (stride & 3)
	{ 
		ae_int16x4 z0, z1, z2, z3;
		ae_int32x2 t0;
		AE_L32X2_IP(t0, castxcc(ae_int32x2, px0), sizeof(ae_int16x4));
		x0 = AE_MOVINT16X4_FROMINT32X2(t0);
		AE_L32X2_IP(t0, castxcc(ae_int32x2, px1), sizeof(ae_int16x4));
		x1 = AE_MOVINT16X4_FROMINT32X2(t0);
		AE_L32X2_IP(t0, castxcc(ae_int32x2, px2), sizeof(ae_int16x4));
		x2 = AE_MOVINT16X4_FROMINT32X2(t0);
		AE_L32X2_IP(t0, castxcc(ae_int32x2, px3), sizeof(ae_int16x4));
		x3 = AE_MOVINT16X4_FROMINT32X2(t0);


		AE_L16X4_IP(tw0102, castxcc(ae_int16x4,ptwd), sizeof(ae_int16x4));
		AE_L16X4_IP(tw0311, castxcc(ae_int16x4,ptwd), sizeof(ae_int16x4));
		AE_L16X4_IP(tw1213, castxcc(ae_int16x4,ptwd), sizeof(ae_int16x4));
		x0 =AE_MULP16X16X4S(x0 ,scale);
		x1 =AE_MULP16X16X4S(x1 ,scale);
		x2 =AE_MULP16X16X4S(x2 ,scale);
		x3 =AE_MULP16X16X4S(x3 ,scale);
		//  Btf 0                         
		DFT4XI2(x0, x1, x2, x3);

		x1 = AE_MULFC16RAS(x1, tw0102);
		x2 = AE_MULFC16RAS(x2, tw0311);
		x3 = AE_MULFC16RAS(x3, tw1213);

		AE_DSEL16X4(z0,z2,x0,x1,dsel);
		AE_DSEL16X4(z1,z3,x2,x3,dsel);

		AE_S16X4X2RNG_IP(z0, z1, castxcc(ae_int16x8, py0), 2 * sizeof(ae_int16x4));
		AE_S16X4X2RNG_IP(z2, z3, castxcc(ae_int16x8, py0), 2 * sizeof(ae_int16x4));
	}
	// update scaling
  {    int a, b; AE_CALCRNG16(a, b, 0, 3); *bexp = 3 - a; (void)b; }
  v[0] <<= 2;
    return shift;
} /* ifft_16x16_stage_first_scl2_DFT4() */
