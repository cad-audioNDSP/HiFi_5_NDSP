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

ALIGN(32) static const int16_t __dft5_tw[] =
{
    (int16_t)0x278E, (int16_t)0x278E, (int16_t)0x278E, (int16_t)0x278E,
    (int16_t)0x8644, (int16_t)0x79BC, (int16_t)0x8644, (int16_t)0x79BC,
    (int16_t)0x9872, (int16_t)0x9872, (int16_t)0x9872, (int16_t)0x9872,
    (int16_t)0xB4C3, (int16_t)0x4B3D, (int16_t)0xB4C3, (int16_t)0x4B3D
};
/* twiddles should be loaded from the table above */
#define DFT5X2(x0, x1, x2, x3, x4, w1, w2, w3, w4)\
{ \
    ae_int16x4 s1, s2, d1, d2;             \
    ae_int16x4 t0, t1, t2, t3;             \
    ae_int16x4 y0, y1, y2, y3;             \
    s1 = AE_ADD16S(x1, x4);                \
    s2 = AE_ADD16S(x2, x3);                \
    d1 = AE_SUB16S(x1, x4);                \
    d2 = AE_SUB16S(x2, x3);                \
    \
    t0 = AE_MULFP16X4RAS(s1, w1);         \
    t1 = AE_MULFP16X4RAS(s2, w3);         \
    t2 = AE_MULFP16X4RAS(s1, w3);         \
    t3 = AE_MULFP16X4RAS(s2, w1);         \
    y0 = AE_ADD16S(x0, AE_ADD16S(t0, t1)); \
    y1 = AE_ADD16S(x0, AE_ADD16S(t2, t3)); \
    \
    t0 = AE_MULFP16X4RAS(d1, w2); \
    t1 = AE_MULFP16X4RAS(d2, w4); \
    t2 = AE_MULFP16X4RAS(d2, w2); \
    t3 = AE_MULFP16X4RAS(d1, w4); \
    y2 = AE_ADD16S(t0, t1);      \
    y3 = AE_SUB16S(t3, t2);      \
    y2 = AE_SEL16_2301(y2, y2);  \
    y3 = AE_SEL16_2301(y3, y3);  \
    \
    x0 = AE_ADD16S(x0, AE_ADD16S(s1, s2)); \
    x1 = AE_ADD16S(y0, y2);               \
    x2 = AE_ADD16S(y1, y3);               \
    x3 = AE_SUB16S(y1, y3);               \
    x4 = AE_SUB16S(y0, y2);               \
}

#define DFT5X2_NOSCALING(x0, x1, x2, x3, x4, w1, w2, w3, w4)\
	{ \
    ae_int16x4 s1, s2, d1, d2;             \
    ae_int16x4 t0, t2;                     \
    ae_int16x4 y0, y1, y2, y3;             \
    AE_ADDANDSUBRNG16RAS_S1(x1,x4);        \
    AE_ADDANDSUBRNG16RAS_S1(x2,x3);        \
    s1=x1;d1=x4; s2=x2; d2=x3;             \
    t0 = AE_MULFD16X16X4RAS(s1,s2, w1,w3); \
    t2 = AE_MULFD16X16X4RAS(s1,s2, w3,w1); \
    y0 = AE_ADD16S(x0, t0);                \
    y1 = AE_ADD16S(x0, t2);                \
    y2 = AE_MULFD16X16X4RAS(d1,d2, w2,w4); \
    y3 = AE_MULFD16X16X4RAS(d1,d2, w4,AE_NEG16S(w2)); \
    y2 = AE_SEL16_2301(y2, y2);           \
    y3 = AE_SEL16_2301(y3, y3);           \
    x0 = AE_ADD16S(x0, AE_ADD16S(s1, s2));\
    AE_ADDANDSUBRNG16RAS_S1(y0, y2);      \
    AE_ADDANDSUBRNG16RAS_S1(y1, y3);      \
    x1=y0; x4=y2; x2=y1; x3=y3;           \
	}

/*
!!!!!!!!!! This is obsolete FFTs stage. Shall be removed.
*  First stage of IFFT 16x16, radix-5, dynamic scaling
*/
int ifft_16x16_stage_first_scl2_DFT5(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{ 
	ALIGN(32) static const int16_t sel_tab[4] = { 0x301, 0x200, 0x705, 0x604 };
	//ALIGN(32) static const int16_t sel_tab[4] = { 0x705, 0x604, 0x301, 0x200 };
	ae_int16x4 dsel = AE_L16X4_I((const ae_int16x4*)sel_tab, 0);
	int i;
	const int N5=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),429496730));//N/5
	const int stride = N5;
	const int shift = XT_MAX(0,3 - *bexp);
	ae_valign va0, va1, va2, va3, va4;
	const ae_int16x4 * restrict px0;
	const ae_int16x4 * restrict px1;
	const ae_int16x4 * restrict px2;
	const ae_int16x4 * restrict px3;
	const ae_int16x4 * restrict px4;
	ae_int16x4 * restrict py0;
	const ae_int16x8 * restrict ptwd;
	const ae_int16x8 * restrict ptwd_dft;
	ae_int16x4 x0, x1, x2, x3, x4;
	ae_int16x4 w1, w2, w3, w4;
	NASSERT_ALIGN16(x);
	NASSERT_ALIGN16(y);
	NASSERT((stride & 1) == 0);

	px0 = (ae_int16x4 *)(x + stride * 2 * 4 + 2);
	px1 = (ae_int16x4 *)(x + stride * 2 * 3 + 2);
	px2 = (ae_int16x4 *)(x + stride * 2 * 2 + 2);
	px3 = (ae_int16x4 *)(x + stride * 2 * 1 + 2);
	px4 = (ae_int16x4 *)(x + stride * 2 * 0 + 2);
	AE_SETCBEGIN0(x); /* circular buffer for px0 */
	AE_SETCEND0(x + stride * 2 * 5);
	AE_LA16X4POS_PC(va0, px0);
	va1 = AE_LA64_PP(px1);
	va2 = AE_LA64_PP(px2);
	va3 = AE_LA64_PP(px3);
	va4 = AE_LA64_PP(px4);

	py0 = (ae_int16x4 *)(y + stride * 2 * 5 - 4);
	ptwd = (const ae_int16x8 *)tw + stride - 1;
	ptwd_dft = (const ae_int16x8 *)__dft5_tw;
	AE_L16X4X2_IP(w1, w2,ptwd_dft, sizeof(ae_int16x8));
	AE_L16X4X2_IP(w3, w4,ptwd_dft, sizeof(ae_int16x8));
	/* Reset RANGE register */
	{ int a, b; AE_CALCRNG16(a, b, 0, 3); (void)b; (void)a;  }
	WUR_AE_SAR(shift * 0x102);

	__Pragma("loop_count min=2");
	for (i = 0; i <(stride >> 1); i++)
	{
		ae_int16x4 y0,y1,y2,y3,y4;
		ae_int16x4 tw1,tw2,tw3,tw4;

		AE_LA16X4_IC(x0, va0, px0);
		AE_LA16X4_IP(x1, va1, px1);
		AE_LA16X4_IP(x2, va2, px2);
		AE_LA16X4_IP(x3, va3, px3);
		AE_LA16X4_IP(x4, va4, px4);

		x0 = AE_SRAA16RS(x0, shift);
		x1 = AE_SRAA16RS(x1, shift);
		x2 = AE_SRAA16RS(x2, shift);
		x3 = AE_SRAA16RS(x3, shift);
		x4 = AE_SRAA16RS(x4, shift);
		DFT5X2_NOSCALING(x0, x1, x2, x3, x4, w1, w2, w3, w4);
		//DFT5X2_SCALING(x0, x1, x2, x3, x4, w1, w2, w3, w4);
		/* Load and unpack twiddles */
		AE_L16X4X2_IP(tw2, tw4, ptwd, -(int)sizeof(ae_int16x8));
		AE_L16X4X2_IP(tw1, tw3, ptwd, -(int)sizeof(ae_int16x8));
		AE_DSEL16X4(tw1,tw2,tw1,tw2,dsel);
		AE_DSEL16X4(tw3,tw4,tw3,tw4,dsel);
		y0=x0;
		y1 = AE_MULFC16RAS(x1, tw1);
		y2 = AE_MULFC16RAS(x2, tw2);
		y3 = AE_MULFC16RAS(x3, tw3);
		y4 = AE_MULFC16RAS(x4, tw4);
		// save with permutation

		AE_S16X4RNG_XP(AE_SEL16_7632(y3, y4), py0, -(int)sizeof(ae_int16x4));
		AE_S16X4RNG_XP(AE_SEL16_7632(y1, y2), py0, -(int)sizeof(ae_int16x4));
		AE_S16X4RNG_XP(AE_SEL16_5432(y4, y0), py0, -(int)sizeof(ae_int16x4));
		AE_S16X4RNG_XP(AE_SEL16_5410(y2, y3), py0, -(int)sizeof(ae_int16x4));
		AE_S16X4RNG_XP(AE_SEL16_5410(y0, y1), py0, -(int)sizeof(ae_int16x4));
	}
	// update block exponent
	{   int a, b; AE_CALCRNG16(a, b, 0, 3); *bexp = 3 - a; (void)b; }
	*v *= 5;
	return shift;
} /* ifft_16x16_stage_first_scl2_DFT5() */
