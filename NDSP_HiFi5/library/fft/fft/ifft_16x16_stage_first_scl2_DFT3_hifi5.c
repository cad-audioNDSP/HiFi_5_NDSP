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
#include "common.h"
#include "NatureDSP_Signal_vector.h"
/* Twiddle factor tables and FFT descriptor structure. */
#include "fft_x16_common.h"
#include "fft_16x16_stages.h"

int ifft_16x16_stage_first_scl2_DFT3(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v_, int tw_step, int *bexp)
#if 0
{
	int i;
	uint32_t *px = (uint32_t*)x, tmp;
	for (i = 1; i < (N + 1) / 2; i++)
	{
		tmp = px[i];
		px[i] = px[N - i];
		px[N - i] = tmp;
	}
	return fft_16x16_stage_first_scl2_DFT3/*fft_16x16_stage_inner_scl2_DFT3*/(tw, x, y, N, v_, tw_step, bexp);
}
#else
{
	int i;
	const int N3 = AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N), 715827883));//N/3
	const int stride = N3;
	const int v = *v_;
	ae_valignx2 va0, va1, va2;
	const ae_int16x4 * restrict px0t = (const ae_int16x4 *)(x + 2 * 2 * stride + 2);
	const ae_int16x4 * restrict px1t = (const ae_int16x4 *)(x + 1 * 2 * stride + 2);
	const ae_int16x4 * restrict px2t = (const ae_int16x4 *)(x + 2);
	ae_int16x4 * restrict py = (ae_int16x4 *)(y + N * 2 - 8);
	ae_int16x4 * restrict ptw = (ae_int16x4 *)(tw + N3 * tw_step * 4 - 4);
	const int shift = XT_MAX(0, 3 - *bexp);
	ae_int16x4 scale_shift;

	ALIGN(32) static const int16_t sel_tab[4] = { 0x301, 0x200, 0x705, 0x604 };
	ae_int16x4 dsel = AE_L16X4_I((const ae_int16x4*)sel_tab, 0);


	NASSERT(shift>-32 && shift<32);
	NASSERT(N % 12 == 0);
	// only in  the first stage!!
	NASSERT(v != stride);
	NASSERT(v == 1); (void)v;
	/* Reset RANGE register */
	WUR_AE_SAR(0);//shift * 0x102); 
	scale_shift = AE_SLAA16S(AE_MOVDA16(0x1000), 3 - shift);

	AE_SETCBEGIN0(x); /*circular buffer*/
	AE_SETCEND0(x + 3 * 2 * stride);

	AE_LA16X4X2POS_PC(va0, castxcc(ae_int16x8, px0t));
	va1 = AE_LA128_PP(px1t);
	va2 = AE_LA128_PP(px2t);

	for (i = 0; i < N3; i += 4)
	{
		ae_int16x4 X0[2], X1[2], X2[2], S[2], D[2], Y0[2], Y1[2], Y2[2], T1[2], T2[2];
		AE_L16X4_XP(X1[1], ptw, -tw_step*sizeof(ae_int16x4));
		AE_L16X4_XP(X0[1], ptw, -tw_step*sizeof(ae_int16x4));
		AE_L16X4_XP(X1[0], ptw, -tw_step*sizeof(ae_int16x4));
		AE_L16X4_XP(X0[0], ptw, -tw_step*sizeof(ae_int16x4));
		AE_DSEL16X4(T1[0], T2[0], X0[0], X1[0], dsel);
		AE_DSEL16X4(T1[1], T2[1], X0[1], X1[1], dsel);

		AE_LA16X4X2_IC(X0[1], X0[0], va0, castxcc(ae_int16x8, px0t));
		AE_LA16X4X2_IP(X1[1], X1[0], va1, castxcc(ae_int16x8, px1t));
		AE_LA16X4X2_IP(X2[1], X2[0], va2, castxcc(ae_int16x8, px2t));
		/*
		DFT3 algorithm:
		x - input complex vector
		y - output complex vector
		y = fft(x)
		y = [ x(1) + x(2)  + x(3);
		x(1) + (x(2) + x(3))*cos(2*pi/3) - 1j*(x(2) - x(3))*sin(2*pi/3);
		x(1) + (x(2) + x(3))*cos(2*pi/3) + 1j*(x(2) - x(3))*sin(2*pi/3) ]
		*/
		X0[0] = AE_SRAA16RS(X0[0], shift);
		X0[1] = AE_SRAA16RS(X0[1], shift);
		S[0] = AE_MULFP16X4RS(X1[0], scale_shift);
		S[1] = AE_MULFP16X4RS(X1[1], scale_shift);
		D[0] = AE_MULFP16X4RS(X2[0], scale_shift);
		D[1] = AE_MULFP16X4RS(X2[1], scale_shift);
		AE_ADDANDSUBRNG16RAS_S2(S[0], D[0]);
		AE_ADDANDSUBRNG16RAS_S2(S[1], D[1]);
		Y0[0] = AE_ADD16S(X0[0], S[0]);
		Y0[1] = AE_ADD16S(X0[1], S[1]);
		S[0] = AE_SRAI16(S[0], 1);
		S[1] = AE_SRAI16(S[1], 1);
		S[0] = AE_SUB16S(X0[0], S[0]);
		S[1] = AE_SUB16S(X0[1], S[1]);
		D[0] = AE_MULFC16RAS(D[0], AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32(0x00006EDA)));
		D[1] = AE_MULFC16RAS(D[1], AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32(0x00006EDA)));
		AE_ADDANDSUBRNG16RAS_S2(S[0], D[0]);
		AE_ADDANDSUBRNG16RAS_S2(S[1], D[1]);
		Y1[0] = AE_MULFC16RAS(D[0], T1[0]);
		Y1[1] = AE_MULFC16RAS(D[1], T1[1]);
		Y2[0] = AE_MULFC16RAS(S[0], T2[0]);
		Y2[1] = AE_MULFC16RAS(S[1], T2[1]);
		AE_S16X4X2RNG_IP(AE_SEL16_5432(Y2[1], Y0[1]), AE_SEL16_7632(Y1[1], Y2[1]), castxcc(ae_int16x8, py), -(int)sizeof(ae_int16x8));
		AE_S16X4X2RNG_IP(AE_SEL16_7632(Y1[0], Y2[0]), AE_SEL16_5410(Y0[1], Y1[1]), castxcc(ae_int16x8, py), -(int)sizeof(ae_int16x8));
		AE_S16X4X2RNG_IP(AE_SEL16_5410(Y0[0], Y1[0]), AE_SEL16_5432(Y2[0], Y0[0]), castxcc(ae_int16x8, py), -(int)sizeof(ae_int16x8));
	}
	*v_ = 3;
	// update range
	{    int a, b;    AE_CALCRNG16(a, b, 0, 3);    *bexp = 3 - a;    (void)b;   }
	return shift;
}
#endif
