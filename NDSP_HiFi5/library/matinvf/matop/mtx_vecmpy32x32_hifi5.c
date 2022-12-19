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
/* Code optimized for HiFi5 core */
#include "NatureDSP_Signal_matop.h"
#include "NatureDSP_types.h"
#include "common.h"

/*-------------------------------------------------------------------------
  Matrix by Vector Multiply
  These functions compute the expression z = 2^lsh * x * y for the matrices 
  x and vector y. 
  NOTE: lsh factor is not relevant for floating point routines.

  Two versions of functions available: regular version (mtx_vecmpy32x32,  
  mtx_vecmpy16x16, mtx_vecmpy8x8, mtx_vecmpy8x16, mtx_vecmpyf) with arbitrary 
  arguments and faster version (mtx_vecmpy32x32_fast, mtx_vecmpy16x16_fast, 
  mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast,  mtx_vecmpyf_fast) that apply 
  some restrictions

  Precision: 
  32x32 32-bit input, 32-bit output
  16x16 16-bit input, 16-bit output
  8x8   8-bit inputs, 8-bit output
  8x16  8/16-bit inputs, 16-bit output
  f     floating point

  Input:
  x[M*N] input matrix,Q31,Q15 or floating point
  y[N]   input vector,Q31,Q15 or floating point
  M      number of rows in matrix x
  N      number of columns in matrix x
  lsh    additional left shift(applied to the fixed-
         point functions only) 
  Output:
  z[M]   output vector,Q31,Q15 or floating point

  Restriction:
  For regular routines (mtx_vecmpy32x32, mtx_vecmpy16x16, mtx_vecmpy8x8,
  mtx_vecmpy8x16,  mtx_vecmpyf)
  x,y,z should not overlap

  For faster routines  (mtx_vecmpy32x32_fast, mtx_vecmpy16x16_fast, 
  mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast, mtx_vecmpyf_fast)
  x,y,z   should not overlap
  x,y     aligned on 16-byte boundary
  N, M    multiples of 4
  lsh     should be in range:
          -31...31 for mtx_vecmpy32x32, mtx_vecmpy32x32_fast
          -15...15 for mtx_vecmpy16x16, mtx_vecmpy16x16_fast, 
                   mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast   
-------------------------------------------------------------------------*/
void mtx_vecmpy32x32 (  int32_t* restrict z,
               const int32_t* restrict x,
               const int32_t* restrict y,
               int M, int N, int lsh)
{
	ae_valignx2 ax0, ax1, ay;
	ae_valign az;
	const ae_int32x4 * restrict px0;
	const ae_int32x4 * restrict px1;
	const ae_int32x4 * restrict py;
	ae_int32x2 * restrict pz;
	ae_int32x2 Y0, Y1, X00, X01, X10, X11;
	xtbool2 bmask0, bmask1;
	ae_f64 B0, B1, B2, B3;
	int m, n;
	if (N <= 0 || M <= 0)    /* exceptional situation */
	{
		for (m = 0; m<M; m++) z[m] = 0;
		return;
	}

	bmask0 = AE_MOVBA2((N & 2)+((N & 2) >> 1)+2*(int)((N&3)==1));  // 2*(N%4) if N%4<2, else 3 
	bmask1 = AE_MOVBA2(((int)((N & 3) == 3))<<1);  // 2 if (N%4)=3, else 0
	az = AE_ZALIGN64();
	px0 = (const ae_int32x4 *)x;
	py = (const ae_int32x4 *)y;
	pz = (ae_int32x2 *)z;
	for (m = 0; m<(M&~1); m += 2)
	{
		B0 = B1 = B2 = B3 = AE_ZERO64();
		py = (const ae_int32x4 *)y;
		px0 = (const ae_int32x4 *)(x + m*N);
		px1 = (const ae_int32x4 *)(x + (m + 1)*N);
		ax0 = AE_LA128_PP(px0);
		ax1 = AE_LA128_PP(px1);
		ay = AE_LA128_PP(py);
		for (n = 0; n<(N >> 2); n++)
		{
			AE_LA32X2X2_IP(Y0, Y1, ay, py);
			AE_LA32X2X2_IP(X00, X01, ax0, px0);
			AE_LA32X2X2_IP(X10, X11, ax1, px1);
			AE_MULAAF2D32RA_HH_LL(B0, B1, X00, X10, Y0, Y0);
			AE_MULAAF2D32RA_HH_LL(B2, B3, X01, X11, Y1, Y1);
		}
		AE_LA32X2X2_IP(Y0, Y1, ay, py);
		AE_LA32X2X2_IP(X00, X01, ax0, px0);
		AE_LA32X2X2_IP(X10, X11, ax1, px1);

		AE_MOVF32X2(Y0, AE_ZERO32(), bmask0);
		AE_MOVF32X2(Y1, AE_ZERO32(), bmask1);
		AE_MULAAF2D32RA_HH_LL(B0, B1, X00, X10, Y0, Y0);
		AE_MULAAF2D32RA_HH_LL(B2, B3, X01, X11, Y1, Y1);

		X00 = AE_TRUNCA32X2F64S(B0+B2, B1+B3, 16 + lsh);
		AE_SA32X2_IP(X00, az, pz);
	}
	AE_SA64POS_FP(az, pz);
	if (M & 1)
	{
		B0 = B1 = AE_ZERO64();
		py = (const ae_int32x4 *)y;
		px0 = (const ae_int32x4 *)(x + m*N);
		ax0 = AE_LA128_PP(px0);
		ay = AE_LA128_PP(py);
		for (n = 0; n<(N >> 2); n++)
		{
			AE_LA32X2X2_IP(Y0, Y1, ay, py);
			AE_LA32X2X2_IP(X00, X01, ax0, px0);
			AE_MULAAF2D32RA_HH_LL(B0, B1, X00, X01, Y0, Y1);
		}
		AE_LA32X2X2_IP(Y0, Y1, ay, py);
		AE_LA32X2X2_IP(X00, X01, ax0, px0);
		AE_MOVF32X2(Y0, AE_ZERO32(), bmask0);
		AE_MOVF32X2(Y1, AE_ZERO32(), bmask1);
		AE_MULAAF2D32RA_HH_LL(B0, B1, X00, X01, Y0, Y1);

		X00 = AE_TRUNCA32X2F64S(B0+B1, B0+B1, 16 + lsh);
		AE_S32_L_IP(X00, castxcc(ae_int32, pz), 4);
	}
}
