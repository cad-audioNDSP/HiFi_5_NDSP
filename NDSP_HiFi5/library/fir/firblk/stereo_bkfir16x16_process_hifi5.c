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
  NatureDSP Signal Processing Library. FIR part
    Real block FIR filter, 16x16-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Filter instance structure. */
#include "stereo_bkfir16x16_common.h"

/*-------------------------------------------------------------------------
  Real FIR filter.
  Computes a real FIR filter (direct-form) using IR stored in vector h. The 
  real data input is stored in vector x. The filter output result is stored 
  in vector y. The filter calculates N output samples using M coefficients 
  and requires last M-1 samples in the delay line.
  NOTE: 
  1. User application is not responsible for management of delay lines
  2. User has an option to set IR externally or copy from original location 
     (i.e. from the slower constant memory). In the first case, user is 
     responsible for right alignment, ordering and zero padding of filter 
     coefficients - usually array is composed from zeroes (left padding), 
     reverted IR and right zero padding.


  Precision: 
  16x16    16-bit data, 16-bit coefficients, 16-bit outputs. Ordinary variant 
           and stereo
  32x16    32-bit data, 16-bit coefficients, 32-bit outputs
  32x32    32-bit data, 32-bit coefficients, 32-bit outputs. Ordinary variant 
           and stereo
  32x32ep  the same as above but using 72-bit accumulator for intermediate 
           computations
  f        floating point. Ordinary variant and stereo

  Input:
  x[N*S]   input samples, Q31, Q15, floating point
  h[M]     filter coefficients in normal order, Q31, Q15, floating point
  hl[M]    for stereo filters: filter coefficients for left channel
  hr[M]    for stereo filters: filter coefficients for right channel
  N        length of sample block, should be a multiple of 4
  M        length of filter, should be a multiple of 4
  extIR    if zero, IR is copied from original location, otherwise not
           but user should keep alignment, order of coefficients 
           and zero padding requirements shown below
  S        1 for ordinary (single channel) filters, 2 - for stereo variant
  
  Output:
  y[N*S]   output samples, Q31, Q15, floating point

  Alignment, ordering and zero padding for external IR  (extIR!=0)
  ------------------------+----------+--------------+--------------+----------------
  Function                |Alignment,|Left zero     |   Coefficient| Right zero 
                          | bytes    |padding, bytes|   order      | padding, bytes
  ------------------------+----------+--------------+--------------+----------------
  bkfir16x16_init         |    16    |      2       |  inverted    |  6
  bkfir32x16_init         |    16    |      2       |  inverted    |  6
  bkfir32x32_init         |    16    |      4       |  inverted    |  12
  bkfir32x32ep_init       |    16    |      4       |  inverted    |  12
  bkfirf_init             |    16    |      4       |  inverted    |  12
  stereo_bkfir16x16_init  |    16    |      2       |  inverted    |  6
  stereo_bkfir32x32_init  |    16    |      4       |  inverted    |  12
  stereo_bkfirf_init      |    16    |      4       |  inverted    |  12
  ------------------------+----------+--------------+--------------+----------------

  Restrictions:
  x, y     should not be overlapping
  x, h     aligned on a 16-bytes boundary
  N, M     multiples of 4 
-------------------------------------------------------------------------*/

#define MAX_BUFFER_SZ ((int)(MAX_ALLOCA_SZ/(sizeof(int16_t)*4)))

void stereo_bkfir16x16_process( stereo_bkfir16x16_handle_t handle,
                         int16_t * restrict  y,
                   const int16_t * restrict  x, int N)
{
    stereo_bkfir16x16_ptr_t stereo_bkfir = (stereo_bkfir16x16_ptr_t)handle;
	const ae_int16x8 * restrict pX;
	      ae_int16x8 * restrict pY;
	      ae_int16x4 * restrict pBxl;
	      ae_int16x4 * restrict pBxr;
	      ae_int16x4 * restrict pByl;
	      ae_int16x4 * restrict pByr;
	ae_int16x4 x0, x1;
	ae_int16x4 xl, xr;
	ae_int16x4 y0, y1;
	ae_int16x4 yl, yr;
    ae_valignx2 yv;
	int n;
	int b;
	// Allocate the four buffers
	ALIGN(16) int16_t bxl[MAX_BUFFER_SZ];
	ALIGN(16) int16_t bxr[MAX_BUFFER_SZ];
	ALIGN(16) int16_t byl[MAX_BUFFER_SZ];
	ALIGN(16) int16_t byr[MAX_BUFFER_SZ];

    static const ALIGN(16) int16_t Sel[2 * 4] = { 0x0706, 0x0504, 0x0302, 0x0100, 0x0705, 0x0301, 0x0604, 0x0200 };
    ae_int16x4 sel1, sel2;
    sel1 = AE_L16X4_I((ae_int16x4*)&Sel, 0);
    sel2 = AE_L16X4_I((ae_int16x4*)&Sel, 4 * sizeof(int16_t));

    NASSERT(stereo_bkfir && y && x);
	NASSERT(stereo_bkfir->magic == STEREO_BKFIR16X16_MAGIC);
    NASSERT_ALIGN(x, 16);
    NASSERT(N % 4 == 0);
    if (N <= 0) return;
    yv = AE_ZALIGN128();

	pX   = (const ae_int16x8 *)(x);
	pY   = (      ae_int16x8 *)(y);

	// Aggregation: @ MAX_BUFFER_SZ stereo-samples
    for (b = 0; b < N; b += MAX_BUFFER_SZ)
    {
        int NN = XT_MIN(N - b, MAX_BUFFER_SZ);
        pBxl = (ae_int16x4 *)(bxl);
        pBxr = (ae_int16x4 *)(bxr);
        pByl = (ae_int16x4 *)(byl);
        pByr = (ae_int16x4 *)(byr);

        // Split x into xl (left channel) and xr (right channel).
        for (n = 0; n < (NN >> 3); n++)
        {
            ae_int16x4 xl0, xl1, xr0, xr1;
            AE_L16X4X2_IP(x0, x1, pX, 8 * sizeof(int16_t));
            AE_DSEL16X4(xl0, xr0, x0, x1, sel1);
            AE_L16X4X2_IP(x0, x1, pX, 8 * sizeof(int16_t));
            AE_DSEL16X4(xl1, xr1, x0, x1, sel1);
            AE_S16X4X2_IP(xl0, xl1, castxcc(ae_int16x8, pBxl), 8 * sizeof(int16_t));
            AE_S16X4X2_IP(xr0, xr1, castxcc(ae_int16x8, pBxr), 8 * sizeof(int16_t));
        }
        if (NN & 4)
        {
            AE_L16X4X2_IP(x0, x1, pX, 8 * sizeof(int16_t));
            AE_DSEL16X4(xl, xr, x0, x1, sel1);
            AE_S16X4_IP(xl, pBxl, sizeof(ae_int16x4));
            AE_S16X4_IP(xr, pBxr, sizeof(ae_int16x4));
        }
		pBxl = (      ae_int16x4 *)bxl;
		pBxr = (      ae_int16x4 *)bxr;

		// Aggregate left channel.
		bkfir16x16_process(stereo_bkfir->bkfir_left, (int16_t *)pByl, (const int16_t *)pBxl, NN);

		// Aggregate right channel.
		bkfir16x16_process(stereo_bkfir->bkfir_right, (int16_t *)pByr, (const int16_t *)pBxr, NN);

		// Merge yl (left channel) and yr (right channel) into y.
		for (n = 0; n < (NN >> 3); n++)
		{
            ae_int16x4 yl0, yl1, yr0, yr1;
            AE_L16X4X2_IP(yl0, yl1, castxcc(ae_int16x8, pByl), 8 * sizeof(int16_t));
            AE_L16X4X2_IP(yr0, yr1, castxcc(ae_int16x8, pByr), 8 * sizeof(int16_t));
            AE_DSEL16X4(y0, y1, yl0, yr0, sel2);
            AE_SA16X4X2_IP(y0, y1, yv, pY);
            AE_DSEL16X4(y0, y1, yl1, yr1, sel2);
            AE_SA16X4X2_IP(y0, y1, yv, pY);
		}
        if (NN & 4)
        {
            AE_L16X4_IP(yl, pByl, sizeof(ae_int16x4));
            AE_L16X4_IP(yr, pByr, sizeof(ae_int16x4));
            AE_DSEL16X4(y0, y1, yl, yr, sel2);
            AE_SA16X4X2_IP(y0, y1, yv, pY);
        }
    }
    AE_SA128POS_FP(yv, pY);
} /* stereo_bkfir16x16_process() */
