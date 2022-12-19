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
    Real block FIR filter, 32x32-bit
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
#include "stereo_bkfir32x32_common.h"

#define AE_DSEL32X2(out0, out1, in0, in1, sel_mask) \
{                                                   \
    ae_int16x4 tmp0, tmp1;                          \
    tmp0 = AE_MOVINT16X4_FROMINT32X2(in0);          \
    tmp1 = AE_MOVINT16X4_FROMINT32X2(in1);          \
    AE_DSEL16X4(tmp0, tmp1, tmp0, tmp1, sel_mask);  \
    out0 = AE_MOVINT32X2_FROMINT16X4(tmp0);         \
    out1 = AE_MOVINT32X2_FROMINT16X4(tmp1);         \
}

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

#define MAX_BUFFER_SZ ((int)(MAX_ALLOCA_SZ/(sizeof(int32_t)*4)))

void stereo_bkfir32x32_process( stereo_bkfir32x32_handle_t handle,
                         int32_t * restrict  y,
                   const int32_t * restrict  x, int N )
{
#if 0
    stereo_bkfir32x32_ptr_t stereo_bkfir;
    stereo_bkfir = (stereo_bkfir32x32_ptr_t)handle;

	const ae_int32x4 * restrict pX;
	      ae_int32x4 * restrict pY;
	      ae_int32x2 * restrict pBxl;
	      ae_int32x2 * restrict pBxr;
	      ae_int32x2 * restrict pByl;
	      ae_int32x2 * restrict pByr;
	ae_int32x2 x0, x1;
    ae_int32x2 xl0, xl1, xr0, xr1;
	ae_int32x2 y0, y1;
    ae_int32x2 yl0, yl1, yr0, yr1;
    ae_valignx2 yv;
	int n;
	int b;
	// Allocate the four buffers
	ALIGN(16) int32_t bxl[MAX_BUFFER_SZ];
	ALIGN(16) int32_t bxr[MAX_BUFFER_SZ];
	ALIGN(16) int32_t byl[MAX_BUFFER_SZ];
	ALIGN(16) int32_t byr[MAX_BUFFER_SZ];

    static const ALIGN(16) int16_t Sel[4] = { 0x0705, 0x0604, 0x0301, 0x0200 };
    ae_int16x4 sel;
    sel = AE_L16X4_I((ae_int16x4*)&Sel, 0);

    NASSERT(stereo_bkfir && y && x);
    NASSERT(stereo_bkfir->magic == STEREO_BKFIR32X32_MAGIC);
    NASSERT_ALIGN(x, 16);
    if (N <= 0) return;
    NASSERT(N % 4 == 0);
	yv = AE_ZALIGN128();

	pX   = (const ae_int32x4 *)(x);
	pY   = (      ae_int32x4 *)(y);

	// Aggregation: @ MAX_BUFFER_SZ stereo-samples
	for (b = 0; b < N; b+=MAX_BUFFER_SZ)
	{
		int NN = XT_MIN(N-b, MAX_BUFFER_SZ);
		pBxl = (      ae_int32x2 *)(bxl);
		pBxr = (      ae_int32x2 *)(bxr);
		pByl = (      ae_int32x2 *)(byl);
		pByr = (      ae_int32x2 *)(byr);

		// Split x into xl (left channel) and xr (right channel).
        __Pragma("loop_count min=1");
        for (n = 0; n < (NN >> 2); n++)
		{
            AE_L32X2X2_IP(x0, x1, pX, 4 * sizeof(int32_t));
            AE_DSEL32X2(xl0, xr0, x0, x1, sel);
            AE_L32X2X2_IP(x0, x1, pX, 4 * sizeof(int32_t));
            AE_DSEL32X2(xl1, xr1, x0, x1, sel);
            AE_S32X2X2_IP(xl0, xl1, castxcc(ae_int32x4, pBxl), 4 * sizeof(int32_t));
            AE_S32X2X2_IP(xr0, xr1, castxcc(ae_int32x4, pBxr), 4 * sizeof(int32_t));
		}
		pBxl = (      ae_int32x2 *)bxl;
		pBxr = (      ae_int32x2 *)bxr;

		// Aggregate left channel.
        bkfir32x32_process(stereo_bkfir->bkfir_left, (int32_t *)pByl, (const int32_t *)pBxl, NN);

		// Aggregate right channel.
        bkfir32x32_process(stereo_bkfir->bkfir_right, (int32_t *)pByr, (const int32_t *)pBxr, NN);

		// Merge yl (left channel) and yr (right channel) into y.
        __Pragma("loop_count min=1");
		for (n = 0; n < (NN >> 2); n++)
		{
            AE_L32X2X2_IP(yl0, yl1, castxcc(ae_int32x4, pByl), 4 * sizeof(int32_t));
            AE_L32X2X2_IP(yr0, yr1, castxcc(ae_int32x4, pByr), 4 * sizeof(int32_t));
            AE_DSEL32X2(y0, y1, yl0, yr0, sel);
            AE_SA32X2X2_IP(y0, y1, yv, pY);
            AE_DSEL32X2(y0, y1, yl1, yr1, sel);
            AE_SA32X2X2_IP(y0, y1, yv, pY);
		}
	}
	AE_SA128POS_FP(yv, pY);
#else
    stereo_bkfir32x32_t * bkfir = (stereo_bkfir32x32_ptr_t)handle;
    const ae_int32x4 *          Cl;
    const ae_int32x4 *          Cr;
          ae_int32x4 * D_wrl;
          ae_int32x4 * D_wrr;
    const ae_int32x4 * pSl0;
    const ae_int32x4 * pSl1;
    const ae_int32x4 * pSl2;
    const ae_int32x4 * pSr0;
    const ae_int32x4 * pSr1;
    const ae_int32x4 * pSr2;
    const ae_int32x4 * X;
          ae_int32x4 * Y;
    ae_valignx2 Y_va;
    ae_f64 ql0, ql1, ql2, ql3, ql4, ql5, ql6, ql7;
    ae_f64 qr0, qr1, qr2, qr3, qr4, qr5, qr6, qr7;
    ae_int32x2 dl0, dl1, dl2, dl3, dl4, dl5;
    ae_int32x2 dr0, dr1, dr2, dr3, dr4, dr5;
    ae_int32x2 cl0, cl1;
    ae_int32x2 cr0, cr1;
    ae_int32x2 x0, x1, x2, x3;
    ae_int32x2 y0, y1, y2, y3;

    int M;
    int m, n;

    static const ALIGN(16) int16_t Sel[4] = { 0x0705, 0x0604, 0x0301, 0x0200 };
    ae_int16x4 sel;
    sel = AE_L16X4_I((ae_int16x4*)&Sel, 0);

    //ASSERT( bkfir && bkfir->magic == MAGIC && y && x );
    if (N <= 0) return;
    M = bkfir->M;
    ASSERT(!(M & 3) && !(N & 3));
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(bkfir->delayLineLeft, 16);
    NASSERT_ALIGN((bkfir->delayLineLeft + bkfir->delayLen), 16);
    NASSERT_ALIGN(bkfir->delayPosLeft, 16);
    NASSERT_ALIGN(bkfir->coefLeft, 16);

    //
    // Setup pointers and circular delay line buffer.
    //

    X = (const ae_int32x4*)x;
    Y = (      ae_int32x4*)y;

    D_wrl = (ae_int32x4*)bkfir->delayPosLeft;
    D_wrr = (ae_int32x4*)bkfir->delayPosRight;

    WUR_AE_CBEGIN0((uintptr_t)(bkfir->delayLineLeft));
    WUR_AE_CEND0  ((uintptr_t)(bkfir->delayLineLeft + bkfir->delayLen));
    WUR_AE_CBEGIN1((uintptr_t)(bkfir->delayLineRight));
    WUR_AE_CEND1  ((uintptr_t)(bkfir->delayLineRight + bkfir->delayLen));

    Y_va = AE_ZALIGN128();

    for (n = 0; n < (N >> 3); n++)
    {
        // Reset the coefficients pointers. Now they look at the taps corresponding
        // to the oldest samples in the delay lines.
        Cl = (const ae_int32x4*)bkfir->coefLeft;
        Cr = (const ae_int32x4*)bkfir->coefRight;

        ql0 = ql1 = ql2 = ql3 = ql4 = ql5 = ql6 = ql7 = AE_ZERO64();
        qr0 = qr1 = qr2 = qr3 = qr4 = qr5 = qr6 = qr7 = AE_ZERO64();
        
        // Jump over 2x4 oldest 32-bit entries with circular address update.
        // Now the pointers are M+4 samples apart from the newest ones.
        pSl0 = D_wrl;AE_ADDCIRC_XC(castxcc(ae_int64, pSl0), 12 * sizeof(int32_t));
        pSl1 = pSl0; AE_ADDCIRC_XC(castxcc(ae_int64, pSl1),  4 * sizeof(int32_t));
        pSl2 = pSl1; AE_ADDCIRC_XC(castxcc(ae_int64, pSl2),  4 * sizeof(int32_t));
        pSr0 = D_wrr;AE_ADDCIRC_XC1(castxcc(ae_int64, pSr0), 12 * sizeof(int32_t));
        pSr1 = pSr0; AE_ADDCIRC_XC1(castxcc(ae_int64, pSr1),  4 * sizeof(int32_t));
        pSr2 = pSr1; AE_ADDCIRC_XC1(castxcc(ae_int64, pSr2),  4 * sizeof(int32_t));

        /*--------------------------------------------------------------------------------*/

        // Load 4 input stereo-samples.
        // Q23 <- Q(23+8) - 8
        AE_L32X2X2_IP(x0, x1, X, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(x2, x3, X, 4 * sizeof(int32_t));
        AE_DSEL32X2(dl0, dr0, x0, x1, sel);
        AE_DSEL32X2(dl1, dr1, x2, x3, sel);
        // Store 2x4 samples to the delay line buffers with circular address update.
        // Q(23+8) <- Q23 + 8
        AE_S32X2X2_XC(dl0, dl1, D_wrl, 4 * sizeof(int32_t));
        AE_S32X2X2_XC1(dr0, dr1, D_wrr, 4 * sizeof(int32_t));

        // Load 4 input stereo-samples.
        // Q23 <- Q(23+8) - 8
        AE_L32X2X2_IP(x0, x1, X, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(x2, x3, X, 4 * sizeof(int32_t));
        AE_DSEL32X2(dl0, dr0, x0, x1, sel);
        AE_DSEL32X2(dl1, dr1, x2, x3, sel);
        // Store 2x4 samples to the delay line buffers with circular address update.
        // Q(23+8) <- Q23 + 8
        AE_S32X2X2_XC(dl0, dl1, D_wrl, 4 * sizeof(int32_t));
        AE_S32X2X2_XC1(dr0, dr1, D_wrr, 4 * sizeof(int32_t));

        /*--------------------------------------------------------------------------------*/

        //
        // Inner loop: process 2x4 taps for 2x4 accumulators on each trip. Totally we 
        // perform 2*(M+4) MACs for each accumulator, 2x4 of which fall on zero taps
        // inserted into the impulse response during initialization.
        //
        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            // Load 2x4 samples from the delay line. Altogether we have 2x8 samples 
            // residing in 2x4 AE registers.
            AE_L32X2X2_XC(dl0, dl1, pSl0, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(dl2, dl3, pSl1, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(dl4, dl5, pSl2, 4 * sizeof(int32_t));
            AE_L32X2X2_XC1(dr0, dr1, pSr0, 4 * sizeof(int32_t));
            AE_L32X2X2_XC1(dr2, dr3, pSr1, 4 * sizeof(int32_t));
            AE_L32X2X2_XC1(dr4, dr5, pSr2, 4 * sizeof(int32_t));
            
            // Load the next 2x4 tap coefficients.
            // Q23 <- Q(23+8) - 8
            AE_L32X2X2_IP(cl0, cl1, Cl, 4 * sizeof(int32_t));
            AE_L32X2X2_IP(cr0, cr1, Cr, 4 * sizeof(int32_t));

            // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
            AE_MULAFD32X2RA_FIR_H(ql0, ql1, dl0, dl1, cl0);
            AE_MULAFD32X2RA_FIR_H(ql2, ql3, dl1, dl2, cl0);
            AE_MULAFD32X2RA_FIR_H(qr0, qr1, dr0, dr1, cr0);
            AE_MULAFD32X2RA_FIR_H(qr2, qr3, dr1, dr2, cr0);
            // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
            AE_MULAFD32X2RA_FIR_H(ql0, ql1, dl1, dl2, cl1);
            AE_MULAFD32X2RA_FIR_H(ql2, ql3, dl2, dl3, cl1);
            AE_MULAFD32X2RA_FIR_H(qr0, qr1, dr1, dr2, cr1);
            AE_MULAFD32X2RA_FIR_H(qr2, qr3, dr2, dr3, cr1);

            // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
            AE_MULAFD32X2RA_FIR_H(ql4, ql5, dl2, dl3, cl0);
            AE_MULAFD32X2RA_FIR_H(ql6, ql7, dl3, dl4, cl0);
            AE_MULAFD32X2RA_FIR_H(qr4, qr5, dr2, dr3, cr0);
            AE_MULAFD32X2RA_FIR_H(qr6, qr7, dr3, dr4, cr0);
            // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
            AE_MULAFD32X2RA_FIR_H(ql4, ql5, dl3, dl4, cl1);
            AE_MULAFD32X2RA_FIR_H(ql6, ql7, dl4, dl5, cl1);
            AE_MULAFD32X2RA_FIR_H(qr4, qr5, dr3, dr4, cr1);
            AE_MULAFD32X2RA_FIR_H(qr6, qr7, dr4, dr5, cr1);
        }

        /*--------------------------------------------------------------------------------*/

        // 2xQ23 <- 2xQ16.47 - 24 w/ rounding and saturation.
        dl0 = AE_ROUND32X2F48SASYM(ql0, ql1);
        dl1 = AE_ROUND32X2F48SASYM(ql2, ql3);
        dr0 = AE_ROUND32X2F48SASYM(qr0, qr1);
        dr1 = AE_ROUND32X2F48SASYM(qr2, qr3);
        AE_DSEL32X2(y0, y1, dl0, dr0, sel);
        AE_DSEL32X2(y2, y3, dl1, dr1, sel);
        // Store 4 filter outputs.
        // 2xQ(23+8) <- 2xQ23 + 8
        AE_SA32X2X2_IP(y0, y1, Y_va, Y);
        AE_SA32X2X2_IP(y2, y3, Y_va, Y);

        // 2xQ23 <- 2xQ16.47 - 24 w/ rounding and saturation.
        dl0 = AE_ROUND32X2F48SASYM(ql4, ql5);
        dl1 = AE_ROUND32X2F48SASYM(ql6, ql7);
        dr0 = AE_ROUND32X2F48SASYM(qr4, qr5);
        dr1 = AE_ROUND32X2F48SASYM(qr6, qr7);
        AE_DSEL32X2(y0, y1, dl0, dr0, sel);
        AE_DSEL32X2(y2, y3, dl1, dr1, sel);
        // Store 4 filter outputs.
        // 2xQ(23+8) <- 2xQ23 + 8
        AE_SA32X2X2_IP(y0, y1, Y_va, Y);
        AE_SA32X2X2_IP(y2, y3, Y_va, Y);
    }

    /*--------------------------------------------------------------------------------*/
    if (N & 4)
    {
        // Reset the coefficients pointers. Now they look at the taps corresponding
        // to the oldest samples in the delay lines.
        Cl = (const ae_int32x4*)bkfir->coefLeft;
        Cr = (const ae_int32x4*)bkfir->coefRight;

        ql0 = ql1 = ql2 = ql3 = qr0 = qr1 = qr2 = qr3 = AE_ZERO64();
        
        // Jump over 2x4 oldest 32-bit entries with circular address update.
        // Now the pointers are M+4 samples apart from the newest ones.
        pSl0 = D_wrl;AE_ADDCIRC_XC(castxcc(ae_int64, pSl0), 12 * sizeof(int32_t));
        pSl1 = pSl0; AE_ADDCIRC_XC(castxcc(ae_int64, pSl1),  4 * sizeof(int32_t));
        pSr0 = D_wrr;AE_ADDCIRC_XC1(castxcc(ae_int64, pSr0), 12 * sizeof(int32_t));
        pSr1 = pSr0; AE_ADDCIRC_XC1(castxcc(ae_int64, pSr1),  4 * sizeof(int32_t));

        // Load 4 input stereo-samples.
        // Q23 <- Q(23+8) - 8
        AE_L32X2X2_IP(x0, x1, X, 4 * sizeof(int32_t));
        AE_L32X2X2_IP(x2, x3, X, 4 * sizeof(int32_t));
        AE_DSEL32X2(dl0, dr0, x0, x1, sel);
        AE_DSEL32X2(dl1, dr1, x2, x3, sel);
        // Store 2x4 samples to the delay line buffers with circular address update.
        // Q(23+8) <- Q23 + 8
        AE_S32X2X2_XC(dl0, dl1, D_wrl, 4 * sizeof(int32_t));
        AE_S32X2X2_XC1(dr0, dr1, D_wrr, 4 * sizeof(int32_t));

        //
        // Inner loop: process 2x4 taps for 2x4 accumulators on each trip. Totally we 
        // perform 2*(M+4) MACs for each accumulator, 2x4 of which fall on zero taps
        // inserted into the impulse response during initialization.
        //
        __Pragma("loop_count min=2");
        for (m = 0; m < (M >> 2) + 1; m++)
        {
            // Load 2x4 samples from the delay line. Altogether we have 2x8 samples 
            // residing in 2x4 AE registers.
            AE_L32X2X2_XC(dl0, dl1, pSl0, 4 * sizeof(int32_t));
            AE_L32X2X2_XC(dl2, dl3, pSl1, 4 * sizeof(int32_t));
            AE_L32X2X2_XC1(dr0, dr1, pSr0, 4 * sizeof(int32_t));
            AE_L32X2X2_XC1(dr2, dr3, pSr1, 4 * sizeof(int32_t));

            // Load the next 2x4 tap coefficients.
            // Q23 <- Q(23+8) - 8
            AE_L32X2X2_IP(cl0, cl1, Cl, 4 * sizeof(int32_t));
            AE_L32X2X2_IP(cr0, cr1, Cr, 4 * sizeof(int32_t));

            // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
            AE_MULAFD32X2RA_FIR_H(ql0, ql1, dl0, dl1, cl0);
            AE_MULAFD32X2RA_FIR_H(ql2, ql3, dl1, dl2, cl0);
            AE_MULAFD32X2RA_FIR_H(qr0, qr1, dr0, dr1, cr0);
            AE_MULAFD32X2RA_FIR_H(qr2, qr3, dr1, dr2, cr0);
            // 2xQ16.47 <- 2xQ16.47 + 2x[ Q23*Q23 + 1 ] + 2x[ Q23*Q23 + 1 ]
            AE_MULAFD32X2RA_FIR_H(ql0, ql1, dl1, dl2, cl1);
            AE_MULAFD32X2RA_FIR_H(ql2, ql3, dl2, dl3, cl1);
            AE_MULAFD32X2RA_FIR_H(qr0, qr1, dr1, dr2, cr1);
            AE_MULAFD32X2RA_FIR_H(qr2, qr3, dr2, dr3, cr1);
        }

        // 2xQ23 <- 2xQ16.47 - 24 w/ rounding and saturation.
        dl0 = AE_ROUND32X2F48SASYM(ql0, ql1);
        dl1 = AE_ROUND32X2F48SASYM(ql2, ql3);
        dr0 = AE_ROUND32X2F48SASYM(qr0, qr1);
        dr1 = AE_ROUND32X2F48SASYM(qr2, qr3);
        AE_DSEL32X2(y0, y1, dl0, dr0, sel);
        AE_DSEL32X2(y2, y3, dl1, dr1, sel);
        // Store 4 filter outputs.
        // 2xQ(23+8) <- 2xQ23 + 8
        AE_SA32X2X2_IP(y0, y1, Y_va, Y);
        AE_SA32X2X2_IP(y2, y3, Y_va, Y);
    }

    AE_SA128POS_FP(Y_va, Y);

    bkfir->delayPosLeft = (int32_t*)D_wrl;
    bkfir->delayPosRight = (int32_t*)D_wrr;
#endif
} /* stereo_bkfir32x32_process() */
