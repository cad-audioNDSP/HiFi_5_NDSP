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
    Real block FIR filter, floating point
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "common_fpu.h"
/* Filter instance structure. */
#include "stereo_bkfirf_common.h"

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

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void,stereo_bkfirf_process,( stereo_bkfirf_handle_t handle,
                         float32_t * restrict  y,
                   const float32_t * restrict  x, int N ))
#elif HAVE_VFPU

#define AE_DSELSX2(out0, out1, in0, in1, sel_mask) \
{                                                  \
    ae_int16x4 tmp0, tmp1;                         \
    tmp0 = AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(in0));         \
    tmp1 = AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(in1));         \
    AE_DSEL16X4(tmp0, tmp1, tmp0, tmp1, sel_mask); \
    out0 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(tmp0));        \
    out1 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(tmp1));        \
}

#define AE_LSX2XC1(d, a, off) \
{                             \
    ae_int32x2 d_;             \
    AE_L32X2_XC1(d_, castxcc(ae_int32x2, a), off); \
    d = AE_MOVXTFLOATX2_FROMINT32X2(d_); \
}

#define XT_LSXC1(d, a, off) \
{                             \
    ae_int32x2 d_;  \
    xtfloatx2 t;             \
    AE_L32_XC1(d_, castxcc(ae_int32, a), off); \
    t = AE_MOVXTFLOATX2_FROMINT32X2(d_); \
    d = LOW_S(t); \
}


#define MAX_BUFFER_SZ ((int)(MAX_ALLOCA_SZ/(sizeof(float32_t)*4)))

/* process block of samples */
void stereo_bkfirf_process( stereo_bkfirf_handle_t handle,
                         float32_t * restrict  y,
                   const float32_t * restrict  x, int N )
{
    stereo_bkfirf_t * bkfir = (stereo_bkfirf_ptr_t)handle;

    const xtfloatx4 * restrict pX;
          xtfloatx4 * restrict pDw;
    const xtfloatx4 * restrict pX0;
    const xtfloatx4 * restrict pX1;
    const xtfloatx4 * restrict pX2;
    const xtfloatx4 * restrict pHl;
    const xtfloatx4 * restrict pHr;
          xtfloatx4 * restrict pY;
    ae_valignx2 aY;
    
    int M;
    int n, m;

    static const ALIGN(16) int16_t Sel[4] = { 0x0705, 0x0604, 0x0301, 0x0200 };
    ae_int16x4 sel;
    sel = AE_L16X4_I((ae_int16x4*)&Sel, 0);

    M = bkfir->M;
    NASSERT(bkfir && bkfir->magic == STEREO_BKFIRF_MAGIC && y && x);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN((bkfir->coefLeft), 16);
    NASSERT_ALIGN((bkfir->coefRight), 16);
    NASSERT_ALIGN((bkfir->delayLine), 16);
    NASSERT_ALIGN((bkfir->delayLine + bkfir->delayLen), 16);
    NASSERT_ALIGN((bkfir->delayPos), 16);
    NASSERT(N % 4 == 0 && M % 4 == 0);
    if (N <= 0) return;

    // Setup pointers and circular delay line buffer.
    pX  = (const xtfloatx4 *)x;
    pY  = (      xtfloatx4 *)y;
    pDw= (      xtfloatx4 *)bkfir->delayPos;
    WUR_AE_CBEGIN0((uintptr_t)(bkfir->delayLine));
    WUR_AE_CEND0  ((uintptr_t)(bkfir->delayLine + bkfir->delayLen));
    aY = AE_ZALIGN128();
    for (n = 0; n < (N >> 3); n++)
    {
        xtfloatx2 y00, y10, y20, y30, y40, y50, y60, y70;
        xtfloatx2 y01, y11, y21, y31, y41, y51, y61, y71;
        xtfloatx2 y0, y1, y2, y3, y4, y5, y6, y7;
        xtfloatx2 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10;
        xtfloatx2 cl0, cl1, cr0, cr1;
        xtfloatx2 c0, c1, c2, c3;     

        //Load 4 stereo samples and save them to the delay line
        AE_LSX2X2_IP(x0, x1, pX, 4 * sizeof(float32_t));
        AE_LSX2X2_IP(x2, x3, pX, 4 * sizeof(float32_t));

        AE_SSX2X2_XC(x0, x1, pDw, 4 * sizeof(float32_t));
        AE_SSX2X2_XC(x2, x3, pDw, 4 * sizeof(float32_t));

        //Load 4 stereo samples and save them to the delay line
        AE_LSX2X2_IP(x0, x1, pX, 4 * sizeof(float32_t));
        AE_LSX2X2_IP(x2, x3, pX, 4 * sizeof(float32_t));

        AE_SSX2X2_XC(x0, x1, pDw, 4 * sizeof(float32_t));
        AE_SSX2X2_XC(x2, x3, pDw, 4 * sizeof(float32_t));

        pHl = (const xtfloatx4 *)bkfir->coefLeft;
        pHr = (const xtfloatx4 *)bkfir->coefRight;

        pX0 = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pX0), 8 * sizeof(int32_t));
        pX1 = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pX1), 16 * sizeof(int32_t));
        pX2 = pDw;
        //preload samples
        AE_LSX2X2_XC(x0, x1, pX2, 4*sizeof(float32_t)); 
        AE_LSX2X2_XC(x2, x3, pX2, 4*sizeof(float32_t));
        AE_LSX2X2_XC(x4, x5, pX2, 4*sizeof(float32_t));
        AE_LSX2X2_XC(x6, x7, pX2, 4*sizeof(float32_t)); 
        AE_LSX2X2_XC(x8, x9, pX2, 4*sizeof(float32_t));
        AE_LSX2XC(x10, castxcc(xtfloatx2, pX2), 4*sizeof(float32_t));            
        
        //first iteration
        {
            AE_LSX2X2_IP(cl0, cl1, pHl, 4*sizeof(float32_t));
            AE_LSX2X2_IP(cr0, cr1, pHr, 4*sizeof(float32_t));

            AE_DSELSX2(c0, c1, cl0, cr0, sel);
            AE_DSELSX2(c2, c3, cl1, cr1, sel);

            MULQ_S(y01, y11, x1, x2, c1);
            MULQ_S(y21, y31, x3, x4, c1);
            MULQ_S(y41, y51, x5, x6, c1);
            MULQ_S(y61, y71, x7, x8, c1);

            AE_LSX2X2_XC(x0, x1, pX0, 4*sizeof(float32_t)); 

            MULQ_S(y00, y10, x2, x3, c2);  
            MULQ_S(y20, y30, x4, x5, c2); 
             

            MADDQ_S(y01, y11, x3, x4, c3);

            AE_LSX2X2_XC(x2, x3, pX0, 4*sizeof(float32_t)); 

            MADDQ_S(y21, y31, x5, x6, c3);

            AE_LSX2X2_XC(x4, x5, pX1, 4*sizeof(float32_t));

            MULQ_S(y40, y50, x6, x7, c2); 
            MULQ_S(y60, y70, x8, x9, c2); 

            MADDQ_S(y41, y51, x7, x8, c3);

            AE_LSX2X2_XC(x6, x7, pX1, 4*sizeof(float32_t));

            MADDQ_S(y61, y71, x9, x10, c3);

            AE_LSX2X2_XC(x8,  x9, pX2, 4*sizeof(float32_t)); 
            AE_LSX2XC(x10, castxcc(xtfloatx2, pX2), 4*sizeof(float32_t)); 
        }

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        { /* 8 cycles unroll=1 MAC bound */
            AE_LSX2X2_IP(cl0, cl1, pHl, 4*sizeof(float32_t));
            AE_LSX2X2_IP(cr0, cr1, pHr, 4*sizeof(float32_t));
            c0 = XT_SEL32_HH_SX2(cl0, cr0);

            MADDQ_S(y00, y10, x0, x1, c0);
            MADDQ_S(y20, y30, x2, x3, c0);
            MADDQ_S(y40, y50, x4, x5, c0);
            MADDQ_S(y60, y70, x6, x7, c0);

            c1 = XT_SEL32_LL_SX2(cl0, cr0);

            MADDQ_S(y01, y11, x1, x2, c1);
            MADDQ_S(y21, y31, x3, x4, c1);
            MADDQ_S(y41, y51, x5, x6, c1);
            MADDQ_S(y61, y71, x7, x8, c1);

            AE_LSX2X2_XC(x0, x1, pX0, 4*sizeof(float32_t)); 

            c2 = XT_SEL32_HH_SX2(cl1, cr1);

            MADDQ_S(y00, y10, x2, x3, c2);  
            MADDQ_S(y20, y30, x4, x5, c2); 

            c3 = XT_SEL32_LL_SX2(cl1, cr1);                

            MADDQ_S(y01, y11, x3, x4, c3);

            AE_LSX2X2_XC(x2, x3, pX0, 4*sizeof(float32_t)); 

            MADDQ_S(y21, y31, x5, x6, c3);

            AE_LSX2X2_XC(x4, x5, pX1, 4*sizeof(float32_t));

            MADDQ_S(y40, y50, x6, x7, c2); 
            MADDQ_S(y60, y70, x8, x9, c2); 

            MADDQ_S(y41, y51, x7, x8, c3);

            AE_LSX2X2_XC(x6, x7, pX1, 4*sizeof(float32_t));

            MADDQ_S(y61, y71, x9, x10, c3);

            AE_LSX2X2_XC(x8,  x9, pX2, 4*sizeof(float32_t)); 
            AE_LSX2XC(x10, castxcc(xtfloatx2, pX2), 4*sizeof(float32_t)); 
        }
        //add accumulators
        ADD_SX2X2(y0, y1, y00, y10, y01, y11);
        ADD_SX2X2(y2, y3, y20, y30, y21, y31);
        ADD_SX2X2(y4, y5, y40, y50, y41, y51);
        ADD_SX2X2(y6, y7, y60, y70, y61, y71);

        /* save computed samples */
        AE_SASX2X2_IP(y0, y1, aY, pY);
        AE_SASX2X2_IP(y2, y3, aY, pY);

        AE_SASX2X2_IP(y4, y5, aY, pY);
        AE_SASX2X2_IP(y6, y7, aY, pY);
    }
    AE_SA128POS_FP(aY, pY);

    //last 4 stereo samples 
    if ( N & 4)
    {
        xtfloatx2 y00, y10, y20, y30;
        xtfloatx2 y01, y11, y21, y31;
        xtfloatx2 y0, y1, y2, y3;
        xtfloatx2 x0, x1, x2, x3, x4, x5, x6, t0;
        xtfloatx2 cl0, cl1, cr0, cr1;
        xtfloatx2 c0, c1, c2, c3;        

        //Load 4 stereo samples and save them to the delay line
        AE_LSX2X2_IP(x0, x1, pX, 4 * sizeof(float32_t));
        AE_LSX2X2_IP(x2, x3, pX, 4 * sizeof(float32_t));

        AE_SSX2X2_XC(x0, x1, pDw, 4 * sizeof(float32_t));
        AE_SSX2X2_XC(x2, x3, pDw, 4 * sizeof(float32_t));

        pHl = (const xtfloatx4 *)bkfir->coefLeft;
        pHr = (const xtfloatx4 *)bkfir->coefRight;

        CONST_SX2X2(y00, y10, 0);
        CONST_SX2X2(y20, y30, 0);
        CONST_SX2X2(y01, y11, 0);
        CONST_SX2X2(y21, y31, 0);

        pX0 = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pX0), 8 * sizeof(int32_t));
        //preload samples
        AE_LSX2X2_XC(x0, x1, pX0, 4*sizeof(float32_t)); 
        AE_LSX2X2_XC(x2, x3, pX0, 4*sizeof(float32_t));
        AE_LSX2X2_XC(x4, x5, pX0, 4*sizeof(float32_t));
        AE_LSX2X2_XC(x6, t0, pX0, 4*sizeof(float32_t));              

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2)+1; m++)
        {
            AE_LSX2X2_IP(cl0, cl1, pHl, 4*sizeof(float32_t));
            AE_LSX2X2_IP(cr0, cr1, pHr, 4*sizeof(float32_t));
            AE_DSELSX2(c0, c1, cl0, cr0, sel);
            AE_DSELSX2(c2, c3, cl1, cr1, sel);

            MADD_SX2X2(y00, y10, x0, x1, c0, c0);
            MADD_SX2X2(y01, y11, x1, x2, c1, c1);
            MADD_SX2X2(y00, y10, x2, x3, c2, c2);
            MADD_SX2X2(y01, y11, x3, x4, c3, c3);

            MADD_SX2X2(y20, y30, x2, x3, c0, c0);
            MADD_SX2X2(y21, y31, x3, x4, c1, c1);
            MADD_SX2X2(y20, y30, x4, x5, c2, c2);
            MADD_SX2X2(y21, y31, x5, x6, c3, c3);

            x0 = x4; x1 = x5; x2 = x6; x3 = t0;
            AE_LSX2X2_XC(x4, x5, pX0, 4*sizeof(float32_t)); 
            AE_LSX2X2_XC(x6, t0, pX0, 4*sizeof(float32_t)); 
        }
        ADD_SX2X2(y0, y1, y00, y10, y01, y11);
        ADD_SX2X2(y2, y3, y20, y30, y21, y31);

        /* save computed samples */
        AE_SASX2X2_IP(y0, y1, aY, pY);
        AE_SASX2X2_IP(y2, y3, aY, pY);
    }
    AE_SA128POS_FP(aY, pY);        
    
    bkfir->delayPos = (float32_t*)pDw;

} /* stereo_bkfirf_process() */
#else

#define MAX_BUFFER_SZ ((int)(MAX_ALLOCA_SZ/(sizeof(float32_t)*4)))

/* process block of samples */
void stereo_bkfirf_process( stereo_bkfirf_handle_t _bkfir,
                         float32_t * restrict  y,
                   const float32_t * restrict  x, int N )
{
    stereo_bkfirf_ptr_t stereo_bkfir;
	stereo_bkfir = (stereo_bkfirf_ptr_t)_bkfir;

	const xtfloat * restrict pX;
	      xtfloat * restrict pY;
	      xtfloat * restrict pBxl;
	      xtfloat * restrict pBxr;
	      xtfloat * restrict pByl;
	      xtfloat * restrict pByr;
	xtfloat xl, xr;
	xtfloat yl, yr;
	int n;
	int b;
	// Allocate the four buffers
	ALIGN(16) float32_t bxl[MAX_BUFFER_SZ];
	ALIGN(16) float32_t bxr[MAX_BUFFER_SZ];
	ALIGN(16) float32_t byl[MAX_BUFFER_SZ];
	ALIGN(16) float32_t byr[MAX_BUFFER_SZ];

    NASSERT(stereo_bkfir && y && x);
	NASSERT(stereo_bkfir->magic == STEREO_BKFIRF_MAGIC);
    NASSERT_ALIGN(x, 8);
    if (N <= 0) return;
    NASSERT(N % 4 == 0);

	pX   = (const xtfloat *)(x);
	pY   = (      xtfloat *)(y);

	// Aggregation: @ MAX_BUFFER_SZ stereo-samples
	for (b = 0; b < N; b+=MAX_BUFFER_SZ)
	{
		int NN = XT_MIN(N-b, MAX_BUFFER_SZ);
		pBxl = (      xtfloat *)(bxl);
		pBxr = (      xtfloat *)(bxr);
		pByl = (      xtfloat *)(byl);
		pByr = (      xtfloat *)(byr);

		// Split x into xl (left channel) and xr (right channel).
		for (n = 0; n < (NN); n++)
		{
			XT_LSIP(xl, pX, sizeof(xtfloat));
			XT_LSIP(xr, pX, sizeof(xtfloat));
			XT_SSIP(xl, pBxl, sizeof(xtfloat));
			XT_SSIP(xr, pBxr, sizeof(xtfloat));
		}
		pBxl = (      xtfloat *)bxl;
		pBxr = (      xtfloat *)bxr;

		// Aggregate left channel.
		bkfirf_process(stereo_bkfir->bkfir_left, (float32_t *)pByl, (const float32_t *)pBxl, NN);

		// Aggregate right channel.
		bkfirf_process(stereo_bkfir->bkfir_right, (float32_t *)pByr, (const float32_t *)pBxr, NN);

		// Merge yl (left channel) and yr (right channel) into y.
		for (n = 0; n < (NN); n++)
		{
			XT_LSIP(yl, pByl, sizeof(xtfloat));
			XT_LSIP(yr, pByr, sizeof(xtfloat));
			XT_SSIP(yl, pY, sizeof(xtfloat));
			XT_SSIP(yr, pY, sizeof(xtfloat));
		}
	}
} /* stereo_bkfirf_process() */
#endif
