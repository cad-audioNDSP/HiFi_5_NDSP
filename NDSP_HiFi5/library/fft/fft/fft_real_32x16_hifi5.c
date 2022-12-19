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
    Real-valued forward FFT: 32-bit data, 16-bit twiddle factors
    C code optimized for HiFi4
  IntegrIT, 2006-2019
*/

/* Signal Processing Library API. */
#include "NatureDSP_Signal_fft.h"
/* Common utility macros. */
#include "common.h"
/* FFT descriptor structure and twiddle factor tables. */
#include "fft_x16_common.h"

ALIGN(16) static const int16_t sel_tab[4] = { 0x705, 0x604, 0x103, 0x002 };

/*
	in-place split part of FFT:
	x[N+2]	input/output
	N		size of FFT
*/
static void splitPart_x2(int32_t *x, const int16_t *twdSplit, int shift, int N)
{
  int i;
  const int step_back = -8;

        ae_int32x2 * restrict p_x0;
        ae_int32x2 * restrict p_x1;
  const ae_int32x2 * restrict p_x0rd;
  const ae_int32x2 * restrict p_x1rd;
  const ae_int32 * restrict p_twd;
  ae_int32x2  vA0, vA1, vB0, vB1, vC0, vC1;
  ae_int16x4  vT;
  ae_int16x4 sel = AE_L16X4_I((ae_int16x4*)sel_tab, 0);

  NASSERT_ALIGN16(x);
  NASSERT((N&7)==0); 

  p_twd = (const ae_int32 *)twdSplit+1;
  p_x0 = (ae_int32x2 *)(x);
  p_x1 = (ae_int32x2 *)(x+N);

  WUR_AE_SAR(1+shift);
  // load data and prepare pointers for pre-increment
  // first and last samples
  vA0 = AE_L32X2_I(p_x0, 0);
  vA0 = AE_SRAA32S(vA0, shift);
  vA1 = AE_SEL32_LH(vA0, vA0);
  AE_ADDANDSUB32S(vB0, vB1, vA0, vA1);
  vA1 = AE_MOVI(0);
  vB0 = AE_SEL32_HH(vB0, vA1);
  vB1 = AE_SEL32_HH(vB1, vA1);
  AE_S32X2_IP(vB0, p_x0, 8);
  AE_S32X2_XP(vB1, p_x1, step_back);

  p_x0rd = (const ae_int32x2 *)p_x0;
  p_x1rd = (const ae_int32x2 *)p_x1;
  {
      // load next data
      AE_L32X2_IP(vA0, p_x0rd, 8);
      AE_L32X2_XP(vA1, p_x1rd, step_back);

      // load twiddle
      AE_L32_XP(vB1, p_twd, 2 * sizeof(int16_t));
      vT = AE_MOVINT16X4_FROMINT32X2(vB1);

      // ADD/SUBB
      AE_ADDANDSUBRNG32(vB0, vB1, vA0, vA1);

      vA0 = AE_SEL32_LH(vB0, vB1);
      vB1 = AE_SEL32_HL(vB0, vB1);

      // do rotation
      vB0 = AE_MULFC32X16RAS_H(vA0, vT);

      // SUM/DIFF
      vC0 = AE_ADDSUB32S(vB1, vB0);
      vC1 = AE_SUBADD32S(vB0, vB1);
      vC1 = AE_NEG32S(vC1);

      AE_S32X2_IP(vC0, p_x0, 8);
      AE_S32X2_XP(vC1, p_x1, step_back);
  }

    __Pragma("no_unroll"); 
    __Pragma( "loop_count min=1" );
    for (i = 2; i < (N>>2); i+=2) 
    {
        /* 5 cycles with unroll=1 */
        ae_int16x4 tmp0, tmp1;
        ae_int32x2 vA00, vA10, vB00, vB10, vC00, vC10;
        ae_int32x2 vA01, vA11, vB01, vB11, vC01, vC11;

        // load next data
        AE_L32X2X2_IP(vA00, vA01, castxcc(ae_int32x4, p_x0rd), 2 * sizeof(complex_fract32));

        AE_L32X2_XP(vA10, p_x1rd, step_back);
        AE_L32X2_XP(vA11, p_x1rd, step_back);

        // load twiddles
        AE_L32X2_XP(vB10, castxcc(ae_int32x2, p_twd), 2 * sizeof(complex_fract16));
        vT = AE_MOVINT16X4_FROMINT32X2(vB10);

        // ADD/SUBB
        AE_ADDANDSUBRNG32(vB00, vB10, vA00, vA10);
        AE_ADDANDSUBRNG32(vB01, vB11, vA01, vA11);
#if 0
        vA00 = AE_SEL32_LH(vB00, vB10);
        vB10 = AE_SEL32_HL(vB00, vB10);
        vA01 = AE_SEL32_LH(vB01, vB11);
        vB11 = AE_SEL32_HL(vB01, vB11);
#else
        AE_DSEL16X4(tmp0, tmp1, AE_MOVINT16X4_FROMINT32X2(vB00), AE_MOVINT16X4_FROMINT32X2(vB10), sel);
        vB10 = AE_MOVINT32X2_FROMINT16X4(tmp0);
        vA00 = AE_MOVINT32X2_FROMINT16X4(tmp1);
        AE_DSEL16X4(tmp0, tmp1, AE_MOVINT16X4_FROMINT32X2(vB01), AE_MOVINT16X4_FROMINT32X2(vB11), sel);
        vB11 = AE_MOVINT32X2_FROMINT16X4(tmp0);
        vA01 = AE_MOVINT32X2_FROMINT16X4(tmp1);
#endif
        // do rotation
        vB00 = AE_MULFC32X16RAS_H(vA00, vT);
        vB01 = AE_MULFC32X16RAS_L(vA01, vT);

        // SUM/DIFF
        vC00 = AE_ADDSUB32S(vB10, vB00);
        vC10 = AE_SUBADD32S(vB00, vB10);
        vC10 = AE_NEG32S(vC10);
        vC01 = AE_ADDSUB32S(vB11, vB01);
        vC11 = AE_SUBADD32S(vB01, vB11);
        vC11 = AE_NEG32S(vC11);

        AE_S32X2X2_IP(vC00, vC01, castxcc(ae_int32x4, p_x0), 2 * sizeof(complex_fract32));
        AE_S32X2_XP(vC10, p_x1, step_back);
        AE_S32X2_XP(vC11, p_x1, step_back);
    }

  // middle sample
  vA0 = AE_L32X2_I(p_x0, 0);
  vA0 = AE_SRAA32S(vA0, shift);
  vB0 = AE_NEG32S(vA0);
  vC0 = AE_SEL32_HL(vA0, vB0);
  AE_S32X2_I(vC0, p_x0, 0);
}

/*-------------------------------------------------------------------------
  FFT on Real Data
  These functions make FFT on real data forming half of spectrum
      Scaling: 
      +-------------------+----------------------------------------+
      |      Function     |           Scaling options              |
      +-------------------+----------------------------------------+
      |  fft_real16x16    |  2 - 16-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  fft_real32x32    |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  fft_real32x16    |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      +-------------------+----------------------------------------+
  NOTES:
  1. Bit-reversal reordering is done here. 
  2. FFT runs in-place so INPUT DATA WILL APPEAR DAMAGED after the call.
  3. Real data FFT function calls fft_cplx() to apply complex FFT of size
     N/2 to input data and then transforms the resulting spectrum.
  4. 32x32,32x16,16x16 FFTs support mixed radix transforms 
  5. N - FFT size

  Precision:
  32x32  32-bit input/outputs, 32-bit twiddles
  32x16  32-bit input/outputs, 16-bit twiddles
  16x16  16-bit input/outputs, 16-bit twiddles

  Input:
  x[N]          input signal
  scalingOpt    scaling option (see table above)
  Output:
  y[(N/2+1)*2]  output spectrum (positive side)

  Returned value: total number of right shifts occurred during scaling 
                  procedure

  Restrictions:
  x,y           should not overlap
  x,y           aligned on a 16-bytes boundary
-------------------------------------------------------------------------*/
int fft_real32x16( int32_t * y, int32_t * x, fft_handle_t h, int scalingOpt )
{
  const int16_t * twdSplit;
  const fft_real_x16_descr_t * h_real;
  const fft_cplx_x16_descr_t * h_cplx;
  int shift, scale=0;
  int N;

  NASSERT_ALIGN16(x);
  NASSERT_ALIGN16(y);
  NASSERT(scalingOpt==3 || scalingOpt==2);

  h_real = (const fft_real_x16_descr_t*)h;
  h_cplx = h_real->cfft_hdl;
  twdSplit = h_real->twd;
  N=(h_cplx->N) << 1;

  scale=fft_cplx32x16(y, x, (fft_handle_t)h_cplx, scalingOpt);/* fft of half-size */

  shift = (scalingOpt==2 || (N&(N-1))!=0) ? 1 : 0;

  splitPart_x2(y, twdSplit, shift, N);

  return scale+shift;
} // fft_real32x16()
