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
    Real-valued inverse FFT: 32-bit data, 16-bit twiddle factors
    C code optimized for HiFi4
  IntegrIT, 2006-2019
*/

/* Signal Processing Library API. */
#include "NatureDSP_Signal_fft.h"
/* Common utility macros. */
#include "common.h"
/* FFT descriptor structure and twiddle factor tables. */
#include "fft_x16_common.h"
#include "fft_32x16_stages.h"

#define SWAP_PTR(_x, _y) {int32_t *tmp = _x; _x = _y ; _y = tmp; } 

ALIGN(16) static const int16_t sel_tab[4] = { 0x705, 0x604, 0x103, 0x002 };

/*
	in-place inverse split part of FFT:
	x[N+2]  input (N+2 samples)/output(N samples)
	N       size of FFT
*/
static void isplitPart_x2(int32_t *x, const int16_t *twdSplit, int shift, int N)
{
  int shiftl, shiftr;
  int i;
  const int step_back = -8;

        ae_int32x2 * restrict p_x0;
        ae_int32x2 * restrict p_x1;
  const ae_int32x2 * restrict p_x0rd;
  const ae_int32x2 * restrict p_x1rd;
  const ae_int32   * restrict p_twd;

  ae_int32x2 vA0, vA1, vB0, vB1, vC0, vC1, scl;
  ae_int16x4 vT;
  ae_int64 T64_0, T64_1;
  ae_int16x4 tmp0, tmp1;
  ae_int16x4 sel = AE_L16X4_I((ae_int16x4*)sel_tab, 0);

  NASSERT_ALIGN16(x);

  shiftl = XT_MAX(0, -shift);
  shiftr = XT_MAX(0,  shift);
  scl = 1 << shiftl;
  WUR_AE_SAR(shiftr);

  p_twd = (const ae_int32 *)twdSplit+1;
  p_x0 = (ae_int32x2 *)(x);
  p_x1 = (ae_int32x2 *)(x+N);

  // first point
  vA0 = AE_L32X2_I(p_x0, 0);
  vA1 = AE_L32X2_I(p_x1, 0);
  vA0 = AE_MULP32X2(vA0, scl);
  vA1 = AE_MULP32X2(vA1, scl);

  AE_ADDANDSUBRNG32(vB0, vB1, vA0, vA1);
  vB0 = AE_SEL32_HH(vB0, vB1);
  vB1 = AE_MOVI(0);

  AE_S32X2_IP(vB0, p_x0, 8);
  AE_S32X2_XP(vB1, p_x1, step_back);

  p_x0rd = (const ae_int32x2 *)p_x0;
  p_x1rd = (const ae_int32x2 *)p_x1;

  {
    // load next data
    AE_L32X2_IP(vA0, p_x0rd, 8);
    AE_L32X2_XP(vA1, p_x1rd, step_back);
    vA0 = AE_MULP32X2(vA0, scl);
    vA1 = AE_MULP32X2(vA1, scl);

    // load twiddle
    AE_L32_XP(vB1, p_twd, 2*sizeof(int16_t));
    vT = AE_MOVINT16X4_FROMINT32X2(vB1);

    // ADD/SUBB
    AE_ADDANDSUBRNG32(vB0, vB1, vA0, vA1);

    vA0 = AE_SEL32_HL(vB1, vB0);
    vB1 = AE_SEL32_LH(vB1, vB0);

    // do rotation
    vB0 = AE_MULFC32X16RAS_H(vA0, vT);

    vC0 = AE_ADDSUB32S(vB1, vB0);
    vC1 = AE_SUBADD32S(vB0, vB1);

    T64_0 = AE_MOVINT64_FROMINT32X2(vC0);
    T64_1 = AE_MOVINT64_FROMINT32X2(vC1);
    AE_S64_IP(T64_0, castxcc(ae_int64, p_x0), 8);
    AE_S64_XP(T64_1, castxcc(ae_int64, p_x1), step_back);
  }

  NASSERT_ALIGN16(p_x0rd); 
  NASSERT_ALIGN8(p_twd);
  NASSERT( (N&7)==0 ); 
  if (shift < 0)
  {
      __Pragma("loop_count min=1");
      for (i = 2; i < (N >> 2); i += 2)
      {
          /* 6 cycles per pipeline stage in steady state with unroll = 1 */
          ae_int32x2 vA00, vA10, vB00, vB10, vC00, vC10;
          ae_int32x2 vA01, vA11, vB01, vB11, vC01, vC11;
          // load data
          AE_L32X2X2_IP(vA00, vA01, castxcc(ae_int32x4, p_x0rd), 2 * sizeof(complex_fract32));
          AE_L32X2_XP(vA10, p_x1rd, step_back);
          AE_L32X2_XP(vA11, p_x1rd, step_back);
          // load twiddles
          AE_L32X2_XP(vB10, castxcc(ae_int32x2, p_twd), 2 * sizeof(complex_fract16));
          vT = AE_MOVINT16X4_FROMINT32X2(vB10);

          /* Prescaling data */
          vA00 = AE_MULP32X2(vA00, scl);
          vA10 = AE_MULP32X2(vA10, scl);
          vA01 = AE_MULP32X2(vA01, scl);
          vA11 = AE_MULP32X2(vA11, scl);

          // ADD/SUBB
          AE_ADDANDSUBRNG32(vB00, vB10, vA00, vA10);
          AE_ADDANDSUBRNG32(vB01, vB11, vA01, vA11);
#if 0
          vA00 = AE_SEL32_HL(vB10, vB00);
          vB10 = AE_SEL32_LH(vB10, vB00);
          vA01 = AE_SEL32_HL(vB11, vB01);
          vB11 = AE_SEL32_LH(vB11, vB01);
#else
          AE_DSEL16X4(tmp0, tmp1, AE_MOVINT16X4_FROMINT32X2(vB10), AE_MOVINT16X4_FROMINT32X2(vB00), sel);
          vA00 = AE_MOVINT32X2_FROMINT16X4(tmp0);
          vB10 = AE_MOVINT32X2_FROMINT16X4(tmp1);
          AE_DSEL16X4(tmp0, tmp1, AE_MOVINT16X4_FROMINT32X2(vB11), AE_MOVINT16X4_FROMINT32X2(vB01), sel);
          vA01 = AE_MOVINT32X2_FROMINT16X4(tmp0);
          vB11 = AE_MOVINT32X2_FROMINT16X4(tmp1);
#endif
          // do rotation
          vB00 = AE_MULFC32X16RAS_H(vA00, vT);
          vB01 = AE_MULFC32X16RAS_L(vA01, vT);

          vC00 = AE_ADDSUB32S(vB10, vB00);
          vC10 = AE_SUBADD32S(vB00, vB10);
          vC01 = AE_ADDSUB32S(vB11, vB01);
          vC11 = AE_SUBADD32S(vB01, vB11);

          AE_S64X2_IP(AE_MOVINT64_FROMINT32X2(vC00), AE_MOVINT64_FROMINT32X2(vC01), castxcc(ae_int64x2, p_x0), 2 * sizeof(complex_fract32));
          AE_S64_XP(AE_MOVINT64_FROMINT32X2(vC10), castxcc(ae_int64, p_x1), step_back);
          AE_S64_XP(AE_MOVINT64_FROMINT32X2(vC11), castxcc(ae_int64, p_x1), step_back);
      }
  }
  else /* if (scl != 1) */
  {
      __Pragma("loop_count min=1");
      __Pragma("no_unroll");
      for (i = 2; i < (N >> 2); i += 2)
      {
          /* 5 cycles per pipeline stage in steady state with unroll = 1 */
          ae_int32x2 vA00, vA10, vB00, vB10, vC00, vC10;
          ae_int32x2 vA01, vA11, vB01, vB11, vC01, vC11;
          // load data
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
          vA00 = AE_SEL32_HL(vB10, vB00);
          vB10 = AE_SEL32_LH(vB10, vB00);
          vA01 = AE_SEL32_HL(vB11, vB01);
          vB11 = AE_SEL32_LH(vB11, vB01);
#else
          AE_DSEL16X4(tmp0, tmp1, AE_MOVINT16X4_FROMINT32X2(vB10), AE_MOVINT16X4_FROMINT32X2(vB00), sel);
          vA00 = AE_MOVINT32X2_FROMINT16X4(tmp0);
          vB10 = AE_MOVINT32X2_FROMINT16X4(tmp1);
          AE_DSEL16X4(tmp0, tmp1, AE_MOVINT16X4_FROMINT32X2(vB11), AE_MOVINT16X4_FROMINT32X2(vB01), sel);
          vA01 = AE_MOVINT32X2_FROMINT16X4(tmp0);
          vB11 = AE_MOVINT32X2_FROMINT16X4(tmp1);
#endif

          // do rotation
          vB00 = AE_MULFC32X16RAS_H(vA00, vT);
          vB01 = AE_MULFC32X16RAS_L(vA01, vT);

          vC00 = AE_ADDSUB32S(vB10, vB00);
          vC10 = AE_SUBADD32S(vB00, vB10);
          vC01 = AE_ADDSUB32S(vB11, vB01);
          vC11 = AE_SUBADD32S(vB01, vB11);

          AE_S64X2_IP(AE_MOVINT64_FROMINT32X2(vC00), AE_MOVINT64_FROMINT32X2(vC01), castxcc(ae_int64x2, p_x0), 2 * sizeof(complex_fract32));
          AE_S64_XP(AE_MOVINT64_FROMINT32X2(vC10), castxcc(ae_int64, p_x1), step_back);
          AE_S64_XP(AE_MOVINT64_FROMINT32X2(vC11), castxcc(ae_int64, p_x1), step_back);
      }
  } /* if (scl != 1) .. else ..   */

  // middle sample
  vA0 = AE_L32X2_I(p_x0rd, 0);
  vA0 = AE_SRAA32(vA0, shift-1);
  vB0 = AE_NEG32S(vA0);
  vC0 = AE_SEL32_HL(vA0, vB0);
  AE_S32X2_I(vC0, p_x0, 0);
}

/*-------------------------------------------------------------------------
  Inverse FFT on Real Data
  These functions make inverse FFT on half spectral data forming real
  data samples.
      Scaling: 
      +-------------------+----------------------------------------+
      |      Function     |           Scaling options              |
      +-------------------+----------------------------------------+
      |  ifft_real16x16   |  2 - 16-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  ifft_real32x32   |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  ifft_real32x16   |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      +-------------------+----------------------------------------+

  NOTES:
  1. Bit-reversing reordering is done here. 
  2. IFFT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after
     the call.
  3. Inverse FFT function for real signal transforms the input spectrum  
     and then calls ifft_cplx() with FFT size set to N/2. 
  4. 32x32,32x16,16x16  FFTs support mixed radix transforms.
  5. N - FFT size

  Precision:
  32x32  32-bit input/outputs, 32-bit twiddles
  32x16  32-bit input/outputs, 16-bit twiddles
  16x16  16-bit input/outputs, 16-bit twiddles

  Input:
  x[(N/2+1)*2]	input spectrum. Real and imaginary data are interleaved  
                and real data goes first. The imaginary part of 0th and
                N/2th input samples should be equal to zero.
  scalingOpt	scaling option (see table above)

  Output:			
  y[N]	        real output signal

  Returned value: total number of right shifts occurred during scaling 
                  procedure

  Restrictions:
  x,y           should not overlap
  x,y           aligned on a 16-bytes boundary
  x[(0)*2+1],
  x[(N/2)*2+1]  should be equal to zero
-------------------------------------------------------------------------*/
int ifft_real32x16( int32_t * y, int32_t * x, fft_handle_t h, int scalingOpt )
{
    int bexp, shift;
    int N;
    int s = 0;
    int v = 1;
    const int isplitScale = 2; // Scaling of isplitPart_x2
    int32_t *pdest = y; 
    const int16_t * twdSplit;
    const fft_real_x16_descr_t * h_real;
    const fft_cplx_x16_descr_t * h_cplx;
    const cint16ptr_t_fft *tw_tab;
    const int *tw_step;
    const ae_int32x4 * restrict pX;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(scalingOpt==3 || scalingOpt==2);

    h_real = (const fft_real_x16_descr_t*)h;
    h_cplx = h_real->cfft_hdl;
    twdSplit = h_real->twd;
    tw_tab = h_cplx->twd;
    tw_step = h_cplx->tw_step;
    const fft_cplx32x16_stage_t*fnstg_Tbl = (const fft_cplx32x16_stage_t*)((scalingOpt == 2) ? h_cplx->fnstages_32x16_s2 : h_cplx->fnstages_32x16_s3);
    N=(h_cplx->N) << 1;

    bexp = 0;
    if (scalingOpt == 2)
    {
        ae_int32x2 x0, x1, x2, x3;
        ae_int16x4 nsa0, nsa1;
        int n;
        pX = (const ae_int32x4 *)x;
        nsa0 = 31; nsa1 = 31;
        NASSERT((N & 3) == 0);
        __Pragma("loop_count min=1");
        for (n = 0; n<(N >> 3); n++)
        {
            AE_L32X2X2_IP(x0, x1, pX, sizeof(ae_int32x4));
            AE_L32X2X2_IP(x2, x3, pX, sizeof(ae_int32x4));
            nsa0 = AE_MIN16(nsa0, AE_NSA32X4(x0, x1));
            nsa1 = AE_MIN16(nsa1, AE_NSA32X4(x2, x3));
        }
        x0 = x1 = AE_L32X2_I((ae_int32x2*)pX, 0);
        nsa0 = AE_MIN16(nsa0, AE_NSA32X4(x0, x1));
        bexp = AE_RMIN16X4(AE_MIN16(nsa0, nsa1));
    }
    shift = isplitScale - bexp;
    isplitPart_x2(x, twdSplit, shift, N);

    /* half-sized complex IFFT */
    N >>= 1;
    bexp = 0;
    while ( tw_tab[s] != NULL )
    {
#if 1
        {
          NASSERT(fnstg_Tbl[0]!=NULL);
          shift += fnstg_Tbl[0](tw_tab[s], x, y, N, &v, tw_step[s], &bexp);
          fnstg_Tbl++;
        }
#else
        {
            int stg_type;
            fft_cplx32x16_stage_t stg_fn;
            stg_type = (int)stg_typeTbl[s];
            NASSERT(stg_type < NUM_FFT_STAGE_TYPES);
            stg_fn = stg_Tbl[stg_type];
            NASSERT(stg_fn != NULL);
            shift += stg_fn(tw_tab[s], x, y, N, &v, tw_step[s], &bexp); 
        }
#endif
        SWAP_PTR(x, y); 
        s++; 
    }

    if (y != pdest)
    {  
        /* Execute the last stage inplace */
        y = x;
    }

    /* Last stage */
#if 1
    {
        NASSERT(fnstg_Tbl[0]!=NULL);
        shift += fnstg_Tbl[0](tw_tab[s], x, y, N, &v, tw_step[s], &bexp);
        fnstg_Tbl++;
    }
#else
    {
        int stg_type;
        fft_cplx32x16_stage_t stg_fn;
        stg_type = (int)stg_typeTbl[s];
        NASSERT(stg_type < NUM_FFT_STAGE_TYPES);
        stg_fn = stg_Tbl[stg_type];
        NASSERT(stg_fn != NULL);
        shift += stg_fn(tw_tab[s], x, y, N, &v, tw_step[s], &bexp);
    }
#endif
    return shift;
} /* ifft_real32x16() */
