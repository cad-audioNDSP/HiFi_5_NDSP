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
    C code optimized for HiFi4
    Integrit, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fft.h"
/* Common utility and macros declarations. */
#include "common.h"

/*-------------------------------------------------------------------------
  Internal stages of functions stereo_fft_cplx32x16_ie(),
  stereo_ifft_cplx32x16_ie()

  Performs all stages of FFT except last one. It is assumed that
  the last stage is radix-2/radix-4 and implemented with
  bit-reversal permutation.
 
  Input/Output:
  x[N]                complex input signal. Real and imaginary data 
                      are interleaved and real data goes first

  Input:
  twd[N*twdstep*3/4]  twiddle factor table of a complex-valued FFT of 
                      size N*twdstep
  N                   FFT size
  twdstep             twiddle step 
  scalingOpt          scaling option

  Returned value:     total number of right shifts occurred during 
                      scaling procedure

  Restrictions:
  x - should not overlap and must be aligned on 16-bytes boundary
-------------------------------------------------------------------------*/
int stereo_fft_cplx32x16_ie_inner(complex_fract32* x, const complex_fract16* twd, int twdstep, int N, int scalingOpt);

/*-------------------------------------------------------------------------
  FFT on Complex Data with Optimized Memory Usage
  These functions make FFT on complex data with optimized memory usage.
  Scaling: 
      +-------------------+----------------------------------------+
      |      Function     |           Scaling options              |
      +-------------------+----------------------------------------+
      |  fft_cplx16x16_ie |  2 - 16-bit dynamic scaling            | 
      |  fft_cplx32x16_ie |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  fft_cplx32x32_ie |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      +-------------------+----------------------------------------+
  NOTES:
  1. Bit-reversing reordering is done here.
  2. FFT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after 
     the call.
  3. FFT of size N may be supplied with constant data
     (twiddle factors) of a larger-sized FFT = N*twdstep.
  4. Stereo FFTs accept inputs/form outputs in the interleaved order:
     left complex sample, right complex sample

  Precision: 
  16x16_ie      16-bit input/outputs, 16-bit twiddles
  32x16_ie      32-bit input/outputs, 16-bit twiddles
  32x32_ie      32-bit input/outputs, 32-bit twiddles
  f_ie          floating point
 
  Input:
  S                     1 for ordinary (single channel) FFT, 2 - for stereo
                        input/outputs
  x[N*S]                complex input signal. Real and imaginary data 
                        are interleaved and real data goes first
  twd[N*twdstep*3/4]    twiddle factor table of a complex-valued FFT of 
                        size N*twdstep
  N                     FFT size
  twdstep               twiddle step 
  scalingOpt            scaling option (see table above), not applicable
                        to the floating point function 
  Output:
  y[N*S]                output spectrum. Real and imaginary data are 
                        interleaved and real data goes first

  Returned value: total number of right shifts occurred during scaling 
                  procedure. Floating point function always return 0.

  Restrictions:
  x,y   should not overlap
  x,y   aligned on 16-bytes boundary
-------------------------------------------------------------------------*/
int stereo_fft_cplx32x16_ie(complex_fract32* y,complex_fract32* x, const complex_fract16* twd, int twdstep, int N, int scalingOpt)
{
    int shift;
    int32_t i,i0,ai;
    ae_int32x4 * restrict p_y0l = (ae_int32x4 *)(y);
    ae_int32x4 * restrict p_y1l = p_y0l + (N>>2);
    ae_int32x4 * restrict p_y2l = p_y1l + (N>>2);
    ae_int32x4 * restrict p_y3l = p_y2l + (N>>2);
    ae_int32x4 * restrict p_xl = (ae_int32x4 *)(x);
    ae_int32x2 vA00, vA10, vA20, vA30;
    ae_int32x2 vA01, vA11, vA21, vA31;
    ae_int32x2 vB00, vB10, vB20, vB30;
    ae_int32x2 vB01, vB11, vB21, vB31;

    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(scalingOpt==3 || scalingOpt==2);
    NASSERT(N>=8 && (N&(N-1))==0);

    /* Perform all stages of FFT except last one */
    shift = stereo_fft_cplx32x16_ie_inner(x, twd, twdstep, N, scalingOpt);

    /* Last stage */
    i = NSA(N);
    ai=((int32_t)0x1)<<i;
    i0=0;

    if ( (i&1)!=0 )
    {
        shift += 1;// Set scaling
        WUR_AE_SAR(1);
        //--------------------------------------------------------------------------
        // last stage is RADIX2 !!!
        //--------------------------------------------------------------------------
        __Pragma("loop_count min=2");
        for (i = 0; i < (N >> 2); i++)
        {
            AE_L32X2X2_IP(vA00, vA01, p_xl, 2 * sizeof(ae_int32x2));
            AE_L32X2X2_IP(vA10, vA11, p_xl, 2 * sizeof(ae_int32x2));
            AE_L32X2X2_IP(vA20, vA21, p_xl, 2 * sizeof(ae_int32x2));
            AE_L32X2X2_IP(vA30, vA31, p_xl, 2 * sizeof(ae_int32x2));

            /* left channel */
            AE_ADDANDSUBRNG32(vB00, vB20, vA00, vA10);
            AE_ADDANDSUBRNG32(vB10, vB30, vA20, vA30);
            /* right channel */
            AE_ADDANDSUBRNG32(vB01, vB21, vA01, vA11);
            AE_ADDANDSUBRNG32(vB11, vB31, vA21, vA31);

            AE_S32X2X2_X(vB00, vB01, (ae_int32x4*)p_y0l, i0);
            AE_S32X2X2_X(vB10, vB11, (ae_int32x4*)p_y1l, i0);
            AE_S32X2X2_X(vB20, vB21, (ae_int32x4*)p_y2l, i0);
            AE_S32X2X2_X(vB30, vB31, (ae_int32x4*)p_y3l, i0);

            i0 = AE_ADDBRBA32(i0, ai);
        }
    }
    else
    {
        shift += 2;// Set scaling
        WUR_AE_SAR(2);
        //--------------------------------------------------------------------------
        // last stage is RADIX4 !!!
        //--------------------------------------------------------------------------
        __Pragma("loop_count min=2"); 
        for (i = 0; i<(N>>2); i++)
        {
            AE_L32X2X2_IP(vA00, vA01, p_xl, 2 * sizeof(ae_int32x2));
            AE_L32X2X2_IP(vA10, vA11, p_xl, 2 * sizeof(ae_int32x2));
            AE_L32X2X2_IP(vA20, vA21, p_xl, 2 * sizeof(ae_int32x2));
            AE_L32X2X2_IP(vA30, vA31, p_xl, 2 * sizeof(ae_int32x2));

            /* left channel */
            AE_ADDANDSUBRNG32(vB00, vB20, vA00, vA20);
            AE_ADDANDSUBRNG32(vB10, vB30, vA10, vA30);
            vB30 = AE_MUL32JS(vB30);
            AE_ADDANDSUB32S(vA00, vA20, vB00, vB10);
            AE_ADDANDSUB32S(vA30, vA10, vB20, vB30);

            /* right channel */
            AE_ADDANDSUBRNG32(vB01, vB21, vA01, vA21);
            AE_ADDANDSUBRNG32(vB11, vB31, vA11, vA31);
            vB31 = AE_MUL32JS(vB31);
            AE_ADDANDSUB32S(vA01, vA21, vB01, vB11);
            AE_ADDANDSUB32S(vA31, vA11, vB21, vB31);

            AE_S32X2X2_X(vA00, vA01, (ae_int32x4*)p_y0l, i0);
            AE_S32X2X2_X(vA10, vA11, (ae_int32x4*)p_y1l, i0);
            AE_S32X2X2_X(vA20, vA21, (ae_int32x4*)p_y2l, i0);
            AE_S32X2X2_X(vA30, vA31, (ae_int32x4*)p_y3l, i0);
            i0 = AE_ADDBRBA32(i0, ai);
        }
    }

    return shift;
} /* stereo_fft_cplx32x16_ie() */
