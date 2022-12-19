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
#include "fft_16x16_stages.h"

#define SWAP_PTR(_x, _y) {complex_fract16 *tmp = _x; _x = _y ; _y = tmp; }


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

extern int stage_inner_DFT4_16x16_ie(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp,  int isStereo);

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
int stereo_fft_cplx16x16_ie(complex_fract16* y, complex_fract16* x, const complex_fract16* twd, int twdstep, int N, int scalingOpt)
{
    ae_int16x4 * restrict px;
    ae_int16x4 * restrict py;
    ae_int32   * restrict ptw;
    int bexp, shift = 0;
    int N4, N8, stridey = 1;
    int i;
    int tw_inc;
    ae_int32x2 t10, t20, t30;
    ae_int16x4 x0, x1, x2, x3;
    ae_int16x4 tw1, tw2, tw3;
    complex_fract16 *pdest = y;

    NASSERT_ALIGN16(twd);
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(x != y);
    NASSERT(scalingOpt == 2);
    NASSERT(N>=8 && (N&(N-1))==0);

    N8 = N >> 3;
    N4 = N >> 2;
    tw_inc = N4 * twdstep*sizeof(complex_fract16);

    /* Calculate the exponent to prescale the input data */
    {
        int n;
        ae_int16x4 x0, x1, x2, x3;
        ae_int16x4 nsa0 = 16, nsa1 = 16;
        px = (ae_int16x4*)x;

        __Pragma("loop_count min=1");
        for (n = 0; n < (N >> 2); n++)
        {
            AE_L16X4X2_IP(x0, x1, castxcc(ae_int16x8, px), sizeof(ae_int16x8));
            AE_L16X4X2_IP(x2, x3, castxcc(ae_int16x8, px), sizeof(ae_int16x8));
            nsa0 = AE_MIN16(nsa0, AE_MIN16(AE_NSA16X4(x0), AE_NSA16X4(x1)));
            nsa1 = AE_MIN16(nsa1, AE_MIN16(AE_NSA16X4(x2), AE_NSA16X4(x3)));
        }
        bexp = AE_RMIN16X4(AE_MIN16(nsa0, nsa1));
    }

    /*
     * First stage of FFT, radix-4
     */
    {
        const int min_shift = 3;
        shift = min_shift - bexp;

        px = (ae_int16x4 *)x;
        py = (ae_int16x4 *)y;
        ptw = (ae_int32 *)twd;
        SetDFT4_Scaling(3);

        __Pragma("loop_count min=1");
        for (i = 0; i < N4; i++)
        {
           
            /* load twiddle factors */
            AE_L32_XP(t10, ptw, tw_inc);
            AE_L32_XP(t20, ptw, tw_inc);
            AE_L32_XP(t30, ptw, twdstep*sizeof(complex_fract16) - 2*tw_inc);
            tw1 = (AE_MOVINT16X4_FROMF32X2(t10));
            tw2 = (AE_MOVINT16X4_FROMF32X2(t20));
            tw3 = (AE_MOVINT16X4_FROMF32X2(t30));

            /* load input data */
            AE_L16X4_XP(x0, px, N4*2*sizeof(complex_fract16));
            AE_L16X4_XP(x1, px, N4*2*sizeof(complex_fract16));
            AE_L16X4_XP(x2, px, N4*2*sizeof(complex_fract16));
            AE_L16X4_XP(x3, px, sizeof(ae_int16x4) - 3*N4*2*sizeof(complex_fract16));

            x0 = AE_SLAA16(x0, bexp);
            x1 = AE_SLAA16(x1, bexp);
            x2 = AE_SLAA16(x2, bexp);
            x3 = AE_SLAA16(x3, bexp);

            /* compute the butterfly */
            DFT4XI2(x0, x1, x2, x3);

            x1 = AE_MULFC16RAS(x1, tw1);
            x2 = AE_MULFC16RAS(x2, tw2);
            x3 = AE_MULFC16RAS(x3, tw3);

            /* store data */
            AE_S16X4RNG_IP(x0, py, sizeof(ae_int16x4));
            AE_S16X4RNG_IP(x1, py, sizeof(ae_int16x4));
            AE_S16X4RNG_IP(x2, py, sizeof(ae_int16x4));
            AE_S16X4RNG_IP(x3, py, sizeof(ae_int16x4));

        }

        {
            int a, b;
            AE_CALCRNG16(a, b, 0, 3); (void)b;
            bexp = 3 - a;
        }

        stridey *= 4;
        SWAP_PTR(x, y);
    }

    /*
     * Next FFT stages except the last, radix-4
     */
    while (stridey < N8)
    {
        twdstep *= 4;
        shift += stage_inner_DFT4_16x16_ie((const int16_t*)twd, (const int16_t*)x, (int16_t*)y, N, &stridey, twdstep, &bexp, 1);
        SWAP_PTR(x, y);
    }

    /*
     * Last FFT stage, radix-2 or radix-4
     */
    if (y != pdest)
    {
        /* Execute the last stage inplace */
        y = x;
    }

    if (stridey == N8)/* radix-8 */
    {
        stridey *= 2; 
        shift += fft_16x16_stage_last_scl2_DFT8(NULL, (const int16_t*)x, (int16_t*)y, N, &stridey, 0, &bexp);
    }
    else/* radix-4 */
    {
        int stage_shift;
        const int min_shift = 2;
    
        stage_shift = XT_MAX(0, min_shift - bexp);
        px = (ae_int16x4 *)x;
        py = (ae_int16x4 *)y;

        SetDFT4_Scaling(stage_shift);
        __Pragma("loop_count min=2");
        for (i = 0; i < N4; i++)
        {
            /* load input data */
            AE_L16X4_XP(x0, px, N4*2*sizeof(complex_fract16));
            AE_L16X4_XP(x1, px, N4*2*sizeof(complex_fract16));
            AE_L16X4_XP(x2, px, N4*2*sizeof(complex_fract16));
            AE_L16X4_XP(x3, px, sizeof(ae_int16x4) - 3*N4*2*sizeof(complex_fract16));

            /* compute the butterfly */
            DFT4XI2(x0, x1, x2, x3);

            /* store data */
            AE_S16X4_XP(x0, py, N4*2*sizeof(complex_fract16));
            AE_S16X4_XP(x1, py, N4*2*sizeof(complex_fract16));
            AE_S16X4_XP(x2, py, N4*2*sizeof(complex_fract16));
            AE_S16X4_XP(x3, py, sizeof(ae_int16x4) - 3*N4*2*sizeof(complex_fract16));
        }
        shift += stage_shift;
    }
    return shift;
} /* stereo_fft_cplx16x16_ie() */
