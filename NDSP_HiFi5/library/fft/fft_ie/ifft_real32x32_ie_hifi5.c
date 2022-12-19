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
#include "NatureDSP_Signal_vector.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "fft_twiddles32x32.h"

#define SWAP_PTR(_x, _y) {complex_fract32 *tmp = _x; _x = _y ; _y = tmp; } 
#define _CONJ32(_x) {_x = AE_SEL32_HL(_x, AE_NEG32S(_x) ); }

/*
The real - to - complex spectrum conversion
MATLAB code:
twd = exp(-2*pi*1j*(0:N/4-1)/N);
a0 = x(1:N/4);
a1 = wrev(x(N/4+2:N/2+1));
b0 = a0+conj(a1);
b1 = (a0-conj(a1))*1j.*conj(twd);
a0 = b0+b1;
a1 = conj(b0-b1);
x = [a0,2*conj(x(N/4+1)),wrev(a1(2:N/4))]; % N/2 complex samples
*/
static int iSpectrConv(complex_fract32 *x, int N, const complex_fract32 *twiddle_table, int twiddle_stride, int scalingOpt, int *bexp)
{
    ae_int32x2 vA0, vA1, vB0, vB1, tw;
    ae_int32x2 * restrict p_x0,
               * restrict p_x1,
               * restrict ptw = (ae_int32x2*)(twiddle_table + twiddle_stride),
               * restrict p_y0,
               * restrict p_y1;
    int n;
    int shift = /*(scalingOpt == 3) ? 2 : */2 - *bexp;
    ae_int32x2 scl;    
    ASSERT(shift>-32 && shift<32);
    int shiftl, shiftr;
    int  N4 = (N + 2) >> 2; 
    ae_int16x4 tmp0, tmp1;
    ALIGN(16) const int16_t sel_tab[4] = { 0x703, 0x602, 0x105, 0x004 };
    ae_int16x4 sel = AE_L16X4_I((ae_int16x4*)sel_tab, 0);

    shiftl = XT_MAX(0, -shift);
    shiftr = XT_MAX(0, shift);
    scl = 1 << shiftl;
    WUR_AE_SAR(shiftr);

    p_x0 = (ae_int32x2 *)(x);
    p_x1 = (ae_int32x2 *)(x + N / 2);
    p_y0 = (ae_int32x2 *)(x);
    p_y1 = (ae_int32x2 *)(x + N / 2 - 1);

    AE_L32X2_IP(vB0, p_x0, sizeof(complex_fract32));
    AE_L32X2_XP(vB1, p_x1, -(int)sizeof(complex_fract32));

    vB0 = AE_SRAA32(vB0, shift);
    vB1 = AE_SRAA32(vB1, shift);
    AE_ADDANDSUB32S(vA0, vA1, vB0, vB1);

    vA0 = AE_SEL32_HH(vA0, vA1);
    AE_S32X2RNG_IP(vA0, p_y0, sizeof(complex_fract32));

    {
        //n=1 for (n = 1; n < N4; n++)
        {
            AE_L32X2_XP(tw, ptw, twiddle_stride * sizeof(complex_fract32));

            AE_L32X2_IP(vB0, p_x0, sizeof(complex_fract32));
            AE_L32X2_XP(vB1, p_x1, -(int)sizeof(complex_fract32));

            vB0 = AE_MULP32X2(vB0, scl);
            vB1 = AE_MULP32X2(vB1, scl);

            // ADD/SUBB
            AE_ADDANDSUBRNG32(vA0, vA1, vB0, vB1);

            /* vB0 = AE_SEL32_HL(vA0, vA1);
            vB1 = AE_SEL32_HL(vA1, vA0); */

            AE_DSEL16X4(tmp0, tmp1,
                AE_MOVINT16X4_FROMINT32X2(vA0),
                AE_MOVINT16X4_FROMINT32X2(vA1), sel);
            vB0 = AE_MOVINT32X2_FROMINT16X4(tmp0);
            vB1 = AE_MOVINT32X2_FROMINT16X4(tmp1);

            tw = AE_MUL32JS(tw);
            vB1 = AE_MULFCJ32RAS(tw, vB1);

            vA0 = AE_SUBADD32S(vB0, vB1);
            vA1 = AE_ADDSUB32S(vB1, vB0);

            AE_S32X2RNG_IP(vA0, p_y0, sizeof(complex_fract32));
            AE_S32X2RNG_XP(vA1, p_y1, -(int)sizeof(complex_fract32));
        }
        if (scalingOpt == 2)
        {
            __Pragma("loop_count min=1");
            for (n = 2; n < N4; n += 2)
            {
                /* 15 cycles unroll = 2 */
                ae_int32x2 tw0, tw1;
                ae_int32x2 vA00, vA10, vA01, vA11, vB00, vB01, vB11, vB10;
                AE_L32X2X2_IP(vB00, vB01, castxcc(ae_int32x4, p_x0), 2 * sizeof(complex_fract32));
                AE_L32X2_XP(tw0, ptw, twiddle_stride * sizeof(complex_fract32));
                AE_L32X2_XP(tw1, ptw, twiddle_stride * sizeof(complex_fract32));

                AE_L32X2_XP(vB10, p_x1, -(int)sizeof(complex_fract32));
                AE_L32X2_XP(vB11, p_x1, -(int)sizeof(complex_fract32));

                vB00 = AE_MULP32X2(vB00, scl);
                vB10 = AE_MULP32X2(vB10, scl);
                vB01 = AE_MULP32X2(vB01, scl);
                vB11 = AE_MULP32X2(vB11, scl);

                AE_ADDANDSUBRNG32(vA00, vA10, vB00, vB10);
                AE_ADDANDSUBRNG32(vA01, vA11, vB01, vB11);

                AE_DSEL16X4(tmp0, tmp1,
                    AE_MOVINT16X4_FROMINT32X2(vA00),
                    AE_MOVINT16X4_FROMINT32X2(vA10), sel);
                vB00 = AE_MOVINT32X2_FROMINT16X4(tmp0);
                vB10 = AE_MOVINT32X2_FROMINT16X4(tmp1);

                AE_DSEL16X4(tmp0, tmp1,
                    AE_MOVINT16X4_FROMINT32X2(vA01),
                    AE_MOVINT16X4_FROMINT32X2(vA11), sel);
                vB01 = AE_MOVINT32X2_FROMINT16X4(tmp0);
                vB11 = AE_MOVINT32X2_FROMINT16X4(tmp1);

                tw0 = AE_MUL32JS(tw0);
                tw1 = AE_MUL32JS(tw1);

                vB10 = AE_MULFCJ32RAS(tw0, vB10);
                vB11 = AE_MULFCJ32RAS(tw1, vB11);

                vA00 = AE_SUBADD32S(vB00, vB10);
                vA10 = AE_ADDSUB32S(vB10, vB00);
                vA01 = AE_SUBADD32S(vB01, vB11);
                vA11 = AE_ADDSUB32S(vB11, vB01);

                AE_S32X2RNG_IP(vA00, p_y0, sizeof(complex_fract32));
                AE_S32X2RNG_IP(vA01, p_y0, sizeof(complex_fract32));
                AE_S32X2RNG_XP(vA10, p_y1, -(int)sizeof(complex_fract32));
                AE_S32X2RNG_XP(vA11, p_y1, -(int)sizeof(complex_fract32));
            }
        }
        else /* if (scalingOpt == 2) */
        {

            __Pragma("loop_count min=1");
            for (n = 2; n < N4; n += 2)
            {
                /* 6 cycles unroll = 1 */
                ae_int32x2 tw0, tw1;
                ae_int32x2 vA00, vA10, vA01, vA11, vB00, vB01, vB11, vB10;
                AE_L32X2X2_IP(vB00, vB01, castxcc(ae_int32x4, p_x0), 2 * sizeof(complex_fract32));
                AE_L32X2_XP(tw0, ptw, twiddle_stride * sizeof(complex_fract32));
                AE_L32X2_XP(tw1, ptw, twiddle_stride * sizeof(complex_fract32));

                AE_L32X2_XP(vB10, p_x1, -(int)sizeof(complex_fract32));
                AE_L32X2_XP(vB11, p_x1, -(int)sizeof(complex_fract32));

                AE_ADDANDSUBRNG32(vA00, vA10, vB00, vB10);
                AE_ADDANDSUBRNG32(vA01, vA11, vB01, vB11);

                /*
                vB00 = AE_SEL32_HL(vA00, vA10);
                vB10 = AE_SEL32_HL(vA10, vA00);
                vB01 = AE_SEL32_HL(vA01, vA11);
                vB11 = AE_SEL32_HL(vA11, vA01);*/
                AE_DSEL16X4(tmp0, tmp1,
                    AE_MOVINT16X4_FROMINT32X2(vA00),
                    AE_MOVINT16X4_FROMINT32X2(vA10), sel);
                vB00 = AE_MOVINT32X2_FROMINT16X4(tmp0);
                vB10 = AE_MOVINT32X2_FROMINT16X4(tmp1);

                AE_DSEL16X4(tmp0, tmp1,
                    AE_MOVINT16X4_FROMINT32X2(vA01),
                    AE_MOVINT16X4_FROMINT32X2(vA11), sel);
                vB01 = AE_MOVINT32X2_FROMINT16X4(tmp0);
                vB11 = AE_MOVINT32X2_FROMINT16X4(tmp1);

                tw0 = AE_MUL32JS(tw0);
                tw1 = AE_MUL32JS(tw1);

                vB10 = AE_MULFCJ32RAS(tw0, vB10);
                vB11 = AE_MULFCJ32RAS(tw1, vB11);

                vA00 = AE_SUBADD32S(vB00, vB10);
                vA10 = AE_ADDSUB32S(vB10, vB00);
                vA01 = AE_SUBADD32S(vB01, vB11);
                vA11 = AE_ADDSUB32S(vB11, vB01);

                AE_S32X2X2RNG_IP(vA00, vA01, castxcc(ae_int32x4, p_y0), 2 * sizeof(complex_fract32));

                AE_S32X2RNG_XP(vA10, p_y1, -(int)sizeof(complex_fract32));
                AE_S32X2RNG_XP(vA11, p_y1, -(int)sizeof(complex_fract32));
            }
        } /* if (scalingOpt == 2) .. else ... */
    } /* if ((N4 & 1) != 0 || (N & 3) != 0) .. else ...  */


    /* 2*conj(x(N/4+1)) */
    vB0 = AE_L32X2_I(p_x0, 0);
    vB0 = AE_SRAA32(vB0, shift - 1);
    _CONJ32(vB0);
    AE_S32X2RNG_I(vB0, p_y0, 0);

    *bexp = 0;
    if (scalingOpt == 2)
    {
        AE_CALCRNG3();
        *bexp = 3 - RUR_AE_SAR();
    }

    return shift;
} /* iSpectrConv */

/* Local version of the ifft_cplx32x32_ie with extenal bexp */
static int __ifft_cplx32x32_ie(complex_fract32* y, complex_fract32* x, const complex_fract32* twd, int twdstep, int N, int scalingOpt, int _bexp)
{
    int bexp = 0;
    int v = 1;
    int shift = 0;

    complex_fract32 *pdest = y;
    int log2N = 30 - NSA(N);

    const fft_cplx32x32_stage_t first_stg_fn = (scalingOpt == 2) ? ifft_stageS2_DFT4_first_32x32 : ifft_stageS3_DFT4_first_32x32;
    const fft_cplx32x32_stage_t stg_fn = (scalingOpt == 2) ? fft_stageS2_DFT4x4_32x32 : fft_stageS3_DFT4x4_32x32;
    const fft_cplx32x32_stage_t last_stg_fn = (log2N & 1) ?
        ((scalingOpt == 2) ? ifft_stageS2_DFT8_last_32x32 : ifft_stageS3_DFT8_last_32x32) :
        ((scalingOpt == 2) ? ifft_stageS2_DFT4_last_32x32 : ifft_stageS3_DFT4_last_32x32);
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT_ALIGN16(twd);
    NASSERT((N&(N - 1)) == 0);
    NASSERT(x != y);
    NASSERT(scalingOpt == 2 || scalingOpt == 3);

    if (scalingOpt == 2)
    {
        bexp = _bexp;
    }

    shift += first_stg_fn((const int32_t*)twd, (int32_t*)x, (int32_t*)y, N, &v, twdstep, &bexp);
    SWAP_PTR(x, y);
    log2N -= 2;
    twdstep *= 4;

    while (log2N >= 4)
    {
        shift += stg_fn((const int32_t*)twd, (int32_t*)x, (int32_t*)y, N, &v, twdstep, &bexp);
        SWAP_PTR(x, y);
        log2N -= 2;
        twdstep *= 4;
    }

    if (y != pdest)
    {
        /* Execute the last stage inplace */
        y = x;
    }

    /* Last stage */
    shift += last_stg_fn(NULL, (int32_t*)x, (int32_t*)y, N, &v, 0, &bexp);
    return shift;
} /* __ifft_cplx32x32_ie() */

/*-------------------------------------------------------------------------
  Inverse FFT on Real Data with Optimized Memory Usage
  These functions make inverse FFT on real data from half of spectrum with
  optimized memory usage.
  Scaling: 
      +-----------------------+--------------------------------------+
      |      Function         |           Scaling options            |
      +-----------------------+--------------------------------------+
      | ifft_real16x16_ie     |  2 - 16-bit dynamic scaling          |
      | ifft_real32x16_ie     |  2 - 32-bit dynamic scaling          |
      |                       |  3 - fixed scaling before each stage |
      | ifft_real32x32_ie     |  2 - 32-bit dynamic scaling          |
      |                       |  3 - fixed scaling before each stage |
      +-----------------------+--------------------------------------+
  NOTES:
  1. Bit-reversing reordering is done here.
  2. INPUT DATA MAY APPEAR DAMAGED after the call.
  3. FFT functions may use input and output buffers for temporal storage
     of intermediate 32-bit data, so FFT functions with 24-bit packed
     I/O (Nx3-byte data) require that the buffers are large enough to 
     keep Nx4-byte data.
  4. FFT of size N may be supplied with constant data (twiddle factors) 
     of a larger-sized FFT = N*twdstep.

  Precision:
  16x16_ie      16-bit input/outputs, 16-bit data, 16-bit twiddles
  32x16_ie      32-bit input/outputs, 32-bit data, 16-bit twiddles
  32x32_ie      32-bit input/outputs, 32-bit data, 32-bit twiddles
  f_ie          floating point

  Input:
  x             input spectrum (positive side). Real and imaginary
                data are interleaved and real data goes first. 
                The imaginary part of 0th and N/2th input samples
                should be equal to zero:
  --------------+----------+-----------------+----------------
  Function      |   Size   |  Allocated Size |       type    |
  --------------+----------+-----------------+----------------
  16x16_ie      |   N/2+1  |      N/2+1      |complex_fract16|
  32x16_ie      |   N/2+1  |      N/2+1      |complex_fract32|
  32x32_ie      |   N/2+1  |      N/2+1      |complex_fract32|
  f_ie          |   N/2+1  |      N/2+1      | complex_float |
  --------------+----------+-----------------+----------------

  twd[2*N*twdstep*3/4]  twiddle factor table of a complex-valued FFT
                        of size N*twdstep
  N                     FFT size
  twdstep               twiddle step
  scalingOpt            scaling option (see table above), not applicable 
                        to the floating point function
  Output:
  y                     output spectrum. Real and imaginary data are 
                        interleaved and real data goes first:
  --------------+----------+-----------------+-----------
  Function      |   Size   |  Allocated Size |  type    |
  --------------+----------+-----------------+-----------
  16x16_ie      |     N    |      N          |  int16_t |
  32x16_ie      |     N    |      N          |  int32_t |
  32x32_ie      |     N    |      N          |  int32_t |
  f_ie          |     N    |      N          | float32_t|
  --------------+----------+-----------------+-----------

  Returned value: total number of right shifts occurred during scaling
  procedure

  Restrictions:
  x,y   should not overlap
  x,y   aligned on 16-bytes boundary
  x[(0)*2+1],
  x[(N/2)*2+1]  should be equal to zero
-------------------------------------------------------------------------*/
int ifft_real32x32_ie(int32_t* y, complex_fract32* x, const complex_fract32* twd, int twdstep, int N, int scalingOpt)
{
    const ae_int32x4 * restrict pX;
    int shift;
    int bexp; 

    NASSERT(scalingOpt==2 || scalingOpt==3); 
    NASSERT(x!=(complex_fract32*)y); 
    NASSERT_ALIGN16(x); 
    NASSERT_ALIGN16(y);
    NASSERT_ALIGN16(twd);


    if (scalingOpt == 2)
    {
        int numComplexSamples = (N>>1) + 1;
        ae_int32x2 x0, x1, x2, x3;
        ae_int16x4 nsa0, nsa1;
        int n;
        pX = (const ae_int32x4 *)x;
        nsa0 = 31; nsa1 = 31;
        NASSERT((N & 3) == 0);
        __Pragma("loop_count min=1");
        for (n = 0; n<(numComplexSamples >> 2); n++)
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
    else
    {
        bexp = 0;
    }

    shift  = iSpectrConv(x, N, twd, 3*twdstep, scalingOpt, &bexp);

    shift += __ifft_cplx32x32_ie((complex_fract32*)y, x, twd, twdstep * 2, N / 2, scalingOpt, bexp);
    return shift;


} /* ifft_real32x32_ie() */
