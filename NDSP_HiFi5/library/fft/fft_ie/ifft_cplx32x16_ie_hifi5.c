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

/* Last stage of computing IFFT */
static int ifft_stage_last_ie( int32_t *x, int32_t *y, int N );

/*-------------------------------------------------------------------------
  Inverse FFT on Complex Data with Optimized Memory Usage
  These functions make inverse FFT on complex data with optimized 
  memory usage.
  Scaling: 
      +-------------------+----------------------------------------+
      |      Function     |           Scaling options              |
      +-------------------+----------------------------------------+
      | ifft_cplx16x16_ie |  2 - 16-bit dynamic scaling            | 
      | ifft_cplx32x16_ie |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      | ifft_cplx32x32_ie |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      +-------------------+----------------------------------------+
  NOTES:
  1. Bit-reversing reordering is done here.
  2. FFT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after 
     the call
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
  S                   1 for ordinary (single channel) FFT, 2 - for stereo
                      input/outputs
  x[N*S]              complex input signal. Real and imaginary data 
                      are interleaved and real data goes first

  twd[N*twdstep*3/4]  twiddle factor table of a complex-valued FFT of 
                      size N*twdstep
  N                   FFT size
  twdstep             twiddle step 
  scalingOpt          scaling option (see table above)

  Output:
  y[N*S]              output spectrum. Real and imaginary data are 
                      interleaved and real data goes first

  Returned value:     total number of right shifts occurred during 
                      scaling procedure

  Restrictions:
  x,y   should not overlap
  x,y   aligned on 16-bytes boundary
-------------------------------------------------------------------------*/
int ifft_cplx32x16_ie(complex_fract32* y,complex_fract32* x, const complex_fract16* twd, int twdstep, int N, int scalingOpt)
{
    int offset_inc, s;
    int shift = 0;
    int tw_step = twdstep;
    int stride = N >> 2;

    ae_int32x4 * restrict px;
    ae_int32x4 * restrict py;

    ae_int32   * restrict p16tw1 = (ae_int32*)(twd);
    ae_int32   * restrict p16tw2 = (ae_int32*)(twd + 1 * twdstep*(N >> 2));
    ae_int32   * restrict p16tw3 = (ae_int32*)(twd + 2 * twdstep*(N >> 2));
    ae_int64   a0, a1, b0, b1, c0, c1, d0, d1;
    ae_int32x2  vB0, vB1, vB2, vB3;

    ae_int32x2 t1, t2, t3;
    ae_f16x4 t1_f16x4;
    ae_f16x4 t2_f16x4;
    ae_f16x4 t3_f16x4;
    int i, j, bexp;
    int log2n = 0;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(scalingOpt == 3 || scalingOpt == 2);
    NASSERT(N >= 16 && (N&(N - 1)) == 0);

    if (scalingOpt == 2)
    {
        const int min_shift = 3;
        int32_t i, n;
        ae_int16x4 nsa0 = 31, nsa1 = 31;
        px = (ae_int32x4 *)x;
        /*
        * autoscaling case
        */
        __Pragma("loop_count min=2");
        for (n = 0; n<(N >> 2); n++)
        {
            /* 2 cycles per pipeline stage in steady state with unroll = 1 */
            ae_int32x2 x0, x1, x2, x3;
            AE_L32X2X2_IP(x0, x1, px, sizeof(ae_int32x4));
            AE_L32X2X2_IP(x2, x3, px, sizeof(ae_int32x4));
            nsa0 = AE_MIN16(nsa0, AE_NSA32X4(x0, x1));
            nsa1 = AE_MIN16(nsa1, AE_NSA32X4(x2, x3));
        }
        bexp = AE_RMIN16X4(AE_MIN16(nsa0, nsa1));

        /* First stage static scaling*/
        {
            ae_int32x2 scl;
            int shiftl, shiftr;

            shift = min_shift - bexp;
            shiftl = XT_MAX(0, -shift);
            shiftr = XT_MAX(0, shift);

            scl = 1 << shiftl;
            WUR_AE_SAR(shiftr);

            offset_inc = tw_step*sizeof(complex_fract16);

            px = (ae_int32x4*)x;
            py = (ae_int32x4*)x;
            __Pragma("loop_count min=1");
            for (i = 0; i < (N >> 3); i++)
            {
                
                ae_int32x2  vA01, vA11, vA21, vA31, vC01, vC11, vC21, vC31;
                ae_int32x2  vA00, vA10, vA20, vA30, vC00, vC10, vC20, vC30;
                /* 10 cycles per pipeline stage in steady state with unroll=1 */

                AE_L64X2_X(d0, d1, (ae_int64x2*)px, 3 * stride*sizeof(ae_int32x2));
                AE_L64X2_X(c0, c1, (ae_int64x2*)px, 2 * stride*sizeof(ae_int32x2));
                AE_L64X2_X(b0, b1, (ae_int64x2*)px, 1 * stride*sizeof(ae_int32x2));
                AE_L64X2_IP(a0, a1, castxcc(ae_int64x2, px), 2 * sizeof(ae_int32x2));

                vA00 = AE_MOVINT32X2_FROMINT64(a0);
                vA10 = AE_MOVINT32X2_FROMINT64(b0);
                vA20 = AE_MOVINT32X2_FROMINT64(c0);
                vA30 = AE_MOVINT32X2_FROMINT64(d0);
                vA01 = AE_MOVINT32X2_FROMINT64(a1);
                vA11 = AE_MOVINT32X2_FROMINT64(b1);
                vA21 = AE_MOVINT32X2_FROMINT64(c1);
                vA31 = AE_MOVINT32X2_FROMINT64(d1);

                vA00 = AE_MULP32X2(vA00, scl);
                vA10 = AE_MULP32X2(vA10, scl);
                vA20 = AE_MULP32X2(vA20, scl);
                vA30 = AE_MULP32X2(vA30, scl);

                AE_ADDANDSUBRNG32(vB0, vB2, vA00, vA20);
                AE_ADDANDSUBRNG32(vB1, vB3, vA10, vA30);
                vB3 = AE_MUL32JS(vB3);
                AE_ADDANDSUB32S(vC00, vC20, vB0, vB1);
                AE_ADDANDSUB32S(vC30, vC10, vB2, vB3);

                AE_L32_XP(t1, p16tw1, offset_inc);
                AE_L32_XP(t2, p16tw2, offset_inc);
                AE_L32_XP(t3, p16tw3, offset_inc);

                t1_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
                t2_f16x4 = AE_MOVF16X4_FROMINT32X2(t2);
                t3_f16x4 = AE_MOVF16X4_FROMINT32X2(t3);
                vC10 = AE_MULFC32X16RAS_L(vC10, t1_f16x4);
                vC20 = AE_MULFC32X16RAS_L(vC20, t2_f16x4);
                vC30 = AE_MULFC32X16RAS_L(vC30, t3_f16x4);
                ///////////
                vA01 = AE_MULP32X2(vA01, scl);
                vA11 = AE_MULP32X2(vA11, scl);
                vA21 = AE_MULP32X2(vA21, scl);
                vA31 = AE_MULP32X2(vA31, scl);

                AE_ADDANDSUBRNG32(vB0, vB2, vA01, vA21);
                AE_ADDANDSUBRNG32(vB1, vB3, vA11, vA31);
                vB3 = AE_MUL32JS(vB3);
                AE_ADDANDSUB32S(vC01, vC21, vB0, vB1);
                AE_ADDANDSUB32S(vC31, vC11, vB2, vB3);

                AE_L32_XP(t1, p16tw1, offset_inc);
                AE_L32_XP(t2, p16tw2, offset_inc);
                AE_L32_XP(t3, p16tw3, offset_inc);

                t1_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
                t2_f16x4 = AE_MOVF16X4_FROMINT32X2(t2);
                t3_f16x4 = AE_MOVF16X4_FROMINT32X2(t3);
                vC11 = AE_MULFC32X16RAS_L(vC11, t1_f16x4);
                vC21 = AE_MULFC32X16RAS_L(vC21, t2_f16x4);
                vC31 = AE_MULFC32X16RAS_L(vC31, t3_f16x4);

                AE_S32X2X2RNG_X(vC30, vC31, py, 3 * stride*sizeof(ae_int32x2));
                AE_S32X2X2RNG_X(vC20, vC21, py, 1 * stride*sizeof(ae_int32x2));
                AE_S32X2X2RNG_X(vC10, vC11, py, 2 * stride*sizeof(ae_int32x2));
                AE_S32X2X2RNG_IP(vC00, vC01, py, 2 * sizeof(ae_int32x2));
            }

            stride >>= 2;
            tw_step <<= 2;
        }
    }
    else
    {
        /* First stage  dynamic scaling*/
        {
            shift = 3;
            WUR_AE_SAR(shift);

            offset_inc = tw_step*sizeof(complex_fract16);

            px = (ae_int32x4*)x;
            py = (ae_int32x4*)x;
            __Pragma("loop_count min=2 factor=2");
            for (i = 0; i < (N >> 3); i++)
            {
                ae_int32x2  vA01, vA11, vA21, vA31, vC01, vC11, vC21, vC31;
                ae_int32x2  vA00, vA10, vA20, vA30, vC00, vC10, vC20, vC30;
                /* 15 cycles per pipeline stage in steady state with unroll=2 */

                AE_L64X2_X(d0, d1, (ae_int64x2*)px, 3 * stride*sizeof(ae_int32x2));
                AE_L64X2_X(c0, c1, (ae_int64x2*)px, 2 * stride*sizeof(ae_int32x2));
                AE_L64X2_X(b0, b1, (ae_int64x2*)px, 1 * stride*sizeof(ae_int32x2));
                AE_L64X2_IP(a0, a1, castxcc(ae_int64x2, px), 2 * sizeof(ae_int32x2));

                vA00 = AE_MOVINT32X2_FROMINT64(a0);
                vA10 = AE_MOVINT32X2_FROMINT64(b0);
                vA20 = AE_MOVINT32X2_FROMINT64(c0);
                vA30 = AE_MOVINT32X2_FROMINT64(d0);
                vA01 = AE_MOVINT32X2_FROMINT64(a1);
                vA11 = AE_MOVINT32X2_FROMINT64(b1);
                vA21 = AE_MOVINT32X2_FROMINT64(c1);
                vA31 = AE_MOVINT32X2_FROMINT64(d1);

                AE_ADDANDSUBRNG32(vB0, vB2, vA00, vA20);
                AE_ADDANDSUBRNG32(vB1, vB3, vA10, vA30);
                vB3 = AE_MUL32JS(vB3);
                AE_ADDANDSUB32S(vC00, vC20, vB0, vB1);
                AE_ADDANDSUB32S(vC30, vC10, vB2, vB3);

                AE_L32_XP(t1, p16tw1, offset_inc);
                AE_L32_XP(t2, p16tw2, offset_inc);
                AE_L32_XP(t3, p16tw3, offset_inc);

                t1_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
                t2_f16x4 = AE_MOVF16X4_FROMINT32X2(t2);
                t3_f16x4 = AE_MOVF16X4_FROMINT32X2(t3);
                vC10 = AE_MULFC32X16RAS_L(vC10, t1_f16x4);
                vC20 = AE_MULFC32X16RAS_L(vC20, t2_f16x4);
                vC30 = AE_MULFC32X16RAS_L(vC30, t3_f16x4);
                ///////////

                AE_ADDANDSUBRNG32(vB0, vB2, vA01, vA21);
                AE_ADDANDSUBRNG32(vB1, vB3, vA11, vA31);
                vB3 = AE_MUL32JS(vB3);
                AE_ADDANDSUB32S(vC01, vC21, vB0, vB1);
                AE_ADDANDSUB32S(vC31, vC11, vB2, vB3);

                AE_L32_XP(t1, p16tw1, offset_inc);
                AE_L32_XP(t2, p16tw2, offset_inc);
                AE_L32_XP(t3, p16tw3, offset_inc);

                t1_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
                t2_f16x4 = AE_MOVF16X4_FROMINT32X2(t2);
                t3_f16x4 = AE_MOVF16X4_FROMINT32X2(t3);
                vC11 = AE_MULFC32X16RAS_L(vC11, t1_f16x4);
                vC21 = AE_MULFC32X16RAS_L(vC21, t2_f16x4);
                vC31 = AE_MULFC32X16RAS_L(vC31, t3_f16x4);

                AE_S32X2X2RNG_X(vC30, vC31, py, 3 * stride*sizeof(ae_int32x2));
                AE_S32X2X2RNG_X(vC20, vC21, py, 1 * stride*sizeof(ae_int32x2));
                AE_S32X2X2RNG_X(vC10, vC11, py, 2 * stride*sizeof(ae_int32x2));
                AE_S32X2X2RNG_IP(vC00, vC01, py, 2 * sizeof(ae_int32x2));
            }
            stride >>= 2;
            tw_step <<= 2;
        }
    }
    /* Intermediate stage */
    while (stride > 1)
    {
        /* set pointers and access stride for twiddle factors */
        p16tw1 = (ae_int32*)(twd);
        p16tw2 = (ae_int32*)(twd + 1 * twdstep*(N >> 2));
        p16tw3 = (ae_int32*)(twd + 2 * twdstep*(N >> 2));
        offset_inc = tw_step*sizeof(complex_fract16);
        if (scalingOpt == 2)
        {
            /* Dynamic scaling */
            AE_CALCRNG3();
            bexp = 3 - RUR_AE_SAR();
            s = XT_MAX(0, 3 - bexp);
        }
        else
        {
            /* Static scaling */
            s = 2;
        }
        shift += s;

        WUR_AE_SAR(s);

        log2n += 2;
        __Pragma("loop_count min=1");
        for (j = 0; j < (stride >> 1); j++)
        {
            /* load twiddle factors */
            AE_L32_XP(t1, p16tw1, offset_inc);
            AE_L32_XP(t2, p16tw2, offset_inc);
            t1 = AE_SEL32_LL(t1, t2);
            t1_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
            AE_L32_XP(t1, p16tw3, offset_inc);
            AE_L32_XP(t2, p16tw1, offset_inc);
            t1 = AE_SEL32_LL(t1, t2);
            t2_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
            AE_L32_XP(t1, p16tw2, offset_inc);
            AE_L32_XP(t2, p16tw3, offset_inc);
            t1 = AE_SEL32_LL(t1, t2);
            t3_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);

            px = (ae_int32x4*)x + j;
            py = (ae_int32x4*)x + j;

            __Pragma("loop_count min=2");
            for (i = 0; i < (1 << log2n); i++)
            {
                ae_int32x2  vA01, vA11, vA21, vA31, vC01, vC11, vC21, vC31;
                ae_int32x2  vA00, vA10, vA20, vA30, vC00, vC10, vC20, vC30;

                /* 6 cycles per pipeline stage in steady state with unroll=1 */

                AE_L32X2X2_XP(vA00, vA01, px, stride*sizeof(ae_int32x2));
                AE_L32X2X2_XP(vA10, vA11, px, stride*sizeof(ae_int32x2));
                AE_L32X2X2_XP(vA20, vA21, px, stride*sizeof(ae_int32x2));
                AE_L32X2X2_XP(vA30, vA31, px, stride*sizeof(ae_int32x2));
                /* butterfly 0 */
                AE_ADDANDSUBRNG32(vB0, vB2, vA00, vA20);
                AE_ADDANDSUBRNG32(vB1, vB3, vA10, vA30);
                vB3 = AE_MUL32JS(vB3);
                AE_ADDANDSUB32S(vC00, vC20, vB0, vB1);
                AE_ADDANDSUB32S(vC30, vC10, vB2, vB3);

                vC10 = AE_MULFC32X16RAS_H(vC10, t1_f16x4);
                vC20 = AE_MULFC32X16RAS_L(vC20, t1_f16x4);
                vC30 = AE_MULFC32X16RAS_H(vC30, t2_f16x4);

                /* butterfly 1 */
                AE_ADDANDSUBRNG32(vB0, vB2, vA01, vA21);
                AE_ADDANDSUBRNG32(vB1, vB3, vA11, vA31);
                vB3 = AE_MUL32JS(vB3);
                AE_ADDANDSUB32S(vC01, vC21, vB0, vB1);
                AE_ADDANDSUB32S(vC31, vC11, vB2, vB3);

                vC11 = AE_MULFC32X16RAS_L(vC11, t2_f16x4);
                vC21 = AE_MULFC32X16RAS_H(vC21, t3_f16x4);
                vC31 = AE_MULFC32X16RAS_L(vC31, t3_f16x4);

                AE_S32X2X2RNG_XP(vC00, vC01, py, stride*sizeof(ae_int32x2));
                AE_S32X2X2RNG_XP(vC20, vC21, py, stride*sizeof(ae_int32x2));
                AE_S32X2X2RNG_XP(vC10, vC11, py, stride*sizeof(ae_int32x2));
                AE_S32X2X2RNG_XP(vC30, vC31, py, stride*sizeof(ae_int32x2));
            } /* for (i = 0; i < (1<<log2n); i++) */
        } /* for (j = 0; j < (stride>>1); j++) */

        stride >>= 2;
        tw_step <<= 2;
    }


    shift += ifft_stage_last_ie((int32_t*)x, (int32_t*)y, N);

    return shift;
} /* ifft_cplx32x16_ie() */

/* Last stage of computing IFFT */
int ifft_stage_last_ie( int32_t *x,
                        int32_t *y,
                        int N)
{
    int32_t i, i0, i1, ai;
    ae_int32x2 * restrict p_y0 = (ae_int32x2 *)(y);
    ae_int32x2 * restrict p_y1 = (p_y0 + (N >> 2));
    ae_int32x2 * restrict p_y2 = (p_y1 + (N >> 2));
    ae_int32x2 * restrict p_y3 = (p_y2 + (N >> 2));
    ae_int32x4 * restrict p_x0 = (ae_int32x4 *)(x);
    int shift;

    i = NSA(N) + 1;
    ai = ((int32_t)0x1) << i;
    i0 = 0;


    if ((i & 1) == 0)
    {
        shift = 1;    //Select scaling
        WUR_AE_SAR(shift);
        //--------------------------------------------------------------------------
        // last stage is RADIX2 !!!
        //--------------------------------------------------------------------------
        __Pragma("loop_count min=1");
        for (i = 0; i < (N >> 3); i++)
        {
            /* 8 cycles per pipeline stage in steady state with unroll=1 */
            ae_int32x2 vA1, vA2, vA3, vA0;
            ae_int32x2 vB1, vB2, vB3, vB0;
            i1 = AE_ADDBRBA32(i0, ai);

            // FFT_BUTTERFLY_R2(i0, shift);
            AE_L32X2X2_IP(vA0, vA1, p_x0, sizeof(ae_int32x4));
            AE_L32X2X2_IP(vA2, vA3, p_x0, sizeof(ae_int32x4));
            AE_ADDANDSUBRNG32(vB0, vB2, vA0, vA1);
            AE_ADDANDSUBRNG32(vB1, vB3, vA2, vA3);

            vB0 = AE_SEL32_LH(vB0, vB0);
            vB1 = AE_SEL32_LH(vB1, vB1);
            vB2 = AE_SEL32_LH(vB2, vB2);
            vB3 = AE_SEL32_LH(vB3, vB3);

            AE_S32X2_X(vB0, p_y0, i0);
            AE_S32X2_X(vB1, p_y1, i0);
            AE_S32X2_X(vB2, p_y2, i0);
            AE_S32X2_X(vB3, p_y3, i0);

            //FFT_BUTTERFLY_R2(i1, shift);
            AE_L32X2X2_IP(vA0, vA1, p_x0, sizeof(ae_int32x4));
            AE_L32X2X2_IP(vA2, vA3, p_x0, sizeof(ae_int32x4));
            AE_ADDANDSUBRNG32(vB0, vB2, vA0, vA1);
            AE_ADDANDSUBRNG32(vB1, vB3, vA2, vA3);

            vB0 = AE_SEL32_LH(vB0, vB0);
            vB1 = AE_SEL32_LH(vB1, vB1);
            vB2 = AE_SEL32_LH(vB2, vB2);
            vB3 = AE_SEL32_LH(vB3, vB3);

            AE_S32X2_X(vB0, p_y0, i1);
            AE_S32X2_X(vB1, p_y1, i1);
            AE_S32X2_X(vB2, p_y2, i1);
            AE_S32X2_X(vB3, p_y3, i1);

            i0 = AE_ADDBRBA32(i1, ai);
        }
    }
    else
    {
        shift = 2;    //Select scaling
        WUR_AE_SAR(shift);
        //--------------------------------------------------------------------------
        // last stage is RADIX4 !!!
        //--------------------------------------------------------------------------
        __Pragma("loop_count min=1");
        for (i = 0; i<(N >> 3); i++)
        {
            /* 8 cycles per pipeline stage in steady state with unroll=1 */
            ae_int32x2 vA1, vA2, vA3, vA0;
            ae_int32x2 vB1, vB2, vB3, vB0;

            //     FFT_BUTTERFLY_R4(i0, shift);
            AE_L32X2X2_IP(vA0, vA1, p_x0, sizeof(ae_int32x4));
            AE_L32X2X2_IP(vA2, vA3, p_x0, sizeof(ae_int32x4));
            AE_ADDANDSUBRNG32(vB0, vB2, vA0, vA2);
            AE_ADDANDSUBRNG32(vB1, vB3, vA1, vA3);
            vB3 = AE_MUL32JS(vB3);
            AE_ADDANDSUB32S(vA0, vA2, vB0, vB1);
            AE_ADDANDSUB32S(vA3, vA1, vB2, vB3);

            vA0 = AE_SEL32_LH(vA0, vA0);
            vA1 = AE_SEL32_LH(vA1, vA1);
            vA2 = AE_SEL32_LH(vA2, vA2);
            vA3 = AE_SEL32_LH(vA3, vA3);

            AE_S32X2_X(vA0, p_y0, i0);
            AE_S32X2_X(vA1, p_y1, i0);
            AE_S32X2_X(vA2, p_y2, i0);
            AE_S32X2_X(vA3, p_y3, i0);
            i0 = AE_ADDBRBA32(i0, ai);

            //     FFT_BUTTERFLY_R4(i0, shift);
            AE_L32X2X2_IP(vA0, vA1, p_x0, sizeof(ae_int32x4));
            AE_L32X2X2_IP(vA2, vA3, p_x0, sizeof(ae_int32x4));
            AE_ADDANDSUBRNG32(vB0, vB2, vA0, vA2);
            AE_ADDANDSUBRNG32(vB1, vB3, vA1, vA3);
            vB3 = AE_MUL32JS(vB3);
            AE_ADDANDSUB32S(vA0, vA2, vB0, vB1);
            AE_ADDANDSUB32S(vA3, vA1, vB2, vB3);

            vA0 = AE_SEL32_LH(vA0, vA0);
            vA1 = AE_SEL32_LH(vA1, vA1);
            vA2 = AE_SEL32_LH(vA2, vA2);
            vA3 = AE_SEL32_LH(vA3, vA3);

            AE_S32X2_X(vA0, p_y0, i0);
            AE_S32X2_X(vA1, p_y1, i0);
            AE_S32X2_X(vA2, p_y2, i0);
            AE_S32X2_X(vA3, p_y3, i0);
            i0 = AE_ADDBRBA32(i0, ai);
        }
    }
    return shift;
} /* ifft_stage_last_ie() */
