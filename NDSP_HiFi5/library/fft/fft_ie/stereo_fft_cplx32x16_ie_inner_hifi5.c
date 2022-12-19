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
int stereo_fft_cplx32x16_ie_inner(complex_fract32* x, const complex_fract16* twd, int twdstep, int N, int scalingOpt)
{
    int offset_inc = twdstep*sizeof(complex_fract16), s;
    int shift = 0;
    int tw_step = twdstep;
    int stride = N>>2;

    ae_int32x4 * restrict px0;
    ae_int32x4 * restrict px1;
    ae_int32x4 * restrict px2;
    ae_int32x4 * restrict px3;
    ae_int32x4 * restrict py0;
    ae_int32x4 * restrict py1;
    ae_int32x4 * restrict py2;
    ae_int32x4 * restrict py3;
    ae_int32x4 * restrict px = (ae_int32x4 *)x;
    ae_int32x4 * restrict py = (ae_int32x4 *)x;

    ae_int32   * restrict p16tw1 = (ae_int32*)(twd);
    ae_int32   * restrict p16tw2 = (ae_int32*)(twd + 1 * twdstep*(N >> 2));
    ae_int32   * restrict p16tw3 = (ae_int32*)(twd + 2 * twdstep*(N >> 2));

    ae_int32x2  vB0, vB1, vB2, vB3;

    ae_int32x2 t1, t2, t3;
    ae_f16x4 t1_f16x4;
    ae_f16x4 t2_f16x4;
    ae_f16x4 t3_f16x4;
    int i, j, bexp;
    int log2n = 0;

    NASSERT_ALIGN16(x);
    NASSERT(scalingOpt==3 || scalingOpt==2);
    NASSERT(N>=16 && (N&(N-1))==0);

    /*
     * autoscaling case
     */
    if ( scalingOpt == 2)
    {
        const int min_shift = 3;
        int32_t n;
        ae_int16x4 nsa0 = 31, nsa1 = 31;
        ae_int32x2 scl;
        int shiftl, shiftr;

        /* Calculate the exponent to prescale the input data */
        __Pragma("loop_count min=2");
        for (n = 0; n<(N >> 1); n++)
        {
            /* 2 cycles per pipeline stage in steady state with unroll = 1 */
            ae_int32x2 x0, x1, x2, x3;
            AE_L32X2X2_IP(x0, x1, px, sizeof(ae_int32x4));
            AE_L32X2X2_IP(x2, x3, px, sizeof(ae_int32x4));
            nsa0 = AE_MIN16(nsa0, AE_NSA32X4(x0, x1));
            nsa1 = AE_MIN16(nsa1, AE_NSA32X4(x2, x3));
        }
        bexp = AE_RMIN16X4(AE_MIN16(nsa0, nsa1));

        /* First stage */
        shift = min_shift - bexp;
        shiftl = XT_MAX(0, -shift);
        shiftr = XT_MAX(0,  shift);
        scl = 1 << shiftl;
        WUR_AE_SAR(shiftr);

        px = (ae_int32x4*)x;

        __Pragma("loop_count min=1");
        for (i = 0; i < (N >> 2); i++)
        {
            ae_int32x2  vA01, vA11, vA21, vA31, vC01, vC11, vC21, vC31;
            ae_int32x2  vA00, vA10, vA20, vA30, vC00, vC10, vC20, vC30;
            /* 10 cycles per pipeline stage in steady state with unroll=1 */

            AE_L32X2X2_X(vA30, vA31, px, 3 * (2*stride)*sizeof(ae_int32x2));
            AE_L32X2X2_X(vA20, vA21, px, 2 * (2*stride)*sizeof(ae_int32x2));
            AE_L32X2X2_X(vA10, vA11, px, 1 * (2*stride)*sizeof(ae_int32x2));
            AE_L32X2X2_IP(vA00, vA01, px, 2 * sizeof(ae_int32x2));

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
            
            vA01 = AE_MULP32X2(vA01, scl);
            vA11 = AE_MULP32X2(vA11, scl);
            vA21 = AE_MULP32X2(vA21, scl);
            vA31 = AE_MULP32X2(vA31, scl);

            AE_ADDANDSUBRNG32(vB0, vB2, vA01, vA21);
            AE_ADDANDSUBRNG32(vB1, vB3, vA11, vA31);
            vB3 = AE_MUL32JS(vB3);
            AE_ADDANDSUB32S(vC01, vC21, vB0, vB1);
            AE_ADDANDSUB32S(vC31, vC11, vB2, vB3);

            vC11 = AE_MULFC32X16RAS_L(vC11, t1_f16x4);
            vC21 = AE_MULFC32X16RAS_L(vC21, t2_f16x4);
            vC31 = AE_MULFC32X16RAS_L(vC31, t3_f16x4);

            AE_S32X2X2RNG_X(vC30, vC31, py, 3 * (2*stride)*sizeof(ae_int32x2));
            AE_S32X2X2RNG_X(vC20, vC21, py, 1 * (2*stride)*sizeof(ae_int32x2));
            AE_S32X2X2RNG_X(vC10, vC11, py, 2 * (2*stride)*sizeof(ae_int32x2));
            AE_S32X2X2RNG_IP(vC00, vC01, py, 2 * sizeof(ae_int32x2));
        }/* for (i = 0; i < (N >> 2); i++) */
        stride>>=2;
        tw_step<<=2;
    }
    else /* if ( scalingOpt == 2) */
    {
        WUR_AE_SAR(3);
        shift += 3;
        /* First stage static scaling */
        __Pragma("loop_count min=1");
        for (i = 0; i < (N >> 2); i++)
        {
            ae_int32x2  vA01, vA11, vA21, vA31, vC01, vC11, vC21, vC31;
            ae_int32x2  vA00, vA10, vA20, vA30, vC00, vC10, vC20, vC30;
            /* 10 cycles per pipeline stage in steady state with unroll=1 */

            AE_L32X2X2_X(vA30, vA31, px, 3 * (2 * stride)*sizeof(ae_int32x2));
            AE_L32X2X2_X(vA20, vA21, px, 2 * (2 * stride)*sizeof(ae_int32x2));
            AE_L32X2X2_X(vA10, vA11, px, 1 * (2 * stride)*sizeof(ae_int32x2));
            AE_L32X2X2_IP(vA00, vA01, px, 2 * sizeof(ae_int32x2));

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

            AE_ADDANDSUBRNG32(vB0, vB2, vA01, vA21);
            AE_ADDANDSUBRNG32(vB1, vB3, vA11, vA31);
            vB3 = AE_MUL32JS(vB3);
            AE_ADDANDSUB32S(vC01, vC21, vB0, vB1);
            AE_ADDANDSUB32S(vC31, vC11, vB2, vB3);

            vC11 = AE_MULFC32X16RAS_L(vC11, t1_f16x4);
            vC21 = AE_MULFC32X16RAS_L(vC21, t2_f16x4);
            vC31 = AE_MULFC32X16RAS_L(vC31, t3_f16x4);

            AE_S32X2X2RNG_X(vC30, vC31, py, 3 * (2 * stride)*sizeof(ae_int32x2));
            AE_S32X2X2RNG_X(vC20, vC21, py, 1 * (2 * stride)*sizeof(ae_int32x2));
            AE_S32X2X2RNG_X(vC10, vC11, py, 2 * (2 * stride)*sizeof(ae_int32x2));
            AE_S32X2X2RNG_IP(vC00, vC01, py, 2 * sizeof(ae_int32x2));
        }/* for (i = 0; i < (N >> 2); i++) */

        stride >>= 2;
        tw_step <<= 2;
    } /*if ( scalingOpt == 2) else ... */

    /* Intermediate stage */
    while( stride > 1 )
    {
        /* set pointers and access stride for twiddle factors */
        p16tw1 = (ae_int32*) (twd);
        p16tw2 = (ae_int32*) (twd+1*twdstep*(N>>2));
        p16tw3 = (ae_int32*) (twd+2*twdstep*(N>>2));
        offset_inc = tw_step*sizeof(complex_fract16);

        AE_CALCRNG3();
        s = RUR_AE_SAR();
        /* Set scaling 2 when scalingOpt = 3  (static scaling mode)*/
        XT_MOVEQZ(s, 2, scalingOpt - 3); 
        WUR_AE_SAR(s);

        shift += s;
        log2n += 2;

        for (j = 0; j < (stride>>1); j++)
        {
            /* load twiddle factors */
            AE_L32_XP( t1, p16tw1, offset_inc );
            AE_L32_XP( t2, p16tw2, offset_inc );
            t1 = AE_SEL32_LL(t1, t2);
            t1_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
            AE_L32_XP( t1, p16tw3, offset_inc );
            AE_L32_XP( t2, p16tw1, offset_inc );
            t1 = AE_SEL32_LL(t1, t2);
            t2_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);
            AE_L32_XP( t1, p16tw2, offset_inc );
            AE_L32_XP( t2, p16tw3, offset_inc );
            t1 = AE_SEL32_LL(t1, t2);
            t3_f16x4 = AE_MOVF16X4_FROMINT32X2(t1);

            px0 = (ae_int32x4*)x + j*2;
            py0 = (ae_int32x4*)x + j*2;
            px1 = px0 + stride;
            py1 = py0 + stride;
            px2 = px1 + stride;
            py2 = py1 + stride;
            px3 = px2 + stride;
            py3 = py2 + stride;

            __Pragma("loop_count min=2");
            for (i = 0; i < (1<<log2n); i++)
            {
                ae_int32x2  vA01, vA11, vA21, vA31, vC01, vC11, vC21, vC31;
                ae_int32x2  vA00, vA10, vA20, vA30, vC00, vC10, vC20, vC30;
                
                AE_L32X2X2_IP(vA00, vA01, px0, sizeof(ae_int32x4));
                AE_L32X2X2_IP(vA10, vA11, px1, sizeof(ae_int32x4));
                AE_L32X2X2_IP(vA20, vA21, px2, sizeof(ae_int32x4));
                AE_L32X2X2_IP(vA30, vA31, px3, sizeof(ae_int32x4));

                /* butterfly 0, left channel */
                AE_ADDANDSUBRNG32(vB0, vB2, vA00, vA20);
                AE_ADDANDSUBRNG32(vB1, vB3, vA10, vA30);
                vB3 = AE_MUL32JS(vB3);
                AE_ADDANDSUB32S(vC00, vC20, vB0, vB1);
                AE_ADDANDSUB32S(vC30, vC10, vB2, vB3);

                vC10 = AE_MULFC32X16RAS_H(vC10, t1_f16x4);
                vC20 = AE_MULFC32X16RAS_L(vC20, t1_f16x4);
                vC30 = AE_MULFC32X16RAS_H(vC30, t2_f16x4);

                /* butterfly 0, right channel */
                AE_ADDANDSUBRNG32(vB0, vB2, vA01, vA21);
                AE_ADDANDSUBRNG32(vB1, vB3, vA11, vA31);
                vB3 = AE_MUL32JS(vB3);
                AE_ADDANDSUB32S(vC01, vC21, vB0, vB1);
                AE_ADDANDSUB32S(vC31, vC11, vB2, vB3);

                vC11 = AE_MULFC32X16RAS_H(vC11, t1_f16x4);
                vC21 = AE_MULFC32X16RAS_L(vC21, t1_f16x4);
                vC31 = AE_MULFC32X16RAS_H(vC31, t2_f16x4);

                AE_S32X2X2RNG_IP(vC00, vC01, py0, sizeof(ae_int32x4));
                AE_S32X2X2RNG_IP(vC20, vC21, py1, sizeof(ae_int32x4));
                AE_S32X2X2RNG_IP(vC10, vC11, py2, sizeof(ae_int32x4));
                AE_S32X2X2RNG_IP(vC30, vC31, py3, sizeof(ae_int32x4));


                AE_L32X2X2_XP(vA00, vA01, px0, (4 * stride * 2 - 2)*sizeof(ae_int32x2));
                AE_L32X2X2_XP(vA10, vA11, px1, (4 * stride * 2 - 2)*sizeof(ae_int32x2));
                AE_L32X2X2_XP(vA20, vA21, px2, (4 * stride * 2 - 2)*sizeof(ae_int32x2));
                AE_L32X2X2_XP(vA30, vA31, px3, (4 * stride * 2 - 2)*sizeof(ae_int32x2));

                /* butterfly 1, left channel */
                AE_ADDANDSUBRNG32(vB0, vB2, vA00, vA20);
                AE_ADDANDSUBRNG32(vB1, vB3, vA10, vA30);
                vB3 = AE_MUL32JS(vB3);
                AE_ADDANDSUB32S(vC00, vC20, vB0, vB1);
                AE_ADDANDSUB32S(vC30, vC10, vB2, vB3);

                vC10 = AE_MULFC32X16RAS_L(vC10, t2_f16x4);
                vC20 = AE_MULFC32X16RAS_H(vC20, t3_f16x4);
                vC30 = AE_MULFC32X16RAS_L(vC30, t3_f16x4);

                /* butterfly 1, right channel */
                AE_ADDANDSUBRNG32(vB0, vB2, vA01, vA21);
                AE_ADDANDSUBRNG32(vB1, vB3, vA11, vA31);
                vB3 = AE_MUL32JS(vB3);
                AE_ADDANDSUB32S(vC01, vC21, vB0, vB1);
                AE_ADDANDSUB32S(vC31, vC11, vB2, vB3);

                vC11 = AE_MULFC32X16RAS_L(vC11, t2_f16x4);
                vC21 = AE_MULFC32X16RAS_H(vC21, t3_f16x4);
                vC31 = AE_MULFC32X16RAS_L(vC31, t3_f16x4);

                AE_S32X2X2RNG_XP(vC00, vC01, py0, (4 * stride * 2 - 2)*sizeof(ae_int32x2));
                AE_S32X2X2RNG_XP(vC20, vC21, py1, (4 * stride * 2 - 2)*sizeof(ae_int32x2));
                AE_S32X2X2RNG_XP(vC10, vC11, py2, (4 * stride * 2 - 2)*sizeof(ae_int32x2));
                AE_S32X2X2RNG_XP(vC30, vC31, py3, (4 * stride * 2 - 2)*sizeof(ae_int32x2));
            }
        }
        stride>>=2;
        tw_step<<=2;
    }
    return shift;
} /* stereo_fft_cplx32x16_ie_inner() */
