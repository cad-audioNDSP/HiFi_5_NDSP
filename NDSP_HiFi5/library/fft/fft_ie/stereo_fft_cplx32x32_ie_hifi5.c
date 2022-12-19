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

#define SWAP_PTR(_x, _y) {complex_fract32 *tmp = _x; _x = _y ; _y = tmp; } 

/* Radix-4 butterfly with normalization   *
 * x0, x1, x2, x3 - input/output samples  *
 * AE_SAR - contains normalization factor */
#define DFT4X1RNG(x0, x1, x2, x3)\
{   \
    ae_int32x2 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG32(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32(s1, d1, x1, x3); \
    d1 = AE_MUL32JS(d1);               \
    AE_ADDANDSUB32S(x0, x2, s0, s1);   \
    AE_ADDANDSUB32S(x3, x1, d0, d1);   \
}




/*------------------------------------------------------------------------------
Radix-4 complex-valued FFT of size 2^n, n=4..12.
This function computes the inner stages.

Notes:
1. Inner means  radix-4 FFT FFT stages excluding first and the last stages.
2. At each stage data buffer is automatically downscaled with appropriate shifts
to avoid overflows if necessary.

Precision:
32-bit  complex input/output data, 32-bit complex twiddle factors

Input/Output:
x[2*2*N]        Interleaved stereo complex input data. Real and imaginary data are
interleaved to form single complex data where real data goes first.
y[2*2*N]        Interleaved stereo complex output data. Real and imaginary data are
interleaved to form single complex data where real data goes first.
twd[2*N*3/4]    Complex twiddle factor table, twiddle data is same for both channel.
*bexp           Common block exponent, that is the minimum number of redundant sign
bits over input (output) data
N               FFT Size

Returned value: total right shift amount applied to dynamically scale the data

Restrictions:
x[],y[],twd[] - must not overlap
x[],y[],twd[] - must be aligned on 16-byte boundary
------------------------------------------------------------------------------*/

static int fft_cplx_inner_32x32_scl2(
    complex_fract32 * restrict y,
    complex_fract32 * restrict x,
    const complex_fract32 *    twd,
    int                        twdstep, 
    int *                      bexp,
    int                        N
    )
{
    const ae_int32x2 *          X0;
    const ae_int32x2 *          X1;
    const ae_int32x2 *          X2;
    const ae_int32x2 *          X3;
    ae_int32x2 * restrict Y0;
    ae_int32x2 * restrict Y1;
    ae_int32x2 * restrict Y2;
    ae_int32x2 * restrict Y3;
    const ae_f32x2   *          TWD;

    int m, n, logN;
    int stride;
    int nsa, shift, shiftSum;
    int shift_ch1, shift_ch2;

    logN = 30 - NSA(N);

    nsa = *bexp;

    shiftSum = 0;

    //----------------------------------------------------------------------------
    // Perform second through the second to the last stages.
    {
        complex_fract32 * X;
        complex_fract32 * Y;

        ae_int32x2 a0, a1, a2, a3;
        ae_int32x2 b0, b1, b2, b3;
        ae_f32x2   c0, c1, c2, c3;

        ae_int32x2 d0, d1, d2, d3;
        ae_int32x2 e0, e1, e2, e3;
        ae_f32x2   f0, f1, f2, f3;

        ae_f32x2 tw1, tw2, tw3;

        X = (((logN + 1) & 2) ? x : y);
        Y = (((logN + 1) & 2) ? y : x);

        shiftSum += (shift = 3 - nsa);

        WUR_AE_SAR(shift);
        twdstep*=4;
        for (stride = N / 16; stride>4; stride /= 4, twdstep *= 4)
        {

            X0 = (ae_int32x2*)((unsigned int)X + 0 * stride * 8 * 2);
            X1 = (ae_int32x2*)((unsigned int)X + 1 * stride * 8 * 2);
            X2 = (ae_int32x2*)((unsigned int)X + 2 * stride * 8 * 2);
            X3 = (ae_int32x2*)((unsigned int)X + 3 * stride * 8 * 2);

            Y0 = (ae_int32x2*)((unsigned int)Y + 0 * stride * 8 * 2);
            Y1 = (ae_int32x2*)((unsigned int)Y + 1 * stride * 8 * 2);
            Y2 = (ae_int32x2*)((unsigned int)Y + 2 * stride * 8 * 2);
            Y3 = (ae_int32x2*)((unsigned int)Y + 3 * stride * 8 * 2);

            __Pragma("loop_count min=1");
            for (m = 0; m*(4 * stride)<N; m++)
            {
                TWD = (ae_f32x2*)twd;

                __Pragma("ymemory( X0 )");
                __Pragma("ymemory( X1 )");
                __Pragma("ymemory( X2 )");
                __Pragma("ymemory( X3 )");
                __Pragma("loop_count min=4, factor=4");
                for (n = 0; n<stride; n++)
                {
                    tw3 = AE_L32X2_I((ae_int32x2*)TWD, 2 * 8);
                    AE_L32X2X2_XP(tw1, tw2, castxcc(ae_int32x4, TWD), +3 * twdstep * 8);

                    //
                    // Group 0
                    //

                    AE_L32X2X2_XP(a0, d0, castxcc(ae_int32x4,X0), +8 * 2);
                    AE_L32X2X2_XP(a1, d1, castxcc(ae_int32x4,X1), +8 * 2);
                    AE_L32X2X2_XP(a2, d2, castxcc(ae_int32x4,X2), +8 * 2);
                    AE_L32X2X2_XP(a3, d3, castxcc(ae_int32x4,X3), +8 * 2);

                    AE_ADDANDSUBRNG32_H(b0, b2, a0, a2);
                    AE_ADDANDSUBRNG32_H(b1, b3, a3, a1);

                    AE_ADDANDSUBRNG32_L(e0, e2, d0, d2);
                    AE_ADDANDSUBRNG32_L(e1, e3, d3, d1);

                    AE_ADDANDSUB32S(c0, c2, b0, b1);
                    AE_ADDANDSUB32JS(c1, c3, b2, b3);

                    AE_ADDANDSUB32S(f0, f2, e0, e1);
                    AE_ADDANDSUB32JS(f1, f3, e2, e3);

                    c1 = AE_MULFC32RAS(c1, tw1);
                    c2 = AE_MULFC32RAS(c2, tw2);
                    c3 = AE_MULFC32RAS(c3, tw3);

                    f1 = AE_MULFC32RAS(f1, tw1);
                    f2 = AE_MULFC32RAS(f2, tw2);
                    f3 = AE_MULFC32RAS(f3, tw3);

                    AE_S32X2X2RNG_XP(c0, f0, castxcc(ae_int32x4,Y0), +8 * 2);
                    AE_S32X2X2RNG_XP(c2, f2, castxcc(ae_int32x4,Y1), +8 * 2);
                    AE_S32X2X2RNG_XP(c1, f1, castxcc(ae_int32x4,Y2), +8 * 2);
                    AE_S32X2X2RNG_XP(c3, f3, castxcc(ae_int32x4,Y3), +8 * 2);

                }
                X0 += 3 * stride * 2;
                X1 += 3 * stride * 2;
                X2 += 3 * stride * 2;
                X3 += 3 * stride * 2;

                Y0 += 3 * stride * 2;
                Y1 += 3 * stride * 2;
                Y2 += 3 * stride * 2;
                Y3 += 3 * stride * 2;
            }

            {
                // Swap input/output buffers between successive stages.
                complex_fract32 * T = X; X = Y; Y = T;
            }

            AE_CALCRNG32(shift_ch1, shift_ch2, 0, 3);
            (void)shift_ch2;
            shiftSum += shift_ch1;

        }
    }

    __Pragma("no_reorder");

    //----------------------------------------------------------------------------
    // Perform the next to last stage.
    {
        ae_int32x2 a0, a1, a2, a3;
        ae_int32x2 b0, b1, b2, b3;
        ae_f32x2   c0, c1, c2, c3;

        ae_int32x2 d0, d1, d2, d3;
        ae_int32x2 e0, e1, e2, e3;
        ae_f32x2   f0, f1, f2, f3;

        TWD = (const ae_f32x2*)((unsigned int)twd + 3 * twdstep * 8);

        X0 = (ae_int32x2*)x;
        Y0 = (ae_int32x2*)y;

        if (stride == 4)
        {
            //
            // Next to last stage for FFT size an even power of two.
            //

            ae_f32x2 tw11, tw12, tw13;
            ae_f32x2 tw21, tw22, tw23;
            ae_f32x2 tw31, tw32, tw33;

            tw21 = ae_f32x2_loadi(TWD, 1 * 8);
            tw31 = ae_f32x2_loadi(TWD, 2 * 8);
            ae_f32x2_loadxp(tw11, TWD, +3 * twdstep * 8);

            tw22 = ae_f32x2_loadi(TWD, 1 * 8);
            tw32 = ae_f32x2_loadi(TWD, 2 * 8);
            ae_f32x2_loadxp(tw12, TWD, +3 * twdstep * 8);

            tw23 = ae_f32x2_loadi(TWD, 1 * 8);
            tw33 = ae_f32x2_loadi(TWD, 2 * 8);
            ae_f32x2_loadxp(tw13, TWD, +3 * twdstep * 8);

            __Pragma("loop_count min=4, factor=4");
            for (n = 0; n<N / (4 * 4); n++)
            {
                //
                // Group 0
                //

                AE_L32X2X2_XP(a0, d0, castxcc(ae_int32x4,X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a1, d1, castxcc(ae_int32x4,X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a2, d2, castxcc(ae_int32x4,X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a3, d3, castxcc(ae_int32x4,X0), -11 * 8 * 2);

                AE_ADDANDSUBRNG32(b0, b2, a0, a2);
                AE_ADDANDSUBRNG32(b1, b3, a3, a1);

                AE_ADDANDSUBRNG32(e0, e2, d0, d2);
                AE_ADDANDSUBRNG32(e1, e3, d3, d1);

                AE_ADDANDSUB32S(c0, c2, b0, b1);
                AE_ADDANDSUB32JS(c1, c3, b2, b3);

                AE_ADDANDSUB32S(f0, f2, e0, e1);
                AE_ADDANDSUB32JS(f1, f3, e2, e3);

                AE_S32X2X2RNG_XP(c0, f0, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c2, f2, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c1, f1, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c3, f3, castxcc(ae_int32x4, Y0), -11 * 8 * 2);

                //
                // Group 1
                //
                AE_L32X2X2_XP(a0, d0, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a1, d1, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a2, d2, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a3, d3, castxcc(ae_int32x4, X0), -11 * 8 * 2);

                AE_ADDANDSUBRNG32(b0, b2, a0, a2);
                AE_ADDANDSUBRNG32(b1, b3, a3, a1);

                AE_ADDANDSUBRNG32(e0, e2, d0, d2);
                AE_ADDANDSUBRNG32(e1, e3, d3, d1);

                AE_ADDANDSUB32S(c0, c2, b0, b1);
                AE_ADDANDSUB32JS(c1, c3, b2, b3);

                AE_ADDANDSUB32S(f0, f2, e0, e1);
                AE_ADDANDSUB32JS(f1, f3, e2, e3);

                c1 = AE_MULFC32RAS(c1, tw11);
                c2 = AE_MULFC32RAS(c2, tw21);
                c3 = AE_MULFC32RAS(c3, tw31);

                f1 = AE_MULFC32RAS(f1, tw11);
                f2 = AE_MULFC32RAS(f2, tw21);
                f3 = AE_MULFC32RAS(f3, tw31);

                AE_S32X2X2RNG_XP(c0, f0, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c2, f2, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c1, f1, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c3, f3, castxcc(ae_int32x4, Y0), -11 * 8 * 2);

                //
                // Group 2
                //
                AE_L32X2X2_XP(a0, d0, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a1, d1, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a2, d2, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a3, d3, castxcc(ae_int32x4, X0), -11 * 8 * 2);

                AE_ADDANDSUBRNG32(b0, b2, a0, a2);
                AE_ADDANDSUBRNG32(b1, b3, a3, a1);

                AE_ADDANDSUBRNG32(e0, e2, d0, d2);
                AE_ADDANDSUBRNG32(e1, e3, d3, d1);

                AE_ADDANDSUB32S(c0, c2, b0, b1);
                AE_ADDANDSUB32JS(c1, c3, b2, b3);

                AE_ADDANDSUB32S(f0, f2, e0, e1);
                AE_ADDANDSUB32JS(f1, f3, e2, e3);

                c1 = AE_MULFC32RAS(c1, tw12);
                c2 = AE_MULFC32RAS(c2, tw22);
                c3 = AE_MULFC32RAS(c3, tw32);

                f1 = AE_MULFC32RAS(f1, tw12);
                f2 = AE_MULFC32RAS(f2, tw22);
                f3 = AE_MULFC32RAS(f3, tw32);

                AE_S32X2X2RNG_XP(c0, f0, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c2, f2, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c1, f1, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c3, f3, castxcc(ae_int32x4, Y0), -11 * 8 * 2);

                //
                // Group 3
                //
                AE_L32X2X2_XP(a0, d0, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a1, d1, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a2, d2, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a3, d3, castxcc(ae_int32x4, X0), +1 * 8 * 2);

                AE_ADDANDSUBRNG32(b0, b2, a0, a2);
                AE_ADDANDSUBRNG32(b1, b3, a3, a1);

                AE_ADDANDSUBRNG32(e0, e2, d0, d2);
                AE_ADDANDSUBRNG32(e1, e3, d3, d1);

                AE_ADDANDSUB32S(c0, c2, b0, b1);
                AE_ADDANDSUB32JS(c1, c3, b2, b3);

                AE_ADDANDSUB32S(f0, f2, e0, e1);
                AE_ADDANDSUB32JS(f1, f3, e2, e3);

                c1 = AE_MULFC32RAS(c1, tw13);
                c2 = AE_MULFC32RAS(c2, tw23);
                c3 = AE_MULFC32RAS(c3, tw33);

                f1 = AE_MULFC32RAS(f1, tw13);
                f2 = AE_MULFC32RAS(f2, tw23);
                f3 = AE_MULFC32RAS(f3, tw33);

                AE_S32X2X2RNG_XP(c0, f0, castxcc(ae_int32x4,Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c2, f2, castxcc(ae_int32x4,Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c1, f1, castxcc(ae_int32x4,Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c3, f3, castxcc(ae_int32x4,Y0), +1 * 8 * 2);
            }

            AE_CALCRNG32(shift_ch1, shift_ch2, 0, 3);

            nsa = 3 - shift_ch1;

        }
        else if (stride == 2)
        {
            //
            // Next to last stage for FFT size an odd power of two.
            //

            ae_f32x2 tw1, tw2, tw3;


            tw1 = ae_f32x2_loadi(TWD, 0 * 8);
            tw2 = ae_f32x2_loadi(TWD, 1 * 8);
            tw3 = ae_f32x2_loadi(TWD, 2 * 8);

            __Pragma("loop_count min=4, factor=4");
            for (n = 0; n<N / (2 * 4); n++)
            {
                //
                // Group 0
                //

                d0 = AE_L32X2_I(X0, +8);
                AE_L32X2_IP(a0, X0, +2 * 8 * 2);

                d1 = AE_L32X2_I(X0, +8);
                AE_L32X2_IP(a1, X0, +2 * 8 * 2);

                d2 = AE_L32X2_I(X0, +8);
                AE_L32X2_IP(a2, X0, +2 * 8 * 2);

                d3 = AE_L32X2_I(X0, +8);
                AE_L32X2_XP(a3, X0, -5 * 8 * 2);

                AE_ADDANDSUBRNG32(b0, b2, a0, a2);
                AE_ADDANDSUBRNG32(b1, b3, a3, a1);

                AE_ADDANDSUBRNG32(e0, e2, d0, d2);
                AE_ADDANDSUBRNG32(e1, e3, d3, d1);

                b3 = AE_MUL32JS(b3);
                e3 = AE_MUL32JS(e3);

                AE_ADDANDSUB32S(c0, c2, b0, b1);
                AE_ADDANDSUB32S(c1, c3, b2, b3);

                AE_ADDANDSUB32S(f0, f2, e0, e1);
                AE_ADDANDSUB32S(f1, f3, e2, e3);

                AE_S32X2RNG_I(f0, Y0, +8);
                AE_S32X2RNG_IP(c0, Y0, +2 * 8 * 2);

                AE_S32X2RNG_I(f2, Y0, +8);
                AE_S32X2RNG_IP(c2, Y0, +2 * 8 * 2);

                AE_S32X2RNG_I(f1, Y0, +8);
                AE_S32X2RNG_IP(c1, Y0, +2 * 8 * 2);

                AE_S32X2RNG_I(f3, Y0, 8);
                AE_S32X2RNG_XP(c3, Y0, -5 * 8 * 2);

                //
                // Group 1
                //

                d0 = AE_L32X2_I(X0, +8);
                AE_L32X2_IP(a0, X0, +2 * 8 * 2);

                d1 = AE_L32X2_I(X0, +8);
                AE_L32X2_IP(a1, X0, +2 * 8 * 2);

                d2 = AE_L32X2_I(X0, +8);
                AE_L32X2_IP(a2, X0, +2 * 8 * 2);

                d3 = AE_L32X2_I(X0, +8);
                AE_L32X2_XP(a3, X0, +1 * 8 * 2);

                AE_ADDANDSUBRNG32(b0, b2, a0, a2);
                AE_ADDANDSUBRNG32(b1, b3, a3, a1);

                AE_ADDANDSUBRNG32(e0, e2, d0, d2);
                AE_ADDANDSUBRNG32(e1, e3, d3, d1);

                b3 = AE_MUL32JS(b3);
                e3 = AE_MUL32JS(e3);

                AE_ADDANDSUB32S(c0, c2, b0, b1);
                AE_ADDANDSUB32S(c1, c3, b2, b3);

                AE_ADDANDSUB32S(f0, f2, e0, e1);
                AE_ADDANDSUB32S(f1, f3, e2, e3);

                c1 = AE_MULFC32RAS(c1, tw1);
                c2 = AE_MULFC32RAS(c2, tw2);
                c3 = AE_MULFC32RAS(c3, tw3);

                f1 = AE_MULFC32RAS(f1, tw1);
                f2 = AE_MULFC32RAS(f2, tw2);
                f3 = AE_MULFC32RAS(f3, tw3);

                AE_S32X2RNG_I(f0, Y0, +8);
                AE_S32X2RNG_IP(c0, Y0, +2 * 8 * 2);

                AE_S32X2RNG_I(f2, Y0, +8);
                AE_S32X2RNG_IP(c2, Y0, +2 * 8 * 2);

                AE_S32X2RNG_I(f1, Y0, +8);
                AE_S32X2RNG_IP(c1, Y0, +2 * 8 * 2);

                AE_S32X2RNG_I(f3, Y0, 8);
                AE_S32X2RNG_XP(c3, Y0, +1 * 8 * 2);
            }

            AE_CALCRNG3();

            nsa = 3 - RUR_AE_SAR();
        }
    }

    *bexp = nsa;

    return (shiftSum);

} // fft_cplx_inner_32x32_scl2()


/*------------------------------------------------------------------------------
Radix-4 complex-valued FFT of size 2^n, n=4..12.
This function computes the inner stages.

Notes:
1. Inner means  radix-4 FFT FFT stages excluding first and the last stages.
2. At each stage data buffer is automatically downscaled with appropriate shifts
to avoid overflows if necessary.

Precision:
32-bit  complex input/output data, 32-bit complex twiddle factors

Input/Output:
x[2*2*N]        Interleaved stereo complex input data. Real and imaginary data are
interleaved to form single complex data where real data goes first.
y[2*2*N]        Interleaved stereo complex output data. Real and imaginary data are
interleaved to form single complex data where real data goes first.
twd[2*N*3/4]    Complex twiddle factor table, twiddle data is same for both channel.
*bexp           Common block exponent, that is the minimum number of redundant sign
bits over input (output) data
N               FFT Size

Returned value: total right shift amount applied to dynamically scale the data

Restrictions:
x[],y[],twd[] - must not overlap
x[],y[],twd[] - must be aligned on 16-byte boundary
------------------------------------------------------------------------------*/

static int fft_cplx_inner_32x32_scl3(
    complex_fract32 * restrict y,
    complex_fract32 * restrict x,
    const complex_fract32 *    twd,
    int                        twdstep,
    int *                      bexp,
    int                        N
    )
{
    const ae_int32x2 *          X0;
    const ae_int32x2 *          X1;
    const ae_int32x2 *          X2;
    const ae_int32x2 *          X3;
    ae_int32x2 * restrict Y0;
    ae_int32x2 * restrict Y1;
    ae_int32x2 * restrict Y2;
    ae_int32x2 * restrict Y3;
    const ae_f32x2   *          TWD;

    int m, n, logN;
    int stride;
    int shift, shiftSum;

    logN = 30 - NSA(N);

    shiftSum = 0;

    //----------------------------------------------------------------------------
    // Perform second through the second to the last stages.
    {
        complex_fract32 * X;
        complex_fract32 * Y;

        ae_int32x2 a0, a1, a2, a3;
        ae_int32x2 b0, b1, b2, b3;
        ae_f32x2   c0, c1, c2, c3;

        ae_int32x2 d0, d1, d2, d3;
        ae_int32x2 e0, e1, e2, e3;
        ae_f32x2   f0, f1, f2, f3;

        ae_f32x2 tw1, tw2, tw3;

        X = (((logN + 1) & 2) ? x : y);
        Y = (((logN + 1) & 2) ? y : x);

        shiftSum += (shift = 2);

        WUR_AE_SAR(shift);
        twdstep *= 4;
        for (stride = N / 16; stride>4; stride /= 4, twdstep *= 4)
        {

            X0 = (ae_int32x2*)((unsigned int)X + 0 * stride * 8 * 2);
            X1 = (ae_int32x2*)((unsigned int)X + 1 * stride * 8 * 2);
            X2 = (ae_int32x2*)((unsigned int)X + 2 * stride * 8 * 2);
            X3 = (ae_int32x2*)((unsigned int)X + 3 * stride * 8 * 2);

            Y0 = (ae_int32x2*)((unsigned int)Y + 0 * stride * 8 * 2);
            Y1 = (ae_int32x2*)((unsigned int)Y + 1 * stride * 8 * 2);
            Y2 = (ae_int32x2*)((unsigned int)Y + 2 * stride * 8 * 2);
            Y3 = (ae_int32x2*)((unsigned int)Y + 3 * stride * 8 * 2);

            __Pragma("loop_count min=1");
            for (m = 0; m*(4 * stride)<N; m++)
            {
                TWD = (ae_f32x2*)twd;

                __Pragma("ymemory( X0 )");
                __Pragma("ymemory( X1 )");
                __Pragma("ymemory( X2 )");
                __Pragma("ymemory( X3 )");
                __Pragma("loop_count min=4, factor=4");
                for (n = 0; n<stride; n++)
                {
                    tw3 = AE_L32X2_I((ae_int32x2*)TWD, 2 * 8);
                    AE_L32X2X2_XP(tw1, tw2, castxcc(ae_int32x4, TWD), +3 * twdstep * 8);

                    //
                    // Group 0
                    //

                    AE_L32X2X2_XP(a0, d0, castxcc(ae_int32x4, X0), +8 * 2);
                    AE_L32X2X2_XP(a1, d1, castxcc(ae_int32x4, X1), +8 * 2);
                    AE_L32X2X2_XP(a2, d2, castxcc(ae_int32x4, X2), +8 * 2);
                    AE_L32X2X2_XP(a3, d3, castxcc(ae_int32x4, X3), +8 * 2);

                    AE_ADDANDSUBRNG32_H(b0, b2, a0, a2);
                    AE_ADDANDSUBRNG32_H(b1, b3, a3, a1);

                    AE_ADDANDSUBRNG32_L(e0, e2, d0, d2);
                    AE_ADDANDSUBRNG32_L(e1, e3, d3, d1);

                    AE_ADDANDSUB32S(c0, c2, b0, b1);
                    AE_ADDANDSUB32JS(c1, c3, b2, b3);

                    AE_ADDANDSUB32S(f0, f2, e0, e1);
                    AE_ADDANDSUB32JS(f1, f3, e2, e3);

                    c1 = AE_MULFC32RAS(c1, tw1);
                    c2 = AE_MULFC32RAS(c2, tw2);
                    c3 = AE_MULFC32RAS(c3, tw3);

                    f1 = AE_MULFC32RAS(f1, tw1);
                    f2 = AE_MULFC32RAS(f2, tw2);
                    f3 = AE_MULFC32RAS(f3, tw3);

                    AE_S32X2X2RNG_XP(c0, f0, castxcc(ae_int32x4, Y0), +8 * 2);
                    AE_S32X2X2RNG_XP(c2, f2, castxcc(ae_int32x4, Y1), +8 * 2);
                    AE_S32X2X2RNG_XP(c1, f1, castxcc(ae_int32x4, Y2), +8 * 2);
                    AE_S32X2X2RNG_XP(c3, f3, castxcc(ae_int32x4, Y3), +8 * 2);

                }
                X0 += 3 * stride * 2;
                X1 += 3 * stride * 2;
                X2 += 3 * stride * 2;
                X3 += 3 * stride * 2;

                Y0 += 3 * stride * 2;
                Y1 += 3 * stride * 2;
                Y2 += 3 * stride * 2;
                Y3 += 3 * stride * 2;
            }

            {
                // Swap input/output buffers between successive stages.
                complex_fract32 * T = X; X = Y; Y = T;
            }
            shiftSum += shift;
        }
    }

    __Pragma("no_reorder");

    //----------------------------------------------------------------------------
    // Perform the next to last stage.
    {
        ae_int32x2 a0, a1, a2, a3;
        ae_int32x2 b0, b1, b2, b3;
        ae_f32x2   c0, c1, c2, c3;

        ae_int32x2 d0, d1, d2, d3;
        ae_int32x2 e0, e1, e2, e3;
        ae_f32x2   f0, f1, f2, f3;

        TWD = (const ae_f32x2*)((unsigned int)twd + 3 * twdstep * 8);

        X0 = (ae_int32x2*)x;
        Y0 = (ae_int32x2*)y;

        if (stride == 4)
        {
            //
            // Next to last stage for FFT size an even power of two.
            //

            ae_f32x2 tw11, tw12, tw13;
            ae_f32x2 tw21, tw22, tw23;
            ae_f32x2 tw31, tw32, tw33;

            tw21 = ae_f32x2_loadi(TWD, 1 * 8);
            tw31 = ae_f32x2_loadi(TWD, 2 * 8);
            ae_f32x2_loadxp(tw11, TWD, +3 * twdstep * 8);

            tw22 = ae_f32x2_loadi(TWD, 1 * 8);
            tw32 = ae_f32x2_loadi(TWD, 2 * 8);
            ae_f32x2_loadxp(tw12, TWD, +3 * twdstep * 8);

            tw23 = ae_f32x2_loadi(TWD, 1 * 8);
            tw33 = ae_f32x2_loadi(TWD, 2 * 8);
            ae_f32x2_loadxp(tw13, TWD, +3 * twdstep * 8);

            __Pragma("loop_count min=4, factor=4");
            for (n = 0; n<N / (4 * 4); n++)
            {
                //
                // Group 0
                //

                AE_L32X2X2_XP(a0, d0, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a1, d1, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a2, d2, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a3, d3, castxcc(ae_int32x4, X0), -11 * 8 * 2);

                AE_ADDANDSUBRNG32(b0, b2, a0, a2);
                AE_ADDANDSUBRNG32(b1, b3, a3, a1);

                AE_ADDANDSUBRNG32(e0, e2, d0, d2);
                AE_ADDANDSUBRNG32(e1, e3, d3, d1);

                AE_ADDANDSUB32S(c0, c2, b0, b1);
                AE_ADDANDSUB32JS(c1, c3, b2, b3);

                AE_ADDANDSUB32S(f0, f2, e0, e1);
                AE_ADDANDSUB32JS(f1, f3, e2, e3);

                AE_S32X2X2RNG_XP(c0, f0, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c2, f2, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c1, f1, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c3, f3, castxcc(ae_int32x4, Y0), -11 * 8 * 2);

                //
                // Group 1
                //
                AE_L32X2X2_XP(a0, d0, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a1, d1, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a2, d2, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a3, d3, castxcc(ae_int32x4, X0), -11 * 8 * 2);

                AE_ADDANDSUBRNG32(b0, b2, a0, a2);
                AE_ADDANDSUBRNG32(b1, b3, a3, a1);

                AE_ADDANDSUBRNG32(e0, e2, d0, d2);
                AE_ADDANDSUBRNG32(e1, e3, d3, d1);

                AE_ADDANDSUB32S(c0, c2, b0, b1);
                AE_ADDANDSUB32JS(c1, c3, b2, b3);

                AE_ADDANDSUB32S(f0, f2, e0, e1);
                AE_ADDANDSUB32JS(f1, f3, e2, e3);

                c1 = AE_MULFC32RAS(c1, tw11);
                c2 = AE_MULFC32RAS(c2, tw21);
                c3 = AE_MULFC32RAS(c3, tw31);

                f1 = AE_MULFC32RAS(f1, tw11);
                f2 = AE_MULFC32RAS(f2, tw21);
                f3 = AE_MULFC32RAS(f3, tw31);

                AE_S32X2X2RNG_XP(c0, f0, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c2, f2, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c1, f1, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c3, f3, castxcc(ae_int32x4, Y0), -11 * 8 * 2);

                //
                // Group 2
                //
                AE_L32X2X2_XP(a0, d0, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a1, d1, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a2, d2, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a3, d3, castxcc(ae_int32x4, X0), -11 * 8 * 2);

                AE_ADDANDSUBRNG32(b0, b2, a0, a2);
                AE_ADDANDSUBRNG32(b1, b3, a3, a1);

                AE_ADDANDSUBRNG32(e0, e2, d0, d2);
                AE_ADDANDSUBRNG32(e1, e3, d3, d1);

                AE_ADDANDSUB32S(c0, c2, b0, b1);
                AE_ADDANDSUB32JS(c1, c3, b2, b3);

                AE_ADDANDSUB32S(f0, f2, e0, e1);
                AE_ADDANDSUB32JS(f1, f3, e2, e3);

                c1 = AE_MULFC32RAS(c1, tw12);
                c2 = AE_MULFC32RAS(c2, tw22);
                c3 = AE_MULFC32RAS(c3, tw32);

                f1 = AE_MULFC32RAS(f1, tw12);
                f2 = AE_MULFC32RAS(f2, tw22);
                f3 = AE_MULFC32RAS(f3, tw32);

                AE_S32X2X2RNG_XP(c0, f0, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c2, f2, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c1, f1, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c3, f3, castxcc(ae_int32x4, Y0), -11 * 8 * 2);

                //
                // Group 3
                //
                AE_L32X2X2_XP(a0, d0, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a1, d1, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a2, d2, castxcc(ae_int32x4, X0), +4 * 8 * 2);
                AE_L32X2X2_XP(a3, d3, castxcc(ae_int32x4, X0), +1 * 8 * 2);

                AE_ADDANDSUBRNG32(b0, b2, a0, a2);
                AE_ADDANDSUBRNG32(b1, b3, a3, a1);

                AE_ADDANDSUBRNG32(e0, e2, d0, d2);
                AE_ADDANDSUBRNG32(e1, e3, d3, d1);

                AE_ADDANDSUB32S(c0, c2, b0, b1);
                AE_ADDANDSUB32JS(c1, c3, b2, b3);

                AE_ADDANDSUB32S(f0, f2, e0, e1);
                AE_ADDANDSUB32JS(f1, f3, e2, e3);

                c1 = AE_MULFC32RAS(c1, tw13);
                c2 = AE_MULFC32RAS(c2, tw23);
                c3 = AE_MULFC32RAS(c3, tw33);

                f1 = AE_MULFC32RAS(f1, tw13);
                f2 = AE_MULFC32RAS(f2, tw23);
                f3 = AE_MULFC32RAS(f3, tw33);

                AE_S32X2X2RNG_XP(c0, f0, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c2, f2, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c1, f1, castxcc(ae_int32x4, Y0), +4 * 8 * 2);
                AE_S32X2X2RNG_XP(c3, f3, castxcc(ae_int32x4, Y0), +1 * 8 * 2);
            }
        }
        else if (stride == 2)
        {
            //
            // Next to last stage for FFT size an odd power of two.
            //

            ae_f32x2 tw1, tw2, tw3;


            tw1 = ae_f32x2_loadi(TWD, 0 * 8);
            tw2 = ae_f32x2_loadi(TWD, 1 * 8);
            tw3 = ae_f32x2_loadi(TWD, 2 * 8);

            __Pragma("loop_count min=4, factor=4");
            for (n = 0; n<N / (2 * 4); n++)
            {
                //
                // Group 0
                //

                d0 = AE_L32X2_I(X0, +8);
                AE_L32X2_IP(a0, X0, +2 * 8 * 2);

                d1 = AE_L32X2_I(X0, +8);
                AE_L32X2_IP(a1, X0, +2 * 8 * 2);

                d2 = AE_L32X2_I(X0, +8);
                AE_L32X2_IP(a2, X0, +2 * 8 * 2);

                d3 = AE_L32X2_I(X0, +8);
                AE_L32X2_XP(a3, X0, -5 * 8 * 2);

                AE_ADDANDSUBRNG32(b0, b2, a0, a2);
                AE_ADDANDSUBRNG32(b1, b3, a3, a1);

                AE_ADDANDSUBRNG32(e0, e2, d0, d2);
                AE_ADDANDSUBRNG32(e1, e3, d3, d1);

                b3 = AE_MUL32JS(b3);
                e3 = AE_MUL32JS(e3);

                AE_ADDANDSUB32S(c0, c2, b0, b1);
                AE_ADDANDSUB32S(c1, c3, b2, b3);

                AE_ADDANDSUB32S(f0, f2, e0, e1);
                AE_ADDANDSUB32S(f1, f3, e2, e3);

                AE_S32X2RNG_I(f0, Y0, +8);
                AE_S32X2RNG_IP(c0, Y0, +2 * 8 * 2);

                AE_S32X2RNG_I(f2, Y0, +8);
                AE_S32X2RNG_IP(c2, Y0, +2 * 8 * 2);

                AE_S32X2RNG_I(f1, Y0, +8);
                AE_S32X2RNG_IP(c1, Y0, +2 * 8 * 2);

                AE_S32X2RNG_I(f3, Y0, 8);
                AE_S32X2RNG_XP(c3, Y0, -5 * 8 * 2);

                //
                // Group 1
                //

                d0 = AE_L32X2_I(X0, +8);
                AE_L32X2_IP(a0, X0, +2 * 8 * 2);

                d1 = AE_L32X2_I(X0, +8);
                AE_L32X2_IP(a1, X0, +2 * 8 * 2);

                d2 = AE_L32X2_I(X0, +8);
                AE_L32X2_IP(a2, X0, +2 * 8 * 2);

                d3 = AE_L32X2_I(X0, +8);
                AE_L32X2_XP(a3, X0, +1 * 8 * 2);

                AE_ADDANDSUBRNG32(b0, b2, a0, a2);
                AE_ADDANDSUBRNG32(b1, b3, a3, a1);

                AE_ADDANDSUBRNG32(e0, e2, d0, d2);
                AE_ADDANDSUBRNG32(e1, e3, d3, d1);

                b3 = AE_MUL32JS(b3);
                e3 = AE_MUL32JS(e3);

                AE_ADDANDSUB32S(c0, c2, b0, b1);
                AE_ADDANDSUB32S(c1, c3, b2, b3);

                AE_ADDANDSUB32S(f0, f2, e0, e1);
                AE_ADDANDSUB32S(f1, f3, e2, e3);

                c1 = AE_MULFC32RAS(c1, tw1);
                c2 = AE_MULFC32RAS(c2, tw2);
                c3 = AE_MULFC32RAS(c3, tw3);

                f1 = AE_MULFC32RAS(f1, tw1);
                f2 = AE_MULFC32RAS(f2, tw2);
                f3 = AE_MULFC32RAS(f3, tw3);

                AE_S32X2RNG_I(f0, Y0, +8);
                AE_S32X2RNG_IP(c0, Y0, +2 * 8 * 2);

                AE_S32X2RNG_I(f2, Y0, +8);
                AE_S32X2RNG_IP(c2, Y0, +2 * 8 * 2);

                AE_S32X2RNG_I(f1, Y0, +8);
                AE_S32X2RNG_IP(c1, Y0, +2 * 8 * 2);

                AE_S32X2RNG_I(f3, Y0, 8);
                AE_S32X2RNG_XP(c3, Y0, +1 * 8 * 2);
            }
        }
    }
    return (shiftSum);
} // fft_cplx_inner_32x32_scl3()

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
int stereo_fft_cplx32x32_ie(complex_fract32* y, complex_fract32* x, const complex_fract32* twd, int twdstep, int N, int scalingOpt)
{

    const ae_int32x2 *          X0;
    const ae_int32x2 *          X1;
    const ae_int32x2 *          X2;
    const ae_int32x2 *          X3;
    ae_int32x2 * restrict Y0;
    ae_int32x2 * restrict Y1;
    ae_int32x2 * restrict Y2;
    ae_int32x2 * restrict Y3;
    const ae_int32x4 *   restrict pX;
    const ae_f32x2   *          TWD;


    int n, logN;
    int isFirstInplace;
    int nsa, shift, shiftSum;
    int minScaling;
    logN = 30 - NSA(N);

    // All the stages from the second to the last are performed out-of-place. For
    // the last stage the reason is the bit reversal permutation, while for other
    // stages this just helps to avoid processor stalls due to memory conflicts.
    // Thus, we choose the first stage to be done either in-place or out-of-place
    // so that the whole FFT would finish in y[].
    isFirstInplace = !((logN + 1) & 2);

    if (scalingOpt == 2)
    {
        /*
        * autoscaling case
        */
        int bexp, shiftl, shiftr;
        ae_int32x2 scl;
        /* Calculate the exponent to prescale the input data */
        complex_fract32 * Y;
        int n;

        ae_int32x2 a0, a1, a2, a3;
        ae_int32x2 b0, b1, b2, b3;
        ae_int32x2 c0, c1, c2, c3;
        ae_int32x2 d0, d1, d2, d3;

        ae_f32x2 tw1, tw2, tw3;
        ae_int32x2 x0, x1, x2, x3;
        ae_int16x4 nsa0, nsa1;
        minScaling  = 0; 
        pX = (const ae_int32x4 *)x;
        nsa0 = 31; nsa1 = 31;
        NASSERT((N & 3) == 0);
        __Pragma("loop_count min=1");
        for (n = 0; n<(N >> 1); n++)
        {
            AE_L32X2X2_IP(x0, x1, pX, sizeof(ae_int32x4));
            AE_L32X2X2_IP(x2, x3, pX, sizeof(ae_int32x4));
            nsa0 = AE_MIN16(nsa0, AE_NSA32X4(x0, x1));
            nsa1 = AE_MIN16(nsa1, AE_NSA32X4(x2, x3));
        }
        bexp = AE_RMIN16X4(AE_MIN16(nsa0, nsa1));

        shiftSum = shift = 3 - bexp;
        shiftl = XT_MAX(0, -shift);
        shiftr = XT_MAX(0, shift);
        scl = 1 << shiftl;
        WUR_AE_SAR(shiftr);
        NASSERT(shift>-32 && shift<4);

        X0 = (ae_int32x2*)((unsigned int)x + 0 * N / 4 * 8 * 2);
        X1 = (ae_int32x2*)((unsigned int)x + 1 * N / 4 * 8 * 2);
        X2 = (ae_int32x2*)((unsigned int)x + 2 * N / 4 * 8 * 2);
        X3 = (ae_int32x2*)((unsigned int)x + 3 * N / 4 * 8 * 2);

        Y = (isFirstInplace ? x : y);

        Y0 = (ae_int32x2*)((unsigned int)Y + 0 * N / 4 * 8 * 2);
        Y1 = (ae_int32x2*)((unsigned int)Y + 1 * N / 4 * 8 * 2);
        Y2 = (ae_int32x2*)((unsigned int)Y + 2 * N / 4 * 8 * 2);
        Y3 = (ae_int32x2*)((unsigned int)Y + 3 * N / 4 * 8 * 2);

        TWD = (ae_f32x2*)twd;

        __Pragma("ymemory( X0 )");
        __Pragma("ymemory( X1 )");
        __Pragma("ymemory( X2 )");
        __Pragma("ymemory( X3 )");
        __Pragma("loop_count min=4, factor=4");

        //----------------------------------------------------------------------------
        // Perform the first stage. We use DIF, all permutations are deferred
        // until the last stage. 
        for (n = 0; n<N / 4; n++)
        {
            AE_L32X2X2_IP(a0, c0, castxcc(ae_int32x4, X0), +16);
            AE_L32X2X2_IP(a1, c1, castxcc(ae_int32x4, X1), +16);
            AE_L32X2X2_IP(a2, c2, castxcc(ae_int32x4, X2), +16);
            AE_L32X2X2_IP(a3, c3, castxcc(ae_int32x4, X3), +16);

            a0 = AE_MULP32X2(a0, scl);
            a1 = AE_MULP32X2(a1, scl);
            a2 = AE_MULP32X2(a2, scl);
            a3 = AE_MULP32X2(a3, scl);
            c0 = AE_MULP32X2(c0, scl);
            c1 = AE_MULP32X2(c1, scl);
            c2 = AE_MULP32X2(c2, scl);
            c3 = AE_MULP32X2(c3, scl);

            ae_f32x2_loadip(tw1, TWD, sizeof(complex_fract32));
            ae_f32x2_loadip(tw2, TWD, sizeof(complex_fract32));
            ae_f32x2_loadxp(tw3, TWD, (3 * twdstep - 2)*sizeof(complex_fract32));

            AE_ADDANDSUBRNG32(b0, b2, a0, a2);
            AE_ADDANDSUBRNG32(b1, b3, a3, a1);

            AE_ADDANDSUBRNG32(d0, d2, c0, c2);
            AE_ADDANDSUBRNG32(d1, d3, c3, c1);

            AE_ADDANDSUB32S(a0, a2, b0, b1);
            AE_ADDANDSUB32JS(a1, a3, b2, b3);

            AE_ADDANDSUB32S(c0, c2, d0, d1);
            AE_ADDANDSUB32JS(c1, c3, d2, d3);

            a1 = AE_MULFC32RAS(a1, tw1);
            a2 = AE_MULFC32RAS(a2, tw2);
            a3 = AE_MULFC32RAS(a3, tw3);

            c1 = AE_MULFC32RAS(c1, tw1);
            c2 = AE_MULFC32RAS(c2, tw2);
            c3 = AE_MULFC32RAS(c3, tw3);

            // Two middle quartiles are swapped to use bit reversal instead of
            // digit reversal at the last stage.

            AE_S32X2X2RNG_IP(a0, c0, castxcc(ae_int32x4, Y0), +16);
            AE_S32X2X2RNG_IP(a2, c2, castxcc(ae_int32x4, Y1), +16);
            AE_S32X2X2RNG_IP(a1, c1, castxcc(ae_int32x4, Y2), +16);
            AE_S32X2X2RNG_IP(a3, c3, castxcc(ae_int32x4, Y3), +16);
        }

        //----------------------------------------------------------------------------
        // Perform second through the next to last stages.
        {
            int shift1, shift2;
            AE_CALCRNG32(shift1, shift2, 0, 3);
            (void)shift2;
            nsa = 3 - shift1;
            shiftSum += fft_cplx_inner_32x32_scl2(x, y, twd, twdstep, &nsa, N);
        }
    }
    else  /* if (scalingOpt == 2) */
    {
        /*
        * static scaling case
        */
        complex_fract32 * Y;
        int n;

        ae_int32x2 a0, a1, a2, a3;
        ae_int32x2 b0, b1, b2, b3;
        ae_int32x2 c0, c1, c2, c3;
        ae_int32x2 d0, d1, d2, d3;

        ae_f32x2 tw1, tw2, tw3;

        minScaling = 2;
        shiftSum = 3;
        WUR_AE_SAR(shiftSum);

        X0 = (ae_int32x2*)((unsigned int)x + 0 * N / 4 * 8 * 2);
        X1 = (ae_int32x2*)((unsigned int)x + 1 * N / 4 * 8 * 2);
        X2 = (ae_int32x2*)((unsigned int)x + 2 * N / 4 * 8 * 2);
        X3 = (ae_int32x2*)((unsigned int)x + 3 * N / 4 * 8 * 2);

        Y = (isFirstInplace ? x : y);

        Y0 = (ae_int32x2*)((unsigned int)Y + 0 * N / 4 * 8 * 2);
        Y1 = (ae_int32x2*)((unsigned int)Y + 1 * N / 4 * 8 * 2);
        Y2 = (ae_int32x2*)((unsigned int)Y + 2 * N / 4 * 8 * 2);
        Y3 = (ae_int32x2*)((unsigned int)Y + 3 * N / 4 * 8 * 2);

        TWD = (ae_f32x2*)twd;

        __Pragma("ymemory( X0 )");
        __Pragma("ymemory( X1 )");
        __Pragma("ymemory( X2 )");
        __Pragma("ymemory( X3 )");
        __Pragma("loop_count min=4, factor=4");

        //----------------------------------------------------------------------------
        // Perform the first stage. We use DIF, all permutations are deferred
        // until the last stage. 
        for (n = 0; n<N / 4; n++)
        {
            AE_L32X2X2_IP(a0, c0, castxcc(ae_int32x4, X0), +16);
            AE_L32X2X2_IP(a1, c1, castxcc(ae_int32x4, X1), +16);
            AE_L32X2X2_IP(a2, c2, castxcc(ae_int32x4, X2), +16);
            AE_L32X2X2_IP(a3, c3, castxcc(ae_int32x4, X3), +16);

            ae_f32x2_loadip(tw1, TWD, sizeof(complex_fract32));
            ae_f32x2_loadip(tw2, TWD, sizeof(complex_fract32));
            ae_f32x2_loadxp(tw3, TWD, (3 * twdstep - 2)*sizeof(complex_fract32));

            AE_ADDANDSUBRNG32(b0, b2, a0, a2);
            AE_ADDANDSUBRNG32(b1, b3, a3, a1);

            AE_ADDANDSUBRNG32(d0, d2, c0, c2);
            AE_ADDANDSUBRNG32(d1, d3, c3, c1);

            AE_ADDANDSUB32S(a0, a2, b0, b1);
            AE_ADDANDSUB32JS(a1, a3, b2, b3);

            AE_ADDANDSUB32S(c0, c2, d0, d1);
            AE_ADDANDSUB32JS(c1, c3, d2, d3);

            a1 = AE_MULFC32RAS(a1, tw1);
            a2 = AE_MULFC32RAS(a2, tw2);
            a3 = AE_MULFC32RAS(a3, tw3);

            c1 = AE_MULFC32RAS(c1, tw1);
            c2 = AE_MULFC32RAS(c2, tw2);
            c3 = AE_MULFC32RAS(c3, tw3);

            // Two middle quartiles are swapped to use bit reversal instead of
            // digit reversal at the last stage.

            AE_S32X2X2RNG_IP(a0, c0, castxcc(ae_int32x4, Y0), +16);
            AE_S32X2X2RNG_IP(a2, c2, castxcc(ae_int32x4, Y1), +16);
            AE_S32X2X2RNG_IP(a1, c1, castxcc(ae_int32x4, Y2), +16);
            AE_S32X2X2RNG_IP(a3, c3, castxcc(ae_int32x4, Y3), +16);
        }
        nsa = 0;
        shiftSum += fft_cplx_inner_32x32_scl3(x, y, twd, twdstep, &nsa, N);
    }       /*if (scalingOpt == 2)*/



    //----------------------------------------------------------------------------
    // Perform the last stage.
    if (!(logN & 1))
    {
        //
        // Last stage for FFT size an even power of two: radix-4.
        //
        ae_int32x2 a0, a1, a2, a3;
        ae_int32x2 b0, b1, b2, b3;
        ae_f32x2   c0, c1, c2, c3;

        ae_int32x2 d0, d1, d2, d3;
        ae_int32x2 e0, e1, e2, e3;
        ae_f32x2   f0, f1, f2, f3;

        unsigned int ix = 0;

        shiftSum += (shift = XT_MAX(minScaling, 2 - nsa));

        WUR_AE_SAR(shift);

        X0 = (ae_int32x2*)x;

        Y0 = (ae_int32x2*)((unsigned int)y + 0 * N / 4 * 8 * 2);
        Y1 = (ae_int32x2*)((unsigned int)y + 1 * N / 4 * 8 * 2);
        Y2 = (ae_int32x2*)((unsigned int)y + 2 * N / 4 * 8 * 2);
        Y3 = (ae_int32x2*)((unsigned int)y + 3 * N / 4 * 8 * 2);

        __Pragma("loop_count min=4, factor=4");
        for (n = 0; n<N / 4; n++)
        {
            AE_L32X2X2_IP(a0, d0, castxcc(ae_int32x4, X0), +8 * 2);
            AE_L32X2X2_IP(a1, d1, castxcc(ae_int32x4, X0), +8 * 2);
            AE_L32X2X2_IP(a2, d2, castxcc(ae_int32x4, X0), +8 * 2);
            AE_L32X2X2_IP(a3, d3, castxcc(ae_int32x4, X0), +8 * 2);

            AE_ADDANDSUBRNG32(b0, b2, a0, a2);
            AE_ADDANDSUBRNG32(b1, b3, a3, a1);

            AE_ADDANDSUBRNG32(e0, e2, d0, d2);
            AE_ADDANDSUBRNG32(e1, e3, d3, d1);

            AE_ADDANDSUB32S(c0, c2, b0, b1);
            AE_ADDANDSUB32JS(c1, c3, b2, b3);

            AE_ADDANDSUB32S(f0, f2, e0, e1);
            AE_ADDANDSUB32JS(f1, f3, e2, e3);

            AE_S32X2X2_X(c0, f0, (ae_int32x4*)Y0, ix);
            AE_S32X2X2_X(c1, f1, (ae_int32x4*)Y1, ix);
            AE_S32X2X2_X(c2, f2, (ae_int32x4*)Y2, ix);
            AE_S32X2X2_X(c3, f3, (ae_int32x4*)Y3, ix);

            ix = AE_ADDBRBA32(ix, 1UL << (32 - (3 + 1 + logN - 2)));
        }
    }
    else
    {
        //
        // Last stage for FFT size an odd power of two: radix-2.
        //
        ae_int32x2 a0, a1, b0, b1;
        ae_int32x2 d0, d1, e0, e1;
        unsigned int ix = 0;
        ASSERT(minScaling == 0 || minScaling == 2);

        minScaling = XT_MAX(minScaling-1, 0); 
        shiftSum += (shift = XT_MAX(minScaling, 1 - nsa));
        WUR_AE_SAR(shift);

        X0 = (ae_int32x2*)x;
        Y0 = (ae_int32x2*)((unsigned int)y + 0 * N / 2 * 8 * 2);
        Y1 = (ae_int32x2*)((unsigned int)y + 1 * N / 2 * 8 * 2);

        __Pragma("loop_count min=16, factor=16");
        for (n = 0; n<N / 2; n++)
        {
            d0 = AE_L32X2_I(X0, +8);
            AE_L32X2_IP(a0, X0, +8 * 2);

            d1 = AE_L32X2_I(X0, +8);
            AE_L32X2_IP(a1, X0, +8 * 2);

            AE_ADDANDSUBRNG32(b0, b1, a0, a1);
            AE_ADDANDSUBRNG32(e0, e1, d0, d1);

            AE_S32X2RNG_X(e0, Y0, ix * 2 + 8);
            AE_S32X2RNG_X(b0, Y0, ix * 2);

            AE_S32X2RNG_X(e1, Y1, ix * 2 + 8);
            AE_S32X2RNG_X(b1, Y1, ix * 2);

            ix = AE_ADDBRBA32(ix, 1UL << (32 - (3 + logN - 1)));
        }
    }
    ASSERT(scalingOpt == 2 || shiftSum == (logN +1 )); 
    return (shiftSum);
} /* stereo_fft_cplx32x32_ie() */
