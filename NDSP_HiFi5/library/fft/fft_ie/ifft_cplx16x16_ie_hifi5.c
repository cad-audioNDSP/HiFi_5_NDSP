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

#define __DFT4XI2(x0, x1, x2, x3, shift)/* output x0, x3, x1, x2*/\
{\
    ae_int16x4 t0, t1, t2, t3;                                        \
    xtbool4 mask = 0x5;                                               \
    x0 = AE_SRAA16RS(x0, shift);                                      \
    x1 = AE_SRAA16RS(x1, shift);                                      \
    x2 = AE_SRAA16RS(x2, shift);                                      \
    x3 = AE_SRAA16RS(x3, shift);                                      \
    t0 = AE_ADD16S(x0, x2); t2 = AE_SUB16S(x0, x2);                   \
    t1 = AE_ADD16S(x1, x3); t3 = AE_SUB16S(x1, x3);                   \
    x0 = t0; x1 = t1; x2 = t2; x3 = t3;                               \
    AE_MOVT16X4(x3, AE_NEG16S(x3), mask); /*  x3 = AE_CONJ16S(x3); */ \
    x3 = AE_SEL16_2301(x3, x3);                                       \
    t0 = AE_ADD16S(x0, x1); t1 = AE_SUB16S(x0, x1);                   \
    t2 = AE_ADD16S(x2, x3); t3 = AE_SUB16S(x2, x3);                   \
    x0 = t0; x1 = t1; x2 = t2; x3 = t3;                               \
    t1 = x3;                                                          \
    t2 = x1;                                                          \
    t3 = x2;                                                          \
    x1 = t1;                                                          \
    x2 = t2;                                                          \
    x3 = t3;                                                          \
}
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

static int stage_first_iDFT4_16x16_ie(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    int i, shift;
    const int R = 4; // stage radix
    const int stride = N / R;
    const int _v = v[0];
    const int tw_inc0 = N / _v / 4 * tw_step* sizeof(complex_fract16);

    ae_int16x4 * restrict _py;
    ae_int16x4 * restrict _px;
    ae_int32   * restrict ptw;

    ae_int32x2 t10, t20, t30, t11, t21, t31;
    //   ae_int16x4 acc16 = AE_MOVINT16X4_FROMF32X2(AE_MOVI(0));


    _py = (ae_int16x4 *)y;
    _px = (ae_int16x4 *)x;
    ptw = (ae_int32 *)((uintptr_t)tw);

    {
        /* Reset RANGE register */
        int a, b;
        AE_CALCRNG16(a, b, 0, 3); (void)b; (void)a;
    }


    shift = 3 - *bexp;
#if __USE_AE_SLAA16__
    /* 9 cycles unroll 1*/
    SetDFT4_Scaling(3);
#else
    /* 10 cycles unroll 1*/
    int r, l;
    l = XT_MAX(-shift, 0);
    r = XT_MAX(shift, 0);
    ae_int16x4 s = 1 << l;
    SetDFT4_Scaling(r);
#endif
    __Pragma("loop_count min=1");
    for (i = 0; i < N / 2 / R; i++)
    {
        ae_int32x2 x0, x1, x2, x3;
        ae_int16x4 _x0, _x1, _x2, _x3;
        ae_int16x4 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5;
        ae_int16x4 tw1, tw2, tw3;
        ae_int16x4 z0, z1, z2, z3;

        AE_L32_XP(t10, ptw, tw_inc0);
        AE_L32_XP(t20, ptw, tw_inc0);
        AE_L32_XP(t30, ptw, tw_step* sizeof(complex_fract16)-2 * tw_inc0);
        AE_L32_XP(t11, ptw, tw_inc0);
        AE_L32_XP(t21, ptw, tw_inc0);
        AE_L32_XP(t31, ptw, tw_step* sizeof(complex_fract16)-2 * tw_inc0);

        tmp0 = (AE_MOVINT16X4_FROMF32X2(t10));
        tmp1 = (AE_MOVINT16X4_FROMF32X2(t20));
        tmp2 = (AE_MOVINT16X4_FROMF32X2(t30));
        tmp3 = (AE_MOVINT16X4_FROMF32X2(t11));
        tmp4 = (AE_MOVINT16X4_FROMF32X2(t21));
        tmp5 = (AE_MOVINT16X4_FROMF32X2(t31));

        tw1 = AE_SEL16_7632(tmp0, tmp3);
        tw2 = AE_SEL16_7632(tmp1, tmp4);
        tw3 = AE_SEL16_7632(tmp2, tmp5);

        AE_L32X2_XP(x0, castxcc(ae_int32x2, _px), stride*sizeof(complex_fract16));
        AE_L32X2_XP(x1, castxcc(ae_int32x2, _px), stride*sizeof(complex_fract16));
        AE_L32X2_XP(x2, castxcc(ae_int32x2, _px), stride*sizeof(complex_fract16));
        AE_L32X2_XP(x3, castxcc(ae_int32x2, _px), sizeof(*_px) - 3 * stride*sizeof(complex_fract16));
        _x0 = AE_MOVINT16X4_FROMF32X2(x0);
        _x1 = AE_MOVINT16X4_FROMF32X2(x1);
        _x2 = AE_MOVINT16X4_FROMF32X2(x2);
        _x3 = AE_MOVINT16X4_FROMF32X2(x3);

#ifdef __USE_AE_SLAA16__
        _x0 = AE_SLAA16(_x0, bexp[0]);
        _x1 = AE_SLAA16(_x1, bexp[0]);
        _x2 = AE_SLAA16(_x2, bexp[0]);
        _x3 = AE_SLAA16(_x3, bexp[0]);
#else
        _x0 = AE_MULP16X16X4S(_x0, s);
        _x1 = AE_MULP16X16X4S(_x1, s);
        _x2 = AE_MULP16X16X4S(_x2, s);
        _x3 = AE_MULP16X16X4S(_x3, s);
#endif

        DFT4XI2(_x0, _x1, _x2, _x3);

        _x1 = AE_MULFC16RAS(_x1, tw1);
        _x2 = AE_MULFC16RAS(_x2, tw2);
        _x3 = AE_MULFC16RAS(_x3, tw3);

        z0 = AE_SEL16_7632(_x0, _x1);
        z1 = AE_SEL16_7632(_x2, _x3);
        z2 = AE_SEL16_5410(_x0, _x1);
        z3 = AE_SEL16_5410(_x2, _x3);

        AE_S16X4RNG_IP(z0, _py, sizeof(*_py));
        AE_S16X4RNG_IP(z1, _py, sizeof(*_py));
        AE_S16X4RNG_IP(z2, _py, sizeof(*_py));
        AE_S16X4RNG_IP(z3, _py, sizeof(*_py));
    }
    {
        int a, b;
        AE_CALCRNG16(a, b, 0, 3); (void)b;
        *bexp = 3 - a;
    }

    *v *= R;
    return shift;
} //stage_first_iDFT4_16x16_ie

static int stage_last_iDFT4_16x16_ie(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    int i, shift;
    const int R = 4; // stage radix
    const int min_shift = 2;
    const int stride = N / R;
    ae_int32x2 * restrict _py;
    ae_int16x4 * restrict _px;

    shift = XT_MAX(0, min_shift - *bexp);
    {
        /* Reset RANGE register */
        int a, b;
        AE_CALCRNG16(a, b, 0, 3); (void)b; (void)a;
    }
    SetDFT4_Scaling(shift);

    _py = (ae_int32x2 *)y;
    _px = (ae_int16x4 *)x;

    // Last phase, without twiddles
    __Pragma("loop_count min=2");
    for (i = 0; i < N / 8; i++)
    {
        ae_int16x4 _x0, _x1, _x2, _x3;

        AE_L16X4_XP(_x0, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x1, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x2, _px, stride*sizeof(complex_fract16));
        AE_L16X4_XP(_x3, _px, sizeof(*_px) - 3 * stride*sizeof(complex_fract16));

        DFT4XI2(_x0, _x1, _x2, _x3);

        AE_S32X2_XP(AE_MOVINT32X2_FROMINT16X4(_x0), _py, stride*sizeof(complex_fract16));
        AE_S32X2_XP(AE_MOVINT32X2_FROMINT16X4(_x1), _py, stride*sizeof(complex_fract16));
        AE_S32X2_XP(AE_MOVINT32X2_FROMINT16X4(_x2), _py, stride*sizeof(complex_fract16));
        AE_S32X2_XP(AE_MOVINT32X2_FROMINT16X4(_x3), _py, sizeof(*_py) - 3 * stride*sizeof(complex_fract16));
    }
    return shift;
} //stage_last_iDFT4_16x16_ie

extern int stage_inner_DFT4_16x16_ie(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp, int isStereo);


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

int ifft_cplx16x16_ie(complex_fract16* y, complex_fract16* x, const complex_fract16* twd, int twdstep, int N, int scalingOpt)
{
    int bexp, shift = 0;
    int v = 1;
    complex_fract16 *pdest = y;
    int log2N = 30 - NSA(N);
    ae_int16x4 * restrict px;

    NASSERT_ALIGN16(twd);
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(x != y);
    NASSERT(scalingOpt == 2);
    NASSERT(N == 128 || N == 256 || N == 512 || N == 1024);

    {
        int i;
        ae_int16x4 acc = AE_MOVINT16X4_FROMINT32X2(AE_MOVI(0)), tmp;
        px = (ae_int16x4*)x;

        __Pragma("loop_count min=2");
        for (i = 0; i < (N >> 1); i++)
        {
            AE_L16X4_IP(tmp, px, sizeof(*px));
            tmp = AE_ABS16S(tmp);
            acc = AE_OR16(acc, tmp);
        }
        acc = AE_OR16(acc, AE_SEL16_5432(acc, acc));
        acc = AE_OR16(acc, AE_SHORTSWAP(acc));

        i = AE_MOVAD16_0(acc);
        bexp = NSA(i) - 16;
        XT_MOVEQZ(bexp, 0, i);
    }

    shift += stage_first_iDFT4_16x16_ie((const int16_t*)twd, (int16_t*)x, (int16_t*)y, N, &v, twdstep, &bexp);
    SWAP_PTR(x, y);
    log2N -= 2;
    twdstep *= 4;

    while (log2N >=4)
    {
        shift += stage_inner_DFT4_16x16_ie((const int16_t*)twd, (int16_t*)x, (int16_t*)y, N, &v, twdstep, &bexp, 0);
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
    if (log2N & 1)
    {
        shift += ifft_16x16_stage_last_scl2_DFT8(NULL, (const int16_t*)x, (int16_t*)y, N, &v, 0, &bexp);
    }
    else
    {
        shift += stage_last_iDFT4_16x16_ie(NULL, (int16_t*)x, (int16_t*)y, N, &v, 0, &bexp);
    }
    return shift;
} /* ifft_cplx16x16_ie() */
