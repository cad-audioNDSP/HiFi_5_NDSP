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
    Reference C code
	Integrit, 2006-2019
*/
#include "NatureDSP_Signal_fft.h"
#include "NatureDSP_Signal_vector.h"
#include "common.h"
/* Twiddle factor tables and FFT descriptor structure. */
#include "fft_x16_common.h"
#include "fft_16x16_stages.h"

#if 0
#include "baseop.h"
static uint16_t _HI(uint32_t x) { return (uint16_t)(x >> 16); }
static uint16_t _LO(uint32_t x) { return (uint16_t)(x); }
static uint32_t _PACK(uint16_t hi, uint16_t lo) { return (((uint32_t)hi) << 16) | (lo); }

static uint32_t _ADD2S(uint32_t x, uint32_t y)
{
  return _PACK(S_add_ss(_HI(x), _HI(y)), S_add_ss(_LO(x), _LO(y)));
}
static uint32_t _SUB2S(uint32_t x, uint32_t y)
{
  return _PACK(S_sub_ss(_HI(x), _HI(y)), S_sub_ss(_LO(x), _LO(y)));
}
/* Two-way bidirectional arithmetic (sign-extending) right shift */
static uint32_t _SHR2R_BIDIR(uint32_t x, int rsh)
{
    int16_t hi = _HI(x), lo = _LO(x);
    if (rsh>0)
    {
        int32_t _hi = (int32_t)hi + (1L<<(rsh-1));
        int32_t _lo = (int32_t)lo + (1L<<(rsh-1));
        _hi >>= rsh; _lo >>= rsh;
        hi = S_sature_l(_hi); lo = S_sature_l(_lo);
    }
    else
    {
        hi <<= -rsh; lo <<= -rsh;
    }
    return _PACK(hi, lo);
}

static void _cmult16x16(uint32_t *result, uint32_t *x, uint32_t *y)
{
    int16_t re = S_extract_l(L_sub_ll(L_mpy_ss(_LO(*x), _LO(*y)), L_mpy_ss(_HI(*x), _HI(*y))));
    int16_t im = S_extract_l(L_add_ll(L_mpy_ss(_LO(*x), _HI(*y)), L_mpy_ss(_HI(*x), _LO(*y))));
    *result = _PACK(im, re);
}

static void _DFT4X1(uint32_t *x0, uint32_t *x1, uint32_t *x2, uint32_t *x3)
{
    uint32_t B0, B1, B2, B3;
    uint32_t C0, C1, C2, C3;

    B0 = _ADD2S(*x0, *x2);
    B1 = _ADD2S(*x1, *x3);
    B2 = _SUB2S(*x0, *x2);
    B3 = _SUB2S(*x1, *x3);
    B3 = _PACK(_LO(B3), S_neg_s((int16_t)_HI(B3)));
    C0 = _ADD2S(B0, B1);
    C1 = _SUB2S(B0, B1);
    C2 = _ADD2S(B2, B3);
    C3 = _SUB2S(B2, B3);

    *x0 = C0;
    *x2 = C1;
    *x3 = C2;
    *x1 = C3;
}

int fft_16x16_stage_inner_scl2_DFT4(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    static const int sar[4]={0,0x102,0x183,0x285};  //scaling for DFT4
    const int R = 4; // stage radix
    const int _v=v[0];
    int i, j;
    const int stride = N / R;
    const int Ninner=stride/_v;

    uint32_t * restrict px = (uint32_t *)x;
    uint32_t * restrict py = (uint32_t *)y;
    uint32_t * restrict ptw = (uint32_t *)tw;
    const ae_int16x8 * restrict pX;
          ae_int16x8 * restrict pY;
    const ae_int32   * restrict pTwd;

    const int shift = XT_MAX(0, 3 - *bexp);
    /* Detect last stage */
    NASSERT(shift>=0 && shift<=3); 
    NASSERT(v[0]!=stride);  // not last
    NASSERT(v[0]!=1);       // not first
    NASSERT((N/R)%4==0);

    for (j = 0; j < _v; j++)
    {
        for (i = 0; i < Ninner; i++)
        {
            uint32_t x0, x1, x2, x3;
            uint32_t t1, t2, t3;
            t1=ptw[3 * tw_step* i+0];
            t2=ptw[3 * tw_step* i+1];
            t3=ptw[3 * tw_step* i+2];
            x0 = px[j + _v * i + 0];
            x1 = px[j + _v * i + 1 * stride];
            x2 = px[j + _v * i + 2 * stride];
            x3 = px[j + _v * i + 3 * stride];

            x0 = _SHR2R_BIDIR(x0, shift);
            x1 = _SHR2R_BIDIR(x1, shift);
            x2 = _SHR2R_BIDIR(x2, shift);
            x3 = _SHR2R_BIDIR(x3, shift);

            _DFT4X1(&x0, &x1, &x2, &x3);

            _cmult16x16(&x1, &x1, &t1);
            _cmult16x16(&x2, &x2, &t2);
            _cmult16x16(&x3, &x3, &t3);
            py[j + 4 * _v * i + 0] = x0;
            py[j + 4 * _v * i + 1 * v[0]] = x1;
            py[j + 4 * _v * i + 2 * v[0]] = x2;
            py[j + 4 * _v * i + 3 * v[0]] = x3;
        }
    }

    *bexp = vec_bexp16(y, 2*N) - 16;
    *v *= R;
    return shift;
}
#endif


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

int fft_16x16_stage_inner_scl2_DFT4(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    static const int sar[4]={0,0x102,0x183,0x285};  //scaling for DFT4
    const int R = 4; // stage radix
    const int _v=v[0];
    int i, j;
    const int stride = N / R;
    const int Ninner=stride/_v;

    const ae_int16x8 * restrict pX;
          ae_int16x8 * restrict pY;
    const ae_int32   * restrict pTwd;

    const int shift = XT_MAX(0, 3 - *bexp);
    /* Detect last stage */
    NASSERT(shift>=0 && shift<=3); 
    NASSERT(_v!=stride);  // not last
    NASSERT(_v!=1);       // not first
    NASSERT(_v%4==0);     // not second after R3/R5 
    NASSERT((N/R)%4==0);
    WUR_AE_SAR(sar[shift]);

    // Ninner is a multiple of 4 - everything is aligned 
    for (j = 0; j < _v; j+=4)
    {
        pX  = (const ae_int16x8 *)(x+j*2);
        pY  = (      ae_int16x8 *)(y+j*2);
        pTwd= (const ae_int32 *)tw;
        for (i = 0; i < Ninner; i++)
        {
            ae_int16x4 x00,x01,x10,x11,x20,x21,x30,x31;
            ae_int16x4 y00,y01,y10,y11,y20,y21,y30,y31;
            ae_int16x4 tw1,tw2,tw3;
            ae_int32x2 t;
            tw2 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(AE_L32_I(pTwd, 1*sizeof(complex_fract16))));
            tw3 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(AE_L32_I(pTwd, 2*sizeof(complex_fract16))));
            AE_L32_XP(t,pTwd, 3 * tw_step*sizeof(complex_fract16));
            tw1 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t));

            AE_L16X4X2_X (x10,x11,pX, 1 * stride*sizeof(complex_fract16));
            AE_L16X4X2_X (x20,x21,pX, 2 * stride*sizeof(complex_fract16));
            AE_L16X4X2_X (x30,x31,pX, 3 * stride*sizeof(complex_fract16));
            AE_L16X4X2_XP(x00,x01,pX,_v * sizeof(complex_fract16));

            DFT4XI2(x00,x10,x20,x30);
            DFT4XI2(x01,x11,x21,x31);
            y00 = x00;
            y10 = AE_MULFC16RAS(x10,tw1);
            y20 = AE_MULFC16RAS(x20,tw2);
            y30 = AE_MULFC16RAS(x30,tw3);
            y01 = x01;
            y11 = AE_MULFC16RAS(x11,tw1);
            y21 = AE_MULFC16RAS(x21,tw2);
            y31 = AE_MULFC16RAS(x31,tw3);
            AE_S16X4X2RNG_XP(y00,y01,pY,_v*sizeof(complex_fract16));
            AE_S16X4X2RNG_XP(y10,y11,pY,_v*sizeof(complex_fract16));
            AE_S16X4X2RNG_XP(y20,y21,pY,_v*sizeof(complex_fract16));
            AE_S16X4X2RNG_XP(y30,y31,pY,_v*sizeof(complex_fract16));
        }
    }
    // update scaling
    {int a, b; AE_CALCRNG16(a, b, 0, 3);*bexp = 3 - a;(void)b; }
    *v *= R;
    return shift;
}
