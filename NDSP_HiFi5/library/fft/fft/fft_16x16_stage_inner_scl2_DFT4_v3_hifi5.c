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

int fft_16x16_stage_inner_scl2_DFT4_v3(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
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

int fft_16x16_stage_inner_scl2_DFT4_v3(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    static const int sar[4]={0,0x102,0x183,0x285};  //scaling for DFT4
    int i;
    const int stride = N>>2;
    const int N12=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),178956971));//N/12
    const int Ninner=N12;
    const ae_int16x4 * restrict pX0=(       ae_int16x4 *)(x);
    const ae_int16x4 * restrict pX1=(       ae_int16x4 *)(x+2*1*stride);
    const ae_int16x4 * restrict pX2=(       ae_int16x4 *)(x+2*2*stride);
    const ae_int16x4 * restrict pX3=(       ae_int16x4 *)(x+2*3*stride);
          ae_int16x8 * restrict pY =(       ae_int16x8 *)y;
    const ae_int32   * restrict pTwd=(const ae_int32   *)tw;
    const int shift = XT_MAX(0, 3 - *bexp);
    /* Detect last stage */
    NASSERT(shift>=0 && shift<=3); 
    NASSERT(3!=stride);  // not last
    NASSERT(N%16==0);
    NASSERT(Ninner%2==0);
    WUR_AE_SAR(sar[shift]);
    for (i = 0; i < Ninner; i+=2)
    {
        ae_int16x4 x001,x023,x045;
        ae_int16x4 x101,x123,x145;
        ae_int16x4 x201,x223,x245;
        ae_int16x4 x301,x323,x345;
        ae_int16x4 tw10,tw11;
        ae_int16x4 tw20,tw21;
        ae_int16x4 tw30,tw31;
        ae_int32x2 t;

        tw20 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(AE_L32_I(pTwd, 1*sizeof(complex_fract16))));
        tw30 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(AE_L32_I(pTwd, 2*sizeof(complex_fract16))));
        AE_L32_XP(t,pTwd, 3 * tw_step*sizeof(complex_fract16));
        tw10 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t));
        tw21 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(AE_L32_I(pTwd, 1*sizeof(complex_fract16))));
        tw31 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(AE_L32_I(pTwd, 2*sizeof(complex_fract16))));
        AE_L32_XP(t,pTwd, 3 * tw_step*sizeof(complex_fract16));
        tw11 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(t));
        // load 3x inputs
        x023=AE_L16X4_I (pX0,1*sizeof(ae_int16x4));
        x045=AE_L16X4_I (pX0,2*sizeof(ae_int16x4));
        AE_L16X4_IP(x001,pX0,3*sizeof(ae_int16x4));
        x123=AE_L16X4_I (pX1,1*sizeof(ae_int16x4));
        x145=AE_L16X4_I (pX1,2*sizeof(ae_int16x4));
        AE_L16X4_IP(x101,pX1,3*sizeof(ae_int16x4));
        x223=AE_L16X4_I (pX2,1*sizeof(ae_int16x4));
        x245=AE_L16X4_I (pX2,2*sizeof(ae_int16x4));
        AE_L16X4_IP(x201,pX2,3*sizeof(ae_int16x4));
        x323=AE_L16X4_I (pX3,1*sizeof(ae_int16x4));
        x345=AE_L16X4_I (pX3,2*sizeof(ae_int16x4));
        AE_L16X4_IP(x301,pX3,3*sizeof(ae_int16x4));
        // DFTs
        DFT4XI2(x001,x101,x201,x301);
        DFT4XI2(x023,x123,x223,x323);
        DFT4XI2(x045,x145,x245,x345);
        x101 = AE_MULFC16RAS(x101,tw10);
        x123 = AE_MULFC16RAS(x123,AE_SEL16_5432(tw10,tw11));
        x145 = AE_MULFC16RAS(x145,tw11);
        x201 = AE_MULFC16RAS(x201,tw20);
        x223 = AE_MULFC16RAS(x223,AE_SEL16_5432(tw20,tw21));
        x245 = AE_MULFC16RAS(x245,tw21);
        x301 = AE_MULFC16RAS(x301,tw30);
        x323 = AE_MULFC16RAS(x323,AE_SEL16_5432(tw30,tw31));
        x345 = AE_MULFC16RAS(x345,tw31);
        // save with permutation
        AE_S16X4X2RNG_IP(x001,AE_SEL16_7632(x023,x101),pY,sizeof(ae_int16x8));
        AE_S16X4X2RNG_IP(AE_SEL16_5432(x101,x123),x201,pY,sizeof(ae_int16x8));
        AE_S16X4X2RNG_IP(AE_SEL16_7632(x223,x301),AE_SEL16_5432(x301,x323),pY,sizeof(ae_int16x8));
        AE_S16X4X2RNG_IP(AE_SEL16_5432(x023,x045),AE_SEL16_5410(x045,x123),pY,sizeof(ae_int16x8));
        AE_S16X4X2RNG_IP(x145,AE_SEL16_5432(x223,x245),pY,sizeof(ae_int16x8));
        AE_S16X4X2RNG_IP(AE_SEL16_5410(x245,x323),x345,pY,sizeof(ae_int16x8));
    }
    // update scaling
    {int a, b; AE_CALCRNG16(a, b, 0, 3);*bexp = 3 - a;(void)b; }
    v[0] = 12;
    return shift;
}
