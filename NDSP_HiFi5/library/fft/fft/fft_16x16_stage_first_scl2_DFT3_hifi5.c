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
    Optimized for HiFi5
	Integrit, 2006-2019
*/
#include "NatureDSP_Signal_fft.h"
#include "common.h"
/* Twiddle factor tables and FFT descriptor structure. */
#include "fft_x16_common.h"
#include "fft_16x16_stages.h"

#if 0
#include "baseop.h"
/* FFT cplx16x16 stages */
inline_ static uint16_t _HI(uint32_t x) { return (uint16_t)(x >> 16); }
inline_ static uint16_t _LO(uint32_t x) { return (uint16_t)(x); }
inline_ static uint32_t _PACK(uint16_t hi, uint16_t lo) { return (((uint32_t)hi) << 16) | (lo); }
inline_ static uint32_t _SHR2(uint32_t x, int rsh)
{
  int16_t hi = _HI(x), lo = _LO(x);
  hi >>= rsh; lo >>= rsh;
  return _PACK(hi, lo);
}
inline_ static uint32_t _SHR2R(uint32_t x, int rsh)
{
  int32_t hi = (int32_t)(int16_t)_HI(x) + (1L<<(rsh-1));
  int32_t lo = (int32_t)(int16_t)_LO(x) + (1L<<(rsh-1));
  hi >>= rsh; lo >>= rsh;
  return _PACK(S_sature_l(hi), S_sature_l(lo));
}
inline_ static uint32_t _ADD2S(uint32_t x, uint32_t y)
{
  return _PACK(S_add_ss(_HI(x), _HI(y)), S_add_ss(_LO(x), _LO(y)));
}
inline_ static uint32_t _SUB2S(uint32_t x, uint32_t y)
{
  return _PACK(S_sub_ss(_HI(x), _HI(y)), S_sub_ss(_LO(x), _LO(y)));
}
/* Two-way bidirectional arithmetic (sign-extending) right shift */
inline_ static uint32_t _SHR2R_BIDIR(uint32_t x, int rsh)
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

inline_ static uint32_t cmult16x16(uint32_t x, uint32_t y)
{
    int32_t re = S_round_l(L_sub_ll(L_mpy_ss(_LO(x), _LO(y)), L_mpy_ss(_HI(x), _HI(y))));
    int32_t im = S_round_l(L_add_ll(L_mpy_ss(_LO(x), _HI(y)), L_mpy_ss(_HI(x), _LO(y))));
    return _PACK(im, re);
}

inline_ static void _cmult16x16(uint32_t *result, uint32_t *x, uint32_t *y)
{
    int16_t re = S_extract_l(L_sub_ll(L_mpy_ss(_LO(*x), _LO(*y)), L_mpy_ss(_HI(*x), _HI(*y))));
    int16_t im = S_extract_l(L_add_ll(L_mpy_ss(_LO(*x), _HI(*y)), L_mpy_ss(_HI(*x), _LO(*y))));
    *result = _PACK(im, re);
}
#endif

int fft_16x16_stage_first_scl2_DFT3(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v_, int tw_step, int *bexp)
#if 0
{
    int i;
    const int stride = N / 3;
    const int v = *v_;
    uint32_t * restrict px = (uint32_t *)x;
    uint32_t * restrict py = (uint32_t *)y;
    uint32_t * restrict ptw = (uint32_t *)tw;
    uint32_t x0, x1, x2, s, d, y0, y1, y2, c;
    uint32_t t1, t2;
    const int shift = XT_MAX(0, 3 - *bexp);

    NASSERT(shift>-32 && shift<32); 
    NASSERT(N%12==0);
    // only in  the first stage!!
    NASSERT(v != stride);
    NASSERT(v == 1);

    for (i = 0; i < N/3; i++)
    {
        t1 = ptw[0];
        t2 = ptw[1];
        ptw += 2*tw_step;
        x0 = px[i + 0];
        x1 = px[i + 1 * stride];
        x2 = px[i + 2 * stride];

        x0 = _SHR2R_BIDIR(x0, shift);
        x1 = _SHR2R_BIDIR(x1, shift);
        x2 = _SHR2R_BIDIR(x2, shift);
        /*
        DFT3 algorithm:
        x - input complex vector
        y - output complex vector
        y = fft(x)
        y = [ x(1) + x(2)  + x(3);
        x(1) + (x(2) + x(3))*cos(2*pi/3) - 1j*(x(2) - x(3))*sin(2*pi/3);
        x(1) + (x(2) + x(3))*cos(2*pi/3) + 1j*(x(2) - x(3))*sin(2*pi/3) ]
        */

        s = _ADD2S(x1, x2);
        y0 = _ADD2S(x0, s);

        s = _SHR2(s, 1); // cos(2*pi/3) = -0.5
        d = _SUB2S(x1, x2);
        c = 0x6EDA0000;  // sin(2*pi/3)*j, Q15 format            
        d = cmult16x16(d, c);

        s = _SUB2S(x0, s);
        y2 = _ADD2S(s, d);
        y1 = _SUB2S(s, d);
        _cmult16x16(&y1, &y1, &t1);
        _cmult16x16(&y2, &y2, &t2);
        py[3*i + 0] = y0;
        py[3*i + 1 * v] = y1;
        py[3*i + 2 * v] = y2;
    }
    *v_ = 3;
    *bexp = vec_bexp16(y, 2*N) - 16;
    return shift;
}
#else
{
    /*  AE_SEL16_7632();
        AE_SEL16_5410();    */
    ALIGN(32) static const int16_t sel_tab[4] = {0x705,0x604,0x301,0x200};
    ae_int16x4 dsel=AE_L16X4_I((const ae_int16x4*)sel_tab,0);

    int i;
    const int N3=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),715827883));//N/3

    const int stride = N3;
    const int v = *v_;
    const ae_int16x4 * restrict px0 = (const ae_int16x4 *)x;
    const ae_int16x4 * restrict px1 = (const ae_int16x4 *)(x+1*2 * stride);
    const ae_int16x4 * restrict px2 = (const ae_int16x4 *)(x+2*2 * stride);
    ae_int16x4 * restrict py = (ae_int16x4 *)y;
    ae_int16x4 * restrict ptw = (ae_int16x4 *)tw;
    const int shift = XT_MAX(0, 3 - *bexp);
    ae_int16x4 scale_shift;

    NASSERT(shift>-32 && shift<32); 
    NASSERT(N%12==0);
    // only in  the first stage!!
    NASSERT(v != stride);
    NASSERT(v == 1); (void)v;
    /* Reset RANGE register */
    WUR_AE_SAR( 0);//shift * 0x102); 
    scale_shift=AE_SLAA16S(AE_MOVDA16(0x1000),3-shift);

    for (i = 0; i < N3; i+=4)
    {
        ae_int16x4 X0[2],X1[2],X2[2],S[2],D[2],Y0[2],Y1[2],Y2[2],T1[2],T2[2];
        AE_L16X4_XP(X0[0],ptw,tw_step*sizeof(ae_int16x4));
        AE_L16X4_XP(X1[0],ptw,tw_step*sizeof(ae_int16x4));
        AE_L16X4_XP(X0[1],ptw,tw_step*sizeof(ae_int16x4));
        AE_L16X4_XP(X1[1],ptw,tw_step*sizeof(ae_int16x4));
        AE_DSEL16X4(T1[0],T2[0],X0[0],X1[0],dsel);
        AE_DSEL16X4(T1[1],T2[1],X0[1],X1[1],dsel);

        AE_L16X4X2_IP(X0[0],X0[1],castxcc(ae_int16x8,px0),sizeof(ae_int16x8));
        AE_L16X4X2_IP(X1[0],X1[1],castxcc(ae_int16x8,px1),sizeof(ae_int16x8));
        AE_L16X4X2_IP(X2[0],X2[1],castxcc(ae_int16x8,px2),sizeof(ae_int16x8));

        /*
        DFT3 algorithm:
        x - input complex vector
        y - output complex vector
        y = fft(x)
        y = [ x(1) + x(2)  + x(3);
        x(1) + (x(2) + x(3))*cos(2*pi/3) - 1j*(x(2) - x(3))*sin(2*pi/3);
        x(1) + (x(2) + x(3))*cos(2*pi/3) + 1j*(x(2) - x(3))*sin(2*pi/3) ]
        */
        X0[0]=AE_SRAA16RS(X0[0],shift);
        X0[1]=AE_SRAA16RS(X0[1],shift);
        S[0]=AE_MULFP16X4RS(X1[0],scale_shift);
        S[1]=AE_MULFP16X4RS(X1[1],scale_shift);
        D[0]=AE_MULFP16X4RS(X2[0],scale_shift);
        D[1]=AE_MULFP16X4RS(X2[1],scale_shift);
        AE_ADDANDSUBRNG16RAS_S2(S[0],D[0]);
        AE_ADDANDSUBRNG16RAS_S2(S[1],D[1]);
        Y0[0]=AE_ADD16S(X0[0],S[0]);
        Y0[1]=AE_ADD16S(X0[1],S[1]);
        S[0]=AE_SRAI16(S[0],1);
        S[1]=AE_SRAI16(S[1],1);
        S[0]=AE_SUB16S(X0[0],S[0]);
        S[1]=AE_SUB16S(X0[1],S[1]);
        D[0]=AE_MULFC16RAS(D[0],AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32(0x00006EDA)));
        D[1]=AE_MULFC16RAS(D[1],AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32(0x00006EDA)));
        AE_ADDANDSUBRNG16RAS_S2(S[0],D[0]); 
        AE_ADDANDSUBRNG16RAS_S2(S[1],D[1]); 
        Y1[0]=AE_MULFC16RAS(D[0],T1[0]);
        Y1[1]=AE_MULFC16RAS(D[1],T1[1]);
        Y2[0]=AE_MULFC16RAS(S[0],T2[0]);
        Y2[1]=AE_MULFC16RAS(S[1],T2[1]);
        AE_S16X4X2RNG_IP(AE_SEL16_7632(Y0[0],Y1[0]),AE_SEL16_7610(Y2[0],Y0[0]),castxcc(ae_int16x8,py),sizeof(ae_int16x8));
        AE_S16X4X2RNG_IP(AE_SEL16_5410(Y1[0],Y2[0]),AE_SEL16_7632(Y0[1],Y1[1]),castxcc(ae_int16x8,py),sizeof(ae_int16x8));
        AE_S16X4X2RNG_IP(AE_SEL16_7610(Y2[1],Y0[1]),AE_SEL16_5410(Y1[1],Y2[1]),castxcc(ae_int16x8,py),sizeof(ae_int16x8));
    }
    *v_ = 3;
    // update range
    {    int a, b;    AE_CALCRNG16(a, b, 0, 3);    *bexp = 3 - a;    (void)b;   }
    return shift;
}
#endif
