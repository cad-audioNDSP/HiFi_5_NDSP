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
static uint16_t _HI(uint32_t x) { return (uint16_t)(x >> 16); }
static uint16_t _LO(uint32_t x) { return (uint16_t)(x); }
static uint32_t _PACK(uint16_t hi, uint16_t lo) { return (((uint32_t)hi) << 16) | (lo); }
static uint32_t _SHR2(uint32_t x, int rsh)
{
  int16_t hi = _HI(x), lo = _LO(x);
  hi >>= rsh; lo >>= rsh;
  return _PACK(hi, lo);
}
static uint32_t _SHR2R(uint32_t x, int rsh)
{
  int32_t hi = (int32_t)(int16_t)_HI(x) + (1L<<(rsh-1));
  int32_t lo = (int32_t)(int16_t)_LO(x) + (1L<<(rsh-1));
  hi >>= rsh; lo >>= rsh;
  return _PACK(S_sature_l(hi), S_sature_l(lo));
}
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

static uint32_t cmult16x16(uint32_t x, uint32_t y)
{
    int32_t re = S_round_l(L_sub_ll(L_mpy_ss(_LO(x), _LO(y)), L_mpy_ss(_HI(x), _HI(y))));
    int32_t im = S_round_l(L_add_ll(L_mpy_ss(_LO(x), _HI(y)), L_mpy_ss(_HI(x), _LO(y))));
    return _PACK(im, re);
}

static void _cmult16x16(uint32_t *result, uint32_t *x, uint32_t *y)
{
    int16_t re = S_extract_l(L_sub_ll(L_mpy_ss(_LO(*x), _LO(*y)), L_mpy_ss(_HI(*x), _HI(*y))));
    int16_t im = S_extract_l(L_add_ll(L_mpy_ss(_LO(*x), _HI(*y)), L_mpy_ss(_HI(*x), _LO(*y))));
    *result = _PACK(im, re);
}
#endif

/*
 *  Intermediate stage of FFT/IFFT 16x16, radix-2, dynamic scaling
 */
int fft_16x16_stage_inner_scl2_DFT2_v3(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
#if 0
{
    int i, j;
    const int stride = N / 2;

    uint32_t * restrict px = (uint32_t *)x;
    uint32_t * restrict py = (uint32_t *)y;
    uint32_t * restrict ptw = (uint32_t *)tw;
    int shift;
    uint32_t x0, x1, y0, y1;
    uint32_t t1;
    int min_shift;
    NASSERT(v[0]==3);
    NASSERT(stride != v[0]); // inner - not last
    NASSERT(N%24==0);
    min_shift = 2;

    shift = min_shift - *bexp; 
    ASSERT(shift>-32 && shift<32);

    for (i = 0; i < N/6; i++)
    {
        t1 = *ptw; ptw += tw_step;
        for (j = 0; j < 3; j++)
        {
            x0 = px[j + 3 * i + 0];
            x1 = px[j + 3 * i + 1 * stride];

            x0 = _SHR2R_BIDIR(x0, shift);
            x1 = _SHR2R_BIDIR(x1, shift);

            y0 = _ADD2S(x0, x1);
            y1 = _SUB2S(x0, x1);

            _cmult16x16(&y1, &y1, &t1);

            py[j + 6 * i + 0] = y0;
            py[j + 6 * i + 3] = y1;
        }
    }

    *bexp = vec_bexp16(y, 2*N) - 16;
    *v = 6;
    return shift;
}
#else
{
    int i;
    const int N12=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),178956971));//N/12
    const int stride = N / 2;
    const int shift = XT_MAX(0, 2 - *bexp);

    ae_int16x4 * restrict px0= (ae_int16x4 *)x;
    ae_int16x4 * restrict px1= (ae_int16x4 *)(x+stride*2);
    ae_int16x8* restrict py = (ae_int16x8*)y;
    ae_int32 * restrict ptw = (ae_int32*)tw;
    NASSERT(v[0]==3);
    NASSERT(stride != v[0]); // inner - not last
    NASSERT(N%24==0);
    NASSERT(shift>-32 && shift<32);
    /* Reset RANGE register */
    {
        int a, b; AE_CALCRNG16(a, b, 0, 3); (void)b; (void)a;
    }
    WUR_AE_SAR( shift*0x102 ); 

    for (i = 0; i < N12; i++)
    {
        ae_int32x2 tmp;
        ae_int16x4 tw0,tw1,tw2;
        ae_int16x4 x001,x023,x045,x101,x123,x145;
        ae_int16x4 y001,y023,y045,y101,y123,y145;
        AE_L16X4_IP(x001,px0,1*sizeof(ae_int16x4));
        AE_L16X4_IP(x023,px0,1*sizeof(ae_int16x4));
        AE_L16X4_IP(x045,px0,1*sizeof(ae_int16x4));
        AE_L16X4_IP(x101,px1,1*sizeof(ae_int16x4));
        AE_L16X4_IP(x123,px1,1*sizeof(ae_int16x4));
        AE_L16X4_IP(x145,px1,1*sizeof(ae_int16x4));

        AE_ADDANDSUBRNG16RAS_S2(x001, x101); 
        AE_ADDANDSUBRNG16RAS_S2(x023, x123); 
        AE_ADDANDSUBRNG16RAS_S2(x045, x145); 

        AE_L32_XP(tmp,ptw,tw_step*sizeof(ae_int32));
        tw0 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(tmp));
        AE_L32_XP(tmp,ptw,tw_step*sizeof(ae_int32));
        tw1 = AE_SHORTSWAP(AE_MOVINT16X4_FROMINT32X2(tmp));
        y001 =x001 ;
        y023 =x023 ;
        y045 =x045 ;
        tw2=AE_SEL16_5410(tw0,tw1);
        y101 = AE_MULFC16RAS(x101, tw0 );
        y123 = AE_MULFC16RAS(x123, tw2 );
        y145 = AE_MULFC16RAS(x145, tw1 );

        AE_S16X4X2RNG_IP(y001, AE_SEL16_7632(y023,y101)                    , py, 2*sizeof(ae_int16x4));
        AE_S16X4X2RNG_IP(AE_SEL16_5432(y101,y123), AE_SEL16_5432(y023,y045), py, 2*sizeof(ae_int16x4));
        AE_S16X4X2RNG_IP(AE_SEL16_5410(y045,y123), y145                    , py, 2*sizeof(ae_int16x4));
    }
    // update scale factor
    {int a, b; AE_CALCRNG16(a, b, 0, 3); *bexp = 3 - a;(void)b;  }
    *v = 6;
    return shift;
}
#endif
