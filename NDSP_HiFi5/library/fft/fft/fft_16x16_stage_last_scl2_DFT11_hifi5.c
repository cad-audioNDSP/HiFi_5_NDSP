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
static uint32_t cmult16x16(uint32_t x, uint32_t y)
{
    int32_t re = S_round_l(L_sub_ll(L_mpy_ss(_LO(x), _LO(y)), L_mpy_ss(_HI(x), _HI(y))));
    int32_t im = S_round_l(L_add_ll(L_mpy_ss(_LO(x), _HI(y)), L_mpy_ss(_HI(x), _LO(y))));
    return _PACK(im, re);
}

static void __dft5x1(uint32_t *x0, uint32_t *x1, uint32_t *x2, uint32_t *x3, uint32_t *x4 )
{
    /*
    w1 =  exp(-1j*2*pi/5);
    w2 =  exp(-1j*2*pi*2/5);
    */
    const uint32_t  real_w1 = 0x0000278EUL;
    const uint32_t jimag_w1 = 0x86440000UL;
    const uint32_t  real_w2 = 0x00009872UL;
    const uint32_t jimag_w2 = 0xB4C30000UL;
    uint32_t y0, y1, y2, y3, y4;
    uint32_t s1, s2, d1, d2;

    /*
    DFT5 algorithm:
    x - input complex vector
    y - output complex vector
    y = fft(x)
    w1 =  exp(-1j*2*pi/5);
    w2 =  exp(-1j*2*pi*2/5);

    y = zeros(5,1);
    s1 = (x1+x4);
    s2 = (x2 + x3);
    d1 = (x1-x4);
    d2 = (x2-x3);

    y(1) = x0 + s1 + s2;
    y(2) = x0 + (s1*real(w1) + s2*real(w2)) + 1j*(d1*imag(w1) + d2*imag(w2));
    y(5) = x0 + (s1*real(w1) + s2*real(w2)) - 1j*(d1*imag(w1) + d2*imag(w2));
    y(3) = x0 + (s1*real(w2) + s2*real(w1)) + 1j*(d1*imag(w2)  - d2*imag(w1));
    y(4) = x0 + (s1*real(w2) + s2*real(w1)) - 1j*(d1*imag(w2)  - d2*imag(w1));
    */

    s1 = _ADD2S(*x1, *x4);
    s2 = _ADD2S(*x2, *x3);
    d1 = _SUB2S(*x1, *x4);
    d2 = _SUB2S(*x2, *x3);

    y0 = _ADD2S(*x0, _ADD2S(s1, s2));
    y1 = _ADD2S(*x0, _ADD2S(cmult16x16(s1, real_w1), cmult16x16(s2, real_w2)));
    y4 = y1;
    y2 = _ADD2S(*x0, _ADD2S(cmult16x16(s1, real_w2), cmult16x16(s2, real_w1)));
    y3 = y2;

    *x1 = _ADD2S(y1, _ADD2S(cmult16x16(d1, jimag_w1), cmult16x16(d2, jimag_w2)));
    *x4 = _SUB2S(y4, _ADD2S(cmult16x16(d1, jimag_w1), cmult16x16(d2, jimag_w2)));
    *x2 = _ADD2S(y2, _SUB2S(cmult16x16(d1, jimag_w2), cmult16x16(d2, jimag_w1)));
    *x3 = _SUB2S(y3, _SUB2S(cmult16x16(d1, jimag_w2), cmult16x16(d2, jimag_w1)));
    *x0 = y0;

}

static void __dft5x2(uint32_t *x0, uint32_t *x1, uint32_t *x2, uint32_t *x3, uint32_t *x4)
{
    __dft5x1(x0,   x1,   x2,   x3,   x4  ); 
    __dft5x1(x0+1, x1+1, x2+1, x3+1, x4+1);
}


static void __dft2x2(uint32_t *x0, uint32_t *x1)
{
    uint32_t y0, y1;
    y0 = _ADD2S(x0[0], x1[0]); 
    y1 = _SUB2S(x0[0], x1[0]); 
    x0[0] = y0; x1[0] = y1; 

    y0 = _ADD2S(x0[1], x1[1]);
    y1 = _SUB2S(x0[1], x1[1]);
    x0[1] = y0; x1[1] = y1;
}


#define DFT5X2(_x0, _x1, _x2, _x3, _x4) \
        __dft5x2((uint32_t*)&_x0, (uint32_t*)&_x1, (uint32_t*)&_x2, (uint32_t*)&_x3, (uint32_t*)&_x4)
#define DFT2xI2(_x0, _x1) \
    __dft2x2((uint32_t*)&_x0, (uint32_t*)&_x1)


static void IDFT10xI2(int16_t *y, int16_t *x)
#if 0
{
    uint64_t *px = (uint64_t *)x;
    uint64_t *py = (uint64_t *)y;
    uint64_t X[10];

    X[0] = px[0];
    X[1] = px[8];
    X[2] = px[6];
    X[3] = px[4];
    X[4] = px[2];
    X[5] = px[5];
    X[6] = px[3];
    X[7] = px[1];
    X[8] = px[9];
    X[9] = px[7];

    DFT5X2(X[0], X[1], X[2], X[3], X[4]);
    DFT5X2(X[5], X[6], X[7], X[8], X[9]);

    DFT2xI2(X[0], X[5]);
    DFT2xI2(X[1], X[6]);
    DFT2xI2(X[2], X[7]);
    DFT2xI2(X[3], X[8]);
    DFT2xI2(X[4], X[9]);

    py[0] = X[0];
    py[6] = X[1];
    py[2] = X[2];
    py[8] = X[3];
    py[4] = X[4];
    py[5] = X[5];
    py[1] = X[6];
    py[7] = X[7];
    py[3] = X[8];
    py[9] = X[9];
}
#else
{
    ae_int16x4 *px = (ae_int16x4 *)x;
    ae_int16x4 *py = (ae_int16x4 *)y;
    ae_int16x4 X[10];

    X[0] = px[0];
    X[1] = px[8];
    X[2] = px[6];
    X[3] = px[4];
    X[4] = px[2];
    X[5] = px[5];
    X[6] = px[3];
    X[7] = px[1];
    X[8] = px[9];
    X[9] = px[7];

    DFT5X2_(X[0], X[1], X[2], X[3], X[4], ((ae_int16x4*)__dft5_tw)[0], ((ae_int16x4*)__dft5_tw)[1],((ae_int16x4*)__dft5_tw)[2],((ae_int16x4*)__dft5_tw)[3]);
    DFT5X2_(X[5], X[6], X[7], X[8], X[9], ((ae_int16x4*)__dft5_tw)[0], ((ae_int16x4*)__dft5_tw)[1],((ae_int16x4*)__dft5_tw)[2],((ae_int16x4*)__dft5_tw)[3]);

    AE_ADDANDSUBRNG16RAS_S2(X[0], X[5]);
    AE_ADDANDSUBRNG16RAS_S2(X[1], X[6]);
    AE_ADDANDSUBRNG16RAS_S2(X[2], X[7]);
    AE_ADDANDSUBRNG16RAS_S2(X[3], X[8]);
    AE_ADDANDSUBRNG16RAS_S2(X[4], X[9]);

    py[0] = X[0];
    py[6] = X[1];
    py[2] = X[2];
    py[8] = X[3];
    py[4] = X[4];
    py[5] = X[5];
    py[1] = X[6];
    py[7] = X[7];
    py[3] = X[8];
    py[9] = X[9];
}
#endif



static void DFT10xI2(int16_t *y, int16_t *x)
#if 0
{
    ae_int16x4 *px = (ae_int16x4 *)x;
    ae_int16x4 *py = (ae_int16x4 *)y;
    ae_int16x4 X[10];

    X[0] = px[0];
    X[1] = px[2];
    X[2] = px[4];
    X[3] = px[6];
    X[4] = px[8];
    X[5] = px[5];
    X[6] = px[7];
    X[7] = px[9];
    X[8] = px[1];
    X[9] = px[3];

    DFT5X2(X[0], X[1], X[2], X[3], X[4]);
    DFT5X2(X[5], X[6], X[7], X[8], X[9]);

    DFT2xI2(X[0], X[5]);
    DFT2xI2(X[1], X[6]);
    DFT2xI2(X[2], X[7]);
    DFT2xI2(X[3], X[8]);
    DFT2xI2(X[4], X[9]);

    py[0] = X[0];
    py[6] = X[1];
    py[2] = X[2];
    py[8] = X[3];
    py[4] = X[4];
    py[5] = X[5];
    py[1] = X[6];
    py[7] = X[7];
    py[3] = X[8];
    py[9] = X[9];
}
#else
{

    ae_int16x4 *px = (ae_int16x4 *)x;
    ae_int16x4 *py = (ae_int16x4 *)y;
    ae_int16x4 X[10];

    X[0] = px[0];
    X[1] = px[2];
    X[2] = px[4];
    X[3] = px[6];
    X[4] = px[8];
    X[5] = px[5];
    X[6] = px[7];
    X[7] = px[9];
    X[8] = px[1];
    X[9] = px[3];

    DFT5X2_(X[0], X[1], X[2], X[3], X[4], ((ae_int16x4*)__dft5_tw)[0], ((ae_int16x4*)__dft5_tw)[1],((ae_int16x4*)__dft5_tw)[2],((ae_int16x4*)__dft5_tw)[3]);
    DFT5X2_(X[5], X[6], X[7], X[8], X[9], ((ae_int16x4*)__dft5_tw)[0], ((ae_int16x4*)__dft5_tw)[1],((ae_int16x4*)__dft5_tw)[2],((ae_int16x4*)__dft5_tw)[3]);
    AE_ADDANDSUBRNG16RAS_S2(X[0], X[5]);
    AE_ADDANDSUBRNG16RAS_S2(X[1], X[6]);
    AE_ADDANDSUBRNG16RAS_S2(X[2], X[7]);
    AE_ADDANDSUBRNG16RAS_S2(X[3], X[8]);
    AE_ADDANDSUBRNG16RAS_S2(X[4], X[9]);

    py[0] = X[0];
    py[6] = X[1];
    py[2] = X[2];
    py[8] = X[3];
    py[4] = X[4];
    py[5] = X[5];
    py[1] = X[6];
    py[7] = X[7];
    py[3] = X[8];
    py[9] = X[9];
}
#endif

ALIGN(32) static const int32_t A_Q16_15[20] =
{

    -3277, 0,
    3130, -10407,
    8638, 6595,
    8327, 6983,
    6784, 8491,
    0, -10868,
    6784, -8491,
    -8327, 6983,
    8638, -6595,
    -3130, -10407,
};

static void _mpy_Q16_15(uint64_t *z, const uint64_t *x, const uint64_t *y)
{
    int16_t y16[4]; 
    int32_t *x32 = (int32_t*)x;
    int16_t *z16 = (int16_t*)z; 

    y16[0] = ((int16_t*)y)[0];
    y16[1] = ((int16_t*)y)[1];
    y16[2] = ((int16_t*)y)[2];
    y16[3] = ((int16_t*)y)[3];
    z16[0] = (int16_t)L_sub_ll(L_mpy_ls(x32[0], y16[0]), L_mpy_ls(x32[1], y16[1])); 
    z16[1] = (int16_t)L_add_ll(L_mpy_ls(x32[0], y16[1]), L_mpy_ls(x32[1], y16[0]));
    z16[2] = (int16_t)L_sub_ll(L_mpy_ls(x32[0], y16[2]), L_mpy_ls(x32[1], y16[3]));
    z16[3] = (int16_t)L_add_ll(L_mpy_ls(x32[0], y16[3]), L_mpy_ls(x32[1], y16[2]));
}

static uint64_t _addx4(uint64_t x, uint64_t y)
{
    uint64_t z;
    int16_t *z16 = (int16_t*)&z;
    int16_t *x16 = (int16_t*)&x;
    int16_t *y16 = (int16_t*)&y;

    z16[0] = S_add_ss(x16[0], y16[0]);
    z16[1] = S_add_ss(x16[1], y16[1]);
    z16[2] = S_add_ss(x16[2], y16[2]);
    z16[3] = S_add_ss(x16[3], y16[3]);

    return z;
}

static  uint64_t _shrx4(uint64_t x, int shift)
{
    uint64_t z;
    int16_t *z16 = (int16_t*)&z;
    int16_t *x16 = (int16_t*)&x;
    int16_t rnd = (1 << shift) >> 1;

    ASSERT(shift < 15); 

    z16[0] = S_shr_s(S_add_ss(x16[0], rnd), (int16_t)shift);
    z16[1] = S_shr_s(S_add_ss(x16[1], rnd), (int16_t)shift);
    z16[2] = S_shr_s(S_add_ss(x16[2], rnd), (int16_t)shift);
    z16[3] = S_shr_s(S_add_ss(x16[3], rnd), (int16_t)shift);

    return z;
}
#endif

// twiddles for DFT5X2_
ALIGN(32) static const int16_t __dft5_tw[] =
{
    (int16_t)0x278E, (int16_t)0x278E, (int16_t)0x278E, (int16_t)0x278E,
    (int16_t)0x8644, (int16_t)0x79BC, (int16_t)0x8644, (int16_t)0x79BC,
    (int16_t)0x9872, (int16_t)0x9872, (int16_t)0x9872, (int16_t)0x9872,
    (int16_t)0xB4C3, (int16_t)0x4B3D, (int16_t)0xB4C3, (int16_t)0x4B3D
};

#define DFT5X2_(x0, x1, x2, x3, x4, w1, w2, w3, w4)\
{                                          \
    ae_int16x4 s1, s2, d1, d2;             \
    ae_int16x4 t0,  t2    ;                \
    ae_int16x4 y0, y1, y2, y3;             \
    AE_ADDANDSUBRNG16RAS_S2(x1,x4);        \
    AE_ADDANDSUBRNG16RAS_S2(x2,x3);        \
    s1=x1;d1=x4;s2=x2;d2=x3;               \
                                           \
    t0 = AE_MULFD16X16X4RAS(s1, s2,w1,w3); \
    t2 = AE_MULFD16X16X4RAS(s1, s2,w3,w1); \
    y0 = AE_ADD16S(x0, t0);                \
    y1 = AE_ADD16S(x0, t2);                \
                                           \
    y2 = AE_MULFD16X16X4RAS(d1,d2, w2,w4); \
    y3 = AE_MULFD16X16X4RAS(d1,d2, w4,AE_NEG16S(w2));          \
    y2 = AE_SEL16_2301(y2, y2);            \
    y3 = AE_SEL16_2301(y3, y3);            \
                                           \
    x0 = AE_ADD16S(x0, AE_ADD16S(s1, s2)); \
    AE_ADDANDSUBRNG16RAS_S2(y0, y2);       \
    AE_ADDANDSUBRNG16RAS_S2(y1, y3);       \
    x1=y0;x4=y2;x2=y1;x3=y3;               \
}


static void DFT11x2_phase (uint64_t *y, uint64_t *x, int stride)
{
    // reordered original twiddles for DFT10->DFT11
    ALIGN(32) static const int16_t A_[40] =
    {
        -3130, -10407,  -3130, -10407,   // 9
        -3277, 0,       -3277, 0,     // 0
        -8327, 6983,    -8327, 6983,     // 7
        3130, -10407,   3130, -10407,    // 1
        8327, 6983,     8327, 6983,      // 3
        8638, -6595,    8638, -6595,     // 8
        6784, -8491,    6784, -8491,     // 6
        8638, 6595,     8638, 6595,      // 2
        0, -10868,      0, -10868,       // 5
        6784, 8491,     6784, 8491       // 4
    };
    const ae_int16x4* restrict pX =(const ae_int16x4*)x;
    ae_int16x4 *restrict pY =(       ae_int16x4 *)y;
    const ae_int16x8* restrict pA = (const ae_int16x8 *)A_;
    const ae_int16x8 *pdft5twd=(const ae_int16x8 *)__dft5_tw;
    int j;
    ae_int16x4 w0,w1,w2,w3;
    for (j=0; j<stride; j++)
    {
        ae_int16x4 X[11], C[10];
        ae_int16x4 y0,a0,a1;
        ae_int16x4 X0,Y0;

        X[10]=AE_L16X4_X(pX,stride*sizeof(ae_int16x4));        AE_L16X4_XP(X[0],pX,2*stride*sizeof(ae_int16x4));
        X[8] =AE_L16X4_X(pX,stride*sizeof(ae_int16x4));        AE_L16X4_XP(X[1],pX,2*stride*sizeof(ae_int16x4));
        X[4] =AE_L16X4_X(pX,stride*sizeof(ae_int16x4));        AE_L16X4_XP(X[2],pX,2*stride*sizeof(ae_int16x4));
        X[7] =AE_L16X4_X(pX,stride*sizeof(ae_int16x4));        AE_L16X4_XP(X[9],pX,2*stride*sizeof(ae_int16x4));
        X[6] =AE_L16X4_X(pX,stride*sizeof(ae_int16x4));        AE_L16X4_XP(X[3],pX,2*stride*sizeof(ae_int16x4));
        AE_L16X4_XP(X[5],pX,(1-10*stride)*sizeof(ae_int16x4));
        y0 = X[0];
        X0=y0;
        y0 = AE_ADD16S(X[0],X[1]);
        y0 = AE_ADD16S(y0,AE_ADD16S(AE_ADD16S(X[2],X[9]),AE_ADD16S(X[4],X[7])));
        y0 = AE_ADD16S(y0,AE_ADD16S(AE_ADD16S(X[6],X[5]),AE_ADD16S(X[8],X[10])));
        y0 = AE_ADD16S(y0,X[3]);

        {
            ae_int16x4 Y[10];

            Y[0] = X[0+1];Y[1] = X[2+1]; Y[2] = X[4+1];Y[3] = X[6+1];
            Y[4] = X[8+1];Y[5] = X[5+1]; Y[6] = X[7+1];Y[7] = X[9+1];
            Y[8] = X[1+1];Y[9] = X[3+1];
            AE_L16X4X2_IP(w0,w1,pdft5twd,sizeof(ae_int16x8));
            AE_L16X4X2_XP(w2,w3,pdft5twd,-(int)sizeof(ae_int16x8));
            DFT5X2_(Y[0], Y[1], Y[2], Y[3], Y[4], w0,w1,w2,w3);
            DFT5X2_(Y[5], Y[6], Y[7], Y[8], Y[9], w0,w1,w2,w3);
            AE_ADDANDSUBRNG16RAS_S2(Y[0], Y[5]);
            AE_ADDANDSUBRNG16RAS_S2(Y[1], Y[6]);
            AE_ADDANDSUBRNG16RAS_S2(Y[2], Y[7]);
            AE_ADDANDSUBRNG16RAS_S2(Y[3], Y[8]);
            AE_ADDANDSUBRNG16RAS_S2(Y[4], Y[9]);

            C[0] = Y[0];    C[6] = Y[1];    C[2] = Y[2];    C[8] = Y[3];
            C[4] = Y[4];    C[5] = Y[5];    C[1] = Y[6];    C[7] = Y[7];
            C[3] = Y[8];    C[9] = Y[9];
        }
        pA = (const ae_int16x8 *)A_;
                   
        Y0=y0;

        AE_L16X4X2_IP(a0,a1,pA,sizeof(ae_int16x8));

        C[9] = AE_MULFC16RAS(C[9],a0);
        C[0] = AE_MULFC16RAS(C[0],a1);

        AE_L16X4X2_IP(a0,a1,pA,sizeof(ae_int16x8));

        C[7] = AE_MULFC16RAS(C[7],a0);
        C[1] = AE_MULFC16RAS(C[1],a1);

        AE_L16X4X2_IP(a0,a1,pA,sizeof(ae_int16x8));

        C[3] = AE_MULFC16RAS(C[3],a0);
        C[8] = AE_MULFC16RAS(C[8],a1);

        AE_L16X4X2_IP(a0,a1,pA,sizeof(ae_int16x8));

        C[6] = AE_MULFC16RAS(C[6],a0);
        C[2] = AE_MULFC16RAS(C[2],a1);

        AE_L16X4X2_IP(a0,a1,pA,sizeof(ae_int16x8));

        C[5] = AE_MULFC16RAS(C[5],a0);
        C[4] = AE_MULFC16RAS(C[4],a1);

        {
            ae_int16x4 X[10];
            

            X[0] = C[0];X[1] = C[8];X[2] = C[6];X[3] = C[4];
            X[4] = C[2];X[5] = C[5];X[6] = C[3];X[7] = C[1];
            X[8] = C[9];X[9] = C[7];


            DFT5X2_(X[0], X[1], X[2], X[3], X[4], w0,w1,w2,w3);
            DFT5X2_(X[5], X[6], X[7], X[8], X[9], w0,w1,w2,w3);

            AE_ADDANDSUBRNG16RAS_S2(X[0], X[5]);
            AE_ADDANDSUBRNG16RAS_S2(X[1], X[6]);
            AE_ADDANDSUBRNG16RAS_S2(X[2], X[7]);
            AE_ADDANDSUBRNG16RAS_S2(X[3], X[8]);
            AE_ADDANDSUBRNG16RAS_S2(X[4], X[9]);

            C[0] = X[0];C[6] = X[1];C[2] = X[2];C[8] = X[3];
            C[4] = X[4];C[5] = X[5];C[1] = X[6];C[7] = X[7];
            C[3] = X[8];C[9] = X[9];
        }
        //        IDFT10xI2((int16_t*)C, (int16_t*)C);

        AE_S16X4_XP(Y0,pY,stride*sizeof(ae_int16x4));
        AE_S16X4_XP(AE_ADD16S(X0,C[9]),pY,stride*sizeof(ae_int16x4));
        AE_S16X4_XP(AE_ADD16S(X0,C[8]),pY,stride*sizeof(ae_int16x4));
        AE_S16X4_XP(AE_ADD16S(X0,C[1]),pY,stride*sizeof(ae_int16x4));
        AE_S16X4_XP(AE_ADD16S(X0,C[7]),pY,stride*sizeof(ae_int16x4));
        AE_S16X4_XP(AE_ADD16S(X0,C[5]),pY,stride*sizeof(ae_int16x4));
        AE_S16X4_XP(AE_ADD16S(X0,C[0]),pY,stride*sizeof(ae_int16x4));
        AE_S16X4_XP(AE_ADD16S(X0,C[2]),pY,stride*sizeof(ae_int16x4));
        AE_S16X4_XP(AE_ADD16S(X0,C[6]),pY,stride*sizeof(ae_int16x4));
        AE_S16X4_XP(AE_ADD16S(X0,C[3]),pY,stride*sizeof(ae_int16x4));
        AE_S16X4_XP(AE_ADD16S(X0,C[4]),pY,(1-10*stride)*sizeof(ae_int16x4));

    }
}




inline_ void DFT11x2_scale (uint64_t* x, int N, int shift)
#if 0
{
    int j;
    NASSERT(shift>=0);
    NASSERT(N%2==0);
    NASSERT_ALIGN16(x);
    for (j=0; j<N; j++)
    {
        x[j]= _shrx4(x[j], shift);
    }
}
#else
{
    const ae_int16x8* restrict pX=(const ae_int16x8*)x;
          ae_int16x8* restrict pY=(      ae_int16x8*)x;
    ae_int16x4 scale=AE_SLAA16S(AE_MOVDA16(0x1),15-shift);
    int j;
    NASSERT(shift>=0);
    NASSERT(N%2==0);
    NASSERT_ALIGN16(x);
    for (j=0; j<(N>>1); j++)
    {
        ae_int16x4 x0,x1;
        AE_L16X4X2_IP(x0,x1,pX,sizeof(ae_int16x8));
        x0=AE_MULFP16X4RS(x0,scale);
        x1=AE_MULFP16X4RS(x1,scale);
        AE_S16X4X2_IP(x0,x1,pY,sizeof(ae_int16x8));
    }
}
#endif

/*
*  Last stage of FFT/IFFT 16x16, radix-11, dynamic scaling
    NOTE: last stage might be inplace!
*/
int fft_16x16_stage_last_scl2_DFT11(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
{
    const int N11=AE_MOVAD32_L(AE_MULFP32X2RAS(AE_MOVDA32(N),195225786));//N/11
    const int stride =N11;
    int shift;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(N%4==0);
    WUR_AE_SAR(0);  // required for subsequent addandsub operations!

    shift = XT_MAX(0,4 - *bexp);
    DFT11x2_scale((uint64_t*)x, (N>>1),shift);
    DFT11x2_phase((uint64_t*)y, (uint64_t*)x, (stride>>1));

    return shift;

} /* fft_16x16_stage_last_scl2_DFT11 */
