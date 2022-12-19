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
    Discrete Cosine Transform, Type II
    C code optimized for HiFi4
   Integrit, 2006-2019
*/

/* Signal Processing Library API. */
#include "NatureDSP_Signal_fft.h"
#include "common.h"
#include "dct2_twd.h"

#define FIRST_STAGE_SCALE 3
#define LAST_STAGE_SCALE 2

/* DSEL16X4 patterns table */
static const int16_t ALIGN(32) dsel_tbl_perm[] = {
    (0<<8)|7, (2<<8)|5, (4<<8)|3, (6<<8)|1,/* Table used in permutations before FFT */
    (5<<8)|7, (4<<8)|6, (1<<8)|3, (0<<8)|2 /* Table used in samples permutation in the 1st stage of FFT */
};

/* kron(exp( -1j*2*pi/8*(1:3)'), [1;1] ) */
static const int16_t ALIGN(32) __fft8_tw1_v2_[] = 
{
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x5A82, (int16_t)0xA57E,
    (int16_t)0x0000, (int16_t)0x8000, (int16_t)0x0000, (int16_t)0x8000,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0xA57E,
};

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

/*
   scaled fft with reordering
   N=16 - size of FFT
   NOTE: y is input and output
*/
void fft16_16x16(int16_t *y, const int16_t *ptwd)
{
    //const int N = 16;
    const ae_int16x8 * restrict p_twd;
    const ae_int16x8 * restrict px0;
          ae_int16x8 * restrict py0;
    ae_int16x4 x0, x1, x2, x3, x4, x5, x6, x7;
    ae_int16x4 z0, z1, z2, z3, z4, z5, z6, z7;
    ae_int16x4 dselTbl;
    ae_int16x4 tw0102, tw0311, tw1213;
    ae_int16x4 tw2122, tw2331, tw3233;

    NASSERT_ALIGN(y, 16);

    p_twd = (const ae_int16x8 *)ptwd;
    px0 = (const ae_int16x8 *)y;
    py0 = (      ae_int16x8 *)y;

    /* Load data */
    AE_L16X4X2_IP(z0, z1, px0, 2*sizeof(ae_int16x4));
    AE_L16X4X2_IP(z2, z3, px0, 2*sizeof(ae_int16x4));
    AE_L16X4X2_IP(z4, z5, px0, 2*sizeof(ae_int16x4));
    AE_L16X4X2_IP(z6, z7, px0, 2*sizeof(ae_int16x4));
    /* Reordering */
    dselTbl = AE_L16X4_I((ae_int16x4 *)dsel_tbl_perm, 0);
    AE_DSEL16X4(x7, x0, z0, z1, dselTbl);
    AE_DSEL16X4(x3, x4, z2, z3, dselTbl);
    AE_DSEL16X4(x6, x1, z4, z5, dselTbl);
    AE_DSEL16X4(x2, x5, z6, z7, dselTbl);

    /* First stage of FFT */
    SetDFT4_Scaling(FIRST_STAGE_SCALE);
    AE_L16X4X2_IP(tw0102, tw0311, p_twd, 2*sizeof(ae_int16x4));
    AE_L16X4X2_IP(tw1213, tw2122, p_twd, 2*sizeof(ae_int16x4));
    AE_L16X4X2_IP(tw2331, tw3233, p_twd, 2*sizeof(ae_int16x4));

    DFT4XI2(x0, x1, x2, x3);
    DFT4XI2(x4, x5, x6, x7);
    x1 = AE_MULFC16RAS(x1, tw0102);
    x2 = AE_MULFC16RAS(x2, tw0311);
    x3 = AE_MULFC16RAS(x3, tw1213);
    x5 = AE_MULFC16RAS(x5, tw2122);
    x6 = AE_MULFC16RAS(x6, tw2331);
    x7 = AE_MULFC16RAS(x7, tw3233);

    dselTbl = AE_L16X4_I((ae_int16x4 *)dsel_tbl_perm, sizeof(ae_int16x4));
    AE_DSEL16X4(z2, z0, x0, x1, dselTbl);
    AE_DSEL16X4(z3, z1, x2, x3, dselTbl);
    AE_DSEL16X4(z6, z4, x4, x5, dselTbl);
    AE_DSEL16X4(z7, z5, x6, x7, dselTbl);

    /* Last stage of FFT */
    SetDFT4_Scaling(LAST_STAGE_SCALE);

    x0 = z0; x4 = z1;
    x1 = z2; x5 = z3;
    x2 = z4; x6 = z5;
    x3 = z6; x7 = z7;
    DFT4XI2(x0, x1, x2, x3);
    DFT4XI2(x4, x5, x6, x7);

    /* Save data */
    AE_S16X4X2_IP(x0, x4, py0, 8*sizeof(int16_t));
    AE_S16X4X2_IP(x1, x5, py0, 8*sizeof(int16_t));
    AE_S16X4X2_IP(x2, x6, py0, 8*sizeof(int16_t));
    AE_S16X4X2_IP(x3, x7, py0, 8*sizeof(int16_t));
}




/*
   scaled fft with reordering
   N=32 - size of FFT
   NOTE: y is input and output
*/
void fft32_16x16(int16_t *y, const int16_t *ptwd)
{
    //const int N = 32;
    const ae_int16x8 * restrict px0;
          ae_int16x8 * restrict py0;
    const ae_int16x8 * restrict p_twd;
    ae_int16x4 x0, x1, x2, x3, x4, x5, x6, x7;
    ae_int16x4 x8, x9, xA, xB, xC, xD, xE, xF;
    ae_int16x4 z0, z1, z2, z3, z4, z5, z6, z7;
    ae_int16x4 z8, z9, zA, zB, zC, zD, zE, zF;
    ae_int16x4 dselTbl;

    NASSERT_ALIGN(y, 16);

    p_twd = (const ae_int16x8 *)ptwd;
    px0 = (const ae_int16x8 *)y;
    py0 = (      ae_int16x8 *)y;
    /* First stage of FFT, radix-4 */
    {
        ae_int16x4 tw0102, tw0311, tw1213, tw2122, tw2331, tw3233;
        ae_int16x4 tw4142, tw4351, tw5253, tw6162, tw6371, tw7273;

        SetDFT4_Scaling(FIRST_STAGE_SCALE); 

        /* Load data */
        AE_L16X4X2_IP(z0, z1, px0, 2*sizeof(ae_int16x4));
        AE_L16X4X2_IP(z2, z3, px0, 2*sizeof(ae_int16x4));
        AE_L16X4X2_IP(z4, z5, px0, 2*sizeof(ae_int16x4));
        AE_L16X4X2_IP(z6, z7, px0, 2*sizeof(ae_int16x4));
        AE_L16X4X2_IP(z8, z9, px0, 2*sizeof(ae_int16x4));
        AE_L16X4X2_IP(zA, zB, px0, 2*sizeof(ae_int16x4));
        AE_L16X4X2_IP(zC, zD, px0, 2*sizeof(ae_int16x4));
        AE_L16X4X2_IP(zE, zF, px0, 2*sizeof(ae_int16x4));
        /* Reordering */
        dselTbl = AE_L16X4_I((ae_int16x4 *)dsel_tbl_perm, 0);
        AE_DSEL16X4(xF, x0, z0, z1, dselTbl);
        AE_DSEL16X4(xB, x4, z2, z3, dselTbl);
        AE_DSEL16X4(x7, x8, z4, z5, dselTbl);
        AE_DSEL16X4(x3, xC, z6, z7, dselTbl);
        AE_DSEL16X4(xE, x1, z8, z9, dselTbl);
        AE_DSEL16X4(xA, x5, zA, zB, dselTbl);
        AE_DSEL16X4(x6, x9, zC, zD, dselTbl);
        AE_DSEL16X4(x2, xD, zE, zF, dselTbl);

        DFT4XI2(x0, x1, x2, x3);
        DFT4XI2(x4, x5, x6, x7);
        DFT4XI2(x8, x9, xA, xB);
        DFT4XI2(xC, xD, xE, xF);

        AE_L16X4X2_IP(tw0102, tw0311, p_twd, 2*sizeof(ae_int16x4));
        AE_L16X4X2_IP(tw1213, tw2122, p_twd, 2*sizeof(ae_int16x4));
        AE_L16X4X2_IP(tw2331, tw3233, p_twd, 2*sizeof(ae_int16x4));
        AE_L16X4X2_IP(tw4142, tw4351, p_twd, 2*sizeof(ae_int16x4));
        AE_L16X4X2_IP(tw5253, tw6162, p_twd, 2*sizeof(ae_int16x4));
        AE_L16X4X2_IP(tw6371, tw7273, p_twd, 2*sizeof(ae_int16x4));

        x1 = AE_MULFC16RAS(x1, tw0102);
        x2 = AE_MULFC16RAS(x2, tw0311);
        x3 = AE_MULFC16RAS(x3, tw1213);
        x5 = AE_MULFC16RAS(x5, tw2122);
        x6 = AE_MULFC16RAS(x6, tw2331);
        x7 = AE_MULFC16RAS(x7, tw3233);
        x9 = AE_MULFC16RAS(x9, tw4142);
        xA = AE_MULFC16RAS(xA, tw4351);
        xB = AE_MULFC16RAS(xB, tw5253);
        xD = AE_MULFC16RAS(xD, tw6162);
        xE = AE_MULFC16RAS(xE, tw6371);
        xF = AE_MULFC16RAS(xF, tw7273);

        dselTbl = AE_L16X4_I((ae_int16x4 *)dsel_tbl_perm, sizeof(ae_int16x4));
        AE_DSEL16X4(z2, z0, x0, x1, dselTbl);
        AE_DSEL16X4(z3, z1, x2, x3, dselTbl);
        AE_DSEL16X4(z6, z4, x4, x5, dselTbl);
        AE_DSEL16X4(z7, z5, x6, x7, dselTbl);
        AE_DSEL16X4(zA, z8, x8, x9, dselTbl);
        AE_DSEL16X4(zB, z9, xA, xB, dselTbl);
        AE_DSEL16X4(zE, zC, xC, xD, dselTbl);
        AE_DSEL16X4(zF, zD, xE, xF, dselTbl);

    }
    /* Last stage of FFT, radix-8 */
    {
        ae_int16x4 s0, s1, s2, s3, s4, s5, s6, s7;
        ae_int16x4 d0, d1, d2, d3, d4, d5, d6, d7;
        ae_int16x4 t1, t2, t3;

        SetDFT4_Scaling(LAST_STAGE_SCALE); 

        x0 = z0; x8 = z1;
        x1 = z2; x9 = z3;
        x2 = z4; xA = z5;
        x3 = z6; xB = z7;
        x4 = z8; xC = z9;
        x5 = zA; xD = zB;
        x6 = zC; xE = zD;
        x7 = zE; xF = zF;

        DFT4XI2(x0, x2, x4, x6); 
        DFT4XI2(x1, x3, x5, x7);
        DFT4XI2(x8, xA, xC, xE);
        DFT4XI2(x9, xB, xD, xF);

        t1 = AE_L16X4_I((ae_int16x4*)__fft8_tw1_v2_,              0);
        t2 = AE_L16X4_I((ae_int16x4*)__fft8_tw1_v2_,     sizeof(t1));
        t3 = AE_L16X4_I((ae_int16x4*)__fft8_tw1_v2_, 2 * sizeof(t1));
        x3 = AE_MULFC16RAS(x3, t1); 
        x5 = AE_MULFC16RAS(x5, t2);
        x7 = AE_MULFC16RAS(x7, t3);
        xB = AE_MULFC16RAS(xB, t1);
        xD = AE_MULFC16RAS(xD, t2);
        xF = AE_MULFC16RAS(xF, t3);

        s0 = x0;    s1 = x2;    s2 = x4;    s3 = x6;
        s4 = x8;    s5 = xA;    s6 = xC;    s7 = xE;
        d0 = x1;    d1 = x3;    d2 = x5;    d3 = x7;
        d4 = x9;    d5 = xB;    d6 = xD;    d7 = xF;

        AE_ADDANDSUBRNG16RAS_S1(s0, d0); 
        AE_ADDANDSUBRNG16RAS_S1(s1, d1);
        AE_ADDANDSUBRNG16RAS_S1(s2, d2);
        AE_ADDANDSUBRNG16RAS_S1(s3, d3);
        AE_ADDANDSUBRNG16RAS_S1(s4, d4);
        AE_ADDANDSUBRNG16RAS_S1(s5, d5);
        AE_ADDANDSUBRNG16RAS_S1(s6, d6);
        AE_ADDANDSUBRNG16RAS_S1(s7, d7);

        x0 = s0;    x4 = d0;    x8 = s4;    xC = d4;
        x1 = s1;    x5 = d1;    x9 = s5;    xD = d5;
        x2 = s2;    x6 = d2;    xA = s6;    xE = d6;
        x3 = s3;    x7 = d3;    xB = s7;    xF = d7;

        /* Save data */
        AE_S16X4X2_IP(x0, x8, py0, 2*sizeof(ae_int16x4));
        AE_S16X4X2_IP(x1, x9, py0, 2*sizeof(ae_int16x4));
        AE_S16X4X2_IP(x2, xA, py0, 2*sizeof(ae_int16x4));
        AE_S16X4X2_IP(x3, xB, py0, 2*sizeof(ae_int16x4));
        AE_S16X4X2_IP(x4, xC, py0, 2*sizeof(ae_int16x4));
        AE_S16X4X2_IP(x5, xD, py0, 2*sizeof(ae_int16x4));
        AE_S16X4X2_IP(x6, xE, py0, 2*sizeof(ae_int16x4));
        AE_S16X4X2_IP(x7, xF, py0, 2*sizeof(ae_int16x4));
    }
}
