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

#ifndef AE_DSEL32X2_HH_LLSWP
/*
   Equal to:
   a = AE_SEL32_HH(c, d)
   b = AE_SEL32_LL(d, c)
*/
#define AE_DSEL32X2_HH_LLSWP(a,b,c,d) \
{\
    ae_int16x4 aa,bb,cc,dd,sel; \
    sel = AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32X2((7L<<24)|(1L<<16)|(6<<8)|(0), (3L<<24)|(5L<<16)|(2<<8)|(4) )); \
    cc = AE_MOVINT16X4_FROMINT32X2(c); \
    dd = AE_MOVINT16X4_FROMINT32X2(d); \
    AE_DSEL16X4(aa,bb,cc,dd,sel); \
    a = AE_MOVINT32X2_FROMINT16X4(aa); \
    b = AE_MOVINT32X2_FROMINT16X4(bb); \
}
#endif

static const int32_t ALIGN(32) __fft8_tw1[] =
{
  (int32_t)0x00000000, (int32_t)0x80000000,
  (int32_t)0xA57D8666, (int32_t)0xA57D8666,
};

#define DFT4X1RNG(x0, x1, x2, x3)\
{   \
    ae_int32x2 s0, s1, d0, d1;         \
    AE_ADDANDSUBRNG32(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32(s1, d1, x1, x3); \
    AE_ADDANDSUB32S(x0, x2, s0, s1);   \
    AE_ADDANDSUB32JS(x3, x1, d0, d1);  \
}

/* radix-4 butterfly with normalization  */
#define DFT4X1RNG_H(x0, x1, x2, x3) \
{ \
    ae_int32x2 s0, s1, d0, d1;           \
    AE_ADDANDSUBRNG32_H(s0, d0, x0, x2); \
    AE_ADDANDSUBRNG32_H(s1, d1, x1, x3); \
    AE_ADDANDSUB32S(x0, x2, s0, s1);     \
    AE_ADDANDSUB32JS(x3, x1, d0, d1);    \
}

/*
   scaled fft with reordering
   N=16 - size of FFT
   NOTE: y is input and output, x - temporary
*/
void fft16_32x32(int32_t *y, int32_t *x, const int32_t *tw)
{
    const ae_int32x4 * restrict px0;
          ae_int32x4 * restrict py0;
    const ae_int32x4 * restrict ptwd;
    ae_int32x2 x0, x1, x2, x3, x4, x5, x6, x7;
    ae_int32x2 x8, x9, xA, xB, xC, xD, xE, xF;
    ae_int32x2 y0, y1, y2, y3, y4, y5, y6, y7;
    ae_int32x2 y8, y9, yA, yB, yC, yD, yE, yF;
    ae_int32x2 tw1, tw2, tw3, tw5, tw6, tw7;
    ae_int32x2 tw9, twA, twB, twD, twE, twF;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);

    px0  = (const ae_int32x4 *)y;
    py0  = (      ae_int32x4 *)y;
    ptwd = (const ae_int32x4 *)tw;

    /*
     * Load data
     */
    AE_L32X2X2_IP(x0, x1, px0, 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(x2, x3, px0, 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(x4, x5, px0, 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(x6, x7, px0, 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(x8, x9, px0, 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(xA, xB, px0, 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(xC, xD, px0, 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(xE, xF, px0, 2*sizeof(ae_int32x2));
    /* Reordering */
    AE_DSEL32X2_HH_LLSWP(y0, y1, x0, x1);
    AE_DSEL32X2_HH_LLSWP(y2, y3, x2, x3);
    AE_DSEL32X2_HH_LLSWP(y4, y5, x4, x5);
    AE_DSEL32X2_HH_LLSWP(y6, y7, x6, x7);
    AE_DSEL32X2_HH_LLSWP(y8, y9, x8, x9);
    AE_DSEL32X2_HH_LLSWP(yA, yB, xA, xB);
    AE_DSEL32X2_HH_LLSWP(yC, yD, xC, xD);
    AE_DSEL32X2_HH_LLSWP(yE, yF, xE, xF);

    /*
     * 1st stage of FFT
     */
    x0 = y0;    x1 = y8;    x2 = yF;    x3 = y7;
    x4 = y2;    x5 = yA;    x6 = yD;    x7 = y5;
    x8 = y4;    x9 = yC;    xA = yB;    xB = y3;
    xC = y6;    xD = yE;    xE = y9;    xF = y1;

    WUR_AE_SAR(3);
    DFT4X1RNG(x0, x1, x2, x3);
    DFT4X1RNG(x4, x5, x6, x7);
    DFT4X1RNG(x8, x9, xA, xB);
    DFT4X1RNG(xC, xD, xE, xF);

    AE_L32X2X2_IP(tw1, tw2, ptwd, 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(tw3, tw5, ptwd, 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(tw6, tw7, ptwd, 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(tw9, twA, ptwd, 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(twB, twD, ptwd, 2*sizeof(ae_int32x2));
    AE_L32X2X2_IP(twE, twF, ptwd, 2*sizeof(ae_int32x2));
    x1 = AE_MULFC32RAS(x1, tw1);
    x2 = AE_MULFC32RAS(x2, tw2);
    x3 = AE_MULFC32RAS(x3, tw3);
    x5 = AE_MULFC32RAS(x5, tw5);
    x6 = AE_MULFC32RAS(x6, tw6);
    x7 = AE_MULFC32RAS(x7, tw7);
    x9 = AE_MULFC32RAS(x9, tw9);
    xA = AE_MULFC32RAS(xA, twA);
    xB = AE_MULFC32RAS(xB, twB);
    xD = AE_MULFC32RAS(xD, twD);
    xE = AE_MULFC32RAS(xE, twE);
    xF = AE_MULFC32RAS(xF, twF);

    /*
     * Last stage of FFT
     */
    y0 = x0;    y1 = x4;    y2 = x8;    y3 = xC;
    y4 = x1;    y5 = x5;    y6 = x9;    y7 = xD;
    y8 = x2;    y9 = x6;    yA = xA;    yB = xE;
    yC = x3;    yD = x7;    yE = xB;    yF = xF;

    WUR_AE_SAR(2);
    DFT4X1RNG(y0, y1, y2, y3);
    DFT4X1RNG(y4, y5, y6, y7);
    DFT4X1RNG(y8, y9, yA, yB);
    DFT4X1RNG(yC, yD, yE, yF);

    /*
     * Save data
     */
    AE_S32X2X2_IP(y0, y4, py0, 2*sizeof(ae_int32x2));
    AE_S32X2X2_IP(y8, yC, py0, 2*sizeof(ae_int32x2));
    AE_S32X2X2_IP(y1, y5, py0, 2*sizeof(ae_int32x2));
    AE_S32X2X2_IP(y9, yD, py0, 2*sizeof(ae_int32x2));
    AE_S32X2X2_IP(y2, y6, py0, 2*sizeof(ae_int32x2));
    AE_S32X2X2_IP(yA, yE, py0, 2*sizeof(ae_int32x2));
    AE_S32X2X2_IP(y3, y7, py0, 2*sizeof(ae_int32x2));
    AE_S32X2X2_IP(yB, yF, py0, 2*sizeof(ae_int32x2));
}




/*
   scaled fft with reordering
   N=32 - size of FFT
   NOTE: y is input and output, x - temporary
*/
void fft32_32x32(int32_t *y, int32_t *x, const int32_t *tw)
{
    #define N 32
    int k;
    ae_int32x4 * restrict px0;
    ae_int32x4 * restrict px1;
    ae_int32x4 * restrict px2;
    ae_int32x4 * restrict px3;
    ae_int32x4 * restrict py0;
    ae_int32x4 * restrict py1;
    const ae_int32x4 * restrict ptwd;
    ae_int32x2 x0, x1, x2, x3, x4, x5, x6, x7;
    ae_int32x2 x8, x9, xA, xB, xC, xD, xE, xF;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);

    /*
     * Reorder data
     */
    {
        ae_int32x2 y0, y1, y2, y3;
        px0  = (ae_int32x4 *)(x);
        px1  = (ae_int32x4 *)(x+2*N-4);
        py0  = (ae_int32x4 *)(y);

    	for (k=0; k<(N>>2); k++) 
        { 
            AE_L32X2X2_IP(x0, x1, py0, 2*sizeof(ae_int32x2));
            AE_L32X2X2_IP(x2, x3, py0, 2*sizeof(ae_int32x2));
            AE_DSEL32X2_HH_LLSWP(y0, y1, x0, x1);
            AE_DSEL32X2_HH_LLSWP(y2, y3, x2, x3);
            AE_S32X2X2_IP(y0, y2, px0, 2*sizeof(ae_int32x2));
            AE_S32X2X2_IP(y3, y1, px1, -2*(int)sizeof(ae_int32x2));
        }
    }
    __Pragma("no_reorder");
    /*
     * The first stage of FFT, radix-4
     */
    {
        #define stride (N/4)
        ae_int32x2 tw1, tw2, tw3, tw5, tw6, tw7;
        ae_int32x2 tw9, twA, twB, twD, twE, twF;

        px0 = (ae_int32x4 *)x;
        px1 = (ae_int32x4 *)px0 + stride/2;
        px2 = (ae_int32x4 *)px1 + stride/2;
        px3 = (ae_int32x4 *)px2 + stride/2;
        py0 = (ae_int32x4 *)y;
        ptwd = (const ae_int32x4 *)tw;

        WUR_AE_SAR(3);

        /* first 4 butterflies */
        AE_L32X2X2_IP(x0, x4, px0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x1, x5, px1, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x2, x6, px2, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x3, x7, px3, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x8, xC, px0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x9, xD, px1, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(xA, xE, px2, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(xB, xF, px3, 2*sizeof(ae_int32x2));

        DFT4X1RNG(x0, x1, x2, x3);
        DFT4X1RNG(x4, x5, x6, x7);
        DFT4X1RNG(x8, x9, xA, xB);
        DFT4X1RNG(xC, xD, xE, xF);

        AE_L32X2X2_IP(tw1, tw2, ptwd, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(tw3, tw5, ptwd, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(tw6, tw7, ptwd, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(tw9, twA, ptwd, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(twB, twD, ptwd, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(twE, twF, ptwd, 2*sizeof(ae_int32x2));
        x1 = AE_MULFC32RAS(x1, tw1);
        x2 = AE_MULFC32RAS(x2, tw2);
        x3 = AE_MULFC32RAS(x3, tw3);
        x5 = AE_MULFC32RAS(x5, tw5);
        x6 = AE_MULFC32RAS(x6, tw6);
        x7 = AE_MULFC32RAS(x7, tw7);
        x9 = AE_MULFC32RAS(x9, tw9);
        xA = AE_MULFC32RAS(xA, twA);
        xB = AE_MULFC32RAS(xB, twB);
        xD = AE_MULFC32RAS(xD, twD);
        xE = AE_MULFC32RAS(xE, twE);
        xF = AE_MULFC32RAS(xF, twF);

        AE_S32X2X2_IP(x0, x1, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x2, x3, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x4, x5, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x6, x7, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x8, x9, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(xA, xB, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(xC, xD, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(xE, xF, py0, 2*sizeof(ae_int32x2));

        /* next 4 butterflies */
        AE_L32X2X2_IP(x0, x4, px0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x1, x5, px1, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x2, x6, px2, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x3, x7, px3, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x8, xC, px0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x9, xD, px1, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(xA, xE, px2, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(xB, xF, px3, 2*sizeof(ae_int32x2));

        DFT4X1RNG(x0, x1, x2, x3);
        DFT4X1RNG(x4, x5, x6, x7);
        DFT4X1RNG(x8, x9, xA, xB);
        DFT4X1RNG(xC, xD, xE, xF);

        AE_L32X2X2_IP(tw1, tw2, ptwd, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(tw3, tw5, ptwd, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(tw6, tw7, ptwd, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(tw9, twA, ptwd, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(twB, twD, ptwd, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(twE, twF, ptwd, 2*sizeof(ae_int32x2));
        x1 = AE_MULFC32RAS(x1, tw1);
        x2 = AE_MULFC32RAS(x2, tw2);
        x3 = AE_MULFC32RAS(x3, tw3);
        x5 = AE_MULFC32RAS(x5, tw5);
        x6 = AE_MULFC32RAS(x6, tw6);
        x7 = AE_MULFC32RAS(x7, tw7);
        x9 = AE_MULFC32RAS(x9, tw9);
        xA = AE_MULFC32RAS(xA, twA);
        xB = AE_MULFC32RAS(xB, twB);
        xD = AE_MULFC32RAS(xD, twD);
        xE = AE_MULFC32RAS(xE, twE);
        xF = AE_MULFC32RAS(xF, twF);

        AE_S32X2X2_IP(x0, x1, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x2, x3, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x4, x5, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x6, x7, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x8, x9, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(xA, xB, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(xC, xD, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(xE, xF, py0, 2*sizeof(ae_int32x2));
        #undef stride
    }
    __Pragma("no_reorder");
    /*
     * Last stage of FFT, radix-8
     */
    {
        #define stride (N/8)
        ae_int32x2 s0, s1, s2, s3, s4, s5, s6, s7;
        ae_int32x2 d0, d1, d2, d3, d4, d5, d6, d7;
        ae_int32x2 tw1, tw2;

        px0 = (ae_int32x4 *)y;
        px1 = px0 + 2*stride;
        py0 = (ae_int32x4 *)y;
        py1 = py0 + 2*stride;
        ptwd = (const ae_int32x4 *)__fft8_tw1;
        AE_L32X2X2_I(tw1, tw2, ptwd, 0);
        AE_MOVSARA7X2(2, 1);

        /* first 2 butterflies */
        AE_L32X2X2_IP(x0, x8, px0, stride*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x1, x9, px0, stride*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x2, xA, px0, stride*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x3, xB, px0, (2-3*stride)*(int)sizeof(ae_int32x2));
        AE_L32X2X2_IP(x4, xC, px1, stride*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x5, xD, px1, stride*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x6, xE, px1, stride*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x7, xF, px1, (2-3*stride)*(int)sizeof(ae_int32x2));

        DFT4X1RNG_H(x0, x2, x4, x6); 
        DFT4X1RNG_H(x1, x3, x5, x7);
        DFT4X1RNG_H(x8, xA, xC, xE); 
        DFT4X1RNG_H(x9, xB, xD, xF);
        x3 = AE_MULFCJ32RAS(x3, tw2);
        x5 = AE_MULFC32RAS(x5, tw1);
        x7 = AE_MULFC32RAS(x7, tw2);
        xB = AE_MULFCJ32RAS(xB, tw2);
        xD = AE_MULFC32RAS(xD, tw1);
        xF = AE_MULFC32RAS(xF, tw2);

        AE_ADDANDSUBRNG32_L(s0, d0, x0, x1);
        AE_ADDANDSUBRNG32_L(d1, s1, x2, x3);
        AE_ADDANDSUBRNG32_L(s2, d2, x4, x5);
        AE_ADDANDSUBRNG32_L(s3, d3, x6, x7);
        AE_ADDANDSUBRNG32_L(s4, d4, x8, x9);
        AE_ADDANDSUBRNG32_L(d5, s5, xA, xB);
        AE_ADDANDSUBRNG32_L(s6, d6, xC, xD);
        AE_ADDANDSUBRNG32_L(s7, d7, xE, xF);
        x0 = s0;    x1 = s1;    x2 = s2;    x3 = s3;
        x4 = d0;    x5 = d1;    x6 = d2;    x7 = d3;
        x8 = s4;    x9 = s5;    xA = s6;    xB = s7;
        xC = d4;    xD = d5;    xE = d6;    xF = d7;

        AE_S32X2X2_IP(x0, x8, py0, stride*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x1, x9, py0, stride*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x2, xA, py0, stride*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x3, xB, py0, (2-3*stride)*(int)sizeof(ae_int32x2));
        AE_S32X2X2_IP(x4, xC, py1, stride*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x5, xD, py1, stride*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x6, xE, py1, stride*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x7, xF, py1, (2-3*stride)*(int)sizeof(ae_int32x2));

        /* next 2 butterflies */
        AE_L32X2X2_IP(x0, x8, px0, stride*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x1, x9, px0, stride*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x2, xA, px0, stride*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x3, xB, px0, (2-3*stride)*(int)sizeof(ae_int32x2));
        AE_L32X2X2_IP(x4, xC, px1, stride*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x5, xD, px1, stride*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x6, xE, px1, stride*sizeof(ae_int32x2));
        AE_L32X2X2_IP(x7, xF, px1, (2-3*stride)*(int)sizeof(ae_int32x2));

        DFT4X1RNG_H(x0, x2, x4, x6); 
        DFT4X1RNG_H(x1, x3, x5, x7);
        DFT4X1RNG_H(x8, xA, xC, xE); 
        DFT4X1RNG_H(x9, xB, xD, xF);
        x3 = AE_MULFCJ32RAS(x3, tw2);
        x5 = AE_MULFC32RAS(x5, tw1);
        x7 = AE_MULFC32RAS(x7, tw2);
        xB = AE_MULFCJ32RAS(xB, tw2);
        xD = AE_MULFC32RAS(xD, tw1);
        xF = AE_MULFC32RAS(xF, tw2);

        AE_ADDANDSUBRNG32_L(s0, d0, x0, x1);
        AE_ADDANDSUBRNG32_L(d1, s1, x2, x3);
        AE_ADDANDSUBRNG32_L(s2, d2, x4, x5);
        AE_ADDANDSUBRNG32_L(s3, d3, x6, x7);
        AE_ADDANDSUBRNG32_L(s4, d4, x8, x9);
        AE_ADDANDSUBRNG32_L(d5, s5, xA, xB);
        AE_ADDANDSUBRNG32_L(s6, d6, xC, xD);
        AE_ADDANDSUBRNG32_L(s7, d7, xE, xF);
        x0 = s0;    x1 = s1;    x2 = s2;    x3 = s3;
        x4 = d0;    x5 = d1;    x6 = d2;    x7 = d3;
        x8 = s4;    x9 = s5;    xA = s6;    xB = s7;
        xC = d4;    xD = d5;    xE = d6;    xF = d7;

        AE_S32X2X2_IP(x0, x8, py0, stride*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x1, x9, py0, stride*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x2, xA, py0, stride*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x3, xB, py0, (2-3*stride)*(int)sizeof(ae_int32x2));
        AE_S32X2X2_IP(x4, xC, py1, stride*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x5, xD, py1, stride*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x6, xE, py1, stride*sizeof(ae_int32x2));
        AE_S32X2X2_IP(x7, xF, py1, (2-3*stride)*(int)sizeof(ae_int32x2));
        #undef stride
    }
#undef N
}
