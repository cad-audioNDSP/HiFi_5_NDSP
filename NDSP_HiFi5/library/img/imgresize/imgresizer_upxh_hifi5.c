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
#include "NatureDSP_types.h"
#include "NatureDSP_Signal_img.h"
#include "common.h"
#include "img_common.h"
#include "img_getCoef_up.h"
#include "imgresizer_bilinear_common.h"

/*    image resizer, upsample 1...2x in horizontal direction */

/* returns size of coefficients */
static size_t getCoefSz(const imgsize_t* in,const imgsize_t* out)
{
    int win=in->width,wout=out->width;
    (void)in,(void)out;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    return img_getCoef_up_alloc(win,wout);
}
/* returns coefficients */
static void getCoef(void* coef, const imgsize_t* in,const imgsize_t* out)
{
    int win=in->width,wout=out->width;
    (void)in,(void)out,(void)coef;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    img_getCoef_up_init(coef,win,wout);
}

static size_t getScratchSize(const imgsize_t* in, const imgsize_t* out)
{
    (void)in, (void)out;
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    /* img size */
    return ((1 + (in->height >> 3))* (((in->width + 7) & ~7) * sizeof(ae_int16x8)));;
}

/* in-place image resize */
static void process(void* pScr, void* pCoef, const void* img, void* imgOut, const imgsize_t* in, const imgsize_t* out, int fast)
{
    int16_t* restrict t = (int16_t*)pScr;
    const int16_t* restrict w = ((img_coefup_t*)pCoef)->coef;
    const int16_t* restrict left = ((img_coefup_t*)pCoef)->left;
    const int16_t* restrict pL;
    const ae_int16* restrict pW;
    const ae_int16x8* restrict pT0;
    const ae_int16x8* restrict pT1;
    const ae_int16x8* restrict pT2;
    const ae_int16x8* restrict pT3;
    const ae_int16x8* restrict pIn;
    ae_int16x4* restrict pOut;
    ae_int16x4 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15;
    ae_int16x4 w0, w1, w2, w3, w4, w5, w6, w7;
    ae_int16x4 y0, y1, y2, y3, y4, y5, y6, y7;
    ae_int16x4 sel;
    ae_int32x2 Y0;
    int m, n, k0, k1, k2, k3, shift;
    int hin = in->height,
        win = in->width,
        wout = out->width,
        stride = out->stride;

    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    NASSERT_ALIGN(img, ALIGNMENT);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->height == out->height && wout<2 * win && wout>win && in->stride == out->stride);

    pW = (const ae_int16*)w;
    pL = left;
    sel = AE_MOVINT16X4_FROMINT64(0x0706050403020100); // SEL7531 + SEL6420
    win = (win + 7) & ~7;
    stride = stride * sizeof(int16_t);
    shift = (win) * sizeof(ae_int16x8);
    pOut = (ae_int16x4*)t;
    for (n = 0; n < (hin & ~7); n += 8)
    {
        /* Interleave input samples from 8 rows and save them to the scratch */
        pIn = (const ae_int16x8*)((uintptr_t)img + stride * n);
        pOut = (ae_int16x4*)((int16_t*)t + win * n);
        //pOut = (ae_int16x4*)((uintptr_t*)t+(shift));
        for (m = 0; m < win; m += 8)
        {
            AE_L16X4X2_XP(x0, x8, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x1, x9, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x2, x10, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x3, x11, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x4, x12, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x5, x13, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x6, x14, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x7, x15, castxcc(ae_int16x8, pIn), -7 * stride + sizeof(ae_int16x8));

            AE_DSEL16X4(y0, y1, x0, x1, sel);
            AE_DSEL16X4(y2, y3, x2, x3, sel);
            AE_DSEL16X4(x0, x2, y0, y2, sel);
            AE_DSEL16X4(x1, x3, y1, y3, sel);

            AE_DSEL16X4(y0, y1, x4, x5, sel);
            AE_DSEL16X4(y2, y3, x6, x7, sel);
            AE_DSEL16X4(x4, x6, y0, y2, sel);
            AE_DSEL16X4(x5, x7, y1, y3, sel);

            AE_DSEL16X4(y0, y1, x8, x9, sel);
            AE_DSEL16X4(y2, y3, x10, x11, sel);
            AE_DSEL16X4(x8, x10, y0, y2, sel);
            AE_DSEL16X4(x9, x11, y1, y3, sel);

            AE_DSEL16X4(y0, y1, x12, x13, sel);
            AE_DSEL16X4(y2, y3, x14, x15, sel);
            AE_DSEL16X4(x12, x14, y0, y2, sel);
            AE_DSEL16X4(x13, x15, y1, y3, sel);

            AE_S16X4X2_IP(x0, x4, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x1, x5, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x2, x6, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x3, x7, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x8, x12, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x9, x13, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x10, x14, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x11, x15, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));

        }


    }
    pIn = (ae_int16x8*)t;
    for (m = 0; m < wout; m += 4)
    {
        pOut = (ae_int16x4*)((int16_t*)img + m);
        AE_L16_IP(w0, pW, sizeof(int16_t));
        AE_L16_IP(w1, pW, sizeof(int16_t));
        AE_L16_IP(w2, pW, sizeof(int16_t));
        AE_L16_IP(w3, pW, sizeof(int16_t));
        AE_L16_IP(w4, pW, sizeof(int16_t));
        AE_L16_IP(w5, pW, sizeof(int16_t));
        AE_L16_IP(w6, pW, sizeof(int16_t));
        AE_L16_IP(w7, pW, sizeof(int16_t));

        k0 = *pL++;  k1 = *pL++;  k2 = *pL++;  k3 = *pL++;

        pT0 = pIn + k0;
        pT1 = pIn + k1;
        pT2 = pIn + k2;
        pT3 = pIn + k3;


        for (n = 0; n < (hin & ~7); n += 8)
        {
            AE_L16X4X2_I(x2, x3, pT0, sizeof(ae_int16x8));
            AE_L16X4X2_XP(x0, x1, pT0, shift);

            AE_L16X4X2_I(x6, x7, pT1, sizeof(ae_int16x8));
            AE_L16X4X2_XP(x4, x5, pT1, shift);

            AE_L16X4X2_I(x10, x11, pT2, sizeof(ae_int16x8));
            AE_L16X4X2_XP(x8, x9, pT2, shift);

            AE_L16X4X2_I(x14, x15, pT3, sizeof(ae_int16x8));
            AE_L16X4X2_XP(x12, x13, pT3, shift);

            y0 = AE_MULFD16X16X4RAS(x0, x2, w0, w1);
            y4 = AE_MULFD16X16X4RAS(x1, x3, w0, w1);

            y1 = AE_MULFD16X16X4RAS(x4, x6, w2, w3);
            y5 = AE_MULFD16X16X4RAS(x5, x7, w2, w3);

            y2 = AE_MULFD16X16X4RAS(x8, x10, w4, w5);
            y6 = AE_MULFD16X16X4RAS(x9, x11, w4, w5);

            y3 = AE_MULFD16X16X4RAS(x12, x14, w6, w7);
            y7 = AE_MULFD16X16X4RAS(x13, x15, w6, w7);

            AE_DSEL16X4(x0, x1, y0, y1, sel);
            AE_DSEL16X4(x2, x3, y2, y3, sel);
            AE_DSEL16X4(y0, y2, x0, x2, sel);
            AE_DSEL16X4(y1, y3, x1, x3, sel);

            AE_DSEL16X4(x0, x1, y4, y5, sel);
            AE_DSEL16X4(x2, x3, y6, y7, sel);
            AE_DSEL16X4(y4, y6, x0, x2, sel);
            AE_DSEL16X4(y5, y7, x1, x3, sel);

            AE_S16X4_XP(y0, pOut, stride);
            AE_S16X4_XP(y1, pOut, stride);
            AE_S16X4_XP(y2, pOut, stride);
            AE_S16X4_XP(y3, pOut, stride);

            AE_S16X4_XP(y4, pOut, stride);
            AE_S16X4_XP(y5, pOut, stride);
            AE_S16X4_XP(y6, pOut, stride);
            AE_S16X4_XP(y7, pOut, stride);
        }
    }

    /* Process the image by 4 output rows per iteration */
    if ((hin & 7) >= 4)
    {
        n = (hin & ~7);
        /* Interleave input samples from 4 rows and save them to the scratch */
        pIn = (const ae_int16x8*)((uintptr_t)img + n * stride);
        pOut = (ae_int16x4*)t;
        for (m = 0; m < win; m += 8)
        {
            AE_L16X4X2_XP(x0, x4, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x1, x5, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x2, x6, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x3, x7, castxcc(ae_int16x8, pIn), -3 * stride + sizeof(ae_int16x8));

            AE_DSEL16X4(y0, y1, x0, x1, sel);
            AE_DSEL16X4(y2, y3, x2, x3, sel);
            AE_DSEL16X4(x0, x2, y0, y2, sel);
            AE_DSEL16X4(x1, x3, y1, y3, sel);

            AE_DSEL16X4(y0, y1, x4, x5, sel);
            AE_DSEL16X4(y2, y3, x6, x7, sel);
            AE_DSEL16X4(x4, x6, y0, y2, sel);
            AE_DSEL16X4(x5, x7, y1, y3, sel);

            AE_S16X4X2_IP(x0, x1, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x2, x3, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x4, x5, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x6, x7, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
        }
        pIn = (const ae_int16x8*)t; /* load input samples from scratch */
        pOut = (ae_int16x4*)((uintptr_t)img + n * stride);
        pW = (const ae_int16*)w;
        /* Process by 4x4 output samples per iteration */
        for (m = 0; m < wout; m += 4)
        {
            /* Load indexes */
            k0 = left[m + 0];
            k1 = left[m + 1];
            k2 = left[m + 2];
            k3 = left[m + 3];
            /* Load 2x4 elems for each of 4 rows */

            x0 = AE_L16X4_X((const ae_int16x4*)pIn, k0 * sizeof(ae_int16x4));
            x1 = AE_L16X4_X((const ae_int16x4*)pIn, k0 * sizeof(ae_int16x4) + sizeof(ae_int16x4));
            x2 = AE_L16X4_X((const ae_int16x4*)pIn, k1 * sizeof(ae_int16x4));
            x3 = AE_L16X4_X((const ae_int16x4*)pIn, k1 * sizeof(ae_int16x4) + sizeof(ae_int16x4));
            x4 = AE_L16X4_X((const ae_int16x4*)pIn, k2 * sizeof(ae_int16x4));
            x5 = AE_L16X4_X((const ae_int16x4*)pIn, k2 * sizeof(ae_int16x4) + sizeof(ae_int16x4));
            x6 = AE_L16X4_X((const ae_int16x4*)pIn, k3 * sizeof(ae_int16x4));
            x7 = AE_L16X4_X((const ae_int16x4*)pIn, k3 * sizeof(ae_int16x4) + sizeof(ae_int16x4));
            /* Load window coefficients */
            AE_L16_IP(w0, pW, sizeof(int16_t));
            AE_L16_IP(w1, pW, sizeof(int16_t));
            AE_L16_IP(w2, pW, sizeof(int16_t));
            AE_L16_IP(w3, pW, sizeof(int16_t));
            AE_L16_IP(w4, pW, sizeof(int16_t));
            AE_L16_IP(w5, pW, sizeof(int16_t));
            AE_L16_IP(w6, pW, sizeof(int16_t));
            AE_L16_IP(w7, pW, sizeof(int16_t));

            y0 = AE_MULFD16X16X4RAS(x0, x1, w0, w1);
            y1 = AE_MULFD16X16X4RAS(x2, x3, w2, w3);
            y2 = AE_MULFD16X16X4RAS(x4, x5, w4, w5);
            y3 = AE_MULFD16X16X4RAS(x6, x7, w6, w7);

            AE_DSEL16X4(x0, x1, y0, y1, sel);
            AE_DSEL16X4(x2, x3, y2, y3, sel);
            AE_DSEL16X4(y0, y2, x0, x2, sel);
            AE_DSEL16X4(y1, y3, x1, x3, sel);

            AE_S16X4_XP(y0, pOut, stride);
            AE_S16X4_XP(y1, pOut, stride);
            AE_S16X4_XP(y2, pOut, stride);
            AE_S16X4_XP(y3, pOut, -3 * stride + sizeof(ae_int16x4));
        }

    }

    /* Process last 0...3 rows */
    for (n = (hin & ~3); n < hin; n++)
    {
        /* Move input samples to the scratch */
        pIn = (const ae_int16x8*)((uintptr_t)img + n * stride);
        pOut = (ae_int16x4*)t;
        for (m = 0; m < win; m += 8)
        {
            AE_L16X4X2_IP(x0, x1, pIn, sizeof(ae_int16x8));
            AE_S16X4X2_IP(x0, x1, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
        }

        pIn = (const ae_int16x8*)t;
        pOut = (ae_int16x4*)((uintptr_t)img + n * stride);
        pW = (const ae_int16*)w;
        /* Process by 1 output sample per iteration */
        for (m = 0; m < wout; m++)
        {
            k0 = left[m];
            AE_L16_IP(w0, pW, sizeof(int16_t));
            AE_L16_IP(w1, pW, sizeof(int16_t));
            x0 = AE_L16_X((const ae_int16*)pIn, k0 * sizeof(int16_t));
            x1 = AE_L16_X((const ae_int16*)pIn, k0 * sizeof(int16_t) + sizeof(int16_t));
            Y0 = AE_MULF16SS_00(x0, w0);
            AE_MULAF16SS_00(Y0, x1, w1);
            y0 = AE_ROUND16X4F32SASYM(Y0, Y0);
            AE_S16_0_IP(y0, castxcc(ae_int16, pOut), sizeof(int16_t));
        }
    }

}

const imgresizer_api_t imgresizer_api_upxh={"imgresizer_api_upxh",getCoefSz,getCoef,getScratchSize,process};
