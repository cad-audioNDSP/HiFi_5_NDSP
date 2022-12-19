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
#include "img_getCoef_dn.h"
#include "imgresizer_bilinear_common.h"

/*    image resizer, downsample 1...2x in horizontal direction */

/* returns size of coefficients */
static size_t getCoefSz(const imgsize_t* in,const imgsize_t* out)
{
    int win=in->width,wout=out->width;
    NASSERT(in->height==out->height &&  win<2*wout && win>wout && in->stride==out->stride);
    (void)in,(void)out;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    return img_getCoef_dn_alloc(win,wout);
}
/* returns coefficients */
static void getCoef(void* coef, const imgsize_t* in,const imgsize_t* out)
{
    int win=in->width,wout=out->width;
    NASSERT(in->height==out->height &&  win<2*wout && win>wout && in->stride==out->stride);
    (void)in,(void)out,(void)coef;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    img_getCoef_dn_init(coef,win,wout);
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
    const int16_t* restrict w = ((img_coefdn_t*)pCoef)->coef;
    const int16_t* restrict left = ((img_coefdn_t*)pCoef)->left;
    const int16_t* restrict pL;
    const ae_int16* restrict pW;
    const ae_int16x8* restrict pT0;
    const ae_int16x8* restrict pT1;
    const ae_int16x8* restrict pT2;
    const ae_int16x8* restrict pT3;
    const ae_int16x8* restrict pIn;
    ae_int16x4* restrict pOut;
    ae_int16x4 x00, x01, x10, x11, x20, x21, x30, x31, x40, x41, x50, x51, x60, x61, x70, x71;
    ae_int16x4 w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15;
    ae_int16x4 y0, y1, y2, y3, y4, y5, y6, y7, t0, t1;
    ae_int16x4 sel;
    ae_int32x2 Y0, Y1;
    ae_int64 A0, A1, A2, A3;
    ae_valign alIn;
    int m, n, k0, k1, k2, k3, shift;
    int hin = in->height,
        win = in->width,
        wout = out->width,
        stride = out->stride;

    NASSERT_ALIGN(img, ALIGNMENT);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->height == out->height && in->width<2 * out->width && in->width>out->width && in->stride == out->stride);
    NASSERT(wout >= 1);
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);

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
            AE_L16X4X2_XP(x00, x01, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x10, x11, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x20, x21, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x30, x31, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x40, x41, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x50, x51, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x60, x61, castxcc(ae_int16x8, pIn), stride);
            AE_L16X4X2_XP(x70, x71, castxcc(ae_int16x8, pIn), -7 * stride + sizeof(ae_int16x8));

            AE_DSEL16X4(y0, y1, x00, x10, sel);
            AE_DSEL16X4(y2, y3, x20, x30, sel);
            AE_DSEL16X4(x00, x20, y0, y2, sel);
            AE_DSEL16X4(x10, x30, y1, y3, sel);

            AE_DSEL16X4(y0, y1, x40, x50, sel);
            AE_DSEL16X4(y2, y3, x60, x70, sel);
            AE_DSEL16X4(x40, x60, y0, y2, sel);
            AE_DSEL16X4(x50, x70, y1, y3, sel);

            AE_DSEL16X4(y0, y1, x01, x11, sel);
            AE_DSEL16X4(y2, y3, x21, x31, sel);
            AE_DSEL16X4(x01, x21, y0, y2, sel);
            AE_DSEL16X4(x11, x31, y1, y3, sel);

            AE_DSEL16X4(y0, y1, x41, x51, sel);
            AE_DSEL16X4(y2, y3, x61, x71, sel);
            AE_DSEL16X4(x41, x61, y0, y2, sel);
            AE_DSEL16X4(x51, x71, y1, y3, sel);

            AE_S16X4X2_IP(x00, x40, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x10, x50, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x20, x60, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x30, x70, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x01, x41, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x11, x51, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x21, x61, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x31, x71, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));

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

        AE_L16_IP(w8, pW, sizeof(int16_t));
        AE_L16_IP(w9, pW, sizeof(int16_t));
        AE_L16_IP(w10, pW, sizeof(int16_t));
        AE_L16_IP(w11, pW, sizeof(int16_t));
        AE_L16_IP(w12, pW, sizeof(int16_t));
        AE_L16_IP(w13, pW, sizeof(int16_t));
        AE_L16_IP(w14, pW, sizeof(int16_t));
        AE_L16_IP(w15, pW, sizeof(int16_t));


        k0 = *pL++;  k1 = *pL++;  k2 = *pL++;  k3 = *pL++;

        pT0 = pIn + k0;
        pT1 = pIn + k1;
        pT2 = pIn + k2;
        pT3 = pIn + k3;


        for (n = 0; n < (hin & ~7); n += 8)
        {
            /* AE_L16X4X2_X is used for better scheduling*/
            AE_L16X4X2_X(x60, x70, pT0, 3 * sizeof(ae_int16x8));
            AE_L16X4X2_X(x40, x50, pT0, 2 * sizeof(ae_int16x8));
            AE_L16X4X2_X(x20, x30, pT0, sizeof(ae_int16x8));
            AE_L16X4X2_XP(x00, x10, pT0, shift);


            t0 = AE_MULFD16X16X4RAS(x00, x20, w0, w1);
            t1 = AE_MULFD16X16X4RAS(x40, x60, w2, w3);
            y0 = AE_ADD16S(t0, t1);

            t0 = AE_MULFD16X16X4RAS(x10, x30, w0, w1);
            t1 = AE_MULFD16X16X4RAS(x50, x70, w2, w3);
            y4 = AE_ADD16S(t0, t1);

            /* AE_L16X4X2_X is used for better scheduling*/
            AE_L16X4X2_X(x60, x70, pT1, 3 * sizeof(ae_int16x8));
            AE_L16X4X2_X(x40, x50, pT1, 2 * sizeof(ae_int16x8));
            AE_L16X4X2_X(x20, x30, pT1, sizeof(ae_int16x8));
            AE_L16X4X2_XP(x00, x10, pT1, shift);


            t0 = AE_MULFD16X16X4RAS(x00, x20, w4, w5);
            t1 = AE_MULFD16X16X4RAS(x40, x60, w6, w7);
            y1 = AE_ADD16S(t0, t1);

            t0 = AE_MULFD16X16X4RAS(x10, x30, w4, w5);
            t1 = AE_MULFD16X16X4RAS(x50, x70, w6, w7);
            y5 = AE_ADD16S(t0, t1);

            /* AE_L16X4X2_X is used for better scheduling*/
            AE_L16X4X2_X(x60, x70, pT2, 3 * sizeof(ae_int16x8));
            AE_L16X4X2_X(x40, x50, pT2, 2 * sizeof(ae_int16x8));
            AE_L16X4X2_X(x20, x30, pT2, sizeof(ae_int16x8));
            AE_L16X4X2_XP(x00, x10, pT2, shift);


            t0 = AE_MULFD16X16X4RAS(x00, x20, w8, w9);
            t1 = AE_MULFD16X16X4RAS(x40, x60, w10, w11);
            y2 = AE_ADD16S(t0, t1);

            t0 = AE_MULFD16X16X4RAS(x10, x30, w8, w9);
            t1 = AE_MULFD16X16X4RAS(x50, x70, w10, w11);
            y6 = AE_ADD16S(t0, t1);

            /* AE_L16X4X2_X is used for better scheduling*/
            AE_L16X4X2_X(x60, x70, pT3, 3 * sizeof(ae_int16x8));
            AE_L16X4X2_X(x40, x50, pT3, 2 * sizeof(ae_int16x8));
            AE_L16X4X2_X(x20, x30, pT3, sizeof(ae_int16x8));
            AE_L16X4X2_XP(x00, x10, pT3, shift);


            t0 = AE_MULFD16X16X4RAS(x00, x20, w12, w13);
            t1 = AE_MULFD16X16X4RAS(x40, x60, w14, w15);
            y3 = AE_ADD16S(t0, t1);

            t0 = AE_MULFD16X16X4RAS(x10, x30, w12, w13);
            t1 = AE_MULFD16X16X4RAS(x50, x70, w14, w15);
            y7 = AE_ADD16S(t0, t1);

            AE_DSEL16X4(x00, x10, y0, y1, sel);
            AE_DSEL16X4(x20, x30, y2, y3, sel);
            AE_DSEL16X4(y0, y2, x00, x20, sel);
            AE_DSEL16X4(y1, y3, x10, x30, sel);

            AE_DSEL16X4(x00, x10, y4, y5, sel);
            AE_DSEL16X4(x20, x30, y6, y7, sel);
            AE_DSEL16X4(y4, y6, x00, x20, sel);
            AE_DSEL16X4(y5, y7, x10, x30, sel);

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
    /* Process last 0...7 rows */
    for (n = (hin & ~7); n < hin; n++)
    {
        /* Move input samples to the scratch */
        pIn = (const ae_int16x8*)((uintptr_t)img + n * stride);
        pOut = (ae_int16x4*)t;
        for (m = 0; m < win; m += 8)
        {
            AE_L16X4X2_IP(x00, x01, pIn, sizeof(ae_int16x8));
            AE_S16X4X2_IP(x00, x01, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
        }

        pIn = (const ae_int16x8*)t; /* load input samples from scratch */
        pOut = (ae_int16x4*)((uintptr_t)img + n * stride);
        pW = (const ae_int16*)w;
        /* Process by 1x4 output samples per iteration */
        for (m = 0; m < wout; m += 4)
        {
            k0 = left[m + 0];
            k1 = left[m + 1];
            k2 = left[m + 2];
            k3 = left[m + 3];
            /* Load by 1x4 input samples for each of 4 output samples */
            pIn = (ae_int16x8*)XT_ADDX2(k0, (uintptr_t)t);
            alIn = AE_LA64_PP(pIn);  AE_LA16X4_IP(x00, alIn, castxcc(ae_int16x4, pIn));
            pIn = (ae_int16x8*)XT_ADDX2(k1, (uintptr_t)t);
            alIn = AE_LA64_PP(pIn);  AE_LA16X4_IP(x10, alIn, castxcc(ae_int16x4, pIn));
            pIn = (ae_int16x8*)XT_ADDX2(k2, (uintptr_t)t);
            alIn = AE_LA64_PP(pIn);  AE_LA16X4_IP(x20, alIn, castxcc(ae_int16x4, pIn));
            pIn = (ae_int16x8*)XT_ADDX2(k3, (uintptr_t)t);
            alIn = AE_LA64_PP(pIn);  AE_LA16X4_IP(x30, alIn, castxcc(ae_int16x4, pIn));
            /* Load window coefficients */

            AE_L16X4X2_IP(w0, w1, castxcc(ae_int16x8, pW), sizeof(ae_int16x8));
            AE_L16X4X2_IP(w2, w3, castxcc(ae_int16x8, pW), sizeof(ae_int16x8));

            AE_MULZAAAA2Q16(A0, A1, x00, x10, w0, w1);
            AE_MULZAAAA2Q16(A2, A3, x20, x30, w2, w3);

            Y0 = AE_TRUNCA32X2F64S(A0, A1, 33);
            Y1 = AE_TRUNCA32X2F64S(A2, A3, 33);
            y0 = AE_ROUND16X4F32SASYM(Y0, Y1);
            AE_S16X4_IP(y0, castxcc(ae_int16x4, pOut), sizeof(ae_int16x4));
        }
    }
}

const imgresizer_api_t imgresizer_api_dnxh={NULL,getCoefSz,getCoef,getScratchSize,process};
