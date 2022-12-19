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
#include "imgresizer_bicubic_common.h"

/*    image resizer, upsample 2x in vertical direction */

/* returns size of coefficients */
static size_t getCoefSz(const imgsize_t* in,const imgsize_t* out)
{
    (void)in,(void)out;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    return 0;
}
/* returns coefficients */
static void getCoef(void* coef, const imgsize_t* in,const imgsize_t* out)
{
    (void)in,(void)out,(void)coef;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
}

static size_t getScratchSize(const imgsize_t* in,const imgsize_t* out)
{
    (void)in,(void)out;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    return 0;
}

/* in-place image resize */
static void process(void *pScr, void* pCoef, const void* img, void*imgOut, const imgsize_t* in, const imgsize_t* out, int fast)
{
    int m, n;
    const ae_int16x8   * restrict pIn;
          ae_int16x8 * restrict pOut;
    int w = in->width,
        hin = in->height,
        hout = out->height,
        stride = out->stride;
    ae_int16x4 s00, s01;
    ae_int16x4 s10, s11;
    ae_int16x4 s20, s21;
    ae_int16x4 s30, s31;
    ae_int16x4 w0, w1, w2, w3;
    ae_int16x4 Y0, Y1;
    ae_int16x4 y0, y1, y2, y3;

    (void)pCoef;
    (void)hin;
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    NASSERT_ALIGN(img, ALIGNMENT);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->width == out->width &&  hout == 2 * hin && in->stride == out->stride);

    /* Set window coefficients */
    w0 = -768; /* -0.0234375, Q15 */
    w1 = 7423; /* 0.2265625, Q15 */
    w2 = 28416; /* 0.8671875, Q15 */
    w3 = -2304; /* -0.0703125, Q15 */
    stride = stride*sizeof(int16_t);

    /* Process the image by 8 columns per iteration */
    /* Process samples in reverse order */
    for (m = 0; m < w; m += 8)
    {
        pOut = (ae_int16x8*)((uintptr_t)img + (hout - 1) * stride + m * sizeof(int16_t));
        pIn = (const ae_int16x8*)((uintptr_t)img + (hin - 1) * stride + m * sizeof(int16_t));

        AE_L16X4X2_XP(s30, s31, pIn, -stride);
        s10 = s30;
        s20 = s30;
        s11 = s31;
        s21 = s31;
        AE_L16X4X2_XP(s00, s01, pIn, -stride);

        Y0 = AE_MULFD16X16X4RAS(s00, s10, w3, w2);
        Y1 = AE_MULFD16X16X4RAS(s20, s30, w1, w0);
        y0 = AE_ADD16S(Y0,Y1);

        Y0 = AE_MULFD16X16X4RAS(s01, s11, w3, w2);
        Y1 = AE_MULFD16X16X4RAS(s21, s31, w1, w0);
        y1 = AE_ADD16S(Y0, Y1);

        AE_S16X4X2_XP(y0, y1, pOut, -stride);

        AE_MOVD16X8(s30, s31, s20, s21);
        AE_MOVD16X8(s20, s21, s10, s11);
        //s30 = s20; s31 = s21;
        //s20 = s10; s21 = s11;
        //s10 = s00; s11 = s01;
        AE_L16X4X2_X(s10, s11, pIn, stride);
        AE_L16X4X2_XP(s00, s01, pIn, -stride);

        Y0 = AE_MULFD16X16X4RAS(s00, s10, w0, w1);
        Y1 = AE_MULFD16X16X4RAS(s20, s30, w2, w3);
        y0 = AE_ADD16S(Y0, Y1);

        Y0 = AE_MULFD16X16X4RAS(s01, s11, w0, w1);
        Y1 = AE_MULFD16X16X4RAS(s21, s31, w2, w3);
        y1 = AE_ADD16S(Y0, Y1);

        Y0 = AE_MULFD16X16X4RAS(s00, s10, w3, w2);
        Y1 = AE_MULFD16X16X4RAS(s20, s30, w1, w0);
        y2 = AE_ADD16S(Y0, Y1);

        Y0 = AE_MULFD16X16X4RAS(s01, s11, w3, w2);
        Y1 = AE_MULFD16X16X4RAS(s21, s31, w1, w0);
        y3 = AE_ADD16S(Y0, Y1);

        AE_S16X4X2_XP(y0, y1, pOut, -stride);
        AE_S16X4X2_XP(y2, y3, pOut, -stride);


        /* Process by 2x8 output samples per iteration */
        __Pragma("loop_count min=1");
        for (n = 0; n < hin - 3; n++)
        {
            AE_L16X4X2_X(s30, s31, pIn, 3*stride);
            AE_L16X4X2_X(s20, s21, pIn, 2*stride);
            AE_L16X4X2_X (s10, s11, pIn, stride);
            AE_L16X4X2_XP(s00, s01, pIn, -stride);

            Y0 = AE_MULFD16X16X4RAS(s00, s10, w0, w1);
            Y1 = AE_MULFD16X16X4RAS(s20, s30, w2, w3);
            y0 = AE_ADD16S(Y0, Y1);

            Y0 = AE_MULFD16X16X4RAS(s01, s11, w0, w1);
            Y1 = AE_MULFD16X16X4RAS(s21, s31, w2, w3);
            y1 = AE_ADD16S(Y0, Y1);

            Y0 = AE_MULFD16X16X4RAS(s00, s10, w3, w2);
            Y1 = AE_MULFD16X16X4RAS(s20, s30, w1, w0);
            y2 = AE_ADD16S(Y0, Y1);

            Y0 = AE_MULFD16X16X4RAS(s01, s11, w3, w2);
            Y1 = AE_MULFD16X16X4RAS(s21, s31, w1, w0);
            y3 = AE_ADD16S(Y0, Y1);

            AE_S16X4X2_XP(y0, y1, pOut, -stride);
            AE_S16X4X2_XP(y2, y3, pOut, -stride);
        }
        AE_L16X4X2_X(s30, s31, pIn, 3 * stride);
        AE_L16X4X2_X(s20, s21, pIn, 2 * stride);
        AE_L16X4X2_X (s10, s11, pIn, stride);
        //s10 = s00; s11 = s01;

        Y0 = AE_MULFD16X16X4RAS(s00, s10, w0, w1);
        Y1 = AE_MULFD16X16X4RAS(s20, s30, w2, w3);
        y0 = AE_ADD16S(Y0, Y1);

        Y0 = AE_MULFD16X16X4RAS(s01, s11, w0, w1);
        Y1 = AE_MULFD16X16X4RAS(s21, s31, w2, w3);
        y1 = AE_ADD16S(Y0, Y1);

        Y0 = AE_MULFD16X16X4RAS(s00, s10, w3, w2);
        Y1 = AE_MULFD16X16X4RAS(s20, s30, w1, w0);
        y2 = AE_ADD16S(Y0, Y1);

        Y0 = AE_MULFD16X16X4RAS(s01, s11, w3, w2);
        Y1 = AE_MULFD16X16X4RAS(s21, s31, w1, w0);
        y3 = AE_ADD16S(Y0, Y1);

        AE_S16X4X2_XP(y0, y1, pOut, -stride);
        AE_S16X4X2_XP(y2, y3, pOut, -stride);

        s30 = s20; s31 = s21;
        s20 = s10; s21 = s11;
        
        Y0 = AE_MULFD16X16X4RAS(s00, s10, w0, w1);
        Y1 = AE_MULFD16X16X4RAS(s20, s30, w2, w3);
        y0 = AE_ADD16S(Y0, Y1);

        Y0 = AE_MULFD16X16X4RAS(s01, s11, w0, w1);
        Y1 = AE_MULFD16X16X4RAS(s21, s31, w2, w3);
        y1 = AE_ADD16S(Y0, Y1);

        AE_S16X4X2_XP(y0, y1, pOut, -stride);
    }
}

const imgresizer_api_t imgresizer_api_up2xv_cubic={NULL,getCoefSz,getCoef,getScratchSize,process};
