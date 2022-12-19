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

/*    image resizer, downsample 2x in vertical direction */

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
    const ae_int16x8 *restrict pIn;
          ae_int16x8 *restrict pOut;
    ae_int16x4 x00, x10, x20, x30, x40, x50, x61, x71;
    ae_int16x4 x01, x11, x21, x31, x41, x51, x60, x70;
    ae_int16x4 w0, w1, w2, w3;
    ae_int16x4 Y0, Y1, Y2, Y3;
    ae_int16x4 y0, y1;
    int m, n;
    int w = in->width,
        hin = in->height,
        hout = out->height,
        stride = out->stride;
    (void)hin;
    (void)pCoef;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    NASSERT_ALIGN(img,ALIGNMENT);
    NASSERT_ALIGN(pScr,ALIGNMENT);
    NASSERT(w==out->width &&  hin==2*hout && stride==out->stride);
    NASSERT(hout >= 2);

    /* Set window coefficients */
    w0 = -384; /* -0.0117188 in Q15 */
    w1 = -1152; /* -0.0351563 in Q15 */
    w2 = 3712; /* 0.1132813 in Q15 */
    w3 = 14208; /* 0.4335938 in Q15 */

    stride = stride * sizeof(int16_t);

    pOut = (      ae_int16x8 *)img;
    pIn  = (const ae_int16x8 *)img;

    /* Process image by 8 columns per iteration */
    for (m=0; m<w; m+=8)
    {
        pIn = (const ae_int16x8*)((int16_t*)img+m);
        pOut = (     ae_int16x8*)((int16_t*)img+m);

        AE_L16X4X2_XP(x30, x31, pIn, stride);
        x00 = x30;
        x10 = x30;
        x20 = x30;

        x01 = x31;
        x11 = x31;
        x21 = x31;

        AE_L16X4X2_X(x40, x41, pIn, 0);
        AE_L16X4X2_X(x50, x51, pIn, stride);
        AE_L16X4X2_X(x60, x61, pIn, 2*stride);
        AE_L16X4X2_X(x70, x71, pIn, 3*stride);

        /* Process by 1x8 samples per iteration */
        __Pragma("loop_count min=1");
        for (n = 0; n < hout - 3; n++)
        {
            Y0 = AE_MULFD16X16X4RAS(x00, x10, w0, w1);
            Y1 = AE_MULFD16X16X4RAS(x20, x30, w2, w3);
            Y2 = AE_MULFD16X16X4RAS(x40, x50, w3, w2);
            Y3 = AE_MULFD16X16X4RAS(x60, x70, w1, w0);
            Y0 = AE_ADD16S(Y0, Y1);
            Y2 = AE_ADD16S(Y2, Y3);
            y0 = AE_ADD16S(Y0, Y2);

            Y0 = AE_MULFD16X16X4RAS(x01, x11, w0, w1);
            Y1 = AE_MULFD16X16X4RAS(x21, x31, w2, w3);
            Y2 = AE_MULFD16X16X4RAS(x41, x51, w3, w2);
            Y3 = AE_MULFD16X16X4RAS(x61, x71, w1, w0);
            Y0 = AE_ADD16S(Y0, Y1);
            Y2 = AE_ADD16S(Y2, Y3);
            y1 = AE_ADD16S(Y0, Y2);

            AE_S16X4X2_XP(y0, y1, pOut, stride);

            AE_MOVD16X8(x00, x01, x20, x21);
            AE_MOVD16X8(x10, x11, x30, x31);
            //x00 = x20;  x01 = x21;
            //x10 = x30;  x11 = x31;

            AE_L16X4X2_XP(x20, x21, pIn, stride);
            AE_L16X4X2_XP(x30, x31, pIn, stride);
            AE_L16X4X2_X(x40, x41, pIn, 0);
            AE_L16X4X2_X(x50, x51, pIn, stride);
            AE_L16X4X2_X(x60, x61, pIn, 2 * stride);
            AE_L16X4X2_X(x70, x71, pIn, 3 * stride);

        }
        /* Process last 1x8 samples */
        Y0 = AE_MULFD16X16X4RAS(x00, x10, w0, w1);
        Y1 = AE_MULFD16X16X4RAS(x20, x30, w2, w3);
        Y2 = AE_MULFD16X16X4RAS(x40, x50, w3, w2);
        Y3 = AE_MULFD16X16X4RAS(x60, x70, w1, w0);
        Y0 = AE_ADD16S(Y0, Y1);
        Y2 = AE_ADD16S(Y2, Y3);
        y0 = AE_ADD16S(Y0, Y2);

        Y0 = AE_MULFD16X16X4RAS(x01, x11, w0, w1);
        Y1 = AE_MULFD16X16X4RAS(x21, x31, w2, w3);
        Y2 = AE_MULFD16X16X4RAS(x41, x51, w3, w2);
        Y3 = AE_MULFD16X16X4RAS(x61, x71, w1, w0);
        Y0 = AE_ADD16S(Y0, Y1);
        Y2 = AE_ADD16S(Y2, Y3);
        y1 = AE_ADD16S(Y0, Y2);

        AE_S16X4X2_XP(y0, y1, pOut, stride);

        x00 = x20;  x01 = x21;
        x10 = x30;  x11 = x31;
        x20 = x40;  x21 = x41;
        x30 = x50;  x31 = x51;
        x40 = x60;  x41 = x61;
        x50 = x70;  x51 = x71;

        AE_L16X4X2_X(x60, x61, pIn, 4*stride);

        x70 = x60; x71 = x61;

        Y0 = AE_MULFD16X16X4RAS(x00, x10, w0, w1);
        Y1 = AE_MULFD16X16X4RAS(x20, x30, w2, w3);
        Y2 = AE_MULFD16X16X4RAS(x40, x50, w3, w2);
        Y3 = AE_MULFD16X16X4RAS(x60, x70, w1, w0);
        Y0 = AE_ADD16S(Y0, Y1);
        Y2 = AE_ADD16S(Y2, Y3);
        y0 = AE_ADD16S(Y0, Y2);

        Y0 = AE_MULFD16X16X4RAS(x01, x11, w0, w1);
        Y1 = AE_MULFD16X16X4RAS(x21, x31, w2, w3);
        Y2 = AE_MULFD16X16X4RAS(x41, x51, w3, w2);
        Y3 = AE_MULFD16X16X4RAS(x61, x71, w1, w0);
        Y0 = AE_ADD16S(Y0, Y1);
        Y2 = AE_ADD16S(Y2, Y3);
        y1 = AE_ADD16S(Y0, Y2);

        AE_S16X4X2_XP(y0, y1, pOut, stride);

        x00 = x20;  x01 = x21;
        x10 = x30;  x11 = x31;
        x20 = x40;  x21 = x41;
        x30 = x50;  x31 = x51;
        x40 = x60;  x41 = x61;
        x50 = x70;  x51 = x71;

        Y0 = AE_MULFD16X16X4RAS(x00, x10, w0, w1);
        Y1 = AE_MULFD16X16X4RAS(x20, x30, w2, w3);
        Y2 = AE_MULFD16X16X4RAS(x40, x50, w3, w2);
        Y3 = AE_MULFD16X16X4RAS(x60, x70, w1, w0);
        Y0 = AE_ADD16S(Y0, Y1);
        Y2 = AE_ADD16S(Y2, Y3);
        y0 = AE_ADD16S(Y0, Y2);

        Y0 = AE_MULFD16X16X4RAS(x01, x11, w0, w1);
        Y1 = AE_MULFD16X16X4RAS(x21, x31, w2, w3);
        Y2 = AE_MULFD16X16X4RAS(x41, x51, w3, w2);
        Y3 = AE_MULFD16X16X4RAS(x61, x71, w1, w0);
        Y0 = AE_ADD16S(Y0, Y1);
        Y2 = AE_ADD16S(Y2, Y3);
        y1 = AE_ADD16S(Y0, Y2);

        AE_S16X4X2_XP(y0, y1, pOut, stride);

    }
}

const imgresizer_api_t imgresizer_api_dn2xv_cubic={NULL,getCoefSz,getCoef,getScratchSize,process};
