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
#include "img_getCoef_up_cubic.h"
#include "imgresizer_bicubic_common.h"

/*    image resizer, upsample 1...2x in vertical direction */

/* returns size of coefficients */
static size_t getCoefSz(const imgsize_t* in,const imgsize_t* out)
{
    int hin=in->height,hout=out->height;
    NASSERT(hout<2*hin && hout>hin);
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    return img_getCoef_up_cubic_alloc(hin,hout);
}
/* returns coefficients */
static void getCoef(void* coef, const imgsize_t* in,const imgsize_t* out)
{
    int hin=in->height,hout=out->height;
    NASSERT(hout<2*hin && hout>hin);
    (void)in,(void)out,(void)coef;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    img_getCoef_up_cubic_init(coef,hin,hout);
}

static size_t getScratchSize(const imgsize_t* in, const imgsize_t* out)
{
    (void)in, (void)out;
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    return (out->stride) * (out->height) * sizeof(int16_t);
}

/* in-place image resize */
static void process(void* pScr, void* pCoef, const void* img, void* imgOut, const imgsize_t* in, const imgsize_t* out, int fast)
{
    int16_t* restrict t = (int16_t*)pScr;
    const int16_t* restrict w = ((img_coefup_cubic_t*)pCoef)->coef;
    const int16_t* restrict up = ((img_coefup_cubic_t*)pCoef)->left;
    const ae_int16x8* restrict pIn;
    ae_int16x8* restrict pOut;
    const ae_int16x8* restrict pW;
    ae_int16x4 x00, x01, x02, x03;
    ae_int16x4 x10, x11, x12, x13;
    ae_int16x4 w0, w1, w2, w3;
    ae_int32x2 Y0, Y1, Y2, Y3, t0, t1;
    ae_int16x4 y0, y1;
    int m, n;
    int16_t ix;

    int wout = out->width,
        hin = in->height,
        hout = out->height,
        stride = out->stride;
    (void)hin;
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    NASSERT_ALIGN(img, ALIGNMENT);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->width == out->width && hout<2 * hin && hout>hin && in->stride == out->stride);
    NASSERT(hout > 2);

    stride = stride * sizeof(int16_t);

    pW = (const ae_int16x8*)(w);

    for (n = 0; n < hout; n++)
    {
        pOut = (ae_int16x8*)((uintptr_t)t + stride * n);
        ix = up[n];
        pIn = (const ae_int16x8*)((uintptr_t)img + stride * ix);

        AE_L16_IP(w0, castxcc(ae_int16, pW), sizeof(int16_t));
        AE_L16_IP(w1, castxcc(ae_int16, pW), sizeof(int16_t));
        AE_L16_IP(w2, castxcc(ae_int16, pW), sizeof(int16_t));
        AE_L16_IP(w3, castxcc(ae_int16, pW), sizeof(int16_t));

        for (m = 0; m < (wout); m += 8)
        {
            AE_L16X4X2_XP(x00, x10, pIn, stride);
            AE_L16X4X2_XP(x01, x11, pIn, stride);
            AE_L16X4X2_XP(x02, x12, pIn, stride);
            AE_L16X4X2_XP(x03, x13, pIn, -3 * stride + sizeof(ae_int16x8));

            AE_MULFD16X16X4WS(Y0, Y1, x00, x01, w0, w1);
            AE_MULFD16X16X4WS(t0, t1, x02, x03, w2, w3);
            Y0 = AE_ADD32(Y0, t0);
            Y1 = AE_ADD32(Y1, t1);

            AE_MULFD16X16X4WS(Y2, Y3, x10, x11, w0, w1);
            AE_MULFD16X16X4WS(t0, t1, x12, x13, w2, w3);
            Y2 = AE_ADD32(Y2, t0);
            Y3 = AE_ADD32(Y3, t1);

            y0 = AE_TRUNCA16X4F32S(Y0, Y1, 1);
            y1 = AE_TRUNCA16X4F32S(Y2, Y3, 1);

            AE_S16X4X2_IP(y0, y1, pOut, sizeof(ae_int16x8));
        }
    }
    /* Move computed samples from the scratch to the result image */
    for (n = 0; n < hout; n++)
    {
        pIn = (const ae_int16x8*)((uintptr_t)t + n * stride);
        pOut = (ae_int16x8*)((uintptr_t)img + n * stride);
        __Pragma("loop_count min=1");
        for (m = 0; m < (wout); m += 8)
        {
            AE_L16X4X2_IP(w0, w1, pIn, sizeof(ae_int16x8));
            AE_S16X4X2_IP(w0, w1, pOut, sizeof(ae_int16x8));
        }
    }

}

const imgresizer_api_t imgresizer_api_upxv_cubic={NULL,getCoefSz,getCoef,getScratchSize,process};
