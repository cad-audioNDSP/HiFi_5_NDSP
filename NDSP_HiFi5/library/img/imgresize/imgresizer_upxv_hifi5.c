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

/*    image resizer, upsample 1...2x in vertical direction */

/* returns size of coefficients */
static size_t getCoefSz(const imgsize_t* in,const imgsize_t* out)
{
    int hin=in->height,hout=out->height;
    NASSERT(hout<2*hin && hout>hin);
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    return img_getCoef_up_alloc(hin,hout);
}
/* returns coefficients */
static void getCoef(void* coef, const imgsize_t* in,const imgsize_t* out)
{
    int hin=in->height,hout=out->height;
    NASSERT(hout<2*hin && hout>hin);
    (void)in,(void)out,(void)coef;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    img_getCoef_up_init(coef,hin,hout);
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
    const int16_t* restrict w = ((img_coefup_t*)pCoef)->coef;
    const int16_t* restrict up = ((img_coefup_t*)pCoef)->left;
    const ae_int16x8* restrict pIn;
    ae_int16x8* restrict pOut;
    const ae_int16x8* restrict pW;
    ae_int16x4 x00, x01;
    ae_int16x4 x10, x11;
    ae_int16x4 w0, w1;
    ae_int16x4 y0, y1;
    int m, n;
    int wout = out->width,
        hout = out->height,
        hin = in->height,
        stride = out->stride;
    int16_t ix;
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
        for (m = 0; m < (wout); m += 8)
        {
            AE_L16X4X2_XP(x00, x01, pIn, stride);
            AE_L16X4X2_XP(x10, x11, pIn, -stride + sizeof(ae_int16x8));


            y0 = AE_MULFD16X16X4RAS(x00, x10, w0, w1);
            y1 = AE_MULFD16X16X4RAS(x01, x11, w0, w1);

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

const imgresizer_api_t imgresizer_api_upxv={"imgresizer_api_upxv",getCoefSz,getCoef,getScratchSize,process};
