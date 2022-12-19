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
#include "img_getCoef_dn_cubic.h"
#include "imgresizer_bicubic_common.h"

/*    image resizer, downsample 1...2x in vertical direction */

/* returns size of coefficients */
static size_t getCoefSz(const imgsize_t* in,const imgsize_t* out)
{
    int hin=in->height,hout=out->height;
    (void)in,(void)out;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    NASSERT(in->width==out->width &&  hin<2*hout && hin>hout && in->stride==out->stride);
    return img_getCoef_dn_cubic_alloc(hin,hout);
}
/* returns coefficients */
static void getCoef(void* coef, const imgsize_t* in,const imgsize_t* out)
{
    int hin=in->height,hout=out->height;
    (void)in,(void)out;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    NASSERT(in->width==out->width &&  hin<2*hout && hin>hout && in->stride==out->stride);
    img_getCoef_dn_cubic_init(coef,hin,hout);
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
    const int16_t* restrict w = ((img_coefdn_cubic_t*)pCoef)->coef;
    const int16_t* restrict up = ((img_coefdn_cubic_t*)pCoef)->left;
    const ae_int16x8* restrict pIn;
    ae_int16x8* restrict pOut;
    const ae_int16x8* restrict pW;
    ae_int16x4 x00, x01, x02, x03, x04, x05, x06, x07;
    ae_int16x4 x10, x11, x12, x13, x14, x15, x16, x17;
    ae_int16x4 w0, w1, w2, w3, w4, w5, w6, w7;
    ae_int32x2 Y0, Y1, Y2, Y3;
    int m, n;
    int wout = out->width,
        hout = out->height,
        stride = out->stride;
    int16_t ix;

    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    NASSERT_ALIGN(img, ALIGNMENT);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->width == wout && in->height<2 * out->height && in->height>out->height && in->stride == out->stride);
    NASSERT(hout >= 1);

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
        AE_L16_IP(w4, castxcc(ae_int16, pW), sizeof(int16_t));
        AE_L16_IP(w5, castxcc(ae_int16, pW), sizeof(int16_t));
        AE_L16_IP(w6, castxcc(ae_int16, pW), sizeof(int16_t));
        AE_L16_IP(w7, castxcc(ae_int16, pW), sizeof(int16_t));

        for (m = 0; m < (wout); m += 8)
        {
            AE_L16X4X2_XP(x00, x10, pIn, stride);
            AE_L16X4X2_XP(x01, x11, pIn, stride);
            AE_L16X4X2_XP(x02, x12, pIn, stride);
            AE_L16X4X2_XP(x03, x13, pIn, stride);
            AE_L16X4X2_XP(x04, x14, pIn, stride);
            AE_L16X4X2_XP(x05, x15, pIn, stride);
            AE_L16X4X2_XP(x06, x16, pIn, stride);
            AE_L16X4X2_XP(x07, x17, pIn, -7 * stride + sizeof(ae_int16x8));

            AE_MUL16X4(Y0, Y1, x00, w0);
            AE_MULA16X4(Y0, Y1, x01, w1);
            AE_MULA16X4(Y0, Y1, x02, w2);
            AE_MULA16X4(Y0, Y1, x03, w3);
            AE_MULA16X4(Y0, Y1, x04, w4);
            AE_MULA16X4(Y0, Y1, x05, w5);
            AE_MULA16X4(Y0, Y1, x06, w6);
            AE_MULA16X4(Y0, Y1, x07, w7);

            AE_MUL16X4(Y2, Y3, x10, w0);
            AE_MULA16X4(Y2, Y3, x11, w1);
            AE_MULA16X4(Y2, Y3, x12, w2);
            AE_MULA16X4(Y2, Y3, x13, w3);
            AE_MULA16X4(Y2, Y3, x14, w4);
            AE_MULA16X4(Y2, Y3, x15, w5);
            AE_MULA16X4(Y2, Y3, x16, w6);
            AE_MULA16X4(Y2, Y3, x17, w7);

            Y0 = AE_ADD32S(Y0, Y0);
            Y1 = AE_ADD32S(Y1, Y1);
            Y2 = AE_ADD32S(Y2, Y2);
            Y3 = AE_ADD32S(Y3, Y3);

            AE_S16X4RA32S_IP(Y0, Y1, castxcc(ae_int16x4, pOut));
            AE_S16X4RA32S_IP(Y2, Y3, castxcc(ae_int16x4, pOut));
        }
    }
    //return;
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

const imgresizer_api_t imgresizer_api_dnxv_cubic={NULL,getCoefSz,getCoef,getScratchSize,process};
