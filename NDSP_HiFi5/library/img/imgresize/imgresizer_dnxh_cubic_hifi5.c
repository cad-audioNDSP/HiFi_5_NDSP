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

/*    image resizer, downsample 1...2x in horizontal direction */

/* returns size of coefficients */
static size_t getCoefSz(const imgsize_t* in,const imgsize_t* out)
{
    int win=in->width,wout=out->width;
    NASSERT(in->height==out->height &&  win<2*wout && win>wout && in->stride==out->stride);
    (void)in,(void)out;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    return img_getCoef_dn_cubic_alloc(win,wout);
}
/* returns coefficients */
static void getCoef(void* coef, const imgsize_t* in,const imgsize_t* out)
{
    int win=in->width,wout=out->width;
    NASSERT(in->height==out->height &&  win<2*wout && win>wout && in->stride==out->stride);
    (void)in,(void)out,(void)coef;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    img_getCoef_dn_cubic_init(coef,win,wout);
}


static size_t getScratchSize(const imgsize_t* in, const imgsize_t* out)
{
    (void)in, (void)out;
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    /* img size */
    return (out->stride) * (out->height) * sizeof(ae_int16);
}

/* in-place image resize */
static void process(void* pScr, void* pCoef, const void* img, void* imgOut, const imgsize_t* in, const imgsize_t* out, int fast)
{
    int16_t* restrict t = (int16_t*)pScr;
    const int16_t* restrict w = ((img_coefdn_cubic_t*)pCoef)->coef;
    const int16_t* restrict left = ((img_coefdn_cubic_t*)pCoef)->left;
    const int16_t* restrict pL;
    const ae_int16x8* restrict pW;
    ae_int16x8* restrict pT0;
    ae_int16x8* restrict pT1;
    ae_int16x8* restrict pT2;
    ae_int16x8* restrict pT3;
    const ae_int16x8* restrict pIn;
    ae_int16x8* restrict pOut;
    ae_int16x4 x00, x01, x10, x11, x20, x21, x30, x31;
    ae_int16x4 w00, w01, w10, w11, w20, w21, w30, w31;
    ae_int16x4 y0;
    ae_int64 A0, A1, A2, A3;
    ae_int32x2 Y0, Y1;
    ae_valignx2 al0;
    int m, n, k0, k1, k2, k3, shift;
    int hout = in->height,
        wout = out->width,
        stride = out->stride;

    NASSERT_ALIGN(img, ALIGNMENT);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->height == out->height && in->width<2 * out->width && in->width>out->width && in->stride == out->stride);
    NASSERT(wout >= 1);
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    wout = (wout + 7) & ~7;
    pW = (const ae_int16x8*)w;
    pL = left;
    stride = stride * sizeof(ae_int16);
    shift = stride - sizeof(ae_int16) * 8;
    pIn = (ae_int16x8*)(img);
    for (m = 0; m < wout; m += 4)
    {

        pOut = (ae_int16x8*)((uintptr_t)t + m * sizeof(ae_int16));

        AE_L16X4X2_IP(w00, w01, pW, sizeof(ae_int16x8));
        AE_L16X4X2_IP(w10, w11, pW, sizeof(ae_int16x8));
        AE_L16X4X2_IP(w20, w21, pW, sizeof(ae_int16x8));
        AE_L16X4X2_IP(w30, w31, pW, sizeof(ae_int16x8));

        k0 = *pL++;  k1 = *pL++;  k2 = *pL++;  k3 = *pL++;
        pT0 = (ae_int16x8*)XT_ADDX2(k0, (uintptr_t)pIn);
        pT1 = (ae_int16x8*)XT_ADDX2(k1, (uintptr_t)pIn);
        pT2 = (ae_int16x8*)XT_ADDX2(k2, (uintptr_t)pIn);
        pT3 = (ae_int16x8*)XT_ADDX2(k3, (uintptr_t)pIn);
        for (n = 0; n < (hout); n++)
        {
            al0 = AE_LA128_PP(pT0);  AE_LA16X4X2_IP(x00, x01, al0, pT0);
            al0 = AE_LA128_PP(pT1);  AE_LA16X4X2_IP(x10, x11, al0, pT1);
            al0 = AE_LA128_PP(pT2);  AE_LA16X4X2_IP(x20, x21, al0, pT2);
            al0 = AE_LA128_PP(pT3);  AE_LA16X4X2_IP(x30, x31, al0, pT3);


            AE_MULZAAAA2Q16(A0, A1, x00, x10, w00, w10);
            AE_MULZAAAA2Q16(A2, A3, x20, x30, w20, w30);
            AE_MULAAAA2Q16(A0, A1, x01, x11, w01, w11);
            AE_MULAAAA2Q16(A2, A3, x21, x31, w21, w31);

            Y0 = AE_TRUNCA32X2F64S(A0, A1, 33);
            Y1 = AE_TRUNCA32X2F64S(A2, A3, 33);

            y0 = AE_ROUND16X4F32SASYM(Y0, Y1);

            AE_S16X4_XP(y0, castxcc(ae_int16x4, pOut), stride);

            pT0 = (ae_int16x8*)((uintptr_t)pT0 + shift);
            pT1 = (ae_int16x8*)((uintptr_t)pT1 + shift);
            pT2 = (ae_int16x8*)((uintptr_t)pT2 + shift);
            pT3 = (ae_int16x8*)((uintptr_t)pT3 + shift);

        }

    }
    /* Move computed samples from the scratch to the result image */
    for (n = 0; n < hout; n++)
    {
        pIn = (ae_int16x8*)((uintptr_t)t + n * stride);
        pOut = (ae_int16x8*)((uintptr_t)img + n * stride);
        __Pragma("loop_count min=1");
        for (m = 0; m < (wout); m += 8)
        {
            AE_L16X4X2_IP(w00, w01, pIn, sizeof(ae_int16x8));
            AE_S16X4X2_IP(w00, w01, pOut, sizeof(ae_int16x8));
        }
    }
}




const imgresizer_api_t imgresizer_api_dnxh_cubic={NULL,getCoefSz,getCoef,getScratchSize,process};
