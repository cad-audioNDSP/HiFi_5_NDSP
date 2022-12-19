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
#include "imgresizer_bilinear_common.h"

/*    image resizer, downsample 2x in horizontal direction */

/* returns size of coefficients */
static size_t getCoefSz(const imgsize_t* in,const imgsize_t* out)
{
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    return 0;
}
/* returns coefficients */
static void getCoef(void* coef, const imgsize_t* in,const imgsize_t* out)
{
    (void)coef;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
}


#if defined(AE_MULQQ8X16CNV) 

static size_t getScratchSize(const imgsize_t* in, const imgsize_t* out)
{
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    return 0;
}

/* in-place image resize */
static void process(void* pScr, void* pCoef, const void* img, void* imgOut, const imgsize_t* in, const imgsize_t* out, int fast)
{
    const int8_t ALIGN(ALIGNMENT) c[8] = { 1,3,3,1,1,3,3,1 };
    ae_int8x8 c8bit = AE_L8X8_I((const ae_int8x8*)c, 0);
    ae_int16x4 S0, S1, S2, S3, S4, Z0, Z1;
    ae_int32x2 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
    ae_valignx2 aX;
    ae_valign aX_64;

    int16_t* restrict pIn;
    int16_t* restrict pIn1;
    int16_t* restrict pOut;
    int16_t* restrict pOut1;
    int m, n, h = in->height, win = in->width, wout = out->width, stride = out->stride;
    (void)in, (void)out, (void)pCoef;
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    NASSERT_ALIGN(img, ALIGNMENT);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->height == out->height && win == 2 * wout && in->stride == out->stride);

    pIn = (int16_t*)img;
    pIn1 = (int16_t*)img + win - 1;
    pOut = (int16_t*)img - 1;
    pOut1 = (int16_t*)img + win;
    for (n = 0; n < h; n++)
    {
        AE_L16_XP(S0, castxcc(ae_int16, pIn), stride * sizeof(int16_t));
        AE_S16_0_XP(S0, castxcc(ae_int16, pOut), stride * sizeof(int16_t));
        AE_L16_XP(S1, castxcc(ae_int16, pIn1), stride * sizeof(int16_t));
        AE_S16_0_XP(S1, castxcc(ae_int16, pOut1), stride * sizeof(int16_t));
    }

    for (n = 0; n < h; n++)
    {
        pIn = ((int16_t*)img) + n * stride - 1;
        pOut = ((int16_t*)img) + n * stride;
        aX_64 = AE_LA64_PP((ae_int16x4*)pIn);
        AE_LA16X4_IP(S0, aX_64, castxcc(ae_int16x4, pIn));
        aX = AE_LA128_PP((ae_int16x8*)pIn);
        // 15 unroll 4
        for (m = 0; m < wout; m += 8)
        {
            AE_LA16X4X2_IP(S1, S2, aX, castxcc(ae_int16x8, pIn));
            AE_LA16X4X2_IP(S3, S4, aX, castxcc(ae_int16x8, pIn));
            AE_MULQQ8X16CNV(Y0, Y1, Y2, Y3, c8bit, S0, S1, S1, S2);
            AE_MULQQ8X16CNV(Y4, Y5, Y6, Y7, c8bit, S2, S3, S3, S4);
            S0 = S4;
            Y0 = AE_SEL32_HH(Y0, Y1);
            Y2 = AE_SEL32_HH(Y2, Y3);
            Y4 = AE_SEL32_HH(Y4, Y5);
            Y6 = AE_SEL32_HH(Y6, Y7);
            Z0 = AE_TRUNCI16X4F32S(Y0, Y2, 13);
            Z1 = AE_TRUNCI16X4F32S(Y4, Y6, 13);
            AE_S16X4X2_IP(Z0, Z1, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
        }
    }
}
#else

static size_t getScratchSize(const imgsize_t* in, const imgsize_t* out)
{
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    //return out->width*sizeof(int16_t);// one row
    return 0;
}

/* in-place image resize */
static void process(void* pScr, void* pCoef, const void* img, void* imgOut, const imgsize_t* in, const imgsize_t* out, int fast)
{
    const int16_t* restrict pIn;
    const int16_t* restrict pIn1;
    int16_t* restrict pOut;
    int16_t* restrict pOut1;
    ae_int16x4 x0, x1, x2, x3, t0, t1, t2, t3, y0, y1;
    ae_int16x4  sel;
    ae_int16x4 Y0, Y1;
    ae_int16x4 _4096, _12288;
    ae_valignx2 al_in0, al_in1;

    int m, n;
    int h = in->height,
        win = in->width,
        wout = out->width,
        stride = out->stride;
    (void)win;
    (void)pCoef;
    NASSERT_ALIGN(img, ALIGNMENT);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->height == out->height && win == 2 * wout && in->stride == out->stride);
    NASSERT(wout >= 2);
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);

    _4096 = 4096;
    _12288 = 12288;

    sel = AE_MOVINT16X4_FROMINT64(0x0706050403020100); // sel 7531+6420

    pIn = (const int16_t*)img;
    pIn1 = (const int16_t*)img+win-1;
    pOut = (int16_t*)img-1;
    pOut1 = (int16_t*)img+win;
    for (n = 0; n < h; n++)
    {
        AE_L16_XP(x0, castxcc(ae_int16, pIn), stride * sizeof(int16_t));
        AE_S16_0_XP(x0, castxcc(ae_int16, pOut), stride * sizeof(int16_t));
        AE_L16_XP(x1, castxcc(ae_int16, pIn1), stride * sizeof(int16_t));
        AE_S16_0_XP(x1, castxcc(ae_int16, pOut1), stride * sizeof(int16_t));
    }


    for (n = 0; n < h; n++)
    {
        pIn = ((int16_t*)img + n * stride-1);
        pOut = ((int16_t*)img + n * stride);

        al_in0 = AE_LA128_PP(pIn);
        pIn1 = (int16_t*)pIn + 2;
        al_in1 = AE_LA128_PP(pIn1);

        __Pragma("loop_count min=1");
        for (m = 0; m < (wout); m+=8)
        {
            AE_LA16X4X2_IP(t0, t2, al_in0, castxcc(ae_int16x8, pIn));
            AE_LA16X4X2_IP(t1, t3, al_in1, castxcc(ae_int16x8, pIn1));

            AE_DSEL16X4(x0, x1, t0, t2, sel);
            AE_DSEL16X4(x2, x3, t1, t3, sel);

            Y0 = AE_MULFD16X16X4RAS(x0, x1, _4096, _12288);
            Y1 = AE_MULFD16X16X4RAS(x2, x3, _12288, _4096);
            y0 = AE_ADD16S(Y0, Y1);

            AE_LA16X4X2_IP(t0, t2, al_in0, castxcc(ae_int16x8, pIn));
            AE_LA16X4X2_IP(t1, t3, al_in1, castxcc(ae_int16x8, pIn1));

            AE_DSEL16X4(x0, x1, t0, t2, sel);
            AE_DSEL16X4(x2, x3, t1, t3, sel);

            Y0 = AE_MULFD16X16X4RAS(x0, x1, _4096, _12288);
            Y1 = AE_MULFD16X16X4RAS(x2, x3, _12288, _4096);
            y1 = AE_ADD16S(Y0, Y1);

            AE_S16X4X2_IP(y0, y1, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));

        }
    }
}

#endif

const imgresizer_api_t imgresizer_api_dn2xh={NULL,getCoefSz,getCoef,getScratchSize,process};
