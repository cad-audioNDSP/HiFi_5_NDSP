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
    const ae_int16x8 * restrict pIn;
          ae_int16x8 * restrict pOut;
    ae_int16x4 s00, s01, s10, s11, w0, w1;
    ae_int16x4 y0, y1, y2, y3;
    int m, n;
    int w = in->width,
        hin = in->height,
        hout = out->height,
        stride = out->stride;

    (void)pCoef;
    (void)hin;
    imgsize_validate(in,2,1);
    imgsize_validate(out,2,1);
    NASSERT_ALIGN(img,ALIGNMENT);
    NASSERT_ALIGN(pScr,ALIGNMENT);
    NASSERT(in->width==out->width &&  hout==2*hin && in->stride==out->stride);
    NASSERT(hin >= 2);

    stride = stride * sizeof(int16_t);
    /* Set window coefficients */
    w0 = 24576; /* 0.75, Q15 */
    w1 = 8192;  /* 0.25, Q15 */

    /* Process the image by 8 columns per iteration */
    /* Process samples in reverse order */
    for (m = 0; m < w; m+=8)
    {
        pOut = (      ae_int16x8 *)((uintptr_t)img + (hout - 1)*stride + m*sizeof(int16_t));
        pIn  = (const ae_int16x8 *)((uintptr_t)img + (hin  - 1)*stride + m*sizeof(int16_t));

        AE_L16X4X2_XP(s01, s11, pIn, -stride);
        AE_S16X4X2_XP(s01, s11, pOut, -stride);

        /* Process by 2x8 output samples per iteration */
        __Pragma("loop_count min=1");
        for (n = 0; n < hin-1; n++)
        {
            AE_MOVD16X8(s00, s10, s01, s11);

            AE_L16X4X2_XP(s01, s11, pIn, -stride);
            y0 = AE_MULFD16X16X4RAS(s00, s01, w0, w1);
            y1 = AE_MULFD16X16X4RAS(s00, s01, w1, w0);
            y2 = AE_MULFD16X16X4RAS(s10, s11, w0, w1);
            y3 = AE_MULFD16X16X4RAS(s10, s11, w1, w0);
            AE_S16X4X2_XP(y0, y2, pOut, -stride);
            AE_S16X4X2_XP(y1, y3, pOut, -stride);
        }
    }
}

const imgresizer_api_t imgresizer_api_up2xv={NULL,getCoefSz,getCoef,getScratchSize,process};
