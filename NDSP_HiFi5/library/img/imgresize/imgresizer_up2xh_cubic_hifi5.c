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

/*    image resizer, upsample 2x in horizontal direction */

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
static void process(void* pScr, void* pCoef, const void* img, void* imgOut, const imgsize_t* in, const imgsize_t* out, int fast)
{

    const ae_int16* restrict pIn3;
    const ae_int16* restrict pIn0;
    ae_int16* restrict pOut;
    int m, n;
    int h = in->height,
        win = in->width,
        wout = out->width,
        stride = out->stride;
    ae_int16x4 x0, x1, x2, x3;
    ae_int16x4 t0, t1, y0, y1;
    ae_int16x4 c0, c1, c2, c3, sel, delay;
    ae_valign alIn0, alIn3, alOut;

    (void)pCoef;
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    NASSERT_ALIGN(img, ALIGNMENT);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->height == out->height && wout == 2 * win && in->stride == out->stride);
    NASSERT(win >= 4);

    /* Set window coefficients */
    c0 = -768;   /* -0.0234375, Q15 */
    c1 = 7424;   /* 0.2265625, Q15 */
    c2 = 28416;  /* 0.8671875, Q15 */
    c3 = -2304;  /* -0.0703125, Q15 */

    sel = AE_MOVINT16X4_FROMINT64(0x0705030106040200); // SEL7362 + SEL5140
    delay = AE_MOVINT16X4_FROMINT64(0x0102070106000506); // SEL1765 + SEL2106

    /* Process the image by 1 row per iteration */
    /* Process samples in reverse order */
    for (n = 0; n < h; n++)
    {
        pIn0 = (const ae_int16*)img + n * stride + win - 1;
        pOut = (ae_int16*)img + n * stride + wout - 1;

        AE_L16_IP(x1, pIn0, -(int)sizeof(int16_t));
        AE_L16_IP(x0, pIn0, -(int)sizeof(int16_t));
        x2 = x3 = x1;

        t0 = AE_MULFD16X16X4RAS(x0, x1, c3, c2);
        t1 = AE_MULFD16X16X4RAS(x2, x3, c1, c0);
        y0 = AE_ADD16S(t0, t1);

        AE_S16_0_IP(y0, pOut, -(int)sizeof(int16_t));

        x3 = x2;
        x2 = x1;
        x1 = x0;
        AE_L16_IP(x0, pIn0, -(int)sizeof(int16_t));

        t0 = AE_MULFD16X16X4RAS(x0, x1, c0, c1);
        t1 = AE_MULFD16X16X4RAS(x2, x3, c2, c3);
        y0 = AE_ADD16S(t0, t1);

        t0 = AE_MULFD16X16X4RAS(x0, x1, c3, c2);
        t1 = AE_MULFD16X16X4RAS(x2, x3, c1, c0);
        y1 = AE_ADD16S(t0, t1);

        AE_S16_0_IP(y0, pOut, -(int)sizeof(int16_t));
        AE_S16_0_IP(y1, pOut, -(int)sizeof(int16_t));

        __Pragma("no_unroll");
        for (m = 0; m < ((win - 3) & 3); m++)
        {
            x3 = x2;
            x2 = x1;
            x1 = x0;
            AE_L16_IP(x0, pIn0, -(int)sizeof(int16_t));

            t0 = AE_MULFD16X16X4RAS(x0, x1, c0, c1);
            t1 = AE_MULFD16X16X4RAS(x2, x3, c2, c3);
            y0 = AE_ADD16S(t0, t1);

            t0 = AE_MULFD16X16X4RAS(x0, x1, c3, c2);
            t1 = AE_MULFD16X16X4RAS(x2, x3, c1, c0);
            y1 = AE_ADD16S(t0, t1);


            AE_S16_0_IP(y0, pOut, -(int)sizeof(int16_t));
            AE_S16_0_IP(y1, pOut, -(int)sizeof(int16_t));
        }

        pIn3 = pIn0 + 3;
        alIn0 = AE_LA64_PP(pIn0);
        alIn3 = AE_LA64_PP(pIn3);
        alOut = AE_ZALIGN64();
        /* Main loop: process by 8 output samples per iteration */
        for (m = 0; m < ((win - 3) / 4); m++)
        {
            AE_LA16X4_RIP(x3, alIn3, castxcc(ae_int16x4, pIn3));
            AE_LA16X4_RIP(x0, alIn0, castxcc(ae_int16x4, pIn0));
            AE_DSEL16X4(x1, x2, x0, x3, delay);

            t0 = AE_MULFD16X16X4RAS(x0, x1, c0, c1);
            t1 = AE_MULFD16X16X4RAS(x2, x3, c2, c3);
            y0 = AE_ADD16S(t0, t1);

            t0 = AE_MULFD16X16X4RAS(x0, x1, c3, c2);
            t1 = AE_MULFD16X16X4RAS(x2, x3, c1, c0);
            y1 = AE_ADD16S(t0, t1);

            AE_DSEL16X4(y0, y1, y0, y1, sel);

            AE_SA16X4_RIP(y0, alOut, castxcc(ae_int16x4, pOut));
            AE_SA16X4_RIP(y1, alOut, castxcc(ae_int16x4, pOut));
        }
        AE_SA64NEG_FP(alOut, pOut);
        x3 = x2;
        x2 = x1;
        x1 = x0;

        t0 = AE_MULFD16X16X4RAS(x0, x1, c0, c1);
        t1 = AE_MULFD16X16X4RAS(x2, x3, c2, c3);
        y0 = AE_ADD16S(t0, t1);

        t0 = AE_MULFD16X16X4RAS(x0, x1, c3, c2);
        t1 = AE_MULFD16X16X4RAS(x2, x3, c1, c0);
        y1 = AE_ADD16S(t0, t1);

        AE_S16_0_IP(y0, pOut, -(int)sizeof(int16_t));
        AE_S16_0_IP(y1, pOut, -(int)sizeof(int16_t));

        x3 = x2;
        x2 = x1;

        t0 = AE_MULFD16X16X4RAS(x0, x1, c0, c1);
        t1 = AE_MULFD16X16X4RAS(x2, x3, c2, c3);
        y0 = AE_ADD16S(t0, t1);

        AE_S16_0_IP(y0, pOut, -(int)sizeof(int16_t));

    }
}

const imgresizer_api_t imgresizer_api_up2xh_cubic={NULL,getCoefSz,getCoef,getScratchSize,process};
