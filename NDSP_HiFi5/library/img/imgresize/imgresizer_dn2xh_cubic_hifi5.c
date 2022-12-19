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

/*    image resizer, downsample 2x in horizontal direction */

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


static size_t getScratchSize(const imgsize_t* in, const imgsize_t* out)
{
    (void)in, (void)out;
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    return 0;
}

/* in-place image resize */
static void process(void* pScr, void* pCoef, const void* img, void* imgOut, const imgsize_t* in, const imgsize_t* out, int fast)
{
    static const int16_t ALIGN(ALIGNMENT) coefTbl[] =
    {
        /* Half of the window coefficients; the second half is symmetric to the first one */
        /* [-0.0117188, -0.0351563, 0.1132813, 0.4335938] in Q15 */
        -384, -1152, 3712, 14208
    };
    const ae_int16x8* restrict pIn;
    const ae_int16x8* restrict pT0;
    const ae_int16x8* restrict pT1;
    const ae_int16x8* restrict pT2;
    const ae_int16x8* restrict pT3;
    const ae_int16* restrict pIn1;
    ae_int16x4* restrict pOut;
    ae_int16x4* restrict pOut1;
    ae_int16x4 x00, x01, x10, x11, x20, x21, x30, x31;
    ae_int16x4 sel;
    ae_int16x4 cf0, cf1;
    ae_int64 A0, A1, A2, A3;
    ae_int32x2 Y0, Y1;
    ae_valignx2 al_in0;
    ae_valignx2 al_in1;
    ae_valignx2 al_in2;
    ae_valignx2 al_in3;
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
    NASSERT(wout > 2);
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);

    cf0 = AE_L16X4_I((const ae_int16x4*)coefTbl, 0);
    cf1 = AE_MOVINT16X4_FROMINT64(AE_L64_I((const ae_int64*)coefTbl, 0));

    int last = (win - 1) & 3;
    sel = AE_MOVINT16X4_FROMINT16(3 - last);
    xtbool4  choose;
    if (last == 3)
        choose = AE_MOVB4(0);
    else if (last == 2)
        choose = AE_MOVB4(1);
    else if (last == 1)
        choose = AE_MOVB4(3);
    else if (last == 0)
        choose = AE_MOVB4(7);
    pIn = (ae_int16x8*)((int16_t*)img);
    pIn1 = (ae_int16*)((int16_t*)img + ((win - 1) & ~3));
    pOut = (ae_int16x4*)((int16_t*)img - 4);
    pOut1 = (ae_int16x4*)((int16_t*)img + ((win - 1) & ~3));
    for (n = 0; n < h; n++)
    {
        AE_L16_XP(x00, castxcc(ae_int16, pIn), stride * sizeof(int16_t));
        AE_S16X4_XP(x00, castxcc(ae_int16x4, pOut), stride * sizeof(int16_t));
        AE_L16X4_XP(x10, castxcc(ae_int16x4, pIn1), stride * sizeof(int16_t));
        x11 = AE_SEL16X4(x10, x10, sel);
        AE_MOVT16X4(x10, x11, choose);
        AE_S16X4_X(x11, castxcc(ae_int16x4, pOut1), sizeof(ae_int16x4));
        AE_S16X4_XP(x10, castxcc(ae_int16x4, pOut1), stride * sizeof(int16_t));
    }

    sel = AE_MOVINT16X4_FROMINT64(0x0604050304020301); // SEL6543 + SEL4321

    /* process image by 1 row per iteration */
    for (n = 0; n < h; n++)
    {
        pOut = (ae_int16x4*)((int16_t*)img + n * stride);

        pT0 = (ae_int16x8*)((int16_t*)img + n * stride - 3);
        pT1 = (ae_int16x8*)((int16_t*)pT0 + 2);
        pT2 = (ae_int16x8*)((int16_t*)pT0 + 4);
        pT3 = (ae_int16x8*)((int16_t*)pT0 + 6);

        al_in0 = AE_LA128_PP(pT0);
        al_in1 = AE_LA128_PP(pT1);
        al_in2 = AE_LA128_PP(pT2);
        al_in3 = AE_LA128_PP(pT3);

        for (m = 0; m < wout; m += 4)
        {
            AE_LA16X4X2_IP(x00, x01, al_in0, pT0);
            AE_LA16X4X2_IP(x10, x11, al_in1, pT1);
            AE_LA16X4X2_IP(x20, x21, al_in2, pT2);
            AE_LA16X4X2_IP(x30, x31, al_in3, pT3);

            AE_MULZAAAA2Q16(A0, A1, x00, x10, cf0, cf0);
            AE_MULZAAAA2Q16(A2, A3, x20, x30, cf0, cf0);
            AE_MULAAAA2Q16(A0, A1, x01, x11, cf1, cf1);
            AE_MULAAAA2Q16(A2, A3, x21, x31, cf1, cf1);

            Y0 = AE_TRUNCA32X2F64S(A0, A1, 33);
            Y1 = AE_TRUNCA32X2F64S(A2, A3, 33);
            AE_S16X4RA32S_IP(Y0, Y1, pOut);

        }
    }
}


const imgresizer_api_t imgresizer_api_dn2xh_cubic={NULL,getCoefSz,getCoef,getScratchSize,process};

