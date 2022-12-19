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

/*-------------------------------------------------------------------------
  Image normalization
  Function normalize the intensity of pixels to the given range

  Image formats:
  gu8    8-bit unsigned grayscale data
  gs8    8-bit signed grayscale data
  gs16   16-bit grayscale data

  Input:
  inImg   input image
  sz      image size
  minInt  min intensity on output (for linear normalization)
  maxInt  max intensity on output (for non-linear normalization)
  tbl[64] tabulated values (for non-linear normalization)
  Input:
  outImg   input image

  Restrictions:
  see general restrictions applied for all images for fast/generic 
  functions
-------------------------------------------------------------------------*/
void imgfastnorm_gu8_nonlinear ( void * restrict outImg, const void * restrict inImg, const imgsize_t* sz, const int16_t* tbl)
{
    const uint8_t * restrict in = (const uint8_t *)inImg;
    uint8_t * restrict out = (uint8_t *)outImg;
    const signed char * restrict pIn;
    ae_int8x8 * restrict pOut;
    int i, j;
    ae_int8x8 t0;
    ae_int16x4 y00, y01, y02, y03;
    ae_int16x4 y10, y11, y12, y13;
    ae_int16x4 x, ix, x0, z, y0, y1;
    ae_int32x2 d0, d1;
    ae_int16x4  y00_, y10_, y01_, y11_;
    ae_int16x4  y02_, y12_, y03_, y13_;
    ae_int16x4 x_, ix_, x0_, z_, y0_, y1_;
    ae_int32x2 d0_, d1_;
    int h = (int)sz->height;
    int w = (int)sz->width;
    int istride = sz->stride;
    NASSERT(inImg != NULL);
    NASSERT(outImg != NULL);
    NASSERT_ALIGN(inImg, ALIGNMENT);
    NASSERT_ALIGN(outImg, ALIGNMENT);
    imgsize_validate(sz, 1, 1);

    for (i = 0; i < h; i++)
    {
        pIn = (const signed char*)(in + i*istride);
        pOut = (ae_int8x8*)(out + i*istride);
        for (j = 0; j < (w >> 3); j++)
        {  
            AE_L8X4U_IP(x, pIn, 4 * sizeof(int8_t));
            ix = AE_SRAI16(AE_SUB16(x, 2), 2);
            AE_MINMAX16(ix, 0, 62);
            x0 = AE_ADD16(AE_SLAI16S(ix, 2), 2);

            ix = AE_ADD16(ix, ix);
            y00 = AE_L16_X((const ae_int16 *)tbl, AE_MOVAD16_3(ix));
            y10 = AE_L16_X((const ae_int16 *)(tbl + 1), AE_MOVAD16_3(ix));
            y01 = AE_L16_X((const ae_int16 *)tbl, AE_MOVAD16_2(ix));
            y11 = AE_L16_X((const ae_int16 *)(tbl + 1), AE_MOVAD16_2(ix));
            y02 = AE_L16_X((const ae_int16 *)tbl, AE_MOVAD16_1(ix));
            y12 = AE_L16_X((const ae_int16 *)(tbl + 1), AE_MOVAD16_1(ix));
            y03 = AE_L16_X((const ae_int16 *)tbl, AE_MOVAD16_0(ix));
            y13 = AE_L16_X((const ae_int16 *)(tbl + 1), AE_MOVAD16_0(ix));

            AE_MOVT16X4(y00, y01, AE_MOVAB4(4));
            AE_MOVT16X4(y02, y03, AE_MOVAB4(4));
            y0 = AE_SEL16_7632(y00, y02);

            AE_MOVT16X4(y10, y11, AE_MOVAB4(4));
            AE_MOVT16X4(y12, y13, AE_MOVAB4(4));
            y1 = AE_SEL16_7632(y10, y12);

            AE_MULF16X4SS(d0, d1, AE_SLAI16S(AE_SUB16(y1, y0), 1), AE_SLAI16S(AE_SUB16(x, x0), 12));
            z = AE_ADD16S(y0, AE_ROUND16X4F32SASYM(d0, d1));

            AE_MINMAX16(z, 0, 255);

            AE_L8X4U_IP(x_, pIn, 4 * sizeof(int8_t));
            ix_ = AE_SRAI16(AE_SUB16(x_, 2), 2);
            AE_MINMAX16(ix_, 0, 62);
            x0_ = AE_ADD16(AE_SLAI16S(ix_, 2), 2);

            ix_ = AE_ADD16(ix_, ix_);
            y00_ = AE_L16_X((const ae_int16 *)tbl, AE_MOVAD16_3(ix_));
            y10_ = AE_L16_X((const ae_int16 *)(tbl + 1), AE_MOVAD16_3(ix_));
            y01_ = AE_L16_X((const ae_int16 *)tbl, AE_MOVAD16_2(ix_));
            y11_ = AE_L16_X((const ae_int16 *)(tbl + 1), AE_MOVAD16_2(ix_));
            y02_ = AE_L16_X((const ae_int16 *)tbl, AE_MOVAD16_1(ix_));
            y12_ = AE_L16_X((const ae_int16 *)(tbl + 1), AE_MOVAD16_1(ix_));
            y03_ = AE_L16_X((const ae_int16 *)tbl, AE_MOVAD16_0(ix_));
            y13_ = AE_L16_X((const ae_int16 *)(tbl + 1), AE_MOVAD16_0(ix_));

            AE_MOVT16X4(y00_, y01_, AE_MOVAB4(4));
            AE_MOVT16X4(y02_, y03_, AE_MOVAB4(4));
            y0_ = AE_SEL16_7632(y00_, y02_);

            AE_MOVT16X4(y10_, y11_, AE_MOVAB4(4));
            AE_MOVT16X4(y12_, y13_, AE_MOVAB4(4));
            y1_ = AE_SEL16_7632(y10_, y12_);

            AE_MULF16X4SS(d0_, d1_, AE_SLAI16S(AE_SUB16(y1_, y0_), 1), AE_SLAI16S(AE_SUB16(x_, x0_), 12));
            z_ = AE_ADD16S(y0_, AE_ROUND16X4F32SASYM(d0_, d1_));

            AE_MINMAX16(z_, 0, 255);
            t0 = AE_SEL8X8I(AE_MOVINT8X8_FROMINT16X4(z), AE_MOVINT8X8_FROMINT16X4(z_), 25); //Extract even
            AE_S8X8_IP(t0, pOut, sizeof(ae_int16x4));
        }

        for (j = (w&~7); j < w; j++)
        {
            uint8_t p = in[i*sz->stride + j];
            x = AE_MOVDA16(p);
            ix = AE_SRAI16(AE_SUB16(x, 2), 2);
            AE_MOVT16X4(ix, 62, AE_LT16(62, ix));
            AE_MOVT16X4(ix, 0, AE_LT16(ix, 0));

            x0 = AE_ADD16(AE_SLAI16S(ix, 2), 2);
            ix = AE_ADD16(ix, ix);
            y0 = AE_L16_X((const ae_int16 *)tbl, AE_MOVAD16_0(ix));
            y1 = AE_L16_X((const ae_int16 *)(tbl + 1), AE_MOVAD16_0(ix));

            AE_MULF16X4SS(d0, d1, AE_SLAI16S(AE_SUB16(y1, y0), 1), AE_SLAI16S(AE_SUB16(x, x0), 12));
            z = AE_ADD16(y0, AE_ROUND16X4F32SASYM(d0, d1));

            AE_MOVT16X4(z, 0, AE_LT16(z, 0));
            AE_MOVT16X4(z, 255, AE_LT16(255, z));

            out[i*sz->stride + j] = AE_MOVAD16_0(z);
        }
    }
}/*imgfastnorm_gs_nonlinear*/
