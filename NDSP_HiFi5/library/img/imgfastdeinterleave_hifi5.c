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
  Image deinterleave
  Functions convert packed images to planar format 

  Image formats:
  8-bit signed or unsigned data
  16-bit signed data

  Input:
  inImg   packed image (RGB come together)
  sz      image size
  Output:
  outImgR
  outImgG
  outImgB planes with R,G,B components

  Restrictions:
  see general restrictions applied for all images for fast/generic 
  functions
-------------------------------------------------------------------------*/
void imgfastdeinterleave  
                      (       void * restrict outImgR, 
                              void * restrict outImgG, 
                              void * restrict outImgB, 
                        const void * restrict inImg, 
                        const imgsize_t* sz)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
    int i,j,w=(int)sz->width,h=(int)sz->height,stride=sz->stride;
    ae_int8x16 * restrict pOutR;
    uint8_t * restrict outR = (uint8_t *)outImgR;
    ae_int8x16 * restrict pOutG;
    uint8_t * restrict outG = (uint8_t *)outImgG;
    ae_int8x16 * restrict pOutB;
    uint8_t * restrict outB = (uint8_t *)outImgB;
    const uint8_t * restrict in = (const uint8_t *)inImg;
    const ae_int8x16 * restrict pIn;
    ae_int8x8 r, g, b;
    ae_int8x8 r_, g_, b_;
    ae_int8x8 _12, _23, _32, _1, _2, _3;
    ae_int8x8 _12_, _23_, _32_, _1_, _2_, _3_;
    ae_int8x8 sel0, sel1, sel2;
    ae_int64 temp64;
    imgsize_validate(sz,1,1);
    NASSERT(outImgR);
    NASSERT(outImgG);
    NASSERT(outImgB);
    NASSERT(inImg  );
    NASSERT_ALIGN(outImgR,ALIGNMENT);
    NASSERT_ALIGN(outImgG,ALIGNMENT);
    NASSERT_ALIGN(outImgB,ALIGNMENT);
    NASSERT_ALIGN(inImg  ,ALIGNMENT);
    temp64 = 0xbead8f7c59462310;
    sel0 = AE_MOVINT8X8_FROMINT64(temp64);
    temp64 = 0xfdeac5b493827160;
    sel1 = AE_MOVINT8X8_FROMINT64(temp64);
    temp64 = 0xfedcba9876543210;
    sel2 = AE_MOVINT8X8_FROMINT64(temp64);
    for (i = 0; i < h; i++)
    {
        pOutR = (ae_int8x16*)(outR + i*stride);
        pOutG = (ae_int8x16*)(outG + i*stride);
        pOutB = (ae_int8x16*)(outB + i*stride);
        pIn = (const ae_int8x16 *)(in + 3 * i*stride);
        for (j = 0; j < (w >> 4); j++)
        {
            AE_L8X8X2_IP(_1,  _2,  pIn, 2*sizeof(ae_int8x8));
            AE_L8X8X2_IP(_3,  _1_, pIn, 2*sizeof(ae_int8x8));
            AE_L8X8X2_IP(_2_, _3_, pIn, 2*sizeof(ae_int8x8));

            AE_DSEL8X8(_23, _32,  _2,  _3, sel0);
            AE_DSEL8X8(_12,   b,  _1, _32, sel1);            
            AE_DSEL8X8(  r,   g, _12, _23, sel2);

            AE_DSEL8X8(_23_, _32_,  _2_,  _3_, sel0);
            AE_DSEL8X8(_12_,   b_,  _1_, _32_, sel1);            
            AE_DSEL8X8(  r_,   g_, _12_, _23_, sel2);

            AE_S8X8X2_IP(r, r_, pOutR, 2*sizeof(ae_int8x8));
            AE_S8X8X2_IP(g, g_, pOutG, 2*sizeof(ae_int8x8));
            AE_S8X8X2_IP(b, b_, pOutB, 2*sizeof(ae_int8x8));
        }
        if (w&8)
        {
            AE_L8X8_IP(_1, castxcc(ae_int8x8, pIn), sizeof(ae_int8x8));
            AE_L8X8_IP(_2, castxcc(ae_int8x8, pIn), sizeof(ae_int8x8));
            AE_L8X8_IP(_3, castxcc(ae_int8x8, pIn), sizeof(ae_int8x8));

            AE_DSEL8X8(_23, _32,  _2,  _3, sel0);
            AE_DSEL8X8(_12,   b,  _1, _32, sel1);            
            AE_DSEL8X8(  r,   g, _12, _23, sel2);

            AE_S8X8_IP(r, castxcc(ae_int8x8, pOutR), sizeof(ae_int8x8));
            AE_S8X8_IP(g, castxcc(ae_int8x8, pOutG), sizeof(ae_int8x8));
            AE_S8X8_IP(b, castxcc(ae_int8x8, pOutB), sizeof(ae_int8x8));
        }
    }
    for (i = 0; i < h; i++)
    {
        for (j = w&(~7); j < w; j++)
        {
            int16_t r = in[i*stride * 3 + j * 3 + 0];
            int16_t g = in[i*stride * 3 + j * 3 + 1];
            int16_t b = in[i*stride * 3 + j * 3 + 2];
            outR[i*stride + j] = r;
            outG[i*stride + j] = g;
            outB[i*stride + j] = b;
        }
    }
#else
    int i,j,w=(int)sz->width,h=(int)sz->height,stride=sz->stride;
    ae_int8x16 * restrict pOutR;
    uint8_t * restrict outR = (uint8_t *)outImgR;
    ae_int8x16 * restrict pOutG;
    uint8_t * restrict outG = (uint8_t *)outImgG;
    ae_int8x16 * restrict pOutB;
    uint8_t * restrict outB = (uint8_t *)outImgB;
    const uint8_t * restrict in = (const uint8_t *)inImg;
    const ae_int8x16 * restrict pIn;
    ae_int8x8 r, g, b;
    ae_int8x8 r_, g_, b_;
    ae_int8x8 _12, _23, _32, _1, _2, _3;
    ae_int8x8 _12_, _23_, _32_, _1_, _2_, _3_;
    ae_int8x8 sel0, sel1, sel2;
    ae_valignx2 al_rw, al_gw, al_bw, al_inw;
    ae_int64 temp64;
    imgsize_validate(sz,1,1);
    NASSERT(outImgR);
    NASSERT(outImgG);
    NASSERT(outImgB);
    NASSERT(inImg  );
    NASSERT_ALIGN(outImgR,ALIGNMENT);
    NASSERT_ALIGN(outImgG,ALIGNMENT);
    NASSERT_ALIGN(outImgB,ALIGNMENT);
    NASSERT_ALIGN(inImg  ,ALIGNMENT);
    temp64 = 0xbead8f7c59462310;
    sel0 = AE_MOVINT8X8_FROMINT64(temp64);
    temp64 = 0xfdeac5b493827160;
    sel1 = AE_MOVINT8X8_FROMINT64(temp64);
    temp64 = 0xfedcba9876543210;
    sel2 = AE_MOVINT8X8_FROMINT64(temp64);
    for (i = 0; i < h; i++)
    {
        pOutR = (ae_int8x16*)(outR + i*stride);
        pOutG = (ae_int8x16*)(outG + i*stride);
        pOutB = (ae_int8x16*)(outB + i*stride);
        pIn = (const ae_int8x16 *)(in + 3 * i*stride);
        for (j = 0; j < (w >> 4); j++)
        {
            AE_L8X8X2_IP(_1,  _2,  pIn, 2*sizeof(ae_int8x8));
            AE_L8X8X2_IP(_3,  _1_, pIn, 2*sizeof(ae_int8x8));
            AE_L8X8X2_IP(_2_, _3_, pIn, 2*sizeof(ae_int8x8));

            AE_DSEL8X8(_23, _32,  _2,  _3, sel0);
            AE_DSEL8X8(_12,   b,  _1, _32, sel1);            
            AE_DSEL8X8(  r,   g, _12, _23, sel2);

            AE_DSEL8X8(_23_, _32_,  _2_,  _3_, sel0);
            AE_DSEL8X8(_12_,   b_,  _1_, _32_, sel1);            
            AE_DSEL8X8(  r_,   g_, _12_, _23_, sel2);

            AE_S8X8X2_IP(r, r_, pOutR, 2*sizeof(ae_int8x8));
            AE_S8X8X2_IP(g, g_, pOutG, 2*sizeof(ae_int8x8));
            AE_S8X8X2_IP(b, b_, pOutB, 2*sizeof(ae_int8x8));
        }
        if (w&15)
        {
            int off = (w&15);
            int off3 = 3*off;
            al_inw = AE_LA128_PP(pIn);
            AE_LAV8X8X2_XP(_1,  _2,  al_inw, pIn, XT_MIN(16, off3));
            AE_LAV8X8X2_XP(_3,  _1_, al_inw, pIn, XT_MIN(16, XT_MAX(0, off3-16)));
            AE_LAV8X8X2_XP(_2_, _3_, al_inw, pIn, XT_MAX(0, off3-32));          

            AE_DSEL8X8(_23, _32,  _2,  _3, sel0);
            AE_DSEL8X8(_12,   b,  _1, _32, sel1);            
            AE_DSEL8X8(  r,   g, _12, _23, sel2);

            AE_DSEL8X8(_23_, _32_,  _2_,  _3_, sel0);
            AE_DSEL8X8(_12_,   b_,  _1_, _32_, sel1);            
            AE_DSEL8X8(  r_,   g_, _12_, _23_, sel2);

            AE_SAV8X8X2_XP(r, r_, al_rw, pOutR, off);
            AE_SAV8X8X2_XP(g, g_, al_gw, pOutG, off);
            AE_SAV8X8X2_XP(b, b_, al_bw, pOutB, off);
            AE_SA128POS_FP(al_rw, pOutR);
            AE_SA128POS_FP(al_gw, pOutG);
            AE_SA128POS_FP(al_bw, pOutB);
        }
    }
#endif
}
