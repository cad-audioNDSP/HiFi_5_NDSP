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
  Image interleave
  Functions convert planar images to packed format 

  Image formats:
  8-bit signed or unsigned data
  16-bit signed data

  Input:
  inImgR
  inImgG
  inImgB  planes with R,G,B components
  sz      image size
  Output:
  outImg  packed image (RGB come together)

  Restrictions:
  see general restrictions applied for all images for fast/generic 
  functions
-------------------------------------------------------------------------*/
void imgfastinterleave(      void * restrict outImg, 
                       const void * restrict inImgR, 
                       const void * restrict inImgG, 
                       const void * restrict inImgB, 
                       const imgsize_t* sz)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
    int i, j, w = (int)sz->width, h = (int)sz->height, stride = sz->stride;
    const ae_int8x16 * restrict pInR;
    const uint8_t * restrict inR = (const uint8_t *)inImgR;
    const ae_int8x16 * restrict pInG;
    const uint8_t * restrict inG = (const uint8_t *)inImgG;
    const ae_int8x16 * restrict pInB;
    const uint8_t * restrict inB = (const uint8_t *)inImgB;
    uint8_t * restrict out = (uint8_t *)outImg;
    ae_int8x16 * restrict pOut;
    ae_int8x8 r, g, b;
    ae_int8x8 r_, g_, b_;
    ae_int8x8 _12, _23, _32, _1, _2, _3;
    ae_int8x8 _12_, _23_, _32_, _1_, _2_, _3_;
    ae_int8x8 sel0, sel1, sel2;
    ae_int64 temp64;
    imgsize_validate(sz, 1, 1);
    NASSERT(outImg);
    NASSERT(inImgR);
    NASSERT(inImgG);
    NASSERT(inImgB);
    NASSERT_ALIGN(outImg, ALIGNMENT);
    NASSERT_ALIGN(inImgR, ALIGNMENT);
    NASSERT_ALIGN(inImgG, ALIGNMENT);
    NASSERT_ALIGN(inImgB, ALIGNMENT);
    temp64 = 0xfb73ea62d951c840;
    sel0 = AE_MOVINT8X8_FROMINT64(temp64);
    temp64 = 0xf9e875d4c362b1a0;
    sel1 = AE_MOVINT8X8_FROMINT64(temp64);
    temp64 = 0x5c726b4af1e938d0;
    sel2 = AE_MOVINT8X8_FROMINT64(temp64);
    for (i = 0; i < h; i++)
    {
        pInR = (const ae_int8x16*)(inR + i*stride);
        pInG = (const ae_int8x16*)(inG + i*stride);
        pInB = (const ae_int8x16*)(inB + i*stride);
        pOut = (ae_int8x16 *)(out + 3 * i*stride);
        for (j = 0; j < (w >> 4); j++)
        {
            AE_L8X8X2_IP(r, r_, pInR, 2*sizeof(ae_int8x8));
            AE_L8X8X2_IP(g, g_, pInG, 2*sizeof(ae_int8x8));
            AE_L8X8X2_IP(b, b_, pInB, 2*sizeof(ae_int8x8));

            AE_DSEL8X8(_12, _23,   r,   g, sel0);
            AE_DSEL8X8( _1, _32, _12,   b, sel1);
            AE_DSEL8X8( _2,  _3, _23, _32, sel2);

            AE_DSEL8X8(_12_, _23_,   r_,   g_, sel0);
            AE_DSEL8X8( _1_, _32_, _12_,   b_, sel1);
            AE_DSEL8X8( _2_,  _3_, _23_, _32_, sel2);
           
            AE_S8X8X2_IP(_1, _2,  pOut, 2*sizeof(ae_int8x8));
            AE_S8X8X2_IP(_3, _1_, pOut, 2*sizeof(ae_int8x8));
            AE_S8X8X2_IP(_2_, _3_,pOut, 2*sizeof(ae_int8x8));
        }
        if (w & 8)
        {
            AE_L8X8_IP(r, castxcc(ae_int8x8, pInR), sizeof(ae_int8x8));
            AE_L8X8_IP(g, castxcc(ae_int8x8, pInG), sizeof(ae_int8x8));
            AE_L8X8_IP(b, castxcc(ae_int8x8, pInB), sizeof(ae_int8x8));

            AE_DSEL8X8(_12, _23,   r,   g, sel0);
            AE_DSEL8X8( _1, _32, _12,   b, sel1);
            AE_DSEL8X8( _2,  _3, _23, _32, sel2);

            AE_S8X8_IP(_1, castxcc(ae_int8x8, pOut), sizeof(ae_int8x8));
            AE_S8X8_IP(_2, castxcc(ae_int8x8, pOut), sizeof(ae_int8x8));
            AE_S8X8_IP(_3, castxcc(ae_int8x8, pOut), sizeof(ae_int8x8));
        }
    }
    for (i = 0; i < h; i++)
    {
        for (j = w&(~7); j < w; j++)
        {
            uint8_t r = ((const uint8_t*)inImgR)[i*stride + j];
            uint8_t g = ((const uint8_t*)inImgG)[i*stride + j];
            uint8_t b = ((const uint8_t*)inImgB)[i*stride + j];
            ((uint8_t*)outImg)[i*stride * 3 + j * 3 + 0] = r;
            ((uint8_t*)outImg)[i*stride * 3 + j * 3 + 1] = g;
            ((uint8_t*)outImg)[i*stride * 3 + j * 3 + 2] = b;
        }
    }
#else
    int i, j, w = (int)sz->width, h = (int)sz->height, stride = sz->stride;
    const ae_int8x16 * restrict pInR;
    const uint8_t * restrict inR = (const uint8_t *)inImgR;
    const ae_int8x16 * restrict pInG;
    const uint8_t * restrict inG = (const uint8_t *)inImgG;
    const ae_int8x16 * restrict pInB;
    const uint8_t * restrict inB = (const uint8_t *)inImgB;
    uint8_t * restrict out = (uint8_t *)outImg;
    ae_int8x16 * restrict pOut;
    ae_int8x8 r, g, b;
    ae_int8x8 r_, g_, b_;
    ae_int8x8 _12, _23, _32, _1, _2, _3;
    ae_int8x8 _12_, _23_, _32_, _1_, _2_, _3_;
    ae_int8x8 sel0, sel1, sel2;
    ae_valignx2 al_rw, al_gw, al_bw, al_outw;
    ae_int64 temp64;
    imgsize_validate(sz, 1, 1);
    NASSERT(outImg);
    NASSERT(inImgR);
    NASSERT(inImgG);
    NASSERT(inImgB);
    NASSERT_ALIGN(outImg, ALIGNMENT);
    NASSERT_ALIGN(inImgR, ALIGNMENT);
    NASSERT_ALIGN(inImgG, ALIGNMENT);
    NASSERT_ALIGN(inImgB, ALIGNMENT);
    temp64 = 0xfb73ea62d951c840;
    sel0 = AE_MOVINT8X8_FROMINT64(temp64);
    temp64 = 0xf9e875d4c362b1a0;
    sel1 = AE_MOVINT8X8_FROMINT64(temp64);
    temp64 = 0x5c726b4af1e938d0;
    sel2 = AE_MOVINT8X8_FROMINT64(temp64);
    al_outw = AE_ZALIGN128();
    for (i = 0; i < h; i++)
    {
        pInR = (const ae_int8x16*)(inR + i*stride);
        pInG = (const ae_int8x16*)(inG + i*stride);
        pInB = (const ae_int8x16*)(inB + i*stride);
        pOut = (ae_int8x16 *)(out + 3 * i*stride);
        for (j = 0; j < (w >> 4); j++)
        {
            AE_L8X8X2_IP(r, r_, pInR, 2*sizeof(ae_int8x8));
            AE_L8X8X2_IP(g, g_, pInG, 2*sizeof(ae_int8x8));
            AE_L8X8X2_IP(b, b_, pInB, 2*sizeof(ae_int8x8));

            AE_DSEL8X8(_12, _23,   r,   g, sel0);
            AE_DSEL8X8( _1, _32, _12,   b, sel1);
            AE_DSEL8X8( _2,  _3, _23, _32, sel2);

            AE_DSEL8X8(_12_, _23_,   r_,   g_, sel0);
            AE_DSEL8X8( _1_, _32_, _12_,   b_, sel1);
            AE_DSEL8X8( _2_,  _3_, _23_, _32_, sel2);

            AE_S8X8X2_IP(_1, _2,  pOut, 2*sizeof(ae_int8x8));
            AE_S8X8X2_IP(_3, _1_, pOut, 2*sizeof(ae_int8x8));
            AE_S8X8X2_IP(_2_, _3_,pOut, 2*sizeof(ae_int8x8));
        }
        if (w & 15)
        {
            int off = (w&15);
            al_rw = AE_LA128_PP(pInR);
            al_gw = AE_LA128_PP(pInG);
            al_bw = AE_LA128_PP(pInB);  
            AE_LAV8X8X2_XP(r, r_, al_rw, pInR, off);
            AE_LAV8X8X2_XP(g, g_, al_gw, pInG, off);
            AE_LAV8X8X2_XP(b, b_, al_bw, pInB, off);
            off = 3*off;
            AE_DSEL8X8(_12, _23,   r,   g, sel0);
            AE_DSEL8X8( _1, _32, _12,   b, sel1);
            AE_DSEL8X8( _2,  _3, _23, _32, sel2);

            AE_DSEL8X8(_12_, _23_,   r_,   g_, sel0);
            AE_DSEL8X8( _1_, _32_, _12_,   b_, sel1);
            AE_DSEL8X8( _2_,  _3_, _23_, _32_, sel2);

            AE_SAV8X8X2_XP(_1, _2,   al_outw, pOut, XT_MIN(16, off));
            AE_SAV8X8X2_XP(_3, _1_,  al_outw, pOut, XT_MIN(16, XT_MAX(0, off-16)));
            AE_SAV8X8X2_XP(_2_, _3_, al_outw, pOut, XT_MAX(0, off-32));
            AE_SA128POS_FP(al_outw, pOut);
        }
    }
#endif
}
