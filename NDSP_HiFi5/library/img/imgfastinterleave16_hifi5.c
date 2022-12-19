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
void imgfastinterleave16(    void * restrict outImg, 
                       const void * restrict inImgR, 
                       const void * restrict inImgG, 
                       const void * restrict inImgB, 
                       const imgsize_t* sz)
{   
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
    int i, j, w = (int)sz->width, h = (int)sz->height, stride = sz->stride;
    const ae_int16x8 * restrict pInR;
    const int16_t * restrict inR = (const int16_t *)inImgR;
    const ae_int16x8 * restrict pInG;
    const int16_t * restrict inG = (const int16_t *)inImgG;
    const ae_int16x8 * restrict pInB;
    const int16_t * restrict inB = (const int16_t *)inImgB;
    int16_t * restrict out = (int16_t *)outImg;
    ae_int16x8 * restrict pOut;
    imgsize_validate(sz, 2, 1);
    ae_int16x4 sel0, sel1, sel2;
    ae_int16x4 r, g, b, _12, _23, _31, _1, _2, _3;
    ae_int16x4 r_, g_, b_, _12_, _23_, _31_, _1_, _2_, _3_;
    ae_int64 temp64;
    NASSERT(outImg);
    NASSERT(inImgR);
    NASSERT(inImgG);
    NASSERT(inImgB);
    NASSERT_ALIGN(outImg, ALIGNMENT);
    NASSERT_ALIGN(inImgR, ALIGNMENT);
    NASSERT_ALIGN(inImgG, ALIGNMENT);
    NASSERT_ALIGN(inImgB, ALIGNMENT);
    temp64 = 0x0705030106040200;
    sel0 = AE_MOVINT16X4_FROMINT64(temp64);
    temp64 = 0x0704060203010500;
    sel1 = AE_MOVINT16X4_FROMINT64(temp64);
    temp64 = 0x0705060103000204;
    sel2 = AE_MOVINT16X4_FROMINT64(temp64);
    for (i = 0; i < h; i++)
    {
        pInR = (const ae_int16x8*)(inR + i*stride);
        pInG = (const ae_int16x8*)(inG + i*stride);
        pInB = (const ae_int16x8*)(inB + i*stride);
        pOut = (ae_int16x8 *)(out + 3 * i*stride);
        for (j = 0; j < (w >> 3); j++)
        {
            AE_L16X4X2_IP(r, r_, pInR, 8 * sizeof(int16_t));
            AE_L16X4X2_IP(g, g_, pInG, 8 * sizeof(int16_t));
            AE_L16X4X2_IP(b, b_, pInB, 8 * sizeof(int16_t));

            AE_DSEL16X4(_12, _23,   r,   g, sel0);
            AE_DSEL16X4( _1, _31, _12,   b, sel1);
            AE_DSEL16X4( _2,  _3, _31, _23, sel2);

            AE_DSEL16X4(_12_, _23_,   r_,   g_, sel0);
            AE_DSEL16X4( _1_, _31_, _12_,   b_, sel1);
            AE_DSEL16X4( _2_,  _3_, _31_, _23_, sel2);

            AE_S16X4X2_IP( _1,  _2, pOut, 2*sizeof(ae_int16x4));
            AE_S16X4X2_IP( _3, _1_, pOut, 2*sizeof(ae_int16x4));
            AE_S16X4X2_IP(_2_, _3_, pOut, 2*sizeof(ae_int16x4));
        }
        if (w&4)
        {          
            AE_L16X4_IP(r, castxcc(ae_int16x4, pInR), 4 * sizeof(int16_t));
            AE_L16X4_IP(g, castxcc(ae_int16x4, pInG), 4 * sizeof(int16_t));
            AE_L16X4_IP(b, castxcc(ae_int16x4, pInB), 4 * sizeof(int16_t));
            AE_DSEL16X4(_12, _23, r, g, sel0);
            AE_DSEL16X4(_1, _31, _12, b, sel1);
            AE_DSEL16X4(_2, _3, _31, _23, sel2);
            AE_S16X4_IP(_1, castxcc(ae_int16x4, pOut), sizeof(ae_int16x4));
            AE_S16X4_IP(_2, castxcc(ae_int16x4, pOut), sizeof(ae_int16x4));
            AE_S16X4_IP(_3, castxcc(ae_int16x4, pOut), sizeof(ae_int16x4));
        }
    }
    for (i = 0; i < h; i++)
    {
        for (j = w&(~3); j < w; j++)
        {
            int16_t r = inR[i*stride + j];
            int16_t g = inG[i*stride + j];
            int16_t b = inB[i*stride + j];
            out[i*stride * 3 + j * 3 + 0] = r;
            out[i*stride * 3 + j * 3 + 1] = g;
            out[i*stride * 3 + j * 3 + 2] = b;
        }
    }
#else
    int i, j, w = (int)sz->width, h = (int)sz->height, stride = sz->stride;
    const ae_int16x8 * restrict pInR;
    const int16_t * restrict inR = (const int16_t *)inImgR;
    const ae_int16x8 * restrict pInG;
    const int16_t * restrict inG = (const int16_t *)inImgG;
    const ae_int16x8 * restrict pInB;
    const int16_t * restrict inB = (const int16_t *)inImgB;
    int16_t * restrict out = (int16_t *)outImg;
    ae_int16x8 * restrict pOut;
    ae_valignx2 al_rw, al_gw, al_bw, al_outw;
    imgsize_validate(sz, 2, 1);
    ae_int16x4 sel0, sel1, sel2;
    ae_int16x4 r, g, b, _12, _23, _31, _1, _2, _3;
    ae_int16x4 r_, g_, b_, _12_, _23_, _31_, _1_, _2_, _3_;
    ae_int64 temp64;
    NASSERT(outImg);
    NASSERT(inImgR);
    NASSERT(inImgG);
    NASSERT(inImgB);
    NASSERT_ALIGN(outImg, ALIGNMENT);
    NASSERT_ALIGN(inImgR, ALIGNMENT);
    NASSERT_ALIGN(inImgG, ALIGNMENT);
    NASSERT_ALIGN(inImgB, ALIGNMENT);
    temp64 = 0x0705030106040200;
    sel0 = AE_MOVINT16X4_FROMINT64(temp64);
    temp64 = 0x0704060203010500;
    sel1 = AE_MOVINT16X4_FROMINT64(temp64);
    temp64 = 0x0705060103000204;
    sel2 = AE_MOVINT16X4_FROMINT64(temp64);
    al_outw = AE_ZALIGN128();
    for (i = 0; i < h; i++)
    {
        pInR = (const ae_int16x8*)(inR + i*stride);
        pInG = (const ae_int16x8*)(inG + i*stride);
        pInB = (const ae_int16x8*)(inB + i*stride);
        pOut = (ae_int16x8 *)(out + 3 * i*stride);
        for (j = 0; j < (w >> 3); j++)
        {
            AE_L16X4X2_IP(r, r_, pInR, 8 * sizeof(int16_t));
            AE_L16X4X2_IP(g, g_, pInG, 8 * sizeof(int16_t));
            AE_L16X4X2_IP(b, b_, pInB, 8 * sizeof(int16_t));

            AE_DSEL16X4(_12, _23,   r,   g, sel0);
            AE_DSEL16X4( _1, _31, _12,   b, sel1);
            AE_DSEL16X4( _2,  _3, _31, _23, sel2);

            AE_DSEL16X4(_12_, _23_,   r_,   g_, sel0);
            AE_DSEL16X4( _1_, _31_, _12_,   b_, sel1);
            AE_DSEL16X4( _2_,  _3_, _31_, _23_, sel2);

            AE_S16X4X2_IP( _1,  _2, pOut, 2*sizeof(ae_int16x4));
            AE_S16X4X2_IP( _3, _1_, pOut, 2*sizeof(ae_int16x4));
            AE_S16X4X2_IP(_2_, _3_, pOut, 2*sizeof(ae_int16x4));
        }
        if (w&7)
        {             
            int off = (w&7)<<1;
            al_rw = AE_LA128_PP(pInR);
            al_gw = AE_LA128_PP(pInG);
            al_bw = AE_LA128_PP(pInB);         
            AE_LAV16X4X2_XP(r, r_, al_rw, pInR, off);
            AE_LAV16X4X2_XP(g, g_, al_gw, pInG, off);
            AE_LAV16X4X2_XP(b, b_, al_bw, pInB, off);
            off = off*3;
            AE_DSEL16X4(_12, _23,   r,   g, sel0);
            AE_DSEL16X4( _1, _31, _12,   b, sel1);
            AE_DSEL16X4( _2,  _3, _31, _23, sel2);

            AE_DSEL16X4(_12_, _23_,   r_,   g_, sel0);
            AE_DSEL16X4( _1_, _31_, _12_,   b_, sel1);
            AE_DSEL16X4( _2_,  _3_, _31_, _23_, sel2);

            AE_SAV16X4X2_XP( _1,  _2, al_outw, pOut, XT_MIN(16, off));
            AE_SAV16X4X2_XP( _3, _1_, al_outw, pOut, XT_MIN(16, XT_MAX(0, off-16)));
            AE_SAV16X4X2_XP(_2_, _3_, al_outw, pOut, XT_MAX(0, off-32));
            AE_SA128POS_FP(al_outw, pOut);
        }
    }

#endif
}
