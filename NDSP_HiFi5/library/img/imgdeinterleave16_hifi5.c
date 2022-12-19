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
void imgdeinterleave16(       void * restrict outImgR, 
                              void * restrict outImgG, 
                              void * restrict outImgB, 
                        const void * restrict inImg  , 
                        const imgsize_t* sz)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
    int i, j, w = (int)sz->width, h = (int)sz->height, stride = sz->stride;
    ae_int16x8 * restrict pOutR;
    int16_t * restrict outR = (int16_t *)outImgR;
    ae_int16x8 * restrict pOutG;
    int16_t * restrict outG = (int16_t *)outImgG;
    ae_int16x8 * restrict pOutB;
    int16_t * restrict outB = (int16_t *)outImgB;
    const ae_int16x8 * restrict pIn;
    const int16_t * restrict in = (const int16_t *)inImg;
    ae_valign al_r, al_g, al_b, al_in;
    ae_valignx2 al_rw, al_gw, al_bw, al_inw;
    ae_int16x4 r, g, b, _1, _2, _3, _rb, _gb, _bg;
    ae_int16x4 r_, g_, b_, _1_, _2_, _3_, _rb_, _gb_, _bg_;
    ae_int16x4 sel0, sel1, sel2;
    ae_int64 temp64;
    imgsize_validate(sz, 2, 0);
    NASSERT(outImgR);
    NASSERT(outImgG);
    NASSERT(outImgB);
    NASSERT(inImg);
    NASSERT_ALIGN(outImgR, 2);
    NASSERT_ALIGN(outImgG, 2);
    NASSERT_ALIGN(outImgB, 2);
    NASSERT_ALIGN(inImg, 2);
    temp64 = 0x0706040301000502;
    sel0 = AE_MOVINT16X4_FROMINT64(temp64);
    temp64 = 0x0704060305000201;
    sel1 = AE_MOVINT16X4_FROMINT64(temp64);
    temp64 = 0x0703060405020001;
    sel2 = AE_MOVINT16X4_FROMINT64(temp64);
    al_rw = AE_ZALIGN128();
    al_gw = AE_ZALIGN128();
    al_bw = AE_ZALIGN128();
    for (i = 0; i < h; i++)
    {
        pOutR = (ae_int16x8*)(outR + i*stride);
        pOutG = (ae_int16x8*)(outG + i*stride);
        pOutB = (ae_int16x8*)(outB + i*stride);
        pIn = (const ae_int16x8 *)(in + 3 * i*stride);
        al_inw = AE_LA128_PP(pIn);
        for (j = 0; j < (w >> 3); j++)
        {
            AE_LA16X4X2_IP(_1, _2,  al_inw, pIn);
            AE_LA16X4X2_IP(_3, _1_, al_inw, pIn);
            AE_LA16X4X2_IP(_2_, _3_,al_inw, pIn);

            AE_DSEL16X4(_rb, _gb,  _1,  _2, sel0);
            AE_DSEL16X4(  r, _bg, _rb,  _3, sel1);
            AE_DSEL16X4(  g,   b, _gb, _bg, sel2);

            AE_DSEL16X4(_rb_, _gb_,  _1_,  _2_, sel0);
            AE_DSEL16X4(  r_, _bg_, _rb_,  _3_, sel1);
            AE_DSEL16X4(  g_,   b_, _gb_, _bg_, sel2);

            AE_SA16X4X2_IP(r, r_, al_rw, pOutR);
            AE_SA16X4X2_IP(g, g_, al_gw, pOutG);
            AE_SA16X4X2_IP(b, b_, al_bw, pOutB);
        }
        AE_SA128POS_FP(al_rw, pOutR);
        AE_SA128POS_FP(al_gw, pOutG);
        AE_SA128POS_FP(al_bw, pOutB);
        if (w&4)
        {
            al_in = AE_LA64_PP(pIn);
            al_r = AE_ZALIGN64();
            al_g = AE_ZALIGN64();
            al_b = AE_ZALIGN64();
            AE_LA16X4_IP(_1, al_in, castxcc(ae_int16x4, pIn));
            AE_LA16X4_IP(_2, al_in, castxcc(ae_int16x4, pIn));
            AE_LA16X4_IP(_3, al_in, castxcc(ae_int16x4, pIn));

            AE_DSEL16X4(_rb, _gb,  _1,  _2, sel0);
            AE_DSEL16X4(  r, _bg, _rb,  _3, sel1);
            AE_DSEL16X4(  g,   b, _gb, _bg, sel2);

            AE_SA16X4_IP(r, al_r, castxcc(ae_int16x4, pOutR));
            AE_SA16X4_IP(g, al_g, castxcc(ae_int16x4, pOutG));
            AE_SA16X4_IP(b, al_b, castxcc(ae_int16x4, pOutB));
            AE_SA64POS_FP(al_r, pOutR);
            AE_SA64POS_FP(al_g, pOutG);
            AE_SA64POS_FP(al_b, pOutB);
        }
        
    }
    for (i = 0; i < h; i++)
    {
        for (j = w&(~3); j < w; j++)
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
    int i, j, w = (int)sz->width, h = (int)sz->height, stride = sz->stride;
    ae_int16x8 * restrict pOutR;
    int16_t * restrict outR = (int16_t *)outImgR;
    ae_int16x8 * restrict pOutG;
    int16_t * restrict outG = (int16_t *)outImgG;
    ae_int16x8 * restrict pOutB;
    int16_t * restrict outB = (int16_t *)outImgB;
    const ae_int16x8 * restrict pIn;
    const int16_t * restrict in = (const int16_t *)inImg;
    ae_valignx2 al_rw, al_gw, al_bw, al_inw;
    ae_int16x4 r, g, b, _1, _2, _3, _rb, _gb, _bg;
    ae_int16x4 r_, g_, b_, _1_, _2_, _3_, _rb_, _gb_, _bg_;
    ae_int16x4 sel0, sel1, sel2;
    ae_int64 temp64;
    imgsize_validate(sz, 2, 0);
    NASSERT(outImgR);
    NASSERT(outImgG);
    NASSERT(outImgB);
    NASSERT(inImg);
    NASSERT_ALIGN(outImgR, 2);
    NASSERT_ALIGN(outImgG, 2);
    NASSERT_ALIGN(outImgB, 2);
    NASSERT_ALIGN(inImg, 2);
    temp64 = 0x0706040301000502;
    sel0 = AE_MOVINT16X4_FROMINT64(temp64);
    temp64 = 0x0704060305000201;
    sel1 = AE_MOVINT16X4_FROMINT64(temp64);
    temp64 = 0x0703060405020001;
    sel2 = AE_MOVINT16X4_FROMINT64(temp64);
    al_rw = AE_ZALIGN128();
    al_gw = AE_ZALIGN128();
    al_bw = AE_ZALIGN128();
    for (i = 0; i < h; i++)
    {
        pOutR = (ae_int16x8*)(outR + i*stride);
        pOutG = (ae_int16x8*)(outG + i*stride);
        pOutB = (ae_int16x8*)(outB + i*stride);
        pIn = (const ae_int16x8 *)(in + 3 * i*stride);
        al_inw = AE_LA128_PP(pIn);
        for (j = 0; j < (w >> 3); j++)
        {
            AE_LA16X4X2_IP(_1, _2,  al_inw, pIn);
            AE_LA16X4X2_IP(_3, _1_, al_inw, pIn);
            AE_LA16X4X2_IP(_2_, _3_,al_inw, pIn);

            AE_DSEL16X4(_rb, _gb,  _1,  _2, sel0);
            AE_DSEL16X4(  r, _bg, _rb,  _3, sel1);
            AE_DSEL16X4(  g,   b, _gb, _bg, sel2);

            AE_DSEL16X4(_rb_, _gb_,  _1_,  _2_, sel0);
            AE_DSEL16X4(  r_, _bg_, _rb_,  _3_, sel1);
            AE_DSEL16X4(  g_,   b_, _gb_, _bg_, sel2);

            AE_SA16X4X2_IP(r, r_, al_rw, pOutR);
            AE_SA16X4X2_IP(g, g_, al_gw, pOutG);
            AE_SA16X4X2_IP(b, b_, al_bw, pOutB);
        }       
        if (w&7)
        {
            int off = (w&7)<<1;
            int off3 = 3*off;
            AE_LAV16X4X2_XP(_1, _2,  al_inw, pIn, XT_MIN(16, off3));
            AE_LAV16X4X2_XP(_3, _1_, al_inw, pIn, XT_MIN(16, XT_MAX(0, off3-16)));
            AE_LAV16X4X2_XP(_2_, _3_,al_inw, pIn, XT_MAX(0, off3-32));          

            AE_DSEL16X4(_rb, _gb,  _1,  _2, sel0);
            AE_DSEL16X4(  r, _bg, _rb,  _3, sel1);
            AE_DSEL16X4(  g,   b, _gb, _bg, sel2);

            AE_DSEL16X4(_rb_, _gb_,  _1_,  _2_, sel0);
            AE_DSEL16X4(  r_, _bg_, _rb_,  _3_, sel1);
            AE_DSEL16X4(  g_,   b_, _gb_, _bg_, sel2);

            AE_SAV16X4X2_XP(r, r_, al_rw, pOutR, off);
            AE_SAV16X4X2_XP(g, g_, al_gw, pOutG, off);
            AE_SAV16X4X2_XP(b, b_, al_bw, pOutB, off);
        }
        AE_SA128POS_FP(al_rw, pOutR);
        AE_SA128POS_FP(al_gw, pOutG);
        AE_SA128POS_FP(al_bw, pOutB);

    }
#endif
}
