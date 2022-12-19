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
  RGB<->YUV conversion 
  Functions convert RGB planar images to YUV planar. 

  Y=c0*R+c1*G+c2*B
  U=c3*R+c4*G+c5*B + bias
  V=c6*R+c7*G+c8*B + bias
  where bias is 128 for 8-bit data and 16384 for 16-bit data

  Image formats:
  8-bit unsigned RGB/YUV
  16-bit signed RGB/YUV

  Input:
  inImgR,
  inImgG,
  inImgB  planes with R,G,B components
  c[13]   transform coefficients, Q29 format
  sz      image size
  Output:
  outImgY
  outImgU
  outImgV planes with Y,U,V components

  Restrictions:
  see general restrictions applied for all images for fast/generic 
  functions
-------------------------------------------------------------------------*/
void imgfastconvert_rgbyuv( void * restrict outImgY, 
                        void * restrict outImgU, 
                        void * restrict outImgV, 
                        const void * restrict inImgR, 
                        const void * restrict inImgG, 
                        const void * restrict inImgB, 
                        const int32_t * restrict c,
                        const imgsize_t* sz)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
    const signed char * restrict pInR;
    const uint8_t * restrict inR = (const uint8_t *)inImgR;
    const signed char * restrict pInG;
    const uint8_t * restrict inG = (const uint8_t *)inImgG;
    const signed char * restrict pInB;
    const uint8_t * restrict inB = (const uint8_t *)inImgB;
    ae_int8x8 * restrict pOutY;
    uint8_t * restrict outY = (uint8_t *)outImgY;
    ae_int8x8 * restrict pOutU;
    uint8_t * restrict outU = (uint8_t *)outImgU;
    ae_int8x8 * restrict pOutV;
    uint8_t * restrict outV = (uint8_t *)outImgV;
    ae_int16x4 r, g, b, y0, u0, v0, y1, u1, v1;//, t0, t1;
    ae_int8x8 t;
    ae_int32x2 dr0, dr1, dg0, dg1, db0, db1, tmp;
    int i, j, w = (int)sz->width, h = (int)sz->height, stride = sz->stride;
    ae_int16x4 c0, c1, c2, c3, c4, c5, c6, c7, c8;        
    
    NASSERT(outImgY);
    NASSERT(outImgU);
    NASSERT(outImgV);
    NASSERT(inImgR);
    NASSERT(inImgG);
    NASSERT(inImgB);
    NASSERT_ALIGN(outImgY, ALIGNMENT);
    NASSERT_ALIGN(outImgU, ALIGNMENT);
    NASSERT_ALIGN(outImgV, ALIGNMENT);
    NASSERT_ALIGN(inImgR, ALIGNMENT);
    NASSERT_ALIGN(inImgG, ALIGNMENT);
    NASSERT_ALIGN(inImgB, ALIGNMENT);
    imgsize_validate(sz, 1, 1);

    tmp = AE_L32_X((const ae_int32 *)c, 0 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c0 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 1 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c1 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 2 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c2 = AE_ROUND16X4F32SASYM(tmp, tmp);  
    tmp = AE_L32_X((const ae_int32 *)c, 3 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c3 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 4 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c4 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 5 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c5 = AE_ROUND16X4F32SASYM(tmp, tmp);   
    tmp = AE_L32_X((const ae_int32 *)c, 6 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c6 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 7 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c7 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 8 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c8 = AE_ROUND16X4F32SASYM(tmp, tmp);
    for (i = 0; i < h; i++)
    {
        pInR = (const signed char*)(inR + i*stride);
        pInG = (const signed char*)(inG + i*stride);
        pInB = (const signed char*)(inB + i*stride);
        pOutY = (ae_int8x8 *)(outY + i*stride);
        pOutU = (ae_int8x8 *)(outU + i*stride);
        pOutV = (ae_int8x8 *)(outV + i*stride);
        for (j = 0; j < (w >> 3); j++)
        {
            AE_L8X4U_IP(r, pInR, 4 * sizeof(uint8_t));
            AE_L8X4U_IP(g, pInG, 4 * sizeof(uint8_t));
            AE_L8X4U_IP(b, pInB, 4 * sizeof(uint8_t));

            AE_MULFD16X16X4WS(dr0, dr1, r, g, c0, c1);
            AE_MULAF16X4SS(dr0, dr1, b, c2);
            y0 = AE_ROUND16X4F32SASYM(dr0, dr1);

            AE_MULFD16X16X4WS(dr0, dr1, r, g, c3, c4);
            AE_MULAF16X4SS(dr0, dr1, b, c5);
            u0 = AE_ROUND16X4F32SASYM(dr0, dr1);
            u0 = AE_ADD16S(u0, (ae_int16x4)128);
            AE_MINMAX16(u0, 0, 255);

            AE_MULFD16X16X4WS(dr0, dr1, r, g, c6, c7);
            AE_MULAF16X4SS(dr0, dr1, b, c8);
            v0 = AE_ROUND16X4F32SASYM(dr0, dr1);
            v0 = AE_ADD16S(v0, (ae_int16x4)128);
            AE_MINMAX16(v0, 0, 255);

            /*second 4*/

            AE_L8X4U_IP(r, pInR, 4 * sizeof(uint8_t));
            AE_L8X4U_IP(g, pInG, 4 * sizeof(uint8_t));
            AE_L8X4U_IP(b, pInB, 4 * sizeof(uint8_t));

            AE_MULFD16X16X4WS(dr0, dr1, r, g, c0, c1);
            AE_MULAF16X4SS(dr0, dr1, b, c2);
            y1 = AE_ROUND16X4F32SASYM(dr0, dr1);

            AE_MULFD16X16X4WS(dr0, dr1, r, g, c3, c4);
            AE_MULAF16X4SS(dr0, dr1, b, c5);
            u1 = AE_ROUND16X4F32SASYM(dr0, dr1);
            u1 = AE_ADD16S(u1, (ae_int16x4)128);
            AE_MINMAX16(u1, 0, 255); 
            
            AE_MULFD16X16X4WS(dr0, dr1, r, g, c6, c7);
            AE_MULAF16X4SS(dr0, dr1, b, c8);
            v1 = AE_ROUND16X4F32SASYM(dr0, dr1);
            v1 = AE_ADD16S(v1, (ae_int16x4)128);
            AE_MINMAX16(v1, 0, 255);  

            /* save results */
            t = AE_SEL8X8I(AE_MOVINT8X8_FROMINT16X4(y0), AE_MOVINT8X8_FROMINT16X4(y1), 25); //Extract even
            AE_S8X8_IP(t, pOutY, sizeof(ae_int8x8));
            t = AE_SEL8X8I(AE_MOVINT8X8_FROMINT16X4(u0), AE_MOVINT8X8_FROMINT16X4(u1), 25); //Extract even
            AE_S8X8_IP(t, pOutU, sizeof(ae_int8x8));
            t = AE_SEL8X8I(AE_MOVINT8X8_FROMINT16X4(v0), AE_MOVINT8X8_FROMINT16X4(v1), 25); //Extract even
            AE_S8X8_IP(t, pOutV, sizeof(ae_int8x8));
        }
    }

    tmp = AE_L32_X((const ae_int32 *)c, 0 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c0 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 1 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c1 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 2 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c2 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 3 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c3 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 4 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c4 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 5 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c5 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 6 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c6 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 7 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c7 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 8 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c8 = AE_ROUND16X4F32SASYM(tmp, tmp);
    for (i = 0; i < h; i++)
    {
        for (j = w&(~7); j<w; j++)
        {
            ae_int16x4 y, u, v;
            r = AE_MOVDA16(inR[i*stride + j]);
            g = AE_MOVDA16(inG[i*stride + j]);
            b = AE_MOVDA16(inB[i*stride + j]);
            r = AE_SLAI16S(r, 7);
            g = AE_SLAI16S(g, 7);
            b = AE_SLAI16S(b, 7);

            AE_MUL16X4(dr0, dr1, r, c0);
            AE_MUL16X4(dg0, dg1, g, c1);
            AE_MUL16X4(db0, db1, b, c2);
            dr0 = AE_ADD32(dr0, dg0);
            db0 = AE_ADD32(db0, 1 << 21);
            dr0 = AE_ADD32(dr0, db0);
            dr1 = AE_ADD32(dr1, dg1);
            db1 = AE_ADD32(db1, 1 << 21);
            dr1 = AE_ADD32(dr1, db1);
            y = AE_SAT16X4(AE_SRAI32(dr0, 15), AE_SRAI32(dr1, 15));

            AE_MUL16X4(dr0, dr1, r, c3);
            AE_MUL16X4(dg0, dg1, g, c4);
            AE_MUL16X4(db0, db1, b, c5);
            dr0 = AE_ADD32(dr0, dg0);
            db0 = AE_ADD32(db0, 1 << 21);
            dr0 = AE_ADD32(dr0, db0);
            dr1 = AE_ADD32(dr1, dg1);
            db1 = AE_ADD32(db1, 1 << 21);
            dr1 = AE_ADD32(dr1, db1);
            u = AE_SAT16X4(AE_SRAI32(dr0, 15), AE_SRAI32(dr1, 15));

            AE_MUL16X4(dr0, dr1, r, c6);
            AE_MUL16X4(dg0, dg1, g, c7);
            AE_MUL16X4(db0, db1, b, c8);
            dr0 = AE_ADD32(dr0, dg0);
            db0 = AE_ADD32(db0, 1 << 21);
            dr0 = AE_ADD32(dr0, db0);
            dr1 = AE_ADD32(dr1, dg1);
            db1 = AE_ADD32(db1, 1 << 21);
            dr1 = AE_ADD32(dr1, db1);
            v = AE_SAT16X4(AE_SRAI32(dr0, 15), AE_SRAI32(dr1, 15));

            u = AE_ADD16S(u, (ae_int16x4)16384);
            v = AE_ADD16S(v, (ae_int16x4)16384);
            AE_MOVT16X4(v, 0, AE_LT16(v, 0));
            AE_MOVT16X4(u, 0, AE_LT16(u, 0));
            y = AE_SRAI16(y, 7);
            u = AE_SRAI16(u, 7);
            v = AE_SRAI16(v, 7);
            outY[i*stride + j] = (uint8_t)AE_MOVAD16_0(y);
            outU[i*stride + j] = (uint8_t)AE_MOVAD16_0(u);
            outV[i*stride + j] = (uint8_t)AE_MOVAD16_0(v);
        }
    }
#else
    const ae_int8x8 * restrict pInR;
    const ae_int8x8 * restrict pInG;
    const ae_int8x8 * restrict pInB;
    const ae_int8x16 * restrict pInR_;
    const ae_int8x16 * restrict pInG_;
    const ae_int8x16 * restrict pInB_;
    const uint8_t * restrict inR = (const uint8_t *)inImgR;    
    const uint8_t * restrict inG = (const uint8_t *)inImgG;   
    const uint8_t * restrict inB = (const uint8_t *)inImgB;
    ae_int8x8 * restrict pOutY;
    ae_int8x8 * restrict pOutU;
    ae_int8x8 * restrict pOutV;
    ae_int8x16 * restrict pOutY_;
    ae_int8x16 * restrict pOutU_;
    ae_int8x16 * restrict pOutV_;
    ae_valignx2 al_r, al_g, al_b, al_u, al_y, al_v;
    uint8_t * restrict outY = (uint8_t *)outImgY;  
    uint8_t * restrict outU = (uint8_t *)outImgU;
    uint8_t * restrict outV = (uint8_t *)outImgV;
    ae_int16x4 y0, u0, v0, y1, u1, v1;
    ae_int16x4 cy, cu, cv;
    //ae_int8x8 t;
    //ae_int32x2 dr0, dr1, dg0, dg1, db0, db1, 
    ae_int32x2 tmp;
    int i, j, w = (int)sz->width, h = (int)sz->height, stride = sz->stride;
    ae_int16x4 c0, c1, c2, c3, c4, c5, c6, c7, c8;        
    
    NASSERT(outImgY);
    NASSERT(outImgU);
    NASSERT(outImgV);
    NASSERT(inImgR);
    NASSERT(inImgG);
    NASSERT(inImgB);
    NASSERT_ALIGN(outImgY, ALIGNMENT);
    NASSERT_ALIGN(outImgU, ALIGNMENT);
    NASSERT_ALIGN(outImgV, ALIGNMENT);
    NASSERT_ALIGN(inImgR, ALIGNMENT);
    NASSERT_ALIGN(inImgG, ALIGNMENT);
    NASSERT_ALIGN(inImgB, ALIGNMENT);
    imgsize_validate(sz, 1, 1);

    tmp = AE_L32_X((const ae_int32 *)c, 0 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c0 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 1 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c1 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 2 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c2 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 3 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c3 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 4 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c4 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 5 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c5 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 6 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c6 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 7 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c7 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 8 * sizeof(int32_t));
    tmp = AE_SLLI32(tmp, 2);
    c8 = AE_ROUND16X4F32SASYM(tmp, tmp);

    cy = c0;
    AE_MOVT16X4(cy, c1, AE_MOVBA4(4));
    AE_MOVT16X4(cy, c2, AE_MOVBA4(2));
    AE_MOVT16X4(cy, 0, AE_MOVBA4(1));

    cu = c3;
    AE_MOVT16X4(cu, c4, AE_MOVBA4(4));
    AE_MOVT16X4(cu, c5, AE_MOVBA4(2));
    AE_MOVT16X4(cu, 0, AE_MOVBA4(1));

    cv = c6;
    AE_MOVT16X4(cv, c7, AE_MOVBA4(4));
    AE_MOVT16X4(cv, c8, AE_MOVBA4(2));
    AE_MOVT16X4(cv, 0, AE_MOVBA4(1));

    ae_int8x8 _0 = AE_MOVINT8X8_FROMINT64(0x0000000000000000);
    ae_int8x8 _128 = AE_MOVINT8X8_FROMINT64(0x8080808080808080);
    ae_int8x8 r, g, b, r_, g_, b_;
    ae_int8x8 _rg7, _rg5, _b07, _b05;

    ae_int8x8 sel0 = AE_MOVINT8X8_FROMINT64(0xfd75ec64b931a820);
    ae_int8x8 sel1 = AE_MOVINT8X8_FROMINT64(0xfdec7564b9a83120);
    ae_int8x8 rgb0, rgb1, rgb2, rgb3;
    ae_int32x2 y01, y23, y45, y67;
    ae_int32x2 u01, u23, u45, u67;
    ae_int32x2 v01, v23, v45, v67;
    ae_int32x2 _2 = (ae_int32x2)2;
    ae_int8x8 y, u, v;

    for (i = 0; i < h; i++)
    {       
        pInR = (const ae_int8x8 *)(inR + i*stride);
        pInG = (const ae_int8x8 *)(inG + i*stride);
        pInB = (const ae_int8x8 *)(inB + i*stride);
        pOutY = (ae_int8x8 *)(outY + i*stride);
        pOutU = (ae_int8x8 *)(outU + i*stride);
        pOutV = (ae_int8x8 *)(outV + i*stride);
        for (j = 0; j<(w>>3); j++)
        {            
            AE_L8X8_IP(r, pInR, sizeof(ae_int8x8));
            AE_L8X8_IP(g, pInG, sizeof(ae_int8x8));
            AE_L8X8_IP(b, pInB, sizeof(ae_int8x8));
            
            AE_DSEL8X8(_rg7, _rg5, r, g, sel0);
            AE_DSEL8X8(_b07, _b05, b, _0, sel0);
            AE_DSEL8X8(rgb1, rgb0, _rg5, _b05, sel1);
            AE_DSEL8X8(rgb3, rgb2, _rg7, _b07, sel1);

            AE_MULUS4O8X16(y01, y23, y45, y67, rgb3, rgb2, rgb1, rgb0, cy, cy);
            AE_MULUS4O8X16(u01, u23, u45, u67, rgb3, rgb2, rgb1, rgb0, cu, cu);
            AE_MULUS4O8X16(v01, v23, v45, v67, rgb3, rgb2, rgb1, rgb0, cv, cv);

            AE_MUL2P32X4S(y01, y23, y01, y23, _2, _2);
            y0 = AE_ROUND16X4F32SASYM(y01, y23);

            AE_MUL2P32X4S(y45, y67, y45, y67, _2, _2);
            y1 = AE_ROUND16X4F32SASYM(y45, y67);

            AE_MUL2P32X4S(u01, u23, u01, u23, _2, _2);
            u0 = AE_ROUND16X4F32SASYM(u01, u23);

            AE_MUL2P32X4S(u45, u67, u45, u67, _2, _2);
            u1 = AE_ROUND16X4F32SASYM(u45, u67);

            AE_MUL2P32X4S(v01, v23, v01, v23, _2, _2);
            v0 = AE_ROUND16X4F32SASYM(v01, v23);

            AE_MUL2P32X4S(v45, v67, v45, v67, _2, _2);
            v1 = AE_ROUND16X4F32SASYM(v45, v67);

            y = AE_SATU8X8X16(y0, y1);
            u = AE_SAT8X8X16(u0, u1);
            v = AE_SAT8X8X16(v0, v1);
            u = AE_ADD8(u, _128);           
            v = AE_ADD8(v, _128);           

            AE_S8X8_IP(y, pOutY, sizeof(ae_int8x8));
            AE_S8X8_IP(u, pOutU, sizeof(ae_int8x8));
            AE_S8X8_IP(v, pOutV, sizeof(ae_int8x8));
        }
        int off = w&7;
        if (off)
        {
            pInR_ = (const ae_int8x16 *)(pInR);
            pInG_ = (const ae_int8x16 *)(pInG);
            pInB_ = (const ae_int8x16 *)(pInB);
            pOutY_ = (ae_int8x16 *)(pOutY);
            pOutU_ = (ae_int8x16 *)(pOutU);
            pOutV_ = (ae_int8x16 *)(pOutV);
            al_r = AE_LA128_PP(pInR_);
            al_g = AE_LA128_PP(pInG_);
            al_b = AE_LA128_PP(pInB_);
            AE_LAV8X8X2_XP(r, r_, al_r, pInR_, off);
            AE_LAV8X8X2_XP(g, g_, al_g, pInG_, off);
            AE_LAV8X8X2_XP(b, b_, al_b, pInB_, off);

            AE_DSEL8X8(_rg7, _rg5, r, g, sel0);
            AE_DSEL8X8(_b07, _b05, b, _0, sel0);
            AE_DSEL8X8(rgb1, rgb0, _rg5, _b05, sel1);
            AE_DSEL8X8(rgb3, rgb2, _rg7, _b07, sel1);

            AE_MULUS4O8X16(y01, y23, y45, y67, rgb3, rgb2, rgb1, rgb0, cy, cy);
            AE_MULUS4O8X16(u01, u23, u45, u67, rgb3, rgb2, rgb1, rgb0, cu, cu);
            AE_MULUS4O8X16(v01, v23, v45, v67, rgb3, rgb2, rgb1, rgb0, cv, cv);

            AE_MUL2P32X4S(y01, y23, y01, y23, _2, _2);
            y0 = AE_ROUND16X4F32SASYM(y01, y23);

            AE_MUL2P32X4S(y45, y67, y45, y67, _2, _2);
            y1 = AE_ROUND16X4F32SASYM(y45, y67);

            AE_MUL2P32X4S(u01, u23, u01, u23, _2, _2);
            u0 = AE_ROUND16X4F32SASYM(u01, u23);

            AE_MUL2P32X4S(u45, u67, u45, u67, _2, _2);
            u1 = AE_ROUND16X4F32SASYM(u45, u67);

            AE_MUL2P32X4S(v01, v23, v01, v23, _2, _2);
            v0 = AE_ROUND16X4F32SASYM(v01, v23);

            AE_MUL2P32X4S(v45, v67, v45, v67, _2, _2);
            v1 = AE_ROUND16X4F32SASYM(v45, v67);

            y = AE_SATU8X8X16(y0, y1);
            u = AE_SAT8X8X16(u0, u1);
            v = AE_SAT8X8X16(v0, v1);
            u = AE_ADD8(u, _128);
            v = AE_ADD8(v, _128);        

            al_u = al_y = al_v = AE_ZALIGN128();
            AE_SAV8X8X2_XP(y, y, al_y, pOutY_, off);
            AE_SAV8X8X2_XP(u, u, al_u, pOutU_, off);
            AE_SAV8X8X2_XP(v, v, al_v, pOutV_, off);

            AE_SA128POS_FP(al_y, pOutY_);
            AE_SA128POS_FP(al_u, pOutU_);
            AE_SA128POS_FP(al_v, pOutV_);
        }
    }

#endif
}
