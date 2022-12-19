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
  Functions convert YUV (BT.601) planar images to RGB planar 

  R=Y+c9*(V-bias)
  G=Y+c10*(U-bias)+c11*(V-bias)
  B=Y+c12*(U-bias)
  where bias is 128 for 8-bit data and 16384 for 16-bit data

  Image formats:
  8-bit unsigned RGB/YUV
  16-bit signed RGB/YUV

  Input:
  inImgY,
  inImgU,
  inImgV  planes with Y,U,V components
  c[13]   transform coefficients, Q29 format
  sz      image size
  Output:
  outImgR
  outImgG
  outImgB planes with R,G,B components

  Restrictions:
  see general restrictions applied for all images for fast/generic 
  functions
-------------------------------------------------------------------------*/
void imgconvert_yuvrgb( void * restrict outImgR, 
                        void * restrict outImgG, 
                        void * restrict outImgB, 
                        const void * restrict inImgY, 
                        const void * restrict inImgU, 
                        const void * restrict inImgV, 
                        const int32_t * restrict c,
                        const imgsize_t* sz)
{
    int i, j, w = (int)sz->width, h = (int)sz->height, stride = sz->stride;
    ae_int16x4 y, u, v, r0, g0, b0, r1, g1, b1;
    ae_int32x2 tmp;
    ae_int32x2 dr0, dr1, dg0, dg1, db0, db1;
    ae_int8x8 t;
    const signed char * restrict pInY;
    const uint8_t * restrict inY = (const uint8_t *)inImgY;
    const signed char * restrict pInU;
    const uint8_t * restrict inU = (const uint8_t *)inImgU;
    const signed char * restrict pInV;
    const uint8_t * restrict inV = (const uint8_t *)inImgV;
    ae_int8x8 * restrict pOutR;
    uint8_t * restrict outR = (uint8_t *)outImgR;
    ae_int8x8 * restrict pOutG;
    uint8_t * restrict outG = (uint8_t *)outImgG;
    ae_int8x8 * restrict pOutB;
    uint8_t * restrict outB = (uint8_t *)outImgB;
    ae_int16x4 c9, c10, c11, c12;
    ae_valign aly, alu, alv, alr, alg, alb;
    NASSERT(outImgR);
    NASSERT(outImgG);
    NASSERT(outImgB);
    NASSERT(inImgY);
    NASSERT(inImgU);
    NASSERT(inImgV);
    NASSERT_ALIGN(outImgR, 1);
    NASSERT_ALIGN(outImgG, 1);
    NASSERT_ALIGN(outImgB, 1);
    NASSERT_ALIGN(inImgY, 1);
    NASSERT_ALIGN(inImgU, 1);
    NASSERT_ALIGN(inImgV, 1);
    imgsize_validate(sz, 1, 0);

    /*R*/
    tmp = AE_L32_X((const ae_int32 *)c, 9 * sizeof(int32_t));
    c9 = AE_ROUND16X4F32SASYM(tmp, tmp);
    alr = AE_ZALIGN64();
    for (i = 0; i < h; i++)
    {
        pInY = (const signed char*)(inY + i*stride);
        pInV = (const signed char*)(inV + i*stride);
        pOutR = (ae_int8x8 *)(outR + i*stride);
        aly = AE_LA64_PP(pInY);
        alv = AE_LA64_PP(pInV);
        for (j = 0; j < (w>>3); j++)
        {          
            AE_LA8X4U_IP(y, aly, pInY);
            AE_LA8X4U_IP(v, alv, pInV);

            v = AE_SUB16(v, 128);

            AE_MULFD16X16X4WS(dr0, dr1, (ae_int16x4)8192, c9, y, v);;
            dr0 = AE_SLLI32S(dr0, 2);
            dr1 = AE_SLLI32S(dr1, 2);
            r0 = AE_ROUND16X4F32SASYM(dr0, dr1);

            AE_MINMAX16(r0, 0, 255);

            ///*Second half*/

            AE_LA8X4U_IP(y, aly, pInY);
            AE_LA8X4U_IP(v, alv, pInV);
            v = AE_SUB16(v, 128);
            u = AE_SUB16(u, 128);

            AE_MULFD16X16X4WS(dr0, dr1, (ae_int16x4)8192, c9, y, v);
            dr0 = AE_SLLI32S(dr0, 2);
            dr1 = AE_SLLI32S(dr1, 2);
            r1 = AE_ROUND16X4F32SASYM(dr0, dr1);

            AE_MINMAX16(r1, 0, 255);

            /* save results */          
            t = AE_SEL8X8I(AE_MOVINT8X8_FROMINT16X4(r0), AE_MOVINT8X8_FROMINT16X4(r1), 25); //Extract even
            AE_SA8X8_IP(t, alr, pOutR);        
        }
        AE_SA64POS_FP(alr, pOutR);
    }
    /*G*/
    tmp = AE_L32_X((const ae_int32 *)c, 10 * sizeof(int32_t));
    c10 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 11 * sizeof(int32_t));
    c11 = AE_ROUND16X4F32SASYM(tmp, tmp);
    alg = AE_ZALIGN64();
    for (i = 0; i < h; i++)
    {
        pInY = (const signed char*)(inY + i*stride);
        pInU = (const signed char*)(inU + i*stride);
        pInV = (const signed char*)(inV + i*stride);
        pOutG = (ae_int8x8 *)(outG + i*stride);
        aly = AE_LA64_PP(pInY);
        alu = AE_LA64_PP(pInU);
        alv = AE_LA64_PP(pInV);
        for (j = 0; j < (w>>3); j++)
        {          
            AE_LA8X4U_IP(y, aly, pInY);
            AE_LA8X4U_IP(u, alu, pInU);
            AE_LA8X4U_IP(v, alv, pInV);

            v = AE_SUB16(v, 128);
            u = AE_SUB16(u, 128);

            AE_MULFD16X16X4WS(dg0, dg1, (ae_int16x4)8192, c10, y, u);
            AE_MULAF16X4SS(dg0, dg1, c11, v);
            dg0 = AE_SLLI32S(dg0, 2);
            dg1 = AE_SLLI32S(dg1, 2);
            g0 = AE_ROUND16X4F32SASYM(dg0, dg1);

            AE_MINMAX16(g0, 0, 255);

            ///*Second half*/

            AE_LA8X4U_IP(y, aly, pInY);
            AE_LA8X4U_IP(u, alu, pInU);
            AE_LA8X4U_IP(v, alv, pInV);
            v = AE_SUB16(v, 128);
            u = AE_SUB16(u, 128);

            AE_MULFD16X16X4WS(dg0, dg1, (ae_int16x4)8192, c10, y, u);
            AE_MULAF16X4SS(dg0, dg1, c11, v);
            dg0 = AE_SLLI32S(dg0, 2);
            dg1 = AE_SLLI32S(dg1, 2);
            g1 = AE_ROUND16X4F32SASYM(dg0, dg1);

            AE_MINMAX16(g1, 0, 255);

            /* save results */          
            t = AE_SEL8X8I(AE_MOVINT8X8_FROMINT16X4(g0), AE_MOVINT8X8_FROMINT16X4(g1), 25); //Extract even
            AE_SA8X8_IP(t, alg, pOutG);         
        }
        AE_SA64POS_FP(alg, pOutG);
    }
    /*B*/
    tmp = AE_L32_X((const ae_int32 *)c, 12 * sizeof(int32_t));
    c12 = AE_ROUND16X4F32SASYM(tmp, tmp);
    alb = AE_ZALIGN64();
    for (i = 0; i < h; i++)
    {
        pInY = (const signed char*)(inY + i*stride);
        pInU = (const signed char*)(inU + i*stride);
        pInV = (const signed char*)(inV + i*stride);
        pOutB = (ae_int8x8 *)(outB + i*stride);
        aly = AE_LA64_PP(pInY);
        alu = AE_LA64_PP(pInU);
        for (j = 0; j < (w>>3); j++)
        {          
            AE_LA8X4U_IP(y, aly, pInY);
            AE_LA8X4U_IP(u, alu, pInU);
            u = AE_SUB16(u, 128);

            AE_MULFD16X16X4WS(db0, db1, (ae_int16x4)8192, c12, y, u);
            db0 = AE_SLLI32S(db0, 2);
            db1 = AE_SLLI32S(db1, 2);
            b0 = AE_ROUND16X4F32SASYM(db0, db1);

            AE_MINMAX16(b0, 0, 255);

            ///*Second half*/

            AE_LA8X4U_IP(y, aly, pInY);
            AE_LA8X4U_IP(u, alu, pInU);
            u = AE_SUB16(u, 128);

            AE_MULFD16X16X4WS(db0, db1, (ae_int16x4)8192, c12, y, u);
            db0 = AE_SLLI32S(db0, 2);
            db1 = AE_SLLI32S(db1, 2);
            b1 = AE_ROUND16X4F32SASYM(db0, db1);

            AE_MINMAX16(b1, 0, 255);

            /* save results */          
            t = AE_SEL8X8I(AE_MOVINT8X8_FROMINT16X4(b0), AE_MOVINT8X8_FROMINT16X4(b1), 25); //Extract even
            AE_SA8X8_IP(t, alb, pOutB);            
        }
        AE_SA64POS_FP(alb, pOutB);
    }

    tmp = AE_L32_X((const ae_int32 *)c, 9 * sizeof(int32_t));
    c9 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 10 * sizeof(int32_t));
    c10 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 11 * sizeof(int32_t));
    c11 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 12 * sizeof(int32_t));
    c12 = AE_ROUND16X4F32SASYM(tmp, tmp);

    for (i = 0; i < h; i++)
    {
        for (j = w&(~7); j < w; j++)
        {
            ae_int16x4 r, g, b;
            y = AE_MOVDA16(((const uint8_t*)inImgY)[i*stride + j]);
            u = AE_MOVDA16(((const uint8_t*)inImgU)[i*stride + j]);
            v = AE_MOVDA16(((const uint8_t*)inImgV)[i*stride + j]);
            v = AE_SUB16(v, 128);
            u = AE_SUB16(u, 128);

            AE_MULFD16X16X4WS(dr0, dr1, (ae_int16x4)8192, c9, y, v);;
            dr0 = AE_SLLI32S(dr0, 2);
            dr1 = AE_SLLI32S(dr1, 2);
            r = AE_ROUND16X4F32SASYM(dr0, dr1);

            AE_MULFD16X16X4WS(db0, db1, (ae_int16x4)8192, c12, y, u);
            db0 = AE_SLLI32S(db0, 2);
            db1 = AE_SLLI32S(db1, 2);
            b = AE_ROUND16X4F32SASYM(db0, db1);

            AE_MULFD16X16X4WS(dg0, dg1, (ae_int16x4)8192, c10, y, u);
            AE_MULAF16X4SS(dg0, dg1, c11, v);
            dg0 = AE_SLLI32S(dg0, 2);
            dg1 = AE_SLLI32S(dg1, 2);
            g = AE_ROUND16X4F32SASYM(dg0, dg1);

            AE_MINMAX16(r, 0, 255);
            AE_MINMAX16(g, 0, 255);
            AE_MINMAX16(b, 0, 255);
            ((uint8_t*)outImgR)[i*stride + j] = (uint8_t)AE_MOVAD16_0(r);
            ((uint8_t*)outImgG)[i*stride + j] = (uint8_t)AE_MOVAD16_0(g);
            ((uint8_t*)outImgB)[i*stride + j] = (uint8_t)AE_MOVAD16_0(b);
        }
    }
}
