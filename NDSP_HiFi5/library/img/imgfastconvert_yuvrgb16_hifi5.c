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
void imgfastconvert_yuvrgb16( 
                        void * restrict outImgR, 
                        void * restrict outImgG, 
                        void * restrict outImgB, 
                        const void * restrict inImgY, 
                        const void * restrict inImgU, 
                        const void * restrict inImgV, 
                        const int32_t * restrict c,
                        const imgsize_t* sz)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
    int i, j, w = (int)sz->width, h = (int)sz->height, stride = sz->stride;
    ae_int16x4 y, u, v;
    ae_int32x2 dy0, dy1, dv0, dv1, du0, du1, tmp;
    ae_int16x4 r, g, b;
    ae_int32x2 dr0, dr1, dg0, dg1, db0, db1;
    const ae_int16x4 * restrict pInY;
    const int16_t * restrict inY = (const int16_t *)inImgY;
    const ae_int16x4 * restrict pInU;
    const int16_t * restrict inU = (const int16_t *)inImgU;
    const ae_int16x4 * restrict pInV;
    const int16_t * restrict inV = (const int16_t *)inImgV;
    ae_int16x4 * restrict pOutR;
    int16_t * restrict outR = (int16_t *)outImgR;
    ae_int16x4 * restrict pOutG;
    int16_t * restrict outG = (int16_t *)outImgG;
    ae_int16x4 * restrict pOutB;
    int16_t * restrict outB = (int16_t *)outImgB;
    ae_int16x4 c9, c10, c11, c12;
    NASSERT(outImgR);
    NASSERT(outImgG);
    NASSERT(outImgB);
    NASSERT(inImgY);
    NASSERT(inImgU);
    NASSERT(inImgV);
    NASSERT_ALIGN(outImgR, ALIGNMENT);
    NASSERT_ALIGN(outImgG, ALIGNMENT);
    NASSERT_ALIGN(outImgB, ALIGNMENT);
    NASSERT_ALIGN(inImgY, ALIGNMENT);
    NASSERT_ALIGN(inImgU, ALIGNMENT);
    NASSERT_ALIGN(inImgV, ALIGNMENT);
    imgsize_validate(sz, 2, 1);
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
        pInY = (const ae_int16x4*)(inY + i*stride);
        pInU = (const ae_int16x4*)(inU + i*stride);
        pInV = (const ae_int16x4*)(inV + i*stride);
        pOutR = (ae_int16x4 *)(outR + i*stride);
        pOutG = (ae_int16x4 *)(outG + i*stride);
        pOutB = (ae_int16x4 *)(outB + i*stride);
        for (j = 0; j < (w>>2); j++)
        {           
            AE_L16X4_IP(y, pInY, 4 * sizeof(int16_t));
            AE_L16X4_IP(u, pInU, 4 * sizeof(int16_t));
            AE_L16X4_IP(v, pInV, 4 * sizeof(int16_t));

            y = AE_MAX16(0, y);
            u = AE_MAX16(0, u);
            v = AE_MAX16(0, v);

            v = AE_SUB16(v, 16384);
            u = AE_SUB16(u, 16384);

            AE_MULFD16X16X4WS(dr0, dr1, (ae_int16x4)8192, c9, y, v);
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

            r = AE_MAX16(0, r);
            g = AE_MAX16(0, g);
            b = AE_MAX16(0, b);
            AE_S16X4_IP(r, pOutR, 4 * sizeof(int16_t));
            AE_S16X4_IP(g, pOutG, 4 * sizeof(int16_t));
            AE_S16X4_IP(b, pOutB, 4 * sizeof(int16_t));
        }
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
        for (j = w&(~3); j < w; j++)
        {
            ae_int16x4 r, g, b;
            y = AE_MOVDA16(inY[i*stride + j]);
            u = AE_MOVDA16(inU[i*stride + j]);
            v = AE_MOVDA16(inV[i*stride + j]);
            AE_MOVT16X4(y, 0, AE_LT16(y, 0));
            AE_MOVT16X4(u, 0, AE_LT16(u, 0));
            AE_MOVT16X4(v, 0, AE_LT16(v, 0));
            v = AE_SUB16(v, 16384);
            u = AE_SUB16(u, 16384);
            AE_MUL16X4(dy0, dy1, 8192, y);
            AE_MUL16X4(dv0, dv1, c9, v);
            dv0 = AE_ADD32(dv0, dy0);
            dv1 = AE_ADD32(dv1, dy1);
            r = AE_SAT16X4(AE_SRAA32RS(dv0, 13), AE_SRAA32RS(dv1, 13));
            AE_MUL16X4(du0, du1, c12, u);
            du0 = AE_ADD32(du0, dy0);
            du1 = AE_ADD32(du1, dy1);
            b = AE_SAT16X4(AE_SRAA32RS(du0, 13), AE_SRAA32RS(du1, 13));
            AE_MUL16X4(du0, du1, c10, u);
            AE_MUL16X4(dv0, dv1, c11, v);
            du0 = AE_ADD32(du0, dy0);
            dv0 = AE_ADD32(dv0, du0);
            du1 = AE_ADD32(du1, dy1);
            dv1 = AE_ADD32(dv1, du1);
            g = AE_SAT16X4(AE_SRAA32RS(dv0, 13), AE_SRAA32RS(dv1, 13));

            AE_MOVT16X4(r, 0, AE_LT16(r, 0));
            AE_MOVT16X4(g, 0, AE_LT16(g, 0));
            AE_MOVT16X4(b, 0, AE_LT16(b, 0));
            outR[i*stride + j] = AE_MOVAD16_0(r);
            outG[i*stride + j] = AE_MOVAD16_0(g);
            outB[i*stride + j] = AE_MOVAD16_0(b);
        }
    }
#else
    int i, j, w = (int)sz->width, h = (int)sz->height, stride = sz->stride;
    ae_int16x4 y0, u0, v0;
    ae_int16x4 y1, u1, v1;
    ae_int32x2 tmp;
    ae_int16x4 r0, g0, b0;
    ae_int16x4 r1, g1, b1;
    ae_valignx2 aly, alu, alv, alr, alg, alb;
    ae_int32x2 dr0, dr1, dg0, dg1, db0, db1;
    const ae_int16x8 * restrict pInY;
    const int16_t * restrict inY = (const int16_t *)inImgY;
    const ae_int16x8 * restrict pInU;
    const int16_t * restrict inU = (const int16_t *)inImgU;
    const ae_int16x8 * restrict pInV;
    const int16_t * restrict inV = (const int16_t *)inImgV;
    ae_int16x8 * restrict pOutR;
    int16_t * restrict outR = (int16_t *)outImgR;
    ae_int16x8 * restrict pOutG;
    int16_t * restrict outG = (int16_t *)outImgG;
    ae_int16x8 * restrict pOutB;
    int16_t * restrict outB = (int16_t *)outImgB;
    ae_int16x4 c9, c10, c11, c12;
    NASSERT(outImgR);
    NASSERT(outImgG);
    NASSERT(outImgB);
    NASSERT(inImgY);
    NASSERT(inImgU);
    NASSERT(inImgV);
    NASSERT_ALIGN(outImgR, ALIGNMENT);
    NASSERT_ALIGN(outImgG, ALIGNMENT);
    NASSERT_ALIGN(outImgB, ALIGNMENT);
    NASSERT_ALIGN(inImgY, ALIGNMENT);
    NASSERT_ALIGN(inImgU, ALIGNMENT);
    NASSERT_ALIGN(inImgV, ALIGNMENT);
    imgsize_validate(sz, 2, 1);
    tmp = AE_L32_X((const ae_int32 *)c, 9 * sizeof(int32_t));
    c9 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 10 * sizeof(int32_t));
    c10 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 11 * sizeof(int32_t));
    c11 = AE_ROUND16X4F32SASYM(tmp, tmp);
    tmp = AE_L32_X((const ae_int32 *)c, 12 * sizeof(int32_t));
    c12 = AE_ROUND16X4F32SASYM(tmp, tmp);
    alr = AE_ZALIGN128();
    alg = AE_ZALIGN128();
    alb = AE_ZALIGN128();
    for (i = 0; i < h; i++)
    {
        pInY = (const ae_int16x8*)(inY + i*stride);
        pInU = (const ae_int16x8*)(inU + i*stride);
        pInV = (const ae_int16x8*)(inV + i*stride);
        pOutR = (ae_int16x8 *)(outR + i*stride);
        pOutG = (ae_int16x8 *)(outG + i*stride);
        pOutB = (ae_int16x8 *)(outB + i*stride);
        for (j = 0; j < (w>>3); j++)
        {           
            AE_L16X4X2_IP(y0, y1, pInY, 8 * sizeof(int16_t));
            AE_L16X4X2_IP(u0, u1, pInU, 8 * sizeof(int16_t));
            AE_L16X4X2_IP(v0, v1, pInV, 8 * sizeof(int16_t));

            y0 = AE_MAX16(0, y0);
            u0 = AE_MAX16(0, u0);
            v0 = AE_MAX16(0, v0);
            v0 = AE_SUB16(v0, 16384);
            u0 = AE_SUB16(u0, 16384);
            AE_MULFD16X16X4WS(dr0, dr1, (ae_int16x4)8192, c9, y0, v0);
            dr0 = AE_SLLI32S(dr0, 2);
            dr1 = AE_SLLI32S(dr1, 2);
            r0 = AE_ROUND16X4F32SASYM(dr0, dr1);
            AE_MULFD16X16X4WS(db0, db1, (ae_int16x4)8192, c12, y0, u0);
            db0 = AE_SLLI32S(db0, 2);
            db1 = AE_SLLI32S(db1, 2);
            b0 = AE_ROUND16X4F32SASYM(db0, db1);
            AE_MULFD16X16X4WS(dg0, dg1, (ae_int16x4)8192, c10, y0, u0);
            AE_MULAF16X4SS(dg0, dg1, c11, v0);
            dg0 = AE_SLLI32S(dg0, 2);
            dg1 = AE_SLLI32S(dg1, 2);
            g0 = AE_ROUND16X4F32SASYM(dg0, dg1);
            r0 = AE_MAX16(0, r0);
            g0 = AE_MAX16(0, g0);
            b0 = AE_MAX16(0, b0);

            y1 = AE_MAX16(0, y1);
            u1 = AE_MAX16(0, u1);
            v1 = AE_MAX16(0, v1);
            v1 = AE_SUB16(v1, 16384);
            u1 = AE_SUB16(u1, 16384);
            AE_MULFD16X16X4WS(dr0, dr1, (ae_int16x4)8192, c9, y1, v1);
            dr0 = AE_SLLI32S(dr0, 2);
            dr1 = AE_SLLI32S(dr1, 2);
            r1 = AE_ROUND16X4F32SASYM(dr0, dr1);
            AE_MULFD16X16X4WS(db0, db1, (ae_int16x4)8192, c12, y1, u1);
            db0 = AE_SLLI32S(db0, 2);
            db1 = AE_SLLI32S(db1, 2);
            b1 = AE_ROUND16X4F32SASYM(db0, db1);
            AE_MULFD16X16X4WS(dg0, dg1, (ae_int16x4)8192, c10, y1, u1);
            AE_MULAF16X4SS(dg0, dg1, c11, v1);
            dg0 = AE_SLLI32S(dg0, 2);
            dg1 = AE_SLLI32S(dg1, 2);
            g1 = AE_ROUND16X4F32SASYM(dg0, dg1);
            r1 = AE_MAX16(0, r1);
            g1 = AE_MAX16(0, g1);
            b1 = AE_MAX16(0, b1);

            AE_S16X4X2_IP(r0, r1, pOutR, 8 * sizeof(int16_t));
            AE_S16X4X2_IP(g0, g1, pOutG, 8 * sizeof(int16_t));
            AE_S16X4X2_IP(b0, b1, pOutB, 8 * sizeof(int16_t));
        }
        if (w & 7)
        {     
            int off = (w & 7) <<1;   
            aly = AE_LA128_PP(pInY);
            alu = AE_LA128_PP(pInU);
            alv = AE_LA128_PP(pInV);   
            AE_LAV16X4X2_XP(y0, y1, aly, pInY, off);
            AE_LAV16X4X2_XP(u0, u1, alu, pInU, off);
            AE_LAV16X4X2_XP(v0, v1, alv, pInV, off);

            y0 = AE_MAX16(0, y0);
            u0 = AE_MAX16(0, u0);
            v0 = AE_MAX16(0, v0);
            v0 = AE_SUB16(v0, 16384);
            u0 = AE_SUB16(u0, 16384);
            AE_MULFD16X16X4WS(dr0, dr1, (ae_int16x4)8192, c9, y0, v0);
            dr0 = AE_SLLI32S(dr0, 2);
            dr1 = AE_SLLI32S(dr1, 2);
            r0 = AE_ROUND16X4F32SASYM(dr0, dr1);
            AE_MULFD16X16X4WS(db0, db1, (ae_int16x4)8192, c12, y0, u0);
            db0 = AE_SLLI32S(db0, 2);
            db1 = AE_SLLI32S(db1, 2);
            b0 = AE_ROUND16X4F32SASYM(db0, db1);
            AE_MULFD16X16X4WS(dg0, dg1, (ae_int16x4)8192, c10, y0, u0);
            AE_MULAF16X4SS(dg0, dg1, c11, v0);
            dg0 = AE_SLLI32S(dg0, 2);
            dg1 = AE_SLLI32S(dg1, 2);
            g0 = AE_ROUND16X4F32SASYM(dg0, dg1);
            r0 = AE_MAX16(0, r0);
            g0 = AE_MAX16(0, g0);
            b0 = AE_MAX16(0, b0);

            y1 = AE_MAX16(0, y1);
            u1 = AE_MAX16(0, u1);
            v1 = AE_MAX16(0, v1);
            v1 = AE_SUB16(v1, 16384);
            u1 = AE_SUB16(u1, 16384);
            AE_MULFD16X16X4WS(dr0, dr1, (ae_int16x4)8192, c9, y1, v1);
            dr0 = AE_SLLI32S(dr0, 2);
            dr1 = AE_SLLI32S(dr1, 2);
            r1 = AE_ROUND16X4F32SASYM(dr0, dr1);
            AE_MULFD16X16X4WS(db0, db1, (ae_int16x4)8192, c12, y1, u1);
            db0 = AE_SLLI32S(db0, 2);
            db1 = AE_SLLI32S(db1, 2);
            b1 = AE_ROUND16X4F32SASYM(db0, db1);
            AE_MULFD16X16X4WS(dg0, dg1, (ae_int16x4)8192, c10, y1, u1);
            AE_MULAF16X4SS(dg0, dg1, c11, v1);
            dg0 = AE_SLLI32S(dg0, 2);
            dg1 = AE_SLLI32S(dg1, 2);
            g1 = AE_ROUND16X4F32SASYM(dg0, dg1);
            r1 = AE_MAX16(0, r1);
            g1 = AE_MAX16(0, g1);
            b1 = AE_MAX16(0, b1);

            AE_SAV16X4X2_XP(r0, r1, alr, pOutR, off);
            AE_SAV16X4X2_XP(g0, g1, alg, pOutG, off);
            AE_SAV16X4X2_XP(b0, b1, alb, pOutB, off);
            AE_SA128POS_FP(alr, pOutR);
            AE_SA128POS_FP(alg, pOutG);
            AE_SA128POS_FP(alb, pOutB);
        }       
    }
#endif
}
