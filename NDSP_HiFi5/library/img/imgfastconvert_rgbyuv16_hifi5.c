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
void imgfastconvert_rgbyuv16( 
                        void * restrict outImgY, 
                        void * restrict outImgU, 
                        void * restrict outImgV, 
                        const void * restrict inImgR, 
                        const void * restrict inImgG, 
                        const void * restrict inImgB, 
                        const int32_t * restrict c,
                        const imgsize_t* sz)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
    const ae_int16x4 * restrict pInR;
    const int16_t * restrict inR = (const int16_t *)inImgR;
    const ae_int16x4 * restrict pInG;
    const int16_t * restrict inG = (const int16_t *)inImgG;
    const ae_int16x4 * restrict pInB;
    const int16_t * restrict inB = (const int16_t *)inImgB;
    ae_int16x4 * restrict pOutY;
    int16_t * restrict outY = (int16_t *)outImgY;
    ae_int16x4 * restrict pOutU;
    int16_t * restrict outU = (int16_t *)outImgU;
    ae_int16x4 * restrict pOutV;
    int16_t * restrict outV = (int16_t *)outImgV;
    ae_int16x4 r, g, b, y, u, v;
    ae_int32x2 dr0, dr1, dg0, dg1, db0, db1, tmp;
    ae_int32x2 dy0, dy1, du0, du1, dv0, dv1;
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
    imgsize_validate(sz, 2, 1);
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
        pInR = (const ae_int16x4*)(inR + i*stride);
        pInG = (const ae_int16x4*)(inG + i*stride);
        pInB = (const ae_int16x4*)(inB + i*stride);
        pOutY = (ae_int16x4 *)(outY + i*stride);
        pOutU = (ae_int16x4 *)(outU + i*stride);
        pOutV = (ae_int16x4 *)(outV + i*stride);
        for (j = 0; j < (w >> 2); j++)
        {           
            AE_L16X4_IP(r, pInR, 4 * sizeof(int16_t));
            AE_L16X4_IP(g, pInG, 4 * sizeof(int16_t));
            AE_L16X4_IP(b, pInB, 4 * sizeof(int16_t));

            r = AE_MAX16(0, r);
            g = AE_MAX16(0, g);
            b = AE_MAX16(0, b);

            AE_MULFD16X16X4WS(dy0, dy1, r, g, c0, c1);
            AE_MULAF16X4SS(dy0, dy1, b, c2);

            y = AE_ROUND16X4F32SASYM(dy0, dy1);
            AE_S16X4_IP(y, pOutY, sizeof(ae_int16x4));

            AE_MULFD16X16X4WS(du0, du1, r, g, c3, c4);
            AE_MULAF16X4SS(du0, du1, b, c5);

            u = AE_ROUND16X4F32SASYM(du0, du1);

            u = AE_ADD16S(u, (ae_int16x4)16384);;
            u = AE_MAX16(0, u);

            AE_S16X4_IP(u, pOutU, sizeof(ae_int16x4));

            AE_MULFD16X16X4WS(dv0, dv1, r, g, c6, c7);
            AE_MULAF16X4SS(dv0, dv1, b, c8);

            v = AE_ROUND16X4F32SASYM(dv0, dv1);

            v = AE_ADD16S(v, (ae_int16x4)16384);
            v = AE_MAX16(0, v);

            AE_S16X4_IP(v, pOutV, sizeof(ae_int16x4));
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
        for (j = w&(~3); j<w; j++)
        {
            ae_int16x4 y, u, v;
            r = AE_MOVDA16(inR[i*stride + j]);
            g = AE_MOVDA16(inG[i*stride + j]);
            b = AE_MOVDA16(inB[i*stride + j]);
            AE_MOVT16X4(r, 0, AE_LT16(r, 0));
            AE_MOVT16X4(g, 0, AE_LT16(g, 0));
            AE_MOVT16X4(b, 0, AE_LT16(b, 0));

            AE_MUL16X4(dr0, dr1, r, c0);
            AE_MUL16X4(dg0, dg1, g, c1);
            AE_MUL16X4(db0, db1, b, c2);
            dr0 = AE_ADD32(dr0, dg0);
            dr0 = AE_ADD32(dr0, db0);
            dr1 = AE_ADD32(dr1, dg1);
            dr1 = AE_ADD32(dr1, db1);
            y = AE_SAT16X4(AE_SRAA32RS(dr0, 15), AE_SRAA32RS(dr1, 15));

            AE_MUL16X4(dr0, dr1, r, c3);
            AE_MUL16X4(dg0, dg1, g, c4);
            AE_MUL16X4(db0, db1, b, c5);
            dr0 = AE_ADD32(dr0, dg0);
            dr0 = AE_ADD32(dr0, db0);
            dr1 = AE_ADD32(dr1, dg1);
            dr1 = AE_ADD32(dr1, db1);
            u = AE_SAT16X4(AE_SRAA32RS(dr0, 15), AE_SRAA32RS(dr1, 15));

            AE_MUL16X4(dr0, dr1, r, c6);
            AE_MUL16X4(dg0, dg1, g, c7);
            AE_MUL16X4(db0, db1, b, c8);
            dr0 = AE_ADD32(dr0, dg0);
            dr0 = AE_ADD32(dr0, db0);
            dr1 = AE_ADD32(dr1, dg1);
            dr1 = AE_ADD32(dr1, db1);
            v = AE_SAT16X4(AE_SRAA32RS(dr0, 15), AE_SRAA32RS(dr1, 15));

            u = AE_ADD16S(u, (ae_int16x4)16384);
            v = AE_ADD16S(v, (ae_int16x4)16384);
            AE_MOVT16X4(v, 0, AE_LT16(v, 0));
            AE_MOVT16X4(u, 0, AE_LT16(u, 0));
            outY[i*stride + j] = AE_MOVAD16_0(y);
            outU[i*stride + j] = AE_MOVAD16_0(u);
            outV[i*stride + j] = AE_MOVAD16_0(v);
        }
    }
#else
    const ae_int16x8 * restrict pInR;
    const int16_t * restrict inR = (const int16_t *)inImgR;
    const ae_int16x8 * restrict pInG;
    const int16_t * restrict inG = (const int16_t *)inImgG;
    const ae_int16x8 * restrict pInB;
    const int16_t * restrict inB = (const int16_t *)inImgB;
    ae_int16x8 * restrict pOutY;
    int16_t * restrict outY = (int16_t *)outImgY;
    ae_int16x8 * restrict pOutU;
    int16_t * restrict outU = (int16_t *)outImgU;
    ae_int16x8 * restrict pOutV;
    int16_t * restrict outV = (int16_t *)outImgV;
    ae_int16x4 r0, g0, b0, y0, u0, v0;
    ae_int16x4 r1, g1, b1, y1, u1, v1;
    ae_valignx2 aly, alu, alv, alr, alg, alb;
    ae_int32x2 tmp;
    ae_int32x2 dy0, dy1, du0, du1, dv0, dv1;
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
    imgsize_validate(sz, 2, 1);
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
        pInR = (const ae_int16x8*)(inR + i*stride);
        pInG = (const ae_int16x8*)(inG + i*stride);
        pInB = (const ae_int16x8*)(inB + i*stride);
        pOutY = (ae_int16x8 *)(outY + i*stride);
        pOutU = (ae_int16x8 *)(outU + i*stride);
        pOutV = (ae_int16x8 *)(outV + i*stride);
        for (j = 0; j < (w >> 3); j++)
        {           
            AE_L16X4X2_IP(r0, r1, pInR, 8 * sizeof(int16_t));
            AE_L16X4X2_IP(g0, g1, pInG, 8 * sizeof(int16_t));
            AE_L16X4X2_IP(b0, b1, pInB, 8 * sizeof(int16_t));

            r0 = AE_MAX16(0, r0);
            g0 = AE_MAX16(0, g0);
            b0 = AE_MAX16(0, b0);

            AE_MULFD16X16X4WS(dy0, dy1, r0, g0, c0, c1);
            AE_MULAF16X4SS(dy0, dy1, b0, c2);
            y0 = AE_ROUND16X4F32SASYM(dy0, dy1);            

            AE_MULFD16X16X4WS(du0, du1, r0, g0, c3, c4);
            AE_MULAF16X4SS(du0, du1, b0, c5);
            u0 = AE_ROUND16X4F32SASYM(du0, du1);
            u0 = AE_ADD16S(u0, (ae_int16x4)16384);
            u0 = AE_MAX16(0, u0);

            AE_MULFD16X16X4WS(dv0, dv1, r0, g0, c6, c7);
            AE_MULAF16X4SS(dv0, dv1, b0, c8);
            v0 = AE_ROUND16X4F32SASYM(dv0, dv1);
            v0 = AE_ADD16S(v0, (ae_int16x4)16384);
            v0 = AE_MAX16(0, v0);

            r1 = AE_MAX16(0, r1);
            g1 = AE_MAX16(0, g1);
            b1 = AE_MAX16(0, b1);

            AE_MULFD16X16X4WS(dy0, dy1, r1, g1, c0, c1);
            AE_MULAF16X4SS(dy0, dy1, b1, c2);
            y1 = AE_ROUND16X4F32SASYM(dy0, dy1);            

            AE_MULFD16X16X4WS(du0, du1, r1, g1, c3, c4);
            AE_MULAF16X4SS(du0, du1, b1, c5);
            u1 = AE_ROUND16X4F32SASYM(du0, du1);
            u1 = AE_ADD16S(u1, (ae_int16x4)16384);
            u1 = AE_MAX16(0, u1);

            AE_MULFD16X16X4WS(dv0, dv1, r1, g1, c6, c7);
            AE_MULAF16X4SS(dv0, dv1, b1, c8);
            v1 = AE_ROUND16X4F32SASYM(dv0, dv1);
            v1 = AE_ADD16S(v1, (ae_int16x4)16384);
            v1 = AE_MAX16(0, v1);

            AE_S16X4X2_IP(y0, y1, pOutY, 2*sizeof(ae_int16x4));
            AE_S16X4X2_IP(u0, u1, pOutU, 2*sizeof(ae_int16x4));
            AE_S16X4X2_IP(v0, v1, pOutV, 2*sizeof(ae_int16x4));
        }
        if (w & 7)
        {    
            int off = (w & 7) << 1;   
            alr = AE_LA128_PP(pInR);
            alg = AE_LA128_PP(pInG);
            alb = AE_LA128_PP(pInB); 
            alv = AE_ZALIGN128();
            alu = AE_ZALIGN128();
            aly = AE_ZALIGN128();   
            AE_LAV16X4X2_XP(r0, r1, alr, pInR, off);
            AE_LAV16X4X2_XP(g0, g1, alg, pInG, off);
            AE_LAV16X4X2_XP(b0, b1, alb, pInB, off);

            r0 = AE_MAX16(0, r0);
            g0 = AE_MAX16(0, g0);
            b0 = AE_MAX16(0, b0);

            AE_MULFD16X16X4WS(dy0, dy1, r0, g0, c0, c1);
            AE_MULAF16X4SS(dy0, dy1, b0, c2);
            y0 = AE_ROUND16X4F32SASYM(dy0, dy1);            

            AE_MULFD16X16X4WS(du0, du1, r0, g0, c3, c4);
            AE_MULAF16X4SS(du0, du1, b0, c5);
            u0 = AE_ROUND16X4F32SASYM(du0, du1);
            u0 = AE_ADD16S(u0, (ae_int16x4)16384);
            u0 = AE_MAX16(0, u0);

            AE_MULFD16X16X4WS(dv0, dv1, r0, g0, c6, c7);
            AE_MULAF16X4SS(dv0, dv1, b0, c8);
            v0 = AE_ROUND16X4F32SASYM(dv0, dv1);
            v0 = AE_ADD16S(v0, (ae_int16x4)16384);
            v0 = AE_MAX16(0, v0);

            r1 = AE_MAX16(0, r1);
            g1 = AE_MAX16(0, g1);
            b1 = AE_MAX16(0, b1);

            AE_MULFD16X16X4WS(dy0, dy1, r1, g1, c0, c1);
            AE_MULAF16X4SS(dy0, dy1, b1, c2);
            y1 = AE_ROUND16X4F32SASYM(dy0, dy1);            

            AE_MULFD16X16X4WS(du0, du1, r1, g1, c3, c4);
            AE_MULAF16X4SS(du0, du1, b1, c5);
            u1 = AE_ROUND16X4F32SASYM(du0, du1);
            u1 = AE_ADD16S(u1, (ae_int16x4)16384);
            u1 = AE_MAX16(0, u1);

            AE_MULFD16X16X4WS(dv0, dv1, r1, g1, c6, c7);
            AE_MULAF16X4SS(dv0, dv1, b1, c8);
            v1 = AE_ROUND16X4F32SASYM(dv0, dv1);
            v1 = AE_ADD16S(v1, (ae_int16x4)16384);
            v1 = AE_MAX16(0, v1);

            AE_SAV16X4X2_XP(y0, y1, aly, pOutY, off);
            AE_SAV16X4X2_XP(u0, u1, alu, pOutU, off);
            AE_SAV16X4X2_XP(v0, v1, alv, pOutV, off);

            AE_SA128POS_FP(alv, pOutV);
            AE_SA128POS_FP(alu, pOutU);
            AE_SA128POS_FP(aly, pOutY);
        }      
    }
#endif
}
