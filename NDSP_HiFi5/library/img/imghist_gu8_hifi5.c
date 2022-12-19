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

/*------------------------------------------------------------------------------
  Histogram of image 
  Functions compute the histogram of input image. First bin corresponds minimum
  intensity (0 for signed 16-bit and unsigned 8-bit data, -128 for signed 8-bit)
  and the last bin corresponds to the maximum intensity (32767 for signed 16-bit,
  255 for unsigned 8-bit data, 127  for signed 8-bit data)

  Image formats:
  gu8    8-bit unsigned grayscale data
  gs8    8-bit signed grayscale data
  gs16   16-bit grayscale data

  Input:
  inImg   input image
  sz      image size
  M       number of bins
  res->h  pointer to the histogram array
  Output:
  res[1]  histogram values, mean and variance

  Restrictions:
  see general restrictions applied for all images for fast/generic 
  functions
--------------------------------------------------------------------------------*/
void imghist_gu8  (imghist_t * restrict res, const void * restrict inImg, const imgsize_t* sz, int M)
{
    ae_int64 var64;
    ae_int16x4 x;
    ae_int32x2 z, t;
    ae_int32x2 d;
    ae_int32x2 d0, d1, d2, d3;
    ae_valign al;
    ae_valignx2 alw;
    uint8_t p;
    int32_t* h = res->h;
    int32_t mean, N, rN, idx;
    int64_t var;
    int32_t rQ31 = 8421504;//(2^8-1)^(-1) Q31
    const signed char * restrict pIn;
    const ae_int8x16 * restrict pIn_;
    const uint8_t * restrict in = (const uint8_t *)inImg;
    int i, j, k;
    int s, sN;
    int w = (int)sz->width;
    int istride = sz->stride;
    NASSERT(h);
    res->M = M;
    res->mean = 0;
    res->var = 0;
    for (k = 0; k<M; k++) h[k] = 0;
    if (M<1) return;
    imgsize_validate(sz, 1, 0);
    NASSERT(inImg);
    NASSERT_ALIGN(inImg, 1);
    mean = 0;
    N = (int)sz->height*(int)sz->width;

    ae_int64 mean64 = 0;
    ae_int64 mean64_1 = 0;
    ae_int64 mean64_2 = 0;
    ae_int64 mean64_3 = 0;
    for (i = 0; i < (int)sz->height; i++)
    {
        ae_int8x8 x0, x1;
        pIn_ = (const ae_int8x16 *)(in + i*istride);
        alw = AE_LA128_PP(pIn_);
        for (j = 0; j < (w >> 4); j++)
        {   /* 2 cycles unroll=1 */
            AE_LA8X8X2_IP(x0, x1, alw, pIn_);
            x0 = AE_MOVINT8X8_FROMINT16X4(AE_XOR16(AE_MOVINT16X4_FROMINT8X8(x0), 0x8080));
            x1 = AE_MOVINT8X8_FROMINT16X4(AE_XOR16(AE_MOVINT16X4_FROMINT8X8(x1), 0x8080));
            AE_MULAAAA2Q16X8(mean64, mean64_1, 1, 1, x0);
            AE_MULAAAA2Q16X8(mean64_2, mean64_3, 1, 1, x1);
        }
        if (w & 8)
        {
            al = AE_LA64_PP(pIn_);
            AE_LA8X8_IP(x0, al, castxcc(ae_int8x8, pIn_));
            x0 = AE_MOVINT8X8_FROMINT16X4(AE_XOR16(AE_MOVINT16X4_FROMINT8X8(x0), 0x8080));
            AE_MULAAAA2Q16X8(mean64, mean64_1, 1, 1, x0);
        }
        for (j = (w&(~7)); j<w; j++)
        {
            p = in[i*sz->stride + j];
            mean += p;
        }
    }
    mean64_2 += mean64_3;
    mean64 += mean64_1;
    mean64 += mean64_2;
    mean64 = mean64 + (ae_int64)mean + (ae_int64)(128*((int)sz->height)*(w&(~7)));

    // reciprocal of N
    {
        ae_int32x2 Z, E, Y;
        int sy;
        Y = AE_MOVDA32(N);
        sy = AE_NSAZ32_L(Y); Y = AE_SLAA32(Y, sy);
        /* first approximation */
        Z = AE_SUB32(AE_MOVDA32(0xBAEC0000), Y);
        /* 4 iterations to achieve 1 LSB for reciprocal */
        E = 0x40000000; AE_MULSFP32X2RAS(E, Y, Z); AE_MULAFP32X2RAS(Z, Z, AE_ADD32(E, E));
        E = 0x40000000; AE_MULSFP32X2RAS(E, Y, Z); AE_MULAFP32X2RAS(Z, Z, AE_ADD32(E, E));
        E = 0x40000000; AE_MULSFP32X2RAS(E, Y, Z); AE_MULAFP32X2RAS(Z, Z, AE_ADD32(E, E));
        E = 0x40000000; AE_MULSFP32X2RAS(E, Y, Z); AE_MULAFP32X2RAS(Z, Z, AE_ADD32(E, E));
        /* */
        rN = AE_MOVAD32_H(Z);
        sN = sy;
    }

    //Divide mean by N
    s = AE_NSA64(mean64);
    mean64 = AE_SLAA64(mean64, s);
    mean64 = AE_SRAA64(mean64, 32);
    mean = (int32_t)((int64_t)mean64);
    mean = AE_MOVAD32_H(AE_MULFP32X2RS(mean, rN));
    mean = AE_MOVAD32_H(AE_SLAA32(mean, sN - 14 - (s - 32)));
    res->mean = mean;

    var64 = 0;

    for (i = 0; i < (((int)sz->height)>>1); i++)
    {
        pIn = (const signed char *)(in + 2*i*istride);
        al = AE_LA64_PP(pIn);
        for (j = 0; j < (w >> 2); j++)
        {            
            AE_LA8X4U_IP(x , al, pIn);
            AE_MUL16X4(d0, d1, x, (M - 1));
            AE_MULF2P32X4RS(d0, d1, d0, d1, rQ31, rQ31);
            d0 = AE_SLAI32(d0, 2);
            d1 = AE_SLAI32(d1, 2);

            /*increment histogam's bins*/
            idx = AE_MOVAD32_H(d0);
            t = AE_L32_X((const ae_int32*)h, idx);
            t = AE_ADD32(t, 1);
            AE_S32_L_X(t, (ae_int32*)h, idx);

            idx = AE_MOVAD32_L(d0);
            t = AE_L32_X((const ae_int32*)h, idx);
            t = AE_ADD32(t, 1);
            AE_S32_L_X(t, (ae_int32*)h, idx);

            idx = AE_MOVAD32_H(d1);
            t = AE_L32_X((const ae_int32*)h, idx);
            t = AE_ADD32(t, 1);
            AE_S32_L_X(t, (ae_int32*)h, idx);

            idx = AE_MOVAD32_L(d1);
            t = AE_L32_X((const ae_int32*)h, idx);
            t = AE_ADD32(t, 1);
            AE_S32_L_X(t, (ae_int32*)h, idx);

            d2 = AE_AND32(AE_MOVINT32X2_FROMINT16X4(x), 0xffff0000);
            d3 = AE_SRLA32(AE_MOVINT32X2_FROMINT16X4(x), -16);

            d = AE_SUB32(mean, d2);    // Q16
            AE_MULAAD32_HH_LL(var64, d, d); // Q32 
            d = AE_SUB32(mean, d3);    // Q16
            AE_MULAAD32_HH_LL(var64, d, d); // Q32 
        }
        pIn = (const signed char *)(in + (2*i+1)*istride);
        al = AE_LA64_PP(pIn);
        for (j = 0; j < (w >> 2); j++)
        {            
            AE_LA8X4U_IP(x , al, pIn);
            AE_MUL16X4(d0, d1, x, (M - 1));
            AE_MULF2P32X4RS(d0, d1, d0, d1, rQ31, rQ31);
            d0 = AE_SLAI32(d0, 2);
            d1 = AE_SLAI32(d1, 2);

            /*increment histogam's bins*/
            idx = AE_MOVAD32_H(d0);
            t = AE_L32_X((const ae_int32*)h, idx);
            t = AE_ADD32(t, 1);
            AE_S32_L_X(t, (ae_int32*)h, idx);

            idx = AE_MOVAD32_L(d0);
            t = AE_L32_X((const ae_int32*)h, idx);
            t = AE_ADD32(t, 1);
            AE_S32_L_X(t, (ae_int32*)h, idx);

            idx = AE_MOVAD32_H(d1);
            t = AE_L32_X((const ae_int32*)h, idx);
            t = AE_ADD32(t, 1);
            AE_S32_L_X(t, (ae_int32*)h, idx);

            idx = AE_MOVAD32_L(d1);
            t = AE_L32_X((const ae_int32*)h, idx);
            t = AE_ADD32(t, 1);
            AE_S32_L_X(t, (ae_int32*)h, idx);

            d2 = AE_AND32(AE_MOVINT32X2_FROMINT16X4(x), 0xffff0000);
            d3 = AE_SRLA32(AE_MOVINT32X2_FROMINT16X4(x), -16);

            d = AE_SUB32(mean, d2);    // Q16
            AE_MULAAD32_HH_LL(var64, d, d); // Q32 
            d = AE_SUB32(mean, d3);    // Q16
            AE_MULAAD32_HH_LL(var64, d, d); // Q32 
        }
    }
    for (i = (((int)sz->height)&(~1)); i < (int)sz->height; i++)
    {
        pIn = (const signed char *)(in + i*istride);
        al = AE_LA64_PP(pIn);
        for (j = 0; j < (w >> 2); j++)
        {            
            AE_LA8X4U_IP(x , al, pIn);
            AE_MUL16X4(d0, d1, x, (M - 1));
            AE_MULF2P32X4RS(d0, d1, d0, d1, rQ31, rQ31);
            d0 = AE_SLAI32(d0, 2);
            d1 = AE_SLAI32(d1, 2);

            /*increment histogam's bins*/
            idx = AE_MOVAD32_H(d0);
            t = AE_L32_X((const ae_int32*)h, idx);
            t = AE_ADD32(t, 1);
            AE_S32_L_X(t, (ae_int32*)h, idx);

            idx = AE_MOVAD32_L(d0);
            t = AE_L32_X((const ae_int32*)h, idx);
            t = AE_ADD32(t, 1);
            AE_S32_L_X(t, (ae_int32*)h, idx);

            idx = AE_MOVAD32_H(d1);
            t = AE_L32_X((const ae_int32*)h, idx);
            t = AE_ADD32(t, 1);
            AE_S32_L_X(t, (ae_int32*)h, idx);

            idx = AE_MOVAD32_L(d1);
            t = AE_L32_X((const ae_int32*)h, idx);
            t = AE_ADD32(t, 1);
            AE_S32_L_X(t, (ae_int32*)h, idx);

            d2 = AE_AND32(AE_MOVINT32X2_FROMINT16X4(x), 0xffff0000);
            d3 = AE_SRLA32(AE_MOVINT32X2_FROMINT16X4(x), -16);

            d = AE_SUB32(mean, d2);    // Q16
            AE_MULAAD32_HH_LL(var64, d, d); // Q32 
            d = AE_SUB32(mean, d3);    // Q16
            AE_MULAAD32_HH_LL(var64, d, d); // Q32 
        }
    }
    if (w&3)
    {
        for (i = 0; i < (int)sz->height; i++)
        {

            for (j = (w&(~3)); j < w; j++)
            {
                p = in[i*istride + j];
                x = AE_MOVDA16(p);
                AE_MUL16X4(d0, d1, x, (M - 1));
                d0 = AE_MULFP32X2RS(d0, rQ31);
                idx = AE_MOVAD32_H(d0) << 2;
                t = AE_L32_X((const ae_int32*)h, idx);
                t = AE_ADD32(t, 1);
                AE_S32_L_X(t, (ae_int32*)h, idx);
                d1 = AE_AND32(AE_MOVINT32X2_FROMINT16X4(x), 0xffff0000);
                d = AE_SUB32(mean, d1);    // Q16
                AE_MULA32_HH(var64, d, d); // Q32 
            }
        }
    }

    s = AE_NSA64(var64);
    var64 = AE_SLAA64(var64, s);
    var64 = AE_SRLI64(var64, 32);
    var = (int64_t)var64;
    z = AE_MULFP32X2RS((int32_t)var, rN);
    s = s - sN + 30;
    z = AE_SRAA32RS(z, s);
    res->var = AE_MOVAD32_H(z);
}
