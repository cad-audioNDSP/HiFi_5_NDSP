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
void imgfastnorm_gu8  ( void * restrict outImg, const void * restrict inImg, const imgsize_t* sz, int minInt, int maxInt)
{
    ae_int32x2 dQ31,cQ23;
    ae_int16x4 cQ7, cQ8;
    const uint8_t * restrict in =(const uint8_t *)inImg;
          uint8_t * restrict out=(      uint8_t *)outImg;
    const signed char * restrict pIn;
    const ae_int8x16 * restrict pIn_;
          ae_int8x8 * restrict pOut;
    int h=(int)sz->height;
    int w=(int)sz->width;
    int istride=sz->stride;
    int i,j;
    int16_t pmax0,pmin0, tmp;
    int16_t pmax, pmin;
    ae_int8x8 pmaxv, pminv;
    NASSERT(inImg!=NULL);
    NASSERT(outImg!=NULL);
    NASSERT_ALIGN(inImg,ALIGNMENT);
    NASSERT_ALIGN(outImg,ALIGNMENT);
    imgsize_validate(sz,1,1);
    minInt=XT_MAX(minInt,0);
    maxInt=XT_MIN(maxInt,255);
    pmax0=0; pmin0=255;
    // find min/max values
    {
        ae_int8x8 p, p1;
        pmaxv = AE_MOVDA8((int8_t)(-128)); 
        pminv = AE_MOVDA8((int8_t)(127));
        for (i = 0; i<h; i++)
        {
            pIn_ = (const ae_int8x16*)(in + i*istride);
            for (j = 0; j<(w>>4); j++)
            {   /*2 cycles unroll=1*/
                AE_L8X8X2_IP(p, p1, pIn_, sizeof(ae_int8x16));
                p = AE_MOVINT8X8_FROMINT16X4(AE_XOR16(AE_MOVINT16X4_FROMINT8X8(p), 0x8080));
                p1 = AE_MOVINT8X8_FROMINT16X4(AE_XOR16(AE_MOVINT16X4_FROMINT8X8(p1), 0x8080));
                pmaxv = AE_MAX8(p, pmaxv);
                pmaxv = AE_MAX8(p1, pmaxv);
                pminv = AE_MIN8(p, pminv);               
                pminv = AE_MIN8(p1, pminv);
            }
            if (w & 8)
            {   
                AE_L8X8_IP(p, castxcc(ae_int8x8, pIn_), sizeof(ae_int8x8));
                p = AE_MOVINT8X8_FROMINT16X4(AE_XOR16(AE_MOVINT16X4_FROMINT8X8(p), 0x8080));
                pmaxv = AE_MAX8(p, pmaxv);
                pminv = AE_MIN8(p, pminv);
            }
            for (j = (w&(~7)); j < w; j++)
            {
                tmp = in[i*istride + j];
                pmin0 = XT_MIN(pmin0, tmp);
                pmax0 = XT_MAX(pmax0, tmp);
            }
        }
    }  
    pmax = (int8_t)(int16_t)AE_MOVINT16_FROMINT8(AE_RMAX8X8(pmaxv));
    pmin = (int8_t)(int16_t)AE_MOVINT16_FROMINT8(AE_RMIN8X8(pminv));
    pmin = XT_MIN(pmin0, pmin+128);
    pmax = XT_MAX(pmax0, pmax+128);

    if (pmax==pmin) 
    {   // force output to minInt
        maxInt=minInt;
        pmax=1; pmin=0;
    }
    
    {
        ae_int64 r=0x80000000ULL;
        ae_int32x2 t=AE_MOVDA32(pmax-pmin);
        xtbool2 ble1=AE_LE32(t,1);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);

        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
        AE_DIV64D32_H(r,t);    AE_DIV64D32_H(r,t);
    
        dQ31=AE_TRUNCA32F64S(AE_AND64(r,0x7fffffff),32);
        AE_MOVT32X2(dQ31,0x7fffffff,ble1);
    }

    cQ23=AE_MULFP32X2RAS((maxInt-minInt)<<23,dQ31);
    cQ7 = AE_ROUND16X4F32SASYM(cQ23, cQ23);
    cQ8 = AE_SLAI16(cQ7, 1);

    ae_int16x4 p0, p1, y0, y1;
    ae_int8x8 t0;
    for (i = 0; i<h; i++)
    {
        pIn = (const signed char*)(in + i*istride);
        pOut = (ae_int8x8*)(out + i*istride);
        for (j = 0; j<(w>>3); j++)
        {
            AE_L8X4F_IP(p0, pIn, 4*sizeof(int8_t));
            p0 = AE_SRLI16(p0, 1);    
            p0 = AE_SUB16(p0, pmin << 7);
            y0 = AE_MULFP16X4RAS(p0, cQ8);
            y0 = AE_ADD16S(y0, (ae_int16x4)minInt);
            AE_MINMAX16(y0, 0, 255);

            AE_L8X4F_IP(p1, pIn, 4 * sizeof(int8_t));
            p1 = AE_SRLI16(p1, 1);
            p1 = AE_SUB16(p1, pmin << 7);
            y1 = AE_MULFP16X4RAS(p1, cQ8);
            y1 = AE_ADD16S(y1, (ae_int16x4)minInt);
            AE_MINMAX16(y1, 0, 255);

            //save result 
            t0 = AE_SEL8X8I(AE_MOVINT8X8_FROMINT16X4(y0), AE_MOVINT8X8_FROMINT16X4(y1), 25); //Extract 8 lsb from each 16 bit value
            AE_S8X8_IP(t0, pOut, sizeof(ae_int16x4));
        }
        for (j = (w&~7); j < w; j++)
        {
            p0 = AE_MOVDA16(in[i*istride + j]); p0 = AE_SLAI16S(p0, 7);
            p0 = AE_SUB16(p0, pmin << 7);
            y0 = AE_MULFP16X4RAS(p0, cQ8);
            y0 = AE_ADD16S(y0, (ae_int16x4)minInt);
            AE_MINMAX16(y0, 0, 255);
            out[i*istride + j] = (uint8_t)AE_MOVAD16_0(y0);
        }
    }
}/*imgfastnorm_gs()*/
