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
#include "imgrotate_common.h"

/*--------------------------------------------------------------------------
rotation by 0,90,180,270 conterclockwise with format conversion
the image is placed beginning from the row PADDING and from the column
ALIGNMENT/2. All paddings are written with black pixels
8-bit gs  -> signed 16-bit 1 pixel per 16-bit word
16-bit gs -> signed 16-bit 1 pixel per 16-bit word
RGB       -> signed 16-bit 1 pixel per 64-bit word
Input:
inImg, inSz  image and its size
outSz        output image size
fillColor    filling color
Output:
outImg       image
--------------------------------------------------------------------------*/
void img_gu8_rot90_16(void* pScr,void* outImg, imgsize_t * restrict outSz, const void* inImg, const imgsize_t * restrict inSz, int fillColor)
{
#if 1
  int i, j, w, h, ostride, istride, h_out;
  const ae_int8  * restrict src0 = (const ae_int8 *)inImg;
  const ae_int8  * restrict src1 = (const ae_int8 *)inImg;
  const ae_int8  * restrict src2 = (const ae_int8 *)inImg;
  const ae_int8  * restrict src3 = (const ae_int8 *)inImg;
  const ae_int8  * restrict src4 = (const ae_int8 *)inImg;
  const ae_int8  * restrict src5 = (const ae_int8 *)inImg;
  const ae_int8  * restrict src6 = (const ae_int8 *)inImg;
  const ae_int8  * restrict src7 = (const ae_int8 *)inImg;
        ae_int16 * restrict dst0 = (      ae_int16 *)outImg;
        ae_int16 * restrict dst1 = (      ae_int16 *)outImg;
        ae_int16 * restrict dst2;
        ae_int16 * restrict dst3;
        ae_int16 * restrict dst4;
        ae_int16 * restrict dst5;
        ae_int16 * restrict dst6;
        ae_int16 * restrict dst7;

  ae_int16x4 y0;
  ae_valignx2 aX, aY;
  (void)pScr;
  NASSERT_ALIGN(pScr, ALIGNMENT);
  NASSERT_ALIGN(outImg, ALIGNMENT);    /* inImg - not aligned, outImg - aligned */
  NASSERT_ALIGN(outImg, 1);

  imgsize_validate(inSz, 1, 0);
  imgsize_validate(outSz, 2, 1);

  static const uint8_t ALIGN(ALIGNMENT) dsel_tbl[] = { 11 | (15 << 4), 3 | (7 << 4), 10 | (14 << 4), 2 | (6 << 4), 9 | (13 << 4), 1 | (5 << 4), 8 | (12 << 4), 0 | (4 << 4),
                                                       11 | (15 << 4), 10 | (14 << 4), 3 | (7 << 4), 2 | (6 << 4), 9 | (13 << 4), 8 | (12 << 4), 1 | (5 << 4), 0 | (4 << 4),
                                                       15 | (11 << 4), 14 | (10 << 4), 13 | (9 << 4), 12 | (8 << 4), 7 | (3 << 4), 6 | (2 << 4), 5 | (1 << 4), 4 | (0 << 4) };
  
  ae_int8x8 dsel0 = AE_L8X8_I((const ae_int8x8*)dsel_tbl, 0);
  ae_int8x8 dsel1 = AE_L8X8_I((const ae_int8x8*)dsel_tbl, sizeof(ae_int8x8));
  ae_int8x8 dsel2 = AE_L8X8_I((const ae_int8x8*)dsel_tbl, 2 * sizeof(ae_int8x8));

  w = inSz->height;
  h = inSz->width;
  ostride = outSz->stride;
  istride = inSz->stride;
  h_out = outSz->height;


  y0 = fillColor;
  /* filling row PADDING */
  dst1 = (ae_int16 *)((uintptr_t)outImg + ((h_out - (h + PADDING))*ostride)*sizeof(int16_t));
  for (j = 0; j<((PADDING*ostride) >> 3); j++)
  {
    AE_S16X4X2_IP(y0, y0, castxcc(ae_int16x8, dst0), 2 * sizeof(ae_int16x4));
    AE_S16X4X2_IP(y0, y0, castxcc(ae_int16x8, dst1), 2 * sizeof(ae_int16x4));
  }
  for (j = 0; j<((PADDING*ostride) & 7); j++)
  {
    AE_S16_0_IP(y0, castxcc(ae_int16, dst0), sizeof(int16_t));
    AE_S16_0_IP(y0, castxcc(ae_int16, dst1), sizeof(int16_t));
  }
  dst0 = (ae_int16 *)((uintptr_t)outImg + ((PADDING + h)*ostride)*sizeof(int16_t));
  for (j = 0; j<(int)(((h_out - 2 * PADDING - h)*ostride) >> 3); j++)
  {
    AE_S16X4X2_IP(y0, y0, castxcc(ae_int16x8, dst0), 2 * sizeof(ae_int16x4));
  }
  for (j = 0; j<((h_out - 2 * PADDING - h) & 7); j++)
  {
    AE_S16_0_IP(y0, castxcc(ae_int16, dst0), sizeof(int16_t));
  }
  __Pragma("no_reorder");
  /* filling column PADDING */
  dst0 = (ae_int16 *)((uintptr_t)outImg + (PADDING * ostride)*sizeof(int16_t));
  dst1 = (ae_int16 *)((uintptr_t)outImg + (PADDING * ostride + ostride - ALIGNMENT / 2)*sizeof(int16_t));
  for (i = 0; i<h; i++)
  {
    AE_S16X4X2_XP(y0, y0, castxcc(ae_int16x8, dst0), ostride * sizeof(int16_t));
    AE_S16X4X2_XP(y0, y0, castxcc(ae_int16x8, dst1), ostride * sizeof(int16_t));
  }
  for (j = 0; j<((ostride - ALIGNMENT - w + 7) >> 3); j++)
  {
    dst0 = (ae_int16 *)((uintptr_t)outImg + (PADDING * ostride + ostride - ALIGNMENT / 2)*sizeof(int16_t)) - 8 - 8*j;
    for (i = 0; i<h; i++)
    {
      AE_S16X4X2_XP(y0, y0, castxcc(ae_int16x8, dst0), ostride * sizeof(ae_int16));
    }
  }

  __Pragma("no_reorder");

  /* rotate by 90 degrees counterclockwise */
  for (i=0; i<(w>>3); i++)
  {
    ae_int8x8 x0, x1, x2, x3;
    ae_int8x8 x4, x5, x6, x7;
    ae_int8x8 y0, y1, y2, y3;
    ae_int8x8 y4, y5, y6, y7;
    ae_int8x8 z0, z1, z2, z3;
    ae_int8x8 z4, z5, z6, z7;
    ae_int16x4 w0, w1, w2, w3;

    src0 = (const ae_int8 *)((uintptr_t)inImg + ((0 + 8 * i) * istride + (h)-(h & 15)) * sizeof(int8_t)) - 16;
    src1 = (const ae_int8 *)((uintptr_t)inImg + ((1 + 8 * i) * istride + (h)-(h & 15)) * sizeof(int8_t)) - 16;

    dst0 = (ae_int16 *)((uintptr_t)outImg + ((PADDING*ostride + (ALIGNMENT / 2)) + ((h & 15) + 0)* ostride)*sizeof(int16_t)) + 8 * i;
    dst1 = (ae_int16 *)((uintptr_t)outImg + ((PADDING*ostride + (ALIGNMENT / 2)) + ((h & 15) + 1)* ostride)*sizeof(int16_t)) + 8 * i;
    dst2 = (ae_int16 *)((uintptr_t)outImg + ((PADDING*ostride + (ALIGNMENT / 2)) + ((h & 15) + 2)* ostride)*sizeof(int16_t)) + 8 * i;
    dst3 = (ae_int16 *)((uintptr_t)outImg + ((PADDING*ostride + (ALIGNMENT / 2)) + ((h & 15) + 3)* ostride)*sizeof(int16_t)) + 8 * i;
    dst4 = (ae_int16 *)((uintptr_t)outImg + ((PADDING*ostride + (ALIGNMENT / 2)) + ((h & 15) + 4)* ostride)*sizeof(int16_t)) + 8 * i;
    dst5 = (ae_int16 *)((uintptr_t)outImg + ((PADDING*ostride + (ALIGNMENT / 2)) + ((h & 15) + 5)* ostride)*sizeof(int16_t)) + 8 * i;
    dst6 = (ae_int16 *)((uintptr_t)outImg + ((PADDING*ostride + (ALIGNMENT / 2)) + ((h & 15) + 6)* ostride)*sizeof(int16_t)) + 8 * i;
    dst7 = (ae_int16 *)((uintptr_t)outImg + ((PADDING*ostride + (ALIGNMENT / 2)) + ((h & 15) + 7)* ostride)*sizeof(int16_t)) + 8 * i;


    for (j=0; j<(h>>4); j++) //48
    {
      aX = AE_LA128_PP(src0); AE_LA8X8X2_IP(z0, x0, aX, castxcc(ae_int8x16, src0)); src0 += 2 * istride - 16;
      aX = AE_LA128_PP(src1); AE_LA8X8X2_IP(z1, x1, aX, castxcc(ae_int8x16, src1)); src1 += 2 * istride - 16;
      aX = AE_LA128_PP(src0); AE_LA8X8X2_IP(z2, x2, aX, castxcc(ae_int8x16, src0)); src0 += 2 * istride - 16;
      aX = AE_LA128_PP(src1); AE_LA8X8X2_IP(z3, x3, aX, castxcc(ae_int8x16, src1)); src1 += 2 * istride - 16;
      aX = AE_LA128_PP(src0); AE_LA8X8X2_IP(z4, x4, aX, castxcc(ae_int8x16, src0)); src0 += 2 * istride - 16;
      aX = AE_LA128_PP(src1); AE_LA8X8X2_IP(z5, x5, aX, castxcc(ae_int8x16, src1)); src1 += 2 * istride - 16;
      aX = AE_LA128_PP(src0); AE_LA8X8X2_IP(z6, x6, aX, castxcc(ae_int8x16, src0)); src0 += -3 * 2 * istride - 32;
      aX = AE_LA128_PP(src1); AE_LA8X8X2_IP(z7, x7, aX, castxcc(ae_int8x16, src1)); src1 += -3 * 2 * istride - 32;

  
      AE_DSEL8X8(y0, y1, x0, x1, dsel0);
      AE_DSEL8X8(y2, y3, x2, x3, dsel0);
      AE_DSEL8X8(y4, y5, x4, x5, dsel0);
      AE_DSEL8X8(y6, y7, x6, x7, dsel0);

      AE_DSEL8X8(x0, x1, y0, y2, dsel1);
      AE_DSEL8X8(x2, x3, y1, y3, dsel1);
      AE_DSEL8X8(x4, x5, y4, y6, dsel1);
      AE_DSEL8X8(x6, x7, y5, y7, dsel1);

      AE_DSEL8X8(y0, y1, x0, x4, dsel2);
      AE_DSEL8X8(y2, y3, x1, x5, dsel2);
      AE_DSEL8X8(y4, y5, x2, x6, dsel2);
      AE_DSEL8X8(y6, y7, x3, x7, dsel2); 
       
      AE_CVTI16X4X2F8U(w0, w1, y6, 7);
      AE_CVTI16X4X2F8U(w2, w3, y7, 7);
      AE_SA16X4X2_IP(w0, w1, aY, castxcc(ae_int16x8, dst0));  AE_SA128POS_FP(aY, dst0); dst0 += 8 * ostride - 8;
      AE_SA16X4X2_IP(w2, w3, aY, castxcc(ae_int16x8, dst1));  AE_SA128POS_FP(aY, dst1); dst1 += 8 * ostride - 8;
      AE_CVTI16X4X2F8U(w0, w1, y4, 7);
      AE_CVTI16X4X2F8U(w2, w3, y5, 7);
      AE_SA16X4X2_IP(w0, w1, aY, castxcc(ae_int16x8, dst2));  AE_SA128POS_FP(aY, dst2); dst2 += 8 * ostride - 8;
      AE_SA16X4X2_IP(w2, w3, aY, castxcc(ae_int16x8, dst3));  AE_SA128POS_FP(aY, dst3); dst3 += 8 * ostride - 8;
      AE_CVTI16X4X2F8U(w0, w1, y2, 7);
      AE_CVTI16X4X2F8U(w2, w3, y3, 7);
      AE_SA16X4X2_IP(w0, w1, aY, castxcc(ae_int16x8, dst4));  AE_SA128POS_FP(aY, dst4); dst4 += 8 * ostride - 8;
      AE_SA16X4X2_IP(w2, w3, aY, castxcc(ae_int16x8, dst5));  AE_SA128POS_FP(aY, dst5); dst5 += 8 * ostride - 8;
      AE_CVTI16X4X2F8U(w0, w1, y0, 7);
      AE_CVTI16X4X2F8U(w2, w3, y1, 7);
      AE_SA16X4X2_IP(w0, w1, aY, castxcc(ae_int16x8, dst6));  AE_SA128POS_FP(aY, dst6); dst6 += 8 * ostride - 8;
      AE_SA16X4X2_IP(w2, w3, aY, castxcc(ae_int16x8, dst7));  AE_SA128POS_FP(aY, dst7); dst7 += 8 * ostride - 8;

      AE_DSEL8X8(y0, y1, z0, z1, dsel0);
      AE_DSEL8X8(y2, y3, z2, z3, dsel0);
      AE_DSEL8X8(y4, y5, z4, z5, dsel0);
      AE_DSEL8X8(y6, y7, z6, z7, dsel0);

      AE_DSEL8X8(x0, x1, y0, y2, dsel1);
      AE_DSEL8X8(x2, x3, y1, y3, dsel1);
      AE_DSEL8X8(x4, x5, y4, y6, dsel1);
      AE_DSEL8X8(x6, x7, y5, y7, dsel1);

      AE_DSEL8X8(y0, y1, x0, x4, dsel2);
      AE_DSEL8X8(y2, y3, x1, x5, dsel2);
      AE_DSEL8X8(y4, y5, x2, x6, dsel2);
      AE_DSEL8X8(y6, y7, x3, x7, dsel2);

      AE_CVTI16X4X2F8U(w0, w1, y6, 7);
      AE_CVTI16X4X2F8U(w2, w3, y7, 7);
      AE_SA16X4X2_IP(w0, w1, aY, castxcc(ae_int16x8, dst0));  AE_SA128POS_FP(aY, dst0); dst0 += 8 * ostride - 8;
      AE_SA16X4X2_IP(w2, w3, aY, castxcc(ae_int16x8, dst1));  AE_SA128POS_FP(aY, dst1); dst1 += 8 * ostride - 8;
      AE_CVTI16X4X2F8U(w0, w1, y4, 7);
      AE_CVTI16X4X2F8U(w2, w3, y5, 7);
      AE_SA16X4X2_IP(w0, w1, aY, castxcc(ae_int16x8, dst2));  AE_SA128POS_FP(aY, dst2); dst2 += 8 * ostride - 8;
      AE_SA16X4X2_IP(w2, w3, aY, castxcc(ae_int16x8, dst3));  AE_SA128POS_FP(aY, dst3); dst3 += 8 * ostride - 8;
      AE_CVTI16X4X2F8U(w0, w1, y2, 7);
      AE_CVTI16X4X2F8U(w2, w3, y3, 7);
      AE_SA16X4X2_IP(w0, w1, aY, castxcc(ae_int16x8, dst4));  AE_SA128POS_FP(aY, dst4); dst4 += 8 * ostride - 8;
      AE_SA16X4X2_IP(w2, w3, aY, castxcc(ae_int16x8, dst5));  AE_SA128POS_FP(aY, dst5); dst5 += 8 * ostride - 8;
      AE_CVTI16X4X2F8U(w0, w1, y0, 7);
      AE_CVTI16X4X2F8U(w2, w3, y1, 7);
      AE_SA16X4X2_IP(w0, w1, aY, castxcc(ae_int16x8, dst6));  AE_SA128POS_FP(aY, dst6); dst6 += 8 * ostride - 8;
      AE_SA16X4X2_IP(w2, w3, aY, castxcc(ae_int16x8, dst7));  AE_SA128POS_FP(aY, dst7); dst7 += 8 * ostride - 8;

    }
  }
  __Pragma("no_reorder");
  for (j = 0; j<(h & 15); j++)
  {
     src0 = (const ae_int8 *)((uintptr_t)inImg + (0 * istride + (h - 1 - j)) * sizeof(int8_t));
     src1 = (const ae_int8 *)((uintptr_t)inImg + (1 * istride + (h - 1 - j)) * sizeof(int8_t));
     src2 = (const ae_int8 *)((uintptr_t)inImg + (2 * istride + (h - 1 - j)) * sizeof(int8_t));
     src3 = (const ae_int8 *)((uintptr_t)inImg + (3 * istride + (h - 1 - j)) * sizeof(int8_t));
     src4 = (const ae_int8 *)((uintptr_t)inImg + (4 * istride + (h - 1 - j)) * sizeof(int8_t));
     src5 = (const ae_int8 *)((uintptr_t)inImg + (5 * istride + (h - 1 - j)) * sizeof(int8_t));
     src6 = (const ae_int8 *)((uintptr_t)inImg + (6 * istride + (h - 1 - j)) * sizeof(int8_t));
     src7 = (const ae_int8 *)((uintptr_t)inImg + (7 * istride + (h - 1 - j)) * sizeof(int8_t));
    
     dst0 = (ae_int16 *)((uintptr_t)outImg + ((PADDING*ostride + (ALIGNMENT / 2)) + j * ostride)*sizeof(int16_t));
     for (i = 0; i<(w >> 3); i++) 
     {
        ae_int8x8 p0, p1, p2, p3;
        ae_int8x8 p4, p5, p6, p7;
        ae_int16x4 w0, w1, w2, w3;
        ae_int16x4 w4, w5, w6, w7;
        AE_L8_XP(p0, castxcc(ae_int8, src0), (int)(sizeof(int8_t) * 8 * istride));
        AE_L8_XP(p1, castxcc(ae_int8, src1), (int)(sizeof(int8_t) * 8 * istride));
        AE_L8_XP(p2, castxcc(ae_int8, src2), (int)(sizeof(int8_t) * 8 * istride));
        AE_L8_XP(p3, castxcc(ae_int8, src3), (int)(sizeof(int8_t) * 8 * istride));
        AE_L8_XP(p4, castxcc(ae_int8, src4), (int)(sizeof(int8_t) * 8 * istride));
        AE_L8_XP(p5, castxcc(ae_int8, src5), (int)(sizeof(int8_t) * 8 * istride));
        AE_L8_XP(p6, castxcc(ae_int8, src6), (int)(sizeof(int8_t) * 8 * istride));
        AE_L8_XP(p7, castxcc(ae_int8, src7), (int)(sizeof(int8_t) * 8 * istride));
        
        AE_CVTI16X4X2F8U(w0, w0, p0, 7);
        AE_CVTI16X4X2F8U(w1, w1, p1, 7);
        AE_CVTI16X4X2F8U(w2, w2, p2, 7);
        AE_CVTI16X4X2F8U(w3, w3, p3, 7);
        AE_CVTI16X4X2F8U(w4, w4, p4, 7);
        AE_CVTI16X4X2F8U(w5, w5, p5, 7);
        AE_CVTI16X4X2F8U(w6, w6, p6, 7);
        AE_CVTI16X4X2F8U(w7, w7, p7, 7);
        AE_S16_0_IP(w0, castxcc(ae_int16, dst0), sizeof(int16_t));
        AE_S16_0_IP(w1, castxcc(ae_int16, dst0), sizeof(int16_t));
        AE_S16_0_IP(w2, castxcc(ae_int16, dst0), sizeof(int16_t));
        AE_S16_0_IP(w3, castxcc(ae_int16, dst0), sizeof(int16_t));
        AE_S16_0_IP(w4, castxcc(ae_int16, dst0), sizeof(int16_t));
        AE_S16_0_IP(w5, castxcc(ae_int16, dst0), sizeof(int16_t));
        AE_S16_0_IP(w6, castxcc(ae_int16, dst0), sizeof(int16_t));
        AE_S16_0_IP(w7, castxcc(ae_int16, dst0), sizeof(int16_t));    
     }
  }
  __Pragma("no_reorder");
  for (i = 0; i<(w & 7); i++)
  {
    src0 = (const ae_int8 *)((uintptr_t)inImg + ((w - 1 - i)* istride + (h)-(h & 7)) * sizeof(int8_t) - 8 * sizeof(int8_t));
    dst0 = (     ae_int16 *)((uintptr_t)outImg + ((PADDING*ostride + (ALIGNMENT / 2)) + ((h & 7) + 0)* ostride + (w - 1 - i))*sizeof(int16_t));
  
    for (j = 0; j<(h >> 3); j++)
    {
      ae_int8x8 x0;
      ae_int16x4 w0, w1;
      ae_valign aX;
      aX = AE_LA64_PP(src0);
      AE_LA8X8_IP(x0, aX, castxcc(ae_int8x8, src0)); src0 -= 16;
      AE_CVTI16X4X2F8U(w1, w0, x0, 7); 
  
      AE_S16_0_XP(w0, castxcc(ae_int16, dst0), ostride*sizeof(int16_t));
      w0 = AE_SEL16_4321(w0, w0);
      AE_S16_0_XP(w0, castxcc(ae_int16, dst0), ostride*sizeof(int16_t));
      w0 = AE_SEL16_4321(w0, w0);
      AE_S16_0_XP(w0, castxcc(ae_int16, dst0), ostride*sizeof(int16_t));
      w0 = AE_SEL16_4321(w0, w0);
      AE_S16_0_XP(w0, castxcc(ae_int16, dst0), ostride*sizeof(int16_t));
      w0 = AE_SEL16_4321(w0, w0);
      AE_S16_0_XP(w1, castxcc(ae_int16, dst0), ostride*sizeof(int16_t));
      w1 = AE_SEL16_4321(w1, w1);
      AE_S16_0_XP(w1, castxcc(ae_int16, dst0), ostride*sizeof(int16_t));
      w1 = AE_SEL16_4321(w1, w1);
      AE_S16_0_XP(w1, castxcc(ae_int16, dst0), ostride*sizeof(int16_t));
      w1 = AE_SEL16_4321(w1, w1);
      AE_S16_0_XP(w1, castxcc(ae_int16, dst0), ostride*sizeof(int16_t));
   
    }
    src0 = (const ae_int8 *)((uintptr_t)inImg + ((w - 1 - i) * istride + (h - 1)) * sizeof(int8_t));
    dst0 = (ae_int16 *)((uintptr_t)outImg + ((PADDING*ostride + (ALIGNMENT / 2)) + w - 1 - i) *sizeof(int16_t));
    for (j = 0; j<(h & 7); j++)
    {
      ae_int8x8 x0;
      ae_int16x4 y0;
      AE_L8_IP(x0, castxcc(ae_int8, src0), -(int)sizeof(int8_t));
      AE_CVTI16X4X2F8U(y0, y0, x0, 7);
      AE_S16_0_XP(y0, castxcc(ae_int16, dst0), ostride*sizeof(int16_t));
    }
  }
  outSz->width=w+ALIGNMENT;
  outSz->height=h+PADDING*2;
#else
  int i,j,w,h,ostride,istride;
  const uint8_t * restrict in =(const uint8_t *)inImg;
        int16_t * restrict out=(      int16_t *)outImg;
  (void)pScr;
  NASSERT_ALIGN(pScr,ALIGNMENT);
  NASSERT_ALIGN(outImg,ALIGNMENT);    /* inImg - not aligned, outImg - aligned */    
  NASSERT_ALIGN(outImg,1);
  imgsize_validate(inSz,1,0);
  imgsize_validate(outSz,2,1);
  ostride = outSz->stride;
  istride = inSz->stride;
  for (i=0; i<PADDING*ostride; i++) *out++=fillColor;
  // rotate by 90 degrees counterclockwise
  w = inSz->height;
  h = inSz->width;
  for (i = 0; i<h; i++)
  {
    for (j = 0; j<ALIGNMENT / 2; j++) *out++ = fillColor;
    for (j = 0; j<w; j++) *out++ = ((int16_t)in[j*istride + (h - 1 - i)]) << 7;
    for (; j<(ostride - ALIGNMENT / 2); j++) *out++ = fillColor;
  }
  for (i=0; i<(int)(outSz->height-(h+PADDING))*ostride; i++) *out++=fillColor;
  outSz->width=w+ALIGNMENT;
  outSz->height=h+PADDING*2;
#endif
} /* img_gu8_rot90_16() */
