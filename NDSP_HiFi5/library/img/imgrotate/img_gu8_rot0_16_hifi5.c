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
void img_gu8_rot0_16(void* pScr,void* outImg, imgsize_t * restrict outSz, const void* inImg, const imgsize_t * restrict inSz, int fillColor)
{
  int i, j, w, h, ostride, istride, h_out;
  const ae_int8x16 * restrict src = (const    ae_int8x16 *)inImg;
        ae_int16x8 * restrict dst0 = (      ae_int16x8 *)outImg;
        ae_int16x8 * restrict dst1 = (      ae_int16x8 *)outImg;

  ae_int8x8 x0, x1;
  ae_int16x4 y0, y1;
  ae_valignx2 aD, aS;
  (void)pScr;
  NASSERT_ALIGN(pScr, ALIGNMENT);
  NASSERT_ALIGN(outImg, ALIGNMENT);    /* inImg - not aligned, outImg - aligned */
  NASSERT_ALIGN(outImg, 1);

  imgsize_validate(inSz, 1, 0);
  imgsize_validate(outSz, 2, 1);
  ostride = outSz->stride;
  istride = inSz->stride;
  w = inSz->width;
  h = inSz->height;
  h_out = outSz->height;

  y0 = fillColor;
  /* filling row PADDING */
  dst1 = (ae_int16x8 *)((uintptr_t)outImg + ((h_out - (h + PADDING))*ostride)*sizeof(int16_t));
  for (j = 0; j<((PADDING*ostride) >> 3); j++)
  {
    AE_S16X4X2_IP(y0, y0, dst0, 2 * sizeof(ae_int16x4));
    AE_S16X4X2_IP(y0, y0, dst1, 2 * sizeof(ae_int16x4));
  }
  for (j = 0; j<((PADDING*ostride) & 7); j++)
  {
    AE_S16_0_IP(y0, castxcc(ae_int16, dst0), sizeof(int16_t));
    AE_S16_0_IP(y0, castxcc(ae_int16, dst1), sizeof(int16_t));
  }
  dst0 = (ae_int16x8 *)((uintptr_t)outImg + ((PADDING + h)*ostride)*sizeof(int16_t));
  for (j = 0; j<(int)(((h_out - 2 * PADDING - h)*ostride) >> 3); j++)
  {
    AE_S16X4X2_IP(y0, y0, dst0, 2 * sizeof(ae_int16x4));
  }
  for (j = 0; j<((h_out - 2 * PADDING - h) & 7); j++)
  {
    AE_S16_0_IP(y0, castxcc(ae_int16, dst0), sizeof(int16_t));
  }
  __Pragma("no_reorder");
  /* filling column PADDING */
  dst0 = (ae_int16x8 *)((uintptr_t)outImg + (PADDING * ostride)*sizeof(int16_t));
  dst1 = (ae_int16x8 *)((uintptr_t)outImg + (PADDING * ostride + ostride - ALIGNMENT / 2)*sizeof(int16_t));
  for (i = 0; i<h; i++)
  {
    AE_S16X4X2_XP(y0, y0, dst0, ostride * sizeof(int16_t));
    AE_S16X4X2_XP(y0, y0, dst1, ostride * sizeof(int16_t));
  }
  for (j = 0; j<((ostride - ALIGNMENT - w + 7) >> 3); j++)
  {
    dst0 = (ae_int16x8 *)((uintptr_t)outImg + (PADDING * ostride + ostride - ALIGNMENT / 2)*sizeof(int16_t)) - 1 - j;
    for (i = 0; i<h; i++)
    {
      AE_S16X4X2_XP(y0, y0, dst0, ostride * sizeof(ae_int16));
    }
  }

  __Pragma("no_reorder");

  // copy
  w=inSz->width;
  h=inSz->height;
  ostride = outSz->stride;
  istride = inSz->stride;

  aD = AE_ZALIGN128();

  for (i=0; i<h; i++)
  { 
    src = (const ae_int8x16 *)((int8_t*)inImg + i*istride);
    dst0 = (      ae_int16x8 *)((int16_t*)outImg +  (PADDING*ostride + (ALIGNMENT / 2)) + i*ostride);
    aS = AE_LA128_PP(src);
    for (j = 0; j<(w >> 4); j++)
    {
      AE_LA8X8X2_IP(x0, x1, aS, src);
      AE_CVTI16X4X2F8U(y0, y1, x0, 7);
      AE_S16X4X2_IP(y0, y1, dst0, 2 * sizeof(ae_int16x4));
      AE_CVTI16X4X2F8U(y0, y1, x1, 7);
      AE_S16X4X2_IP(y0, y1, dst0, 2 * sizeof(ae_int16x4));
    }
  }
  __Pragma("no_reorder");
  for (i = 0; i<(w & 15); i++)
  {
    src  = (const ae_int8x16 *)((int8_t*)inImg + w - 1 - i);
    dst0 = (      ae_int16x8 *)((int16_t*)outImg + (PADDING*ostride + (ALIGNMENT / 2)) + w - 1 - i);
    for (j = 0; j<h; j++)
    {
      ae_int8x8 p;
      AE_L8_XP(p, castxcc(ae_int8, src), sizeof(int8_t)*istride);
      AE_CVTI16X4X2F8U(y0, y0, p, 7);
      AE_S16_0_XP(y0, castxcc(ae_int16, dst0), sizeof(int16_t)*ostride);
    }
  }

  outSz->width=w+ALIGNMENT;
  outSz->height=h+PADDING*2;
} /* img_gu8_rot0_16() */