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
void imgfast_gs16_rot180_16(void* pScr, void* outImg, imgsize_t * restrict outSz, const void* inImg, const imgsize_t * restrict inSz, int fillColor)
{
#if 1
  static const int16_t ALIGN(ALIGNMENT) dsel_tbl[] = { 0 | (4 << 8), 1 | (5 << 8), 2 | (6 << 8), 3 | (7 << 8) };

  int i, j, w, h, ostride, istride, h_out;
  const ae_int16x8 * restrict src = (const ae_int16x8 *)inImg;
        ae_int16x8 * restrict dst0 = (      ae_int16x8 *)outImg;
        ae_int16x8 * restrict dst1 = (      ae_int16x8 *)outImg;
  (void)pScr;
  ae_int16x4 x0, x1;
  ae_valignx2 aD;
  NASSERT_ALIGN(pScr,ALIGNMENT);
  NASSERT_ALIGN(outImg,ALIGNMENT);
  NASSERT_ALIGN(inImg ,ALIGNMENT);
  ae_int16x4 dsel0 = AE_L16X4_I((const ae_int16x4*)dsel_tbl, 0);

  imgsize_validate(inSz,2,1);
  imgsize_validate(outSz,2,1);
  ostride = outSz->stride;
  istride = inSz->stride;
  w = inSz->width;
  h = inSz->height;
  h_out = outSz->height;
  x0 = fillColor;
  aD = AE_ZALIGN128();
  /* filling row PADDING */
  dst1 = (ae_int16x8 *)((uintptr_t)outImg + ((h_out - (h + PADDING))*ostride)*sizeof(int16_t));
  for (j = 0; j<((PADDING*ostride) >> 3); j++)
  {
    AE_S16X4X2_IP(x0, x0, dst0, 2 * sizeof(ae_int16x4));
    AE_S16X4X2_IP(x0, x0, dst1, 2 * sizeof(ae_int16x4));
  }
  for (j = 0; j<((PADDING*ostride) & 7); j++)
  {
    AE_S16_0_IP(x0, castxcc(ae_int16, dst0), sizeof(int16_t));
    AE_S16_0_IP(x0, castxcc(ae_int16, dst1), sizeof(int16_t));
  }
  dst0 = (ae_int16x8 *)((uintptr_t)outImg + ((PADDING + h)*ostride)*sizeof(int16_t));
  for (j = 0; j<(int)(((h_out - 2 * PADDING - h)*ostride) >> 3); j++)
  {
    AE_S16X4X2_IP(x0, x0, dst0, 2 * sizeof(ae_int16x4));
  }
  for (j = 0; j<((h_out - 2 * PADDING - h) & 7); j++)
  {
    AE_S16_0_IP(x0, castxcc(ae_int16, dst0), sizeof(int16_t));
  }
  __Pragma("no_reorder");

  /* filling column PADDING */
  dst0 = (ae_int16x8 *)((uintptr_t)outImg + (PADDING * ostride)*sizeof(int16_t));
  dst1 = (ae_int16x8 *)((uintptr_t)outImg + (PADDING * ostride + ostride - ALIGNMENT / 2)*sizeof(int16_t));
  for (i = 0; i<h; i++)
  {
    AE_S16X4X2_XP(x0, x0, dst0, ostride * sizeof(int16_t));
    AE_S16X4X2_XP(x0, x0, dst1, ostride * sizeof(int16_t));
  }
  for (j = 0; j<((ostride - ALIGNMENT - w + 7) >> 3); j++)
  {
    dst0 = (ae_int16x8 *)((uintptr_t)outImg + (PADDING * ostride + ostride - ALIGNMENT / 2)*sizeof(int16_t)) - 1 - j;
    for (i = 0; i<h; i++)
    {
      AE_S16X4X2_XP(x0, x0, dst0, ostride * sizeof(ae_int16));
    }
  }

  __Pragma("no_reorder");
  // flip by 180 degrees
  src = (const ae_int16x8 *)((uintptr_t)inImg + ((h - 1)*istride + (w>>3)*8)*sizeof(int16_t)) - 1;
  aD = AE_ZALIGN128();
  dst0 = (ae_int16x8 *)((uintptr_t)outImg + (PADDING*ostride + (ALIGNMENT / 2) + (w & 7))*sizeof(int16_t));

  for (i=0; i<h; i++)
  {    
    for (j=0; j<(w>>3); j++) //1
    {
      AE_L16X4X2_IP(x0, x1, src, -2 * (int)sizeof(ae_int16x4));
      AE_DSEL16X4(x0, x1, x0, x1, dsel0);
      AE_SA16X4X2_IP(x1, x0, aD, dst0);
    }
    AE_SA128POS_FP(aD, dst0);
    src = (const ae_int16x8 *)((uintptr_t)src- (istride-(w>>3)*8)*sizeof(int16_t));
    dst0 = (      ae_int16x8 *)((uintptr_t)dst0+ (ostride-(w>>3)*8)*sizeof(int16_t));
  }  
  __Pragma("no_reorder");
  for (i=0; i<(w&7); i++)
  {
    src  = (const ae_int16x8 *)((int16_t*)inImg + (h - 1)*istride + w - 1 - i);
    dst0 = (      ae_int16x8 *)((int16_t*)outImg + (PADDING*ostride + (ALIGNMENT / 2)) + i);

    for (j=0; j<h; j++)
    {
      ae_int16x4 p;
      AE_L16_XP(p,castxcc(ae_int16,src),-(int)(sizeof(int16_t)*istride));
      AE_S16_0_XP(p,castxcc(ae_int16,dst0),sizeof(int16_t)*ostride);
    }
  }

  outSz->width = w + ALIGNMENT;
  outSz->height = h + PADDING * 2;
#else
  int i,j,w,h,ostride,istride;
  const int16_t * restrict in =(const int16_t *)inImg;
        int16_t * restrict out=(      int16_t *)outImg;
  (void)pScr;
  NASSERT_ALIGN(pScr,ALIGNMENT);
  NASSERT_ALIGN(outImg,ALIGNMENT);
  NASSERT_ALIGN(inImg ,ALIGNMENT);
  imgsize_validate(inSz,2,1);
  imgsize_validate(outSz,2,1);
  ostride = outSz->stride;
  istride = inSz->stride;
  for (i=0; i<PADDING*ostride; i++) *out++=fillColor;
  // flip by 180 degrees
  w = inSz->width;
  h = inSz->height;
  in += (h - 1)*istride;
  for (i = 0; i<h; i++, in -= istride)
  {
    for (j = 0; j<ALIGNMENT / 2; j++) *out++ = fillColor;
    for (j = 0; j<w; j++) *out++ = in[w - 1 - j];
    for (; j<(ostride - ALIGNMENT / 2); j++) *out++ = fillColor;
  }
  for (i = 0; i<(int)(outSz->height - (h + PADDING))*ostride; i++) *out++ = fillColor;
  outSz->width = w + ALIGNMENT;
  outSz->height = h + PADDING * 2;
#endif
} /* imgfast_gs16_rot180_16() */
