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
rotation by 0,90,180,270 with no format conversion
Input:
inImg, inSz  image and its size
outSz        output image size
Output:
outImg       image
--------------------------------------------------------------------------*/
void imgfast_gu8_rot180 (void* pScr, void* outImg, const imgsize_t * restrict outSz, const void* inImg, const imgsize_t * restrict inSz)
{
#if 1
  int i, j, w, h, ostride, istride;
  const ae_int8x16 * restrict src = (const ae_int8x16 *)inImg;
        ae_int8x16 * restrict dst = (      ae_int8x16 *)outImg;

  (void)pScr;
  NASSERT_ALIGN(pScr, ALIGNMENT);
  NASSERT_ALIGN(outImg, ALIGNMENT);
  NASSERT_ALIGN(inImg, ALIGNMENT);

  static const uint8_t ALIGN(ALIGNMENT) dsel_tbl[] = { 0 | (8 << 4), 1 | (9 << 4), 2 | (10 << 4), 3 | (11 << 4), 4 | (12 << 4), 5 | (13 << 4), 6 | (14 << 4), 7 | (15 << 4) };
  ae_int8x8 dsel0 = AE_L8X8_I((const ae_int8x8*)dsel_tbl, 0);
  imgsize_validate(inSz, 1, 1);
  imgsize_validate(outSz, 1, 1);
  w = inSz->width;
  h = inSz->height;
  ostride = outSz->stride;
  istride = inSz->stride;
  ae_valignx2 aD;
  // flip by 180 degrees
  src = (const ae_int8x16 *)((uintptr_t)inImg + ((h - 1)*istride + (w >> 4) * 16)*sizeof(int8_t)) - 1;
  dst = (      ae_int8x16 *)((uintptr_t)outImg + (w & 15)*sizeof(int8_t));
  aD = AE_ZALIGN128();
  for (i = 0; i<h; i++)
  {
    for (j=0; j<(w>>4); j++)
    { 
      ae_int8x8 x0, x1;
      AE_L8X8X2_IP(x0, x1, src, -2 * (int)sizeof(ae_int8x8));
      AE_DSEL8X8(x0, x1, x0, x1, dsel0);
      AE_SA8X8X2_IP(x1, x0, aD, dst);
    }
    AE_SA128POS_FP(aD, dst);

    src = (const ae_int8x16 *)((uintptr_t)src - (istride - (w&~15))*sizeof(uint8_t));
    dst = (      ae_int8x16 *)((uintptr_t)dst + (ostride - (w&~15))*sizeof(uint8_t));
  }
  __Pragma("no_reorder");
  for (i = 0; i<(w&15); i++) 
  {
    src = (const ae_int8x16 *)((int8_t*)inImg + (h - 1)*istride + w - 1 - i);
    dst=(        ae_int8x16 *)((int8_t*)outImg+i) ;

    for (j=0; j<h; j++)
    {
      ae_int8x8 x0;
      AE_L8_XP(x0, castxcc(ae_int8, src), -(istride)*(int)sizeof(uint8_t));
      AE_S8_0_XP(x0, castxcc(ae_int8, dst), (ostride)*sizeof(uint8_t));
    }
  }

#else
  int i, j, w, h, ostride, istride;
  const uint8_t * restrict in = (const uint8_t *)inImg;
  uint8_t * restrict out = (uint8_t *)outImg;
  (void)pScr;
  NASSERT_ALIGN(pScr, ALIGNMENT);
  NASSERT_ALIGN(outImg, ALIGNMENT);
  NASSERT_ALIGN(inImg, ALIGNMENT);
  imgsize_validate(inSz, 1, 1);
  imgsize_validate(outSz, 1, 1);
  w = inSz->width;
  h = inSz->height;
  ostride = outSz->stride;
  istride = inSz->stride;
  in += (h - 1)*istride;
  for (i = 0; i<h; i++, in -= istride)
  {
    for (j = 0; j<w; j++) 
    {
      *out++ = in[w - 1 - j];    
    }
    out += ostride - w;
  }
#endif
} /* imgfast_gu8_rot180() */
