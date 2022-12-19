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
void imgfast_gu8_rot0 (void* pScr, void* outImg, const imgsize_t * restrict outSz, const void* inImg, const imgsize_t * restrict inSz)
{
#if 1
  int i, j, w, h, ostride, istride;
  const ae_int8x16 * restrict src = (const ae_int8x16 *)inImg;
        ae_int8x16 * restrict dst = (      ae_int8x16 *)outImg;

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

 
  for (i = 0; i<h; i++)
  {
    for (j=0; j<(w>>4); j++)
    { 
      ae_int8x8 x0, x1;
      AE_L8X8X2_IP(x0, x1, src, 2 * sizeof(ae_int8x8));
      AE_S8X8X2_IP(x0, x1, dst, 2 * sizeof(ae_int8x8));
    }
    src = (const ae_int8x16 *)((uintptr_t)src + (istride - (w&~15))*sizeof(uint8_t));
    dst = (      ae_int8x16 *)((uintptr_t)dst + (ostride - (w&~15))*sizeof(uint8_t));
  }
  __Pragma("no_reorder");
  for (i = 0; i<(w&15); i++) 
  {
    src = (const ae_int8x16 *)((uint8_t*)inImg + w - 1 - i);
    dst = (      ae_int8x16 *)((uint8_t*)outImg + w - 1 - i);
    for (j=0; j<h; j++)
    {
      ae_int8x8 x0;
      AE_L8_XP(x0, castxcc(ae_int8, src), (istride)*sizeof(uint8_t));
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
  for (i = 0; i<h; i++, in += istride)
  {
    for (j = 0; j<w; j++) *out++ = in[j];
    out += ostride - w;
  }
#endif
} /* imgfast_gu8_rot0() */
