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
 out-of-place cropping with conversion from signed 16-bit to destination format
    Input:
    left,up    left/up corner
    inImg,inSz input image in intermediate format
    outSz      output size
    Output:
    outImg     output image
 --------------------------------------------------------------------------*/
void imgfast_gu8_crop(void *pScr, void* outImg, const imgsize_t * restrict outSz, const void* inImg, const imgsize_t * restrict inSz, int left, int up)
{
  const ae_int16x8 * restrict src = (const ae_int16x8 *)inImg;
        ae_int8x16 * restrict dst = (      ae_int8x16  *)outImg;

  int i, j, w, h, ostride, istride;
  ae_valignx2 aS;
  ae_int16x4 x0, x1, x2, x3; 
  ae_int8x8  y0, y1;
  ae_int16x4 _64;
  (void)pScr;
  NASSERT_ALIGN(pScr,ALIGNMENT);
  w=outSz->width;
  h=outSz->height;
  ostride=outSz->stride;
  istride=inSz->stride;
  NASSERT_ALIGN(outImg,ALIGNMENT);
  NASSERT_ALIGN(inImg,ALIGNMENT);
  imgsize_validate(inSz,2,1);
  imgsize_validate(outSz,1,1);
  _64 = 64;
 
  src = (const ae_int16x8 *)((uintptr_t)inImg + (up* istride + left) * sizeof(int16_t));
  for (i = 0; i<h; i++)
  {
    aS = AE_LA128_PP(src);
    for (j = 0; j<(w >> 4); j++) //15/4
    {
      AE_LA16X4X2_IP(x0, x1, aS, src);
      AE_LA16X4X2_IP(x2, x3, aS, src);

      x0 = AE_MAX16(x0, AE_ZERO16());
      x1 = AE_MAX16(x1, AE_ZERO16());
      x2 = AE_MAX16(x2, AE_ZERO16());
      x3 = AE_MAX16(x3, AE_ZERO16());

      x0 = AE_SRAI16R(x0, 7);
      x1 = AE_SRAI16R(x1, 7);
      x2 = AE_SRAI16R(x2, 7);
      x3 = AE_SRAI16R(x3, 7);

      y0 = AE_SEL8X8I(AE_MOVINT8X8_FROMINT16X4(x0), AE_MOVINT8X8_FROMINT16X4(x1), 25); //AE_SELI_8B_EXTRACT_EVEN
      y1 = AE_SEL8X8I(AE_MOVINT8X8_FROMINT16X4(x2), AE_MOVINT8X8_FROMINT16X4(x3), 25); //AE_SELI_8B_EXTRACT_EVEN
      AE_S8X8X2_IP(y0, y1, dst, 2 * sizeof(ae_int8x8));
    }
    src = (const ae_int16x8 *)((uintptr_t)src + (istride - (w >> 4) * 16)*sizeof(int16_t));
    dst = (      ae_int8x16 *)((uintptr_t)dst + (ostride - (w >> 4) * 16)*sizeof(int8_t));

  }
  __Pragma("no_reorder");
  
  for (i = 0; i<(w & 15); i++)
  {
    src = (const ae_int16x8 *)((uintptr_t)inImg + (up* istride + left + (w >> 4) * 16 + i) * sizeof(int16_t));
    dst = (      ae_int8x16  *)((int8_t*)outImg + (w >> 4) * 16 + i);
    for (j = 0; j<h; j++)
    {
      AE_L16_XP(x0, castxcc(ae_int16, src), (sizeof(int16_t)*istride));
      x0 = AE_MAX16(x0, AE_ZERO16());
      x0 = AE_SRAI16R(x0, 7);
      y0 = AE_MOVINT8X8_FROMINT16X4(x0);
      AE_S8_0_XP(y0, castxcc(ae_int8, dst), ostride);
    }
  }
} /* imgfast_gu8_crop() */
