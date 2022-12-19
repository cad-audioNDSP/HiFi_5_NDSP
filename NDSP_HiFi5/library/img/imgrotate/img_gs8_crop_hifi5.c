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
void img_gs8_crop(void *pScr, void* outImg, const imgsize_t * restrict outSz, const void* inImg, const imgsize_t * restrict inSz, int left, int up)
{
#if 1
  const ae_int16x8 * restrict src = (const ae_int16x8 *)inImg;
        ae_int8x16 * restrict dst = (      ae_int8x16 *)outImg;
  int i, j, w, h, ostride, istride;
  ae_valignx2 aD;
  ae_int16x4 x0, x1, x2, x3;
  ae_int8x8  y0, y1;
  ae_int16x4 _64;
  NASSERT_ALIGN(pScr,ALIGNMENT);
  imgsize_validate(inSz, 2, 1);
  imgsize_validate(outSz, 1, 0);
  NASSERT_ALIGN(outImg, 1);        /* output non-aligned, input aligned */
  NASSERT_ALIGN(inImg, ALIGNMENT);  w = outSz->width;
  _64 = 64;
  h = outSz->height;
  ostride = outSz->stride;
  istride = inSz->stride;
  src = (const ae_int16x8 *)((uintptr_t)inImg + (up* istride + left) * sizeof(int16_t));
  for (i = 0; i<h; i++)
  {
    int s_;
    s_ = ((uintptr_t)src & 15) >> 1;
    if (s_) s_ = 8 - s_;
    for (j = 0; j < s_; j++)
    {
      AE_L16_IP(x0, castxcc(ae_int16, src), sizeof(int16_t));
      x0 = AE_MAX16(x0, AE_ZERO16());
      x0 = AE_SRAI16R(x0, 7);
      y0 = AE_MOVINT8X8_FROMINT16X4(x0);
      y0 = AE_MOVINT8X8_FROMINT32X2(AE_XOR32(AE_MOVINT32X2_FROMINT8X8(y0), 0x80808080));
      AE_S8_0_XP(y0, castxcc(ae_int8, dst), sizeof(int8_t));
    }

    aD = AE_ZALIGN128();
    for (j = 0; j<((w - s_) >> 4); j++)  //9/2
    {
      AE_L16X4X2_IP(x0, x1, src, 2 * sizeof(ae_int16x4));
      AE_L16X4X2_IP(x2, x3, src, 2 * sizeof(ae_int16x4));

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
      y0 = AE_MOVINT8X8_FROMINT32X2(AE_XOR32(AE_MOVINT32X2_FROMINT8X8(y0), 0x80808080));
      y1 = AE_MOVINT8X8_FROMINT32X2(AE_XOR32(AE_MOVINT32X2_FROMINT8X8(y1), 0x80808080));

      AE_SA8X8X2_IP(y0, y1, aD, dst);
    }
    AE_SA128POS_FP(aD, dst);
    for (j = 0; j < ((w - s_) & 15); j++)
    {
      AE_L16_IP(x0, castxcc(ae_int16, src), sizeof(int16_t));
      x0 = AE_MAX16(x0, AE_ZERO16());
      x0 = AE_SRAI16R(x0, 7);
      y0 = AE_MOVINT8X8_FROMINT16X4(x0);
      y0 = AE_MOVINT8X8_FROMINT32X2(AE_XOR32(AE_MOVINT32X2_FROMINT8X8(y0), 0x80808080));
      AE_S8_0_XP(y0, castxcc(ae_int8, dst), sizeof(int8_t));
    }     
    src = (const ae_int16x8 *)((uintptr_t)src + (istride - w)*sizeof(int16_t));
    dst = (      ae_int8x16 *)((uintptr_t)dst + (ostride - w)*sizeof(int8_t));
  }
#else
    const int16_t* restrict in;
          uint8_t* restrict out;
    int i,j,w,h,ostride,istride;
    (void)pScr;
    uint8_t tmp;
    NASSERT_ALIGN(pScr,ALIGNMENT);
    w=outSz->width;
    h=outSz->height;
    ostride=outSz->stride;
    istride=inSz->stride;
    imgsize_validate(inSz,2,1);
    imgsize_validate(outSz,1,0);
    NASSERT_ALIGN(outImg,1);        /* output non-aligned, input aligned */
    NASSERT_ALIGN(inImg,ALIGNMENT);

    in =(const int16_t*)inImg;
    out=(      uint8_t*)outImg;
    in+=up*istride+left;
    for (i=0; i<h; i++,in+=istride)
    {
        for (j=0; j<w; j++)
        {
          tmp = (uint8_t)(S_shr_s(S_add_ss(XT_MAX((in[j]), 0), (1 << 6)), 7));
          *out++ = (int8_t)(tmp^128);
        }
        out+=ostride-w;
    }
#endif
} /* img_gs8_crop() */
