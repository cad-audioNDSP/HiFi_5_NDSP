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
void imgfast_gs16_crop(void *pScr, void* outImg, const imgsize_t * restrict outSz, const void* inImg, const imgsize_t * restrict inSz, int left, int up)
{
#if 1
  const ae_int16x8 * restrict src = (const ae_int16x8 *)inImg;
        ae_int16x8 * restrict dst = (      ae_int16x8 *)outImg;
  int i, j, w, h, ostride, istride;
  ae_valignx2 aS;
  ae_int16x4 x0, x1;
  NASSERT_ALIGN(pScr,ALIGNMENT);
  NASSERT_ALIGN(outImg,ALIGNMENT);
  NASSERT_ALIGN(inImg ,ALIGNMENT);

  w = outSz->width;
  h = outSz->height;
  ostride = outSz->stride;
  istride = inSz->stride;
  imgsize_validate(inSz,2,1);
  imgsize_validate(outSz,2,1);

  src = (const ae_int16x8 *)((uintptr_t)inImg + (up* istride + left) * sizeof(int16_t));
  for (i=0; i<h; i++)
  {    
    aS = AE_LA128_PP(src);
    for (j=0; j<(w>>3); j++) //2/2
    {
      AE_LA16X4X2_IP(x0, x1, aS, src);
      x0 = AE_MAX16(x0, AE_ZERO16());
      x1 = AE_MAX16(x1, AE_ZERO16());
      AE_S16X4X2_IP(x0, x1, dst, 2 * sizeof(ae_int16x4));
    }
    src = (const ae_int16x8 *)((uintptr_t)src + (istride - (w >> 3) * 8)*sizeof(int16_t));
    dst = (      ae_int16x8 *)((uintptr_t)dst + (ostride - (w >> 3) * 8)*sizeof(int16_t));
  }
  for (i=0; i<(w&7); i++)
  {
    src = (const ae_int16x8 *)((uintptr_t)inImg + (up* istride + left + (w>>3)*8 + i) * sizeof(int16_t));
    dst = (        ae_int16x8 *)((int16_t*)outImg  + (w>>3)*8 + i) ;
    for (j=0; j<h; j++)
    {
      ae_int16x4 p;
      AE_L16_XP(p,castxcc(ae_int16,src),(sizeof(int16_t)*istride));
      p = AE_MAX16(p, AE_ZERO16());
      AE_S16_0_XP(p,castxcc(ae_int16,dst),sizeof(int16_t)*ostride);
    }
  }
#else
    const int16_t* restrict in;
          int16_t* restrict out;
    int i,j,w,h,ostride,istride;
    (void)pScr;
    NASSERT_ALIGN(pScr,ALIGNMENT);
    w=outSz->width;
    h=outSz->height;
    ostride=outSz->stride;
    istride=inSz->stride;
    NASSERT_ALIGN(outImg,ALIGNMENT);
    NASSERT_ALIGN(inImg,ALIGNMENT);
    imgsize_validate(inSz,2,1);
    imgsize_validate(outSz,2,1);

    in =(const int16_t*)inImg;
    out=(      int16_t*)outImg;
    in+=up*istride+left;
    for (i=0; i<h; i++,in+=istride)
    {
        for (j=0; j<w      ; j++) *out++=XT_MAX(0,in[j]);
        out+=ostride-w;
    }
#endif
} /* imgfast_gs16_crop() */
