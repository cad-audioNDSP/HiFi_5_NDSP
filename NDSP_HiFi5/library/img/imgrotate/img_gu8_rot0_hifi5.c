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
void img_gu8_rot0 (void* pScr, void* outImg, const imgsize_t * restrict outSz, const void* inImg, const imgsize_t * restrict inSz)
{
#if 1
  const ae_int8x16 * restrict src = (const ae_int8x16 *)inImg;
        ae_int8x16 * restrict dst = (      ae_int8x16 *)outImg;

  int i,j,w,h,istride,ostride;
  int sh;
  (void)pScr;
  ae_int8x8 x0, x1;
  ae_valignx2 aD, aS;
  NASSERT_ALIGN(pScr, ALIGNMENT);
  imgsize_validate(inSz, 1, 0);
  imgsize_validate(outSz, 1, 0);
  NASSERT_ALIGN(outImg, 1);   /* both input and outputs nonaligned */
  NASSERT_ALIGN(inImg, 1);

  // copy
  w=inSz->width;
  h=inSz->height;
  ostride = outSz->stride;
  istride = inSz->stride;

  aD = AE_ZALIGN128();
  if (w >= 16)
  {
    for (i=0; i<h; i++)
    { 
      src = (const ae_int8x16 *)((int8_t*)inImg + i*istride);
      dst = (      ae_int8x16 *)((int8_t*)outImg + i*ostride);
      aS = AE_LA128_PP(src);
      /* align input */
      sh = ((uintptr_t)src & 15) ; 
    
      AE_LA8X8X2_IP(x0, x1, aS, src);
      AE_SA8X8X2_IP(x0, x1, aD, dst);
      AE_SA128POS_FP(aD, dst);
      src = (const ae_int8x16*)((int8_t*)src - sh);
      dst = (      ae_int8x16*)((int8_t*)dst - sh);
      for (j = 0; j<((w >> 4) - 1); j++)
      {
        AE_L8X8X2_IP(x0, x1, src, 2 * sizeof(ae_int8x8));
        AE_SA8X8X2_IP(x0, x1, aD, dst);
      }
      AE_SA128POS_FP(aD, dst);
      for (j = 0; j < sh; j++)
      {
        ae_int8x8 p;
        AE_L8_IP(p, castxcc(ae_int8, src), sizeof(int8_t));
        AE_S8_0_IP(p, castxcc(ae_int8, dst), sizeof(int8_t));
      }
    }
  }  
  __Pragma("no_reorder");
  for (i = 0; i<(w & 15); i++)
  {
    src = (const ae_int8x16 *)((int8_t*)inImg + w - 1 - i);
    dst = (      ae_int8x16 *)((int8_t*)outImg + w - 1 - i);
    for (j = 0; j<h; j++)
    {
      ae_int8x8 p;
      AE_L8_XP(p, castxcc(ae_int8, src), sizeof(int8_t)*istride);
      AE_S8_0_XP(p, castxcc(ae_int8, dst), sizeof(int8_t)*ostride);
    }
  }

#else
  int i,j,w,h, istride,ostride;
  const uint8_t * restrict in =(const uint8_t *)inImg;
        uint8_t * restrict out=(      uint8_t *)outImg;
  (void)pScr;
  NASSERT_ALIGN(pScr,ALIGNMENT);
  imgsize_validate(inSz ,1,0);
  imgsize_validate(outSz,1,0);
  NASSERT_ALIGN(outImg,1);   /* both input and outputs nonaligned */ 
  NASSERT_ALIGN(inImg, 1);
  // copy
  w = inSz->width;
  h = inSz->height;
  istride = inSz->stride;
  ostride = outSz->stride;
  for (i = 0; i<h; i++, in += istride)
  {
    for (j = 0; j<w; j++) *out++ = in[j];
    out += ostride - w;
  }
#endif
} /* img_gu8_rot0() */
