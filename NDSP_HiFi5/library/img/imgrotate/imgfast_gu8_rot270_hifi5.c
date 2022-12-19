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
void imgfast_gu8_rot270 (void* pScr, void* outImg, const imgsize_t * restrict outSz, const void* inImg, const imgsize_t * restrict inSz)
{
#if 1
  const ae_int8x16 * restrict src0;
  const ae_int8x16 * restrict src1;
  const ae_int8x16 * restrict src2;
  const ae_int8x16 * restrict src3;
  const ae_int8x16 * restrict src4;
  const ae_int8x16 * restrict src5;
  const ae_int8x16 * restrict src6;
  const ae_int8x16 * restrict src7;
        ae_int8x8 * restrict dst0;

 
  int i, j, w, h, ostride, istride;
  ae_int8x8 x0, x1, x2, x3;
  ae_int8x8 x4, x5, x6, x7;
  ae_int8x8 y0, y1, y2, y3;
  ae_int8x8 y4, y5, y6, y7;
  ae_int8x8 z0, z1, z2, z3;
  ae_int8x8 z4, z5, z6, z7;
  NASSERT_ALIGN(pScr,ALIGNMENT);
  NASSERT_ALIGN(outImg,ALIGNMENT);
  NASSERT_ALIGN(inImg ,ALIGNMENT);

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

  for (i=0; i<(w>>3); i++)
  {
    src0 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 8 * i - 0) * istride) * sizeof(int8_t));
    src1 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 8 * i - 1) * istride) * sizeof(int8_t));
    src2 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 8 * i - 2) * istride) * sizeof(int8_t));
    src3 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 8 * i - 3) * istride) * sizeof(int8_t));
    src4 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 8 * i - 4) * istride) * sizeof(int8_t));
    src5 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 8 * i - 5) * istride) * sizeof(int8_t));
    src6 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 8 * i - 6) * istride) * sizeof(int8_t));
    src7 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 8 * i - 7) * istride) * sizeof(int8_t));


    dst0 = (ae_int8x8 *)((uintptr_t)outImg + (0 * ostride) * sizeof(int8_t)) + i;
    
    for (j=0; j<(h>>4); j++) 
    {
      AE_L8X8X2_XP(x0, z0, src0, 2 * (int)sizeof(ae_int8x8));
      AE_L8X8X2_XP(x1, z1, src1, 2 * (int)sizeof(ae_int8x8));
      AE_L8X8X2_XP(x2, z2, src2, 2 * (int)sizeof(ae_int8x8));
      AE_L8X8X2_XP(x3, z3, src3, 2 * (int)sizeof(ae_int8x8));
      AE_L8X8X2_XP(x4, z4, src4, 2 * (int)sizeof(ae_int8x8));
      AE_L8X8X2_XP(x5, z5, src5, 2 * (int)sizeof(ae_int8x8));
      AE_L8X8X2_XP(x6, z6, src6, 2 * (int)sizeof(ae_int8x8));
      AE_L8X8X2_XP(x7, z7, src7, 2 * (int)sizeof(ae_int8x8));        

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
       
      AE_S8X8_XP(y1, dst0, ostride*sizeof(int8_t));
      AE_S8X8_XP(y0, dst0, ostride*sizeof(int8_t));
      AE_S8X8_XP(y3, dst0, ostride*sizeof(int8_t));
      AE_S8X8_XP(y2, dst0, ostride*sizeof(int8_t));
      AE_S8X8_XP(y5, dst0, ostride*sizeof(int8_t));
      AE_S8X8_XP(y4, dst0, ostride*sizeof(int8_t));
      AE_S8X8_XP(y7, dst0, ostride*sizeof(int8_t));
      AE_S8X8_XP(y6, dst0, ostride*sizeof(int8_t));

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

      AE_S8X8_XP(y1, dst0, ostride*sizeof(int8_t));
      AE_S8X8_XP(y0, dst0, ostride*sizeof(int8_t));
      AE_S8X8_XP(y3, dst0, ostride*sizeof(int8_t));
      AE_S8X8_XP(y2, dst0, ostride*sizeof(int8_t));
      AE_S8X8_XP(y5, dst0, ostride*sizeof(int8_t));
      AE_S8X8_XP(y4, dst0, ostride*sizeof(int8_t));
      AE_S8X8_XP(y7, dst0, ostride*sizeof(int8_t));
      AE_S8X8_XP(y6, dst0, ostride*sizeof(int8_t));

    }
  }
  __Pragma("no_reorder");
  for (j = 0; j<(h & 15); j++)
  {
     src0 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 0) * istride + (h - 1 - j)) * sizeof(int8_t));
     src1 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 1) * istride + (h - 1 - j)) * sizeof(int8_t));
     src2 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 2) * istride + (h - 1 - j)) * sizeof(int8_t));
     src3 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 3) * istride + (h - 1 - j)) * sizeof(int8_t));
     src4 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 4) * istride + (h - 1 - j)) * sizeof(int8_t));
     src5 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 5) * istride + (h - 1 - j)) * sizeof(int8_t));
     src6 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 6) * istride + (h - 1 - j)) * sizeof(int8_t));
     src7 = (const ae_int8x16 *)((uintptr_t)inImg + ((w - 1 - 7) * istride + (h - 1 - j)) * sizeof(int8_t));
    
     dst0 = (ae_int8x8 *)((uintptr_t)outImg + ((h - 1 - j) * ostride)*sizeof(int8_t));
     for (i = 0; i<(w >> 3); i++) 
     {
       ae_int8x8 p0, p1, p2, p3;
       ae_int8x8 p4, p5, p6, p7;
       AE_L8_XP(p0, castxcc(ae_int8, src0), -(int)(sizeof(int8_t) * 8 * istride));
       AE_L8_XP(p1, castxcc(ae_int8, src1), -(int)(sizeof(int8_t) * 8 * istride));
       AE_L8_XP(p2, castxcc(ae_int8, src2), -(int)(sizeof(int8_t) * 8 * istride));
       AE_L8_XP(p3, castxcc(ae_int8, src3), -(int)(sizeof(int8_t) * 8 * istride));
       AE_L8_XP(p4, castxcc(ae_int8, src4), -(int)(sizeof(int8_t) * 8 * istride));
       AE_L8_XP(p5, castxcc(ae_int8, src5), -(int)(sizeof(int8_t) * 8 * istride));
       AE_L8_XP(p6, castxcc(ae_int8, src6), -(int)(sizeof(int8_t) * 8 * istride));
       AE_L8_XP(p7, castxcc(ae_int8, src7), -(int)(sizeof(int8_t) * 8 * istride));
    
    
       AE_S8_0_IP(p0, castxcc(ae_int8, dst0), sizeof(int8_t));
       AE_S8_0_IP(p1, castxcc(ae_int8, dst0), sizeof(int8_t));
       AE_S8_0_IP(p2, castxcc(ae_int8, dst0), sizeof(int8_t));
       AE_S8_0_IP(p3, castxcc(ae_int8, dst0), sizeof(int8_t));
       AE_S8_0_IP(p4, castxcc(ae_int8, dst0), sizeof(int8_t));
       AE_S8_0_IP(p5, castxcc(ae_int8, dst0), sizeof(int8_t));
       AE_S8_0_IP(p6, castxcc(ae_int8, dst0), sizeof(int8_t));
       AE_S8_0_IP(p7, castxcc(ae_int8, dst0), sizeof(int8_t));
    
     }
  }
  __Pragma("no_reorder");
  for (i = 0; i<(w & 7); i++)
  {

    src0 = (const ae_int8x16 *)((uintptr_t)inImg + (i * istride) * sizeof(int8_t) );
    dst0 = (      ae_int8x8 *)((uintptr_t)outImg + ((7)* ostride + (w - 1 - i))*sizeof(int8_t) );
    for (j = 0; j<(h >> 3); j++)
    {
      AE_L8X8_XP(x0, castxcc(ae_int8x8, src0), 1 * (int)sizeof(ae_int8x8));
      AE_S8_0_XP(x0, castxcc(ae_int8, dst0), -ostride*sizeof(int8_t));
      x0 = AE_SEL8X8I(x0, x0, 19); //AE_SELI_8B_ROTATE_RIGHT_1
      AE_S8_0_XP(x0, castxcc(ae_int8, dst0), -ostride*sizeof(int8_t));
      x0 = AE_SEL8X8I(x0, x0, 19); //AE_SELI_8B_ROTATE_RIGHT_1
      AE_S8_0_XP(x0, castxcc(ae_int8, dst0), -ostride*sizeof(int8_t));
      x0 = AE_SEL8X8I(x0, x0, 19); //AE_SELI_8B_ROTATE_RIGHT_1
      AE_S8_0_XP(x0, castxcc(ae_int8, dst0), -ostride*sizeof(int8_t));
      x0 = AE_SEL8X8I(x0, x0, 19); //AE_SELI_8B_ROTATE_RIGHT_1
      AE_S8_0_XP(x0, castxcc(ae_int8, dst0), -ostride*sizeof(int8_t));
      x0 = AE_SEL8X8I(x0, x0, 19); //AE_SELI_8B_ROTATE_RIGHT_1
      AE_S8_0_XP(x0, castxcc(ae_int8, dst0), -ostride*sizeof(int8_t));
      x0 = AE_SEL8X8I(x0, x0, 19); //AE_SELI_8B_ROTATE_RIGHT_1
      AE_S8_0_XP(x0, castxcc(ae_int8, dst0), -ostride*sizeof(int8_t));
      x0 = AE_SEL8X8I(x0, x0, 19); //AE_SELI_8B_ROTATE_RIGHT_1
      AE_S8_0_XP(x0, castxcc(ae_int8, dst0), 15*ostride*sizeof(int8_t));
   
    }
    src0 = (const ae_int8x16 *)((uintptr_t)inImg + (i* istride + (h - 1)) * sizeof(int8_t));
    dst0 = (      ae_int8x8 *)((uintptr_t)outImg +  ((h - 1)* ostride + (w - 1 - i)) *sizeof(int8_t));
    for (j = 0; j<(h & 7); j++)
    {
      AE_L8_IP(x0, castxcc(ae_int8, src0), -(int)sizeof(int8_t));
      AE_S8_0_XP(x0, castxcc(ae_int8, dst0), -ostride*sizeof(int8_t));
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
  w = inSz->height;
  h = inSz->width;
  ostride = outSz->stride;
  istride = inSz->stride;
  for (i = 0; i<h; i++)
  {
    for (j = 0; j<w; j++)
    {
      out[i*ostride + j] = in[(w - 1 - j)*istride + i];

    }
  }
#endif
} /* imgfast_gu8_rot270() */
