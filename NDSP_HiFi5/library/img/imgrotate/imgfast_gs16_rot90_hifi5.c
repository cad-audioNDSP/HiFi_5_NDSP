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
void imgfast_gs16_rot90 (void* pScr, void* outImg, const imgsize_t * restrict outSz, const void* inImg, const imgsize_t * restrict inSz)
{
#if 1
  const ae_int16x8 * restrict src0;
  const ae_int16x8 * restrict src1;
  const ae_int16x8 * restrict src2;
  const ae_int16x8 * restrict src3;
  const ae_int16x8 * restrict src4;
  const ae_int16x8 * restrict src5;
  const ae_int16x8 * restrict src6;
  const ae_int16x8 * restrict src7;
        ae_int16x8 * restrict dst0;
        ae_int16x8 * restrict dst1;
        ae_int16x8 * restrict dst2;
        ae_int16x8 * restrict dst3;
        ae_int16x8 * restrict dst4;
        ae_int16x8 * restrict dst5;
        ae_int16x8 * restrict dst6;
        ae_int16x8 * restrict dst7;

  int i, j, w, h, ostride, istride;
  ae_int16x4 x0, x1, x2, x3;
  ae_int16x4 x4, x5, x6, x7;
  ae_int16x4 y0, y1, y2, y3;
  ae_int16x4 y4, y5, y6, y7;
  ae_int16x4 z0, z1, z2, z3;
  ae_int16x4 z4, z5, z6, z7;
  NASSERT_ALIGN(pScr,ALIGNMENT);
  NASSERT_ALIGN(outImg,ALIGNMENT);
  NASSERT_ALIGN(inImg ,ALIGNMENT);
  static const int16_t ALIGN(ALIGNMENT) dsel_tbl[] = { 5 | (7 << 8), 1 | (3 << 8), 4 | (6 << 8), 0 | (2 << 8),
                                                       7 | (5 << 8), 6 | (4 << 8), 3 | (1 << 8), 2 | (0 << 8)};

  ae_int16x4 dsel0 = AE_L16X4_I((const ae_int16x4*)dsel_tbl, 0);
  ae_int16x4 dsel1 = AE_L16X4_I((const ae_int16x4*)dsel_tbl, sizeof(ae_int16x4));


  w = inSz->height;
  h = inSz->width;
  ostride = outSz->stride;
  istride = inSz->stride;

  for (i=0; i<(w>>3); i++)
  {
    src0 = (const ae_int16x8 *)((uintptr_t)inImg + ((0 + 8 * i) * istride + (h)-(h & 7)) * sizeof(int16_t)) - 1;
    src1 = (const ae_int16x8 *)((uintptr_t)inImg + ((1 + 8 * i) * istride + (h)-(h & 7)) * sizeof(int16_t)) - 1;
    src2 = (const ae_int16x8 *)((uintptr_t)inImg + ((2 + 8 * i) * istride + (h)-(h & 7)) * sizeof(int16_t)) - 1;
    src3 = (const ae_int16x8 *)((uintptr_t)inImg + ((3 + 8 * i) * istride + (h)-(h & 7)) * sizeof(int16_t)) - 1;
    src4 = (const ae_int16x8 *)((uintptr_t)inImg + ((4 + 8 * i) * istride + (h)-(h & 7)) * sizeof(int16_t)) - 1;
    src5 = (const ae_int16x8 *)((uintptr_t)inImg + ((5 + 8 * i) * istride + (h)-(h & 7)) * sizeof(int16_t)) - 1;
    src6 = (const ae_int16x8 *)((uintptr_t)inImg + ((6 + 8 * i) * istride + (h)-(h & 7)) * sizeof(int16_t)) - 1;
    src7 = (const ae_int16x8 *)((uintptr_t)inImg + ((7 + 8 * i) * istride + (h)-(h & 7)) * sizeof(int16_t)) - 1;


    dst0 = (ae_int16x8 *)((uintptr_t)outImg + (((h & 7) + 0)* ostride)*sizeof(int16_t)) + i;
    dst4 = (ae_int16x8 *)((uintptr_t)outImg + (((h & 7) + 4)* ostride)*sizeof(int16_t)) + i;

    __Pragma("loop_count min =1 ")
    for (j=0; j<(h>>3); j++) //8
    {
      AE_L16X4X2_IP(x0, y0, src0, -2 * (int)sizeof(ae_int16x4));
      AE_L16X4X2_IP(x1, y1, src1, -2 * (int)sizeof(ae_int16x4));
      AE_L16X4X2_IP(x2, y2, src2, -2 * (int)sizeof(ae_int16x4));
      AE_L16X4X2_IP(x3, y3, src3, -2 * (int)sizeof(ae_int16x4));
      

     AE_L16X4X2_IP(x4, y4, src4, -2 * (int)sizeof(ae_int16x4));
     AE_L16X4X2_IP(x5, y5, src5, -2 * (int)sizeof(ae_int16x4));
     AE_L16X4X2_IP(x6, y6, src6, -2 * (int)sizeof(ae_int16x4));
     AE_L16X4X2_IP(x7, y7, src7, -2 * (int)sizeof(ae_int16x4));

      AE_DSEL16X4(z0, z2, x0, x1, dsel0);
      AE_DSEL16X4(z1, z3, x2, x3, dsel0);
      AE_DSEL16X4(x0, x1, z2, z3, dsel1);
      AE_DSEL16X4(x2, x3, z0, z1, dsel1);

      AE_DSEL16X4(z0, z2, x4, x5, dsel0);
      AE_DSEL16X4(z1, z3, x6, x7, dsel0);
      AE_DSEL16X4(x4, x5, z2, z3, dsel1);
      AE_DSEL16X4(x6, x7, z0, z1, dsel1);

      AE_S16X4X2_XP(x0, x4, dst4,  ostride*sizeof(int16_t));
      AE_S16X4X2_XP(x1, x5, dst4,  ostride*sizeof(int16_t));
      AE_S16X4X2_XP(x2, x6, dst4,  ostride*sizeof(int16_t));
      AE_S16X4X2_XP(x3, x7, dst4, 5 * ostride*sizeof(int16_t)); 

      AE_DSEL16X4(z4, z6, y0, y1, dsel0);
      AE_DSEL16X4(z5, z7, y2, y3, dsel0);
      AE_DSEL16X4(y0, y1, z6, z7, dsel1);
      AE_DSEL16X4(y2, y3, z4, z5, dsel1);

     
      AE_DSEL16X4(z4, z6, y4, y5, dsel0);
      AE_DSEL16X4(z5, z7, y6, y7, dsel0);
      AE_DSEL16X4(y4, y5, z6, z7, dsel1);
      AE_DSEL16X4(y6, y7, z4, z5, dsel1);


      AE_S16X4X2_XP(y0, y4, dst0, ostride*sizeof(int16_t));
      AE_S16X4X2_XP(y1, y5, dst0, ostride*sizeof(int16_t));
      AE_S16X4X2_XP(y2, y6, dst0, ostride*sizeof(int16_t));
      AE_S16X4X2_XP(y3, y7, dst0, 5 * ostride*sizeof(int16_t));

    }
  }
  __Pragma("no_reorder");
  for (j = 0; j<(h&7); j++)
  {
    src0 = (const ae_int16x8 *)((uintptr_t)inImg + (0 * istride + (h - 1 - j)) * sizeof(int16_t));
    src1 = (const ae_int16x8 *)((uintptr_t)inImg + (1 * istride + (h - 1 - j)) * sizeof(int16_t));
    src2 = (const ae_int16x8 *)((uintptr_t)inImg + (2 * istride + (h - 1 - j)) * sizeof(int16_t));
    src3 = (const ae_int16x8 *)((uintptr_t)inImg + (3 * istride + (h - 1 - j)) * sizeof(int16_t));
    src4 = (const ae_int16x8 *)((uintptr_t)inImg + (4 * istride + (h - 1 - j)) * sizeof(int16_t));
    src5 = (const ae_int16x8 *)((uintptr_t)inImg + (5 * istride + (h - 1 - j)) * sizeof(int16_t));
    src6 = (const ae_int16x8 *)((uintptr_t)inImg + (6 * istride + (h - 1 - j)) * sizeof(int16_t));
    src7 = (const ae_int16x8 *)((uintptr_t)inImg + (7 * istride + (h - 1 - j)) * sizeof(int16_t));

    dst0 = (ae_int16x8 *)((uintptr_t)outImg + (j * ostride)*sizeof(int16_t));  
    for (i = 0; i<(w >> 3); i++) //4
    {
      ae_int16x4 p0, p1, p2, p3;
      ae_int16x4 p4, p5, p6, p7;
      AE_L16_XP(p0, castxcc(ae_int16, src0), (int)(sizeof(int16_t) * 8 * istride));
      AE_L16_XP(p1, castxcc(ae_int16, src1), (int)(sizeof(int16_t) * 8 * istride));
      AE_L16_XP(p2, castxcc(ae_int16, src2), (int)(sizeof(int16_t) * 8 * istride));
      AE_L16_XP(p3, castxcc(ae_int16, src3), (int)(sizeof(int16_t) * 8 * istride));
      AE_L16_XP(p4, castxcc(ae_int16, src4), (int)(sizeof(int16_t) * 8 * istride));
      AE_L16_XP(p5, castxcc(ae_int16, src5), (int)(sizeof(int16_t) * 8 * istride));
      AE_L16_XP(p6, castxcc(ae_int16, src6), (int)(sizeof(int16_t) * 8 * istride));
      AE_L16_XP(p7, castxcc(ae_int16, src7), (int)(sizeof(int16_t) * 8 * istride));


      AE_S16_0_IP(p0, castxcc(ae_int16, dst0), sizeof(int16_t));
      AE_S16_0_IP(p1, castxcc(ae_int16, dst0), sizeof(int16_t));
      AE_S16_0_IP(p2, castxcc(ae_int16, dst0), sizeof(int16_t));
      AE_S16_0_IP(p3, castxcc(ae_int16, dst0), sizeof(int16_t));
      AE_S16_0_IP(p4, castxcc(ae_int16, dst0), sizeof(int16_t));
      AE_S16_0_IP(p5, castxcc(ae_int16, dst0), sizeof(int16_t));
      AE_S16_0_IP(p6, castxcc(ae_int16, dst0), sizeof(int16_t));
      AE_S16_0_IP(p7, castxcc(ae_int16, dst0), sizeof(int16_t));

    }
  }
  __Pragma("no_reorder");
  for (i=0; i<(w&7); i++)
  {
    
    src0 = (const ae_int16x8 *)((uintptr_t)inImg + ((w - 1 - i)* istride + (h)-(h & 7)) * sizeof(int16_t)) - 1;
    dst0 = (ae_int16x8 *)((uintptr_t)outImg + (((h & 7) + 0)* ostride  + (w - 1 - i))*sizeof(int16_t));
    dst1 = (ae_int16x8 *)((uintptr_t)outImg + (((h & 7) + 1)* ostride  + (w - 1 - i))*sizeof(int16_t));
    dst2 = (ae_int16x8 *)((uintptr_t)outImg + (((h & 7) + 2)* ostride  + (w - 1 - i))*sizeof(int16_t));
    dst3 = (ae_int16x8 *)((uintptr_t)outImg + (((h & 7) + 3)* ostride  + (w - 1 - i))*sizeof(int16_t));
    dst4 = (ae_int16x8 *)((uintptr_t)outImg + (((h & 7) + 4)* ostride + (w - 1 - i))*sizeof(int16_t));
    dst5 = (ae_int16x8 *)((uintptr_t)outImg + (((h & 7) + 5)* ostride + (w - 1 - i))*sizeof(int16_t));
    dst6 = (ae_int16x8 *)((uintptr_t)outImg + (((h & 7) + 6)* ostride + (w - 1 - i))*sizeof(int16_t));
    dst7 = (ae_int16x8 *)((uintptr_t)outImg + (((h & 7) + 7)* ostride + (w - 1 - i))*sizeof(int16_t));

    for (j = 0; j<(h>>3); j++)
    {
      AE_L16X4X2_XP(x0, x1, castxcc(ae_int16x8, src0), -2 * (int)sizeof(ae_int16x4));
      AE_S16_0_XP(x0, castxcc(ae_int16, dst4), 8 * ostride*sizeof(int16_t));
      x0 = AE_SEL16_4321(x0, x0);
      AE_S16_0_XP(x0, castxcc(ae_int16, dst5), 8 * ostride*sizeof(int16_t));
      x0 = AE_SEL16_4321(x0, x0);
      AE_S16_0_XP(x0, castxcc(ae_int16, dst6), 8 * ostride*sizeof(int16_t));
      x0 = AE_SEL16_4321(x0, x0);
      AE_S16_0_XP(x0, castxcc(ae_int16, dst7), 8 * ostride*sizeof(int16_t));

      AE_S16_0_XP(x1, castxcc(ae_int16, dst0), 8 * ostride*sizeof(int16_t));
      x1 = AE_SEL16_4321(x1, x1);
      AE_S16_0_XP(x1, castxcc(ae_int16, dst1), 8 * ostride*sizeof(int16_t));
      x1 = AE_SEL16_4321(x1, x1);
      AE_S16_0_XP(x1, castxcc(ae_int16, dst2), 8 * ostride*sizeof(int16_t));
      x1 = AE_SEL16_4321(x1, x1);
      AE_S16_0_XP(x1, castxcc(ae_int16, dst3), 8 * ostride*sizeof(int16_t));

    }
    src0 = (const ae_int16x8 *)((uintptr_t)inImg + ((w - 1 - i) * istride + (h - 1)) * sizeof(int16_t));
    dst0 = (      ae_int16x8 *)((uintptr_t)outImg + (w - 1 - i) *sizeof(int16_t));
    for (j = 0; j<(h & 7); j++)
    {
      AE_L16_IP(x0, castxcc(ae_int16, src0), -(int)sizeof(int16_t));
      AE_S16_0_XP(x0, castxcc(ae_int16, dst0), ostride*sizeof(int16_t));
    }
  }
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

  // copy
  w = inSz->height;
  h = inSz->width;
  ostride = outSz->stride;
  istride = inSz->stride;
  for (i = 0; i<h; i++)
  {
    for (j = 0; j<w; j++)    out[i*ostride + j] = in[j*istride + (h - 1 - i)];
  }
  NASSERT(outSz->width ==w);
  NASSERT(outSz->height==h);
#endif
} /* imgfast_gs16_rot90() */
