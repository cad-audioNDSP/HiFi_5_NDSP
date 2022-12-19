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
#include "img_common.h"
#include "imgcopy_common.h"

#define MAX(a,b)        ((a)>(b) ? (a) : (b))
#define MIN(a,b)        ((a)<(b) ? (a) : (b))

/*-------------------------------------------------------------------------
image copying (8bit non generic->16-bit aligned temporary)
Input:
in  - input size
out - output size
inImg - input imgae
Output:
outImg - output image
Temporary:
pScr
-------------------------------------------------------------------------*/
static void copy(void* pScr, void* restrict outImg,const void* restrict inImg,const imgsize_t* restrict in,const imgsize_t* restrict out)
{
    const uint8_t * restrict src;
          int16_t * restrict dst;
    int m, n,
        w = (int)in->width,
        h = (int)in->height,
        istride = in->stride,
        ostride = out->stride;
    ae_valignx2 alIn;
    
    (void)pScr;
    NASSERT(out->width >= in->width);
    NASSERT(out->height >= in->height);
    imgsize_validate(in, 1, 0);
    imgsize_validate(out, 2, 1);
    w = (w + 15) & ~15;
    for (n = 0; n < h; n++)
    {
        src  = ((const uint8_t *)inImg + n*istride);
        dst  = ((      int16_t *)outImg + n*ostride);
        alIn = AE_LA128_PP(src);
        for (m = 0; m < (w >> 4); ++m)
        {
            ae_int8x8 s0, s1;
            ae_int16x4 p0, p1, p2, p3;
            AE_LA8X8X2_IP(s0, s1, alIn, castxcc(ae_int8x16, src));
            AE_CVTI16X4X2F8U(p0, p1, s0, 7);
            AE_CVTI16X4X2F8U(p2, p3, s1, 7);
            AE_S16X4X2_IP(p0, p1, castxcc(ae_int16x8, dst), sizeof(ae_int16x8));
            AE_S16X4X2_IP(p2, p3, castxcc(ae_int16x8, dst), sizeof(ae_int16x8));
        }
    }
}

static size_t getScratchSize(const imgsize_t* in,const imgsize_t* out)
{
    (void)in; (void)out;
    imgsize_validate(in,  1, 0);
    imgsize_validate(out, 2, 1);
    return 0;
}
const imgcopy_api imgcopy_api_gu8_temp={getScratchSize,copy};
