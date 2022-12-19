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
image copying (16-bit aligned temporary->8bit aligned)
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
#if 0
{
    int m,n,w,h,
        wout=(int)out->width, win=(int)out->width,
        hout=(int)out->height,hin=(int)out->height,
        istride=in->stride,ostride=out->stride;
    NASSERT(wout <= win );
    NASSERT(hout <= hin );
    h=hout;
    w=wout;
    imgsize_validate(in,  2, 1);
    imgsize_validate(out, 2, 1);
    for (n=0; n<h; n++)
    for (m=0; m<w; m++)
    {
        int16_t p=((const int16_t*)inImg)[m+n*istride];
        p=MAX(0,p);
        ((int16_t*)outImg)[m+n*ostride]=(int16_t)p;
    }
}
#else
{
    const ae_int16x8 * restrict pIn;
          ae_int16x8 * restrict pOut;
    int m,n,w,h,
        wout=(int)out->width, win=(int)out->width,
        hout=(int)out->height,hin=(int)out->height,
        istride=in->stride,ostride=out->stride;
    NASSERT(wout <= win );
    NASSERT(hout <= hin );
    (void)win,(void)hin;
    h=hout;
    w=wout;
    imgsize_validate(in,  2, 1);
    imgsize_validate(out, 2, 1);
    NASSERT(w<=ostride);
    if (w<=0 || h<=0) return;
    pIn = (const ae_int16x8*)(inImg);
    pOut = (ae_int16x8*)(outImg);
    NASSERT_ALIGN(pIn,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(pOut,HIFI_SIMD_WIDTH);
    for (n = 0; n < h; n++)
    {
        pIn = (const ae_int16x8*)((int16_t*)inImg + n * istride);
        pOut = (ae_int16x8*)((int16_t*)outImg + n * ostride);
        for (m = 0; m < (w >> 3); m++)
        {
            ae_int16x4 p0, p1;
            AE_L16X4X2_IP(p0, p1, pIn, sizeof(ae_int16x8));
            p0 = AE_MAX16(p0, 0);
            p1 = AE_MAX16(p1, 0);
            AE_S16X4X2_IP(p0, p1, pOut, sizeof(ae_int16x8));
        }
#if defined(AE_LAV16X4X2_XP)&defined(AE_SAV8X8X2_XP)
        if (w & 7)
        {
            ae_valignx2 alIn, alOut;
            ae_int16x4 p0, p1;
            alIn = AE_LA128_PP(pIn);
            alOut = AE_ZALIGN128();
            AE_LAV16X4X2_XP(p0, p1, alIn, pIn, sizeof(ae_int16) * (w & 7));
            p0 = AE_MAX16(p0, 0);
            p1 = AE_MAX16(p1, 0);
            AE_SAV16X4X2_XP(p0, p1, alOut, pOut, sizeof(ae_int16) * (w & 7));
            AE_SA128POS_FP(alOut, pOut);
        }

#else
        if (w & 7)
        {
            for (m = 0; m < (w & 7); ++m)
            {
                ae_int16x4 p0;
                AE_L16_IP(p0, castxcc(ae_int16,pIn), sizeof(ae_int16));
                p0 = AE_MAX16(p0, 0);
                AE_S16_0_IP(p0, castxcc(ae_int16, pOut), sizeof(ae_int16));
            }
        }
#endif
    }
}
#endif


static size_t getScratchSize(const imgsize_t* in,const imgsize_t* out)
{
    (void)in; (void)out;
    imgsize_validate(in,  2, 1);
    imgsize_validate(out, 2, 1);
    return 0;
}
const imgcopy_api imgcopy_api_temp_gs16fast={getScratchSize,copy};
