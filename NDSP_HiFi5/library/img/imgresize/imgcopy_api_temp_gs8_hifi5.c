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
image copying (16-bit aligned temporary->8bit generic)
Input:
in  - input size
out - output size
inImg - input image
Output:
outImg - output image
Temporary:
pScr
-------------------------------------------------------------------------*/
static void copy(void* pScr, void* restrict outImg,const void* restrict inImg,const imgsize_t* restrict in,const imgsize_t* restrict out)
{
    const int16_t* restrict pIn;
    int8_t* restrict pOut;
    ae_valignx2 alOut;
    ae_int16x4 _256 = 256;
    ae_int16x4 _16384 = 16384;
    ae_int16x4 _0 = 0;
    ae_int16x4 _m256 = -256;
    int m, n, w, h,
        wout = (int)out->width,
        win = (int)out->width,
        hout = (int)out->height,
        hin = (int)out->height,
        istride = in->stride,
        ostride = out->stride;
    (void)hin, (void)win;
    NASSERT(wout <= win);
    NASSERT(hout <= hin);
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 1, 0);
    h = hout;
    alOut = AE_ZALIGN128();
    w = wout;

#if defined(AE_LAV16X4X2_XP)&defined(AE_SAV8X8X2_XP)
    int last0, last1;
    last0 = (w & 15) > 7 ? 8 : w & 15;
    last1 = (w & 15) > 7 ? (w & 15) - 8 : 0;
#endif

    for (n = 0; n < h; n++)
    {
        pIn = ((const int16_t*)inImg + n * istride);
        pOut = ((int8_t*)outImg + n * ostride);
       

        for (m = 0; m < (w >> 4); m++)
        {
            ae_int16x4 p0, p1, p2, p3;
            ae_int8x8 y0, y1;
            AE_L16X4X2_IP(p0, p1, castxcc(ae_int16x8, pIn), sizeof(ae_int16x8));
            AE_L16X4X2_IP(p2, p3, castxcc(ae_int16x8, pIn), sizeof(ae_int16x8));
        
            p0 = AE_MAX16(p0, _0);
            p1 = AE_MAX16(p1, _0);
            p2 = AE_MAX16(p2, _0);
            p3 = AE_MAX16(p3, _0);
       
            p0 = AE_MULFD16X16X4RAS(p0, _16384, _256, _m256);
            p1 = AE_MULFD16X16X4RAS(p1, _16384, _256, _m256);
            p2 = AE_MULFD16X16X4RAS(p2, _16384, _256, _m256);
            p3 = AE_MULFD16X16X4RAS(p3, _16384, _256, _m256);
       
            y0 = AE_SAT8X8X16(p0, p1);
            y1 = AE_SAT8X8X16(p2, p3);
       
            AE_SA8X8X2_IP(y0, y1, alOut, castxcc(ae_int8x16, pOut));
        }
       
#if defined(AE_LAV16X4X2_XP)&defined(AE_SAV8X8X2_XP)
        if (w & 15)
        {
            ae_int16x4 p0, p1, p2, p3;
            ae_int8x8 y0, y1;
            ae_valignx2 alIn;

            alIn = AE_LA128_PP(pIn);

            AE_LAV16X4X2_XP(p0, p1, alIn, castxcc(ae_int16x8, pIn), last0 * sizeof(ae_int16));
            AE_LAV16X4X2_XP(p2, p3, alIn, castxcc(ae_int16x8, pIn), last1 * sizeof(ae_int16));

            p0 = AE_MULFD16X16X4RAS(p0, _16384, _256, _m256);
            p1 = AE_MULFD16X16X4RAS(p1, _16384, _256, _m256);
            p2 = AE_MULFD16X16X4RAS(p2, _16384, _256, _m256);
            p3 = AE_MULFD16X16X4RAS(p3, _16384, _256, _m256);

            y0 = AE_SAT8X8X16(p0, p1);
            y1 = AE_SAT8X8X16(p2, p3);

            AE_SAV8X8X2_XP(y0, y1, alOut, castxcc(ae_int8x16, pOut), (w & 15) * sizeof(ae_int8));
        }
        AE_SA128POS_FP(alOut, pOut);
#else
        AE_SA128POS_FP(alOut, pOut);
        for (m = 0; m < (w&15) ; m++)
        {
            int16_t p = *pIn++;
            p -= 128*128;
            p >>= 7;
            p = XT_MIN(128, XT_MAX(-128, p));
            *pOut++ = (int8_t)p;
        }
#endif
    }
}

static size_t getScratchSize(const imgsize_t* in,const imgsize_t* out)
{
    (void)in; (void)out;
    imgsize_validate(in,  2, 1);
    imgsize_validate(out, 1, 0);
    return 0;
}
const imgcopy_api imgcopy_api_temp_gs8={getScratchSize,copy};
