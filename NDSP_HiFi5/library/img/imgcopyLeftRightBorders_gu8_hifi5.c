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
#include "imgpad_common.h"


/*------------------------------------------------
    Functon copies left/right edges
    Input:
    params - parameters of padding
    inImg  - input image
    Output:
    left[h] - left edge extended with lu and 
              lb points
    right[h]- right edge extended with ru and 
              rb points
------------------------------------------------*/
void imgcopyLeftRightBorders_gu8(void *left, void *right, const void * restrict inImg, const imgpad_params_t* params)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
    const ae_int8 * restrict in=(const ae_int8*)inImg;
    const uint8_t * restrict in_=(const uint8_t*)inImg;
    int n,h=params->out.height,y,y0=params->y,hin=params->in.height,istride=params->in.stride,win=params->in.width;
    uint8_t lu=in_[0];
    uint8_t ru=in_[win-1];
    uint8_t lb=in_[istride*(hin-1)];
    uint8_t rb=in_[istride*(hin-1)+win-1];
    uint8_t* restrict outl = (uint8_t*)left;
    uint8_t* restrict outr = (uint8_t*)right;
    ae_int8x8 xl, xr;
    ae_int8x16 * restrict pOutl = (ae_int8x16 *)outl;
    ae_int8x16 * restrict pOutr = (ae_int8x16 *)outr;
    ae_valign all, alr;
    ae_valignx2 allw, alrw;
    NASSERT(left );
    NASSERT(right);
    NASSERT_ALIGN(left ,ALIGNMENT);
    NASSERT_ALIGN(right,ALIGNMENT);
    y=XT_MAX(0,XT_MIN(y0,h));
    n = 0;
    if (params->fill==IMGPAD_EDGE)
    {
        if (y>0)
        {
            xl = AE_L8_I(in, 0);
            xr = AE_L8_X(in, win-1);
            for (n = 0; n < (y >> 4); n++)
            {
                AE_S8X8X2_IP(xl, xl, pOutl, 2*sizeof(ae_int8x8));
                AE_S8X8X2_IP(xr, xr, pOutr, 2*sizeof(ae_int8x8));
            }
            if (y & 8)
            {
                AE_S8X8_IP(xl, castxcc(ae_int8x8, pOutl), sizeof(ae_int8x8));
                AE_S8X8_IP(xr, castxcc(ae_int8x8, pOutr), sizeof(ae_int8x8));
            }
            for (n = (y&(~7)); n < y; n++)
            {
                outl[n] = lu;
                outr[n] = ru;
            }
        }
        int l, r;
        l = XT_MAX(0, y);
        r = XT_MIN(h, (y0 + hin));
        for (n=l; n<r; n++)
        {
            outl[n] = in_[istride*(n - y0)];
            outr[n] = in_[istride*(n - y0) + win - 1];
        }
        l = XT_MAX(l, r);
        if (l < h)
        {
            all = alr = AE_ZALIGN64();
            allw = alrw = AE_ZALIGN128();
            pOutl = (ae_int8x16 *)(outl + l);
            pOutr = (ae_int8x16 *)(outr + l);
            xl = AE_L8_X(in, istride*(hin-1));
            xr = AE_L8_X(in, istride*(hin-1)+win-1);           
            for (n = 0; n < ((h - l) >> 4); n++)
            {
                AE_SA8X8X2_IP(xl, xl, allw, pOutl);
                AE_SA8X8X2_IP(xr, xr, alrw, pOutr);
            }
            AE_SA128POS_FP(allw, pOutl);
            AE_SA128POS_FP(alrw, pOutr);
            if ((h-l) & 8)
            {
                AE_SA8X8_IP(xl, all, castxcc(ae_int8x8, pOutl));
                AE_SA8X8_IP(xr, alr, castxcc(ae_int8x8, pOutr));
                AE_SA64POS_FP(all, pOutl);
                AE_SA64POS_FP(alr, pOutr);
            }
            
            for (n = l + ((h - l)&(~7)); n < h; n++)
            {
                outl[n] = lb;
                outr[n] = rb;
            }
        }
    }
    else
    {
        uint8_t fill=(uint8_t)XT_MIN(255,XT_MAX(0,params->fill));
        const ae_int8 * tmp = (ae_int8*)(&fill);
        xl = AE_L8_I(tmp, 0);
        xr = xl;
        for (n = 0; n < (h >> 4); n++)
        {
            AE_S8X8X2_IP(xl, xl, pOutl, 2*sizeof(ae_int8x8));
            AE_S8X8X2_IP(xr, xr, pOutr, 2*sizeof(ae_int8x8));
        }
        if (h & 8)
        {
            AE_S8X8_IP(xl, castxcc(ae_int8x8, pOutl), sizeof(ae_int8x8));
            AE_S8X8_IP(xr, castxcc(ae_int8x8, pOutr), sizeof(ae_int8x8));
        }
        for (n = (h&(~7)); n < h; n++)
        {
            outl[n] = fill;
            outr[n] = fill;
        }
    }
#else
    const ae_int8 * restrict in=(const ae_int8*)inImg;
    const uint8_t * restrict in_=(const uint8_t*)inImg;
    int n,h=params->out.height,y,y0=params->y,hin=params->in.height,istride=params->in.stride,win=params->in.width;
    uint8_t* restrict outl = (uint8_t*)left;
    uint8_t* restrict outr = (uint8_t*)right;
    ae_int8x8 xl, xr;
    ae_int8x16 * restrict pOutl = (ae_int8x16 *)outl;
    ae_int8x16 * restrict pOutr = (ae_int8x16 *)outr;
    ae_valignx2 allw, alrw;
    NASSERT(left );
    NASSERT(right);
    NASSERT_ALIGN(left ,ALIGNMENT);
    NASSERT_ALIGN(right,ALIGNMENT);
    y=XT_MAX(0,XT_MIN(y0,h));
    n = 0;
    allw = alrw = AE_ZALIGN128();
    if (params->fill==IMGPAD_EDGE)
    {
        if (y>0)
        {
            xl = AE_L8_I(in, 0);
            xr = AE_L8_X(in, win-1);
            for (n = 0; n < (y >> 4); n++)
            {
                AE_S8X8X2_IP(xl, xl, pOutl, 2*sizeof(ae_int8x8));
                AE_S8X8X2_IP(xr, xr, pOutr, 2*sizeof(ae_int8x8));
            }
            if (y & 15)
            {
                int off = y&15;
                AE_SAV8X8X2_XP(xl, xl, allw, pOutl, off);
                AE_SAV8X8X2_XP(xr, xr, alrw, pOutr, off);
                AE_SA128POS_FP(allw, pOutl);
                AE_SA128POS_FP(alrw, pOutr);
            }         
        }
        int l, r;
        l = XT_MAX(0, y);
        r = XT_MIN(h, (y0 + hin));
        for (n=l; n<r; n++)
        {
            outl[n] = in_[istride*(n - y0)];
            outr[n] = in_[istride*(n - y0) + win - 1];
        }
        l = XT_MAX(l, r);
        if (l < h)
        {
            allw = alrw = AE_ZALIGN128();
            pOutl = (ae_int8x16 *)(outl + l);
            pOutr = (ae_int8x16 *)(outr + l);
            xl = AE_L8_X(in, istride*(hin-1));
            xr = AE_L8_X(in, istride*(hin-1)+win-1);           
            for (n = 0; n < ((h - l) >> 4); n++)
            {
                AE_SA8X8X2_IP(xl, xl, allw, pOutl);
                AE_SA8X8X2_IP(xr, xr, alrw, pOutr);
            }
            if ((h-l) & 15)
            {
                int off = (h-l)&15;
                AE_SAV8X8X2_XP(xl, xl, allw, pOutl, off);
                AE_SAV8X8X2_XP(xr, xr, alrw, pOutr, off);
            }
            AE_SA128POS_FP(allw, pOutl);
            AE_SA128POS_FP(alrw, pOutr);
        }
    }
    else
    {
        uint8_t fill=(uint8_t)XT_MIN(255,XT_MAX(0,params->fill));
        const ae_int8 * tmp = (ae_int8*)(&fill);
        xl = AE_L8_I(tmp, 0);
        xr = xl;
        for (n = 0; n < (h >> 4); n++)
        {
            AE_S8X8X2_IP(xl, xl, pOutl, 2*sizeof(ae_int8x8));
            AE_S8X8X2_IP(xr, xr, pOutr, 2*sizeof(ae_int8x8));
        }
        if (h & 15)
        {
            int off = h&15;
            allw = alrw = AE_ZALIGN128();
            AE_SAV8X8X2_XP(xl, xl, allw, pOutl, off);
            AE_SAV8X8X2_XP(xr, xr, alrw, pOutr, off);
            AE_SA128POS_FP(allw, pOutl);
            AE_SA128POS_FP(alrw, pOutr);
        }        
    }
#endif
}
