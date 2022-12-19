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
    Functon copies upper/bottom edges
    Input:
    params - parameters of padding
    inImg  - input image
    Output:
    up[h] -  upper edge extended with lu and 
             ru points
    btm[h]-  bottom edge extended with lb and 
             rb points
------------------------------------------------*/
void imgcopyUpBtmBorders_gu8(void *up, void *btm, const void * restrict inImg, const imgpad_params_t* params)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
    const uint8_t* restrict in=(const uint8_t*)inImg;
    const ae_int8* restrict in_=(const ae_int8*)inImg;
    int n,w=params->out.width,x,x0=params->x,win=params->in.width,istride=params->in.stride,hin=params->in.height;
    uint8_t lu=in[0];
    uint8_t ru=in[win-1];
    uint8_t lb=in[istride*(hin-1)];
    uint8_t rb=in[istride*(hin-1)+win-1];
    uint8_t* restrict outu = (uint8_t*)up;
    uint8_t* restrict outb = (uint8_t*)btm;
    const ae_int8x16 * restrict pInu;
    const ae_int8x16 * restrict pInb;
    ae_int8x16 * restrict pOutu = (ae_int8x16 *)outu;
    ae_int8x16 * restrict pOutb = (ae_int8x16 *)outb;
    ae_valign aliu, alib, alou, alob;
    ae_valignx2 aliuw, alibw, alouw, alobw;
    ae_int8x8 xu, xb;
    ae_int8x8 xu_, xb_;
    NASSERT(up );
    NASSERT(btm);
    NASSERT_ALIGN(up,ALIGNMENT);
    NASSERT_ALIGN(btm,ALIGNMENT);
    x=XT_MAX(0,XT_MIN(x0,w));
    n = 0;
    if (params->fill==IMGPAD_EDGE)
    {
        if (x>0)
        {
            xu = AE_L8_I(in_, 0);
            xb = AE_L8_X(in_, istride*(hin-1));
            for (n = 0; n < (x >> 4); n++)
            {
                AE_S8X8X2_IP(xu, xu, pOutu, 2 * sizeof(ae_int8x8));
                AE_S8X8X2_IP(xb, xb, pOutb, 2 * sizeof(ae_int8x8));
            }
            if (x & 8)
            {
                AE_S8X8_IP(xu, castxcc(ae_int8x8, pOutu), sizeof(ae_int8x8));
                AE_S8X8_IP(xb, castxcc(ae_int8x8, pOutb), sizeof(ae_int8x8));
            }
            for (n = (x&(~7)); n < x; n++)
            {
                outu[n] = lu;
                outb[n] = lb;
            }
        }
        int l, r;
        l = XT_MAX(0, x);
        r = XT_MIN(w, (x0 + win));

        if (r > l)
        {         
            pOutu = (ae_int8x16 *)(outu + l);
            pOutb = (ae_int8x16 *)(outb + l);
            pInu = (const ae_int8x16 *)(in + (l - x0));
            pInb = (const ae_int8x16 *)(in + (l - x0) + (hin - 1)*istride);

            alouw = AE_ZALIGN128();
            alobw = AE_ZALIGN128();
            aliuw = AE_LA128_PP(pInu);            
            alibw = AE_LA128_PP(pInb);
            for (n = 0; n < ((r - l) >> 4); n++)
            {
                AE_LA8X8X2_IP(xu, xu_, aliuw, pInu);
                AE_LA8X8X2_IP(xb, xb_, alibw, pInb);
                AE_SA8X8X2_IP(xu, xu_, alouw, pOutu);               
                AE_SA8X8X2_IP(xb, xb_, alobw, pOutb);
            }
            AE_SA128POS_FP(alouw, pOutu);
            AE_SA128POS_FP(alobw, pOutb);
            if((r - l)& 8)
            {
                alou = AE_ZALIGN64();
                aliu = AE_LA64_PP(pInu);
                alob = AE_ZALIGN64();
                alib = AE_LA64_PP(pInb);
                AE_LA8X8_IP(xu, aliu, castxcc(ae_int8x8, pInu ));                
                AE_LA8X8_IP(xb, alib, castxcc(ae_int8x8, pInb ));
                AE_SA8X8_IP(xu, alou, castxcc(ae_int8x8, pOutu));
                AE_SA8X8_IP(xb, alob, castxcc(ae_int8x8, pOutb));
                AE_SA64POS_FP(alou, pOutu);
                AE_SA64POS_FP(alob, pOutb);
            }
            for (n = l + ((r - l)&(~7)); n < r; n++)
            {
                outu[n] = in[(n - x0)];
                outb[n] = in[(n - x0) + (hin - 1)*istride];
            }
        }

        l = XT_MAX(l, r);
        if (l < w)
        {
            alobw = alouw = AE_ZALIGN128();
            pOutu = (ae_int8x16 *)(outu + l);
            pOutb = (ae_int8x16 *)(outb + l);
            xu = AE_L8_X(in_, win-1);
            xb = AE_L8_X(in_, istride*(hin-1) + win-1);
            for (n = 0; n < ((w - l) >> 4); n++)
            {
                AE_SA8X8X2_IP(xu, xu, alouw, pOutu);
                AE_SA8X8X2_IP(xb, xb, alobw, pOutb);
            }
            AE_SA128POS_FP(alouw, pOutu);
            AE_SA128POS_FP(alobw, pOutb);
            if ((w-l) & 8)
            {
                alob = alou = AE_ZALIGN64();
                AE_SA8X8_IP(xu, alou, castxcc(ae_int8x8, pOutu));
                AE_SA8X8_IP(xb, alob, castxcc(ae_int8x8, pOutb));
                AE_SA64POS_FP(alou, pOutu);
                AE_SA64POS_FP(alob, pOutb);
            }            
            for (n = l + ((w - l)&(~7)); n < w; n++)
            {
                outu[n] = ru;
                outb[n] = rb;
            }
        }
    }
    else
    {
        uint8_t fill=(uint8_t)XT_MIN(255,XT_MAX(0,params->fill));
        const ae_int8 * tmp = (ae_int8*)(&fill);
        xu = AE_L8_I(tmp, 0);
        xb = xu;
        for (n = 0; n < (w >> 4); n++)
        {
            AE_S8X8X2_IP(xu, xu, pOutu, 2*sizeof(ae_int8x8));
            AE_S8X8X2_IP(xb, xb, pOutb, 2*sizeof(ae_int8x8));
        }
        if (w & 8)
        {
            AE_S8X8_IP(xu, castxcc(ae_int8x8, pOutu), sizeof(ae_int8x8));
            AE_S8X8_IP(xb, castxcc(ae_int8x8, pOutb), sizeof(ae_int8x8));
        }
        for (n = (w&(~7)); n < w; n++)
        {
            outu[n] = fill;
            outb[n] = fill;
        }
    }
#else
    const uint8_t* restrict in=(const uint8_t*)inImg;
    const ae_int8* restrict in_=(const ae_int8*)inImg;
    int n,w=params->out.width,x,x0=params->x,win=params->in.width,istride=params->in.stride,hin=params->in.height;
    //uint8_t lu=in[0];
    //uint8_t ru=in[win-1];
    //uint8_t lb=in[istride*(hin-1)];
    //uint8_t rb=in[istride*(hin-1)+win-1];
    uint8_t* restrict outu = (uint8_t*)up;
    uint8_t* restrict outb = (uint8_t*)btm;
    const ae_int8x16 * restrict pInu;
    const ae_int8x16 * restrict pInb;
    ae_int8x16 * restrict pOutu = (ae_int8x16 *)outu;
    ae_int8x16 * restrict pOutb = (ae_int8x16 *)outb;
    //ae_valign aliu, alib, alou, alob;
    ae_valignx2 aliuw, alibw, alouw, alobw;
    ae_int8x8 xu, xb;
    ae_int8x8 xu_, xb_;
    NASSERT(up );
    NASSERT(btm);
    NASSERT_ALIGN(up,ALIGNMENT);
    NASSERT_ALIGN(btm,ALIGNMENT);
    x=XT_MAX(0,XT_MIN(x0,w));
    n = 0;
    alouw = AE_ZALIGN128();
    alobw = AE_ZALIGN128();
    if (params->fill==IMGPAD_EDGE)
    {
        if (x>0)
        {
            xu = AE_L8_I(in_, 0);
            xb = AE_L8_X(in_, istride*(hin-1));
            for (n = 0; n < (x >> 4); n++)
            {
                AE_S8X8X2_IP(xu, xu, pOutu, 2 * sizeof(ae_int8x8));
                AE_S8X8X2_IP(xb, xb, pOutb, 2 * sizeof(ae_int8x8));
            }
            if (x & 15)
            {
                int off = x & 15;
                AE_SAV8X8X2_XP(xu, xu, alouw, pOutu, off);
                AE_SAV8X8X2_XP(xb, xb, alobw, pOutb, off);
                AE_SA128POS_FP(alouw, pOutu);
                AE_SA128POS_FP(alobw, pOutb);
            }            
            
        }
        int l, r;
        l = XT_MAX(0, x);
        r = XT_MIN(w, (x0 + win));

        if (r > l)
        {         
            pOutu = (ae_int8x16 *)(outu + l);
            pOutb = (ae_int8x16 *)(outb + l);
            pInu = (const ae_int8x16 *)(in + (l - x0));
            pInb = (const ae_int8x16 *)(in + (l - x0) + (hin - 1)*istride);

            alouw = AE_ZALIGN128();
            alobw = AE_ZALIGN128();
            aliuw = AE_LA128_PP(pInu);            
            alibw = AE_LA128_PP(pInb);
            for (n = 0; n < ((r - l) >> 4); n++)
            {
                AE_LA8X8X2_IP(xu, xu_, aliuw, pInu);
                AE_LA8X8X2_IP(xb, xb_, alibw, pInb);
                AE_SA8X8X2_IP(xu, xu_, alouw, pOutu);               
                AE_SA8X8X2_IP(xb, xb_, alobw, pOutb);
            }
            if((r - l)& 15)
            {
                int off = (r-l) & 15;
                AE_LAV8X8X2_XP(xu, xu_, aliuw,  pInu, off);
                AE_LAV8X8X2_XP(xb, xb_, alibw,  pInb, off);
                AE_SAV8X8X2_XP(xu, xu_, alouw, pOutu, off);               
                AE_SAV8X8X2_XP(xb, xb_, alobw, pOutb, off);
            }
            AE_SA128POS_FP(alouw, pOutu);
            AE_SA128POS_FP(alobw, pOutb);
        }

        l = XT_MAX(l, r);
        if (l < w)
        {
            alobw = alouw = AE_ZALIGN128();
            pOutu = (ae_int8x16 *)(outu + l);
            pOutb = (ae_int8x16 *)(outb + l);
            xu = AE_L8_X(in_, win-1);
            xb = AE_L8_X(in_, istride*(hin-1) + win-1);
            for (n = 0; n < ((w - l) >> 4); n++)
            {
                AE_SA8X8X2_IP(xu, xu, alouw, pOutu);
                AE_SA8X8X2_IP(xb, xb, alobw, pOutb);
            }  
            if ((w-l) & 15)
            {
                int off = (w-l) & 15;
                AE_SAV8X8X2_XP(xu, xu, alouw, pOutu, off);
                AE_SAV8X8X2_XP(xb, xb, alobw, pOutb, off);
            }            
            AE_SA128POS_FP(alouw, pOutu);
            AE_SA128POS_FP(alobw, pOutb);
        }
    }
    else
    {
        uint8_t fill=(uint8_t)XT_MIN(255,XT_MAX(0,params->fill));
        const ae_int8 * tmp = (ae_int8*)(&fill);
        xu = AE_L8_I(tmp, 0);
        xb = xu;
        for (n = 0; n < (w >> 4); n++)
        {
            AE_S8X8X2_IP(xu, xu, pOutu, 2*sizeof(ae_int8x8));
            AE_S8X8X2_IP(xb, xb, pOutb, 2*sizeof(ae_int8x8));
        }
        if (w & 15)
        {
            int off = w & 15;
            alobw = alouw = AE_ZALIGN128();
            AE_SAV8X8X2_XP(xu, xu, alouw, pOutu, off);
            AE_SAV8X8X2_XP(xb, xb, alobw, pOutb, off);
            AE_SA128POS_FP(alouw, pOutu);
            AE_SA128POS_FP(alobw, pOutb);
        }            
        
    }
#endif
}
