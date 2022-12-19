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
void imgcopyUpBtmBorders_gs16(void *up, void *btm, const void * restrict inImg, const imgpad_params_t* params)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
    const int16_t* restrict in=(const int16_t*)inImg;
    int n,w=params->out.width,x,x0=params->x,win=params->in.width,istride=params->in.stride,hin=params->in.height;
    int16_t lu=in[0];
    int16_t ru=in[win-1];
    int16_t lb=in[istride*(hin-1)];
    int16_t rb=in[istride*(hin-1)+win-1];
    int16_t* restrict outu = (int16_t*)up;
    int16_t* restrict outb = (int16_t*)btm;
    const ae_int16x8 * restrict pInu;
    const ae_int16x8 * restrict pInb;
    ae_int16x8 * restrict pOutu = (ae_int16x8 *)outu;
    ae_int16x8 * restrict pOutb = (ae_int16x8 *)outb;
    ae_valign aliu, alib, alou, alob;
    ae_valignx2 aliuw, alibw, alouw, alobw;
    ae_int16x4 xu, xb, xu_, xb_;
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
            xu = AE_MOVDA16(lu);
            xb = AE_MOVDA16(lb);
            for (n = 0; n < (x >> 3); n++)
            {
                AE_S16X4X2_IP(xu, xu, pOutu, 8 * sizeof(int16_t));
                AE_S16X4X2_IP(xb, xb, pOutb, 8 * sizeof(int16_t));
            }
            if (x & 4)
            {
                AE_S16X4_IP(xu, castxcc(ae_int16x4, pOutu), 4 * sizeof(int16_t));
                AE_S16X4_IP(xb, castxcc(ae_int16x4, pOutb), 4 * sizeof(int16_t));
            }
            for (n = (x&(~3)); n < x; n++)
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
            alobw = alouw = AE_ZALIGN128();
            pOutu = (ae_int16x8 *)(outu + l);
            pOutb = (ae_int16x8 *)(outb + l);
            pInu = (const ae_int16x8 *)(in + (l - x0));
            pInb = (const ae_int16x8 *)(in + (l - x0) + (hin - 1)*istride);
            aliuw = AE_LA128_PP(pInu);
            alibw = AE_LA128_PP(pInb);
            for (n = 0; n<((r-l)>>3); n++)
            {
                AE_LA16X4X2_IP(xu, xu_, aliuw, pInu);
                AE_LA16X4X2_IP(xb, xb_, alibw, pInb);
                AE_SA16X4X2_IP(xu, xu_, alouw, pOutu);
                AE_SA16X4X2_IP(xb, xb_, alobw, pOutb);
            }
            AE_SA128POS_FP(alobw, pOutb);
            AE_SA128POS_FP(alouw, pOutu);
            if ((r-l)&4)
            {
                aliu = AE_LA64_PP(pInu);
                alib = AE_LA64_PP(pInb);
                alob = alou = AE_ZALIGN64();
                AE_LA16X4_IP(xu, aliu, castxcc(ae_int16x4, pInu ));
                AE_LA16X4_IP(xb, alib, castxcc(ae_int16x4, pInb ));
                AE_SA16X4_IP(xu, alou, castxcc(ae_int16x4, pOutu));
                AE_SA16X4_IP(xb, alob, castxcc(ae_int16x4, pOutb));
                AE_SA64POS_FP(alob, pOutb);
                AE_SA64POS_FP(alou, pOutu);
            }
            
            for (n = l + ((r - l)&(~3)); n<r; n++)
            {
                outu[n] = in[(n - x0)];
                outb[n] = in[(n - x0) + (hin - 1)*istride];
            }
        }
        l = XT_MAX(l, r);
        if (l < w)
        {
            alobw = alouw = AE_ZALIGN128();
            pOutu = (ae_int16x8 *)(outu + l);
            pOutb = (ae_int16x8 *)(outb + l);
            xu = AE_MOVDA16(ru);
            xb = AE_MOVDA16(rb);
            for (n = 0; n < ((w-l)>>3); n++)
            {
                AE_SA16X4X2_IP(xu, xu, alouw, pOutu);
                AE_SA16X4X2_IP(xb, xb, alobw, pOutb);
            }
            AE_SA128POS_FP(alobw, pOutb);
            AE_SA128POS_FP(alouw, pOutu);
            if ((w-l)&4)
            {
                alob = alou = AE_ZALIGN64();
                AE_SA16X4_IP(xu, alou, castxcc(ae_int16x4, pOutu));
                AE_SA16X4_IP(xb, alob, castxcc(ae_int16x4, pOutb));
                AE_SA64POS_FP(alob, pOutb);
                AE_SA64POS_FP(alou, pOutu);
            }
            
            for (n = l+ ((w-l)&(~3)); n < w; n++)
            {
                outu[n] = ru;
                outb[n] = rb;
            }
        }

    }
    else
    {
        int16_t fill = (int16_t)XT_MIN(32767, XT_MAX(0, params->fill));
        xu = AE_MOVDA16(fill);
        xb = AE_MOVDA16(fill);
        for (n = 0; n < (w >> 3); n++)
        {
            AE_S16X4X2_IP(xu, xu,  pOutu, 8 * sizeof(int16_t));
            AE_S16X4X2_IP(xb, xb,  pOutb, 8 * sizeof(int16_t));
        }
        if (w & 4)
        {
            AE_S16X4_IP(xu, castxcc(ae_int16x4, pOutu), 4 * sizeof(int16_t));
            AE_S16X4_IP(xb, castxcc(ae_int16x4, pOutb), 4 * sizeof(int16_t));
        }
        for (n = (w&(~3)); n < w; n++)
        {
            outu[n] = fill;
            outb[n] = fill;
        }
    }
#else
    const int16_t* restrict in=(const int16_t*)inImg;
    int n,w=params->out.width,x,x0=params->x,win=params->in.width,istride=params->in.stride,hin=params->in.height;
    int16_t lu=in[0];
    int16_t ru=in[win-1];
    int16_t lb=in[istride*(hin-1)];
    int16_t rb=in[istride*(hin-1)+win-1];
    int16_t* restrict outu = (int16_t*)up;
    int16_t* restrict outb = (int16_t*)btm;
    const ae_int16x8 * restrict pInu;
    const ae_int16x8 * restrict pInb;
    ae_int16x8 * restrict pOutu = (ae_int16x8 *)outu;
    ae_int16x8 * restrict pOutb = (ae_int16x8 *)outb;
    //ae_valign aliu, alib, alou, alob;
    ae_valignx2 aliuw, alibw, alouw, alobw;
    ae_int16x4 xu, xb, xu_, xb_;
    NASSERT(up );
    NASSERT(btm);
    NASSERT_ALIGN(up,ALIGNMENT);
    NASSERT_ALIGN(btm,ALIGNMENT);
    x=XT_MAX(0,XT_MIN(x0,w));
    n = 0;
    alobw = alouw = AE_ZALIGN128();
    if (params->fill==IMGPAD_EDGE)
    {
        if (x>0)
        {
            xu = AE_MOVDA16(lu);
            xb = AE_MOVDA16(lb);
            for (n = 0; n < (x >> 3); n++)
            {
                AE_S16X4X2_IP(xu, xu, pOutu, 8 * sizeof(int16_t));
                AE_S16X4X2_IP(xb, xb, pOutb, 8 * sizeof(int16_t));
            }
            if (x & 7)
            {
                int off = (x&7)<<1;
                AE_SAV16X4X2_XP(xu, xu, alouw, pOutu, off);
                AE_SAV16X4X2_XP(xb, xb, alobw, pOutb, off);
                AE_SA128POS_FP(alobw, pOutb);
                AE_SA128POS_FP(alouw, pOutu);
            }        
        }
        int l, r;
        l = XT_MAX(0, x);
        r = XT_MIN(w, (x0 + win));
        if (r > l)
        {
            alobw = alouw = AE_ZALIGN128();
            pOutu = (ae_int16x8 *)(outu + l);
            pOutb = (ae_int16x8 *)(outb + l);
            pInu = (const ae_int16x8 *)(in + (l - x0));
            pInb = (const ae_int16x8 *)(in + (l - x0) + (hin - 1)*istride);
            aliuw = AE_LA128_PP(pInu);
            alibw = AE_LA128_PP(pInb);
            for (n = 0; n<((r-l)>>3); n++)
            {
                AE_LA16X4X2_IP(xu, xu_, aliuw, pInu);
                AE_LA16X4X2_IP(xb, xb_, alibw, pInb);
                AE_SA16X4X2_IP(xu, xu_, alouw, pOutu);
                AE_SA16X4X2_IP(xb, xb_, alobw, pOutb);
            }
            if ((r-l)&7)
            {
                int off = ((r-l)&7)<<1;
                AE_LAV16X4X2_XP(xu, xu_, aliuw, pInu , off);
                AE_LAV16X4X2_XP(xb, xb_, alibw, pInb , off);
                AE_SAV16X4X2_XP(xu, xu_, alouw, pOutu, off);
                AE_SAV16X4X2_XP(xb, xb_, alobw, pOutb, off);
            }
            AE_SA128POS_FP(alobw, pOutb);
            AE_SA128POS_FP(alouw, pOutu);
        }
        l = XT_MAX(l, r);
        if (l < w)
        {
            alobw = alouw = AE_ZALIGN128();
            pOutu = (ae_int16x8 *)(outu + l);
            pOutb = (ae_int16x8 *)(outb + l);
            xu = AE_MOVDA16(ru);
            xb = AE_MOVDA16(rb);
            for (n = 0; n < ((w-l)>>3); n++)
            {
                AE_SA16X4X2_IP(xu, xu, alouw, pOutu);
                AE_SA16X4X2_IP(xb, xb, alobw, pOutb);
            }
            if ((w-l)&7)
            {
                int off = ((w-l)&7)<<1;
                AE_SAV16X4X2_XP(xu, xu, alouw, pOutu, off);
                AE_SAV16X4X2_XP(xb, xb, alobw, pOutb, off);
            }
            AE_SA128POS_FP(alobw, pOutb);
            AE_SA128POS_FP(alouw, pOutu);
        }
    }
    else
    {
        int16_t fill = (int16_t)XT_MIN(32767, XT_MAX(0, params->fill));
        xu = AE_MOVDA16(fill);
        xb = AE_MOVDA16(fill);
        for (n = 0; n < (w >> 3); n++)
        {
            AE_S16X4X2_IP(xu, xu,  pOutu, 8 * sizeof(int16_t));
            AE_S16X4X2_IP(xb, xb,  pOutb, 8 * sizeof(int16_t));
        }
        if ((w)&7)
        {
            int off = ((w)&7)<<1;
            alobw = alouw = AE_ZALIGN128();
            AE_SAV16X4X2_XP(xu, xu, alouw, pOutu, off);
            AE_SAV16X4X2_XP(xb, xb, alobw, pOutb, off);
            AE_SA128POS_FP(alobw, pOutb);
            AE_SA128POS_FP(alouw, pOutu);
        }      
    }
#endif
}
