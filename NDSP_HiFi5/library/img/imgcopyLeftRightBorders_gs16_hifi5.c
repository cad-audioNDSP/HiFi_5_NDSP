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
void imgcopyLeftRightBorders_gs16(void *left, void *right, const void * restrict inImg, const imgpad_params_t* params)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
    const int16_t* restrict in=(const int16_t*)inImg;
    int16_t* restrict outl = (int16_t*)left;
    int16_t* restrict outr = (int16_t*)right;
    ae_int16x8 * restrict pOutl = (ae_int16x8 *)outl;
    ae_int16x8 * restrict pOutr = (ae_int16x8 *)outr;
    int n,h=params->out.height,y,y0=params->y,hin=params->in.height,istride=params->in.stride,win=params->in.width;
    int16_t lu=in[0];
    int16_t ru=in[win-1];
    int16_t lb=in[istride*(hin-1)];
    int16_t rb=in[istride*(hin-1)+win-1];
    ae_int16x4 xl, xr;
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
            xl = AE_MOVDA16(lu);
            xr = AE_MOVDA16(ru);
            for (n = 0; n < (y >> 3); n++)
            {
                AE_S16X4X2_IP(xl, xl,  pOutl, 8 * sizeof(int16_t));
                AE_S16X4X2_IP(xr, xr,  pOutr, 8 * sizeof(int16_t));
            }
            if (y & 4)
            {
                AE_S16X4_IP(xl, castxcc(ae_int16x4, pOutl), 4 * sizeof(int16_t));
                AE_S16X4_IP(xr, castxcc(ae_int16x4, pOutr), 4 * sizeof(int16_t));
            }
            for (n = (y&(~3)); n < y; n++)
            {
                outl[n] = lu;
                outr[n] = ru;
            }
        }
        int l, r;
        l = XT_MAX(0, y);
        r = XT_MIN(h, (y0 + hin));
        for(n=l; n<r; n++)
        {
            outl[n]=in[istride*(n-y0)];
            outr[n]=in[istride*(n-y0)+win-1]; 
        }

        l = XT_MAX(l, r);
        if (l < h)
        {
            allw = alrw = AE_ZALIGN128();
            pOutl = (ae_int16x8 *)(outl + l);
            pOutr = (ae_int16x8 *)(outr + l);
            xl = AE_MOVDA16(lb);
            xr = AE_MOVDA16(rb);
            for (n = 0; n < ((h - l) >> 3); n++)
            {
                AE_SA16X4X2_IP(xl, xl, allw, pOutl);
                AE_SA16X4X2_IP(xr, xr, alrw, pOutr);
            }
            AE_SA128POS_FP(allw, pOutl);
            AE_SA128POS_FP(alrw, pOutr);
            if ((h-l)&4)
            {
                all = alr = AE_ZALIGN64();
                AE_SA16X4_IP(xl, all, castxcc(ae_int16x4, pOutl));
                AE_SA16X4_IP(xr, alr, castxcc(ae_int16x4, pOutr));
                AE_SA64POS_FP(all, pOutl);
                AE_SA64POS_FP(alr, pOutr);
            }
            
            for (n = l + ((h - l)&(~3)); n < h; n++)
            {
                outl[n] = lb;
                outr[n] = rb;
            }
        }
    }
    else
    {
        int16_t fill=(int16_t)XT_MIN(32767,XT_MAX(0,params->fill));
        xl = AE_MOVDA16(fill);
        xr = AE_MOVDA16(fill);
        for (n = 0; n < (h >> 3); n++)
        {
            AE_S16X4X2_IP(xl, xl, pOutl, 8 * sizeof(int16_t));
            AE_S16X4X2_IP(xr, xr, pOutr, 8 * sizeof(int16_t));
        }
        if (h & 4)
        {
            AE_S16X4_IP(xl, castxcc(ae_int16x4, pOutl), 4 * sizeof(int16_t));
            AE_S16X4_IP(xr, castxcc(ae_int16x4, pOutr), 4 * sizeof(int16_t));
        }
        for (n = (h&(~3)); n < h; n++)
        {
            outl[n] = fill;
            outr[n] = fill;
        }
    }
#else
    const int16_t* restrict in=(const int16_t*)inImg;
    int16_t* restrict outl = (int16_t*)left;
    int16_t* restrict outr = (int16_t*)right;
    ae_int16x8 * restrict pOutl = (ae_int16x8 *)outl;
    ae_int16x8 * restrict pOutr = (ae_int16x8 *)outr;
    int n,h=params->out.height,y,y0=params->y,hin=params->in.height,istride=params->in.stride,win=params->in.width;
    int16_t lu=in[0];
    int16_t ru=in[win-1];
    int16_t lb=in[istride*(hin-1)];
    int16_t rb=in[istride*(hin-1)+win-1];
    ae_int16x4 xl, xr;
    //ae_valign all, alr;
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
            xl = AE_MOVDA16(lu);
            xr = AE_MOVDA16(ru);
            for (n = 0; n < (y >> 3); n++)
            {
                AE_S16X4X2_IP(xl, xl,  pOutl, 8 * sizeof(int16_t));
                AE_S16X4X2_IP(xr, xr,  pOutr, 8 * sizeof(int16_t));
            }
            if (y & 7)
            {
                int off = (y&7)<<1;
                AE_SAV16X4X2_XP(xl, xl, allw, pOutl, off);
                AE_SAV16X4X2_XP(xr, xr, alrw, pOutr, off);
                AE_SA128POS_FP(allw, pOutl);
                AE_SA128POS_FP(alrw, pOutr);
            }  
        }
        int l, r;
        l = XT_MAX(0, y);
        r = XT_MIN(h, (y0 + hin));
        for(n=l; n<r; n++)
        {
            outl[n]=in[istride*(n-y0)];
            outr[n]=in[istride*(n-y0)+win-1]; 
        }

        l = XT_MAX(l, r);
        if (l < h)
        {
            allw = alrw = AE_ZALIGN128();
            pOutl = (ae_int16x8 *)(outl + l);
            pOutr = (ae_int16x8 *)(outr + l);
            xl = AE_MOVDA16(lb);
            xr = AE_MOVDA16(rb);
            for (n = 0; n < ((h - l) >> 3); n++)
            {
                AE_SA16X4X2_IP(xl, xl, allw, pOutl);
                AE_SA16X4X2_IP(xr, xr, alrw, pOutr);
            }
            if ((h-l)&7)
            {
                int off = ((h-l)&7)<<1;
                AE_SAV16X4X2_XP(xl, xl, allw, pOutl, off);
                AE_SAV16X4X2_XP(xr, xr, alrw, pOutr, off);
            }
            AE_SA128POS_FP(allw, pOutl);
            AE_SA128POS_FP(alrw, pOutr);
        }
    }
    else
    {
        int16_t fill=(int16_t)XT_MIN(32767,XT_MAX(0,params->fill));
        xl = AE_MOVDA16(fill);
        xr = AE_MOVDA16(fill);
        for (n = 0; n < (h >> 3); n++)
        {
            AE_S16X4X2_IP(xl, xl, pOutl, 8 * sizeof(int16_t));
            AE_S16X4X2_IP(xr, xr, pOutr, 8 * sizeof(int16_t));
        }
        if (h & 7)
        {
            int off = (h&7)<<1;
            allw = alrw = AE_ZALIGN128();
            AE_SAV16X4X2_XP(xl, xl, allw, pOutl, off);
            AE_SAV16X4X2_XP(xr, xr, alrw, pOutr, off);
            AE_SA128POS_FP(allw, pOutl);
            AE_SA128POS_FP(alrw, pOutr);
        }
        
    }
#endif
}
