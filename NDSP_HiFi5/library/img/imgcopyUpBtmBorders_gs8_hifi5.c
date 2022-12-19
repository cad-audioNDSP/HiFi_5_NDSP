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
void imgcopyUpBtmBorders_gs8(void *up, void *btm, const void * restrict inImg, const imgpad_params_t* params)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
    int n,w=params->out.width;
    int8_t* restrict outu = (int8_t*)up;
    int8_t* restrict outb = (int8_t*)btm;
    ae_int8x16 * restrict pOutu = (ae_int8x16 *)outu;
    ae_int8x16 * restrict pOutb = (ae_int8x16 *)outb;
    ae_int8x8 xu, xb;
    NASSERT(up );
    NASSERT(btm);
    NASSERT_ALIGN(up,ALIGNMENT);
    NASSERT_ALIGN(btm,ALIGNMENT);
    if (params->fill==IMGPAD_EDGE)
    {
        imgcopyUpBtmBorders_gu8(up, btm, inImg, params);
    }
    else
    {
        int8_t fill=(int8_t)XT_MIN(127,XT_MAX(-128,params->fill));
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
    int n,w=params->out.width;
    ae_int8x16 * restrict pOutu = (ae_int8x16 *)up;
    ae_int8x16 * restrict pOutb = (ae_int8x16 *)btm;
    ae_valignx2 alouw, alobw;
    ae_int8x8 xu, xb;
    NASSERT(up );
    NASSERT(btm);
    NASSERT_ALIGN(up,ALIGNMENT);
    NASSERT_ALIGN(btm,ALIGNMENT);
    alouw = AE_ZALIGN128();
    alobw = AE_ZALIGN128();
    if (params->fill==IMGPAD_EDGE)
    {
        imgcopyUpBtmBorders_gu8(up, btm, inImg, params);
    }
    else
    {
        int8_t fill=(int8_t)XT_MIN(127,XT_MAX(-128,params->fill));
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
