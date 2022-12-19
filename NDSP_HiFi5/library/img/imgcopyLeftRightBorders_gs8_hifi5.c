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
void imgcopyLeftRightBorders_gs8(void *left, void *right, const void * restrict inImg, const imgpad_params_t* params)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
    int n,h=params->out.height;
    int8_t* restrict outl = (int8_t*)left;
    int8_t* restrict outr = (int8_t*)right;
    ae_int8x8 xl, xr;
    ae_int8x16 * restrict pOutl = (ae_int8x16 *)outl;
    ae_int8x16 * restrict pOutr = (ae_int8x16 *)outr;
    NASSERT(left );
    NASSERT(right);
    NASSERT_ALIGN(left ,ALIGNMENT);
    NASSERT_ALIGN(right,ALIGNMENT);
    if (params->fill==IMGPAD_EDGE)
    {
        imgcopyLeftRightBorders_gu8(left, right,  inImg, params);
    }
    else
    {
        int8_t fill=(int8_t)XT_MIN(127,XT_MAX(-128,params->fill));
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
    int n,h=params->out.height;
    ae_int8x8 xl, xr;
    ae_int8x16 * restrict pOutl = (ae_int8x16 *)left;
    ae_int8x16 * restrict pOutr = (ae_int8x16 *)right;
    ae_valignx2 allw, alrw;
    NASSERT(left );
    NASSERT(right);
    NASSERT_ALIGN(left ,ALIGNMENT);
    NASSERT_ALIGN(right,ALIGNMENT);
    allw = alrw = AE_ZALIGN128();
    if (params->fill==IMGPAD_EDGE)
    {
        imgcopyLeftRightBorders_gu8(left, right,  inImg, params);
    }
    else
    {
        int8_t fill=(int8_t)XT_MIN(127,XT_MAX(-128,params->fill));
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
