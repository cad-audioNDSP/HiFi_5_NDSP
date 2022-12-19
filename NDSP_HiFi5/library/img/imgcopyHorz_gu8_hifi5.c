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
    Functon copies horizontal edges to rectange
    Input:
    edge[h] - vertical edge 
    w,h,ostride - copied region
    Output:
    outImg   - output image
------------------------------------------------*/
void imgcopyHorz_gu8(void* outImg, const void* edge, int w, int h, int ostride)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
          uint8_t * restrict out=(uint8_t * )outImg;
    const uint8_t * restrict in =(const uint8_t * )edge;
    int m,n;
    ae_int8x16 * restrict pOut;
    const ae_int8x16 * restrict pIn = (const ae_int8x16 *)in;
    ae_int8x8 x0, x1;
    ae_valign al_out;
    ae_valignx2 al_outw;
    uint8_t tmp;
    NASSERT_ALIGN(edge,ALIGNMENT);
    if (!(ostride & 15) && !((uintptr_t)out & 15))
    {
        for (m = 0; m<(w >> 4); m++)
        {
            AE_L8X8X2_IP(x0, x1, pIn, 2* sizeof(ae_int8x8));
            pOut = (ae_int8x16 *)(out + 16 * m);
            for (n = 0; n<h; n++)
            {
                AE_S8X8X2_XP(x0, x1, pOut, sizeof(uint8_t)*ostride);
            }
        }
        if (w & 8)
        {
            AE_L8X8_IP(x0, castxcc(ae_int8x8, pIn), sizeof(ae_int8x8));
            pOut = (ae_int8x16 *)(out + (w&~15));
            for (n = 0; n<h; n++)
            {
                AE_S8X8_XP(x0, castxcc(ae_int8x8, pOut), sizeof(uint8_t)*ostride);
            }
        }
        for (m = (w&(~7)); m < w; m++)
        {
            tmp = in[m];
            for (n = 0; n < h; n++)
            {
                out[n*ostride + m] = tmp;
            }
        }
    }
    else
    {
        al_out = AE_ZALIGN64();
        al_outw = AE_ZALIGN128();
        for (n = 0; n < h; n++)
        {
            pOut = (ae_int8x16 *)(out + n*ostride);
            pIn = (const ae_int8x16 *)in;
            for (m = 0; m < (w >> 4); m++)
            {
                AE_L8X8X2_IP(x0, x1, pIn, 2*sizeof(ae_int8x8));
                AE_SA8X8X2_IP(x0, x1, al_outw, pOut);
            }
            AE_SA128POS_FP(al_outw, pOut);
            if (w & 8)
            {
                AE_L8X8_IP(x0, castxcc(ae_int8x8, pIn), sizeof(ae_int8x8));
                AE_SA8X8_IP(x0, al_out, castxcc(ae_int8x8, pOut));
                AE_SA64POS_FP(al_out, pOut);
            }
        }
        for (m = (w&(~7)); m < w; m++)
        {
            tmp = in[m];
            for (n = 0; n < h; n++)
            {
                out[n*ostride + m] = tmp;
            }
        }
    }
#else
    uint8_t * restrict out=(uint8_t * )outImg;
    const uint8_t * restrict in =(const uint8_t * )edge;
    int m,n;
    ae_int8x16 * restrict pOut;
    const ae_int8x16 * restrict pIn = (const ae_int8x16 *)in;
    ae_int8x8 x0, x1;
    ae_valignx2 al_outw, al_inw;
    NASSERT_ALIGN(edge,ALIGNMENT);
    if (!(ostride & 15) && !((uintptr_t)out & 15))
#if 0
    {
        al_outw = AE_ZALIGN128();
        for (n = 0; n < h; n++)
        {
            pOut = (ae_int8x16 *)(out + n*ostride);
            pIn = (const ae_int8x16 *)in;
            for (m = 0; m < (w >> 4); m++)
            {
                AE_L8X8X2_IP(x0, x1, pIn,  2*sizeof(ae_int8x8));
                AE_S8X8X2_IP(x0, x1, pOut, 2*sizeof(ae_int8x8));
            }            
            if (w & 15)
            {
                int off = w & 15;
                al_inw = AE_LA128_PP(pIn);
                AE_LAV8X8X2_XP(x0, x1, al_inw,  pIn,  off);
                AE_SAV8X8X2_XP(x0, x1, al_outw, pOut, off);
                AE_SA128POS_FP(al_outw, pOut);
            }            
        }
    }
#else
    {
        uint8_t tmp;
        for (m = 0; m<(w >> 4); m++)
        {
            AE_L8X8X2_IP(x0, x1, pIn, 2* sizeof(ae_int8x8));
            pOut = (ae_int8x16 *)(out + 16 * m);
            for (n = 0; n<h; n++)
            {
                AE_S8X8X2_XP(x0, x1, pOut, sizeof(uint8_t)*ostride);
            }
        }
        if (w & 8)
        {
            AE_L8X8_IP(x0, castxcc(ae_int8x8, pIn), sizeof(ae_int8x8));
            pOut = (ae_int8x16 *)(out + (w&~15));
            for (n = 0; n<h; n++)
            {
                AE_S8X8_XP(x0, castxcc(ae_int8x8, pOut), sizeof(uint8_t)*ostride);
            }
        }
        for (m = (w&(~7)); m < w; m++)
        {
            tmp = in[m];
            for (n = 0; n < h; n++)
            {
                out[n*ostride + m] = tmp;
            }
        }
    }
#endif
    else
    {
        al_outw = AE_ZALIGN128();
        for (n = 0; n < h; n++)
        {
            pOut = (ae_int8x16 *)(out + n*ostride);
            pIn = (const ae_int8x16 *)in;
            for (m = 0; m < (w >> 4); m++)
            {
                AE_L8X8X2_IP(x0, x1, pIn, 2*sizeof(ae_int8x8));
                AE_SA8X8X2_IP(x0, x1, al_outw, pOut);
            }            
            if (w & 15)
            {
                int off = w & 15;
                al_inw = AE_LA128_PP(pIn);
                AE_LAV8X8X2_XP(x0, x1, al_inw,  pIn,  off);
                AE_SAV8X8X2_XP(x0, x1, al_outw, pOut, off);
            }
            AE_SA128POS_FP(al_outw, pOut);
        }
    }
#endif
}
