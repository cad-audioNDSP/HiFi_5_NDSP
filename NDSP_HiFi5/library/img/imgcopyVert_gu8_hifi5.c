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
    Functon copies vertical edges to rectange
    Input:
    edge[h] - vertical edge 
    w,h,ostride - copied region
    Output:
    outImg   - output image
------------------------------------------------*/
void imgcopyVert_gu8(void* outImg, const void* edge, int w, int h, int ostride)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
          uint8_t * restrict out=(uint8_t * )outImg;
    const ae_int8 * restrict in =(const ae_int8 * )edge;
    int m,n;
    ae_int8x8 x;
    ae_int8x16 * restrict pOut;
    ae_int8 * restrict out_;
    ae_valign al_out;
    ae_valignx2 al_outw;
    NASSERT(edge );
    NASSERT_ALIGN(edge,ALIGNMENT);
    if (w <= 0 || h <= 0) return;
    if (!(ostride & 15) && !((uintptr_t)out & 15))
    {
        for (n = 0; n < h; n++)
        {
            AE_L8_IP(x, in, sizeof(uint8_t));
            pOut = (ae_int8x16 *)(out + n*ostride);
            for (m = 0; m < (w >> 4); m++)
            {
                AE_S8X8X2_IP(x, x, pOut, 2 * sizeof(ae_int8x8));
            }
            if (w & 8)
            {
                AE_S8X8_IP(x, castxcc(ae_int8x8, pOut), sizeof(ae_int8x8));
            }
            out_ = (ae_int8 *)(pOut);
            for (m = (w&(~7)); m < w; m++)
            {
                AE_S8_0_IP(x, out_, sizeof(uint8_t));
            }
        }
    }
    else
    {
        al_out = AE_ZALIGN64();
        al_outw = AE_ZALIGN128();
        for (n = 0; n < h; n++)
        {
            AE_L8_IP(x, in, sizeof(uint8_t));
            pOut = (ae_int8x16 *)(out + n*ostride);
            for (m = 0; m < (w >> 4); m++)
            {
                AE_SA8X8X2_IP(x, x, al_outw, pOut);
            }
            AE_SA128POS_FP(al_outw, pOut);
            if (w & 8)
            {
                AE_SA8X8_IP(x, al_out, castxcc(ae_int8x8, pOut));
                AE_SA64POS_FP(al_out, pOut);
            }
            out_ = (ae_int8 *)(pOut);
            for (m = (w&(~7)); m < w; m++)
            {
                AE_S8_0_IP(x, out_, sizeof(uint8_t));
            }
        }
    }
#else
    uint8_t * restrict out=(uint8_t * )outImg;
    const ae_int8 * restrict in =(const ae_int8 * )edge;
    int m,n;
    ae_int8x8 x;
    ae_int8x16 * restrict pOut;
    ae_valignx2 al_outw;
    NASSERT(edge );
    NASSERT_ALIGN(edge,ALIGNMENT);
    if (w <= 0 || h <= 0) return;
    if (!(ostride & 15) && !((uintptr_t)out & 15))
    {
        al_outw = AE_ZALIGN128();
        for (n = 0; n < h; n++)
        {
            AE_L8_IP(x, in, sizeof(uint8_t));
            pOut = (ae_int8x16 *)(out + n*ostride);
            for (m = 0; m < (w >> 4); m++)
            {
                AE_S8X8X2_IP(x, x, pOut, 2 * sizeof(ae_int8x8));
            }
            if (w & 15)
            {
                AE_SAV8X8X2_XP(x, x, al_outw, pOut, (w & 15));
                AE_SA128POS_FP(al_outw, pOut);
            }
        }
    }
    else
    {
        al_outw = AE_ZALIGN128();
        for (n = 0; n < h; n++)
        {
            AE_L8_IP(x, in, sizeof(uint8_t));
            pOut = (ae_int8x16 *)(out + n*ostride);
            for (m = 0; m < (w >> 4); m++)
            {
                AE_SA8X8X2_IP(x, x, al_outw, pOut);
            }           
            if (w & 15)
            {
                AE_SAV8X8X2_XP(x, x, al_outw, pOut, (w & 15));
            }
            AE_SA128POS_FP(al_outw, pOut);
        }
    }
#endif
}
