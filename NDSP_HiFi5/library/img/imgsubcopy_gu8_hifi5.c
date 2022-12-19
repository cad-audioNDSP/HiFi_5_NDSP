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
    Functon copies one subimage to another subimage
    Input:
    in - input image
    width,height - size
    istride - stride of input image
    ostride - stride of output image
    Output:
    out - output image
------------------------------------------------*/
void imgsubcopy_gu8(void* outImg, const void* inImg, int width,int height, int istride, int ostride)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
          uint8_t * restrict out=(uint8_t * )outImg;
    const uint8_t * restrict in =(const uint8_t * )inImg;
    ae_int8x16 * restrict pOut;
    const ae_int8x16 * restrict pIn;
    ae_valign al_in, al_out;
    ae_valignx2 al_inw, al_outw;
    ae_int8x8 x0, x1;
    int m,n;
    NASSERT(outImg );
    NASSERT(inImg );
    if ((width <= 0) || (height <= 0)) return;
    //perfect alignment
    if (!(istride & 15) && !(ostride & 15) && !((uintptr_t)in & 15) && !((uintptr_t)out & 15))
    {
        for (n = 0; n < height; n++)
        {
            pIn = (const ae_int8x16 *)(in + n*istride);
            pOut = (ae_int8x16 *)(out + n*ostride);
            for (m = 0; m < (width >> 4); m++)
            {
                AE_L8X8X2_IP(x0, x1, pIn,  2*sizeof(ae_int8x8));
                AE_S8X8X2_IP(x0, x1, pOut, 2*sizeof(ae_int8x8));
            }
            if (width & 8)
            {
                AE_L8X8_IP(x0, castxcc(ae_int8x8, pIn ),sizeof(ae_int8x8));
                AE_S8X8_IP(x0, castxcc(ae_int8x8, pOut),sizeof(ae_int8x8));
            }
        } 
    }
    else
    {
        al_out = AE_ZALIGN64();
        al_outw = AE_ZALIGN128();
        for (n = 0; n < height; n++)
        {
            pIn = (const ae_int8x16 *)(in + n*istride);
            pOut = (ae_int8x16 *)(out + n*ostride);
            al_inw = AE_LA128_PP(pIn);
            for (m = 0; m < (width >> 4); m++)
            {
                AE_LA8X8X2_IP(x0, x1, al_inw,  pIn);
                AE_SA8X8X2_IP(x0, x1, al_outw, pOut);
            }
            AE_SA128POS_FP(al_outw, pOut);
            if (width & 8)
            {
                al_in = AE_LA64_PP(pIn);
                AE_LA8X8_IP(x0, al_in,  castxcc(ae_int8x8, pIn ));
                AE_SA8X8_IP(x0, al_out, castxcc(ae_int8x8, pOut));
                AE_SA64POS_FP(al_out, pOut);
            }           
        }
    }
    for (n = 0; n < height; n++)
    {
        for (m = width&(~7); m < width; m++)
        {
            out[n*ostride + m] = in[n*istride + m];
        }
    }
#else
    uint8_t * restrict out=(uint8_t * )outImg;
    const uint8_t * restrict in =(const uint8_t * )inImg;
    ae_int8x16 * restrict pOut;
    const ae_int8x16 * restrict pIn;
    ae_valignx2 al_inw, al_outw;
    ae_int8x8 x0, x1;
    int m,n;
    NASSERT(outImg );
    NASSERT(inImg );
    if ((width <= 0) || (height <= 0)) return;
    //perfect alignment
    if (!(istride & 15) && !(ostride & 15) && !((uintptr_t)in & 15) && !((uintptr_t)out & 15))
    {
        al_outw = AE_ZALIGN128();
        for (n = 0; n < height; n++)
        {
            pIn = (const ae_int8x16 *)(in + n*istride);
            pOut = (ae_int8x16 *)(out + n*ostride);
            for (m = 0; m < (width >> 4); m++)
            {
                AE_L8X8X2_IP(x0, x1, pIn,  2*sizeof(ae_int8x8));
                AE_S8X8X2_IP(x0, x1, pOut, 2*sizeof(ae_int8x8));
            }
            if (width & 15)
            {
                al_inw = AE_LA128_PP(pIn);
                AE_LAV8X8X2_XP(x0, x1, al_inw,  pIn , (width & 15));
                AE_SAV8X8X2_XP(x0, x1, al_outw, pOut, (width & 15));
                AE_SA128POS_FP(al_outw, pOut);
            }           
        } 
    }
    else
    {
        al_outw = AE_ZALIGN128();
        for (n = 0; n < height; n++)
        {
            pIn = (const ae_int8x16 *)(in + n*istride);
            pOut = (ae_int8x16 *)(out + n*ostride);
            al_inw = AE_LA128_PP(pIn);
            for (m = 0; m < (width >> 4); m++)
            {
                AE_LA8X8X2_IP(x0, x1, al_inw,  pIn);
                AE_SA8X8X2_IP(x0, x1, al_outw, pOut);
            }            
            if (width & 15)
            {
                AE_LAV8X8X2_XP(x0, x1, al_inw,  pIn , (width & 15));
                AE_SAV8X8X2_XP(x0, x1, al_outw, pOut, (width & 15));
            } 
            AE_SA128POS_FP(al_outw, pOut);          
        }
    }
#endif
}
