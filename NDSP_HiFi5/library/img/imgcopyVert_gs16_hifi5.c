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
void imgcopyVert_gs16(void* outImg, const void* edge, int w, int h, int ostride)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
          int16_t * restrict out=(int16_t * )outImg;
    const ae_int16 * restrict in =(const ae_int16 * )edge;
    ae_int16x8 * restrict pOut;
    ae_int16 * restrict out_;
    ae_int16x4 x;
    ae_valign al_out;
    ae_valignx2 al_outw;
    int m,n;
    NASSERT(edge );
    NASSERT_ALIGN(edge,ALIGNMENT);
    //aligned version
    if (!(ostride & 7) && !((uintptr_t)out & 15))
    {
        for (n = 0; n < h; n++)
        {
            AE_L16_IP(x, in, sizeof(int16_t));
            pOut = (ae_int16x8 *)(out + n*ostride);
            for (m = 0; m < (w >> 3); m++)
            {
                AE_S16X4X2_IP(x, x, pOut, 8 * sizeof(int16_t));
            }
            if (w & 4)
            {
                AE_S16X4_IP(x, castxcc(ae_int16x4, pOut), 4 * sizeof(int16_t));
            }
            out_ = (ae_int16 *)(pOut);
            for (m = (w&(~3)); m < w; m++)
            {
                AE_S16_0_IP(x, out_, sizeof(int16_t));
            }
        }
    }
    else //not aligned
    {
        al_out = AE_ZALIGN64();
        al_outw = AE_ZALIGN128();
        for (n = 0; n < h; n++)
        {
            AE_L16_IP(x, in, sizeof(int16_t));
            pOut = (ae_int16x8 *)(out + n*ostride);
            for (m = 0; m < (w >> 3); m++)
            {
                AE_SA16X4X2_IP(x, x, al_outw, pOut);
            }
            AE_SA128POS_FP(al_outw, pOut);
            if (w & 4)
            {
                AE_SA16X4_IP(x, al_out, castxcc(ae_int16x4, pOut));
                AE_SA64POS_FP(al_out, pOut);
            }
            out_ = (ae_int16 *)(out + n*ostride + (w&(~3)));
            for (m = (w&(~3)); m < w; m++)
            {
                AE_S16_0_IP(x, out_, sizeof(int16_t));
            }
        }
    }
#else
    int16_t * restrict out=(int16_t * )outImg;
    const ae_int16 * restrict in =(const ae_int16 * )edge;
    ae_int16x8 * restrict pOut;
    ae_int16x4 x;
    ae_valignx2 al_outw;
    int m,n;
    NASSERT(edge );
    NASSERT_ALIGN(edge,ALIGNMENT);
    //aligned version
    if (!(ostride & 7) && !((uintptr_t)out & 15))
    {
        al_outw = AE_ZALIGN128();
        for (n = 0; n < h; n++)
        {
            AE_L16_IP(x, in, sizeof(int16_t));
            pOut = (ae_int16x8 *)(out + n*ostride);
            for (m = 0; m < (w >> 3); m++)
            {
                AE_S16X4X2_IP(x, x, pOut, 8 * sizeof(int16_t));
            }
            if (w & 7)
            {
                AE_SAV16X4X2_XP(x, x, al_outw, pOut, ((w&7)<<1));
                AE_SA128POS_FP(al_outw, pOut);
            }           
        }
    }
    else //not aligned
    {
        al_outw = AE_ZALIGN128();
        for (n = 0; n < h; n++)
        {
            AE_L16_IP(x, in, sizeof(int16_t));
            pOut = (ae_int16x8 *)(out + n*ostride);
            for (m = 0; m < (w >> 3); m++)
            {
                AE_SA16X4X2_IP(x, x, al_outw, pOut);
            }            
            if (w & 7)
            {
                AE_SAV16X4X2_XP(x, x, al_outw, pOut, ((w&7)<<1));
            }
            AE_SA128POS_FP(al_outw, pOut);
        }
    }
#endif
}
