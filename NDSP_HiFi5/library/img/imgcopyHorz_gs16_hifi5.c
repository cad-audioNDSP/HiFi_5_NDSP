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
void imgcopyHorz_gs16(void* restrict outImg, const void* edge, int w, int h, int ostride)
{
#if !(XCHAL_HAVE_HIFI5_NN_MAC)
          int16_t * restrict out=(int16_t * )outImg;
    const int16_t * restrict in =(const int16_t * )edge;
    int m,n;
    ae_int16x8 * restrict pOut;
    const ae_int16x8 * restrict pIn = (const ae_int16x8 *)in;
    ae_int16x4 x0, x1;
    ae_valign al_out;
    ae_valignx2 al_outw;
    int16_t tmp;
    NASSERT(edge );
    NASSERT_ALIGN(edge,ALIGNMENT);
    if ((!(ostride & 7)) && (!((uintptr_t)out & 15)))
    {        
        for (m=0; m<(w>>3); m++) 
        {
            AE_L16X4X2_IP(x0, x1, pIn, 8*sizeof(int16_t)); 
            pOut = (ae_int16x8 *)(out + 8*m);
            for (n = 0; n<h; n++)
            {
                AE_S16X4X2_XP(x0, x1, pOut, sizeof(int16_t)*ostride);
            }
        } 
        if (w &4) 
        {
            AE_L16X4_IP(x0, castxcc(ae_int16x4, pIn), 4*sizeof(int16_t)); 
            pOut = (ae_int16x8 *)(out + (w&~7));
            for (n = 0; n<h; n++)
            {
                AE_S16X4_XP(x0, castxcc(ae_int16x4, pOut), sizeof(int16_t)*ostride);
            }
        }      
    }
    else //not aligned
    {
        al_outw = AE_ZALIGN128();
        al_out = AE_ZALIGN64();
        for (n = 0; n < h; n++)
        {
            pOut = (ae_int16x8 *)(out + n*ostride);
            pIn = (const ae_int16x8 *)in;
            for (m = 0; m < (w >> 3); m++)
            {
                AE_L16X4X2_IP(x0, x1, pIn, 8 * sizeof(int16_t));
                AE_SA16X4X2_IP(x0, x1, al_outw, pOut);
            }
            AE_SA128POS_FP(al_outw, pOut);
            if (w&4)
            {
                AE_L16X4_IP(x0, castxcc(ae_int16x4, pIn), 4 * sizeof(int16_t));
                AE_SA16X4_IP(x0, al_out, castxcc(ae_int16x4, pOut));
                AE_SA64POS_FP(al_out, pOut);
            }
        }
    }
    for (m = (w&(~3)); m < w; m++)
    {
        tmp = in[m];
        for (n = 0; n < h; n++)
        {
            out[n*ostride + m] = tmp;
        }
    }
#else
              int16_t * restrict out=(int16_t * )outImg;
    const int16_t * restrict in =(const int16_t * )edge;
    int m,n;
    ae_int16x8 * restrict pOut;
    const ae_int16x8 * restrict pIn = (const ae_int16x8 *)in;
    ae_int16x4 x0, x1;
    ae_valignx2 al_outw, al_inw;

    NASSERT(edge );
    NASSERT_ALIGN(edge,ALIGNMENT);
    if ((!(ostride & 7)) && (!((uintptr_t)out & 15)))
    {
#if 0        
        al_outw = AE_ZALIGN128();
        for (n = 0; n < h; n++)
        {
            pOut = (ae_int16x8 *)(out + n*ostride);
            pIn = (const ae_int16x8 *)in;
            for (m = 0; m < (w >> 3); m++)
            {
                AE_L16X4X2_IP(x0, x1, pIn, 8 * sizeof(int16_t));
                AE_S16X4X2_IP(x0, x1, pOut,  8 * sizeof(int16_t));
            }           
            if (w&7)
            {
                int off = (w&7)<<1;
                al_inw = AE_LA128_PP(pIn);
                AE_LAV16X4X2_XP(x0, x1, al_inw, pIn, off);
                AE_SAV16X4X2_XP(x0, x1, al_outw, pOut, off);
                AE_SA128POS_FP(al_outw, pOut);
            }
            
        } 
#else
        int16_t tmp;
        for (m=0; m<(w>>3); m++) 
        {
            AE_L16X4X2_IP(x0, x1, pIn, 8*sizeof(int16_t)); 
            pOut = (ae_int16x8 *)(out + 8*m);
            for (n = 0; n<h; n++)
            {
                AE_S16X4X2_XP(x0, x1, pOut, sizeof(int16_t)*ostride);
            }
        } 
        if (w &4) 
        {
            AE_L16X4_IP(x0, castxcc(ae_int16x4, pIn), 4*sizeof(int16_t)); 
            pOut = (ae_int16x8 *)(out + (w&~7));
            for (n = 0; n<h; n++)
            {
                AE_S16X4_XP(x0, castxcc(ae_int16x4, pOut), sizeof(int16_t)*ostride);
            }
        }
        for (m = (w&(~3)); m < w; m++)
        {
            tmp = in[m];
            for (n = 0; n < h; n++)
            {
                out[n*ostride + m] = tmp;
            }
        }      
#endif     
    }
    else //not aligned
    {
        al_outw = AE_ZALIGN128();
        for (n = 0; n < h; n++)
        {
            pOut = (ae_int16x8 *)(out + n*ostride);
            pIn = (const ae_int16x8 *)in;
            for (m = 0; m < (w >> 3); m++)
            {
                AE_L16X4X2_IP(x0, x1, pIn, 8 * sizeof(int16_t));
                AE_SA16X4X2_IP(x0, x1, al_outw, pOut);
            }           
            if (w&7)
            {
                int off = (w&7)<<1;
                al_inw = AE_LA128_PP(pIn);
                AE_LAV16X4X2_XP(x0, x1, al_inw, pIn, off);
                AE_SAV16X4X2_XP(x0, x1, al_outw, pOut, off);
            }
            AE_SA128POS_FP(al_outw, pOut);
        }
    }
    
#endif
}
