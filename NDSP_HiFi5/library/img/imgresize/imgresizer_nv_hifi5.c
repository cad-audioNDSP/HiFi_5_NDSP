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
#include "imgresizer_nearest_common.h"

/*    image resizer, nearest, vertical direction */

static void getIdx(int16_t* ix, int M, int N)
#if 1
{
    ae_int16x4* pIx=(ae_int16x4*)ix;
    ae_int32x2 x0,x1,inc,y0,y1,e,z;
    int32_t x;
    int n,s,s2N1M;
    NASSERT_ALIGN(ix,ALIGNMENT);
    x=N; s=NSA(x); x<<=s;
    s2N1M=NSA((2*N-1)*M);
    s=(31-s)+s2N1M;
    // compute 1/N
    /* first approximation */
    z=AE_SUB32(0xBAEC0000,x);
    /* 4 iterations to achieve 1 LSB for reciprocal */
    e=0x40000000; AE_MULSFP32X2RAS(e,x,z); AE_MULAFP32X2RAS(z,z,AE_ADD32(e,e));
    e=0x40000000; AE_MULSFP32X2RAS(e,x,z); AE_MULAFP32X2RAS(z,z,AE_ADD32(e,e));
    e=0x40000000; AE_MULSFP32X2RAS(e,x,z); AE_MULAFP32X2RAS(z,z,AE_ADD32(e,e));
    e=0x40000000; AE_MULSFP32X2RAS(e,x,z); AE_MULAFP32X2RAS(z,z,AE_ADD32(e,e));
    x0=AE_MOVDA32X2(1,3);
    x1=AE_MOVDA32X2(5,7);
    x0=AE_SLAA32S(AE_MULP32X2(x0,M),s2N1M);
    x1=AE_SLAA32S(AE_MULP32X2(x1,M),s2N1M);
    inc=AE_SLAA32S(M<<3,s2N1M);
    for (n=0; n<((N+3)&~3); n+=4)
    {
        ae_int16x4 res;
        y0=AE_SRAA32S(AE_MULFP32X2RAS(x0,z),s);
        y1=AE_SRAA32S(AE_MULFP32X2RAS(x1,z),s);
        y0=AE_MIN32(M-1,AE_MAX32(0,y0));
        y1=AE_MIN32(M-1,AE_MAX32(0,y1));
        res=AE_SAT16X4(y0,y1);
        AE_S16X4_IP(res,pIx,sizeof(ae_int16x4));

        x0=AE_ADD32(x0,inc);
        x1=AE_ADD32(x1,inc);
    }
}
#else
{
    int m,n;
    for (n=0; n<N; n++)
    {
        float32_t u;
        u=((n+0.5f)*M)/N;
        m=(int)floorf(u);
        m=MIN(M-1,MAX(0,m));
        ix[n]=m;
    }
}
#endif

/* returns size of coefficients */
static size_t getCoefSz(const imgsize_t* in,const imgsize_t* out)
{
    int hin=in->height,hout=out->height;
    (void)hin,(void)hout;
    (void)in,(void)out;
    imgsize_validate(in,2,0);
    imgsize_validate(out,2,0);
    NASSERT(in->width == out->width);
    return sizeof(int16_t)*hout;
}
/* returns coefficients */
static void getCoef(void* coef, const imgsize_t* in,const imgsize_t* out)
{
    int hin=in->height,hout=out->height;
    (void)in,(void)out;
    imgsize_validate(in,2,0);
    imgsize_validate(out,2,0);
    NASSERT(in->width == out->width);
    getIdx((int16_t*)coef,hin,hout);
}

static size_t getScratchSize(const imgsize_t* in,const imgsize_t* out)
{
    (void)in,(void)out;
    imgsize_validate(in,2,0);
    imgsize_validate(out,2,0);
    return 0;
}
static void process_gs16(void *pScr, void* pCoef, const void* imgIn, void*imgOut, const imgsize_t* in, const imgsize_t* out);
static void process_fast_gs(void *pScr, void* pCoef, const void* imgIn, void*imgOut, const imgsize_t* in, const imgsize_t* out);
static void process_gs(void *pScr, void* pCoef, const void* imgIn, void*imgOut, const imgsize_t* in, const imgsize_t* out);

/* not in-place image resize */
static void process(void *pScr, void* pCoef, const void* imgIn, void*imgOut, const imgsize_t* in, const imgsize_t* out, int fast)
{
    const ae_int16x8 * restrict pIn;
          ae_int16x8 * restrict pOut;
    const int16_t * restrict ix = (const int16_t*)pCoef;
    ae_int16x4 x0,x1;
    ae_valignx2 al_out;
    int m, n,
        wout = out->width,
        hin = in->height,
        hout = out->height,
        ostride = out->stride,
        istride = in->stride;

    if (fast == 2)
    {
        return process_fast_gs(pScr, pCoef, imgIn, imgOut, in, out);
    }
    if (fast == 3)
    {
        return process_gs(pScr, pCoef, imgIn, imgOut, in, out);
    }
    if (fast == 1)
    {
        return process_gs16(pScr, pCoef, imgIn, imgOut, in, out);
    }

    (void)hin;
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    NASSERT_ALIGN(imgIn, ALIGNMENT);
    NASSERT_ALIGN(imgOut, ALIGNMENT);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->width == wout);// && in->stride == out->stride);

    ostride = ostride * sizeof(int16_t);
    istride = istride * sizeof(int16_t);

    al_out = AE_ZALIGN128();

    if (hin < hout)
    {
        /* Upsampling */
        /* Go from the last pixels to the first in the column */
        for (n = hout - 1; n >= 0; n--)
        {                    
            pIn = (const ae_int16x8 *)((uintptr_t)imgIn + (ix[n] * istride));
            pOut = (ae_int16x8 *)((uintptr_t)imgOut + n*ostride);
            for (m = 0; m < (wout >> 3); m++)
            {
                AE_L16X4X2_IP(x0, x1, pIn , 8 * sizeof(int16_t));
                AE_S16X4X2_IP(x0, x1, pOut, 8 * sizeof(int16_t));
            }
#if  defined(AE_LAV16X4X2_XP)& defined(AE_SAV8X8X2_XP)
            if (wout & 7)
            {
                ae_valignx2 al_in;
                al_in = AE_LA128_PP(pIn);
                AE_LAV16X4X2_XP(x0, x1, al_in, pIn, (wout & 7) * sizeof(int16_t));
                AE_SAV16X4X2_XP(x0, x1, al_out, pOut, (wout & 7) * sizeof(int16_t));
                AE_SA128POS_FP(al_out, pOut);
            }
#else
            if (wout & 4)
            {

                AE_L16X4_IP(x0, castxcc(ae_int16x4, pIn ), 4 * sizeof(int16_t));
                AE_S16X4_IP(x0, castxcc(ae_int16x4, pOut), 4 * sizeof(int16_t));
            }
            if (wout & 3)
            {
                int16_t* pout = (int16_t*)pOut;
                AE_L16X4_IP(x0, castxcc(ae_int16x4, pIn), 4 * sizeof(int16_t));
                pout[0] = AE_MOVAD16_3(x0);
                if (wout & 2)
                {
                    pout[1] = AE_MOVAD16_2(x0);
                }
                if ((wout & 3) == 3)
                {
                    pout[2] = AE_MOVAD16_1(x0);
                }
            }
#endif
        }
    }
    else
    {
        /* Downsampling */
        /* Go from the first pixels to the last in the column */
        for (n = 0; n < hout; n++)
        {
            pIn = (const ae_int16x8*)((uintptr_t)imgIn + (ix[n] * istride));
            pOut = (ae_int16x8*)((uintptr_t)imgOut + n * ostride);
            for (m = 0; m < (wout >> 3); m++)
            {
                AE_L16X4X2_IP(x0, x1, pIn, 8 * sizeof(int16_t));
                AE_S16X4X2_IP(x0, x1, pOut, 8 * sizeof(int16_t));
            }
#if  defined(AE_LAV16X4X2_XP)& defined(AE_SAV8X8X2_XP)
            if (wout & 7)
            {
                ae_valignx2 al_in;
                al_in = AE_LA128_PP(pIn);
                AE_LAV16X4X2_XP(x0, x1, al_in, pIn, (wout & 7) * sizeof(int16_t));
                AE_SAV16X4X2_XP(x0, x1, al_out, pOut, (wout & 7) * sizeof(int16_t));
                AE_SA128POS_FP(al_out, pOut);
            }
#else
            if (wout & 4)
            {

                AE_L16X4_IP(x0, castxcc(ae_int16x4, pIn), 4 * sizeof(int16_t));
                AE_S16X4_IP(x0, castxcc(ae_int16x4, pOut), 4 * sizeof(int16_t));
            }
            if (wout & 3)
            {
                int16_t* pout = (int16_t*)pOut;
                AE_L16X4_IP(x0, castxcc(ae_int16x4, pIn), 4 * sizeof(int16_t));
                pout[0] = AE_MOVAD16_3(x0);
                if (wout & 2)
                {
                    pout[1] = AE_MOVAD16_2(x0);
                }
                if ((wout & 3) == 3)
                {
                    pout[2] = AE_MOVAD16_1(x0);
                }
            }
#endif
        }
    }

}

/* not in-place image resize, specialised for non aligned 16 bit images */
static void process_gs16(void *pScr, void* pCoef, const void* imgIn, void*imgOut, const imgsize_t* in, const imgsize_t* out)
{
    const ae_int16x8 * restrict pIn;
    ae_int16x8 * restrict pOut;
    const int16_t * restrict ix = (const int16_t*)pCoef;
    ae_int16x4 x0,x1;
    ae_valignx2 al_in, al_out;
    int m, n,
        wout = out->width,
        hin = in->height,
        hout = out->height,
        ostride = out->stride,
        istride = in->stride;

    (void)hin;
    imgsize_validate(in, 2, 0);
    imgsize_validate(out, 2, 0);
    NASSERT_ALIGN(imgIn, 2);
    NASSERT_ALIGN(imgOut, 2);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->width == wout);// && in->stride == out->stride);

    ostride = ostride * sizeof(int16_t);
    istride = istride * sizeof(int16_t);
    al_out = AE_ZALIGN128();
    if (hin < hout)
    {
        /* Upsampling */
        /* Go from the last pixels to the first in the column */
        for (n = hout - 1; n >= 0; n--)
        {                    
            pIn = (const ae_int16x8 *)((uintptr_t)imgIn + (ix[n] * istride));
            pOut = (ae_int16x8 *)((uintptr_t)imgOut + n*ostride);
            al_in = AE_LA128_PP(pIn);
            for (m = 0; m < (wout >> 3); m++)
            {
                AE_LA16X4X2_IP(x0, x1, al_in, pIn);
                AE_SA16X4X2_IP(x0, x1, al_out, pOut);
            }
#if  defined(AE_LAV16X4X2_XP)& defined(AE_SAV8X8X2_XP)
            if (wout & 7)
            {
                AE_LAV16X4X2_XP(x0, x1, al_in, pIn, (wout & 7) * sizeof(ae_int16));
                AE_SAV16X4X2_XP(x0, x1, al_out, pOut, (wout & 7) * sizeof(ae_int16));
            }
            AE_SA128POS_FP(al_out, pOut);
#else
            AE_SA128POS_FP(al_out, pOut);
            if (wout & 7)
            {
                ae_valign al_in64;
                al_in64 = AE_LA64_PP(castxcc(ae_int16x4, pIn));
                if (wout & 4)
                {
                    ae_valign al_out64 = AE_ZALIGN64();
                    AE_LA16X4_IP(x0, al_in64, castxcc(ae_int16x4, pIn));
                    AE_SA16X4_IP(x0, al_out64,castxcc(ae_int16x4, pOut));
                    AE_SA64POS_FP(al_out64, pOut);
                }
                if (wout & 3)
                {
                    int16_t* pout = (int16_t*)pOut;
                    AE_LA16X4_IP(x0, al_in64, castxcc(ae_int16x4, pIn));
                    pout[0] = AE_MOVAD16_3(x0);
                    if (wout & 2)
                    {
                        pout[1] = AE_MOVAD16_2(x0);
                    }
                    if ((wout & 3) == 3)
                    {
                        pout[2] = AE_MOVAD16_1(x0);
                    }
                }
            }
#endif
        }
    }
    else
    {
        /* Downsampling */
        /* Go from the first pixels to the last in the column */
        for (n = 0; n < hout; n++)
        {
            pIn = (const ae_int16x8*)((uintptr_t)imgIn + (ix[n] * istride));
            pOut = (ae_int16x8*)((uintptr_t)imgOut + n * ostride);
            al_in = AE_LA128_PP(pIn);
            for (m = 0; m < (wout >> 3); m++)
            {
                AE_LA16X4X2_IP(x0, x1, al_in, pIn);
                AE_SA16X4X2_IP(x0, x1, al_out, pOut);
            }
#if  defined(AE_LAV16X4X2_XP)& defined(AE_SAV8X8X2_XP)
            if (wout & 7)
            {
                AE_LAV16X4X2_XP(x0, x1, al_in, pIn, (wout & 7) * sizeof(ae_int16));
                AE_SAV16X4X2_XP(x0, x1, al_out, pOut, (wout & 7) * sizeof(ae_int16));
            }
            AE_SA128POS_FP(al_out, pOut);
#else
            AE_SA128POS_FP(al_out, pOut);
            if (wout & 7)
            {
                ae_valign al_in64;
                al_in64 = AE_LA64_PP(castxcc(ae_int16x4, pIn));
                if (wout & 4)
                {
                    ae_valign al_out64 = AE_ZALIGN64();
                    AE_LA16X4_IP(x0, al_in64, castxcc(ae_int16x4, pIn));
                    AE_SA16X4_IP(x0, al_out64, castxcc(ae_int16x4, pOut));
                    AE_SA64POS_FP(al_out64, pOut);
                }
                if (wout & 3)
                {
                    int16_t* pout = (int16_t*)pOut;
                    AE_LA16X4_IP(x0, al_in64, castxcc(ae_int16x4, pIn));
                    pout[0] = AE_MOVAD16_3(x0);
                    if (wout & 2)
                    {
                        pout[1] = AE_MOVAD16_2(x0);
                    }
                    if ((wout & 3) == 3)
                    {
                        pout[2] = AE_MOVAD16_1(x0);
                    }
                }
            }
#endif
        }
    }
}

/* 8 bit specific process for fast */
static void process_fast_gs(void *pScr, void* pCoef, const void* imgIn, void*imgOut, const imgsize_t* in, const imgsize_t* out)
{
    const ae_int8x16 * restrict pIn;
    ae_int8x16 * restrict pOut;
    const int16_t * restrict ix = (const int16_t*)pCoef;
    ae_int8x8 x0,x1;
    int m, n,
        wout = out->width,
        hin = in->height,
        hout = out->height,
        ostride = out->stride,
        istride = in->stride;

    

    (void)hin;
    imgsize_validate(in, 1, 1);
    imgsize_validate(out, 1, 1);
    NASSERT_ALIGN(imgIn, ALIGNMENT);
    NASSERT_ALIGN(imgOut, ALIGNMENT);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->width == wout);

    if (hin < hout)
    {
        /* Upsampling */
        /* Go from the last pixels to the first in the column */
        for (n = hout - 1; n >= 0; n--)
        {
            pIn = (const ae_int8x16 *)((uintptr_t)imgIn + (ix[n] * istride));
            pOut = (ae_int8x16 *)((uintptr_t)imgOut + n*ostride);
            for (m = 0; m < (wout >> 4); m++)
            {
                AE_L8X8X2_IP(x0, x1, pIn, sizeof(ae_int8x16));
                AE_S8X8X2_IP(x0, x1, pOut, sizeof(ae_int8x16));
            }
#if  defined(AE_LAV16X4X2_XP)& defined(AE_SAV8X8X2_XP)
            if (wout & 15)
            {
                ae_valignx2 al_in, al_out;
                al_in = AE_LA128_PP(pIn);
                al_out = AE_ZALIGN128();

                AE_LAV8X8X2_XP(x0, x1, al_in, pIn, (wout & 15) * sizeof(int8_t));
                AE_SAV8X8X2_XP(x0, x1, al_out, pOut, (wout & 15) * sizeof(int8_t));
                AE_SA128POS_FP(al_out, pOut);
            }
#else
            if (wout & 8)
            {
                AE_L8X8_IP(x0, castxcc(ae_int8x8, pIn),  sizeof(ae_int8x8));
                AE_S8X8_IP(x0, castxcc(ae_int8x8, pOut), sizeof(ae_int8x8));
            }
            if (wout & 7)
            {
                uint8_t* pout = (uint8_t*)pOut;
                uint8_t* pin = (uint8_t*)pIn;
                for (m = 0; m < (wout & 7); m++)
                {
                    pout[m] = pin[m];
                }
            }
#endif

        }
    }
    else
    {
        /* Downsampling */
        /* Go from the first pixels to the last in the column */
        for (n = 0; n < hout; n++)
        {
            pIn = (const ae_int8x16*)((uintptr_t)imgIn + (ix[n] * istride));
            pOut = (ae_int8x16*)((uintptr_t)imgOut + n * ostride);
            for (m = 0; m < (wout >> 4); m++)
            {
                AE_L8X8X2_IP(x0, x1, pIn, sizeof(ae_int8x16));
                AE_S8X8X2_IP(x0, x1, pOut, sizeof(ae_int8x16));
            }
#if  defined(AE_LAV16X4X2_XP)& defined(AE_SAV8X8X2_XP)
            if (wout & 15)
            {
                ae_valignx2 al_in, al_out;
                al_in = AE_LA128_PP(pIn);
                al_out = AE_ZALIGN128();

                AE_LAV8X8X2_XP(x0, x1, al_in, pIn, (wout & 15) * sizeof(int8_t));
                AE_SAV8X8X2_XP(x0, x1, al_out, pOut, (wout & 15) * sizeof(int8_t));
                AE_SA128POS_FP(al_out, pOut);
            }
#else
            if (wout & 8)
            {
                AE_L8X8_IP(x0, castxcc(ae_int8x8, pIn), sizeof(ae_int8x8));
                AE_S8X8_IP(x0, castxcc(ae_int8x8, pOut), sizeof(ae_int8x8));
            }
            if (wout & 7)
            {
                uint8_t* pout = (uint8_t*)pOut;
                uint8_t* pin = (uint8_t*)pIn;
                for (m = 0; m < (wout & 7); m++)
                {
                    pout[m] = pin[m];
                }
            }
#endif
        }
    }
}

/* 8 bit specific process */
static void process_gs(void *pScr, void* pCoef, const void* imgIn, void*imgOut, const imgsize_t* in, const imgsize_t* out)
{
    const ae_int8x16 * restrict pIn;
    ae_int8x16* restrict pOut;
    ae_valignx2 al_in, al_out;
    const int16_t * restrict ix = (const int16_t*)pCoef;
    ae_int8x8 x0,x1;

    int m, n,
        wout = out->width,
        hin = in->height,
        hout = out->height,
        ostride = out->stride,
        istride = in->stride;


    (void)hin;
    imgsize_validate(in, 1, 0);
    imgsize_validate(out, 1, 0);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->width == wout);
    al_out = AE_ZALIGN128();

    if (hin < hout)
    {
        /* Go from the last pixels to the first in the column */
        for (n = hout - 1; n >= 0; n--)
        {
            pIn = (const ae_int8x16*)((uintptr_t)imgIn + (ix[n] * istride));
            pOut = (ae_int8x16*)((uintptr_t)imgOut + n * ostride);

            al_in = AE_LA128_PP(pIn);
            for (m = 0; m < ((wout) >> 4); m++)
            {
                AE_LA8X8X2_IP(x0, x1, al_in, pIn);
                AE_SA8X8X2_IP(x0, x1, al_out, pOut);
            }
#if  defined(AE_LAV16X4X2_XP)& defined(AE_SAV8X8X2_XP)
            if ((wout) & 15)
            {
                AE_LAV8X8X2_XP(x0, x1, al_in, pIn, (wout & 15) * sizeof(int8_t));
                AE_SAV8X8X2_XP(x0, x1, al_out, pOut, (wout & 15) * sizeof(int8_t));
            }
            AE_SA128POS_FP(al_out, pOut);
#else
            AE_SA128POS_FP(al_out, pOut);
            if (wout & 8)
            {
                ae_valign al_in64;
                ae_valign al_out64 = AE_ZALIGN64();
                al_in64 = AE_LA64_PP(castxcc(ae_int16x4, pIn));
                AE_LA8X8_IP(x0, al_in64, castxcc(ae_int8x8, pIn));
                AE_SA8X8_IP(x0, al_out64, castxcc(ae_int8x8, pOut));
                AE_SA64POS_FP(al_out64, pOut);
            }
            if (wout & 7)
            {
                uint8_t* pout = (uint8_t*)pOut;
                uint8_t* pin = (uint8_t*)pIn;
                for (m = 0; m < (wout & 7); m++)
                {
                    pout[m] = pin[m];
                }
            }
#endif
        }
    } 
    else
    {
        /* Go from the last pixels to the first in the column */
        for ( n = 0; n < hout; n++)
        {
            pIn = (const ae_int8x16*)((uintptr_t)imgIn + (ix[n] * istride));
            pOut = (ae_int8x16*)((uintptr_t)imgOut + n * ostride);

            al_in = AE_LA128_PP(pIn);
            for (m = 0; m < ((wout) >> 4); m++)
            {
                AE_LA8X8X2_IP(x0, x1, al_in, pIn);
                AE_SA8X8X2_IP(x0, x1, al_out, pOut);
            }
#if  defined(AE_LAV16X4X2_XP)& defined(AE_SAV8X8X2_XP)
            if ((wout) & 15)
            {
                AE_LAV8X8X2_XP(x0, x1, al_in, pIn, (wout & 15) * sizeof(int8_t));
                AE_SAV8X8X2_XP(x0, x1, al_out, pOut, (wout & 15) * sizeof(int8_t));
            }
            AE_SA128POS_FP(al_out, pOut);
#else
            AE_SA128POS_FP(al_out, pOut);
            if (wout & 8)
            {
                ae_valign al_in64;
                ae_valign al_out64 = AE_ZALIGN64();
                al_in64 = AE_LA64_PP(castxcc(ae_int16x4, pIn));
                AE_LA8X8_IP(x0, al_in64, castxcc(ae_int8x8, pIn));
                AE_SA8X8_IP(x0, al_out64, castxcc(ae_int8x8, pOut));
                AE_SA64POS_FP(al_out64, pOut);
            }
            if (wout & 7)
            {
                uint8_t* pout = (uint8_t*)pOut;
                uint8_t* pin = (uint8_t*)pIn;
                for (m = 0; m < (wout & 7); m++)
                {
                    pout[m] = pin[m];
                }
            }
#endif
        }
    }

}

const imgresizer_api_t imgresizer_api_nv={NULL,getCoefSz,getCoef,getScratchSize,process};
