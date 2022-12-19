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

/*    image resizer, nearest method, horizontal direction */

#define MAX(a,b) (a>b?a:b)

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
    for (n=0; n<((N+7)&~7); n+=4)
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
    for (n=0; n<N; n++)
    {
        float32_t u;
        u=((n+0.5f)*M)/N;
        x=L_mpy_ll(((2*n+1)*M)<<s2N1M,z);
        m=(int)floorf(u);
        m=MIN(M-1,MAX(0,m));
        ix[n]=m;
    }
}
#endif

/* returns size of coefficients */
static size_t getCoefSz(const imgsize_t* in,const imgsize_t* out)
{
    int win=in->width,wout=out->width;
    NASSERT(in->height == out->height);
    (void)in,(void)out;
    (void)win,(void)wout;
    imgsize_validate(in,2,0);
    imgsize_validate(out,2,0);
    return sizeof(int16_t)*((wout+3)&~3);   // multiple of 4
}
/* returns coefficients */
static void getCoef(void* coef, const imgsize_t* in,const imgsize_t* out)
{
    int win=in->width,wout=out->width;
    NASSERT(in->height == out->height);
    (void)in,(void)out,(void)coef;
    imgsize_validate(in,2,0);
    imgsize_validate(out,2,0);
    getIdx((int16_t*)coef,win,wout);
}

static size_t getScratchSize(const imgsize_t* in,const imgsize_t* out)
{
    (void)in,(void)out;
    imgsize_validate(in,2,0);
    imgsize_validate(out,2,0);
    /* img size + 8 aligned rows*/
    return ((in->height >> 3))* (((in->width + 7) & ~7) * sizeof(ae_int16x8));
}

static void process_gs16(void *pScr, void* pCoef, const void* imgIn, void*imgOut, const imgsize_t* in, const imgsize_t* out);
static void process_fast_gs(void *pScr, void* pCoef, const void* imgIn, void*imgOut, const imgsize_t* in, const imgsize_t* out);
static void process_gs(void *pScr, void* pCoef, const void* imgIn, void*imgOut, const imgsize_t* in, const imgsize_t* out);

/* not in-place image resize */
static void process(void *pScr, void* pCoef, const void* imgIn, void*imgOut, const imgsize_t* in, const imgsize_t* out, int fast)
{
    int16_t* restrict t;
    const int16_t* restrict ix;
    const ae_int16x8* restrict pIn;
    const ae_int16x8* restrict pIn1;
    const ae_int16x8* restrict pT0;
    const ae_int16x8* restrict pT1;
    const ae_int16x8* restrict pT2;
    const ae_int16x8* restrict pT3;
    const ae_int16x8* restrict pT4;
    const ae_int16x8* restrict pT5;
    const ae_int16x8* restrict pT6;
    const ae_int16x8* restrict pT7;
    ae_int16x8* restrict pOut;
    ae_int16x8* restrict pOut1;
    const ae_int16x8* restrict pIdx;

    ae_int16x4 x0, x1, x2, x3, x4, x5, x6, x7;
    ae_int16x4 x0_, x1_, x2_, x3_, x4_, x5_, x6_, x7_;
    ae_int16x4 t0, t1, t2, t3, t0_, t1_, t2_, t3_;
    ae_int16x4 sel0,sel1;
    ae_int16x4 idx0, idx1, ofs;
    int m, n, win, wout, hout, ostride, istride, ostridex2, istridex2;

    (void)win;

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
        //if (win * 2 == wout) return;
        return process_gs16(pScr, pCoef, imgIn, imgOut, in, out);
    }

    t = (int16_t*)pScr;
    ix = (const int16_t*)pCoef;
    win = in->width;
    wout = out->width;
    hout = out->height;
    ostride = out->stride;
    istride = in->stride;
    imgsize_validate(in, 2, 1);
    imgsize_validate(out, 2, 1);
    NASSERT_ALIGN(imgIn, ALIGNMENT);
    NASSERT_ALIGN(imgOut, ALIGNMENT);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->height == out->height);// && in->stride == out->stride);

    if (win == wout * 2)
    {
        sel0 = AE_MOVINT16X4_FROMINT64(0x0700050003000100); // SEL 7531 

        /* 2x decimation */
        for (n = 0; n < (hout & ~1); n+=2)
        {
            pIn = (const ae_int16x8*)((int16_t*)imgIn + n * istride);
            pIn1 = (const ae_int16x8*)((int16_t*)imgIn + (n+1) * istride);
            pOut = (ae_int16x8*)((int16_t*)imgOut + n * ostride);
            pOut1 = (ae_int16x8*)((int16_t*)imgOut + (n+1) * ostride);

            for (m = 0; m < (win >> 4); m++)
            {
                AE_L16X4X2_IP(x0, x1, pIn, sizeof(ae_int16x8));
                AE_L16X4X2_IP(x0_, x1_, pIn1, sizeof(ae_int16x8));
                AE_L16X4X2_IP(x2, x3, pIn, sizeof(ae_int16x8));
                AE_L16X4X2_IP(x2_, x3_, pIn1, sizeof(ae_int16x8));
                AE_DSEL16X4(t0,t2,x0,x1,sel0); 
                AE_DSEL16X4(t0_,t2_,x0_,x1_,sel0); 
                AE_DSEL16X4(t1,t3,x2,x3,sel0);  
                AE_DSEL16X4(t1_,t3_,x2_,x3_,sel0);  
                AE_S16X4X2_IP(t0, t1, pOut, sizeof(ae_int16x8));
                AE_S16X4X2_IP(t0_, t1_, pOut1, sizeof(ae_int16x8));
            }
            if (win & 8)
            {
                AE_L16X4X2_IP(x0, x1, pIn, sizeof(ae_int16x8));
                AE_L16X4X2_IP(x0_, x1_, pIn1, sizeof(ae_int16x8));
                t0 = AE_SEL16_7531(x0, x1);
                t0_ = AE_SEL16_7531(x0_, x1_);
                AE_S16X4_IP(t0, castxcc(ae_int16x4, pOut), sizeof(ae_int16x4));
                AE_S16X4_IP(t0_, castxcc(ae_int16x4, pOut1), sizeof(ae_int16x4));
            }
            if (win & 4)
            {
                AE_L16X4_IP(x0, castxcc(ae_int16x4, pIn), sizeof(ae_int16x4));
                AE_L16X4_IP(x0_, castxcc(ae_int16x4, pIn1), sizeof(ae_int16x4));
                t0 = AE_SEL16_7531(x0, x0);
                t0_ = AE_SEL16_7531(x0_, x0_);
                t0 = AE_SHORTSWAP(t0);
                t0_ = AE_SHORTSWAP(t0_);
                AE_S32_L_IP(AE_MOVINT32X2_FROMINT16X4(t0), castxcc(ae_int32, pOut), 2 * sizeof(int16_t));
                AE_S32_L_IP(AE_MOVINT32X2_FROMINT16X4(t0_), castxcc(ae_int32, pOut1), 2 * sizeof(int16_t));
            }
            if (win & 2)
            {
                AE_L16X4_IP(x0, castxcc(ae_int16x4, pIn), sizeof(ae_int16x4));
                AE_L16X4_IP(x0_, castxcc(ae_int16x4, pIn1), sizeof(ae_int16x4));
                x0 = AE_SHORTSWAP(x0);
                x0_ = AE_SHORTSWAP(x0_);
                AE_S16_0_IP(x0, castxcc(ae_int16, pOut), sizeof(int16_t));
                AE_S16_0_IP(x0_, castxcc(ae_int16, pOut1), sizeof(int16_t));
            }
        }
        if (hout&1)
        {
            pIn = (const ae_int16x8*)((int16_t*)imgIn + n * istride);
            pOut = (ae_int16x8*)((int16_t*)imgOut + n * ostride);

            for (m = 0; m < (win >> 4); m++)
            {
                AE_L16X4X2_IP(x0, x1, pIn, sizeof(ae_int16x8));
                AE_L16X4X2_IP(x2, x3, pIn, sizeof(ae_int16x8));
                AE_DSEL16X4(t0, t2, x0, x1, sel0);
                AE_DSEL16X4(t1, t3, x2, x3, sel0);
                AE_S16X4X2_IP(t0, t1, pOut, sizeof(ae_int16x8));
            }
            if (win & 8)
            {
                AE_L16X4X2_IP(x0, x1, pIn, sizeof(ae_int16x8));
                t0 = AE_SEL16_7531(x0, x1);
                AE_S16X4_IP(t0, castxcc(ae_int16x4, pOut), sizeof(ae_int16x4));
            }
            if (win & 4)
            {
                AE_L16X4_IP(x0, castxcc(ae_int16x4, pIn), sizeof(ae_int16x4));
                t0 = AE_SEL16_7531(x0, x0);
                t0 = AE_SHORTSWAP(t0);
                AE_S32_L_IP(AE_MOVINT32X2_FROMINT16X4(t0), castxcc(ae_int32, pOut), 2 * sizeof(int16_t));
            }
            if (win & 2)
            {
                AE_L16X4_IP(x0, castxcc(ae_int16x4, pIn), sizeof(ae_int16x4));
                x0 = AE_SHORTSWAP(x0);
                AE_S16_0_IP(x0, castxcc(ae_int16, pOut), sizeof(int16_t));
            }
        }

        return;
    }
    else if (win * 2 == wout)
    {

        /* 2x interpolation */
        sel0 = AE_MOVINT16X4_FROMINT64(0x0703070306020602);
        sel1 = AE_MOVINT16X4_FROMINT64(0x0501050104000400);
        for (n = 0; n < hout; n++)
        {

            pIn = (const ae_int16x8*)((int16_t*)imgIn + n * istride + win - 1);
            pOut = (ae_int16x8*)((int16_t*)imgOut + n * ostride + wout - 2);

            for (m = 0; m < (win & 7); m++)
            {
                AE_L16_IP(x0, castxcc(ae_int16, pIn), -(int)sizeof(int16_t));
                AE_S32_L_IP(AE_MOVINT32X2_FROMINT16X4(x0), castxcc(ae_int32, pOut), -2 * (int)sizeof(int16_t));
            }
            pIn = (const ae_int16x8*)((int16_t*)pIn - 7);
            pOut = (ae_int16x8*)((int16_t*)pOut - 6);
            for (m = 0; m < (win >> 3); m++)
            {
                AE_L16X4X2_XP(x0, x1, pIn, -(int)sizeof(ae_int16x8));
                AE_DSEL16X4(t0, t2, x0, x1, sel0);
                AE_DSEL16X4(t1, t3, x0, x1, sel1);
                AE_S16X4X2_XP(t2, t3, pOut, -(int)sizeof(ae_int16x8));
                AE_S16X4X2_XP(t0, t1, pOut, -(int)sizeof(ae_int16x8));
            }
        }
        return;
    }
    
    sel0 = AE_MOVINT16X4_FROMINT64(0x0706050403020100); // SEL7531 + SEL6420
    istridex2 = istride * sizeof(int16_t);
    ostridex2 = ostride * sizeof(int16_t);

    int shift = ((win+7)&~7) * sizeof(ae_int16x8);
    for (n = 0; n < (hout & ~7); n += 8)
    {
        /* Interleave input samples from 8 rows and save them to the scratch */
        pIn = (const ae_int16x8*)((uintptr_t)imgIn + istridex2 * n);
        pOut = (ae_int16x8*)((int16_t*)t + ((win+7)&~7) * n);
        //pOut = (ae_int16x4*)((uintptr_t*)t+(shift));
        for (m = 0; m < win; m += 8)
        {
            AE_L16X4X2_XP(x0, x0_, castxcc(ae_int16x8, pIn), istridex2);
            AE_L16X4X2_XP(x1, x1_, castxcc(ae_int16x8, pIn), istridex2);
            AE_L16X4X2_XP(x2, x2_, castxcc(ae_int16x8, pIn), istridex2);
            AE_L16X4X2_XP(x3, x3_, castxcc(ae_int16x8, pIn), istridex2);
            AE_L16X4X2_XP(x4, x4_, castxcc(ae_int16x8, pIn), istridex2);
            AE_L16X4X2_XP(x5, x5_, castxcc(ae_int16x8, pIn), istridex2);
            AE_L16X4X2_XP(x6, x6_, castxcc(ae_int16x8, pIn), istridex2);
            AE_L16X4X2_XP(x7, x7_, castxcc(ae_int16x8, pIn), -7 * istridex2 + sizeof(ae_int16x8));

            AE_DSEL16X4(t0, t1, x0, x1, sel0);
            AE_DSEL16X4(t2, t3, x2, x3, sel0);
            AE_DSEL16X4(x0, x2, t0, t2, sel0);
            AE_DSEL16X4(x1, x3, t1, t3, sel0);

            AE_DSEL16X4(t0, t1, x4, x5, sel0);
            AE_DSEL16X4(t2, t3, x6, x7, sel0);
            AE_DSEL16X4(x4, x6, t0, t2, sel0);
            AE_DSEL16X4(x5, x7, t1, t3, sel0);

            AE_DSEL16X4(t0, t1, x0_, x1_, sel0);
            AE_DSEL16X4(t2, t3, x2_, x3_, sel0);
            AE_DSEL16X4(x0_, x2_, t0, t2, sel0);
            AE_DSEL16X4(x1_, x3_, t1, t3, sel0);

            AE_DSEL16X4(t0, t1, x4_, x5_, sel0);
            AE_DSEL16X4(t2, t3, x6_, x7_, sel0);
            AE_DSEL16X4(x4_, x6_, t0, t2, sel0);
            AE_DSEL16X4(x5_, x7_, t1, t3, sel0);

            AE_S16X4X2_IP(x0, x4, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x1, x5, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x2, x6, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x3, x7, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x0_, x4_, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x1_, x5_, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x2_, x6_, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
            AE_S16X4X2_IP(x3_, x7_, castxcc(ae_int16x8, pOut), sizeof(ae_int16x8));
        }
    }
    pIn = (ae_int16x8*)t;
    pIdx = (const ae_int16x8*)(ix);
    for (m = 0; m < (wout&~7); m += 8)
    {
        int k0, k1, k2, k3, k4, k5, k6, k7;
        pOut = (ae_int16x8*)((int16_t*)imgOut + m);

        AE_L16X4X2_IP(idx0, idx1, pIdx, sizeof(ae_int16x8));
        ofs = AE_SLAI16S(idx0, 4); /* offset = idx*sizeof(ae_int16x8) */
        k0 = AE_MOVAD16_3(ofs);
        k1 = AE_MOVAD16_2(ofs);
        k2 = AE_MOVAD16_1(ofs);
        k3 = AE_MOVAD16_0(ofs);
        ofs = AE_SLAI16S(idx1, 4); /* offset = idx*sizeof(ae_int16x8) */
        k4 = AE_MOVAD16_3(ofs);
        k5 = AE_MOVAD16_2(ofs);
        k6 = AE_MOVAD16_1(ofs);
        k7 = AE_MOVAD16_0(ofs);

        pT0 = (ae_int16x8*)((uintptr_t)pIn + k0);
        pT1 = (ae_int16x8*)((uintptr_t)pIn + k1);
        pT2 = (ae_int16x8*)((uintptr_t)pIn + k2);
        pT3 = (ae_int16x8*)((uintptr_t)pIn + k3);
        pT4 = (ae_int16x8*)((uintptr_t)pIn + k4);
        pT5 = (ae_int16x8*)((uintptr_t)pIn + k5);
        pT6 = (ae_int16x8*)((uintptr_t)pIn + k6);
        pT7 = (ae_int16x8*)((uintptr_t)pIn + k7);


        for (n = 0; n < (hout>>3) ; n++)
        {

            AE_L16X4X2_XP(x0, x4, pT0, shift);
            AE_L16X4X2_XP(x1, x5, pT1, shift);
            AE_L16X4X2_XP(x2, x6, pT2, shift);
            AE_L16X4X2_XP(x3, x7, pT3, shift);
  
            AE_L16X4X2_XP(x0_, x4_, pT4, shift);
            AE_L16X4X2_XP(x1_, x5_, pT5, shift);
            AE_L16X4X2_XP(x2_, x6_, pT6, shift);
            AE_L16X4X2_XP(x3_, x7_, pT7, shift);

            AE_DSEL16X4(t0, t1, x0, x1, sel0);
            AE_DSEL16X4(t2, t3, x2, x3, sel0);
            AE_DSEL16X4(x0, x2, t0, t2, sel0);
            AE_DSEL16X4(x1, x3, t1, t3, sel0);

            AE_DSEL16X4(t0, t1, x4, x5, sel0);
            AE_DSEL16X4(t2, t3, x6, x7, sel0);
            AE_DSEL16X4(x4, x6, t0, t2, sel0);
            AE_DSEL16X4(x5, x7, t1, t3, sel0);

            AE_DSEL16X4(t0, t1, x0_, x1_, sel0);
            AE_DSEL16X4(t2, t3, x2_, x3_, sel0);
            AE_DSEL16X4(x0_, x2_, t0, t2, sel0);
            AE_DSEL16X4(x1_, x3_, t1, t3, sel0);

            AE_DSEL16X4(t0, t1, x4_, x5_, sel0);
            AE_DSEL16X4(t2, t3, x6_, x7_, sel0);
            AE_DSEL16X4(x4_, x6_, t0, t2, sel0);
            AE_DSEL16X4(x5_, x7_, t1, t3, sel0);

            AE_S16X4X2_XP(x0, x0_, pOut, ostridex2);
            AE_S16X4X2_XP(x1, x1_, pOut, ostridex2);
            AE_S16X4X2_XP(x2, x2_, pOut, ostridex2);
            AE_S16X4X2_XP(x3, x3_, pOut, ostridex2);
            AE_S16X4X2_XP(x4, x4_, pOut, ostridex2);
            AE_S16X4X2_XP(x5, x5_, pOut, ostridex2);
            AE_S16X4X2_XP(x6, x6_, pOut, ostridex2);
            AE_S16X4X2_XP(x7, x7_, pOut, ostridex2);
        }
    }
    m = wout & ~7;
    if  (wout & 4)
    {
        int k0, k1, k2, k3;
        pOut = (ae_int16x8*)((int16_t*)imgOut + m);

        AE_L16X4_IP(idx0, castxcc(ae_int16x4, pIdx), sizeof(ae_int16x4));
        ofs = AE_SLAI16S(idx0, 4); /* offset = idx*sizeof(ae_int16x8) */
        k0 = AE_MOVAD16_3(ofs);
        k1 = AE_MOVAD16_2(ofs);
        k2 = AE_MOVAD16_1(ofs);
        k3 = AE_MOVAD16_0(ofs);

        pT0 = (ae_int16x8*)((uintptr_t)pIn + k0);
        pT1 = (ae_int16x8*)((uintptr_t)pIn + k1);
        pT2 = (ae_int16x8*)((uintptr_t)pIn + k2);
        pT3 = (ae_int16x8*)((uintptr_t)pIn + k3);

        for (n = 0; n < (hout >> 3); n++)
        {

            AE_L16X4X2_XP(x0, x4, pT0, shift);
            AE_L16X4X2_XP(x1, x5, pT1, shift);
            AE_L16X4X2_XP(x2, x6, pT2, shift);
            AE_L16X4X2_XP(x3, x7, pT3, shift);

            AE_DSEL16X4(t0, t1, x0, x1, sel0);
            AE_DSEL16X4(t2, t3, x2, x3, sel0);
            AE_DSEL16X4(x0, x2, t0, t2, sel0);
            AE_DSEL16X4(x1, x3, t1, t3, sel0);

            AE_DSEL16X4(t0, t1, x4, x5, sel0);
            AE_DSEL16X4(t2, t3, x6, x7, sel0);
            AE_DSEL16X4(x4, x6, t0, t2, sel0);
            AE_DSEL16X4(x5, x7, t1, t3, sel0);

            AE_DSEL16X4(t0, t1, x0_, x1_, sel0);
            AE_DSEL16X4(t2, t3, x2_, x3_, sel0);
            AE_DSEL16X4(x0_, x2_, t0, t2, sel0);
            AE_DSEL16X4(x1_, x3_, t1, t3, sel0);

            AE_DSEL16X4(t0, t1, x4_, x5_, sel0);
            AE_DSEL16X4(t2, t3, x6_, x7_, sel0);
            AE_DSEL16X4(x4_, x6_, t0, t2, sel0);
            AE_DSEL16X4(x5_, x7_, t1, t3, sel0);

            AE_S16X4_XP(x0, castxcc(ae_int16x4,pOut), ostridex2);
            AE_S16X4_XP(x1, castxcc(ae_int16x4,pOut), ostridex2);
            AE_S16X4_XP(x2, castxcc(ae_int16x4,pOut), ostridex2);
            AE_S16X4_XP(x3, castxcc(ae_int16x4,pOut), ostridex2);
            AE_S16X4_XP(x4, castxcc(ae_int16x4,pOut), ostridex2);
            AE_S16X4_XP(x5, castxcc(ae_int16x4,pOut), ostridex2);
            AE_S16X4_XP(x6, castxcc(ae_int16x4,pOut), ostridex2);
            AE_S16X4_XP(x7, castxcc(ae_int16x4,pOut), ostridex2);
        }

        m += 4;
    }
    if (wout & 3)
    {
        int k0, k1, k2, k3;
        pOut = (ae_int16x8*)((int16_t*)imgOut + m);

        AE_L16X4_IP(idx0, castxcc(ae_int16x4, pIdx), sizeof(ae_int16x4));
        ofs = AE_SLAI16S(idx0, 4); /* offset = idx*sizeof(ae_int16x8) */
        k0 = AE_MOVAD16_3(ofs);
        k1 = AE_MOVAD16_2(ofs);
        k2 = AE_MOVAD16_1(ofs);
        k3 = AE_MOVAD16_0(ofs);

        pT0 = (ae_int16x8*)((uintptr_t)pIn + k0);
        pT1 = (ae_int16x8*)((uintptr_t)pIn + k1);
        pT2 = (ae_int16x8*)((uintptr_t)pIn + k2);
        pT3 = (ae_int16x8*)((uintptr_t)pIn + k3);
        int16_t* pout = (int16_t*)pOut;
        for (n = 0; n < (hout >> 3); n++)
        {
            AE_L16X4X2_XP(x0, x4, pT0, shift);
            AE_L16X4X2_XP(x1, x5, pT1, shift);
            AE_L16X4X2_XP(x2, x6, pT2, shift);
            AE_L16X4X2_XP(x3, x7, pT3, shift);

            AE_DSEL16X4(t0, t1, x0, x1, sel0);
            AE_DSEL16X4(t2, t3, x2, x3, sel0);
            AE_DSEL16X4(x0, x2, t0, t2, sel0);
            AE_DSEL16X4(x1, x3, t1, t3, sel0);

            AE_DSEL16X4(t0, t1, x4, x5, sel0);
            AE_DSEL16X4(t2, t3, x6, x7, sel0);
            AE_DSEL16X4(x4, x6, t0, t2, sel0);
            AE_DSEL16X4(x5, x7, t1, t3, sel0);

            pout[0 * ostride + 0] = AE_MOVAD16_3(x0);
            pout[1 * ostride + 0] = AE_MOVAD16_3(x1);
            pout[2 * ostride + 0] = AE_MOVAD16_3(x2);
            pout[3 * ostride + 0] = AE_MOVAD16_3(x3);
            pout[4 * ostride + 0] = AE_MOVAD16_3(x4);
            pout[5 * ostride + 0] = AE_MOVAD16_3(x5);
            pout[6 * ostride + 0] = AE_MOVAD16_3(x6);
            pout[7 * ostride + 0] = AE_MOVAD16_3(x7);
            if (wout & 2)
            {
                pout[0 * ostride + 1] = AE_MOVAD16_2(x0);
                pout[1 * ostride + 1] = AE_MOVAD16_2(x1);
                pout[2 * ostride + 1] = AE_MOVAD16_2(x2);
                pout[3 * ostride + 1] = AE_MOVAD16_2(x3);
                pout[4 * ostride + 1] = AE_MOVAD16_2(x4);
                pout[5 * ostride + 1] = AE_MOVAD16_2(x5);
                pout[6 * ostride + 1] = AE_MOVAD16_2(x6);
                pout[7 * ostride + 1] = AE_MOVAD16_2(x7);
            }
            if ((wout & 3) == 3)
            {
                pout[0 * ostride + 2] = AE_MOVAD16_1(x0);
                pout[1 * ostride + 2] = AE_MOVAD16_1(x1);
                pout[2 * ostride + 2] = AE_MOVAD16_1(x2);
                pout[3 * ostride + 2] = AE_MOVAD16_1(x3);
                pout[4 * ostride + 2] = AE_MOVAD16_1(x4);
                pout[5 * ostride + 2] = AE_MOVAD16_1(x5);
                pout[6 * ostride + 2] = AE_MOVAD16_1(x6);
                pout[7 * ostride + 2] = AE_MOVAD16_1(x7);
            }
            pout += 8 * ostride;
        }
    }

     /* Process last 0...7 rows */
    sel0 = AE_MOVINT16X4_FROMINT64(0x0505040401010000);
    sel1 = AE_MOVINT16X4_FROMINT64(0x0505010104040000);
    for (n = (hout & ~7); n < hout; n++)
    {
        pIn = (const ae_int16x8*)((uintptr_t)imgIn + n * istridex2);
        
        pOut = (ae_int16x8*)(t);
    
        __Pragma("loop_count min=1");
        for (m = 0; m < win; m += 8)
        {
            AE_L16X4X2_IP(x0, x1, pIn, sizeof(ae_int16x8));
            AE_S16X4X2_IP(x0, x1, pOut, sizeof(ae_int16x8));
        }
    
        pIn = (const ae_int16x8*)(t);
        pOut = (ae_int16x8*)((uintptr_t)imgOut + n * ostridex2);
        
        pIdx = (const ae_int16x8*)(ix);
        for (m = 0; m < (wout >> 3); m++)
        {
            /* Load 1x4 elems for the 8 columns */
            AE_L16X4X2_IP(idx0, idx1, pIdx, sizeof(ae_int16x8));
    
            ofs = AE_SLAI16S(idx0, 1); /* offset = idx*sizeof(int16_t) */
            x0 = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_3(ofs));
            x1 = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_2(ofs));
            x2 = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_1(ofs));
            x3 = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_0(ofs));
    
            ofs = AE_SLAI16S(idx1, 1); /* offset = idx*sizeof(int16_t) */
            x0_ = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_3(ofs));
            x1_ = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_2(ofs));
            x2_ = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_1(ofs));
            x3_ = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_0(ofs));
    
            AE_DSEL16X4(x0,t0,x0,x1,sel1);
            AE_DSEL16X4(x2,t1,x2,x3,sel1);
            AE_DSEL16X4(x0,t2,x0,x2,sel0);
            
            AE_DSEL16X4(x0_, t0, x0_, x1_, sel1);
            AE_DSEL16X4(x2_, t1, x2_, x3_, sel1);
            AE_DSEL16X4(x0_, t2, x0_, x2_, sel0);
    
            AE_S16X4X2_IP(x0, x0_, pOut, sizeof(ae_int16x8));
        }
        if (wout & 4)
        {
            AE_L16X4_IP(idx0, castxcc(ae_int16x4, pIdx), sizeof(ae_int16x4));

            ofs = AE_SLAI16S(idx0, 1); /* offset = idx*sizeof(int16_t) */
            x0 = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_3(ofs));
            x1 = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_2(ofs));
            x2 = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_1(ofs));
            x3 = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_0(ofs));

            AE_DSEL16X4(x0, t0, x0, x1, sel1);
            AE_DSEL16X4(x2, t1, x2, x3, sel1);
            AE_DSEL16X4(x0, t2, x0, x2, sel0);

            AE_S16X4_IP(x0, castxcc(ae_int16x4,pOut), sizeof(ae_int16x4));

        }
        if (wout & 3)
        {

            int16_t* pout = (int16_t*)pOut;

            AE_L16X4_IP(idx0, castxcc(ae_int16x4, pIdx), sizeof(ae_int16x4));

            ofs = AE_SLAI16S(idx0, 1); /* offset = idx*sizeof(int16_t) */
            x0 = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_3(ofs));
            x1 = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_2(ofs));
            x2 = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_1(ofs));
            x3 = AE_L16_X((const ae_int16*)pIn, AE_MOVAD16_0(ofs));

            AE_DSEL16X4(x0, t0, x0, x1, sel1);
            AE_DSEL16X4(x2, t1, x2, x3, sel1);
            AE_DSEL16X4(x0, t2, x0, x2, sel0);

            pout[0 * ostride + 0] = AE_MOVAD16_3(x0);
            if (wout & 2)
            {
                pout[0 * ostride + 1] = AE_MOVAD16_2(x0);
            }
            if ((wout & 3) == 3)
            {
                pout[0 * ostride + 2] = AE_MOVAD16_1(x0);
            }
        }
    }
}

/* not in-place image resize */
static void process_gs16(void *pScr, void* pCoef, const void* imgIn, void*imgOut, const imgsize_t* in, const imgsize_t* out)
{
    int16_t *t = (int16_t *)pScr;
    const int16_t * restrict  ix;
    const ae_int16x8* restrict pIn;
    const ae_int16x8* restrict pT;
    const ae_int16x8* restrict pIn0;
    const ae_int16x8* restrict pIn1;
    const ae_int16x8* restrict pIn2;
    const ae_int16x8* restrict pIn3;
    ae_int16x8* restrict pOut;
    ae_int16x8* restrict pOut0;
    ae_int16x8* restrict pOut1;
    ae_int16x8* restrict pOut2;
    ae_int16x8* restrict pOut3;

    ae_valignx2 al_in, al_out;
    ae_valignx2 al_in0, al_in1, al_in2, al_in3;
    ae_valignx2 al_out0, al_out1, al_out2, al_out3;

    ae_int16x4 x0, x1, x2, x3;
    ae_int16x4 x0_, x1_, x2_, x3_;
    ae_int16x4 t0, t1, t2, t3;
    ae_int16x4 sel;
    ae_int16x4 idx0, idx1, ofs;
    int m, n,
        win = in->width,
        wout = out->width,
        hout = out->height,
        ostride = out->stride,
        istride = in->stride;
    

    (void)win;
    imgsize_validate(in, 2, 0);
    imgsize_validate(out, 2, 0);
    NASSERT_ALIGN(imgIn, 2);
    NASSERT_ALIGN(imgOut, 2);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->height == out->height);

#if ! (defined(AE_LAV16X4X2_XP)&defined(AE_SAV8X8X2_XP))
    int16_t* restrict pin;
    int16_t* restrict pout;
#endif
    
    al_out = AE_ZALIGN128();

    if (win == wout * 2)
    {
#if defined(AE_LAV16X4X2_XP)&defined(AE_SAV8X8X2_XP)
        int residue0, residue1;
        residue0 = (win & 15) > 7 ? 8 : win & 15;
        residue1 = ((win & 15) > 7) ? ((win & 15) - 8) : 0;
#endif 
        /* 2x decimation */
        for (n = 0; n < hout; n++)
        {
            pIn = (const ae_int16x8 *)((int16_t *)imgIn + n*istride);
            pOut = (ae_int16x8 *)((int16_t *)imgOut + n*ostride);
            al_in = AE_LA128_PP(pIn);
            
            for (m = 0; m < (win >> 4); m++)
            {
                AE_LA16X4X2_IP(x0, x1, al_in, pIn);
                AE_LA16X4X2_IP(x2, x3, al_in, pIn);
                t0 = AE_SEL16_7531(x0, x1);
                t1 = AE_SEL16_7531(x2, x3);
                AE_SA16X4X2_IP(t0, t1, al_out, pOut);
            }
            
            if (win & 15)
            {
#if defined(AE_LAV16X4X2_XP)&defined(AE_SAV8X8X2_XP)
                AE_LAV16X4X2_XP(x0, x1, al_in, pIn, residue0 * sizeof(int16_t));
                AE_LAV16X4X2_XP(x2, x3, al_in, pIn, residue1 * sizeof(int16_t));
                t0 = AE_SEL16_7531(x0, x1);
                t1 = AE_SEL16_7531(x2, x3);
                AE_SAV16X4X2_XP(t0, t1, al_out, pOut, (wout & 7) * sizeof(int16_t));
#else
                pin = (int16_t*)pIn;
                pout = (int16_t*)pOut;
                for (m = 0; m < (wout & 7); m++)
                {
                    pout[m] = pin[2 * m];
                }
#endif
            }
            AE_SA128POS_FP(al_out, pOut);
        }
        return;
    }
    else if (win * 2 == wout)
    {

        /* 2x interpolation */
#if defined(AE_LAV16X4X2_XP)&defined(AE_SAV8X8X2_XP)
        int residue0, residue1;
        residue0 = (wout & 15) > 7 ? 8 : wout & 15;
        residue1 = ((wout & 15) > 7) ? ((wout & 15) - 8) : 0;
#endif
        sel = AE_MOVINT16X4_FROMINT64(0x0705070506040604);
        for (n = 0; n < hout; n++)
        {
            pIn = (const ae_int16x8*)((int16_t*)imgIn + n * istride);
            pOut = (ae_int16x8*)((int16_t*)imgOut + n * ostride);
            al_in = AE_LA128_PP(pIn);
            for (m = 0; m < (wout >> 4); m++)
            {
                AE_LA16X4X2_IP(x0, x1, al_in, pIn);
                
                AE_DSEL16X4(t0,t1,x0,x0,sel);
                AE_DSEL16X4(t2,t3,x1,x1,sel);

                AE_SA16X4X2_IP(t0, t1, al_out, pOut);
                AE_SA16X4X2_IP(t2, t3, al_out, pOut);
            }

            if (win & 7)
            {
#if defined(AE_LAV16X4X2_XP)&defined(AE_SAV8X8X2_XP)
                AE_LAV16X4X2_XP(x0, x1, al_in, pIn, (win & 7) * sizeof(int16_t));

                AE_DSEL16X4(t0, t1, x0, x0, sel);
                AE_DSEL16X4(t2, t3, x1, x1, sel);

                AE_SAV16X4X2_XP(t0, t1, al_out, pOut, residue0 * sizeof(int16_t));
                AE_SAV16X4X2_XP(t2, t3, al_out, pOut, residue1 * sizeof(int16_t));
#else
                pin = (int16_t*)pIn;
                pout = (int16_t*)pOut;
                for (m = 0; m < (win & 7); m++)
                {
                    pout[2*m]   = pin[m];
                    pout[2*m+1] = pin[m];
                }
#endif
            }
            AE_SA128POS_FP(al_out, pOut);
        }
        return;
    }
    al_out0 = AE_ZALIGN128();
    al_out1 = AE_ZALIGN128();
    al_out2 = AE_ZALIGN128();
    al_out3 = AE_ZALIGN128();
    sel = AE_MOVINT16X4_FROMINT64(0x0706050403020100); // SEL7531 + SEL6420
    /* Process the image by 4 rows per iteration */
    for (n = 0; n < (hout & ~3); n += 4)
    {
        /* Interleave input samples from 4 rows and save them to the scratch */
        pIn0 = (const ae_int16x8 *)((int16_t *)imgIn + (n+0)*istride);
        pIn1 = (const ae_int16x8 *)((int16_t *)imgIn + (n+1)*istride);
        pIn2 = (const ae_int16x8 *)((int16_t *)imgIn + (n+2)*istride);
        pIn3 = (const ae_int16x8 *)((int16_t *)imgIn + (n+3)*istride);
        al_in0 = AE_LA128_PP(pIn0);
        al_in1 = AE_LA128_PP(pIn1);
        al_in2 = AE_LA128_PP(pIn2);
        al_in3 = AE_LA128_PP(pIn3);
        pOut = (ae_int16x8 *)(t);
        __Pragma("loop_count min=1");
        for (m = 0; m < win; m += 8)
        {
            AE_LA16X4X2_IP(x0, x0_, al_in0, pIn0);
            AE_LA16X4X2_IP(x1, x1_, al_in1, pIn1);
            AE_LA16X4X2_IP(x2, x2_, al_in2, pIn2);
            AE_LA16X4X2_IP(x3, x3_, al_in3, pIn3);

            AE_DSEL16X4(t0, t1, x0, x1, sel);
            AE_DSEL16X4(t2, t3, x2, x3, sel);

            AE_DSEL16X4(x0, x2, t0, t2, sel);
            AE_DSEL16X4(x1, x3, t1, t3, sel);

            AE_S16X4X2_IP(x0, x1, pOut, sizeof(ae_int16x8));
            AE_S16X4X2_IP(x2, x3, pOut, sizeof(ae_int16x8));

            AE_DSEL16X4(t0, t1, x0_, x1_, sel);
            AE_DSEL16X4(t2, t3, x2_, x3_, sel);

            AE_DSEL16X4(x0_, x2_, t0, t2, sel);
            AE_DSEL16X4(x1_, x3_, t1, t3, sel);

            AE_S16X4X2_IP(x0_, x1_, pOut, sizeof(ae_int16x8));
            AE_S16X4X2_IP(x2_, x3_, pOut, sizeof(ae_int16x8));

        }

        pT = (const ae_int16x8 *)(t);
        pOut0 = (ae_int16x8 *)((int16_t *)imgOut + (n+0)*ostride);
        pOut1 = (ae_int16x8 *)((int16_t *)imgOut + (n+1)*ostride);
        pOut2 = (ae_int16x8 *)((int16_t *)imgOut + (n+2)*ostride);
        pOut3 = (ae_int16x8 *)((int16_t *)imgOut + (n+3)*ostride);
        ix = (const int16_t*)pCoef;
        for (m = 0; m < (wout >> 3); m++)
        {

            /* Load 4x4 elems for the 4 columns */
            AE_L16X4X2_IP(idx0, idx1, castxcc(ae_int16x8,ix), sizeof(ae_int16x8));
            ofs = AE_SLAI16S(idx0, 3); /* offset = idx*sizeof(ae_int16x4) */
            x0 = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_3(ofs));
            x1 = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_2(ofs));
            x2 = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_1(ofs));
            x3 = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_0(ofs));

            ofs = AE_SLAI16S(idx1, 3); /* offset = idx*sizeof(ae_int16x4) */
            x0_ = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_3(ofs));
            x1_ = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_2(ofs));
            x2_ = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_1(ofs));
            x3_ = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_0(ofs));


            AE_DSEL16X4(t0, t1, x0, x1, sel);
            AE_DSEL16X4(t2, t3, x2, x3, sel);

            AE_DSEL16X4(x0, x2, t0, t2, sel);
            AE_DSEL16X4(x1, x3, t1, t3, sel);

            AE_DSEL16X4(t0, t1, x0_, x1_, sel);
            AE_DSEL16X4(t2, t3, x2_, x3_, sel);

            AE_DSEL16X4(x0_, x2_, t0, t2, sel);
            AE_DSEL16X4(x1_, x3_, t1, t3, sel);

            AE_SA16X4X2_IP(x0, x0_, al_out0, pOut0);
            AE_SA16X4X2_IP(x1, x1_, al_out1, pOut1);
            AE_SA16X4X2_IP(x2, x2_, al_out2, pOut2);
            AE_SA16X4X2_IP(x3, x3_, al_out3, pOut3);
        }

        if (wout & 7)
        {
#if defined(AE_LAV16X4X2_XP)&defined(AE_SAV8X8X2_XP)
            //al_id = AE_LA128_PP(pIdx);
            /* Load 4x4 elems for the 4 columns */
            AE_L16X4X2_IP(idx0, idx1, castxcc(ae_int16x8, ix), sizeof(ae_int16x8));
            ofs = AE_SLAI16S(idx0, 3); /* offset = idx*sizeof(ae_int16x4) */
            x0 = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_3(ofs));
            x1 = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_2(ofs));
            x2 = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_1(ofs));
            x3 = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_0(ofs));

            ofs = AE_SLAI16S(idx1, 3); /* offset = idx*sizeof(ae_int16x4) */
            x0_ = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_3(ofs));
            x1_ = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_2(ofs));
            x2_ = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_1(ofs));
            x3_ = AE_L16X4_X((const ae_int16x4*)pT, AE_MOVAD16_0(ofs));


            AE_DSEL16X4(t0, t1, x0, x1, sel);
            AE_DSEL16X4(t2, t3, x2, x3, sel);

            AE_DSEL16X4(x0, x2, t0, t2, sel);
            AE_DSEL16X4(x1, x3, t1, t3, sel);

            AE_DSEL16X4(t0, t1, x0_, x1_, sel);
            AE_DSEL16X4(t2, t3, x2_, x3_, sel);

            AE_DSEL16X4(x0_, x2_, t0, t2, sel);
            AE_DSEL16X4(x1_, x3_, t1, t3, sel);


            AE_SAV16X4X2_XP(x0, x0_, al_out0, pOut0, (wout&7)*sizeof(int16_t));
            AE_SAV16X4X2_XP(x1, x1_, al_out1, pOut1, (wout&7)*sizeof(int16_t));
            AE_SAV16X4X2_XP(x2, x2_, al_out2, pOut2, (wout&7)*sizeof(int16_t));
            AE_SAV16X4X2_XP(x3, x3_, al_out3, pOut3, (wout&7)*sizeof(int16_t));
#else 
            pin = (int16_t*)imgIn + (n) * istride;
            pout = (int16_t*)imgOut + (n) * ostride;
            for (m = (wout & (~7)); m < wout; m++)
            {
                pout[0*ostride+m] = pin[0*istride+*ix];
                pout[1*ostride+m] = pin[1*istride+*ix];
                pout[2*ostride+m] = pin[2*istride+*ix];
                pout[3*ostride+m] = pin[3*istride+*ix++];
            }
            
#endif
        }
        AE_SA128POS_FP(al_out0, pOut0);
        AE_SA128POS_FP(al_out1, pOut1);
        AE_SA128POS_FP(al_out2, pOut2);
        AE_SA128POS_FP(al_out3, pOut3);
    }
    /* Process last 0...3 rows */
    for (n = (hout & ~3); n < hout; n++)
    {
        pIn = (const ae_int16x8 *)((int16_t *)imgIn + n*istride);
        al_in = AE_LA128_PP(pIn);
        pOut = (ae_int16x8 *)(t);

        __Pragma("loop_count min=1");
        for (m = 0; m < win; m += 8)
        {
            AE_LA16X4X2_IP(x0, x1, al_in, pIn);
            AE_S16X4X2_IP(x0, x1, pOut, sizeof(ae_int16x8));
        }

        pT = (const ae_int16x8 *)(t);
        pOut = (ae_int16x8 *)((int16_t *)imgOut + n*ostride);
        al_out = AE_ZALIGN128();
        ix = (const int16_t*)pCoef;
        for (m = 0; m < (wout >> 3); m++)
        {
            /* Load 1x4 elems for the 4 columns */
            AE_L16X4X2_IP(idx0, idx1, castxcc(ae_int16x8, ix), sizeof(ae_int16x8));

            ofs = AE_SLAI16S(idx0, 1); /* offset = idx*sizeof(int16_t) */
            x0 = AE_L16_X((const ae_int16 *)pT, AE_MOVAD16_3(ofs));
            x1 = AE_L16_X((const ae_int16 *)pT, AE_MOVAD16_2(ofs));
            x2 = AE_L16_X((const ae_int16 *)pT, AE_MOVAD16_1(ofs));
            x3 = AE_L16_X((const ae_int16 *)pT, AE_MOVAD16_0(ofs));

            ofs = AE_SLAI16S(idx1, 1); /* offset = idx*sizeof(int16_t) */
            x0_ = AE_L16_X((const ae_int16*)pT, AE_MOVAD16_3(ofs));
            x1_ = AE_L16_X((const ae_int16*)pT, AE_MOVAD16_2(ofs));
            x2_ = AE_L16_X((const ae_int16*)pT, AE_MOVAD16_1(ofs));
            x3_ = AE_L16_X((const ae_int16*)pT, AE_MOVAD16_0(ofs));

            x0 = AE_SEL16_5140(x0, x1);
            x2 = AE_SEL16_5140(x2, x3);
            x0 = AE_SEL16_5410(x0, x2);

            x0_ = AE_SEL16_5140(x0_, x1_);
            x2_ = AE_SEL16_5140(x2_, x3_);
            x0_ = AE_SEL16_5410(x0_, x2_);

            AE_SA16X4X2_IP(x0, x0_, al_out, pOut);
        }
        AE_SA128POS_FP(al_out, pOut);
        if (wout & 7)
        {
#if defined(AE_LAV16X4X2_XP)&defined(AE_SAV8X8X2_XP)
            ae_valignx2 al_id = AE_LA128_PP(castxcc(ae_int16x8, ix));

            AE_LAV16X4X2_XP(idx0, idx1, al_id, castxcc(ae_int16x8, ix), (wout & 7) * sizeof(int16_t));

            ofs = AE_SLAI16S(idx0, 1); /* offset = idx*sizeof(int16_t) */
            x0 = AE_L16_X((const ae_int16*)pT, AE_MOVAD16_3(ofs));
            x1 = AE_L16_X((const ae_int16*)pT, AE_MOVAD16_2(ofs));
            x2 = AE_L16_X((const ae_int16*)pT, AE_MOVAD16_1(ofs));
            x3 = AE_L16_X((const ae_int16*)pT, AE_MOVAD16_0(ofs));

            ofs = AE_SLAI16S(idx1, 1); /* offset = idx*sizeof(int16_t) */
            x0_ = AE_L16_X((const ae_int16*)pT, AE_MOVAD16_3(ofs));
            x1_ = AE_L16_X((const ae_int16*)pT, AE_MOVAD16_2(ofs));
            x2_ = AE_L16_X((const ae_int16*)pT, AE_MOVAD16_1(ofs));
            x3_ = AE_L16_X((const ae_int16*)pT, AE_MOVAD16_0(ofs));

            x0 = AE_SEL16_5140(x0, x1);
            x2 = AE_SEL16_5140(x2, x3);
            x0 = AE_SEL16_5410(x0, x2);

            x0_ = AE_SEL16_5140(x0_, x1_);
            x2_ = AE_SEL16_5140(x2_, x3_);
            x0_ = AE_SEL16_5410(x0_, x2_);

            AE_SAV16X4X2_XP(x0, x0_, al_out, pOut, (wout & 7) * sizeof(int16_t));
            AE_SA128POS_FP(al_out, pOut);
#else
            pin = (int16_t*)imgIn + n * istride;
            pout = ((int16_t*)imgOut) + n *ostride;
            for (m = (wout & (~7)); m < wout; m++)
            {
                pout[m] = pin[*ix++];
            }
#endif
        }
        
    }
}


/* 8 bit fast process */
static void process_fast_gs(void *pScr, void* pCoef, const void* imgIn, void*imgOut, const imgsize_t* in, const imgsize_t* out)
{
    const uint8_t* restrict pin;
    uint8_t* restrict pout;
    int16_t* restrict t = (int16_t*)pScr;
    const int16_t* restrict ix = (const int16_t*)pCoef;
    const ae_int8x16* restrict pIn;
    ae_int8x16* restrict pOut;
    const ae_int16x8* restrict pIdx;

    ae_valignx2 al_out;

    ae_int8x8 x00, x10, x20, x30, x40, x50, x60, x70;
    ae_int8x8 x01, x11, x21, x31, x41, x51, x61, x71;

    ae_int8x8 x0, x1, x2, x3;
    ae_int8x8 t0, t1, t2, t3, t4, t5, t6, t7, sel, interliver;
    ae_int8x8 t0_1, t1_1, t2_1, t3_1, t4_1, t5_1, t6_1, t7_1;
    int m, n,
        win = in->width,
        wout = out->width,
        hout = out->height,
        ostride = out->stride,
        istride = in->stride;


    (void)win;
    imgsize_validate(in, 1, 1);
    imgsize_validate(out, 1, 1);
    NASSERT_ALIGN(imgIn, ALIGNMENT);
    NASSERT_ALIGN(imgOut, ALIGNMENT);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->height == out->height);

    al_out = AE_ZALIGN128();

    if (win == wout * 2)
    {

        /* 2x decimation */
#if defined(AE_LAV16X4X2_XP)&defined(AE_SAV8X8X2_XP)
        int residue0, residue1;
        residue0 = (win & 31) > 15 ? 16 : win & 31;
        residue1 = ((win & 31) > 15) ? ((win & 31) - 16) : 0;
#endif
        for (n = 0; n < hout; n++)
        {
            pIn = (const ae_int8x16*)((uint8_t*)imgIn + n * istride);
            pOut = (ae_int8x16*)((uint8_t*)imgOut + n * ostride);

            for (m = 0; m < (wout >> 4); m++)
            {
                AE_L8X8X2_IP(x0, x1, pIn, sizeof(ae_int8x16));
                AE_L8X8X2_IP(x2, x3, pIn, sizeof(ae_int8x16));

                t0 = AE_SEL8X8I(x0, x1, 25);
                t1 = AE_SEL8X8I(x2, x3, 25);

                AE_S8X8X2_IP(t0, t1, pOut, sizeof(ae_int8x16));
            }
            if (wout & 15)
            {
#if defined(AE_LAV16X4X2_XP)&defined(AE_SAV8X8X2_XP)
                ae_valignx2 al_in, al_out;
                al_in = AE_LA128_PP(pIn);
                al_out = AE_ZALIGN128();
                AE_LAV8X8X2_XP(x0, x1, al_in, pIn, residue0 * sizeof(int8_t));
                AE_LAV8X8X2_XP(x2, x3, al_in, pIn, residue1 * sizeof(int8_t));

                t0 = AE_SEL8X8I(x0, x1, 25);
                t1 = AE_SEL8X8I(x2, x3, 25);

                AE_SAV8X8X2_XP(t0, t1, al_out, pOut, (wout & 15) * sizeof(int8_t));
                AE_SA128POS_FP(al_out, pOut);
#else
                uint8_t* pin = (uint8_t*)pIn;
                uint8_t* pout = (uint8_t*)pOut;
                for (m = 0; m < (wout & 15); m++)
                {
                    pout[m] = pin[2 * m];
                }
#endif
            }

        }
        return;
    }
    else if (win * 2 == wout)
    {
        sel = AE_MOVINT8X8_FROMINT64(0x7373626251514040);
#if defined(AE_LAV16X4X2_XP)&defined(AE_SAV8X8X2_XP)
        int residue0, residue1;
        residue0 = (wout & 31) > 15 ? 16 : wout & 31;
        residue1 = ((wout & 31) > 15) ? ((wout & 31) - 16) : 0;
#endif
        /* 2x interpolation */
        for (n = 0; n < hout; n++)
        {
            pIn = (const ae_int8x16*)((uint8_t*)imgIn + n * istride);
            pOut = (ae_int8x16*)((uint8_t*)imgOut + n * ostride);            

            for (m = 0; m < (win >> 4); m++)
            {
                AE_L8X8X2_IP(x0, x1, pIn, sizeof(ae_int8x16));

                AE_DSEL8X8(t0, t1, x0, x0, sel);
                AE_DSEL8X8(t2, t3, x1, x1, sel);

                AE_S8X8X2_IP(t0, t1, pOut, sizeof(ae_int8x16));
                AE_S8X8X2_IP(t2, t3, pOut, sizeof(ae_int8x16));
            }
            if (win & 15)
            {
#if defined(AE_LAV16X4X2_XP)&defined(AE_SAV8X8X2_XP)
                ae_valignx2 al_in, al_out;
                al_in = AE_LA128_PP(pIn);
                al_out = AE_ZALIGN128();
                AE_LAV8X8X2_XP(x0, x1, al_in, pIn, (win & 15) * sizeof(int8_t));

                AE_DSEL8X8(t0, t1, x0, x0, sel);
                AE_DSEL8X8(t2, t3, x1, x1, sel);

                AE_SAV8X8X2_XP(t0, t1, al_out, pOut, residue0 * sizeof(int8_t));
                AE_SAV8X8X2_XP(t2, t3, al_out, pOut, residue1 * sizeof(int8_t));
                AE_SA128POS_FP(al_out, pOut);
#else
                uint8_t* pin = (uint8_t*)pIn;
                uint8_t* pout = (uint8_t*)pOut;
                for (m = 0; m < (win & 15); m++)
                {
                    pout[2*m  ] = pin[m];
                    pout[2*m+1] = pin[m];
                }

#endif
            }


        }
        return;
    }

    interliver = AE_MOVINT8X8_FROMINT64(0xFEDCBA9876543210);

    /* Process the image by 8 rows per iteration */
    for (n = 0; n < (hout & ~7); n += 8)
    {

        /* Interleave input samples from 8 rows and save them to the scratch */
        pIn = (const ae_int8x16*)((uint8_t*)imgIn + n * istride);

        pOut = (ae_int8x16*)t;

        for (m = 0; m < win; m += 16)
        {
            AE_L8X8X2_XP(x00, x01, pIn, istride);
            AE_L8X8X2_XP(x10, x11, pIn, istride);
            AE_L8X8X2_XP(x20, x21, pIn, istride);
            AE_L8X8X2_XP(x30, x31, pIn, istride);
            AE_L8X8X2_XP(x40, x41, pIn, istride);
            AE_L8X8X2_XP(x50, x51, pIn, istride);
            AE_L8X8X2_XP(x60, x61, pIn, istride);
            AE_L8X8X2_XP(x70, x71, pIn, -7*istride + sizeof(ae_int8x16));


            AE_DSEL8X8(t0, t1, x00, x10, interliver);
            AE_DSEL8X8(t2, t3, x20, x30, interliver);
            AE_DSEL8X8(t4, t5, x40, x50, interliver);
            AE_DSEL8X8(t6, t7, x60, x70, interliver);

            AE_DSEL8X8(t0_1, t1_1, t0, t2, interliver);
            AE_DSEL8X8(t2_1, t3_1, t1, t3, interliver);
            AE_DSEL8X8(t4_1, t5_1, t4, t6, interliver);
            AE_DSEL8X8(t6_1, t7_1, t5, t7, interliver);
                                         
            AE_DSEL8X8(x00, x40, t0_1, t4_1, interliver);
            AE_DSEL8X8(x20, x60, t1_1, t5_1, interliver);
            AE_DSEL8X8(x10, x50, t2_1, t6_1, interliver);
            AE_DSEL8X8(x30, x70, t3_1, t7_1, interliver);

            AE_DSEL8X8(t0, t1, x01, x11, interliver);
            AE_DSEL8X8(t2, t3, x21, x31, interliver);
            AE_DSEL8X8(t4, t5, x41, x51, interliver);
            AE_DSEL8X8(t6, t7, x61, x71, interliver);

            AE_DSEL8X8(t0_1, t1_1, t0, t2, interliver);
            AE_DSEL8X8(t2_1, t3_1, t1, t3, interliver);
            AE_DSEL8X8(t4_1, t5_1, t4, t6, interliver);
            AE_DSEL8X8(t6_1, t7_1, t5, t7, interliver);

            AE_DSEL8X8(x01, x41, t0_1, t4_1, interliver);
            AE_DSEL8X8(x21, x61, t1_1, t5_1, interliver);
            AE_DSEL8X8(x11, x51, t2_1, t6_1, interliver);
            AE_DSEL8X8(x31, x71, t3_1, t7_1, interliver);

            AE_S8X8X2_IP(x00, x10, pOut, sizeof(ae_int8x16));
            AE_S8X8X2_IP(x20, x30, pOut, sizeof(ae_int8x16));
            AE_S8X8X2_IP(x40, x50, pOut, sizeof(ae_int8x16));
            AE_S8X8X2_IP(x60, x70, pOut, sizeof(ae_int8x16));
            AE_S8X8X2_IP(x01, x11, pOut, sizeof(ae_int8x16));
            AE_S8X8X2_IP(x21, x31, pOut, sizeof(ae_int8x16));
            AE_S8X8X2_IP(x41, x51, pOut, sizeof(ae_int8x16));
            AE_S8X8X2_IP(x61, x71, pOut, sizeof(ae_int8x16));

        }

        pIdx = (const ae_int16x8*)(ix);
        pOut = (ae_int8x16*)((uint8_t*)imgOut + (n * ostride));
        pIn = (const ae_int8x16*)t;

        for (m = 0; m < (wout >> 4); m++)
        {
            ae_int16x4 idx0, idx1, ofs;

            AE_L16X4X2_IP(idx0,idx1, pIdx, sizeof(ae_int16x8));
            ofs = AE_SLAI16S(idx0, 3);
            x00 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_3(ofs));
            x10 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_2(ofs));
            x20 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_1(ofs));
            x30 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_0(ofs));

            ofs = AE_SLAI16S(idx1, 3);
            x40 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_3(ofs));
            x50 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_2(ofs));
            x60 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_1(ofs));
            x70 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_0(ofs));


            AE_L16X4X2_IP(idx0, idx1, pIdx, sizeof(ae_int16x8));
            ofs = AE_SLAI16S(idx0, 3);
            x01 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_3(ofs));
            x11 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_2(ofs));
            x21 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_1(ofs));
            x31 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_0(ofs));

            ofs = AE_SLAI16S(idx1, 3);
            x41 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_3(ofs));
            x51 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_2(ofs));
            x61 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_1(ofs));
            x71 = AE_L8X8_X((const ae_int8x8*)pIn, AE_MOVAD16_0(ofs));

            AE_DSEL8X8(t0, t1, x00, x10, interliver);
            AE_DSEL8X8(t2, t3, x20, x30, interliver);
            AE_DSEL8X8(t4, t5, x40, x50, interliver);
            AE_DSEL8X8(t6, t7, x60, x70, interliver);

            AE_DSEL8X8(t0_1, t1_1, t0, t2, interliver);
            AE_DSEL8X8(t2_1, t3_1, t1, t3, interliver);
            AE_DSEL8X8(t4_1, t5_1, t4, t6, interliver);
            AE_DSEL8X8(t6_1, t7_1, t5, t7, interliver);

            AE_DSEL8X8(x00, x40, t0_1, t4_1, interliver);
            AE_DSEL8X8(x20, x60, t1_1, t5_1, interliver);
            AE_DSEL8X8(x10, x50, t2_1, t6_1, interliver);
            AE_DSEL8X8(x30, x70, t3_1, t7_1, interliver);

            AE_DSEL8X8(t0, t1, x01, x11, interliver);
            AE_DSEL8X8(t2, t3, x21, x31, interliver);
            AE_DSEL8X8(t4, t5, x41, x51, interliver);
            AE_DSEL8X8(t6, t7, x61, x71, interliver);

            AE_DSEL8X8(t0_1, t1_1, t0, t2, interliver);
            AE_DSEL8X8(t2_1, t3_1, t1, t3, interliver);
            AE_DSEL8X8(t4_1, t5_1, t4, t6, interliver);
            AE_DSEL8X8(t6_1, t7_1, t5, t7, interliver);

            AE_DSEL8X8(x01, x41, t0_1, t4_1, interliver);
            AE_DSEL8X8(x21, x61, t1_1, t5_1, interliver);
            AE_DSEL8X8(x11, x51, t2_1, t6_1, interliver);
            AE_DSEL8X8(x31, x71, t3_1, t7_1, interliver);
               
            AE_S8X8X2_XP(x00, x01, pOut, ostride);
            AE_S8X8X2_XP(x10, x11, pOut, ostride);
            AE_S8X8X2_XP(x20, x21, pOut, ostride);
            AE_S8X8X2_XP(x30, x31, pOut, ostride);
            AE_S8X8X2_XP(x40, x41, pOut, ostride);
            AE_S8X8X2_XP(x50, x51, pOut, ostride);
            AE_S8X8X2_XP(x60, x61, pOut, ostride);
            AE_S8X8X2_XP(x70, x71, pOut, -7*ostride+sizeof(ae_int8x16));
        }

        if (wout & 15)
        {
            int i;
            for (i = 0; i < 8; i++)
            {
                pin = (uint8_t*)imgIn + (n + i) * istride;
                pout = ((uint8_t*)imgOut) + (n + i) * ostride;
                for (m = (wout & (~15)); m < wout; m++)
                {
                    pout[m] = pin[ix[m]];
                }
            }
        }
    }
    /* Process last 0...7 rows */
    for (n = (hout & ~7); n < hout; n++)
    {
        pin = ((uint8_t*)imgIn) + n * istride;
        pout = ((uint8_t*)imgOut) + n * ostride;
        for (m = 0; m < wout; m++)
        {
            pout[m] = pin[ix[m]];
        }
    }
}

/* 8 bit process */
static void process_gs(void *pScr, void* pCoef, const void* imgIn, void*imgOut, const imgsize_t* in, const imgsize_t* out)
{
    const uint8_t * restrict pin;
    uint8_t * restrict pout;
    int16_t * restrict t = (int16_t *)pScr;
    const int16_t * restrict ix = (const int16_t *)pCoef;
    const ae_int8x16 * restrict pIn;
    const ae_int8x8 * restrict pInx8;
    const ae_int8x16 * restrict pIn0;
    const ae_int8x16 * restrict pIn1;
    const ae_int8x16 * restrict pIn2;
    const ae_int8x16 * restrict pIn3;
    ae_valignx2 al_in, al_out;
    ae_valignx2 al_in0, al_in1, al_in2, al_in3;
    ae_valign al_ind, al_out0, al_out1, al_out2, al_out3;


    ae_int8x16 * restrict pOut;
    ae_int8x8 * restrict pOut0;
    ae_int8x8 * restrict pOut1;
    ae_int8x8 * restrict pOut2;
    ae_int8x8 * restrict pOut3;
    const ae_int16x4 * restrict pIdx;

    ae_int8x8 x00, x10, x20, x30;
    ae_int8x8 x01, x11, x21, x31;

    ae_int8x8 x0, x1, x2, x3;
    ae_int8x8 x0_, x1_, x2_, x3_;
    ae_int8x8 t0, t1, t2, t3, sel, selHL;
    ae_int8x8 sel0, sel1, sel2, sel3;
    ae_int8x8 t0_0, t1_0, t2_0, t3_0;
    ae_int8x8 t0_1, t1_1, t2_1, t3_1;
    int m, n,
        win = in->width,
        wout = out->width,
        hout = out->height,
        ostride = out->stride,
        istride = in->stride;

    (void)win;
    imgsize_validate(in, 1, 0);
    imgsize_validate(out, 1, 0);
    NASSERT_ALIGN(pScr, ALIGNMENT);
    NASSERT(in->height == out->height);

    if (win == wout * 2)
    {
#if defined(AE_LAV16X4X2_XP)& defined(AE_SAV8X8X2_XP)
        int residue0, residue1;
        residue0 = (win & 31) > 15 ? 16 : win & 31;
        residue1 = ((win&31)>15)?((win & 31) - 16):0;
#endif
        /* 2x decimation */
        for (n = 0; n < hout; n++)
        {
            pIn = (const ae_int8x16 *)((uint8_t *)imgIn + n*istride);
            pOut = (ae_int8x16 *)((uint8_t *)imgOut + n*ostride);

            al_in = AE_LA128_PP(pIn);
            al_out = AE_ZALIGN128();
            for(m = 0; m < (wout >> 4); m++)
            {
                AE_LA8X8X2_IP(x0, x1, al_in, pIn);
                AE_LA8X8X2_IP(x2, x3, al_in, pIn);

                t0 = AE_SEL8X8I(x0, x1, 25);
                t1 = AE_SEL8X8I(x2, x3, 25);

                AE_SA8X8X2_IP(t0, t1, al_out, pOut);
            }
            if (wout &15)
            {
#if defined(AE_LAV16X4X2_XP)& defined(AE_SAV8X8X2_XP)
                AE_LAV8X8X2_XP(x0, x1, al_in, pIn, residue0  * sizeof(int8_t));
                AE_LAV8X8X2_XP(x2, x3, al_in, pIn, residue1  * sizeof(int8_t));

                t0 = AE_SEL8X8I(x0, x1, 25);
                t1 = AE_SEL8X8I(x2, x3, 25);

                AE_SAV8X8X2_XP(t0, t1, al_out, pOut, (wout & 15) * sizeof(int8_t));

#else

                uint8_t* pin = (uint8_t*)pIn;
                uint8_t* pout = (uint8_t*)pOut;
                for (m = 0; m < (wout & 15); m++)
                {
                    pout[m] = pin[2 * m];
                }

#endif
            }
            AE_SA128POS_FP(al_out, pOut);
        }
        return;
    }
    else if (win * 2 == wout)
    {
        sel = AE_MOVINT8X8_FROMINT64(0x7373626251514040);

#if defined(AE_LAV16X4X2_XP)& defined(AE_SAV8X8X2_XP)
        int residue0, residue1;
        residue0 = (wout & 31)>15?16: wout & 31;
        residue1 = ((wout & 31) > 15) ? ((wout & 31) - 16) : 0;
#endif
        /* 2x interpolation */
        for (n = 0; n < hout; n++)
        {
            pIn = (const ae_int8x16*)((uint8_t*)imgIn + n * istride);
            pOut = (ae_int8x16*)((uint8_t*)imgOut + n * ostride);

            al_in = AE_LA128_PP(pIn);
            al_out = AE_ZALIGN128();

            for (m = 0; m < (win >> 4); m++)
            {
                AE_LA8X8X2_IP(x0, x1, al_in, pIn);

                AE_DSEL8X8(t0, t1, x0, x0, sel);
                AE_DSEL8X8(t2, t3, x1, x1, sel);

                AE_SA8X8X2_IP(t0, t1, al_out, pOut);
                AE_SA8X8X2_IP(t2, t3, al_out, pOut);
            }
            if (win & 15)
            {
#if defined(AE_LAV16X4X2_XP)& defined(AE_SAV8X8X2_XP)
                AE_LAV8X8X2_XP(x0, x1, al_in, pIn, (win & 15) * sizeof(int8_t));

                AE_DSEL8X8(t0, t1, x0, x0, sel);
                AE_DSEL8X8(t2, t3, x1, x1, sel);

                AE_SAV8X8X2_XP(t0, t1, al_out, pOut, residue0 * sizeof(int8_t));
                AE_SAV8X8X2_XP(t2, t3, al_out, pOut, residue1 * sizeof(int8_t));
#else
                uint8_t* pin = (uint8_t*)pIn;
                uint8_t* pout = (uint8_t*)pOut;
                for (m = 0; m < (win & 15); m++)
                {
                    pout[2 * m] = pin[m];
                    pout[2 * m + 1] = pin[m];
                }
#endif
            }
            AE_SA128POS_FP(al_out, pOut);

        }
        return;
    }

    //ABCDEF FEDCBA98
    sel0 = AE_MOVINT8X8_FROMINT64(0xFD750000EC640000); // DSEL x0 x1
    sel1 = AE_MOVINT8X8_FROMINT64(0x0000FD750000EC64); // DSEL x2 x3

    sel2 = AE_MOVINT8X8_FROMINT64(0xB9310000A8200000); // DSEL x0 x1 part 2
    sel3 = AE_MOVINT8X8_FROMINT64(0x0000B9310000A820); // DSEL x2 x3 part 2

    sel = AE_MOVINT8X8_FROMINT64(0x0F0E05040B0A0100); // SEL  t00 t01

    selHL = AE_MOVINT8X8_FROMINT64(0xFBEAD9C8FBEAD9C8); // DSEL HL t00 t01

    /* Process the image by 4 rows per iteration */
    for (n = 0; n < (hout & ~3); n += 4)
    {

        /* Interleave input samples from 4 rows and save them to the scratch */
        pIn0 = (const ae_int8x16 *)((uint8_t *)imgIn + (n + 0)*istride);
        pIn1 = (const ae_int8x16 *)((uint8_t *)imgIn + (n + 1)*istride);
        pIn2 = (const ae_int8x16 *)((uint8_t *)imgIn + (n + 2)*istride);
        pIn3 = (const ae_int8x16 *)((uint8_t *)imgIn + (n + 3)*istride);

        al_in0 = AE_LA128_PP(pIn0);
        al_in1 = AE_LA128_PP(pIn1);
        al_in2 = AE_LA128_PP(pIn2);
        al_in3 = AE_LA128_PP(pIn3);

        pOut  = (ae_int8x16 *)t;
 

        for (m = 0; m < win; m+=16)
        {

            AE_LA8X8X2_IP(x00, x01, al_in0, pIn0);
            AE_LA8X8X2_IP(x10, x11, al_in1, pIn1);
            AE_LA8X8X2_IP(x20, x21, al_in2, pIn2);
            AE_LA8X8X2_IP(x30, x31, al_in3, pIn3);


            AE_DSEL8X8(t0, t1, x00, x10, sel0);
            AE_DSEL8X8(t2, t3, x00, x10, sel2);
            
            AE_DSEL8X8(t0_0, t1_0, x20, x30, sel1);
            AE_DSEL8X8(t2_0, t3_0, x20, x30, sel3);

            t0 = AE_SEL8X8(t0, t0_0, sel);
            t1 = AE_SEL8X8(t1, t1_0, sel);
            t2 = AE_SEL8X8(t2, t2_0, sel);
            t3 = AE_SEL8X8(t3, t3_0, sel);

            AE_DSEL8X8(t0_0, t0_1, t0, t0, selHL);
            AE_DSEL8X8(t1_0, t1_1, t1, t1, selHL);
            AE_DSEL8X8(t2_0, t2_1, t2, t2, selHL);
            AE_DSEL8X8(t3_0, t3_1, t3, t3, selHL);

            AE_S8X8X2_IP(t0_0, t0_1, pOut, sizeof(ae_int8x16));
            AE_S8X8X2_IP(t1_0, t1_1, pOut, sizeof(ae_int8x16));
            AE_S8X8X2_IP(t2_0, t2_1, pOut, sizeof(ae_int8x16));
            AE_S8X8X2_IP(t3_0, t3_1, pOut, sizeof(ae_int8x16));

            AE_DSEL8X8(t0, t1, x01, x11, sel0);
            AE_DSEL8X8(t2, t3, x01, x11, sel2);

            AE_DSEL8X8(t0_0, t1_0, x21, x31, sel1);
            AE_DSEL8X8(t2_0, t3_0, x21, x31, sel3);

            t0 = AE_SEL8X8(t0, t0_0, sel);
            t1 = AE_SEL8X8(t1, t1_0, sel);
            t2 = AE_SEL8X8(t2, t2_0, sel);
            t3 = AE_SEL8X8(t3, t3_0, sel);

            AE_DSEL8X8(t0_0, t0_1, t0, t0, selHL);
            AE_DSEL8X8(t1_0, t1_1, t1, t1, selHL);
            AE_DSEL8X8(t2_0, t2_1, t2, t2, selHL);
            AE_DSEL8X8(t3_0, t3_1, t3, t3, selHL);

            AE_S8X8X2_IP(t0_0, t0_1, pOut, sizeof(ae_int8x16));
            AE_S8X8X2_IP(t1_0, t1_1, pOut, sizeof(ae_int8x16));
            AE_S8X8X2_IP(t2_0, t2_1, pOut, sizeof(ae_int8x16));
            AE_S8X8X2_IP(t3_0, t3_1, pOut, sizeof(ae_int8x16));
        }

        pIdx = (const ae_int16x4 *)(ix);
        al_ind = AE_LA64_PP(pIdx);
        pOut0 = (ae_int8x8 *)((uint8_t *)imgOut + (n + 0)*ostride);
        pOut1 = (ae_int8x8 *)((uint8_t *)imgOut + (n + 1)*ostride);
        pOut2 = (ae_int8x8 *)((uint8_t *)imgOut + (n + 2)*ostride);
        pOut3 = (ae_int8x8 *)((uint8_t *)imgOut + (n + 3)*ostride);
        al_out0 = AE_ZALIGN64();
        al_out1 = AE_ZALIGN64();
        al_out2 = AE_ZALIGN64();
        al_out3 = AE_ZALIGN64();
        pInx8 = (const ae_int8x8 *)t;

        for (m = 0; m < (wout >> 3); m++)
        {
            ae_int16x4 idx, ofs;

            /* Load 4x8 elems for the 8 columns */
            AE_LA16X4_IP(idx, al_ind, pIdx);
            ofs = AE_SLAI16S(idx, 3); 

            x0  = AE_L8X8_X(pInx8, AE_MOVAD16_3(ofs));
            x1  = AE_L8X8_X(pInx8, AE_MOVAD16_2(ofs));
            x2  = AE_L8X8_X(pInx8, AE_MOVAD16_1(ofs));
            x3  = AE_L8X8_X(pInx8, AE_MOVAD16_0(ofs));

            AE_LA16X4_IP(idx, al_ind, pIdx);
            ofs = AE_SLAI16S(idx, 3);
            x0_ = AE_L8X8_X(pInx8, AE_MOVAD16_3(ofs));
            x1_ = AE_L8X8_X(pInx8, AE_MOVAD16_2(ofs));
            x2_ = AE_L8X8_X(pInx8, AE_MOVAD16_1(ofs));
            x3_ = AE_L8X8_X(pInx8, AE_MOVAD16_0(ofs));


            t0_0 = AE_SEL8X8I(x0, x1,20);
            t2_0 = AE_SEL8X8I(x0_, x1_,20);
            t1_0 = AE_SEL8X8I(x2, x3, 20);
            t3_0 = AE_SEL8X8I(x2_, x3_,20);

            t0_1 = AE_SEL8X8I(t0_0, t1_0, 9);
            t1_1 = AE_SEL8X8I(t0_0, t1_0, 11);
            t2_1 = AE_SEL8X8I(t2_0, t3_0, 9);
            t3_1 = AE_SEL8X8I(t2_0, t3_0, 11);

            t0 = AE_SEL8X8I(t0_1, t2_1, 1);
            t1 = AE_SEL8X8I(t0_1, t2_1, 3);
            t2 = AE_SEL8X8I(t1_1, t3_1, 1);
            t3 = AE_SEL8X8I(t1_1, t3_1, 3);

            AE_SA8X8_IP(t0, al_out0, pOut0);
            AE_SA8X8_IP(t1, al_out1, pOut1);
            AE_SA8X8_IP(t2, al_out2, pOut2);
            AE_SA8X8_IP(t3, al_out3, pOut3);
        }

        AE_SA64POS_FP(al_out0, pOut0);
        AE_SA64POS_FP(al_out1, pOut1);
        AE_SA64POS_FP(al_out2, pOut2);
        AE_SA64POS_FP(al_out3, pOut3);

        if (wout & 7)
        {
            int i;
            for (i = 0; i < 4; i++)
            {
                pin = (uint8_t *)imgIn + (n + i)*istride;
                pout = ((uint8_t *)imgOut) + (n + i)*ostride;
                for (m = (wout&(~7)); m < wout; m++)
                {
                    pout[m] = pin[ix[m]];
                }
            }
        }
    }
    /* Process last 0...3 rows */
    for (n = (hout & ~3); n<hout; n++)
    {
        pin = ((uint8_t *)imgIn) + n*istride;
        pout = ((uint8_t *)imgOut) + n*ostride;
        for (m = 0; m<wout; m++)
        {
            pout[m] = pin[ix[m]];
        }
    }
}

const imgresizer_api_t imgresizer_api_nh={NULL,getCoefSz,getCoef,getScratchSize,process};
