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
/*
  NatureDSP Signal Processing Library. FIR filters
    Gerenalized Cross-Correlation Phase Transform (GCC-PHAT)
    Real floating point input/output data
    C code optimized for HiFi5 core with VFPU/SFPU
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "common_fpu.h"
/* GCC-PHAT common declarations. */
#include "gccphat_common.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
#include "NatureDSP_Signal_fft.h"
#include <float.h>
#include "inff_tbl.h"

#define ALIGN_SIZE     (HIFI_SIMD_WIDTH)
#define ALIGN_PAD      (ALIGN_SIZE-1)
#define ALIGN_PTR(p)   (void*)(((uintptr_t)(p)+ALIGN_PAD)&~ALIGN_PAD)
#define sz_f32         sizeof(float32_t)
#define sz_cf32        sizeof(complex_float)

/*-------------------------------------------------------------------------
  Generalized Cross-Correlation Phase Transform (GCC-PHAT)
  Estimate the generalized cross-correlation of input vectors x and y to define
  the time-shift between them:
  NOTE:
    GCC-PHAT implementation performs a few FFT transforms of size N rounded up
    to the next power of 2. That is, these functions are utilized most efficiently
    when the input argument N is already a power of 2.

  Precision: 
  32x32   real 32-bit input/output data
  f       real floating point input/output data

  Input:
  N       input/output vectors size
  x[N]    input vector x
  y[N]    input vector y
  Output:
  z[N]    GCC-PHAT results vector, Q31 for the fixed-point variant
  Returned value:
          Fixed-point variant function returns the sum of bi-directional right 
          shift amounts applied throughout the inverse FFT transform.
  Temporary:
  pScr    scratch memory area of size specified by the corresponding scratch
          allocation function (in bytes).

  Restrictions:
  x,y,z   must not overlap and must be aligned by 16-bytes
  pScr    must be aligned by 16-bytes
  N       multiple of 32, 64<=N<=320
-------------------------------------------------------------------------*/
#if (!HAVE_VFPU && !HAVE_FPU) 
DISCARD_FUN(void, gccphatf, (void* pScr, 
                             float32_t * z, 
                       const float32_t * x, 
                       const float32_t * y, 
                       int N));

/* Returns: size of scratch memory area, in bytes. */
size_t gccphatf_getScratchSize(int N)
{
    (void)N;
    return 0;
}
#elif HAVE_VFPU
void gccphatf(void* pScr, float32_t * z, const float32_t * x, const float32_t * y, int N)
{
    const xtfloatx4* restrict pX;
    const xtfloatx4* restrict pY;
          xtfloatx4* restrict pZ;

    float32_t *a_x, *a_y, *a_z; /* Time-domain representation of vectors x, y and z. */
    complex_float *a_X, *a_Y, *a_Z; /* Frequency-domain representation of vectors x, y and z. */
    const complex_float * twdTbl = (complex_float*)gccphatf_twd_tbl;
    int logNup, Nup;
    int twdStp, n;
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(z, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT(0==(N%32) && 64<=N && N<=320);
    NASSERT(N<=GCCPHATF_TWD_TBL_SZ);
    /* Round N up to the next power of 2 and compute the step size for the FFT twiddles. */
    logNup = 31-NSA(N-1); Nup = 1<<logNup;
    twdStp = 1<<(GCCPHATF_TWD_TBL_LOG_SZ-logNup);
    /* Partition the scratch memory area */
    {
        void *p = pScr;
        a_x = (float32_t    *)p; p = ALIGN_PTR(a_x + Nup);
        a_X = (complex_float*)p; p = ALIGN_PTR(a_X + Nup/2+1);
        a_Y = (complex_float*)p; p = ALIGN_PTR(a_Y + Nup/2+1);
#ifdef _DEBUG
        NASSERT((int8_t*)p-(int8_t*)pScr <= (int)gccphatf_getScratchSize(N));
#endif
        /* Partially re-use the allocated arrays. */
        a_y = a_x; a_Z = a_X;
        a_z = N==Nup ? z : a_x;
    }
    /* Copy data from the input argument x to scratch array, to feel free for
     * the FFT routine to distort them. Then, pad the FFT window with zeros. */
    NASSERT(N%32==0 && Nup%32==0);
    pX=(xtfloatx4*)x;
    pZ=(xtfloatx4*)a_x;
    __Pragma("loop_count min=4")
    for ( n=0; n<(N>>3); n++ ) 
    {
        xtfloatx2 t0,t1,t2,t3;
        AE_LSX2X2_IP(t0,t1,pX,sizeof(xtfloatx4));
        AE_LSX2X2_IP(t2,t3,pX,sizeof(xtfloatx4));
        AE_SSX2X2_IP(t0,t1,pZ,sizeof(xtfloatx4));
        AE_SSX2X2_IP(t2,t3,pZ,sizeof(xtfloatx4));
    }
    for ( ; n<(Nup>>3); n++ ) 
    {
        AE_SSX2X2_IP(XT_CONST_S(0),XT_CONST_S(0),pZ,sizeof(xtfloatx4));
        AE_SSX2X2_IP(XT_CONST_S(0),XT_CONST_S(0),pZ,sizeof(xtfloatx4));
    }
    /* X <- fft(x) */
    fft_realf_ie(a_X, a_x, twdTbl, twdStp, Nup);
    /* Copy data from the argument y to the FFT window and pad it with zeros. */
    NASSERT(N%32==0 && Nup%32==0);
    pX=(xtfloatx4*)y;
    pZ=(xtfloatx4*)a_y;
    __Pragma("loop_count min=4")
    for ( n=0; n<(N>>3); n++ ) 
    {
        xtfloatx2 t0,t1,t2,t3;
        AE_LSX2X2_IP(t0,t1,pX,sizeof(xtfloatx4));
        AE_LSX2X2_IP(t2,t3,pX,sizeof(xtfloatx4));
        AE_SSX2X2_IP(t0,t1,pZ,sizeof(xtfloatx4));
        AE_SSX2X2_IP(t2,t3,pZ,sizeof(xtfloatx4));
    }
    for ( ; n<(Nup>>3); n++ ) 
    {
        AE_SSX2X2_IP(XT_CONST_S(0),XT_CONST_S(0),pZ,sizeof(xtfloatx4));
        AE_SSX2X2_IP(XT_CONST_S(0),XT_CONST_S(0),pZ,sizeof(xtfloatx4));
    }
    /* Y <- fft(y); */
    fft_realf_ie(a_Y, a_y, twdTbl, twdStp, Nup);
    /* Z <- X.*conj(Y) */
    /* Z <- exp(1j*angle(Z)) */
    pX=(const xtfloatx4*)a_X;
    pY=(const xtfloatx4*)a_Y;
    pZ=(      xtfloatx4*)a_Z;
#if 0
    __Pragma("loop_count min=4,factor=4")
    for (n=0; n<(Nup>>2); n++)
    {
        int aset1;
        xtbool2 bdenorm,binf,bsetone0,bsetone1;
        xtfloatx2 complexone=AE_SEL32_HH_SX2(XT_CONST_S(1),XT_CONST_S(0));
        xtfloatx2 x0,x1,y0,y1;
        xtfloatx2 z0,z1,f,g;
        AE_LSX2X2_IP(x0,x1,pX,sizeof(xtfloatx4));
        AE_LSX2X2_IP(y0,y1,pY,sizeof(xtfloatx4));
        MULCCONJ_SX2(z0,z1,x0,x1,y0,y1);
        MULCCONJ_SX2(y0,y1,z0,z1,z0,z1);
        f=AE_SEL32_HH_SX2(y0,y1);
        // set to 1 on very small of very big arguments
        bdenorm=XT_OLT_SX2(f,FLT_MIN);
        binf   =XT_UEQ_SX2(f,plusInff.f);
        aset1 = AE_MOVAB2(bdenorm)|AE_MOVAB2(binf);
        bsetone0 =  AE_MOVBA1X2(aset1>>1,aset1>>1);
        bsetone1 =  AE_MOVBA1X2(aset1 &1,aset1&1 );
        g = XT_RSQRT_SX2(f);
        x0=AE_SEL32_HH_SX2(g,g);
        x1=AE_SEL32_LL_SX2(g,g);
        MUL_SX2X2(z0,z1,z0,z1,x0,x1);
        XT_MOVT_SX2(z0,complexone,bsetone0);
        XT_MOVT_SX2(z1,complexone,bsetone1);
        AE_SSX2X2_IP(z0,z1,pZ,sizeof(xtfloatx4));
    }
#else
    for (n=0; n<(Nup>>4); n++)
    {
        xtfloatx2 x0,x1,y0,y1;
        xtfloatx2 z0,z1;
        AE_LSX2X2_IP(x0,x1,pX,sizeof(xtfloatx4));
        AE_LSX2X2_IP(y0,y1,pY,sizeof(xtfloatx4));
        MULCCONJ_SX2(z0,z1,x0,x1,y0,y1);
        AE_SSX2X2_IP(z0,z1,pZ,sizeof(xtfloatx4));
        AE_LSX2X2_IP(x0,x1,pX,sizeof(xtfloatx4));
        AE_LSX2X2_IP(y0,y1,pY,sizeof(xtfloatx4));
        MULCCONJ_SX2(z0,z1,x0,x1,y0,y1);
        AE_SSX2X2_IP(z0,z1,pZ,sizeof(xtfloatx4));
        AE_LSX2X2_IP(x0,x1,pX,sizeof(xtfloatx4));
        AE_LSX2X2_IP(y0,y1,pY,sizeof(xtfloatx4));
        MULCCONJ_SX2(z0,z1,x0,x1,y0,y1);
        AE_SSX2X2_IP(z0,z1,pZ,sizeof(xtfloatx4));
        AE_LSX2X2_IP(x0,x1,pX,sizeof(xtfloatx4));
        AE_LSX2X2_IP(y0,y1,pY,sizeof(xtfloatx4));
        MULCCONJ_SX2(z0,z1,x0,x1,y0,y1);
        AE_SSX2X2_IP(z0,z1,pZ,sizeof(xtfloatx4));
    }
    __Pragma("no_reorder")
    pX=(const xtfloatx4*)a_Z;
    pZ=(      xtfloatx4*)a_Z;
    __Pragma("loop_count min=4, factor=4")
    for (n=0; n<(Nup>>2); n++)
    {
        int aset1;
        xtbool2 bdenorm,binf,bsetone0,bsetone1;
        xtfloatx2 complexone=AE_SEL32_HH_SX2(XT_CONST_S(1),XT_CONST_S(0));
        xtfloatx2 x0,x1,y0,y1;
        xtfloatx2 z0,z1,f,g;
        AE_LSX2X2_IP(z0,z1,pX,sizeof(xtfloatx4));
        MULCCONJ_SX2(y0,y1,z0,z1,z0,z1);
        f=AE_SEL32_HH_SX2(y0,y1);
        // set to 1 on very small of very big arguments
        bdenorm=XT_OLT_SX2(f,FLT_MIN);
        binf   =XT_UEQ_SX2(f,plusInff.f);
        aset1 = AE_MOVAB2(bdenorm)|AE_MOVAB2(binf);
        bsetone0 =  AE_MOVBA1X2(aset1>>1,aset1>>1);
        bsetone1 =  AE_MOVBA1X2(aset1 &1,aset1&1 );
        g = XT_RSQRT_SX2(f);
        x0=AE_SEL32_HH_SX2(g,g);
        x1=AE_SEL32_LL_SX2(g,g);
        MUL_SX2X2(z0,z1,z0,z1,x0,x1);
        XT_MOVT_SX2(z0,complexone,bsetone0);
        XT_MOVT_SX2(z1,complexone,bsetone1);
        AE_SSX2X2_IP(z0,z1,pZ,sizeof(xtfloatx4));
    }
#endif
    // last odd sample
    {
        int aset1;
        xtbool2 binf,bdenorm,bsetone;
        xtfloatx2 complexone=AE_SEL32_HH_SX2(XT_CONST_S(1),XT_CONST_S(0));
        xtfloatx2 x0,y0,z0,f,g;
        x0=XT_LSX2I((const xtfloatx2*)(a_X+Nup/2),0);
        y0=XT_LSX2I((const xtfloatx2*)(a_Y+Nup/2),0);
        z0=MULCCONJ_S(x0,y0);
        y0=MULCCONJ_S(z0,z0);
        f=AE_SEL32_HH_SX2(y0,y0);
        bdenorm=XT_OLT_S(f,FLT_MIN);
        binf   =XT_UEQ_S(f,plusInff.f);
        aset1 = AE_MOVAB2(bdenorm)|AE_MOVAB2(binf);
        bsetone =  AE_MOVBA1X2(aset1>>1,aset1>>1);
        g=XT_RSQRT_S(f); 
        x0=AE_SEL32_HH_SX2(g,g);
        z0=MUL_SX2(z0,x0);
        XT_MOVT_SX2(z0,complexone,bsetone);
        XT_SSX2I(z0,(xtfloatx2*)pZ,0);
    }
    /* z <- ifft(Z) */
    ifft_realf_ie(a_z, a_Z, twdTbl, twdStp, Nup);
    /* Copy results to the output argument z, if necessary. */
    if (z!=a_z) {
        pX=(const xtfloatx4*)a_z;
        pZ=(      xtfloatx4*)z;
        __Pragma("loop_count min=4")
        for ( n=0; n<(N>>3); n++ ) 
        {
            xtfloatx2 t0,t1,t2,t3;
            AE_LSX2X2_IP(t0,t1,pX,sizeof(xtfloatx4));
            AE_LSX2X2_IP(t2,t3,pX,sizeof(xtfloatx4));
            AE_SSX2X2_IP(t0,t1,pZ,sizeof(xtfloatx4));
            AE_SSX2X2_IP(t2,t3,pZ,sizeof(xtfloatx4));
        }
    }
} /* gccphatf() */

/* Returns: size of scratch memory area, in bytes. */
size_t gccphatf_getScratchSize(int N)
{
    int Nup = 1<<(31-NSA(N-1)); /* Round N up to the next power of 2. */
    NASSERT(0==(N%32) && 64<=N && N<=320);
    return Nup*sz_f32 +                    /* a_x[Nup]     */
           (Nup/2+1)*sz_cf32 + ALIGN_PAD + /* a_X[Nup/2+1] */
           (Nup/2+1)*sz_cf32 + ALIGN_PAD;  /* a_Y[Nup/2+1] */
}
#else
// code for SFPU
void gccphatf(void* pScr, float32_t * z, const float32_t * x, const float32_t * y, int N)
{
    const ae_int32x4* restrict pX;
    const ae_int32x4* restrict pY;
          ae_int32x4* restrict pZ;

    float32_t *a_x, *a_y, *a_z; /* Time-domain representation of vectors x, y and z. */
    complex_float *a_X, *a_Y, *a_Z; /* Frequency-domain representation of vectors x, y and z. */
    const complex_float * twdTbl = (complex_float*)gccphatf_twd_tbl;
    int logNup, Nup;
    int twdStp, n;
    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(z, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT(0==(N%32) && 64<=N && N<=320);
    NASSERT(N<=GCCPHATF_TWD_TBL_SZ);
    /* Round N up to the next power of 2 and compute the step size for the FFT twiddles. */
    logNup = 31-NSA(N-1); Nup = 1<<logNup;
    twdStp = 1<<(GCCPHATF_TWD_TBL_LOG_SZ-logNup);
    /* Partition the scratch memory area */
    {
        void *p = pScr;
        a_x = (float32_t    *)p; p = ALIGN_PTR(a_x + Nup);
        a_X = (complex_float*)p; p = ALIGN_PTR(a_X + Nup/2+1);
        a_Y = (complex_float*)p; p = ALIGN_PTR(a_Y + Nup/2+1);
#ifdef _DEBUG
        NASSERT((int8_t*)p-(int8_t*)pScr <= (int)gccphatf_getScratchSize(N));
#endif
        /* Partially re-use the allocated arrays. */
        a_y = a_x; a_Z = a_X;
        a_z = N==Nup ? z : a_x;
    }
    /* Copy data from the input argument x to scratch array, to feel free for
     * the FFT routine to distort them. Then, pad the FFT window with zeros. */
    NASSERT(N%32==0 && Nup%32==0);
    pX=(ae_int32x4*)x;
    pZ=(ae_int32x4*)a_x;
    for ( n=0; n<(N>>3); n++ ) 
    {
        ae_int32x2 t0,t1,t2,t3;
        AE_L32X2X2_IP(t0,t1,pX,sizeof(ae_int32x4));
        AE_L32X2X2_IP(t2,t3,pX,sizeof(ae_int32x4));
        AE_S32X2X2_IP(t0,t1,pZ,sizeof(ae_int32x4));
        AE_S32X2X2_IP(t2,t3,pZ,sizeof(ae_int32x4));
    }
    for ( ; n<(Nup>>3); n++ ) 
    {
        AE_S32X2X2_IP(0,0,pZ,sizeof(ae_int32x4));
        AE_S32X2X2_IP(0,0,pZ,sizeof(ae_int32x4));
    }
    /* X <- fft(x) */
    fft_realf_ie(a_X, a_x, twdTbl, twdStp, Nup);
    /* Copy data from the argument y to the FFT window and pad it with zeros. */
    NASSERT(N%32==0 && Nup%32==0);
    pX=(ae_int32x4*)y;
    pZ=(ae_int32x4*)a_y;
    for ( n=0; n<(N>>3); n++ ) 
    {
        ae_int32x2 t0,t1,t2,t3;
        AE_L32X2X2_IP(t0,t1,pX,sizeof(ae_int32x4));
        AE_L32X2X2_IP(t2,t3,pX,sizeof(ae_int32x4));
        AE_S32X2X2_IP(t0,t1,pZ,sizeof(ae_int32x4));
        AE_S32X2X2_IP(t2,t3,pZ,sizeof(ae_int32x4));
    }
    for ( ; n<(Nup>>3); n++ ) 
    {
        AE_S32X2X2_IP(0,0,pZ,sizeof(ae_int32x4));
        AE_S32X2X2_IP(0,0,pZ,sizeof(ae_int32x4));
    }
    /* Y <- fft(y); */
    fft_realf_ie(a_Y, a_y, twdTbl, twdStp, Nup);
    /* Z <- X.*conj(Y) */
    /* Z <- exp(1j*angle(Z)) */
    pX=(const ae_int32x4*)a_X;
    pY=(const ae_int32x4*)a_Y;
    pZ=(ae_int32x4*)a_Z;
    for (n=0; n<Nup/2+1; n++)
    {
        xtbool binf,bdenorm,bsetone;
        xtfloat f,zr,zi,g,xr,xi,yr,yi;
        XT_LSIP(xr,castxcc(xtfloat,pX),sizeof(xtfloat));
        XT_LSIP(xi,castxcc(xtfloat,pX),sizeof(xtfloat));
        XT_LSIP(yr,castxcc(xtfloat,pY),sizeof(xtfloat));
        XT_LSIP(yi,castxcc(xtfloat,pY),sizeof(xtfloat));
        zr=XT_MUL_S(xr,yr); XT_MADD_S(zr,xi,yi);
        zi=XT_MUL_S(xi,yr); XT_MSUB_S(zi,xr,yi);

        f = XT_MUL_S(zr,zr);
        XT_MADD_S(f,zi,zi);
        bdenorm=XT_OLT_S(f,FLT_MIN);
        binf   =XT_UEQ_S(f,plusInff.f);
        g=XT_RSQRT_S(f); 
        zr=XT_MUL_S(zr,g); 
        zi=XT_MUL_S(zi,g);
        bsetone=AE_MOVBA(AE_MOVAB(bdenorm)|AE_MOVAB(binf));
        XT_MOVT_S(zr,XT_CONST_S(1),bsetone);
        XT_MOVT_S(zi,XT_CONST_S(0),bsetone);
        XT_SSIP(zr,castxcc(xtfloat,pZ),sizeof(xtfloat));
        XT_SSIP(zi,castxcc(xtfloat,pZ),sizeof(xtfloat));
    }
    /* z <- ifft(Z) */
    ifft_realf_ie(a_z, a_Z, twdTbl, twdStp, Nup);
    /* Copy results to the output argument z, if necessary. */
    if (z!=a_z) {
        pX=(ae_int32x4*)a_z;
        pZ=(ae_int32x4*)z;
        for ( n=0; n<(N>>3); n++ ) 
        {
            ae_int32x2 t0,t1,t2,t3;
            AE_L32X2X2_IP(t0,t1,pX,sizeof(ae_int32x4));
            AE_L32X2X2_IP(t2,t3,pX,sizeof(ae_int32x4));
            AE_S32X2X2_IP(t0,t1,pZ,sizeof(ae_int32x4));
            AE_S32X2X2_IP(t2,t3,pZ,sizeof(ae_int32x4));
        }
    }
} /* gccphatf() */

/* Returns: size of scratch memory area, in bytes. */
size_t gccphatf_getScratchSize(int N)
{
    int Nup = 1<<(31-NSA(N-1)); /* Round N up to the next power of 2. */
    NASSERT(0==(N%32) && 64<=N && N<=320);
    return Nup*sz_f32 +                    /* a_x[Nup]     */
           (Nup/2+1)*sz_cf32 + ALIGN_PAD + /* a_X[Nup/2+1] */
           (Nup/2+1)*sz_cf32 + ALIGN_PAD;  /* a_Y[Nup/2+1] */
}
#endif
