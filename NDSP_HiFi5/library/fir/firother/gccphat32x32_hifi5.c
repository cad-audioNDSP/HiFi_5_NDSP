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
    Real 32-bit input/output data
    C code optimized for HiFi5 core
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
#include "NatureDSP_Signal_fft.h"

#define ALIGN_SIZE     (HIFI_SIMD_WIDTH)
#define ALIGN_PAD      (ALIGN_SIZE-1)
#define ALIGN_PTR(p)   (void*)(((uintptr_t)(p)+ALIGN_PAD)&~ALIGN_PAD)
#define sz_i16         sizeof(int16_t)
#define sz_i32         sizeof(int32_t)
#define sz_ci32        sizeof(complex_fract32)

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
int gccphat32x32 (void* pScr, int32_t * z, const int32_t * x, const int32_t * y, int N)
{
    const ae_int32x4* restrict pX;
    const ae_int32x4* restrict pY;
          ae_int32x4* restrict pZ;
          ae_int32x2* restrict pF;
          ae_int16x4* restrict pE;
          ae_int16x4* restrict pG;

    int32_t *a_x, *a_y, *a_z; /* Time-domain representation of vectors x, y and z. */
    complex_fract32 *a_X, *a_Y, *a_Z; /* Frequency-domain representation of vectors x, y and z. */
    int32_t *a_f, *a_g; /* Fractional part of |Z|.^2, 1./|Z| */
    int16_t *a_ef, *a_eg; /* Exponential part of |Z|.^2, 1./|Z| */
    typedef const fft_handle_t chandle;
    chandle hfft[2][4] = {{rfft32_64 , rfft32_128 , rfft32_256 , rfft32_512 },
                          {rifft32_64, rifft32_128, rifft32_256, rifft32_512}};
    int logNup, Nup;
    int ez, n;
    NASSERT_ALIGN(pScr, ALIGN_SIZE);
    NASSERT_ALIGN(z, ALIGN_SIZE);
    NASSERT_ALIGN(x, ALIGN_SIZE);
    NASSERT_ALIGN(y, ALIGN_SIZE);
    NASSERT(0==(N%32) && 64<=N && N<=320);
    /* Round N up to the next power of 2 and compute the step size for the FFT twiddles. */
    logNup = 31-NSA(N-1); Nup = 1<<logNup;
    /* Partition the scratch memory area */
    {
        void *p = pScr;
        a_x  = (int32_t        *)p; p = ALIGN_PTR(a_x  + Nup);
        a_X  = (complex_fract32*)p; p = ALIGN_PTR(a_X  + Nup/2+1);
        a_Y  = (complex_fract32*)p; p = ALIGN_PTR(a_Y  + Nup/2+1);
        a_ef = (int16_t        *)p; p = ALIGN_PTR(a_ef + Nup/2+1);
#ifdef _DEBUG
        NASSERT((int8_t*)p-(int8_t*)pScr <= (int)gccphat32x32_getScratchSize(N));
#endif
        /* Partially re-use the allocated arrays. */
        a_f = a_y = a_x; a_Z = a_X; a_g = (int32_t*)a_Y; 
        a_eg = (int16_t*)ALIGN_PTR(a_g + Nup/2+1);
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
    /* X <- fft(x) w/ 32-bit dynamic scaling. Total shift amount returned by the
     * FFT routine is discarded because here we are concerned only with the phase
     * of a spectrum bin. */
    fft_real32x32((int32_t*)a_X, a_x, hfft[0][logNup-6], 2);
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
    fft_real32x32((int32_t*)a_Y, a_y, hfft[0][logNup-6], 2);
    /* Compute pointwise product of complex input vectors X[len] and Y[len], with the factor Y being conjugated:
    Z <- X.*conj(Y);
    We should care most for the phase of the product, so multiplicands may be subject to appropriate scaling.
    Besides of pointwise multiplication, compute the squared absolute value for each product, and normalize
    it for the forthcoming rsqrt computation:
    f*2^ef <- |Z|^2, 0.25<=f<1.0, ef odd.
    Note that any of the input vectors X[len] and Y[len] may be reused for the resulting vector Z[len]! */
    pX=(const ae_int32x4*)a_X;
    pY=(const ae_int32x4*)a_Y;
    pZ=(      ae_int32x4*)a_Z;
    pF=(      ae_int32x2*)a_f;
    __Pragma("loop_count min=4,factor=4")
    for ( n=0; n<(Nup>>2); n++) 
    {
        ae_int64 wr0,wr1,wi0,wi1;
        ae_int32x2 x0,x1,y0,y1,w0,w1;
        ae_int16x4 nsa_;
        int nsa0,nsa1;
        AE_L32X2X2_IP(x0,x1,pX,sizeof(ae_int32x4));
        AE_L32X2X2_IP(y0,y1,pY,sizeof(ae_int32x4));
        nsa_=AE_NSA32X4(x0,x1);
        nsa_=AE_MIN16(AE_SEL16_2301(nsa_,nsa_),nsa_);
        AE_CVTI32X4F16(w0,w1,AE_NEG16S(nsa_),0);
        x0=AE_SRAV32RS(x0,w0);
        x1=AE_SRAV32RS(x1,w1);
        nsa_=AE_NSA32X4(y0,y1);
        nsa_=AE_MIN16(AE_SEL16_2301(nsa_,nsa_),nsa_);
        AE_CVTI32X4F16(w0,w1,AE_NEG16S(nsa_),0);
        y0=AE_SRAV32RS(y0,w0);
        y1=AE_SRAV32RS(y1,w1);
        AE_MULZAAF2D32RA_HH_LL(wr0,wr1,x0,x1,y0,y1);
        AE_MULZASF2D32RA_HL_LH(wi0,wi1,y0,y1,x0,x1);
        x0=AE_TRUNCI32X2F64S(wr0,wi0,15);
        x1=AE_TRUNCI32X2F64S(wr1,wi1,15);
        AE_S32X2X2_IP(x0,x1,pZ,sizeof(ae_int32x4));
        wr0=AE_MULZAAD32S_HH_LL(x0,x0);
        wr1=AE_MULZAAD32S_HH_LL(x1,x1);
        nsa0=AE_NSA64(wr0)&~1;
        nsa1=AE_NSA64(wr1)&~1;
        wr0=AE_SLAA64(wr0,nsa0);
        wr1=AE_SLAA64(wr1,nsa1);
        AE_S32X2_IP(AE_ROUND32X2F64SASYM(wr0,wr1),pF,sizeof(ae_int32x2));
        a_ef[2*n+0]=nsa0-3;
        a_ef[2*n+1]=nsa1-3;
    }
    __Pragma("no_reorder")
    /* g <- 1./sqrt(f);  */
    /* simplified rsqrt32x32: assumes positive input and less accurate last iteration */
    pX=(const ae_int32x4*)a_f;
    pZ=(      ae_int32x4*)a_g;
    pE=(      ae_int16x4*)a_eg;
    __Pragma("loop_count min=2,factor=2")
    for (n = 0; n<(Nup >> 3); n++)
    {
        ae_int32x2 x0,x1,z0,z1,e0,e1,y0,y1,t0,t1;
        ae_int16x4 nsax;
        AE_L32X2X2_IP(x0, x1, pX,sizeof(ae_int32x4));
        nsax = AE_NSA32X4(x0, x1);   

        x0 = AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(AE_AND16(nsax, ~1))));
        x1 = AE_SRAV32RS(x1, AE_SEXT32X2D16_10(AE_NEG16S(AE_AND16(nsax, ~1))));
        nsax = AE_SRAI16(nsax, 1);
        nsax = AE_ADD16(nsax, 1);
        AE_S16X4_IP(nsax, pE, sizeof(ae_int16x4));
    
        y0 = AE_SUB32(0x80000000, AE_SRAI32(x0, 1)); /* Q30 */
        y1 = AE_SUB32(0x80000000, AE_SRAI32(x1, 1)); /* Q30 */
    
        /*1-st iteration*/
        AE_MULF2P32X4RAS(z0, z1, y0, y1, y0, y1);
        t0 = t1 = 0x20000000;
        AE_MULSF2P32X4RAS(t0, t1, x0, x1, z0, z1); e0 = t0; e1 = t1;
        AE_MULF2P32X4RAS(z0, z1, y0, y1, e0, e1);
        AE_MULA2P32X4(y0, y1, AE_MOVDA32(2), AE_MOVDA32(2), z0, z1); 
    
        AE_MULF2P32X4RAS(z0, z1, y0, y1, y0, y1);
        t0 = t1 = 0x20000000;
        AE_MULSF2P32X4RAS(t0, t1, x0, x1, z0, z1); e0 = t0; e1 = t1;
        AE_MULF2P32X4RAS(z0, z1, y0, y1, e0, e1);
        AE_MULA2P32X4(y0, y1, AE_MOVDA32(2), AE_MOVDA32(2), z0, z1);
    
        AE_MULF2P32X4RAS(z0, z1, y0, y1, y0, y1);
        t0 = t1 = 0x20000000;
        AE_MULSF2P32X4RAS(t0, t1, x0, x1, z0, z1); e0 = t0; e1 = t1;
        AE_MULF2P32X4RAS(z0, z1, y0, y1, e0, e1);
        AE_MULA2P32X4(y0, y1, AE_MOVDA32(2), AE_MOVDA32(2), z0, z1);
        AE_S32X2X2_IP(y0, y1, pZ, sizeof(ae_int32x4));
    }
/* Normalize each complex element of the input vector X[N] by its absolute value, in-place:
   Z <- exp(1j*angle(Z));
   For zero input elements the result is set to (1,0).
   Note that the input vector X[len] may be reused for the resulting vector Z[len]. */
    __Pragma("no_reorder")
    pX=(const ae_int32x4*)a_g;
    pY=(const ae_int32x4*)a_Z;
    pZ=(      ae_int32x4*)a_Z;
    pG=(      ae_int16x4*)a_eg;
    pE=(      ae_int16x4*)a_ef;
    __Pragma("loop_count min=2,factor=2")
    for ( n=0; n<(Nup>>3); n++) 
    {
        ae_int32x2 complexone=AE_MOVDA32X2(0x7fffffff,0);
        xtbool2 bzero0,bzero1,bzero2,bzero3;
        ae_int32x2 h0,h1,y0,y1,y2,y3;
        ae_int64 zr0,zi0,zr1,zi1,zr2,zi2,zr3,zi3;
        ae_int16x4 ef,eg,e;
        /* Multiply the reciprocal absolute value by sqrt(2) to account for the
         * oddness of ef. 
         * Q(30-eg-ef/2) = Q30*Q(31-eg-ef/2) - 31 w/ asym. rounding */
        AE_L32X2X2_IP(h0,h1,pX,sizeof(ae_int32x4));
        AE_L32X2X2_IP(y0,y1,pY,sizeof(ae_int32x4));
        AE_L32X2X2_IP(y2,y3,pY,sizeof(ae_int32x4));
        bzero0=AE_EQ32(AE_OR32(y0,AE_SEL32_LH(y0,y0)),0);
        bzero1=AE_EQ32(AE_OR32(y1,AE_SEL32_LH(y1,y1)),0);
        bzero2=AE_EQ32(AE_OR32(y2,AE_SEL32_LH(y2,y2)),0);
        bzero3=AE_EQ32(AE_OR32(y3,AE_SEL32_LH(y3,y3)),0);
        AE_MULF2P32X4RAS(h0,h1,h0,h1,(int32_t)1518500250,(int32_t)1518500250);
        AE_MUL32X2S_HH_LL(zr0,zi0, y0, AE_SEL32_HH(h0,h0));
        AE_MUL32X2S_HH_LL(zr1,zi1, y1, AE_SEL32_LL(h0,h0));
        AE_MUL32X2S_HH_LL(zr2,zi2, y2, AE_SEL32_HH(h1,h1));
        AE_MUL32X2S_HH_LL(zr3,zi3, y3, AE_SEL32_LL(h1,h1));
        AE_L16X4_IP(ef,pE,sizeof(ae_int16x4));
        AE_L16X4_IP(eg,pG,sizeof(ae_int16x4));
        e=AE_ADD16(3,eg);
        e=AE_ADD16(e,AE_SRAI16(ef,1));
        y0=AE_TRUNCA32X2F64S(zr0,zi0,AE_MOVAD16_3(e));
        y1=AE_TRUNCA32X2F64S(zr1,zi1,AE_MOVAD16_2(e));
        y2=AE_TRUNCA32X2F64S(zr2,zi2,AE_MOVAD16_1(e));
        y3=AE_TRUNCA32X2F64S(zr3,zi3,AE_MOVAD16_0(e));
        AE_MOVT32X2(y0,complexone,bzero0);
        AE_MOVT32X2(y1,complexone,bzero1);
        AE_MOVT32X2(y2,complexone,bzero2);
        AE_MOVT32X2(y3,complexone,bzero3);
        AE_S32X2X2_IP(y0,y1,pZ,sizeof(ae_int32x4));
        AE_S32X2X2_IP(y2,y3,pZ,sizeof(ae_int32x4));
    } 
    //-----------------------------------------------------------
    // last odd iteration
    //-----------------------------------------------------------
    {
        ae_int64 wr0,wi0;
        ae_int32x2 h0,x0,y0,z0,e0,t0,w0,w1;
        ae_int16x4 nsa_;
        ae_int16x4 nsax;
        xtbool2 bzero0;
        ae_int64 zr0,zi0;
        int nsa0;
        ae_int32x2 complexone=AE_MOVDA32X2(0x7fffffff,0);

        x0=AE_L32X2_I((const ae_int32x2*)(a_X+Nup/2),0);
        y0=AE_L32X2_I((const ae_int32x2*)(a_Y+Nup/2),0);
        nsa_=AE_NSA32X4(x0,x0);
        nsa_=AE_MIN16(AE_SEL16_2301(nsa_,nsa_),nsa_);
        AE_CVTI32X4F16(w0,w1,AE_NEG16S(nsa_),0);
        x0=AE_SRAV32RS(x0,w0);
        nsa_=AE_NSA32X4(y0,y0);
        nsa_=AE_MIN16(AE_SEL16_2301(nsa_,nsa_),nsa_);
        AE_CVTI32X4F16(w0,w1,AE_NEG16S(nsa_),0);
        y0=AE_SRAV32RS(y0,w0);
        wr0=AE_MULZAAFD32RA_HH_LL(x0,y0);
        wi0=AE_MULZASFD32RA_HL_LH(y0,x0);
        x0=AE_TRUNCI32X2F64S(wr0,wi0,15);
        h0=x0;
        wr0=AE_MULZAAD32S_HH_LL(x0,x0);
        nsa0=AE_NSA64(wr0)&~1;
        wr0=AE_SLAA64(wr0,nsa0);
        // rsqrt
        x0=AE_ROUND32X2F64SASYM(wr0,wr0);
        nsax = AE_NSA32X4(x0, x0);   
        x0 = AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(AE_AND16(nsax, ~1))));
        nsax = AE_SRAI16(nsax, 1);
        nsax = AE_ADD16(nsax, 1);
        y0 = AE_SUB32(0x80000000, AE_SRAI32(x0, 1)); /* Q30 */
        z0=AE_MULFP32X2RAS(y0, y0);t0 = 0x20000000; AE_MULSFP32X2RAS(t0, x0, z0); e0 = t0; 
        z0=AE_MULFP32X2RAS(y0, e0);AE_MULAP32X2(y0,AE_MOVDA32(2), z0); 
        z0=AE_MULFP32X2RAS(y0, y0);t0 = 0x20000000; AE_MULSFP32X2RAS(t0, x0, z0); e0 = t0; 
        z0=AE_MULFP32X2RAS(y0, e0);AE_MULAP32X2(y0,AE_MOVDA32(2), z0); 
        z0=AE_MULFP32X2RAS(y0, y0);t0 = 0x20000000; AE_MULSFP32X2RAS(t0, x0, z0); e0 = t0; 
        z0=AE_MULFP32X2RAS(y0, e0);AE_MULAP32X2(y0,AE_MOVDA32(2), z0); 
        /* angle */
        y0=AE_MULFP32X2RAS(y0,(int32_t)1518500250);
        z0=h0;
        bzero0=AE_EQ32(AE_OR32(z0,AE_SEL32_LH(z0,z0)),0);
        AE_MUL32X2S_HH_LL(zr0,zi0, z0, y0);
        z0=AE_TRUNCA32X2F64S(zr0,zi0,3+AE_MOVAD16_3(nsax)+((nsa0-3)>>1));
        AE_MOVT32X2(z0,complexone,bzero0);
        AE_S32X2_I(z0,(ae_int32x2*)(a_Z+Nup/2),0);
    } 


    /* z <- ifft(Z). Here the dynamic scaling improves the resulting SINAD by only a
     * few dB, so it may be sacrifyed for the sake of cycles performance if the fixed
     * scaling option makes the IFFT transform significantly faster. */
    ez = ifft_real32x32(a_z, (int32_t*)a_Z, hfft[1][logNup-6], 2);
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
    return ez;
} /* gccphat32x32() */

/* Returns: size of scratch memory area, in bytes. */
size_t gccphat32x32_getScratchSize(int N)
{
    int Nup = 1<<(31-NSA(N-1)); /* Round N up to the next power of 2. */
    NASSERT(0==(N%32) && 64<=N && N<=320);
    return Nup*sz_i32 +                    /* a_x [Nup]     */
           (Nup/2+1)*sz_ci32 + ALIGN_PAD + /* a_X [Nup/2+1] */
           (Nup/2+1)*sz_ci32 + ALIGN_PAD + /* a_Y [Nup/2+1] */
           (Nup/2+1)*sz_i16  + ALIGN_PAD;  /* a_ef[Nup/2+1] */
}
