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
 * Power spectrum, 16x32 precision
 * Code optimized for HiFi5 core
 */
#include "NatureDSP_Signal_fft.h"
#include "common.h"

static const int32_t ALIGN(32) polysqrt[]=
                                  {-1564282316, 1610612736, -139612678, // Q29
                                    893795525, -832682052, // Q31
                                    1239415238};           // last coefficient+1, Q31

inline_ void sqrtQ31(ae_int32x2* px)
{
    ae_int32x2 x0=*px,r,y,r2,d,t;
    /* first compute 0.5/sqrt(x), Q31 via polynomial */
    d=AE_SUB32(x0,1610612736);
    r=polysqrt[0]; // q29
    r=AE_MULFP32X2RAS(d,r);
    t=polysqrt[2]; AE_MULAFP32X2RAS(t,d,r); r=t;
    r=AE_SLLI32S(r,2); 
    t=polysqrt[3]; AE_MULAFP32X2RAS(t,d,r); r=t;
    t=polysqrt[4]; AE_MULAFP32X2RAS(t,d,r); r=t;
    t=polysqrt[5]; AE_MULAFP32X2RAS(t,d,r); r=t;
    /* next iterate  0.5/sqrt(x) */
    r2=AE_MULFP32X2RAS(r,r);
    d =AE_MULFP32X2RAS(x0,r2);
    d =AE_SUB32((1<<30),AE_ADD32(d,d));
    AE_MULAFP32X2RAS(r,d,r);
    /* last iteration for sqrt(x) is done in higher precision */
    y=AE_SLLI32(AE_MULFP32X2RAS(x0,r),1);
    /* truncate, not saturate here!! */
    d=AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(AE_SLLI64(AE_MUL32_HH(y,r),5)),AE_MOVINT32X2_FROMINT64(AE_SLLI64(AE_MUL32_LL(y,r),5)));
    d=AE_MULFP32X2RAS(d,y);
    AE_MULSFP32X2RAS(y,d,0x08000000);
    *px=y;
}

inline_ void sqrtQ31_sx2x2(ae_int32x2* px0,ae_int32x2* px1)
{
    const ae_int32* pPoly=(const ae_int32*)polysqrt;
    ae_int32x2 x0=*px0,x1=*px1,r0,r1,y0,y1,r20,r21,d0,d1,t0,t1;
    ae_int64 a0,a1;
    /* first compute 0.5/sqrt(x), Q31 via polynomial */
#if 0
    d0=AE_SUB32(x0,1610612736); 
    d1=AE_SUB32(x1,1610612736);
    r0=r1=polysqrt[0]; // q29
    AE_MULF2P32X4RAS(r0,r1,d0,d1,r0,r1);
    t0=t1=polysqrt[2]; AE_MULAF2P32X4RAS(t0,t1,d0,d1,r0,r1); AE_MOVD32X4(r0,r1,t0,t1);
    r0=AE_SLLI32S(r0,2);  r1=AE_SLLI32S(r1,2); 
    t0=t1=polysqrt[3]; AE_MULAF2P32X4RAS(t0,t1,d0,d1,r0,r1); AE_MOVD32X4(r0,r1,t0,t1);
    t0=t1=polysqrt[4]; AE_MULAF2P32X4RAS(t0,t1,d0,d1,r0,r1); AE_MOVD32X4(r0,r1,t0,t1);
    t0=t1=polysqrt[5]; AE_MULAF2P32X4RAS(t0,t1,d0,d1,r0,r1); AE_MOVD32X4(r0,r1,t0,t1);
#else
    r0=AE_L32_I(pPoly,1*sizeof(int32_t));
    d0=AE_SUB32(x0,r0); 
    d1=AE_SUB32(x1,r0);
    AE_L32_XP(r0,pPoly,2*sizeof(int32_t));r1=r0;
    AE_MULF2P32X4RAS(r0,r1,d0,d1,r0,r1);
    AE_L32_IP(t0,pPoly,sizeof(int32_t)); t1=t0;
    AE_MULAF2P32X4RAS(t0,t1,d0,d1,r0,r1); AE_MOVD32X4(r0,r1,t0,t1);
    r0=AE_SLLI32S(r0,2);  r1=AE_SLLI32S(r1,2); 
    AE_L32_IP(t0,pPoly,sizeof(int32_t)); t1=t0;
    AE_MULAF2P32X4RAS(t0,t1,d0,d1,r0,r1); AE_MOVD32X4(r0,r1,t0,t1);
    AE_L32_IP(t0,pPoly,sizeof(int32_t)); t1=t0;
    AE_MULAF2P32X4RAS(t0,t1,d0,d1,r0,r1); AE_MOVD32X4(r0,r1,t0,t1);
    AE_L32_IP(t0,pPoly,sizeof(int32_t)); t1=t0;
    AE_MULAF2P32X4RAS(t0,t1,d0,d1,r0,r1); AE_MOVD32X4(r0,r1,t0,t1);
#endif
    /* next iterate  0.5/sqrt(x) */
    AE_MULF2P32X4RAS(r20,r21,r0,r1,r0,r1);
    AE_MULF2P32X4RAS(d0,d1,x0,x1,r20,r21);
    AE_MOVD32X4(t0,t1,1<<30,1<<30);
    AE_MULS2P32X4(t0,t1,d0,d1,2,2);
    AE_MULAF2P32X4RAS(r0,r1,t0,t1,r0,r1);
    /* last iteration for sqrt(x) is done in higher precision */
    AE_MULF2P32X4RAS(y0,y1,x0,x1,r0,r1);
    y0=AE_SLLI32(y0,1); 
    y1=AE_SLLI32(y1,1);
    /* truncate, not saturate here!! */
    AE_MUL32X2S_HH_LL(a0,a1,y0,r0);d0=AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(AE_SLLI64(a0,5)),AE_MOVINT32X2_FROMINT64(AE_SLLI64(a1,5)));
    AE_MUL32X2S_HH_LL(a0,a1,y1,r1);d1=AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(AE_SLLI64(a0,5)),AE_MOVINT32X2_FROMINT64(AE_SLLI64(a1,5)));
    AE_MULF2P32X4RAS(d0,d1,d0,d1,y0,y1);
    AE_MULAF2P32X4RAS(y0,y1,d0,d1,0xF8000000,0xF8000000);
    *px0=y0; *px1=y1;
}
/*-------------------------------------------------------------------------
  Power spectrum
  These functions compute a normalized power spectrum from the output signal 
  generated by an FFT function. The N argument specifies the size of the FFT 
  and must be a power of 2. The  mode argument is used to specify the type 
  of FFT function used to generate the x array. If the x array has been 
  generated from a frequency-domain complex input signal (output of complex 
  FFT function), the mode argument must be set to 0. Otherwise the mode 
  argument must be set to 1 to signify that the x array has been generated 
  from a frequency-domain real input signal (output of real FFT function).
  The block_exponent argument is used to control the normalization of the 
  power spectrum. It will usually be set to the block_exponent that is 
  returned by corresponding FFT functions.  If the input array was 
  generated by some other means, then the value specified for the 
  block_exponent argument will depend upon how the FFT was calculated. 
  If the function used to calculate the FFT did not scale the intermediate 
  results at any of the stages of the computation, then set block_exponent 
  to zero; if the FFT function scaled the intermediate results at each 
  stage of the computation, then set block_exponent to -1; otherwise set 
  block_exponent to the sum of negated base-2 logarithm of all scaling 
  factors applied to data at intermediate FFT stages. This value will be 
  in the range 0 to log2(N).
  fft_spectrum functions write the power spectrum to the output array y. 
  If mode is set to 0, then the length of the power spectrum will be N. If 
  mode is set to 1, then the length of the power spectrum will be (N/2+1)

  Precision:
  16x32   16-bit inputs, 32-bit outputs
  32x32   32-bit inputs/outputs
  f       floating point inputs/outputs. Requires VFPU/SFPU core option


  Input:
  for mode==0
  x[N]           input spectrum . Real and imaginary
                 data are interleaved and real data goes first:
  for mode==1
  x[N/2+1]       input spectrum (positive side). Real and imaginary
                 data are interleaved and real data goes first:
  block_exponent power spectrum normalization control
  N              FFT size
  mode           power spectrum mode:
                 0 – complex signal
                 1 – real signal

  twdstep               twiddle step
  scalingOpt            scaling option (see table above), not applicable 
                        to the floating point function
  Output:
  for mode==0
  y[N]           output power spectrum
  for mode==1
  y[N/2+1]       output power spectrum
  Returned value:  none

  Restrictions:
  x,y   should not overlap
  x,y   aligned on 16-bytes boundary
-------------------------------------------------------------------------*/
void fft_spectrum16x32( fract32 * restrict y, const complex_fract16 * restrict x, 
                         int N, int block_exponent, int mode )
{
    ae_valignx2 aX;
    ae_valignx2 aY;
    /*
    * Accuracy: 2 ULP
    */
    int n, logN;
    NASSERT( x );
    NASSERT( y );
    logN = 30 - NSA( N );
    if ( N<2 || 0 != (N&(N-1)) || block_exponent < -1 || block_exponent > logN || ( mode != 0 && mode != 1 ) ) return;

    /* Negative block exponent value corresponds to static scaling, i.e. the spectrum
    * has been divided by N in the course of FFT computation. */
    if ( block_exponent<0 ) block_exponent = logN;
    N=( mode ? N/2+1 : N );
    aX=AE_LA128_PP(x);
    aY=AE_ZALIGN128();
    for ( n=0; n<(N>>2); n++)
    {
        ae_int16x4 x01,x23;
        ae_int32x2 a01,a23,w01,w23;
        ae_int16x4 nsa0123,sh0123;
        AE_LA16X4X2_IP(x01,x23,aX,castxcc(ae_int16x8,x));
        a01=AE_MULZAA2D16SS_HH_LL(x01,x01);
        a23=AE_MULZAA2D16SS_HH_LL(x23,x23);
        /* Even normalization is required by sqrt algorithm. */
        nsa0123=AE_NSA32X4(a01,a23); nsa0123=AE_AND16(nsa0123,~1);
        AE_CVTI32X4F16(w01,w23,AE_NEG16S(nsa0123),0);
        w01=AE_SRAV32RS(a01,w01);
        w23=AE_SRAV32RS(a23,w23);
        sqrtQ31_sx2x2(&w01,&w23);
        sh0123=AE_SUB16(AE_SRAI16(nsa0123,1),AE_MOVDA16(mode + block_exponent - logN + 1));
        // multiply by sqrt(0.5) and shift to the right position
        AE_MULF2P32X4RAS(w01,w23,w01,w23,1518500250,1518500250);
        AE_CVTI32X4F16(a01,a23,sh0123,0);
        w01=AE_SRAV32RS(w01,a01);
        w23=AE_SRAV32RS(w23,a23);
        AE_SA32X2X2_IP(w01,w23,aY,castxcc(ae_int32x4,y));
    }
    AE_SA128POS_FP(aY,y);
    __Pragma("loop_count max=3")
    for ( n=0; n<(N&3); n++ )
    {
        ae_int16x4 x0;
        ae_int32x2 a;
        ae_int16x4 nsa;
        AE_L32_IP(a,castxcc(ae_int32,x),sizeof(complex_fract16)); x0=AE_MOVINT16X4_FROMINT32X2(a);
        a=AE_MULZAA2D16SS_HH_LL(x0,x0);
        /* Even normalization is required by sqrt algorithm. */
        nsa=AE_AND16(AE_NSA32X4(a,a),~1);
        a=AE_SRAV32RS(a,AE_SEXT32X2D16_32(AE_NEG16S(nsa)));
        sqrtQ31(&a);
        nsa=AE_SUB16(AE_SRAI16(nsa,1),AE_MOVDA16(mode + block_exponent - logN + 1));
        // multiply by sqrt(0.5) and shift to the right position
        a=AE_MULFP32X2RAS(a,1518500250);
        a=AE_SRAV32RS(a,AE_SEXT32X2D16_32(nsa));
        AE_S32_L_IP(a,castxcc(ae_int32,y),sizeof(int32_t));
    }
}
