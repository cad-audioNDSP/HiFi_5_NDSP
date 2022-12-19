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
  NatureDSP Signal Processing Library. Audio processing part
    Compute first order regression coefficients
    32-bit fixed-point variant
    C code optimized for HiFi5 core
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_audio.h"

#define ALIGN_SIZE       (HIFI_SIMD_WIDTH)

/*-------------------------------------------------------------------------
  Compute first order regression coefficients as defined in the Hidden
  Markov Models Toolkit (HTK) [1].
  For a set of 2*N+1 frames (each of M static coefficients) specified via 
  the input argument c[(2*N+1)*M], functions compute M delta coefficients
  and store results to the output argument d[M].
  The user may specify any window size N, but the computation is most 
  efficient when N equals 2.
  [1] S. Young, G. Evermann, M. Gales, T. Hain, D. Kershaw, X. Liu, G. Moore,
      J. Odell, D. Ollason, D. Povey, V. Valtchev, P. Woodland,
      The HTK Book (for HTK version 3.4), Chapter 5.9, Eq. (5.16)
      Cambridge University Engineering Department, 2009. 
      http://htk.eng.cam.ac.uk/docs/docs.shtml
      http://www1.icsi.berkeley.edu/Speech/docs/HTKBook/node65_mn.html
  Precision: 
  32x32         32-bit fixed-point input/output data
  f             Floating-point input/output data; requires VFPU/SFPU core option
  Input:
  M             Static coefficients frame size
  N             Window size
  c[(2*N+1)*M]  Static coefficient frames. Coefficients that belong to the
                first frame are stored to c[0..M-1], coefficients for the
                second frame - to c[M..2*M-1] and so forth.
  Output:
  d[M]          Delta coefficients. The number of fractional positions for
                the fixed-point variant is the same as for the input 
                argument c[(2*N+1)*M].
  Restrictions:
  c[],d[]       Must not overlap and must be aligned by 16-bytes
  M             12..30
-------------------------------------------------------------------------*/
#if 0
DISCARD_FUN(void, htkdelta32x32, (int32_t * d,  const int32_t * c, int M, int N));
#else
void htkdelta32x32(int32_t * d, const int32_t * c, int M, int N )
{
    const ae_int32x4* restrict pC0;
    const ae_int32x4* restrict pC1;
    const ae_int32x4* restrict pC3;
    const ae_int32x4* restrict pC4;
    ae_int32x2 s2,r_frac;
    int r_exp;
    int m, n;
    static const short ALIGN(16) dsel_ind[4] = { 1797, 1540, 769, 512 };
    ae_int16x4 dsel; 
    dsel=AE_L16X4_I((const ae_int16x4*)dsel_ind,0);

    NASSERT_ALIGN(d, ALIGN_SIZE);
    NASSERT_ALIGN(c, ALIGN_SIZE);
    NASSERT(12<=M && M<=30);
    /*
     * MATLAB reference:
     *   if N>0
     *      cn = reshape(c,M,2*N+1)';
     *      d = (-N:N)*cn/(2*sum((1:N).^2));
     *   else
     *      d = zeros(1,M);
     *   end
     */
    NASSERT_ALIGN(d, ALIGN_SIZE);
    NASSERT_ALIGN(c, ALIGN_SIZE);
    NASSERT(12<=M && M<=30);
    if (N<=0) 
    {
        for ( m=0; m<M; m++ ) 
        {
            d[m] = 0;
        } /* m */
        return;
    } /* (N<=0) */
    /* special case for N==2 */
    if(N==2)
    {
        ae_valignx2 aC1,aC3;
        pC0=(const ae_int32x4*)(c+0*M);
        pC1=(const ae_int32x4*)(c+1*M);
        pC3=(const ae_int32x4*)(c+3*M);
        pC4=(const ae_int32x4*)(c+4*M);
        aC1=AE_LA128_PP(pC1);
        aC3=AE_LA128_PP(pC3);
        for ( m=0; m<(M>>2); m++) 
        {
            ae_int32x2 t0,t1,c00,c01,c10,c11,c30,c31,c40,c41;
            ae_int64 a0,a1,a2,a3;
            ae_ep aep;
            ae_int64 f0,f1,f2,f3;

            AE_L32X2X2_IP (c00,c01,    castxcc(ae_int32x4,pC0),sizeof(ae_int32x4));
            AE_LA32X2X2_IP(c10,c11,aC1,castxcc(ae_int32x4,pC1));
            AE_LA32X2X2_IP(c30,c31,aC3,castxcc(ae_int32x4,pC3));
            AE_L32X2X2_IP (c40,c41,    castxcc(ae_int32x4,pC4),sizeof(ae_int32x4));
#if 0
            ae_int32x2 c4000h,c4000l;
            ae_int32x2 c4101h,c4101l;
            ae_int32x2 c3010h,c3010l;
            ae_int32x2 c3111h,c3111l;
            ae_int16x4 t4000h,t4000l;
            ae_int16x4 t4101h,t4101l;
            ae_int16x4 t3010h,t3010l;
            ae_int16x4 t3111h,t3111l;
            AE_DSEL16X4(t4000h,t4000l,AE_MOVINT16X4_FROMINT32X2(c40),AE_MOVINT16X4_FROMINT32X2(c00),  dsel);
            AE_DSEL16X4(t4101h,t4101l,AE_MOVINT16X4_FROMINT32X2(c41),AE_MOVINT16X4_FROMINT32X2(c01),  dsel);
            AE_DSEL16X4(t3010h,t3010l,AE_MOVINT16X4_FROMINT32X2(c30),AE_MOVINT16X4_FROMINT32X2(c10),  dsel);
            AE_DSEL16X4(t3111h,t3111l,AE_MOVINT16X4_FROMINT32X2(c31),AE_MOVINT16X4_FROMINT32X2(c11),  dsel);
            c4000h=AE_MOVINT32X2_FROMINT16X4(t4000h); c4000l=AE_MOVINT32X2_FROMINT16X4(t4000l);
            c4101h=AE_MOVINT32X2_FROMINT16X4(t4101h); c4101l=AE_MOVINT32X2_FROMINT16X4(t4101l);
            c3010h=AE_MOVINT32X2_FROMINT16X4(t3010h); c3010l=AE_MOVINT32X2_FROMINT16X4(t3010l);
            c3111h=AE_MOVINT32X2_FROMINT16X4(t3111h); c3111l=AE_MOVINT32X2_FROMINT16X4(t3111l);
            f0=AE_MULZASD32_HL_LH(c4000h,2);
            f1=AE_MULZASD32_HL_LH(c4000l,2);
            AE_MULASD32_HL_LH(f0,c3010h,1);
            AE_MULASD32_HL_LH(f1,c3010l,1);
            f2=AE_MULZASD32_HL_LH(c4101h,2);
            f3=AE_MULZASD32_HL_LH(c4101l,2);
            AE_MULASD32_HL_LH(f2,c3111h,1);
            AE_MULASD32_HL_LH(f3,c3111l,1);
#else
            f0=AE_MULZASD32_HL_LH(AE_SEL32_HH(c40,c00),2);
            f1=AE_MULZASD32_HL_LH(AE_SEL32_LL(c40,c00),2);
            AE_MULASD32_HL_LH(f0,AE_SEL32_HH(c30,c10),1);
            AE_MULASD32_HL_LH(f1,AE_SEL32_LL(c30,c10),1);
            f2=AE_MULZASD32_HL_LH(AE_SEL32_HH(c41,c01),2);
            f3=AE_MULZASD32_HL_LH(AE_SEL32_LL(c41,c01),2);
            AE_MULASD32_HL_LH(f2,AE_SEL32_HH(c31,c11),1);
            AE_MULASD32_HL_LH(f3,AE_SEL32_LL(c31,c11),1);
#endif
            AE_MUL32USEP_LL(aep,a0,AE_MOVINT32X2_FROMINT64(f0),1717986918);
            AE_MUL32USEP_LL(aep,a1,AE_MOVINT32X2_FROMINT64(f1),1717986918);
            AE_MUL32USEP_LL(aep,a2,AE_MOVINT32X2_FROMINT64(f2),1717986918);
            AE_MUL32USEP_LL(aep,a3,AE_MOVINT32X2_FROMINT64(f3),1717986918);
            a0=AE_SRAI72(aep,a0,32);
            a1=AE_SRAI72(aep,a1,32);
            a2=AE_SRAI72(aep,a2,32);
            a3=AE_SRAI72(aep,a3,32);
            AE_MULA32S_HL(a0,AE_MOVINT32X2_FROMINT64(f0),1717986918);
            AE_MULA32S_HL(a1,AE_MOVINT32X2_FROMINT64(f1),1717986918);
            AE_MULA32S_HL(a2,AE_MOVINT32X2_FROMINT64(f2),1717986918);
            AE_MULA32S_HL(a3,AE_MOVINT32X2_FROMINT64(f3),1717986918);
            t0=AE_TRUNCA32X2F64S(a0,a1,30);
            t1=AE_TRUNCA32X2F64S(a2,a3,30);
            AE_S32X2X2_IP(t0,t1,castxcc(ae_int32x4,d),sizeof(ae_int32x4));
        } 
        if(M&2)
        {
            ae_valign aC1,aC3;
            ae_int32x2 t,c0,c1,c3,c4;
            ae_int64 a0,a1;
            ae_ep aep;
            ae_int64 f0,f1;
            aC1=AE_LA64_PP(pC1);
            aC3=AE_LA64_PP(pC3);

            AE_L32X2_IP (c0,    castxcc(ae_int32x2,pC0),sizeof(ae_int32x2));
            AE_LA32X2_IP(c1,aC1,castxcc(ae_int32x2,pC1));
            AE_LA32X2_IP(c3,aC3,castxcc(ae_int32x2,pC3));
            AE_L32X2_IP (c4,    castxcc(ae_int32x2,pC4),sizeof(ae_int32x2));
            f0=AE_MULZASD32_HL_LH(AE_SEL32_HH(c4,c0),2);
            f1=AE_MULZASD32_HL_LH(AE_SEL32_LL(c4,c0),2);
            AE_MULASD32_HL_LH(f0,AE_SEL32_HH(c3,c1),1);
            AE_MULASD32_HL_LH(f1,AE_SEL32_LL(c3,c1),1);

            AE_MUL32USEP_LL(aep,a0,AE_MOVINT32X2_FROMINT64(f0),1717986918);
            AE_MUL32USEP_LL(aep,a1,AE_MOVINT32X2_FROMINT64(f1),1717986918);
            a0=AE_SRAI72(aep,a0,32);
            a1=AE_SRAI72(aep,a1,32);
            AE_MULA32S_HL(a0,AE_MOVINT32X2_FROMINT64(f0),1717986918);
            AE_MULA32S_HL(a1,AE_MOVINT32X2_FROMINT64(f1),1717986918);
            t=AE_TRUNCA32X2F64S(a0,a1,30);
            AE_S32X2_IP(t,castxcc(ae_int32x2,d),sizeof(ae_int32x2));
            m+=2;
        } 
        if(M&1)
        {
            ae_int32x2 c0,c1,c3,c4,t;
            ae_int64 a;
            ae_ep aep;
            ae_int64 f;
            AE_L32_IP (c0,castxcc(ae_int32,pC0),sizeof(ae_int32x2));
            AE_L32_IP (c1,castxcc(ae_int32,pC1),sizeof(ae_int32x2));
            AE_L32_IP (c3,castxcc(ae_int32,pC3),sizeof(ae_int32x2));
            AE_L32_IP (c4,castxcc(ae_int32,pC4),sizeof(ae_int32x2));
            f=AE_MULZASD32_HL_LH(AE_SEL32_HH(c4,c0),2);
            AE_MULASD32_HL_LH(f,AE_SEL32_HH(c3,c1),1);

            AE_MUL32USEP_LL(aep,a,AE_MOVINT32X2_FROMINT64(f),1717986918);
            a=AE_SRAI72(aep,a,32);
            AE_MULA32S_HL(a,AE_MOVINT32X2_FROMINT64(f),1717986918);
            t=AE_TRUNCA32X2F64S(a,a,30);
            AE_S32_L_IP(t,castxcc(ae_int32,d),sizeof(int32_t));
        } 
        return;
    }
    s2 = AE_MOVDA32(N*(N+1)*(2*N+1)); /* Q0 */
    /* compute reciprocal Q(62-r_exp) = Q62/Q0 - r_exp; r_exp<=30 */
    {
        int exponent;
        ae_int32x2 e, y;
        exponent = AE_NSA32_L(s2);
        s2 = AE_SLAA32(s2, exponent);// x in 0.5..1
        r_exp = exponent + 1;
        /* first approximation */
        y = AE_SUB32((int32_t)0xBAEC0000,s2);
        /* 4 iterations to achieve 1 LSB */
        e=0x40000000;AE_MULSFP32X2RAS(e,s2,y); e=AE_ADD32(e,e); AE_MULAFP32X2RAS(y,y,e);
        e=0x40000000;AE_MULSFP32X2RAS(e,s2,y); e=AE_ADD32(e,e); AE_MULAFP32X2RAS(y,y,e);
        e=0x40000000;AE_MULSFP32X2RAS(e,s2,y); e=AE_ADD32(e,e); AE_MULAFP32X2RAS(y,y,e);
        e=0x40000000;AE_MULSFP32X2RAS(e,s2,y); e=AE_ADD32(e,e); AE_MULAFP32X2RAS(y,y,e);
        r_frac = y;
    }
    for ( m=0; m<(M>>2); m++) 
    {
        ae_int32x2 t,t0,t1;
        ae_int64 a0,a1,a2,a3;
        ae_ep aep;
        ae_int64 f0,f1,f2,f3;
        t=N;
        f0=f1=f2=f3=0;
        for ( n=0; n<N; n++ ) 
        {
            ae_int16x4 t1000h,t1000l,t1101h,t1101l;
            ae_int32x2 c1000h,c1000l,c1101h,c1101l;
            ae_int32x2 c00,c01,c10,c11;
            ae_valignx2 aC0,aC1;
            pC0=(const ae_int32x4*)(c+n*M);
            pC1=(const ae_int32x4*)(c+(2*N-n)*M);
            aC0=AE_LA128_PP(pC0);
            aC1=AE_LA128_PP(pC1);
            AE_LA32X2X2_IP(c00,c01,aC0,pC0);
            AE_LA32X2X2_IP(c10,c11,aC1,pC1);
            AE_DSEL16X4(t1000h,t1000l,AE_MOVINT16X4_FROMINT32X2(c10),AE_MOVINT16X4_FROMINT32X2(c00),  dsel);
            AE_DSEL16X4(t1101h,t1101l,AE_MOVINT16X4_FROMINT32X2(c11),AE_MOVINT16X4_FROMINT32X2(c01),  dsel);
            c1000h=AE_MOVINT32X2_FROMINT16X4(t1000h); c1000l=AE_MOVINT32X2_FROMINT16X4(t1000l);
            c1101h=AE_MOVINT32X2_FROMINT16X4(t1101h); c1101l=AE_MOVINT32X2_FROMINT16X4(t1101l);
            AE_MULASD32_HL_LH(f0,c1000h,t);
            AE_MULASD32_HL_LH(f1,c1000l,t);
            AE_MULASD32_HL_LH(f2,c1101h,t);
            AE_MULASD32_HL_LH(f3,c1101l,t);
            t= AE_SUB32(t,1);
        } 
        f0=AE_ADD64(f0,AE_SLAI64(f0,1));
        f1=AE_ADD64(f1,AE_SLAI64(f1,1));
        f2=AE_ADD64(f2,AE_SLAI64(f2,1));
        f3=AE_ADD64(f3,AE_SLAI64(f3,1));
        AE_MUL32USEP_LL(aep,a0,AE_MOVINT32X2_FROMINT64(f0),r_frac);
        AE_MUL32USEP_LL(aep,a1,AE_MOVINT32X2_FROMINT64(f1),r_frac);
        AE_MUL32USEP_LL(aep,a2,AE_MOVINT32X2_FROMINT64(f2),r_frac);
        AE_MUL32USEP_LL(aep,a3,AE_MOVINT32X2_FROMINT64(f3),r_frac);
        a0=AE_SRAI72(aep,a0,32);
        a1=AE_SRAI72(aep,a1,32);
        a2=AE_SRAI72(aep,a2,32);
        a3=AE_SRAI72(aep,a3,32);
        AE_MULA32S_HL(a0,AE_MOVINT32X2_FROMINT64(f0),r_frac);
        AE_MULA32S_HL(a1,AE_MOVINT32X2_FROMINT64(f1),r_frac);
        AE_MULA32S_HL(a2,AE_MOVINT32X2_FROMINT64(f2),r_frac);
        AE_MULA32S_HL(a3,AE_MOVINT32X2_FROMINT64(f3),r_frac);
        t0=AE_TRUNCA32X2F64S(a0,a1,r_exp+2);
        t1=AE_TRUNCA32X2F64S(a2,a3,r_exp+2);
        AE_S32X2X2_IP(t0,t1,castxcc(ae_int32x4,d),sizeof(ae_int32x4));
        c+=4;
    } 
    if (M&2)
    {
        ae_int32x2 t,t0;
        ae_int64 a0,a1;
        ae_ep aep;
        ae_int64 f0,f1;
        t=N;
        f0=f1=0;
        for ( n=0; n<N; n++ ) 
        {
            ae_int32x2 c00,c10;
            ae_valign aC0,aC1;
            pC0=(const ae_int32x4*)(c+n*M);
            pC1=(const ae_int32x4*)(c+(2*N-n)*M);
            aC0=AE_LA64_PP(pC0);
            aC1=AE_LA64_PP(pC1);
            AE_LA32X2_IP(c00,aC0,castxcc(ae_int32x2,pC0));
            AE_LA32X2_IP(c10,aC1,castxcc(ae_int32x2,pC1));
            AE_MULASD32_HL_LH(f0,AE_SEL32_HH(c10,c00),t);
            AE_MULASD32_HL_LH(f1,AE_SEL32_LL(c10,c00),t);
            t= AE_SUB32(t,1);
        } 
        f0=AE_ADD64(f0,AE_SLAI64(f0,1));
        f1=AE_ADD64(f1,AE_SLAI64(f1,1));
        AE_MUL32USEP_LL(aep,a0,AE_MOVINT32X2_FROMINT64(f0),r_frac);
        AE_MUL32USEP_LL(aep,a1,AE_MOVINT32X2_FROMINT64(f1),r_frac);
        a0=AE_SRAI72(aep,a0,32);
        a1=AE_SRAI72(aep,a1,32);
        AE_MULA32S_HL(a0,AE_MOVINT32X2_FROMINT64(f0),r_frac);
        AE_MULA32S_HL(a1,AE_MOVINT32X2_FROMINT64(f1),r_frac);
        t0=AE_TRUNCA32X2F64S(a0,a1,r_exp+2);
        AE_S32X2_IP(t0,castxcc(ae_int32x2,d),sizeof(ae_int32x2));
        c+=2;
        m+=2;
    }
    if (M&1)
    {
        ae_int32x2 t;
        ae_int64 a0;
        ae_ep aep;
        ae_int64 f0;
        t=N;
        f0=0;
        for ( n=0; n<N; n++ ) 
        {
            ae_int32x2 c00,c10;
            pC0=(const ae_int32x4*)(c+n*M);
            pC1=(const ae_int32x4*)(c+(2*N-n)*M);
            AE_L32_IP(c00,castxcc(ae_int32,pC0),sizeof(ae_int32));
            AE_L32_IP(c10,castxcc(ae_int32,pC1),sizeof(ae_int32));
            AE_MULASD32_HL_LH(f0,AE_SEL32_HH(c10,c00),t);
            t= AE_SUB32(t,1);
        } 
        f0=AE_ADD64(f0,AE_SLAI64(f0,1));
        AE_MUL32USEP_LL(aep,a0,AE_MOVINT32X2_FROMINT64(f0),r_frac);
        a0=AE_SRAI72(aep,a0,32);
        AE_MULA32S_HL(a0,AE_MOVINT32X2_FROMINT64(f0),r_frac);
        t=AE_TRUNCA32X2F64S(a0,a0,r_exp+2);
        AE_S32_L_IP(t,castxcc(ae_int32,d),sizeof(ae_int32));
    }
} /* htkdelta_32x32() */
#endif
