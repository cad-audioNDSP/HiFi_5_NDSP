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
    Floating-point variant
    C code optimized for HiFi5 core with VFPU/SFPU option
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "common_fpu.h"
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

#if (!HAVE_VFPU && !HAVE_FPU) 
DISCARD_FUN(void, htkdeltaf, (float32_t * d, const float32_t * c, int M, int N));
#elif HAVE_VFPU
void htkdeltaf(float32_t * d, const float32_t * c, int M, int N)
{
    const xtfloatx2 * restrict pC0;
    const xtfloatx2 * restrict pC1;
    const xtfloatx2 * restrict pC3;
    const xtfloatx2 * restrict pC4;
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
    float32_t r;
    int m, n;
    if (N<=0) 
    {
        for ( m=0; m<M; m++ ) 
        {
            d[m] = XT_CONST_S(0);
        } /* m */
        return;
    } /* (N<=0) */
    if (N==2)
    {
        // specially optimized case for N==2
        ae_valignx2 aC1,aC3;
        pC0=(const xtfloatx2*)(c    );
        pC1=(const xtfloatx2*)(c+1*M);
        pC3=(const xtfloatx2*)(c+3*M);
        pC4=(const xtfloatx2*)(c+4*M);
        aC1=AE_LA128_PP(pC1);
        aC3=AE_LA128_PP(pC3);
        r = 1.0f/10.f;
        for ( m=0; m<(M>>2); m++) 
        {
            xtfloatx2 f00,f01,f10,f11;
            xtfloatx2 c00,c01,c40,c41;
            AE_LSX2X2_IP (c00,c01,    castxcc(xtfloatx4,pC0),sizeof(xtfloatx4));
            AE_LASX2X2_IP(f00,f01,aC1,castxcc(xtfloatx4,pC1));
            AE_LASX2X2_IP(f10,f11,aC3,castxcc(xtfloatx4,pC3));
            AE_LSX2X2_IP (c40,c41,    castxcc(xtfloatx4,pC4),sizeof(xtfloatx4));
            MADDQ_S(f00,f01,c00,c01,XT_CONST_S(2));
            MADDQ_S(f10,f11,c40,c41,XT_CONST_S(2));
            f10=XT_SUB_SX2(f10,f00);
            f11=XT_SUB_SX2(f11,f01);
            MULQ_S(f10,f11,f10,f11,r);
            AE_SSX2X2_IP(f10,f11,castxcc(xtfloatx4,d),sizeof(xtfloatx4));
        }
        if (M&2)
        {
            ae_valign aC1,aC3;
            xtfloatx2 f00,f10;
            xtfloatx2 c00,c40;
            aC1=AE_LA64_PP(pC1);
            aC3=AE_LA64_PP(pC3);

            XT_LSX2IP (c00,    castxcc(xtfloatx2,pC0),sizeof(xtfloatx2));
            XT_LASX2IP(f00,aC1,castxcc(xtfloatx2,pC1));
            XT_LASX2IP(f10,aC3,castxcc(xtfloatx2,pC3));
            XT_LSX2IP (c40,    castxcc(xtfloatx2,pC4),sizeof(xtfloatx2));
            XT_MADD_SX2(f00,c00,XT_CONST_S(2));
            XT_MADD_SX2(f10,c40,XT_CONST_S(2));
            f10=XT_SUB_SX2(f10,f00);
            f10=XT_MUL_SX2(f10,r);
            XT_SSX2IP(f10,castxcc(xtfloatx2,d),sizeof(xtfloatx2));
        }
        if (M&1)
        {
            xtfloat  f0,f1;
            xtfloat c0,c4;
            XT_LSIP(c0,castxcc(xtfloat,pC0),sizeof(xtfloat));
            XT_LSIP(c4,castxcc(xtfloat,pC4),sizeof(xtfloat));
            XT_LSIP(f0,castxcc(xtfloat,pC1),sizeof(xtfloat));
            XT_LSIP(f1,castxcc(xtfloat,pC3),sizeof(xtfloat));
            XT_MADD_S(f0,2.f,c0);
            XT_MADD_S(f1,2.f,c4);
            *d++ = XT_MUL_S(XT_SUB_S(f1,f0),r);
        }

        return;
    }
    /* s2 <- sum((1:N).^2);
     r <- 1/(2*s2)
     We use the following formula for the sum of first N natural numbers squared:
       s2 == N*(N+1)*(2*N+1)/6
     See Ken Ward's Mathematics Pages:
      trans4mind.com/personal_development/mathematics/series/sumNaturalSquares.htm */
    r = XT_MUL_S(3.f,XT_RECIP_S((float32_t)(N*(N+1)*(2*N+1))));
    for ( m=0; m<(M&~3); m+=4) 
    {
        xtfloatx2 f02,f13,f46,f57,theta;
        f46=f57=f02=f13=XT_CONST_S(0);
        theta=N;
        for ( n=0; n<N; n++ ) 
        {
            ae_valignx2 aC1,aC0;
            xtfloatx2 c02,c13,c46,c57;
            pC0=(const xtfloatx2*)&c[n*M];
            pC1=(const xtfloatx2*)&c[(2*N-n)*M];
            aC0=AE_LA128_PP(pC0);
            aC1=AE_LA128_PP(pC1);
            AE_LASX2X2_IP(c02,c46,aC0,castxcc(xtfloatx4,pC0));
            AE_LASX2X2_IP(c13,c57,aC1,castxcc(xtfloatx4,pC1));
            MADDQ_S(f02,f46,c02,c46,theta);
            MADDQ_S(f13,f57,c13,c57,theta);
            theta=XT_SUB_SX2(theta,XT_CONST_S(1));
        }
        f13=XT_SUB_SX2(f13,f02);
        f57=XT_SUB_SX2(f57,f46);
        MULQ_S(f13,f57,f13,f57,r);
        AE_SSX2X2_IP(f13,f57,castxcc(xtfloatx4,d),sizeof(xtfloatx4));
        c+=4;
    } 
    if (M&2)
    {
        xtfloatx2 f02,f13,theta;
        f02=f13=XT_CONST_S(0);
        theta=N;
        for ( n=0; n<N; n++ ) 
        {
            ae_valign aC1,aC0;
            xtfloatx2 c02,c13;
            pC0=(const xtfloatx2*)&c[n*M];
            pC1=(const xtfloatx2*)&c[(2*N-n)*M];
            aC0=AE_LA64_PP(pC0);
            aC1=AE_LA64_PP(pC1);
            XT_LASX2IP(c02,aC0,pC0);
            XT_LASX2IP(c13,aC1,pC1);
            XT_MADD_SX2(f02,theta,c02);
            XT_MADD_SX2(f13,theta,c13);
            theta=XT_SUB_SX2(theta,XT_CONST_S(1));
        }
        XT_SSX2IP(XT_MUL_SX2(XT_SUB_SX2(f13,f02),r),castxcc(xtfloatx2,d),sizeof(xtfloatx2));
        c+=2;
    } 
    if (M&1)
    {
        float32_t f0,f1,theta;
        f0=f1=XT_CONST_S(0);
        theta=N;
        for ( n=0; n<N; n++ ) 
        {
            XT_MADD_S(f0,theta,c[n*M]        );
            XT_MADD_S(f1,theta,c[(2*N-n)*M]  );
            theta=XT_SUB(theta,XT_CONST_S(1));
        }
        *d++ = XT_MUL_S(f1-f0,r);
    }
} /* htkdeltaf() */
#else
// scalar FPU
void htkdeltaf(float32_t * d, const float32_t * c, int M, int N)
{
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
    float32_t r;
    int m, n;
    if (N<=0) 
    {
        for ( m=0; m<M; m++ ) 
        {
            d[m] = XT_CONST_S(0);
        } /* m */
        return;
    } /* (N<=0) */
    /* s2 <- sum((1:N).^2);
     r <- 1/(2*s2)
     We use the following formula for the sum of first N natural numbers squared:
       s2 == N*(N+1)*(2*N+1)/6
     See Ken Ward's Mathematics Pages:
       trans4mind.com/personal_development/mathematics/series/sumNaturalSquares.htm */
    r = XT_MUL_S(3.f,XT_RECIP_S((float32_t)(N*(N+1)*(2*N+1))));
    for ( m=0; m<(M&~1); m+=2) 
    {
        float32_t f0,f1,f2,f3,theta;
        f0=f1=f2=f3=XT_CONST_S(0);
        theta=N;
        for ( n=0; n<N; n++ ) 
        {
            XT_MADD_S(f0,theta,c[n*M+m]        );
            XT_MADD_S(f1,theta,c[(2*N-n)*M+m]  );
            XT_MADD_S(f2,theta,c[n*M+m+1]      );
            XT_MADD_S(f3,theta,c[(2*N-n)*M+m+1]);
            theta=XT_SUB(theta,XT_CONST_S(1));
        }
        *d++ = XT_MUL_S(f1-f0,r);
        *d++ = XT_MUL_S(f3-f2,r);
    } 
    if (M&1)
    {
        float32_t f0,f1,theta;
        f0=f1=XT_CONST_S(0);
        theta=N;
        for ( n=0; n<N; n++ ) 
        {
            XT_MADD_S(f0,theta,c[n*M+m]        );
            XT_MADD_S(f1,theta,c[(2*N-n)*M+m]  );
            theta=XT_SUB(theta,XT_CONST_S(1));
        }
        *d++ = XT_MUL_S(f1-f0,r);
    }
} 
#endif
