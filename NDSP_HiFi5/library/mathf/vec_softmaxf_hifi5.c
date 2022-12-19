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
  NatureDSP Signal Processing Library. Vector matematics
    Softmax
    Code optimized for HiFi5 core
  IntegrIT, 2006-2019
*/
#include "NatureDSP_Signal_math.h"
#include "NatureDSP_types.h"
#include "common.h"
#include "common_fpu.h"
#include "inff_tbl.h"
#include "nanf_tbl.h"
#include "expf_tbl.h"
#define SMALLER_CODESIZE 0
/*-------------------------------------------------------------------------
  Softmax
  The function computes the softmax (normalized exponential function) of 
  input data. 32-bit fixed-point functions accept inputs in Q6.25 and form 
  outputs in Q16.15 format. 

  Precision:
  32x32  32-bit inputs, 32-bit output. Accuracy: 2 LSB (see Note below)
  f      floating point input, floating point output

  Note: Accuracy of function may depend on amount of data and their 
  distribution. Given accuracy is achieved for N=2 for any pair of data 
  from input domain.

  Input:
  x[N]   input data, Q6.25 or floating point
  N      length of vectors
  Output:
  y[N]   result, Q16.15 or floating point

  Restriction:
  x,y should not overlap

-------------------------------------------------------------------------*/
#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(void,vec_softmaxf,(float32_t * y, const float32_t * x,int N))
#elif HAVE_VFPU
void vec_softmaxf(float32_t* y, const float32_t* x, int N)
{
    const xtfloatx2* restrict pX = (const xtfloatx2*)x;
    xtfloatx2* restrict pY = (xtfloatx2*)y;
    int n;
    xtfloatx2 xmax, ysum, t,t0,t1;
    if (N <= 0) return;

    /* compute maximum of x */
    {
        ae_valignx2 aX;
        xmax = minusInff.f;
        aX = AE_LA128_PP(pX);
#if SMALLER_CODESIZE
        __Pragma("no_unroll");
#endif
        for (n = 0; n < (N & ~3); n += 4)
        {
            AE_LASX2X2_IP(t0, t1, aX, castxcc(xtfloatx4, pX));
            xmax = MAX_SX2(xmax, MAX_SX2(t0, t1));
        }

#if SMALLER_CODESIZE
        __Pragma("no_unroll");
#endif
        __Pragma("loop_count min=0, max=3")
        for (n = 0; n < (N & 3); ++n)
        {
            ae_int32x2 t;
            AE_L32_IP(t, castxcc(ae_int32, pX), sizeof(ae_int32));
            xmax = MAX_SX2(xmax, AE_MOVXTFLOATX2_FROMINT32X2(t));
        }
        t = AE_SEL32_LH_SX2(xmax, xmax);
        xmax = MAX_SX2(xmax, t);
    }


    /* subtract maximum of x from input data */
    {
        ae_valignx2 aX,aY;
        pX = (const xtfloatx2*)x;
        aX=AE_LA128_PP(pX);
        aY = AE_ZALIGN128();
#if SMALLER_CODESIZE
        __Pragma("no_unroll");
#endif
        for (n = 0; n < (N & ~3); n += 4)
        {
            AE_LASX2X2_IP(t0, t1, aX, castxcc(xtfloatx4, pX));
            SUB_SX2X2(t0, t1, t0, t1, xmax, xmax);
            AE_SASX2X2_IP(t0, t1, aY, castxcc(xtfloatx4, pY));
        }
        AE_SA128POS_FP(aY, pY);
        __Pragma("no_unroll");
        for (n=0; n<(N&3); ++n)
        {
            ae_int32x2 tmp;
            AE_L32_IP(tmp, castxcc(ae_int32, pX), sizeof(xtfloat));
            t = XT_AE_MOVXTFLOATX2_FROMINT32X2(tmp);
            t = XT_SUB_SX2(t, xmax);
            tmp = XT_AE_MOVINT32X2_FROMXTFLOATX2(t);
            AE_S32_L_IP(tmp, castxcc(ae_int32, pY), sizeof(xtfloat));
        }
    }


#if SMALLER_CODESIZE
    /* compute exp() */
    vec_antilognf(y, y, N);
    /* sum results */
    pY = (xtfloatx2*)y;
    ysum = XT_CONST_S(0);
    {
        ae_valignx2 aY;
        aY = AE_LA128_PP(pY);
        xtfloatx2 t0, t1;
        xtfloatx2 ysum1 = CONST_SX2(0);
        xtfloatx2 ysum2 = CONST_SX2(0);
        xtfloatx2 ysum3 = CONST_SX2(0);
        __Pragma("no_unroll");
        for (n = 0; n < (N >> 3); n++)
        {
            AE_LASX2X2_IP(t0, t1, aY, castxcc(xtfloatx4,pY)); ADD_SX2X2(ysum, ysum1, ysum, ysum1, t0, t1);
            AE_LASX2X2_IP(t0, t1, aY, castxcc(xtfloatx4,pY)); ADD_SX2X2(ysum2, ysum3, ysum2, ysum3, t0, t1);
        }
        if (N & 4)
        {
            AE_LASX2X2_IP(t0, t1, aY, castxcc(xtfloatx4, pY)); ADD_SX2X2(ysum2, ysum3, ysum2, ysum3, t0, t1);
        }
        if (N & 2)
        {
            ae_valign aY;
            aY = AE_LA64_PP(pY);
            XT_LASX2IP(t, aY, pY); ysum1 = XT_ADD_SX2(ysum1, t);
        }
        ADD_SX2X2(ysum, ysum2, ysum, ysum2, ysum1, ysum3);
        ysum = XT_ADD_SX2(ysum, ysum2);
    }
    ysum = ADD_HL_LH_S(ysum, ysum);
    if (N & 1)
    {
        t = XT_LSI((const xtfloat*)pY, 0);
        ysum = XT_ADD_SX2(ysum, t);
    }
    ysum = XT_SEL32_HH_SX2(ysum, ysum);
#else
    const ae_int32* restrict TBL = (ae_int32*)expftbl_Q30;
    {
        const xtfloatx4* X = (xtfloatx4*)y;
        const xtfloatx4* X1 = (xtfloatx4*)y;
       
        ae_valignx2 X_va, X1_va, Y_va;
        xtfloatx2 x0, x1, x2, x3, y0, y1, z0, z1;
        ae_int32x2 tb0, tb1, tb2, tb3, tb4, tb5, tb6;
        ae_int32x2 u0, u1, n0, n1;
        ae_int32x2 e0_0, e0_1, e1_0, e1_1, e0, e1;
        ae_int64 wh0, wl0, wh1, wl1;
        ae_f32x2 f0, f1;
        xtbool2 b0_nan, b1_nan;
        xtfloatx2 ysum0, ysum1;
        int n;

        CONST_SX2X2(ysum0, ysum1, 0);
        Y_va = AE_ZALIGN128();
        X_va = AE_LA128_PP(X);
        X1_va = AE_LA128_PP(X1);
        pY = (xtfloatx2*)y;

        for (n = 0; n < (N >> 2); n++)
        {
            AE_LASX2X2_IP(x0, x1, X_va, X);
            /* scale input by 1/ln(2) and convert to Q31 */
            u0 = XT_TRUNC_SX2(x0, 24);
            u1 = XT_TRUNC_SX2(x1, 24);
            AE_MUL32X2S_HH_LL(wh0, wl0, u0, invln2_Q30);
            AE_MUL32X2S_HH_LL(wh1, wl1, u1, invln2_Q30);
            e0 = AE_TRUNCA32X2F64S(wh0, wl0, -22);
            e1 = AE_TRUNCA32X2F64S(wh1, wl1, -22);

            wh0 = AE_SRLI64(wh0, 23); wl0 = AE_SRLI64(wl0, 23);
            wh1 = AE_SRLI64(wh1, 23); wl1 = AE_SRLI64(wl1, 23);
            u0 = AE_SEL32_LL(AE_MOVINT32X2_FROMINT64(wh0), AE_MOVINT32X2_FROMINT64(wl0));
            u1 = AE_SEL32_LL(AE_MOVINT32X2_FROMINT64(wh1), AE_MOVINT32X2_FROMINT64(wl1));

            u0 = AE_AND32(u0, AE_MOVINT32X2_FROMINT32(0x7FFFFFFF));
            u1 = AE_AND32(u1, AE_MOVINT32X2_FROMINT32(0x7FFFFFFF));

            tb0 = AE_L32_I(TBL, 0 * 4);
            tb1 = AE_L32_I(TBL, 1 * 4);
            tb2 = AE_L32_I(TBL, 2 * 4);
            tb3 = AE_L32_I(TBL, 3 * 4);
            tb4 = AE_L32_I(TBL, 4 * 4);
            tb5 = AE_L32_I(TBL, 5 * 4);
            tb6 = AE_L32_I(TBL, 6 * 4);

            AE_MOVD32X4(f0, f1, tb1, tb2);
            AE_MULAF2P32X4RS(tb1, f0, u0, u1, tb0, tb0);
            AE_MULAF2P32X4RS(tb2, f1, u0, u1, tb1, f0);
            AE_MOVD32X4(n0, n1, tb3, tb4);
            AE_MULAF2P32X4RS(tb3, n0, u0, u1, tb2, f1);
            AE_MULAF2P32X4RS(tb4, n1, u0, u1, tb3, n0);
            AE_MOVD32X4(f0, f1, tb5, tb6);
            AE_MULAF2P32X4RS(tb5, f0, u0, u1, tb4, n1);
            AE_MULAF2P32X4RS(tb6, f1, u0, u1, tb5, f0);

            x0 = XT_FLOAT_SX2(tb6, 30);
            x1 = XT_FLOAT_SX2(f1, 30);

            e0_0 = AE_SRAI32(e0, 1);
            e0_1 = AE_SRAI32(e1, 1);
            e1_0 = AE_SUB32S(e0, e0_0);
            e1_1 = AE_SUB32S(e1, e0_1);
            y0 = FLOATEXP_SX2(AE_MOVINT8X8_FROMINT32X2(e0_0));
            y1 = FLOATEXP_SX2(AE_MOVINT8X8_FROMINT32X2(e0_1));
            z0 = FLOATEXP_SX2(AE_MOVINT8X8_FROMINT32X2(e1_0));
            z1 = FLOATEXP_SX2(AE_MOVINT8X8_FROMINT32X2(e1_1));

            AE_LASX2X2_IP(x2, x3, X1_va, X1);
            b0_nan = XT_UN_SX2(x2, x2);
            b1_nan = XT_UN_SX2(x3, x3);
            XT_MOVT_SX2(x0, x2, b0_nan);
            XT_MOVT_SX2(x1, x3, b1_nan);

            MUL_SX2X2(y0, y1, y0, y1, z0, z1);
            MUL_SX2X2(y0, y1, x0, x1, y0, y1);

            AE_SASX2X2_IP(y0, y1, Y_va, castxcc(xtfloatx4,pY));
            ADD_SX2X2(ysum0, ysum1, ysum0, ysum1, y0, y1);
        }
        AE_SA128POS_FP(Y_va, pY);
        ysum = ADD_SX2(ysum0, ysum1);
        ysum = ADD_HL_LH_S(ysum, ysum);
        if ((N&3) > 0)
        {
            xtfloatx2 x0, x1, e0, e1;
            int32_t ALIGN(32) scratch[4];
            xtfloatx4* pScr;
            pScr = (xtfloatx4*)scratch;
            AE_SSX2X2_IP(0, 0, pScr, 0);
            __Pragma("no_unroll");
            for (n = 0; n < (N&3); n++)
            {
                xtfloat t;
                AE_LSIP(t, castxcc(xtfloat, X), sizeof(xtfloat));
                AE_SSIP(t, castxcc(xtfloat, pScr), sizeof(xtfloat));
            }
            pScr = (xtfloatx4*)scratch;
            AE_LSX2X2_IP(x0, x1, pScr, 0);

            /* scale input by 1/ln(2) and convert to Q31 */
            u0 = XT_TRUNC_SX2(x0, 24);
            u1 = XT_TRUNC_SX2(x1, 24);
            AE_MUL32X2S_HH_LL(wh0, wl0, u0, invln2_Q30);
            AE_MUL32X2S_HH_LL(wh1, wl1, u1, invln2_Q30);
            e0 = AE_TRUNCA32X2F64S(wh0, wl0, -22);
            e1 = AE_TRUNCA32X2F64S(wh1, wl1, -22);

            wh0 = AE_SRLI64(wh0, 23);
            wl0 = AE_SRLI64(wl0, 23);
            wh1 = AE_SRLI64(wh1, 23);
            wl1 = AE_SRLI64(wl1, 23);

            u0 = AE_SEL32_LL(AE_MOVINT32X2_FROMINT64(wh0), AE_MOVINT32X2_FROMINT64(wl0));
            u1 = AE_SEL32_LL(AE_MOVINT32X2_FROMINT64(wh1), AE_MOVINT32X2_FROMINT64(wl1));

            u0 = AE_AND32(u0, AE_MOVINT32X2_FROMINT32(0x7FFFFFFF));
            u1 = AE_AND32(u1, AE_MOVINT32X2_FROMINT32(0x7FFFFFFF));

            tb0 = AE_L32_I(TBL, 0 * 4);
            tb1 = AE_L32_I(TBL, 1 * 4);
            tb2 = AE_L32_I(TBL, 2 * 4);
            tb3 = AE_L32_I(TBL, 3 * 4);
            tb4 = AE_L32_I(TBL, 4 * 4);
            tb5 = AE_L32_I(TBL, 5 * 4);
            tb6 = AE_L32_I(TBL, 6 * 4);

            AE_MOVD32X4(f0, f1, tb1, tb2);
            AE_MULAF2P32X4RS(tb1, f0, u0, u1, tb0, tb0);
            AE_MULAF2P32X4RS(tb2, f1, u0, u1, tb1, f0);
            AE_MOVD32X4(n0, n1, tb3, tb4);
            AE_MULAF2P32X4RS(tb3, n0, u0, u1, tb2, f1);
            AE_MULAF2P32X4RS(tb4, n1, u0, u1, tb3, n0);
            AE_MOVD32X4(f0, f1, tb5, tb6);
            AE_MULAF2P32X4RS(tb5, f0, u0, u1, tb4, n1);
            AE_MULAF2P32X4RS(tb6, f1, u0, u1, tb5, f0);

            x0 = XT_FLOAT_SX2(tb6, 30);
            x1 = XT_FLOAT_SX2(f1, 30);

            e0_0 = AE_SRAI32(e0, 1);
            e0_1 = AE_SRAI32(e1, 1);
            e1_0 = AE_SUB32(e0, e0_0);
            e1_1 = AE_SUB32(e1, e0_1);
            y0 = FLOATEXP_SX2(AE_MOVINT8X8_FROMINT32X2(e0_0));
            y1 = FLOATEXP_SX2(AE_MOVINT8X8_FROMINT32X2(e0_1));
            z0 = FLOATEXP_SX2(AE_MOVINT8X8_FROMINT32X2(e1_0));
            z1 = FLOATEXP_SX2(AE_MOVINT8X8_FROMINT32X2(e1_1));

            AE_LASX2X2_IP(x2, x3, X1_va, X1);
            b0_nan = XT_UN_SX2(x2, x2);
            b1_nan = XT_UN_SX2(x3, x3);
            XT_MOVT_SX2(x0, x2, b0_nan);
            XT_MOVT_SX2(x1, x3, b1_nan);

            MUL_SX2X2(y0, y1, y0, y1, z0, z1);
            MUL_SX2X2(y0, y1, x0, x1, y0, y1);
            AE_SSX2X2_I(y0, y1, pScr, 0 * sizeof(xtfloatx4));
            __Pragma("no_unroll")
            for (n = 0; n < (N&3); n++)
            {
                xtfloat t;
                AE_LSIP(t, castxcc(xtfloat, pScr), sizeof(xtfloat));
                ysum = XT_ADD_SX2(ysum, t);
                AE_SSIP(t, castxcc(xtfloat, pY), sizeof(xtfloat));
            }
        }
        ysum = XT_SEL32_HH_SX2(ysum, ysum);
    } 
#endif



    /* normalize output */
    ysum = XT_RECIP_S(ysum);

    {
        ae_valignx2 aY;
        xtfloatx2 t0,t1;
        int M = N;
        pX = (xtfloatx2*)y;
        pY = (xtfloatx2*)y;


        if (M > 3)
        {
            for (n = 0; n < (int)((uintptr_t)pX & 15); n += 4)
            {
                ae_int32x2 tmp;
                AE_L32_IP(tmp, castxcc(ae_int32, pX), sizeof(xtfloat));
                t = XT_AE_MOVXTFLOATX2_FROMINT32X2(tmp);
                t = XT_MUL_SX2(t, ysum);
                tmp = XT_AE_MOVINT32X2_FROMXTFLOATX2(t);
                AE_S32_L_IP(tmp, castxcc(ae_int32, pY), sizeof(xtfloat));
                M--;
            }

            aY = AE_ZALIGN128(); 
            for (n = 0; n < (M & ~3); n += 4)
            {
                AE_LSX2X2_IP(t0, t1, castxcc(xtfloatx4, pX), sizeof(xtfloatx4));
                MULQ_S(t0, t1, t0, t1, ysum);
                AE_SASX2X2_IP(t0, t1, aY, castxcc(xtfloatx4, pY));
            }
            AE_SA128POS_FP(aY, pY);
        }
        
        __Pragma("no_unroll");
        __Pragma("loop_count min=0, max=3");
        for (n=0; n<(M&3); ++n)
        {
            xtfloat t;
            AE_LSIP(t, castxcc(xtfloat, pX), sizeof(xtfloat));
            t = MUL_S(t, ysum);
            AE_SSIP(t, castxcc(xtfloat, pY), sizeof(xtfloat));
        }
    }
}
#else
// code for scalar FPU
void vec_softmaxf    (float32_t * y, const float32_t * x,int N)
{
    const xtfloat* restrict pX=(const xtfloat*)x;
          xtfloat* restrict pY=(      xtfloat*)y;
    int n;
    xtfloat xmax,ysum;
    if (N<0) return;
    /* compute maximum of x */
    xmax=minusInff.f;
    for (n=0; n<N; n++)
    {
        xtfloat t;
        XT_LSIP(t,pX,sizeof(xtfloat));
        XT_MOVT_S(xmax,t,XT_OLT_S(xmax,t));
    }
    /* subtract maximum of x from input data */
    pX=(const xtfloat*)x;
    for (n=0; n<N; n++)
    {
        xtfloat t;
        XT_LSIP(t,pX,sizeof(xtfloat));
        t=XT_SUB_S(t,xmax);
        XT_SSIP(t,pY,sizeof(xtfloat));
    }
    /* compute exp() */
    vec_antilognf(y,y,N);
    /* sum results */
    pY=(xtfloat*)y;
    ysum=XT_CONST_S(0);
    for (n=0; n<N; n++)
    {
        xtfloat t;
        XT_LSIP(t,pY,sizeof(xtfloat));
        ysum=XT_ADD_S(ysum,t);
    }
    /* normalize output */
    ysum=XT_RECIP_S(ysum);
    __Pragma("no_reorder")
    pX=(xtfloat*)y;
    pY=(xtfloat*)y;
    for (n=0; n<N; n++) 
    {
        xtfloat t;
        XT_LSIP(t,pX,sizeof(xtfloat));
        t=XT_MUL_S(t,ysum);
        XT_SSIP(t,pY,sizeof(xtfloat));
    }
} /* vec_softmaxf() */
#endif
