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
  NatureDSP Signal Processing Library. FIR part
  2D Convolution  
  IntegrIT, 2006-2020
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
#include "common_fpu.h"


/*-------------------------------------------------------------------------
  2D convolution
  Functions compute the two-dimensional convolution of input matrix x[M][N]
  and y[P][Q] and store the result in matrix z[M+P-1][N+Q-1]
  Additional parameter rsh allows to contro l fixed point representation of 
  output data.
  Two versions of functions available: 
  - generic version with _gen_ suffix. 
    These functions work with arbitrary arguments.
  - fast version with no _gen_ suffix. 
    These functions expose some additional restrictions on argument


  Precision: 
  8x8      8-bit coefficients, 8-bit data, 8-bit output, Q7
  8x16     8-bit coefficients Q7, 16-bit data, 16-bit output, Q15
  16x16    16-bit coefficients, 16-bit data, 16-bit output, Q15
  f        single precision floating point data
  fp16     half precision floating point data

  Input:
  x[M][N]   input data Q15, Q7, floating point
  y[P][Q]   input data Q15, Q7, floating point
  M         number of rows in the matrix x
  N         number of columns in the matrix x
  P         number of rows in the matrix y
  Q         number of columns in the matrix y
  rsh       additional right shift (for fixed point API only)

  Output:
  z	[M+P-1][N+Q-1] output data, Q(7-rsh), Q(15-rsh)

  Temporary:
  pScr     scratch data. Should have size at least as requested by 
           corresponding scratch allocation function

  Restrictions:
  For regular routines:
  x,y,z        should not overlap
  pScr         aligned on a 16-bytes boundary
  P, Q	       >0

  For fast routines:
  x,y,z        should not overlap
  x,y,z,pScr   aligned on a 16-bytes boundary
  P, Q	       >0 and multiplies of 8
-------------------------------------------------------------------------*/

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void,conv2d_gen_11x7f,(void* pScr, float32_t *z, const float32_t * x, const float32_t * y, int P, int Q))
size_t conv2d_gen_11x7f_getScratchSize (int P, int Q) 
{ 
    (void)P,(void)Q;
    return 0;
}
#elif HAVE_VFPU

void conv2d_gen_11x7f     (void* pScr, float32_t *z, const float32_t * x, const float32_t * y, int P, int Q)
{
#define M 11
#define N 7
    const xtfloatx4* restrict pY;
    const xtfloatx4* restrict pY0;
    const xtfloatx4* restrict pY1;
    const xtfloatx4* restrict pT_read;
    xtfloatx4* restrict pT;

    xtfloatx4* restrict pW;
    xtfloatx4* restrict pZ0;
    const xtfloat* pX = (const xtfloat*)x;

    ae_valignx2 al0, al1, alZ;

    int i, j, mmin, mmax, xstart, ystart, m;
    xtfloatx2 Y01, Y12, Y23, Y34, Y45, Y56, Y67, Y78, Y89, Y910;
    xtfloatx2 W0, W1, W2, W3, W4, W5, W6;
    xtfloatx2 dummy;
    xtfloatx2 S0, S1, S2, S3;

    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0);
    NASSERT(Q >= 0);
    if (P <= 0 || Q <= 0) return;


    alZ = AE_ZALIGN128();

    pW = (xtfloatx4*)pScr + (M-1) * 4;

    // preload X and save in scratch
    {
        xtfloatx2 w0, w1, w2, w3, w4, w5, w6;
        ae_int32x2 v_i;

        for (m = 0; m < M; ++m)
        {
            AE_L32_IP(v_i, castxcc(const ae_int32, pX), sizeof(xtfloat));
            w0 = AE_MOVXTFLOATX2_FROMINT32X2(v_i);
            AE_L32_IP(v_i, castxcc(const ae_int32, pX), sizeof(xtfloat));
            w1 = AE_MOVXTFLOATX2_FROMINT32X2(v_i);
            AE_L32_IP(v_i, castxcc(const ae_int32, pX), sizeof(xtfloat));
            w2 = AE_MOVXTFLOATX2_FROMINT32X2(v_i);
            AE_L32_IP(v_i, castxcc(const ae_int32, pX), sizeof(xtfloat));
            w3 = AE_MOVXTFLOATX2_FROMINT32X2(v_i);
            AE_L32_IP(v_i, castxcc(const ae_int32, pX), sizeof(xtfloat));
            w4 = AE_MOVXTFLOATX2_FROMINT32X2(v_i);
            AE_L32_IP(v_i, castxcc(const ae_int32, pX), sizeof(xtfloat));
            w5 = AE_MOVXTFLOATX2_FROMINT32X2(v_i);
            AE_L32_IP(v_i, castxcc(const ae_int32, pX), sizeof(xtfloat));
            w6 = AE_MOVXTFLOATX2_FROMINT32X2(v_i);

            AE_SSX2X2_IP(w6, w5, pW, sizeof(xtfloatx4));
            AE_SSX2X2_IP(w4, w3, pW, sizeof(xtfloatx4));
            AE_SSX2X2_IP(w2, w1, pW, sizeof(xtfloatx4));
            AE_SSX2X2_XP(w0, CONST_SX2(0), pW, -7 * (int)sizeof(xtfloatx4));
        }
    }

    __Pragma("loop_count min=1");
    for (i = 0; i < P + M - 1; i += 1)
    {
        mmin = i < M - 1 ? M - 1 - i : 0;
        mmax = i >= P ? M - 1 + P - i : M;
        xstart = i < M - 1 ? (mmin) : 0;
        ystart = i < M - 1 ? 0 : i - (M - 1);

        pW = (xtfloatx4*)pScr + 4 * xstart;
        pY0 = (xtfloatx4*)(y + ystart * Q);

        pZ0 = (xtfloatx4*)(z + i * (Q + N - 1));

        if (mmax - mmin == 1)
        {
            Y01 = CONST_SX2(0);
            Y12 = CONST_SX2(0);
            Y23 = CONST_SX2(0);
            Y34 = CONST_SX2(0);
            Y45 = CONST_SX2(0);

            al0 = AE_LA128_PP(pY0);
            AE_LASX2X2_IP(Y67, Y89, al0, pY0);
            Y56 = AE_SEL32_LH_SX2(Y45,Y67);
            pY0 = (xtfloatx4*)((xtfloat*)(pY0)-4);

            pY1 = (xtfloatx4*)((xtfloat*)pY0 + 1);
            al0 = AE_LA128_PP(pY0);
            al1 = AE_LA128_PP(pY1);

            AE_LSX2X2_IP(W0, W1, pW, sizeof(xtfloatx4));
            AE_LSX2X2_IP(W2, W3, pW, sizeof(xtfloatx4));
            AE_LSX2X2_IP(W4, W5, pW, sizeof(xtfloatx4));
            AE_LSX2X2_IP(W6, dummy, pW, sizeof(xtfloatx4));

            for (j=0; j < (Q & ~3); j += 4)
            {
                AE_LASX2X2_IP(Y67, Y89, al0, pY0);
                AE_LASX2X2_IP(Y78, Y910, al1, pY1);

                MULQ_S (S0, S1, Y01, Y23, W0);
                MADDQ_S(S0, S1, Y12, Y34, W1);
                MADDQ_S(S0, S1, Y23, Y45, W2);
                MADDQ_S(S0, S1, Y34, Y56, W3);
                MADDQ_S(S0, S1, Y45, Y67, W4);
                MADDQ_S(S0, S1, Y56, Y78, W5);
                MADDQ_S(S0, S1, Y67, Y89, W6);

                AE_SASX2X2_IP(S0, S1, alZ, pZ0);

                Y01 = Y45; Y12 = Y56;
                Y23 = Y67; Y34 = Y78;
                Y45 = Y89; Y56 = Y910;
            }
            AE_SA128POS_FP(alZ, pZ0);
            if (Q & 3)
            {
                AE_LASX2X2_IP(Y67, Y89, al0, pY0);
                AE_LASX2X2_IP(Y78, Y910, al1, pY1);

                MULQ_S(S0, S1, Y01, Y23, W0);
                MADDQ_S(S0, S1, Y12, Y34, W1);
                MADDQ_S(S0, S1, Y23, Y45, W2);
                MADDQ_S(S0, S1, Y34, Y56, W3);
                MADDQ_S(S0, S1, Y45, Y67, W4);
                MADDQ_S(S0, S1, Y56, Y78, W5);
                MADDQ_S(S0, S1, Y67, Y89, W6);

                for (j= (Q & ~3); j < Q; ++j)
                {
                    xtfloatx2 tmp;

                    tmp = AE_SEL32_LH_SX2(S0, S1); S1 = AE_SEL32_LH_SX2(S1, S0); S0 = tmp;

                    Y01 = Y12; Y12 = Y23; Y23 = Y34; Y34 = Y45; Y45 = Y56; Y56 = Y67; Y67 = Y78; Y78 = Y89; Y89 = Y910;

                    AE_S32_L_IP(AE_MOVINT32X2_FROMXTFLOATX2(S1), castxcc(ae_int32, pZ0), sizeof(ae_int32));
                }
            }
            //tail 4 el
            {
                Y910 = CONJC_SX2(0);
                Y89 = CONST_SX2(0);
                Y78 = CONST_SX2(0);
                Y67 = CONST_SX2(0); 
                Y56 = AE_SEL32_HL_SX2(Y56, CONST_SX2(0));

                MULQ_S(S0, S1, Y01, Y23, W0);
                MADDQ_S(S0, S1, Y12, Y34, W1);
                MADDQ_S(S0, S1, Y23, Y45, W2);
                MADDQ_S(S0, S1, Y34, Y56, W3);
                MADDQ_S(S0, S1, Y45, Y67, W4);
                MADDQ_S(S0, S1, Y56, Y78, W5);
                MADDQ_S(S0, S1, Y67, Y89, W6);

                AE_SASX2X2_IP(S0, S1, alZ, pZ0);

                Y01 = Y45; Y12 = Y56;
                Y23 = Y67; Y34 = Y78;
                Y45 = Y89; Y56 = Y910;
            }
            AE_SA128POS_FP(alZ, pZ0);
            //tail 2 el
            {
                xtfloatx2 tmp;

                MULQ_S(S0, S1, Y01, Y23, W0);
                MADDQ_S(S0, S1, Y12, Y34, W1);
                MADDQ_S(S0, S1, Y23, Y45, W2);
                MADDQ_S(S0, S1, Y34, Y56, W3);
                MADDQ_S(S0, S1, Y45, Y67, W4);
                MADDQ_S(S0, S1, Y56, Y78, W5);
                MADDQ_S(S0, S1, Y67, Y89, W6);

                tmp = AE_SEL32_LH_SX2(S0, S1); S1 = AE_SEL32_LH_SX2(S1, S0); S0 = tmp;
                AE_S32_L_IP(AE_MOVINT32X2_FROMXTFLOATX2(S1), castxcc(ae_int32, pZ0), sizeof(ae_int32));

                tmp = AE_SEL32_LH_SX2(S0, S1); S1 = AE_SEL32_LH_SX2(S1, S0); S0 = tmp;
                AE_S32_L_IP(AE_MOVINT32X2_FROMXTFLOATX2(S1), castxcc(ae_int32, pZ0), sizeof(ae_int32));
            }
        }
        else
        {
            //first iteration
            pY = pY0;
            {
                pT = (xtfloatx4*)pScr + 4 * 11;

                Y01 = CONST_SX2(0);
                Y12 = CONST_SX2(0);
                Y23 = CONST_SX2(0);
                Y34 = CONST_SX2(0);
                Y45 = CONST_SX2(0);

                al0 = AE_LA128_PP(pY0);
                AE_LASX2X2_IP(Y67, Y89, al0, pY0);
                Y56 = AE_SEL32_LH_SX2(Y45, Y67);
                pY0 = (xtfloatx4*)((xtfloat*)(pY0)-4);

                pY1 = (xtfloatx4*)((xtfloat*)pY0 + 1);
                al0 = AE_LA128_PP(pY0);
                al1 = AE_LA128_PP(pY1);

                AE_LSX2X2_IP(W0, W1, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W2, W3, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W4, W5, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W6, dummy, pW, sizeof(xtfloatx4));

                for (j=0; j < (Q & ~3); j += 4)
                {
                    AE_LASX2X2_IP(Y67, Y89, al0, pY0);
                    AE_LASX2X2_IP(Y78, Y910, al1, pY1);

                    MULQ_S(S0, S1, Y01, Y23, W0);
                    MULQ_S(S2, S3, Y12, Y34, W1);
                    MADDQ_S(S0, S1, Y23, Y45, W2);
                    MADDQ_S(S2, S3, Y34, Y56, W3);
                    MADDQ_S(S0, S1, Y45, Y67, W4);
                    MADDQ_S(S2, S3, Y56, Y78, W5);
                    MADDQ_S(S0, S1, Y67, Y89, W6);
                    ADD_SX2X2(S0, S1, S0, S1, S2, S3);

                    AE_SSX2X2_IP(S0, S1, pT, sizeof(xtfloatx4));

                    Y01 = Y45; Y12 = Y56;
                    Y23 = Y67; Y34 = Y78;
                    Y45 = Y89; Y56 = Y910;
                }
                // tail
                if (Q & 3)
                {
                    AE_LASX2X2_IP(Y67, Y89, al0, pY0);
                    AE_LASX2X2_IP(Y78, Y910, al1, pY1);

                    MULQ_S(S0, S1, Y01, Y23, W0);
                    MULQ_S(S2, S3, Y12, Y34, W1);
                    MADDQ_S(S0, S1, Y23, Y45, W2);
                    MADDQ_S(S2, S3, Y34, Y56, W3);
                    MADDQ_S(S0, S1, Y45, Y67, W4);
                    MADDQ_S(S2, S3, Y56, Y78, W5);
                    MADDQ_S(S0, S1, Y67, Y89, W6);
                    ADD_SX2X2(S0, S1, S0, S1, S2, S3);

                    AE_SSX2X2_IP(S0, S1, pT, sizeof(xtfloatx4));
                    for (j = Q - (Q & 3); j < Q; ++j)
                    {
                        Y01 = Y12; Y12 = Y23; Y23 = Y34; Y34 = Y45; Y45 = Y56; Y56 = Y67; Y67 = Y78; Y78 = Y89; Y89 = Y910;
                    }
                }
                //last N-1 elements
                {
                    Y910 = CONJC_SX2(0);
                    Y89 = CONST_SX2(0);
                    Y78 = CONST_SX2(0);
                    Y67 = CONST_SX2(0);
                    Y56 = AE_SEL32_HL_SX2(Y56, CONST_SX2(0));

                    MULQ_S(S0, S1, Y01, Y23, W0);
                    MULQ_S(S2, S3, Y12, Y34, W1);
                    MADDQ_S(S0, S1, Y23, Y45, W2);
                    MADDQ_S(S2, S3, Y34, Y56, W3);
                    MADDQ_S(S0, S1, Y45, Y67, W4);
                    MADDQ_S(S2, S3, Y56, Y78, W5);
                    MADDQ_S(S0, S1, Y67, Y89, W6);
                    ADD_SX2X2(S0, S1, S0, S1, S2, S3);

                    AE_SSX2X2_IP(S0, S1, pT, sizeof(xtfloatx4));

                    Y01 = Y45; Y12 = Y56;
                    Y23 = Y67; Y34 = Y78;
                    Y45 = Y89; Y56 = Y910;

                    MULQ_S(S0, S1, Y01, Y23, W0);
                    MULQ_S(S2, S3, Y12, Y34, W1);
                    MADDQ_S(S0, S1, Y23, Y45, W2);
                    MADDQ_S(S2, S3, Y34, Y56, W3);
                    MADDQ_S(S0, S1, Y45, Y67, W4);
                    MADDQ_S(S2, S3, Y56, Y78, W5);
                    MADDQ_S(S0, S1, Y67, Y89, W6);
                    ADD_SX2X2(S0, S1, S0, S1, S2, S3);

                    AE_SSX2X2_IP(S0, S1, pT, sizeof(xtfloatx4));
                }
            }
            pY = (xtfloatx4*)((xtfloat*)pY + Q);
            pY0 = pY;
            //next iterations 
            __Pragma("loop_count min=0, max=9");
            for (m = mmin + 1; m < mmax - 1; ++m)
            {
                pT_read = pT = (xtfloatx4*)pScr + 4 * 11;

                Y01 = CONST_SX2(0);
                Y12 = CONST_SX2(0);
                Y23 = CONST_SX2(0);
                Y34 = CONST_SX2(0);
                Y45 = CONST_SX2(0);

                al0 = AE_LA128_PP(pY0);
                AE_LASX2X2_IP(Y67, Y89, al0, pY0);
                Y56 = AE_SEL32_LH_SX2(Y45, Y67);
                pY0 = (xtfloatx4*)((xtfloat*)(pY0)-4);

                pY1 = (xtfloatx4*)((xtfloat*)pY0 + 1);
                al0 = AE_LA128_PP(pY0);
                al1 = AE_LA128_PP(pY1);

                AE_LSX2X2_IP(W0, W1, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W2, W3, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W4, W5, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W6, dummy, pW, sizeof(xtfloatx4));

                for (j=0; j < (Q & ~3); j += 4)
                {
                    AE_LSX2X2_IP(S0, S1, pT_read, sizeof(xtfloatx4));

                    AE_LASX2X2_IP(Y67, Y89, al0, pY0);
                    AE_LASX2X2_IP(Y78, Y910, al1, pY1);

                    MULQ_S(S2, S3, Y01, Y23, W0);
                    MADDQ_S(S0, S1, Y12, Y34, W1);
                    MADDQ_S(S2, S3, Y23, Y45, W2);
                    MADDQ_S(S0, S1, Y34, Y56, W3);
                    MADDQ_S(S2, S3, Y45, Y67, W4);
                    MADDQ_S(S0, S1, Y56, Y78, W5);
                    MADDQ_S(S2, S3, Y67, Y89, W6);
                    ADD_SX2X2(S0, S1, S0, S1, S2, S3);

                    AE_SSX2X2_IP(S0, S1, pT, sizeof(xtfloatx4));

                    Y01 = Y45; Y12 = Y56;
                    Y23 = Y67; Y34 = Y78;
                    Y45 = Y89; Y56 = Y910;
                }
                // tail
                if (Q & 3)
                {
                    AE_LSX2X2_IP(S0, S1, pT_read, sizeof(xtfloatx4));

                    AE_LASX2X2_IP(Y67, Y89, al0, pY0);
                    AE_LASX2X2_IP(Y78, Y910, al1, pY1);

                    MULQ_S(S2, S3, Y01, Y23, W0);
                    MADDQ_S(S0, S1, Y12, Y34, W1);
                    MADDQ_S(S2, S3, Y23, Y45, W2);
                    MADDQ_S(S0, S1, Y34, Y56, W3);
                    MADDQ_S(S2, S3, Y45, Y67, W4);
                    MADDQ_S(S0, S1, Y56, Y78, W5);
                    MADDQ_S(S2, S3, Y67, Y89, W6);
                    ADD_SX2X2(S0, S1, S0, S1, S2, S3);

                    AE_SSX2X2_IP(S0, S1, pT, sizeof(xtfloatx4));

                    for (j = Q - (Q & 3); j < Q; ++j)
                    {
                        Y01 = Y12; Y12 = Y23; Y23 = Y34; Y34 = Y45; Y45 = Y56; Y56 = Y67; Y67 = Y78; Y78 = Y89; Y89 = Y910;
                    }
                }
                //last N-1 elements
                {
                    AE_LSX2X2_IP(S0, S1, pT_read, sizeof(xtfloatx4));

                    Y910 = CONJC_SX2(0);
                    Y89 = CONST_SX2(0);
                    Y78 = CONST_SX2(0);
                    Y67 = CONST_SX2(0);
                    Y56 = AE_SEL32_HL_SX2(Y56, CONST_SX2(0));

                    MULQ_S(S2, S3, Y01, Y23, W0);
                    MADDQ_S(S0, S1, Y12, Y34, W1);
                    MADDQ_S(S2, S3, Y23, Y45, W2);
                    MADDQ_S(S0, S1, Y34, Y56, W3);
                    MADDQ_S(S2, S3, Y45, Y67, W4);
                    MADDQ_S(S0, S1, Y56, Y78, W5);
                    MADDQ_S(S2, S3, Y67, Y89, W6);
                    ADD_SX2X2(S0, S1, S0, S1, S2, S3);

                    AE_SSX2X2_IP(S0, S1, pT, sizeof(xtfloatx4));

                    Y01 = Y45; Y12 = Y56;
                    Y23 = Y67; Y34 = Y78;
                    Y45 = Y89; Y56 = Y910;

                    AE_LSX2X2_IP(S0, S1, pT_read, sizeof(xtfloatx4));

                    MULQ_S(S2, S3, Y01, Y23, W0);
                    MADDQ_S(S0, S1, Y12, Y34, W1);
                    MADDQ_S(S2, S3, Y23, Y45, W2);
                    MADDQ_S(S0, S1, Y34, Y56, W3);
                    MADDQ_S(S2, S3, Y45, Y67, W4);
                    MADDQ_S(S0, S1, Y56, Y78, W5);
                    MADDQ_S(S2, S3, Y67, Y89, W6);
                    ADD_SX2X2(S0, S1, S0, S1, S2, S3);

                    AE_SSX2X2_IP(S0, S1, pT, sizeof(xtfloatx4));
                }

                pY = (xtfloatx4*)((xtfloat*)pY + Q);
                pY0 = pY;
            }

            //last iterations
            {
                pT_read = (xtfloatx4*)pScr + 4 * 11;

                Y01 = CONST_SX2(0);
                Y12 = CONST_SX2(0);
                Y23 = CONST_SX2(0);
                Y34 = CONST_SX2(0);
                Y45 = CONST_SX2(0);

                al0 = AE_LA128_PP(pY0);
                AE_LASX2X2_IP(Y67, Y89, al0, pY0);
                Y56 = AE_SEL32_LH_SX2(Y45, Y67);
                pY0 = (xtfloatx4*)((xtfloat*)(pY0)-4);

                pY1 = (xtfloatx4*)((xtfloat*)pY0 + 1);
                al0 = AE_LA128_PP(pY0);
                al1 = AE_LA128_PP(pY1);

                AE_LSX2X2_IP(W0, W1, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W2, W3, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W4, W5, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W6, dummy, pW, sizeof(xtfloatx4));

                for (j=0; j < (Q & ~3); j += 4)
                {
                    AE_LSX2X2_IP(S0, S1, pT_read, sizeof(xtfloatx4));

                    AE_LASX2X2_IP(Y67, Y89, al0, pY0);
                    AE_LASX2X2_IP(Y78, Y910, al1, pY1);

                    MULQ_S(S2, S3, Y01, Y23, W0);
                    MADDQ_S(S0, S1, Y12, Y34, W1);
                    MADDQ_S(S2, S3, Y23, Y45, W2);
                    MADDQ_S(S0, S1, Y34, Y56, W3);
                    MADDQ_S(S2, S3, Y45, Y67, W4);
                    MADDQ_S(S0, S1, Y56, Y78, W5);
                    MADDQ_S(S2, S3, Y67, Y89, W6);
                    ADD_SX2X2(S0, S1, S0, S1, S2, S3);

                    AE_SASX2X2_IP(S0, S1, alZ, pZ0);

                    Y01 = Y45; Y12 = Y56;
                    Y23 = Y67; Y34 = Y78;
                    Y45 = Y89; Y56 = Y910;
                }
                AE_SA128POS_FP(alZ, pZ0);
                // tail
                if (Q & 3)
                {
                    AE_LSX2X2_IP(S0, S1, pT_read, sizeof(xtfloatx4));

                    AE_LASX2X2_IP(Y67, Y89, al0, pY0);
                    AE_LASX2X2_IP(Y78, Y910, al1, pY1);

                    MULQ_S(S2, S3, Y01, Y23, W0);
                    MADDQ_S(S0, S1, Y12, Y34, W1);
                    MADDQ_S(S2, S3, Y23, Y45, W2);
                    MADDQ_S(S0, S1, Y34, Y56, W3);
                    MADDQ_S(S2, S3, Y45, Y67, W4);
                    MADDQ_S(S0, S1, Y56, Y78, W5);
                    MADDQ_S(S2, S3, Y67, Y89, W6);
                    ADD_SX2X2(S0, S1, S0, S1, S2, S3);

                    for (; j < Q; ++j)
                    {
                        xtfloatx2 tmp;

                        tmp = AE_SEL32_LH_SX2(S0, S1); S1 = AE_SEL32_LH_SX2(S1, S0); S0 = tmp;

                        Y01 = Y12; Y12 = Y23; Y23 = Y34; Y34 = Y45; Y45 = Y56; Y56 = Y67; Y67 = Y78; Y78 = Y89; Y89 = Y910;

                        AE_S32_L_IP(AE_MOVINT32X2_FROMXTFLOATX2(S1), castxcc(ae_int32, pZ0), sizeof(ae_int32));
                    }
                }
                //last N-1 elements
                {
                    AE_LSX2X2_IP(S0, S1, pT_read, sizeof(xtfloatx4));

                    xtfloatx2 tmp;

                    Y910 = CONJC_SX2(0);
                    Y89 = CONST_SX2(0);
                    Y78 = CONST_SX2(0);
                    Y67 = CONST_SX2(0);
                    Y56 = AE_SEL32_HL_SX2(Y56, CONST_SX2(0));

                    MULQ_S(S2, S3, Y01, Y23, W0);
                    MADDQ_S(S0, S1, Y12, Y34, W1);
                    MADDQ_S(S2, S3, Y23, Y45, W2);
                    MADDQ_S(S0, S1, Y34, Y56, W3);
                    MADDQ_S(S2, S3, Y45, Y67, W4);
                    MADDQ_S(S0, S1, Y56, Y78, W5);
                    MADDQ_S(S2, S3, Y67, Y89, W6);
                    ADD_SX2X2(S0, S1, S0, S1, S2, S3);

                    AE_SASX2X2_IP(S0, S1, alZ, pZ0);
                    AE_SA128POS_FP(alZ, pZ0);

                    Y01 = Y45; Y12 = Y56;
                    Y23 = Y67; Y34 = Y78;
                    Y45 = Y89; Y56 = Y910;

                    AE_LSX2X2_IP(S0, S1, pT_read, sizeof(xtfloatx4));

                    MULQ_S(S2, S3, Y01, Y23, W0);
                    MADDQ_S(S0, S1, Y12, Y34, W1);
                    MADDQ_S(S2, S3, Y23, Y45, W2);
                    MADDQ_S(S0, S1, Y34, Y56, W3);
                    MADDQ_S(S2, S3, Y45, Y67, W4);
                    MADDQ_S(S0, S1, Y56, Y78, W5);
                    MADDQ_S(S2, S3, Y67, Y89, W6);
                    ADD_SX2X2(S0, S1, S0, S1, S2, S3);

                    tmp = AE_SEL32_LH_SX2(S0, S1); S1 = AE_SEL32_LH_SX2(S1, S0); S0 = tmp;
                    AE_S32_L_IP(AE_MOVINT32X2_FROMXTFLOATX2(S1), castxcc(ae_int32, pZ0), sizeof(ae_int32));

                    tmp = AE_SEL32_LH_SX2(S0, S1); S1 = AE_SEL32_LH_SX2(S1, S0); S0 = tmp;
                    AE_S32_L_IP(AE_MOVINT32X2_FROMXTFLOATX2(S1), castxcc(ae_int32, pZ0), sizeof(ae_int32));
                }
            }

        }
    }
#undef M
#undef N
}

size_t conv2d_gen_11x7f_getScratchSize(int P, int Q)
{
    const int M = 11;
    (void)P;
    return 4 * M * sizeof(xtfloatx4) + (Q + 16) * sizeof(xtfloat);
} // MxN=11x7


#else 

void conv2d_gen_11x7f(void* pScr, float32_t* z, const float32_t* x, const float32_t* y, int P, int Q)
{
#define M 11
#define N 7
    const xtfloat* restrict pY;
    const xtfloat* restrict pY0;
    const xtfloat* restrict pT_read;
    xtfloat* restrict pT;
    xtfloat* restrict pZ0;
    xtfloat* restrict pW;

    xtfloat y0_0, y0_1, y0_2, y0_3, y0_4, y0_5, y0_6;
    xtfloat s0;
    xtfloat w0, w1, w2, w3, w4, w5, w6;
    int i, j, m, mmin, mmax, xstart, ystart;

    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0);
    NASSERT(Q >= 0);
    if (P <= 0 || Q <= 0) return;

    pW = (xtfloat*)pScr + (M - 1) * N;
    // preload X and save in scratch
    {
        for (i = 0; i < M; ++i)
        {
            XT_LSIP(w6, x, sizeof(xtfloat));
            XT_LSIP(w5, x, sizeof(xtfloat));
            XT_LSIP(w4, x, sizeof(xtfloat));
            XT_LSIP(w3, x, sizeof(xtfloat));
            XT_LSIP(w2, x, sizeof(xtfloat));
            XT_LSIP(w1, x, sizeof(xtfloat));
            XT_LSIP(w0, x, sizeof(xtfloat));

            XT_SSIP(w0, pW, sizeof(xtfloat));
            XT_SSIP(w1, pW, sizeof(xtfloat));
            XT_SSIP(w2, pW, sizeof(xtfloat));
            XT_SSIP(w3, pW, sizeof(xtfloat));
            XT_SSIP(w4, pW, sizeof(xtfloat));
            XT_SSIP(w5, pW, sizeof(xtfloat));
            XT_SSXP(w6, pW, -13 * (int)sizeof(xtfloat));
        }
    }

    for (i = 0; i < P + M - 1; i += 1)
    {
        mmin = i < M - 1 ? M - 1 - i : 0;
        mmax = i >= P ? M - 1 + P - i : M;
        xstart = i < M - 1 ? (mmin) : 0;
        ystart = i < M - 1 ? 0 : i - (M - 1);

        pW = (xtfloat*)pScr + N * xstart;
        pY0 = (xtfloat*)(y + ystart * Q);

        pZ0 = (xtfloat*)(z + i * (Q + N - 1));

        if (mmax - mmin == 1)
        {

            y0_0 = y0_1 = y0_2 = y0_3 = y0_4 = y0_5 = XT_CONST_S(0);
            XT_LSIP(w0, pW, sizeof(xtfloat));
            XT_LSIP(w1, pW, sizeof(xtfloat));
            XT_LSIP(w2, pW, sizeof(xtfloat));
            XT_LSIP(w3, pW, sizeof(xtfloat));
            XT_LSIP(w4, pW, sizeof(xtfloat));
            XT_LSIP(w5, pW, sizeof(xtfloat));
            XT_LSIP(w6, pW, sizeof(xtfloat));

            for (j = 0; j < Q; ++j)
            {
                XT_LSIP(y0_6, pY0, sizeof(xtfloat));

                s0 = XT_MUL_S(w0, y0_0);
                XT_MADD_S(s0, w1, y0_1);
                XT_MADD_S(s0, w2, y0_2);
                XT_MADD_S(s0, w3, y0_3);
                XT_MADD_S(s0, w4, y0_4);
                XT_MADD_S(s0, w5, y0_5);
                XT_MADD_S(s0, w6, y0_6);
                XT_SSIP(s0, pZ0, sizeof(xtfloat));

                y0_0 = y0_1; y0_1 = y0_2; y0_2 = y0_3; 
                y0_3 = y0_4; y0_4 = y0_5; y0_5 = y0_6;
            }
            for (; j < Q + N - 1; ++j)
            {
                y0_6 = XT_CONST_S(0);
                s0 = XT_MUL_S(w0, y0_0);
                XT_MADD_S(s0, w1, y0_1);
                XT_MADD_S(s0, w2, y0_2);
                XT_MADD_S(s0, w3, y0_3);
                XT_MADD_S(s0, w4, y0_4);
                XT_MADD_S(s0, w5, y0_5);
                XT_MADD_S(s0, w6, y0_6);
                XT_SSIP(s0, pZ0, sizeof(xtfloat));
                y0_0 = y0_1; y0_1 = y0_2; y0_2 = y0_3;
                y0_3 = y0_4; y0_4 = y0_5; y0_5 = y0_6;
            }
        }
        else
        {
            pY = pY0 + Q;
            {
                pT = (xtfloat*)pScr + M * N;

                y0_0 = y0_1 = y0_2 = y0_3 = y0_4 = y0_5 = XT_CONST_S(0);
                XT_LSIP(w0, pW, sizeof(xtfloat));
                XT_LSIP(w1, pW, sizeof(xtfloat));
                XT_LSIP(w2, pW, sizeof(xtfloat));
                XT_LSIP(w3, pW, sizeof(xtfloat));
                XT_LSIP(w4, pW, sizeof(xtfloat));
                XT_LSIP(w5, pW, sizeof(xtfloat));
                XT_LSIP(w6, pW, sizeof(xtfloat));

                for (j = 0; j < Q; ++j)
                {
                    XT_LSIP(y0_6, pY0, sizeof(xtfloat));

                    s0 = XT_MUL_S(w0, y0_0);
                    XT_MADD_S(s0, w1, y0_1);
                    XT_MADD_S(s0, w2, y0_2);
                    XT_MADD_S(s0, w3, y0_3);
                    XT_MADD_S(s0, w4, y0_4);
                    XT_MADD_S(s0, w5, y0_5);
                    XT_MADD_S(s0, w6, y0_6);
                    XT_SSIP(s0, pT, sizeof(xtfloat));

                    y0_0 = y0_1; y0_1 = y0_2; y0_2 = y0_3;
                    y0_3 = y0_4; y0_4 = y0_5; y0_5 = y0_6;
                }
                for (; j < Q + N - 1; ++j)
                {
                    y0_6 = XT_CONST_S(0);
                    s0 = XT_MUL_S(w0, y0_0);
                    XT_MADD_S(s0, w1, y0_1);
                    XT_MADD_S(s0, w2, y0_2);
                    XT_MADD_S(s0, w3, y0_3);
                    XT_MADD_S(s0, w4, y0_4);
                    XT_MADD_S(s0, w5, y0_5);
                    XT_MADD_S(s0, w6, y0_6);
                    XT_SSIP(s0, pT, sizeof(xtfloat));
                    
                    y0_0 = y0_1; y0_1 = y0_2; y0_2 = y0_3;
                    y0_3 = y0_4; y0_4 = y0_5; y0_5 = y0_6;
                }
            }

            for (m = mmin + 1; m < mmax - 1; ++m)
            {
                pY0 = pY;
                pY += Q;
                pT_read = pT = (xtfloat*)pScr + M * N;

                y0_0 = y0_1 = y0_2 = y0_3 = y0_4 = y0_5 = XT_CONST_S(0);
                XT_LSIP(w0, pW, sizeof(xtfloat));
                XT_LSIP(w1, pW, sizeof(xtfloat));
                XT_LSIP(w2, pW, sizeof(xtfloat));
                XT_LSIP(w3, pW, sizeof(xtfloat));
                XT_LSIP(w4, pW, sizeof(xtfloat));
                XT_LSIP(w5, pW, sizeof(xtfloat));
                XT_LSIP(w6, pW, sizeof(xtfloat));

                for (j = 0; j < Q; ++j)
                {
                    XT_LSIP(y0_6, pY0, sizeof(xtfloat));
                    XT_LSIP(s0, pT_read, sizeof(xtfloat));

                    XT_MADD_S(s0, w0, y0_0);
                    XT_MADD_S(s0, w1, y0_1);
                    XT_MADD_S(s0, w2, y0_2);
                    XT_MADD_S(s0, w3, y0_3);
                    XT_MADD_S(s0, w4, y0_4);
                    XT_MADD_S(s0, w5, y0_5);
                    XT_MADD_S(s0, w6, y0_6);
                    XT_SSIP(s0, pT, sizeof(xtfloat));

                    y0_0 = y0_1; y0_1 = y0_2; y0_2 = y0_3;
                    y0_3 = y0_4; y0_4 = y0_5; y0_5 = y0_6;
                }
                for (; j < Q + N - 1; ++j)
                {
                    XT_LSIP(s0, pT_read, sizeof(xtfloat));
                    y0_6 = XT_CONST_S(0);
                    XT_MADD_S(s0, w0, y0_0);
                    XT_MADD_S(s0, w1, y0_1);
                    XT_MADD_S(s0, w2, y0_2);
                    XT_MADD_S(s0, w3, y0_3);
                    XT_MADD_S(s0, w4, y0_4);
                    XT_MADD_S(s0, w5, y0_5);
                    XT_MADD_S(s0, w6, y0_6);
                    XT_SSIP(s0, pT, sizeof(xtfloat));

                    y0_0 = y0_1; y0_1 = y0_2; y0_2 = y0_3;
                    y0_3 = y0_4; y0_4 = y0_5; y0_5 = y0_6;
                }
            }

            {
                pY0 = pY;

                pT_read = (xtfloat*)pScr + M * N;

                y0_0 = y0_1 = y0_2 = y0_3 = y0_4 = y0_5 = XT_CONST_S(0);
                XT_LSIP(w0, pW, sizeof(xtfloat));
                XT_LSIP(w1, pW, sizeof(xtfloat));
                XT_LSIP(w2, pW, sizeof(xtfloat));
                XT_LSIP(w3, pW, sizeof(xtfloat));
                XT_LSIP(w4, pW, sizeof(xtfloat));
                XT_LSIP(w5, pW, sizeof(xtfloat));
                XT_LSIP(w6, pW, sizeof(xtfloat));

                for (j = 0; j < Q; ++j)
                {
                    XT_LSIP(y0_6, pY0, sizeof(xtfloat));
                    XT_LSIP(s0, pT_read, sizeof(xtfloat));

                    XT_MADD_S(s0, w0, y0_0);
                    XT_MADD_S(s0, w1, y0_1);
                    XT_MADD_S(s0, w2, y0_2);
                    XT_MADD_S(s0, w3, y0_3);
                    XT_MADD_S(s0, w4, y0_4);
                    XT_MADD_S(s0, w5, y0_5);
                    XT_MADD_S(s0, w6, y0_6);
                    XT_SSIP(s0, pZ0, sizeof(xtfloat));

                    y0_0 = y0_1; y0_1 = y0_2; y0_2 = y0_3;
                    y0_3 = y0_4; y0_4 = y0_5; y0_5 = y0_6;
                }
                for (; j < Q + N - 1; ++j)
                {
                    XT_LSIP(s0, pT_read, sizeof(xtfloat));
                    y0_6 = XT_CONST_S(0);
                    XT_MADD_S(s0, w0, y0_0);
                    XT_MADD_S(s0, w1, y0_1);
                    XT_MADD_S(s0, w2, y0_2);
                    XT_MADD_S(s0, w3, y0_3);
                    XT_MADD_S(s0, w4, y0_4);
                    XT_MADD_S(s0, w5, y0_5);
                    XT_MADD_S(s0, w6, y0_6);
                    XT_SSIP(s0, pZ0, sizeof(xtfloat));

                    y0_0 = y0_1; y0_1 = y0_2; y0_2 = y0_3;
                    y0_3 = y0_4; y0_4 = y0_5; y0_5 = y0_6;
                }
            }

        }

    }
#undef M
#undef N
}

// scratch allocatation functions. return required scratch size in bytes
size_t conv2d_gen_11x7f_getScratchSize(int P, int Q)
{
    (void)P;
    const int M = 11, N = 7;
    return M * N * sizeof(xtfloat) + (Q + N - 1) * sizeof(xtfloat);
}
#endif
