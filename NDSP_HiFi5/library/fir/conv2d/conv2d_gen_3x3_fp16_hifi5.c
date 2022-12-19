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
  IntegrIT, 2006-2018
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
#include "common.h"
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

#if HAVE_HPFPU

void conv2d_gen_3x3_fp16(void* pScr, float16_t* z, const float16_t* x, const float16_t* y, int P, int Q)
{
#define M 3
#define N 3
    const xthalfx8 * restrict pY0;
    const xthalfx8 * restrict pY1;
    const xthalfx8 * restrict pY2;
    xthalfx8* restrict pW;
    xthalfx8* restrict pZ0;

    const xthalf* pX = (const xthalf*)x;

    ae_valignx2 al0, al1, al2, alZ;

    int i, j, mmin, mmax, xstart, ystart;
    xthalfx4 Y00, Y01, Y02, Y03, Y04;
    xthalfx4 Y10, Y11, Y12, Y13, Y14;
    xthalfx4 Y20, Y21, Y22, Y23, Y24;
    xthalfx4 W0, W1, W2;
    xthalfx4 dummy;
    xthalfx4 S0, S1, S2, S3, S4, S5;
    xthalfx4 S6, S7, S8, S9, S10, S11;
    ae_int16x4 sel, t0, t1, t2, t3;
    xthalfx4 SEL0, SEL1, SEL2, SEL3;
    xthalfx4 SEL0_, SEL1_, SEL2_, SEL3_;

    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0);
    NASSERT(Q >= 0);
    if (P <= 0 || Q <= 0) return;


    alZ = AE_ZALIGN128();

    pW = (xthalfx8*)pScr;

    // preload X and save in scratch
    {
        xthalfx4 w0, w1, w2, w3, w4, w5, w6, w7, w8;
        ae_int16x4 v_i;
        
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
        w0 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
        w1 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
        w2 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
        w3 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
        w4 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
        w5 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
        w6 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
        w7 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
        w8 = AE_MOVXTHALFX4_FROMINT16X4(v_i);

        AE_SHX4X2_IP(w8, w7, pW, sizeof(xthalfx8));
        AE_SHX4X2_IP(w6, ZERO_HX4(), pW, sizeof(xthalfx8));
        AE_SHX4X2_IP(w5, w4, pW, sizeof(xthalfx8));
        AE_SHX4X2_IP(w3, ZERO_HX4(), pW, sizeof(xthalfx8));
        AE_SHX4X2_IP(w2, w1, pW, sizeof(xthalfx8));
        AE_SHX4X2_IP(w0, ZERO_HX4(), pW, sizeof(xthalfx8));
    }

    sel = AE_MOVINT16X4_FROMINT64(0x0504040303020201); // sel 5432 4321

    __Pragma("loop_count min=1");
    for (i = 0; i < P + M - 1; i += 1)
    {
        mmin = i < M - 1 ? M - 1 - i : 0;
        mmax = i >= P ? M - 1 + P - i : M;
        xstart = i < M - 1 ? (mmin) : 0;
        ystart = i < M - 1 ? 0 : i - (M - 1);

        pW  = (xthalfx8*)pScr + 2*xstart;
        pY0 = (xthalfx8*)(y + ystart * Q);
        pY1 = (xthalfx8*)(y + (ystart + 1) * Q);
        pY2 = (xthalfx8*)(y + (ystart + 2) * Q);

        pZ0 = (xthalfx8*)(z + i * (Q + N - 1));

        // all 3 row needed
        if (mmax - mmin == 3)
        {
            Y00 = ZERO_HX4();;
            Y10 = ZERO_HX4();;
            Y20 = ZERO_HX4();;

            al0 = AE_LA128_PP(pY0);
            al1 = AE_LA128_PP(pY1);
            al2 = AE_LA128_PP(pY2);


            for (j = 0; j < (Q & ~15); j += 16)
            {
                AE_LAHX4X2_IP(Y01, Y02, al0, pY0); AE_LAHX4X2_IP(Y03, Y04, al0, pY0);
                AE_LAHX4X2_IP(Y11, Y12, al1, pY1); AE_LAHX4X2_IP(Y13, Y14, al1, pY1);
                AE_LAHX4X2_IP(Y21, Y22, al2, pY2); AE_LAHX4X2_IP(Y23, Y24, al2, pY2);

                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_IP(W2, dummy, pW, sizeof(xthalfx8));

                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y02), AE_MOVF16X4_FROMHALFX4(Y03), sel); SEL0_ = AE_MOVHALFX4_FROMF16X4(t2); SEL1_ = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y03), AE_MOVF16X4_FROMHALFX4(Y04), sel); SEL2_ = AE_MOVHALFX4_FROMF16X4(t2); SEL3_ = AE_MOVHALFX4_FROMF16X4(t3);

                MULQ_H(S0, S1, SEL0, SEL2, W0); MULQ_H(S2, S3, SEL1, SEL3, W1); MULQ_H(S4, S5, Y01, Y02, W2);
                MULQ_H(S6, S7, SEL0_, SEL2_, W0); MULQ_H(S8, S9, SEL1_, SEL3_, W1); MULQ_H(S10, S11, Y03, Y04, W2);

                Y00 = Y04;

                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_IP(W2, dummy, pW, sizeof(xthalfx8));


                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y10), AE_MOVF16X4_FROMHALFX4(Y11), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y11), AE_MOVF16X4_FROMHALFX4(Y12), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y12), AE_MOVF16X4_FROMHALFX4(Y13), sel); SEL0_ = AE_MOVHALFX4_FROMF16X4(t2); SEL1_ = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y13), AE_MOVF16X4_FROMHALFX4(Y14), sel); SEL2_ = AE_MOVHALFX4_FROMF16X4(t2); SEL3_ = AE_MOVHALFX4_FROMF16X4(t3);

                MADDQ_H(S0, S1, SEL0, SEL2, W0);   MADDQ_H(S2, S3, SEL1, SEL3, W1);   MADDQ_H(S4, S5, Y11, Y12, W2);
                MADDQ_H(S6, S7, SEL0_, SEL2_, W0); MADDQ_H(S8, S9, SEL1_, SEL3_, W1); MADDQ_H(S10, S11, Y13, Y14, W2);

                Y10 = Y14;
        
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_XP(W2, dummy, pW, -5 * (int)sizeof(xthalfx8));


                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y20), AE_MOVF16X4_FROMHALFX4(Y21), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y21), AE_MOVF16X4_FROMHALFX4(Y22), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y22), AE_MOVF16X4_FROMHALFX4(Y23), sel); SEL0_ = AE_MOVHALFX4_FROMF16X4(t2); SEL1_ = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y23), AE_MOVF16X4_FROMHALFX4(Y24), sel); SEL2_ = AE_MOVHALFX4_FROMF16X4(t2); SEL3_ = AE_MOVHALFX4_FROMF16X4(t3);

                MADDQ_H(S0, S1, SEL0, SEL2, W0);   MADDQ_H(S2, S3, SEL1, SEL3, W1);   MADDQ_H(S4, S5, Y21, Y22, W2);
                MADDQ_H(S6, S7, SEL0_, SEL2_, W0); MADDQ_H(S8, S9, SEL1_, SEL3_, W1); MADDQ_H(S10, S11, Y23, Y24, W2);

                Y20 = Y24;

                ADD_HX4X2(S0, S1, S0, S1, S2, S3); ADD_HX4X2(S0, S1, S0, S1, S4, S5);
                ADD_HX4X2(S6, S7, S6, S7, S8, S9); ADD_HX4X2(S6, S7, S6, S7, S10, S11);

                AE_SAHX4X2_IP(S0, S1, alZ, pZ0);  AE_SAHX4X2_IP(S6, S7, alZ, pZ0);
            }
            if ((Q & ~7)>(Q & ~15))
            {
                AE_LAHX4X2_IP(Y01, Y02, al0, pY0);
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_IP(W2, dummy, pW, sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MULQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y01, Y02, W2);
                Y00 = Y02;

                AE_LAHX4X2_IP(Y11, Y12, al1, pY1);
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_IP(W2, dummy, pW, sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y10), AE_MOVF16X4_FROMHALFX4(Y11), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y11), AE_MOVF16X4_FROMHALFX4(Y12), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MULQ_H(S2, S3, SEL0, SEL2, W0);
                MADDQ_H(S2, S3, SEL1, SEL3, W1);
                MADDQ_H(S2, S3, Y11, Y12, W2);
                Y10 = Y12;

                AE_LAHX4X2_IP(Y21, Y22, al2, pY2);
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_XP(W2, dummy, pW, -5 * (int)sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y20), AE_MOVF16X4_FROMHALFX4(Y21), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y21), AE_MOVF16X4_FROMHALFX4(Y22), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MULQ_H(S4, S5, SEL0, SEL2, W0);
                MADDQ_H(S4, S5, SEL1, SEL3, W1);
                MADDQ_H(S4, S5, Y21, Y22, W2);
                Y20 = Y22;

                ADD_HX4X2(S0, S1, S0, S1, S2, S3);
                ADD_HX4X2(S0, S1, S0, S1, S4, S5);

                AE_SAHX4X2_IP(S0, S1, alZ, pZ0);
            }
       
            AE_SA128POS_FP(alZ, pZ0);
            // tail
            if (Q&7)
            {

                AE_LAHX4X2_IP(Y01, Y02, al0, pY0);
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_IP(W2, dummy, pW, sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MULQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y01, Y02, W2);

                AE_LAHX4X2_IP(Y11, Y12, al1, pY1);
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_IP(W2, dummy, pW, sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y10), AE_MOVF16X4_FROMHALFX4(Y11), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y11), AE_MOVF16X4_FROMHALFX4(Y12), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MADDQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y11, Y12, W2);

                AE_LAHX4X2_IP(Y21, Y22, al2, pY2);
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_XP(W2, dummy, pW, -5 * (int)sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y20), AE_MOVF16X4_FROMHALFX4(Y21), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y21), AE_MOVF16X4_FROMHALFX4(Y22), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MADDQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y21, Y22, W2);

                for (j = Q - (Q&7); j < Q; ++j)
                {
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(S0), AE_MOVF16X4_FROMHALFX4(S1), AE_MOVINT16X4_FROMINT64(0x0602050104000307));
                    S0 = AE_MOVHALFX4_FROMF16X4(t0); S1 = AE_MOVHALFX4_FROMF16X4(t1);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    Y00 = AE_MOVHALFX4_FROMF16X4(t0); Y01 = AE_MOVHALFX4_FROMF16X4(t1);
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    Y01 = AE_MOVHALFX4_FROMF16X4(t0); Y02 = AE_MOVHALFX4_FROMF16X4(t1);
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y10), AE_MOVF16X4_FROMHALFX4(Y11), AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    Y10 = AE_MOVHALFX4_FROMF16X4(t0); Y11 = AE_MOVHALFX4_FROMF16X4(t1);
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y11), AE_MOVF16X4_FROMHALFX4(Y12), AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    Y11 = AE_MOVHALFX4_FROMF16X4(t0); Y12 = AE_MOVHALFX4_FROMF16X4(t1);
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y20), AE_MOVF16X4_FROMHALFX4(Y21), AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    Y20 = AE_MOVHALFX4_FROMF16X4(t0); Y21 = AE_MOVHALFX4_FROMF16X4(t1);
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y21), AE_MOVF16X4_FROMHALFX4(Y22), AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    Y21 = AE_MOVHALFX4_FROMF16X4(t0); Y22 = AE_MOVHALFX4_FROMF16X4(t1);

                    AE_S16_0_IP(AE_MOVINT16X4_FROMXTHALFX4(S1), castxcc(ae_int16, pZ0), sizeof(ae_int16));
                }
            }
            //last N-1 elements
            {
                Y01 = ZERO_HX4(); Y02 = ZERO_HX4();
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_IP(W2, dummy, pW, sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MULQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y01, Y02, W2);

                Y11 = ZERO_HX4(); Y12 = ZERO_HX4();
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_IP(W2, dummy, pW, sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y10), AE_MOVF16X4_FROMHALFX4(Y11), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y11), AE_MOVF16X4_FROMHALFX4(Y12), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MADDQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y11, Y12, W2);

                Y21 = ZERO_HX4(); Y22 = ZERO_HX4();
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_XP(W2, dummy, pW, -5 * (int)sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y20), AE_MOVF16X4_FROMHALFX4(Y21), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y21), AE_MOVF16X4_FROMHALFX4(Y22), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MADDQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y21, Y22, W2);

                for (j=Q; j < Q + N - 1; ++j)
                {
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(S0), AE_MOVF16X4_FROMHALFX4(S1), AE_MOVINT16X4_FROMINT64(0x0602050104000307));
                    S0 = AE_MOVHALFX4_FROMF16X4(t0); S1 = AE_MOVHALFX4_FROMF16X4(t1);

                    AE_S16_0_IP(AE_MOVINT16X4_FROMXTHALFX4(S1), castxcc(ae_int16, pZ0), sizeof(ae_int16));
                }
            }
        }
        // if only one Y row need to 
        else if (mmax - mmin == 1)
        {
            Y00 = ZERO_HX4();

            al0 = AE_LA128_PP(pY0);

            for (j = 0; j < (Q & ~7); j += 8)
            {

                AE_LAHX4X2_IP(Y01, Y02, al0, pY0);
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_XP(W2, dummy, pW, -(int)sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MULQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y01, Y02, W2);
                Y00 = Y02;

                AE_SAHX4X2_IP(S0, S1, alZ, pZ0);

            }
            AE_SA128POS_FP(alZ, pZ0);
            // tail
            if (j < Q)
            {

                AE_LAHX4X2_IP(Y01, Y02, al0, pY0);
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_XP(W2, dummy, pW, -(int)sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MULQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y01, Y02, W2);


                for (; j < Q; ++j)
                {
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(S0), AE_MOVF16X4_FROMHALFX4(S1), AE_MOVINT16X4_FROMINT64(0x0602050104000307));
                    S0 = AE_MOVHALFX4_FROMF16X4(t0); S1 = AE_MOVHALFX4_FROMF16X4(t1);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    Y00 = AE_MOVHALFX4_FROMF16X4(t0); Y01 = AE_MOVHALFX4_FROMF16X4(t1);
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    Y01 = AE_MOVHALFX4_FROMF16X4(t0); Y02 = AE_MOVHALFX4_FROMF16X4(t1);

                    AE_S16_0_IP(AE_MOVINT16X4_FROMXTHALFX4(S1), castxcc(ae_int16, pZ0), sizeof(ae_int16));
                }
            }
            //last N-1 elements
            {
                Y01 = ZERO_HX4(); Y02 = ZERO_HX4();
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_XP(W2, dummy, pW, -(int)sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MULQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y01, Y02, W2);

                for (; j < Q + N - 1; ++j)
                {
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(S0), AE_MOVF16X4_FROMHALFX4(S1), AE_MOVINT16X4_FROMINT64(0x0602050104000307));
                    S0 = AE_MOVHALFX4_FROMF16X4(t0); S1 = AE_MOVHALFX4_FROMF16X4(t1);

                    AE_S16_0_IP(AE_MOVINT16X4_FROMXTHALFX4(S1), castxcc(ae_int16, pZ0), sizeof(ae_int16));
                }
            }
        }
        else if (mmax - mmin == 2)
        {
            Y00 = ZERO_HX4();;
            Y10 = ZERO_HX4();;

            al0 = AE_LA128_PP(pY0);
            al1 = AE_LA128_PP(pY1);


            for (j = 0; j < (Q & ~7); j += 8)
            {

                AE_LAHX4X2_IP(Y01, Y02, al0, pY0);
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_IP(W2, dummy, pW, sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MULQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y01, Y02, W2);
                Y00 = Y02;

                AE_LAHX4X2_IP(Y11, Y12, al1, pY1);
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_XP(W2, dummy, pW, -3*(int)sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y10), AE_MOVF16X4_FROMHALFX4(Y11), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y11), AE_MOVF16X4_FROMHALFX4(Y12), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MADDQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y11, Y12, W2);
                Y10 = Y12;

                AE_SAHX4X2_IP(S0, S1, alZ, pZ0);

            }
            AE_SA128POS_FP(alZ, pZ0);
            // tail
            if (j < Q)
            {

                AE_LAHX4X2_IP(Y01, Y02, al0, pY0);
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_IP(W2, dummy, pW, sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MULQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y01, Y02, W2);

                AE_LAHX4X2_IP(Y11, Y12, al1, pY1);
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_XP(W2, dummy, pW, -3 * (int)sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y10), AE_MOVF16X4_FROMHALFX4(Y11), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y11), AE_MOVF16X4_FROMHALFX4(Y12), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MADDQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y11, Y12, W2);

                for (; j < Q; ++j)
                {
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(S0), AE_MOVF16X4_FROMHALFX4(S1), AE_MOVINT16X4_FROMINT64(0x0602050104000307));
                    S0 = AE_MOVHALFX4_FROMF16X4(t0); S1 = AE_MOVHALFX4_FROMF16X4(t1);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    Y00 = AE_MOVHALFX4_FROMF16X4(t0); Y01 = AE_MOVHALFX4_FROMF16X4(t1);
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    Y01 = AE_MOVHALFX4_FROMF16X4(t0); Y02 = AE_MOVHALFX4_FROMF16X4(t1);
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y10), AE_MOVF16X4_FROMHALFX4(Y11), AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    Y10 = AE_MOVHALFX4_FROMF16X4(t0); Y11 = AE_MOVHALFX4_FROMF16X4(t1);
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y11), AE_MOVF16X4_FROMHALFX4(Y12), AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    Y11 = AE_MOVHALFX4_FROMF16X4(t0); Y12 = AE_MOVHALFX4_FROMF16X4(t1);

                    AE_S16_0_IP(AE_MOVINT16X4_FROMXTHALFX4(S1), castxcc(ae_int16, pZ0), sizeof(ae_int16));
                }
            }
            //last N-1 elements
            {
                Y01 = ZERO_HX4(); Y02 = ZERO_HX4();
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_IP(W2, dummy, pW, sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MULQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y01, Y02, W2);

                Y11 = ZERO_HX4(); Y12 = ZERO_HX4();
                AE_LHX4X2_IP(W0, W1, pW, sizeof(xthalfx8));
                AE_LHX4X2_XP(W2, dummy, pW, -3 * (int)sizeof(xthalfx8));
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y10), AE_MOVF16X4_FROMHALFX4(Y11), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t2); SEL1 = AE_MOVHALFX4_FROMF16X4(t3);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y11), AE_MOVF16X4_FROMHALFX4(Y12), sel); SEL2 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
                MADDQ_H(S0, S1, SEL0, SEL2, W0);
                MADDQ_H(S0, S1, SEL1, SEL3, W1);
                MADDQ_H(S0, S1, Y11, Y12, W2);

                for (; j < Q + N - 1; ++j)
                {
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(S0), AE_MOVF16X4_FROMHALFX4(S1), AE_MOVINT16X4_FROMINT64(0x0602050104000307));
                    S0 = AE_MOVHALFX4_FROMF16X4(t0); S1 = AE_MOVHALFX4_FROMF16X4(t1);

                    AE_S16_0_IP(AE_MOVINT16X4_FROMXTHALFX4(S1), castxcc(ae_int16, pZ0), sizeof(ae_int16));
                }
            }
        }
    }
#undef M
#undef N
}
size_t conv2d_gen_3x3_fp16_getScratchSize(int P, int Q)
{
    (void)P, (void)Q;
    return 3*4*sizeof(xthalfx4);
}
#else
DISCARD_FUN(void,conv2d_gen_3x3_fp16,(void* pScr, float16_t *z, const float16_t * x, const float16_t * y, int P, int Q))
size_t conv2d_gen_3x3_fp16_getScratchSize (int P, int Q) 
{ 
    (void)P,(void)Q;
    return 0;
}
#endif

