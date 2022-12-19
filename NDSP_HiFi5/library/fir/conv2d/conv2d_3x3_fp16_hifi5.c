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
void conv2d_3x3_fp16(void* pScr, float16_t *z, const float16_t * x, const float16_t * y, int P, int Q)
{
#   define M 3
#   define N 3
    int i, j;
    const xthalfx8* restrict pY0;
    const xthalfx8* restrict pY1;
    const xthalfx8* restrict pY2;
    const xthalfx8* restrict pY3;
    xthalfx8* restrict pW = (xthalfx8*) pScr;
    xthalfx8* restrict pZ;
    xthalfx8* restrict pZ1;

    ae_valignx2 alZ;
    ae_valignx2 alZ1;

    xthalfx4 S00, S10, S20, S30;
    xthalfx4 S01, S21, S11, S31;
    xthalfx4 Y00, Y01, Y02;
    xthalfx4 Y10, Y11, Y12;
    xthalfx4 Y20, Y21, Y22;
    xthalfx4 Y30, Y31, Y32;
    xthalfx4 SEL0, SEL1, SEL2, SEL3;
    xthalfx4 SEL0_, SEL1_, SEL2_, SEL3_;
    ae_int16x4 t0, t1, t2, t3;

    xthalfx4 W00, W01, W02;
    xthalfx4 W10, W11, W12;
    xthalfx4 W20, W21, W22;

    ae_int32x2 tmp;
    ae_int16x4 sel;
    // 15x ae64 regs

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    NASSERT(pScr);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT_ALIGN(z, 16);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0 && P % 8 == 0);
    NASSERT(Q >= 0 && Q % 8 == 0);

    if (P <= 0 || Q <= 0) return;

    //preload X
    {
        ae_int16x4 v_i;		
        const xtfloatx4* pX;
        pX = (xtfloatx4*)(x);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
	    W00 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
	    W01 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
	    W02 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
	    W10 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
	    W11 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
	    W12 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
	    W20 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
	    W21 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(xthalf));
	    W22 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
    }
    AE_SHX4X2_I(W22, W21, pW, 0*sizeof(xthalfx8));
    AE_SHX4X2_I(W12, W11, pW, 1*sizeof(xthalfx8));
    AE_SHX4X2_I(W02, W01, pW, 2*sizeof(xthalfx8));
    AE_SHX4X2_I(W20, W10, pW, 3*sizeof(xthalfx8));
    AE_SHX4X(W00, (xthalfx4*)pW, 8*sizeof(xthalfx4));

    alZ  = AE_ZALIGN128();
    alZ1 = AE_ZALIGN128();
    sel = AE_MOVINT16X4_FROMINT64(0x0504040303020201); // sel 5432 4321
    
                                                       
                                                       
    // Processing of convolution
    //first 2 rows 
    {
        pY2 = (xthalfx8*)(y);
        pY3 = (xthalfx8*)((float16_t*)pY2 + Q);

        Y20 = CONST_HX4(0);
        Y30 = CONST_HX4(0);

        pZ = (xthalfx8*)(z);
        pZ1 = (xthalfx8*)(z +  (N + Q - 1));

        __Pragma("loop_count min=1");
        for (j = 0; j < Q; j += 8)
        {
            AE_LHX4X2_IP(Y21, Y22, pY2, sizeof(xthalfx8));
            AE_LHX4X2_IP(Y31, Y32, pY3, sizeof(xthalfx8));

            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y20), AE_MOVF16X4_FROMHALFX4(Y21), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t0); SEL2 = AE_MOVHALFX4_FROMF16X4(t1);
            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y21), AE_MOVF16X4_FROMHALFX4(Y22), sel); SEL1 = AE_MOVHALFX4_FROMF16X4(t0); SEL3 = AE_MOVHALFX4_FROMF16X4(t1);

            MUL_HX4X2(S20, S30, W12, W12, SEL0, SEL1);
            MUL_HX4X2(S21, S31, W11, W11, SEL2, SEL3);
            MUL_HX4X2(S00, S10, W02, W02, SEL0, SEL1);
            MUL_HX4X2(S01, S11, W01, W01, SEL2, SEL3);

            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y30), AE_MOVF16X4_FROMHALFX4(Y31), sel); SEL0_ = AE_MOVHALFX4_FROMF16X4(t0); SEL2_ = AE_MOVHALFX4_FROMF16X4(t1);
            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y31), AE_MOVF16X4_FROMHALFX4(Y32), sel); SEL1_ = AE_MOVHALFX4_FROMF16X4(t0); SEL3_ = AE_MOVHALFX4_FROMF16X4(t1);

            MADD_HX4X2(S20, S30, W02, W02, SEL0_, SEL1_);
            MADD_HX4X2(S21, S31, W01, W01, SEL2_, SEL3_);
            MADD_HX4X2(S00, S10, W00, W00, Y21, Y22);
            ADD_HX4X2(S00, S10, S00, S10, S01, S11);

            MADD_HX4X2(S20, S30, W10, W10, Y21, Y22);
            MADD_HX4X2(S21, S31, W00, W00, Y31, Y32);

            ADD_HX4X2(S20, S30, S20, S30, S21, S31);

            AE_SAHX4X2_IP(S00, S10, alZ, pZ);
            AE_SAHX4X2_IP(S20, S30, alZ1, pZ1);

            Y20 = Y22;
            Y30 = Y32;
        }

        AE_SA128POS_FP(alZ, pZ);
        AE_SA128POS_FP(alZ1, pZ1);

        //last 2
        {

            Y21 = CONST_HX4(0);
            Y31 = CONST_HX4(0);

            MUL_HX4X2( S00, S20, W02, W12, AE_SELH_5432(Y20, Y21), AE_SELH_5432(Y20, Y21));
            MADD_HX4X2(S00, S20, W01, W11, AE_SELH_4321(Y20, Y21), AE_SELH_4321(Y20, Y21));

            MADD_HX4(S20, W02, AE_SELH_5432(Y30, Y31));
            MADD_HX4(S20, W01, AE_SELH_4321(Y30, Y31));


            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(S00), AE_MOVF16X4_FROMHALFX4(S20), AE_MOVINT16X4_FROMINT64(0x0602070304000501));
            S00 = AE_MOVHALFX4_FROMF16X4(t0);
            S20 = AE_MOVHALFX4_FROMF16X4(t1);
            tmp = AE_MOVINT32X2_FROMINT16X4(AE_MOVINT16X4_FROMXTHALFX4(S00));
            AE_S32_H_I(tmp, (ae_int32*)pZ, 0);
            tmp = AE_MOVINT32X2_FROMINT16X4(AE_MOVINT16X4_FROMXTHALFX4(S20));
            AE_S32_H_I(tmp, (ae_int32*)pZ1, 0);

        }

    }
    __Pragma("loop_count min=1");
    for (i = 2; i < P; i+=2)
    {
        pY0 = (xthalfx8*)(y + (i - 2) * Q);
        pY1 = (xthalfx8*)((float16_t*)pY0 + Q);
        pY2 = (xthalfx8*)((float16_t*)pY1 + Q);
        pY3 = (xthalfx8*)((float16_t*)pY2 + Q);

        Y00 = CONST_HX4(0);
        Y10 = CONST_HX4(0);
        Y20 = CONST_HX4(0);
        Y30 = CONST_HX4(0);

        pZ  = (xthalfx8*)(z + (i + 0) * (N + Q - 1));
        pZ1 = (xthalfx8*)(z + (i + 1) * (N + Q - 1));

        __Pragma("loop_count min=1");
        for (j = 0; j < Q; j += 8)
        {

            AE_LHX4X2_IP(Y01, Y02, pY0, sizeof(xthalfx8));
            AE_LHX4X2_IP(Y11, Y12, pY1, sizeof(xthalfx8));

            // W22 W21
            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y10), AE_MOVF16X4_FROMHALFX4(Y11), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t0); SEL2 = AE_MOVHALFX4_FROMF16X4(t1);
            AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), sel); SEL1 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
            Y00 = Y02;Y10 = Y12;
            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y11), AE_MOVF16X4_FROMHALFX4(Y12), sel); SEL0_ = AE_MOVHALFX4_FROMF16X4(t0);SEL2_ = AE_MOVHALFX4_FROMF16X4(t1);
            AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), sel); SEL1_ = AE_MOVHALFX4_FROMF16X4(t2);SEL3_ = AE_MOVHALFX4_FROMF16X4(t3);
            AE_LHX4X2_IP(W22, W21, pW, sizeof(xthalfx8));
            MUL_HX4X2(S00, S20, W22, W22, SEL1, SEL0);
            MUL_HX4X2(S01, S21, W22, W22, SEL1_, SEL0_);
            MUL_HX4X2(S10, S30, W21, W21, SEL3, SEL2);
            MUL_HX4X2(S11, S31, W21, W21, SEL3_, SEL2_);

            AE_LHX4X2_IP(Y21, Y22, pY2, sizeof(xthalfx8));
            AE_LHX4X2_IP(Y31, Y32, pY3, sizeof(xthalfx8));
            // W12 W11
            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y20), AE_MOVF16X4_FROMHALFX4(Y21), sel); SEL1 = AE_MOVHALFX4_FROMF16X4(t0);SEL3 = AE_MOVHALFX4_FROMF16X4(t1);
            Y20 = Y22;
            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y21), AE_MOVF16X4_FROMHALFX4(Y22), sel); SEL1_ = AE_MOVHALFX4_FROMF16X4(t0);SEL3_ = AE_MOVHALFX4_FROMF16X4(t1);
            AE_LHX4X2_IP(W12, W11, pW, sizeof(xthalfx8));
            MADD_HX4X2(S00, S20, W12, W12, SEL0, SEL1);
            MADD_HX4X2(S01, S21, W12, W12, SEL0_, SEL1_);
            MADD_HX4X2(S10, S30, W11, W11, SEL2, SEL3);
            MADD_HX4X2(S11, S31, W11, W11, SEL2_, SEL3_);
            //W02 W01
            AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y30), AE_MOVF16X4_FROMHALFX4(Y31), sel);SEL0 = AE_MOVHALFX4_FROMF16X4(t2);SEL2 = AE_MOVHALFX4_FROMF16X4(t3);
            Y30 = Y32;
            AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y31), AE_MOVF16X4_FROMHALFX4(Y32), sel);SEL0_ = AE_MOVHALFX4_FROMF16X4(t2);SEL2_ = AE_MOVHALFX4_FROMF16X4(t3);
            AE_LHX4X2_IP(W02, W01, pW, sizeof(xthalfx8));
            MADD_HX4X2(S00, S20, W02, W02, SEL1, SEL0);
            MADD_HX4X2(S01, S21, W02, W02, SEL1_, SEL0_);
            MADD_HX4X2(S10, S30, W01, W01, SEL3, SEL2);
            MADD_HX4X2(S11, S31, W01, W01, SEL3_, SEL2_);

            //w20 w10
            AE_LHX4X2_IP(W20, W10, pW, sizeof(xthalfx8));
            MADD_HX4X2(S00, S20, W20, W20, Y01, Y11);
            MADD_HX4X2(S01, S21, W20, W20, Y02, Y12);
            MADD_HX4X2(S10, S30, W10, W10, Y11, Y21);
            MADD_HX4X2(S11, S31, W10, W10, Y12, Y22);

            //w00
            AE_LHX4XP(W00, castxcc(xthalfx4, pW), -4 * (int)sizeof(xthalfx8));
            MADD_HX4X2(S00, S20, W00, W00, Y21, Y31);
            MADD_HX4X2(S01, S21, W00, W00, Y22, Y32);

            ADD_HX4X2(S00, S20, S00, S20, S10, S30);
            ADD_HX4X2(S01, S21, S01, S21, S11, S31);

            AE_SAHX4X2_IP(S00, S01, alZ, pZ);
            AE_SAHX4X2_IP(S20, S21, alZ1, pZ1);
        }
        AE_LHX4X2_IP(W22, W21, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W12, W11, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W02, W01, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W20, W10, pW, -3*(int)sizeof(xthalfx8));

        AE_SA128POS_FP(alZ, pZ);
        AE_SA128POS_FP(alZ1, pZ1);
        //last 2 
        {
            Y01 = CONST_HX4(0);
            Y11 = CONST_HX4(0);
            Y21 = CONST_HX4(0);
            Y31 = CONST_HX4(0);

            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y10), AE_MOVF16X4_FROMHALFX4(Y11), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t0); SEL2 = AE_MOVHALFX4_FROMF16X4(t1);
            AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y20), AE_MOVF16X4_FROMHALFX4(Y21), sel); SEL1 = AE_MOVHALFX4_FROMF16X4(t2); SEL3 = AE_MOVHALFX4_FROMF16X4(t3);
            AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), sel); SEL0_ = AE_MOVHALFX4_FROMF16X4(t2); SEL2_ = AE_MOVHALFX4_FROMF16X4(t3);
            AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(Y30), AE_MOVF16X4_FROMHALFX4(Y31), sel); SEL1_ = AE_MOVHALFX4_FROMF16X4(t2); SEL3_ = AE_MOVHALFX4_FROMF16X4(t3);

            MUL_HX4X2( S00, S20, W22, W22, SEL0_, SEL0);
            MUL_HX4X2( S10, S30, W21, W21, SEL2_, SEL2);
            MADD_HX4X2(S00, S20, W12, W12, SEL0, SEL1);
            MADD_HX4X2(S10, S30, W11, W11, SEL2, SEL3);

            MADD_HX4X2(S00, S20, W02, W02, SEL1, SEL1_);
            MADD_HX4X2(S10, S30, W01, W01, SEL3, SEL3_);

            ADD_HX4X2(S00, S20, S00, S20, S10, S30);

            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(S00), AE_MOVF16X4_FROMHALFX4(S20), AE_MOVINT16X4_FROMINT64(0x0602070304000501)); 
            S00 = AE_MOVHALFX4_FROMF16X4(t0); 
            S20 = AE_MOVHALFX4_FROMF16X4(t1);

            tmp = AE_MOVINT32X2_FROMINT16X4(AE_MOVINT16X4_FROMXTHALFX4(S00));
            AE_S32_H_I(tmp, (ae_int32*)pZ, 0);
            tmp = AE_MOVINT32X2_FROMINT16X4(AE_MOVINT16X4_FROMXTHALFX4(S20));
            AE_S32_H_I(tmp, (ae_int32*)pZ1, 0);
        }
    }
    //last 2 rows

    {
        pY0 = (xthalfx8*)(y + (P - 2) * Q);
        pY1 = (xthalfx8*)((float16_t*)pY0 + Q);


        Y00 = CONST_HX4(0);
        Y10 = CONST_HX4(0);

        pZ = (xthalfx8*)(z + (P) * (N + Q - 1));
        pZ1 = (xthalfx8*)(z + (P + 1) * (N + Q - 1));

        __Pragma("loop_count min=1");
        for (j = 0; j < Q; j += 8)
        {

            AE_LHX4X2_IP(Y01, Y02, pY0, sizeof(xthalfx8));
            AE_LHX4X2_IP(Y11, Y12, pY1, sizeof(xthalfx8));

            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y10), AE_MOVF16X4_FROMHALFX4(Y11), sel); SEL0 = AE_MOVHALFX4_FROMF16X4(t0); SEL2 = AE_MOVHALFX4_FROMF16X4(t1);
            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y11), AE_MOVF16X4_FROMHALFX4(Y12), sel); SEL1 = AE_MOVHALFX4_FROMF16X4(t0); SEL3 = AE_MOVHALFX4_FROMF16X4(t1);

            MUL_HX4X2(S00, S10, W12, W12, SEL0, SEL1);
            MUL_HX4X2(S01, S11, W11, W11, SEL2, SEL3);
            MUL_HX4X2(S20, S30, W22, W22, SEL0, SEL1);
            MADD_HX4X2(S20, S30, W21, W21, SEL2, SEL3);

            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), sel); SEL0_ = AE_MOVHALFX4_FROMF16X4(t0); SEL2_ = AE_MOVHALFX4_FROMF16X4(t1);
            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), sel); SEL1_ = AE_MOVHALFX4_FROMF16X4(t0); SEL3_ = AE_MOVHALFX4_FROMF16X4(t1);

            MADD_HX4X2(S00, S10, W22, W22, SEL0_, SEL1_);
            MADD_HX4X2(S01, S11, W21, W21, SEL2_, SEL3_);
           
            MADD_HX4X2(S00, S10, W20, W20, Y01, Y02);
            MADD_HX4X2(S01, S11, W10, W10, Y11, Y12);

            MADD_HX4X2(S20, S30, W20, W20, Y11, Y12);

            ADD_HX4X2(S00, S10, S00, S10, S01, S11);

            AE_SAHX4X2_IP(S00, S10, alZ, pZ);
            AE_SAHX4X2_IP(S20, S30, alZ1, pZ1);

            Y00 = Y02;
            Y10 = Y12;

        }

        AE_SA128POS_FP(alZ, pZ);
        AE_SA128POS_FP(alZ1, pZ1);


        {

            Y01 = CONST_HX4(0);
            Y11 = CONST_HX4(0);

            MUL_HX4X2 (S00, S20, W22, W22, AE_SELH_5432(Y00, Y01), AE_SELH_5432(Y10, Y11));
            MADD_HX4X2(S00, S20, W21, W21, AE_SELH_4321(Y00, Y01), AE_SELH_4321(Y10, Y11));

            MADD_HX4(S00, W12, AE_SELH_5432(Y10, Y11));
            MADD_HX4(S00, W11, AE_SELH_4321(Y10, Y11));

            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(S00), AE_MOVF16X4_FROMHALFX4(S20), AE_MOVINT16X4_FROMINT64(0x0602070304000501));
            S00 = AE_MOVHALFX4_FROMF16X4(t0);
            S20 = AE_MOVHALFX4_FROMF16X4(t1);
            tmp = AE_MOVINT32X2_FROMINT16X4(AE_MOVINT16X4_FROMXTHALFX4(S00));
            AE_S32_H_I(tmp, (ae_int32*)pZ, 0);
            tmp = AE_MOVINT32X2_FROMINT16X4(AE_MOVINT16X4_FROMXTHALFX4(S20));
            AE_S32_H_I(tmp, (ae_int32*)pZ1, 0);
        }

    }
#   undef M
#   undef N
}

size_t conv2d_3x3_fp16_getScratchSize(int P, int Q)
{
	return 9*sizeof(xthalfx4);
} 
#elif 0 // code using MULCNV
#include "float16.h"
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )
#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
void conv2d_3x3_fp16(void* pScr, float16_t* z, const float16_t* x, const float16_t* y, int P, int Q)
{
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT_ALIGN(z, 16);
    NASSERT_ALIGN(pScr, 16);
    int i, j;
    int M = 3, N = 3;
    if (N <= 0 || M <= 0 || P <= 0 || Q <= 0) return;
    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    int stride = ((Q+4) + 7) & ~7;

    const xthalfx8* restrict pY0;
    const xthalfx8* restrict pY1;
    const xthalfx8* restrict pY2;
    xthalfx8* restrict  pZ = (xthalfx8*)z;
    ae_valignx2 alY0, alY1, alY2, alZ;
    alZ = AE_ZALIGN128();
    
    ae_int16x4 zero = 0;
    ae_int16x8* pS = (ae_int16x8*)pScr;
    ae_int16x8* pR = (ae_int16x8*)y;
    //zero first 2 row
    for (j = 0; j < (stride>>2); j++)
    {
        AE_S16X4X2_IP(zero, zero, pS, sizeof(ae_int16x8));
    }
    for (i = 0; i < P; i++)
    {
        AE_S16_0_IP(zero, castxcc(ae_int16, pS), sizeof(ae_int16));
        AE_S16_0_IP(zero, castxcc(ae_int16, pS), sizeof(ae_int16));
        alY0 = AE_LA128_PP(pR);
        for (j = 0; j < (Q & ~7); j+=8)
        {
            ae_int16x4 y0, y1;
            AE_LA16X4X2_IP(y0, y1, alY0, pR);
            AE_SA16X4X2_IP(y0, y1, alZ, pS);
        }
        AE_SA128POS_FP(alZ, pS);
        for (j = 0; j < (Q&7); j++)
        {
            ae_int16x4 y0;
            AE_L16_IP(y0, castxcc(ae_int16,pR), sizeof(ae_int16));
            AE_S16_0_IP(y0, castxcc(ae_int16, pS), sizeof(ae_int16));
        }
        for (j = Q+2; j < stride; j++)
        {
            AE_S16_0_IP(zero, castxcc(ae_int16, pS), sizeof(ae_int16));
        }
    }
    //zero last 2 row
    for (j = 0; j < (stride >> 2); j++)
    {
        AE_S16X4X2_IP(zero, zero, pS, sizeof(ae_int16x8));
    }

    xthalfx4 w0, w1, w2;
    {
        ae_int16x4 t0, t1, t2;
        ae_int16x4 mask0;
        ae_int16x4 x0123, x4567, x8;
        AE_L16X4X2_IP(x0123, x4567, castxcc(ae_int16x8, x), sizeof(xthalfx8));
        AE_L16_IP(x8, castxcc(ae_int16, x), 0);
        t0 = AE_SEL16X4(x8, x4567,    AE_MOVF16X4_FROMF64(0x0004000000010000)); // 4010
        t1 = AE_SEL16X4(x4567, x0123, AE_MOVF16X4_FROMF64(0x0006000700000000)); // 6700
        t2 = AE_SEL16X4(x0123, x0123, AE_MOVF16X4_FROMF64(0x0001000200030000)); // 1230
        mask0 = AE_MOVF16X4_FROMINT64(0xFFFFFFFFFFFF0000);
        w0 = AE_MOVHALFX4_FROMF16X4(AE_AND16(t0, mask0));
        w1 = AE_MOVHALFX4_FROMF16X4(AE_AND16(t1, mask0));
        w2 = AE_MOVHALFX4_FROMF16X4(AE_AND16(t2, mask0));
    }


    ae_int16x4 sel = AE_MOVINT16X4_FROMINT64(0x0706050403020100); // SEL7531 + SEL6420

    for (i = 0; i < M + P - 1; i++)
    {
        xthalfx4 y00, y01, y02, y10, y11, y12, y20, y21, y22;
        pY0 = (xthalfx8*)((uintptr_t)pScr + sizeof(xthalf) * stride * i      );
        pY1 = (xthalfx8*)((uintptr_t)pScr + sizeof(xthalf) * stride * (i + 1));
        pY2 = (xthalfx8*)((uintptr_t)pScr + sizeof(xthalf) * stride * (i + 2));
        
        AE_LHX4IP(y00, castxcc(xthalfx4,pY0), sizeof(xthalfx4));
        AE_LHX4IP(y10, castxcc(xthalfx4,pY1), sizeof(xthalfx4));
        AE_LHX4IP(y20, castxcc(xthalfx4,pY2), sizeof(xthalfx4));

        alY0 = AE_LA128_PP(pY0);
        alY1 = AE_LA128_PP(pY1);
        alY2 = AE_LA128_PP(pY2);

        for (j = 0; j < (Q + 2)/8; j++)
        {
            xthalfx4 z0, z1, z2, z3;
            xthalfx4 r0, r1, r2, r3;
            ae_int16x4 t0, t1, t2, t3;
            ae_int16x4 k0, k1, k2, k3;

            AE_LAHX4X2_IP(y01, y02, alY0, pY0);
            AE_LAHX4X2_IP(y11, y12, alY1, pY1);
            AE_LAHX4X2_IP(y21, y22, alY2, pY2);

            MULCNVH_HX4X2(z0, z1, w0, y00, y01);
            MULCNVL_HX4X2(z2, z3, w0, y00, y01);
            MULACNVH_HX4X2(z0, z1, w1, y10, y11);
            MULACNVL_HX4X2(z2, z3, w1, y10, y11);
            MULACNVH_HX4X2(z0, z1, w2, y20, y21);
            MULACNVL_HX4X2(z2, z3, w2, y20, y21);

            MULCNVH_HX4X2(r0, r1, w0, y01, y02);
            MULCNVL_HX4X2(r2, r3, w0, y01, y02);
            MULACNVH_HX4X2(r0, r1, w1, y11, y12);
            MULACNVL_HX4X2(r2, r3, w1, y11, y12);
            MULACNVH_HX4X2(r0, r1, w2, y21, y22);
            MULACNVL_HX4X2(r2, r3, w2, y21, y22);

            y00 = y02;
            y10 = y12;
            y20 = y22;

            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
            AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
            AE_DSEL16X4(t0, t2, t0, t2, sel);
            AE_DSEL16X4(t1, t3, t1, t3, sel);
            z0 = AE_MOVHALFX4_FROMF16X4(t0);
            z1 = AE_MOVHALFX4_FROMF16X4(t2);
            z2 = AE_MOVHALFX4_FROMF16X4(t1);

            AE_DSEL16X4(k0, k1, AE_MOVF16X4_FROMHALFX4(r0), AE_MOVF16X4_FROMHALFX4(r1), sel);
            AE_DSEL16X4(k2, k3, AE_MOVF16X4_FROMHALFX4(r2), AE_MOVF16X4_FROMHALFX4(r3), sel);
            AE_DSEL16X4(k0, k2, k0, k2, sel);
            AE_DSEL16X4(k1, k3, k1, k3, sel);
            r0 = AE_MOVHALFX4_FROMF16X4(k0);
            r1 = AE_MOVHALFX4_FROMF16X4(k2);
            r2 = AE_MOVHALFX4_FROMF16X4(k1);

            ADD_HX4X2(z0, r0, z0, r0, z1, r1);
            ADD_HX4X2(z0, r0, z0, r0, z2, r2);
            AE_SAHX4X2_IP(z0, r0, alZ, pZ);

        }
        if ((Q + 2) & 7) 
        {
            xthalfx4 z0, z1, z2, z3;
            xthalfx4 r0, r1, r2, r3;
            ae_int16x4 t0, t1, t2, t3;
            ae_int16x4 k0, k1, k2, k3;

            AE_LAHX4X2_IP(y01, y02, alY0, pY0);
            AE_LAHX4X2_IP(y11, y12, alY1, pY1);
            AE_LAHX4X2_IP(y21, y22, alY2, pY2);

            MULCNVH_HX4X2(z0, z1, w0, y00, y01);
            MULCNVL_HX4X2(z2, z3, w0, y00, y01);
            MULACNVH_HX4X2(z0, z1, w1, y10, y11);
            MULACNVL_HX4X2(z2, z3, w1, y10, y11);
            MULACNVH_HX4X2(z0, z1, w2, y20, y21);
            MULACNVL_HX4X2(z2, z3, w2, y20, y21);

            MULCNVH_HX4X2(r0, r1, w0, y01, y02);
            MULCNVL_HX4X2(r2, r3, w0, y01, y02);
            MULACNVH_HX4X2(r0, r1, w1, y11, y12);
            MULACNVL_HX4X2(r2, r3, w1, y11, y12);
            MULACNVH_HX4X2(r0, r1, w2, y21, y22);
            MULACNVL_HX4X2(r2, r3, w2, y21, y22);

            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
            AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
            AE_DSEL16X4(t0, t2, t0, t2, sel);
            AE_DSEL16X4(t1, t3, t1, t3, sel);
            z0 = AE_MOVHALFX4_FROMF16X4(t0);
            z1 = AE_MOVHALFX4_FROMF16X4(t2);
            z2 = AE_MOVHALFX4_FROMF16X4(t1);

            AE_DSEL16X4(k0, k1, AE_MOVF16X4_FROMHALFX4(r0), AE_MOVF16X4_FROMHALFX4(r1), sel);
            AE_DSEL16X4(k2, k3, AE_MOVF16X4_FROMHALFX4(r2), AE_MOVF16X4_FROMHALFX4(r3), sel);
            AE_DSEL16X4(k0, k2, k0, k2, sel);
            AE_DSEL16X4(k1, k3, k1, k3, sel);
            r0 = AE_MOVHALFX4_FROMF16X4(k0);
            r1 = AE_MOVHALFX4_FROMF16X4(k2);
            r2 = AE_MOVHALFX4_FROMF16X4(k1);

            ADD_HX4X2(z0, r0, z0, r0, z1, r1);
            ADD_HX4X2(z0, r0, z0, r0, z2, r2);
            AE_SAVHX4X2_XP(z0, r0, alZ, pZ, ((Q + 2) & 7)*sizeof(xthalf));
        }
    }
    AE_SA128POS_FP(alZ, pZ);
}
size_t conv2d_3x3_fp16_getScratchSize(int P, int Q)
{
    return (P+4)*(((Q+4) + 7) & ~7)*sizeof(float16_t);
}

#elif HAVE_SFPU & HAVE_HPFPU
#include "float16.h"
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )
#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
void conv2d_3x3_fp16(void* pScr, float16_t* z, const float16_t* x, const float16_t* y, int P, int Q)
{
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT_ALIGN(z, 16);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0 && P % 8 == 0);
    NASSERT(Q >= 0 && Q % 8 == 0);
    conv2d_fp16(pScr, z, x, y, 3, 3, P, Q);
}
void conv2d_fp16(void* pScr, float16_t* z, const float16_t* x, const float16_t* y, int M, int N, int P, int Q)
{
    int i, j, m, n, n0, n1, m0, m1;
    if (N <= 0 || M <= 0 || P <= 0 || Q <= 0) return;
    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    for (i = 0; i < M + P - 1; i++)
        for (j = 0; j < N + Q - 1; j++)
        {
            float16_t s;
            m0 = MAX(i - P + 1, 0);
            m1 = MIN(i + 1, M);
            n0 = MAX(j - Q + 1, 0);
            n1 = MIN(j + 1, N);
            s = 0;
            for (n = n0; n < n1; n++)
                for (m = m0; m < m1; m++)
                {
                    s = fma_f16(x[m * N + n], y[(i - m) * Q + (j - n)], s);
                }
            z[i * (N + Q - 1) + j] = s;
        }
}/* conv2df() */

size_t conv2d_3x3_fp16_getScratchSize(int M, int N, int P, int Q)
{
    return 0;
}

#else
DISCARD_FUN(void,conv2d_3x3_fp16,(void* pScr, float16_t *z, const float16_t * x, const float16_t * y, int P, int Q))
size_t conv2d_3x3_fp16_getScratchSize (int P, int Q) 
{ 
    (void)P,(void)Q;
    return 0;
}
#endif

