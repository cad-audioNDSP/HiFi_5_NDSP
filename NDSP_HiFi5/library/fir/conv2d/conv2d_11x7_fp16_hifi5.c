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

#if !HAVE_HPFPU
DISCARD_FUN(void,conv2d_11x7_fp16,(void* pScr, float16_t *z, const float16_t * x, const float16_t * y, int P, int Q))
size_t conv2d_11x7_fp16_getScratchSize (int P, int Q) 
{ 
    (void)P,(void)Q;
    return 0;
}
#elif HAVE_HPFPU
void conv2d_11x7_fp16(void* pScr, float16_t* z, const float16_t* x, const float16_t* y, int P, int Q)
{
#   define M 11
#   define N 7
    int i, j, m, mmin, mmax, ystart, xstart;


    const xthalfx8* restrict pY;
    const xthalfx8* restrict pY0;

    // weights
    const xthalfx8* restrict pW;
    
    // rezult
    xthalfx8* restrict pZ;

    // temp scratch mem
    xthalfx8* restrict pT;
    const xthalfx8* restrict pTread;

    xthalfx4 Y00, Y01, Y02, Y03;
    xthalfx4 Y10, Y11, Y12, Y13;
    xthalfx4 Y20, Y21, Y22, Y23;


    xthalfx4 z0, z1, z2, z3;
    xthalfx4 r0, r1, r2, r3;
    ae_int16x4 t0, t1, t2, t3;
    ae_int16x4 k0, k1, k2, k3;
    int stride = Q * sizeof(xthalf);
    xthalfx4 W0123, W456z;
    ae_int16x4 sel = AE_MOVINT16X4_FROMINT64(0x0706050403020100); // SEL7531 + SEL6420

    ae_valignx2 alZ;

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

    //load X and replicate, prepare zero row 
    {
        xthalfx4 W0, W1, W2, W3, W4, W5, W6;
        xthalfx4* pWw = (xthalfx4*)pScr+2*M - 2;
        const xthalfx4* pX;
        pX = (xthalfx4*)(x);
        for (i = 0; i < 11; ++i)
        {
            ae_int16x4 v_i;
            xthalfx4 w0123, w456z;

            AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(ae_int16));
            W6 = AE_MOVXTHALFX4_FROMINT16X4(v_i);

            AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(ae_int16));
            W5 = AE_MOVXTHALFX4_FROMINT16X4(v_i);

            AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(ae_int16));
            W4 = AE_MOVXTHALFX4_FROMINT16X4(v_i);

            AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(ae_int16));
            W3 = AE_MOVXTHALFX4_FROMINT16X4(v_i);

            AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(ae_int16));
            W2 = AE_MOVXTHALFX4_FROMINT16X4(v_i);

            AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(ae_int16));
            W1 = AE_MOVXTHALFX4_FROMINT16X4(v_i);

            AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(ae_int16));
            W0 = AE_MOVXTHALFX4_FROMINT16X4(v_i);

            w0123 = AE_SELH_7362(W0, W1);
            w0123 = AE_SELH_7610(w0123, W2);
            w0123 = AE_MOVHALFX4_FROMF16X4(AE_SEL16X4(AE_MOVF16X4_FROMHALFX4(w0123),
                AE_MOVF16X4_FROMHALFX4(W3),
                AE_MOVF16X4_FROMINT64(0x0007000600050003)));

            w456z = AE_SELH_7362(W4, W5);
            w456z = AE_SELH_7610(w456z, W6);
            w456z = AE_MOVHALFX4_FROMF16X4(AE_SEL16X4(AE_MOVF16X4_FROMHALFX4(w456z),
                AE_MOVF16X4_FROMHALFX4(CONST_HX4(0)),
                AE_MOVF16X4_FROMINT64(0x0007000600050003)));

            AE_SHX4X2_IP(w0123, w456z, castxcc(xthalfx8, pWw), -(int)sizeof(xthalfx8));

        }

        pWw = (xthalfx4*)(pScr)+ 2*M;
        for (i = 0; i < Q >> 2; ++i)
        {
            AE_SHX4IP(CONST_HX4(0), pWw, 8);
        }
    }

    alZ = AE_ZALIGN128();

    // Processing of convolution

    __Pragma("loop_count min=1");
    for (i = 0; i < P + M - 1; i += 1)
    {

        mmin = i < M - 1 ? M - 1 - i: 0;
        mmax = i >= P ? M - 1 + P - i : M;
        xstart = i < M - 1 ? (mmin) : 0;
        ystart = i < M - 1 ? 0 : i - (M - 1);

        pW = (xthalfx8*)pScr + xstart;
        pY = (xthalfx8*)(y + ystart * Q);

        // if only one Y row need to process
        if (mmax-mmin == 1)
        {
            AE_LHX4X2_IP(W0123, W456z, pW, sizeof(xthalfx8));
            //setup pointers and init variables
            CONST_HX4X2(Y00, Y01, 0);
            pY0 = pY;
            pZ = (xthalfx8*)(z + (i + 0) * (N + Q - 1));

            __Pragma("loop_count min=1");
            for (j = 0; j < Q; j += 8)
            {
                AE_LHX4X2_IP(Y02, Y03, pY0, sizeof(xthalfx8));

                MULCNVL_HX4X2(z0, z1, W0123, Y00, Y01);
                MULCNVH_HX4X2(z2, z3, W0123, Y01, Y02);
                MULACNVL_HX4X2(z0, z1, W456z, Y01, Y02);
                MULACNVH_HX4X2(z2, z3, W456z, Y02, Y02);

                MULCNVL_HX4X2(r0, r1, W0123, Y01, Y02);
                MULCNVH_HX4X2(r2, r3, W0123, Y02, Y03);
                MULACNVL_HX4X2(r0, r1, W456z, Y02, Y03);
                MULACNVH_HX4X2(r2, r3, W456z, Y03, Y03);

                Y00 = Y02; Y01 = Y03;

                AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                AE_DSEL16X4(t0, t2, t0, t2, sel);
                AE_DSEL16X4(t1, t3, t1, t3, sel);
                z0 = AE_MOVHALFX4_FROMF16X4(t0);
                z1 = AE_MOVHALFX4_FROMF16X4(t2);
                z2 = AE_MOVHALFX4_FROMF16X4(t1);
                z3 = AE_MOVHALFX4_FROMF16X4(t3);

                AE_DSEL16X4(k0, k1, AE_MOVF16X4_FROMHALFX4(r0), AE_MOVF16X4_FROMHALFX4(r1), sel);
                AE_DSEL16X4(k2, k3, AE_MOVF16X4_FROMHALFX4(r2), AE_MOVF16X4_FROMHALFX4(r3), sel);
                AE_DSEL16X4(k0, k2, k0, k2, sel);
                AE_DSEL16X4(k1, k3, k1, k3, sel);
                r0 = AE_MOVHALFX4_FROMF16X4(k0);
                r1 = AE_MOVHALFX4_FROMF16X4(k2);
                r2 = AE_MOVHALFX4_FROMF16X4(k1);
                r3 = AE_MOVHALFX4_FROMF16X4(k3);

                ADD_HX4X2(z0, z1, z0, z1, z2, z3);
                ADD_HX4X2(r0, r1, r0, r1, r2, r3);
                ADD_HX4X2(z0, r0, z0, r0, z1, r1);

                AE_SAHX4X2_IP(z0, r0, alZ, pZ);

            }
            //last 6
            {
                CONST_HX4X2(Y02, Y03, 0);

                MULCNVL_HX4X2(z0, z1, W0123, Y00, Y01);
                MULCNVH_HX4X2(z2, z3, W0123, Y01, Y02);
                MULACNVL_HX4X2(z0, z1, W456z, Y01, Y02);
                MULACNVH_HX4X2(z2, z3, W456z, Y02, Y02);

                MULCNVL_HX4X2(r0, r1, W0123, Y01, Y02);
                MULCNVH_HX4X2(r2, r3, W0123, Y02, Y03);
                MULACNVL_HX4X2(r0, r1, W456z, Y02, Y03);
                MULACNVH_HX4X2(r2, r3, W456z, Y03, Y03);

                Y00 = Y02; Y01 = Y03;

                AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                AE_DSEL16X4(t0, t2, t0, t2, sel);
                AE_DSEL16X4(t1, t3, t1, t3, sel);
                z0 = AE_MOVHALFX4_FROMF16X4(t0);
                z1 = AE_MOVHALFX4_FROMF16X4(t2);
                z2 = AE_MOVHALFX4_FROMF16X4(t1);
                z3 = AE_MOVHALFX4_FROMF16X4(t3);

                AE_DSEL16X4(k0, k1, AE_MOVF16X4_FROMHALFX4(r0), AE_MOVF16X4_FROMHALFX4(r1), sel);
                AE_DSEL16X4(k2, k3, AE_MOVF16X4_FROMHALFX4(r2), AE_MOVF16X4_FROMHALFX4(r3), sel);
                AE_DSEL16X4(k0, k2, k0, k2, sel);
                AE_DSEL16X4(k1, k3, k1, k3, sel);
                r0 = AE_MOVHALFX4_FROMF16X4(k0);
                r1 = AE_MOVHALFX4_FROMF16X4(k2);
                r2 = AE_MOVHALFX4_FROMF16X4(k1);
                r3 = AE_MOVHALFX4_FROMF16X4(k3);

                ADD_HX4X2(z0, z1, z0, z1, z2, z3);
                ADD_HX4X2(r0, r1, r0, r1, r2, r3);
                ADD_HX4X2(z0, r0, z0, r0, z1, r1);

                AE_SAVHX4X2_XP(z0, r0, alZ, pZ, 6 * sizeof(xthalf));
            }
            AE_SA128POS_FP(alZ, pZ);
        }
        else
        {
            // if not all 11 rows of Y need to process
            if (mmax-mmin!=M)
            {
                //first itteration
                {
                    //setup pointers and init variables
                    CONST_HX4X2(Y00, Y01, 0);
                    pY0 = pY;
                    pY = (xthalfx8*)((float16_t*)pY + Q);
                    pT = (xthalfx8*)((xthalfx4*)(pScr)+2 * M + (Q / 2));

                    AE_LHX4X2_IP(W0123, W456z, pW, sizeof(xthalfx8));
                    __Pragma("loop_count min=1");
                    for (j = 0; j < Q; j += 8)
                    {
                        AE_LHX4X2_IP(Y02, Y03, pY0, sizeof(xthalfx8));

                        MULCNVL_HX4X2(z0, z1, W0123, Y00, Y01);
                        MULCNVH_HX4X2(z2, z3, W0123, Y01, Y02);
                        MULACNVL_HX4X2(z0, z1, W456z, Y01, Y02);
                        MULACNVH_HX4X2(z2, z3, W456z, Y02, Y02);

                        MULCNVL_HX4X2(r0, r1, W0123, Y01, Y02);
                        MULCNVH_HX4X2(r2, r3, W0123, Y02, Y03);
                        MULACNVL_HX4X2(r0, r1, W456z, Y02, Y03);
                        MULACNVH_HX4X2(r2, r3, W456z, Y03, Y03);

                        Y00 = Y02; Y01 = Y03;

                        AE_SHX4X2_IP(z0, z1, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r0, r1, pT, sizeof(xthalfx8));

                        AE_SHX4X2_IP(z2, z3, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r2, r3, pT, sizeof(xthalfx8));
                    }
                    //last 6
                    {
                        CONST_HX4X2(Y02, Y03, 0);

                        MULCNVL_HX4X2(z0, z1, W0123, Y00, Y01);
                        MULCNVH_HX4X2(z2, z3, W0123, Y01, Y02);
                        MULACNVL_HX4X2(z0, z1, W456z, Y01, Y02);
                        MULACNVH_HX4X2(z2, z3, W456z, Y02, Y02);

                        MULCNVL_HX4X2(r0, r1, W0123, Y01, Y02);
                        MULCNVH_HX4X2(r2, r3, W0123, Y02, Y03);
                        MULACNVL_HX4X2(r0, r1, W456z, Y02, Y03);
                        MULACNVH_HX4X2(r2, r3, W456z, Y03, Y03);

                        Y00 = Y02; Y01 = Y03;

                        AE_SHX4X2_IP(z0, z1, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r0, r1, pT, sizeof(xthalfx8));

                        AE_SHX4X2_IP(z2, z3, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r2, r3, pT, sizeof(xthalfx8));
                    }

                }
                //main inner loops
                __Pragma("loop_count min=0, max=9");
                for (m = mmin + 1; m < mmax - 1; ++m)
                {
                    //setup pointers and init variables
                    CONST_HX4X2(Y00, Y01, 0);
                    pY0 = pY;
                    pY = (xthalfx8*)((float16_t*)pY + Q);

                    pT = (xthalfx8*)((xthalfx4*)(pScr)+2 * M + (Q / 2));
                    pTread = pT;

                    AE_LHX4X2_IP(W0123, W456z, pW, sizeof(xthalfx8));
                    for (j = 0; j < Q; j += 8)
                    {
                        AE_LHX4X2_IP(Y02, Y03, pY0, sizeof(xthalfx8));

                        AE_LHX4X2_IP(z0, z1, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r0, r1, pTread, sizeof(xthalfx8));

                        AE_LHX4X2_IP(z2, z3, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r2, r3, pTread, sizeof(xthalfx8));

                        MULACNVL_HX4X2(z0, z1, W0123, Y00, Y01);
                        MULACNVH_HX4X2(z2, z3, W0123, Y01, Y02);
                        MULACNVL_HX4X2(z0, z1, W456z, Y01, Y02);
                        MULACNVH_HX4X2(z2, z3, W456z, Y02, Y02);

                        MULACNVL_HX4X2(r0, r1, W0123, Y01, Y02);
                        MULACNVH_HX4X2(r2, r3, W0123, Y02, Y03);
                        MULACNVL_HX4X2(r0, r1, W456z, Y02, Y03);
                        MULACNVH_HX4X2(r2, r3, W456z, Y03, Y03);

                        Y00 = Y02; Y01 = Y03;

                        AE_SHX4X2_IP(z0, z1, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r0, r1, pT, sizeof(xthalfx8));

                        AE_SHX4X2_IP(z2, z3, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r2, r3, pT, sizeof(xthalfx8));
                    }

                    //last 6
                    {
                        CONST_HX4X2(Y02, Y03, 0);

                        AE_LHX4X2_IP(z0, z1, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r0, r1, pTread, sizeof(xthalfx8));

                        AE_LHX4X2_IP(z2, z3, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r2, r3, pTread, sizeof(xthalfx8));

                        MULACNVL_HX4X2(z0, z1, W0123, Y00, Y01);
                        MULACNVH_HX4X2(z2, z3, W0123, Y01, Y02);
                        MULACNVL_HX4X2(z0, z1, W456z, Y01, Y02);
                        MULACNVH_HX4X2(z2, z3, W456z, Y02, Y02);

                        MULACNVL_HX4X2(r0, r1, W0123, Y01, Y02);
                        MULACNVH_HX4X2(r2, r3, W0123, Y02, Y03);
                        MULACNVL_HX4X2(r0, r1, W456z, Y02, Y03);
                        MULACNVH_HX4X2(r2, r3, W456z, Y03, Y03);

                        Y00 = Y02; Y01 = Y03;

                        AE_SHX4X2_IP(z0, z1, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r0, r1, pT, sizeof(xthalfx8));

                        AE_SHX4X2_IP(z2, z3, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r2, r3, pT, sizeof(xthalfx8));
                    }


                }

                // last itteration
                {

                    AE_LHX4X2_IP(W0123, W456z, pW, sizeof(xthalfx8));

                    //setup pointers and init variables
                    CONST_HX4X2(Y00, Y01, 0);
                    pY0 = pY;
                    pZ = (xthalfx8*)(z + (i + 0) * (N + Q - 1));
                    pTread = (xthalfx8*)((xthalfx4*)(pScr)+2 * M + (Q / 2));

                    __Pragma("loop_count min=1");
                    for (j = 0; j < Q; j += 8)
                    {
                        AE_LHX4X2_IP(Y02, Y03, pY0, sizeof(xthalfx8));

                        AE_LHX4X2_IP(z0, z1, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r0, r1, pTread, sizeof(xthalfx8));

                        AE_LHX4X2_IP(z2, z3, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r2, r3, pTread, sizeof(xthalfx8));

                        MULACNVL_HX4X2(z0, z1, W0123, Y00, Y01);
                        MULACNVH_HX4X2(z2, z3, W0123, Y01, Y02);
                        MULACNVL_HX4X2(z0, z1, W456z, Y01, Y02);
                        MULACNVH_HX4X2(z2, z3, W456z, Y02, Y02);

                        MULACNVL_HX4X2(r0, r1, W0123, Y01, Y02);
                        MULACNVH_HX4X2(r2, r3, W0123, Y02, Y03);
                        MULACNVL_HX4X2(r0, r1, W456z, Y02, Y03);
                        MULACNVH_HX4X2(r2, r3, W456z, Y03, Y03);

                        Y00 = Y02; Y01 = Y03;

                        AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                        AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                        AE_DSEL16X4(t0, t2, t0, t2, sel);
                        AE_DSEL16X4(t1, t3, t1, t3, sel);
                        z0 = AE_MOVHALFX4_FROMF16X4(t0);
                        z1 = AE_MOVHALFX4_FROMF16X4(t2);
                        z2 = AE_MOVHALFX4_FROMF16X4(t1);
                        z3 = AE_MOVHALFX4_FROMF16X4(t3);

                        AE_DSEL16X4(k0, k1, AE_MOVF16X4_FROMHALFX4(r0), AE_MOVF16X4_FROMHALFX4(r1), sel);
                        AE_DSEL16X4(k2, k3, AE_MOVF16X4_FROMHALFX4(r2), AE_MOVF16X4_FROMHALFX4(r3), sel);
                        AE_DSEL16X4(k0, k2, k0, k2, sel);
                        AE_DSEL16X4(k1, k3, k1, k3, sel);
                        r0 = AE_MOVHALFX4_FROMF16X4(k0);
                        r1 = AE_MOVHALFX4_FROMF16X4(k2);
                        r2 = AE_MOVHALFX4_FROMF16X4(k1);
                        r3 = AE_MOVHALFX4_FROMF16X4(k3);

                        ADD_HX4X2(z0, z1, z0, z1, z2, z3);
                        ADD_HX4X2(r0, r1, r0, r1, r2, r3);
                        ADD_HX4X2(z0, r0, z0, r0, z1, r1);

                        AE_SAHX4X2_IP(z0, r0, alZ, pZ);

                    }
                    //last 6
                    {
                        CONST_HX4X2(Y02, Y03, 0);

                        AE_LHX4X2_IP(z0, z1, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r0, r1, pTread, sizeof(xthalfx8));

                        AE_LHX4X2_IP(z2, z3, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r2, r3, pTread, sizeof(xthalfx8));

                        MULACNVL_HX4X2(z0, z1, W0123, Y00, Y01);
                        MULACNVH_HX4X2(z2, z3, W0123, Y01, Y02);
                        MULACNVL_HX4X2(z0, z1, W456z, Y01, Y02);
                        MULACNVH_HX4X2(z2, z3, W456z, Y02, Y02);

                        MULACNVL_HX4X2(r0, r1, W0123, Y01, Y02);
                        MULACNVH_HX4X2(r2, r3, W0123, Y02, Y03);
                        MULACNVL_HX4X2(r0, r1, W456z, Y02, Y03);
                        MULACNVH_HX4X2(r2, r3, W456z, Y03, Y03);

                        Y00 = Y02; Y01 = Y03;

                        AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                        AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                        AE_DSEL16X4(t0, t2, t0, t2, sel);
                        AE_DSEL16X4(t1, t3, t1, t3, sel);
                        z0 = AE_MOVHALFX4_FROMF16X4(t0);
                        z1 = AE_MOVHALFX4_FROMF16X4(t2);
                        z2 = AE_MOVHALFX4_FROMF16X4(t1);
                        z3 = AE_MOVHALFX4_FROMF16X4(t3);

                        AE_DSEL16X4(k0, k1, AE_MOVF16X4_FROMHALFX4(r0), AE_MOVF16X4_FROMHALFX4(r1), sel);
                        AE_DSEL16X4(k2, k3, AE_MOVF16X4_FROMHALFX4(r2), AE_MOVF16X4_FROMHALFX4(r3), sel);
                        AE_DSEL16X4(k0, k2, k0, k2, sel);
                        AE_DSEL16X4(k1, k3, k1, k3, sel);
                        r0 = AE_MOVHALFX4_FROMF16X4(k0);
                        r1 = AE_MOVHALFX4_FROMF16X4(k2);
                        r2 = AE_MOVHALFX4_FROMF16X4(k1);
                        r3 = AE_MOVHALFX4_FROMF16X4(k3);

                        ADD_HX4X2(z0, z1, z0, z1, z2, z3);
                        ADD_HX4X2(r0, r1, r0, r1, r2, r3);
                        ADD_HX4X2(z0, r0, z0, r0, z1, r1);

                        AE_SAVHX4X2_XP(z0, r0, alZ, pZ, 6 * sizeof(xthalf));
                    }
                    AE_SA128POS_FP(alZ, pZ);
                }
            }
            else
            {
                // MAIN PARTS OF COMPLEXITY

                //first 3 iter
                {
                    //setup pointers and init variables
                    CONST_HX4X2(Y00, Y01, 0);
                    CONST_HX4X2(Y10, Y11, 0);
                    CONST_HX4X2(Y20, Y21, 0);
                    pY0 = pY;
                    pY = (xthalfx8*)((float16_t*)pY + Q*3);

                    pT = (xthalfx8*)((xthalfx4*)(pScr)+2 * M + (Q / 2));
                    pTread = pT;

                    for (j = 0; j < Q; j += 8)
                    {
                        AE_LHX4X2_IP(W0123, W456z, pW, sizeof(xthalfx8));
                        AE_LHX4X2_XP(Y02, Y03, pY0, stride);
                        MULCNVL_HX4X2(z0, z1, W0123, Y00, Y01);
                        MULCNVH_HX4X2(z2, z3, W0123, Y01, Y02);
                        MULACNVL_HX4X2(z0, z1, W456z, Y01, Y02);
                        MULACNVH_HX4X2(z2, z3, W456z, Y02, Y02);

                        MULCNVL_HX4X2(r0, r1, W0123, Y01, Y02);
                        MULCNVH_HX4X2(r2, r3, W0123, Y02, Y03);
                        MULACNVL_HX4X2(r0, r1, W456z, Y02, Y03);
                        MULACNVH_HX4X2(r2, r3, W456z, Y03, Y03);
                        Y00 = Y02; Y01 = Y03;

                        AE_LHX4X2_IP(W0123, W456z, pW, sizeof(xthalfx8));
                        AE_LHX4X2_XP(Y12, Y13, pY0, stride);
                        MULACNVL_HX4X2(z0, z1, W0123, Y10, Y11);
                        MULACNVH_HX4X2(z2, z3, W0123, Y11, Y12);
                        MULACNVL_HX4X2(z0, z1, W456z, Y11, Y12);
                        MULACNVH_HX4X2(z2, z3, W456z, Y12, Y12);

                        MULACNVL_HX4X2(r0, r1, W0123, Y11, Y12);
                        MULACNVH_HX4X2(r2, r3, W0123, Y12, Y13);
                        MULACNVL_HX4X2(r0, r1, W456z, Y12, Y13);
                        MULACNVH_HX4X2(r2, r3, W456z, Y13, Y13);
                        Y10 = Y12; Y11 = Y13;


                        AE_LHX4X2_XP(W0123, W456z, pW, -2 * (int)sizeof(xthalfx8));
                        AE_LHX4X2_XP(Y22, Y23, pY0, -2 * stride + (int)sizeof(xthalfx8));
                        MULACNVL_HX4X2(z0, z1, W0123, Y20, Y21);
                        MULACNVH_HX4X2(z2, z3, W0123, Y21, Y22);
                        MULACNVL_HX4X2(z0, z1, W456z, Y21, Y22);
                        MULACNVH_HX4X2(z2, z3, W456z, Y22, Y22);

                        MULACNVL_HX4X2(r0, r1, W0123, Y21, Y22);
                        MULACNVH_HX4X2(r2, r3, W0123, Y22, Y23);
                        MULACNVL_HX4X2(r0, r1, W456z, Y22, Y23);
                        MULACNVH_HX4X2(r2, r3, W456z, Y23, Y23);
                        Y20 = Y22; Y21 = Y23;


                        AE_SHX4X2_IP(z0, z1, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r0, r1, pT, sizeof(xthalfx8));

                        AE_SHX4X2_IP(z2, z3, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r2, r3, pT, sizeof(xthalfx8));
                    }

                    //last 6
                    {
                        CONST_HX4X2(Y02, Y03, 0);
                        CONST_HX4X2(Y12, Y13, 0);
                        CONST_HX4X2(Y22, Y23, 0);

                        AE_LHX4X2_IP(W0123, W456z, pW, sizeof(xthalfx8));
                        MULCNVL_HX4X2(z0, z1, W0123, Y00, Y01);
                        MULCNVH_HX4X2(z2, z3, W0123, Y01, Y02);
                        MULACNVL_HX4X2(z0, z1, W456z, Y01, Y02);
                        MULACNVH_HX4X2(z2, z3, W456z, Y02, Y02);

                        MULCNVL_HX4X2(r0, r1, W0123, Y01, Y02);
                        MULCNVH_HX4X2(r2, r3, W0123, Y02, Y03);
                        MULACNVL_HX4X2(r0, r1, W456z, Y02, Y03);
                        MULACNVH_HX4X2(r2, r3, W456z, Y03, Y03);

                        AE_LHX4X2_IP(W0123, W456z, pW, sizeof(xthalfx8));
                        MULACNVL_HX4X2(z0, z1, W0123, Y10, Y11);
                        MULACNVH_HX4X2(z2, z3, W0123, Y11, Y12);
                        MULACNVL_HX4X2(z0, z1, W456z, Y11, Y12);
                        MULACNVH_HX4X2(z2, z3, W456z, Y12, Y12);

                        MULACNVL_HX4X2(r0, r1, W0123, Y11, Y12);
                        MULACNVH_HX4X2(r2, r3, W0123, Y12, Y13);
                        MULACNVL_HX4X2(r0, r1, W456z, Y12, Y13);
                        MULACNVH_HX4X2(r2, r3, W456z, Y13, Y13);

                        AE_LHX4X2_XP(W0123, W456z, pW, -2 * (int)sizeof(xthalfx8));
                        MULACNVL_HX4X2(z0, z1, W0123, Y20, Y21);
                        MULACNVH_HX4X2(z2, z3, W0123, Y21, Y22);
                        MULACNVL_HX4X2(z0, z1, W456z, Y21, Y22);
                        MULACNVH_HX4X2(z2, z3, W456z, Y22, Y22);

                        MULACNVL_HX4X2(r0, r1, W0123, Y21, Y22);
                        MULACNVH_HX4X2(r2, r3, W0123, Y22, Y23);
                        MULACNVL_HX4X2(r0, r1, W456z, Y22, Y23);
                        MULACNVH_HX4X2(r2, r3, W456z, Y23, Y23);

                        AE_SHX4X2_IP(z0, z1, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r0, r1, pT, sizeof(xthalfx8));

                        AE_SHX4X2_IP(z2, z3, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r2, r3, pT, sizeof(xthalfx8));
                    }
                }

                pW += 3;
                //second 3*2 iter
                for (m=0; m<3; ++m)
                {
                    //setup pointers and init variables
                    CONST_HX4X2(Y00, Y01, 0);
                    CONST_HX4X2(Y10, Y11, 0);
                    pY0 = pY;
                    pY = (xthalfx8*)((float16_t*)pY + Q * 2);

                    pT = (xthalfx8*)((xthalfx4*)(pScr)+2 * M + (Q / 2));
                    pTread = pT;

                    for (j = 0; j < Q; j += 8)
                    {
                        AE_LHX4X2_IP(z0, z1, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r0, r1, pTread, sizeof(xthalfx8));

                        AE_LHX4X2_IP(z2, z3, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r2, r3, pTread, sizeof(xthalfx8));

                        AE_LHX4X2_IP(W0123, W456z, pW, sizeof(xthalfx8));
                        AE_LHX4X2_XP(Y02, Y03, pY0, stride);
                        MULACNVL_HX4X2(z0, z1, W0123, Y00, Y01);
                        MULACNVH_HX4X2(z2, z3, W0123, Y01, Y02);
                        MULACNVL_HX4X2(z0, z1, W456z, Y01, Y02);
                        MULACNVH_HX4X2(z2, z3, W456z, Y02, Y02);

                        MULACNVL_HX4X2(r0, r1, W0123, Y01, Y02);
                        MULACNVH_HX4X2(r2, r3, W0123, Y02, Y03);
                        MULACNVL_HX4X2(r0, r1, W456z, Y02, Y03);
                        MULACNVH_HX4X2(r2, r3, W456z, Y03, Y03);
                        Y00 = Y02; Y01 = Y03;

                        AE_LHX4X2_IP(W0123, W456z, pW, -(int)sizeof(xthalfx8));
                        AE_LHX4X2_XP(Y12, Y13, pY0, -stride+(int)sizeof(xthalfx8));
                        MULACNVL_HX4X2(z0, z1, W0123, Y10, Y11);
                        MULACNVH_HX4X2(z2, z3, W0123, Y11, Y12);
                        MULACNVL_HX4X2(z0, z1, W456z, Y11, Y12);
                        MULACNVH_HX4X2(z2, z3, W456z, Y12, Y12);

                        MULACNVL_HX4X2(r0, r1, W0123, Y11, Y12);
                        MULACNVH_HX4X2(r2, r3, W0123, Y12, Y13);
                        MULACNVL_HX4X2(r0, r1, W456z, Y12, Y13);
                        MULACNVH_HX4X2(r2, r3, W456z, Y13, Y13);
                        Y10 = Y12; Y11 = Y13;

                        AE_SHX4X2_IP(z0, z1, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r0, r1, pT, sizeof(xthalfx8));

                        AE_SHX4X2_IP(z2, z3, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r2, r3, pT, sizeof(xthalfx8));
                    }

                    //last 6
                    {
                        AE_LHX4X2_IP(z0, z1, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r0, r1, pTread, sizeof(xthalfx8));

                        AE_LHX4X2_IP(z2, z3, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r2, r3, pTread, sizeof(xthalfx8));

                        CONST_HX4X2(Y02, Y03, 0);
                        CONST_HX4X2(Y12, Y13, 0);

                        AE_LHX4X2_IP(W0123, W456z, pW, sizeof(xthalfx8));
                        MULACNVL_HX4X2(z0, z1, W0123, Y00, Y01);
                        MULACNVH_HX4X2(z2, z3, W0123, Y01, Y02);
                        MULACNVL_HX4X2(z0, z1, W456z, Y01, Y02);
                        MULACNVH_HX4X2(z2, z3, W456z, Y02, Y02);

                        MULACNVL_HX4X2(r0, r1, W0123, Y01, Y02);
                        MULACNVH_HX4X2(r2, r3, W0123, Y02, Y03);
                        MULACNVL_HX4X2(r0, r1, W456z, Y02, Y03);
                        MULACNVH_HX4X2(r2, r3, W456z, Y03, Y03);

                        AE_LHX4X2_IP(W0123, W456z, pW, -(int)sizeof(xthalfx8));
                        MULACNVL_HX4X2(z0, z1, W0123, Y10, Y11);
                        MULACNVH_HX4X2(z2, z3, W0123, Y11, Y12);
                        MULACNVL_HX4X2(z0, z1, W456z, Y11, Y12);
                        MULACNVH_HX4X2(z2, z3, W456z, Y12, Y12);

                        MULACNVL_HX4X2(r0, r1, W0123, Y11, Y12);
                        MULACNVH_HX4X2(r2, r3, W0123, Y12, Y13);
                        MULACNVL_HX4X2(r0, r1, W456z, Y12, Y13);
                        MULACNVH_HX4X2(r2, r3, W456z, Y13, Y13);


                        AE_SHX4X2_IP(z0, z1, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r0, r1, pT, sizeof(xthalfx8));

                        AE_SHX4X2_IP(z2, z3, pT, sizeof(xthalfx8));
                        AE_SHX4X2_IP(r2, r3, pT, sizeof(xthalfx8));
                    }

                    pW += 2;
                }

                // last 2 itteration
                {
                    //setup pointers and init variables
                    CONST_HX4X2(Y00, Y01, 0);
                    CONST_HX4X2(Y10, Y11, 0);
                    pY0 = pY;
                    pZ = (xthalfx8*)(z + (i + 0) * (N + Q - 1));
                    pTread = (xthalfx8*)((xthalfx4*)(pScr)+2 * M + (Q / 2));

                    __Pragma("loop_count min=1");
                    for (j = 0; j < Q; j += 8)
                    {
                        

                        AE_LHX4X2_IP(z0, z1, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r0, r1, pTread, sizeof(xthalfx8));

                        AE_LHX4X2_IP(z2, z3, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r2, r3, pTread, sizeof(xthalfx8));


                        AE_LHX4X2_IP(W0123, W456z, pW, sizeof(xthalfx8));
                        AE_LHX4X2_XP(Y02, Y03, pY0, stride);
                        MULACNVL_HX4X2(z0, z1, W0123, Y00, Y01);
                        MULACNVH_HX4X2(z2, z3, W0123, Y01, Y02);
                        MULACNVL_HX4X2(z0, z1, W456z, Y01, Y02);
                        MULACNVH_HX4X2(z2, z3, W456z, Y02, Y02);

                        MULACNVL_HX4X2(r0, r1, W0123, Y01, Y02);
                        MULACNVH_HX4X2(r2, r3, W0123, Y02, Y03);
                        MULACNVL_HX4X2(r0, r1, W456z, Y02, Y03);
                        MULACNVH_HX4X2(r2, r3, W456z, Y03, Y03);

                        Y00 = Y02; Y01 = Y03;

                        AE_LHX4X2_IP(W0123, W456z, pW, -(int)sizeof(xthalfx8));
                        AE_LHX4X2_XP(Y12, Y13, pY0, -stride+sizeof(xthalfx8));
                        MULACNVL_HX4X2(z0, z1, W0123, Y10, Y11);
                        MULACNVH_HX4X2(z2, z3, W0123, Y11, Y12);
                        MULACNVL_HX4X2(z0, z1, W456z, Y11, Y12);
                        MULACNVH_HX4X2(z2, z3, W456z, Y12, Y12);

                        MULACNVL_HX4X2(r0, r1, W0123, Y11, Y12);
                        MULACNVH_HX4X2(r2, r3, W0123, Y12, Y13);
                        MULACNVL_HX4X2(r0, r1, W456z, Y12, Y13);
                        MULACNVH_HX4X2(r2, r3, W456z, Y13, Y13);
                        Y10 = Y12; Y11 = Y13;


                        AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                        AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                        AE_DSEL16X4(t0, t2, t0, t2, sel);
                        AE_DSEL16X4(t1, t3, t1, t3, sel);
                        z0 = AE_MOVHALFX4_FROMF16X4(t0);
                        z1 = AE_MOVHALFX4_FROMF16X4(t2);
                        z2 = AE_MOVHALFX4_FROMF16X4(t1);
                        z3 = AE_MOVHALFX4_FROMF16X4(t3);

                        AE_DSEL16X4(k0, k1, AE_MOVF16X4_FROMHALFX4(r0), AE_MOVF16X4_FROMHALFX4(r1), sel);
                        AE_DSEL16X4(k2, k3, AE_MOVF16X4_FROMHALFX4(r2), AE_MOVF16X4_FROMHALFX4(r3), sel);
                        AE_DSEL16X4(k0, k2, k0, k2, sel);
                        AE_DSEL16X4(k1, k3, k1, k3, sel);
                        r0 = AE_MOVHALFX4_FROMF16X4(k0);
                        r1 = AE_MOVHALFX4_FROMF16X4(k2);
                        r2 = AE_MOVHALFX4_FROMF16X4(k1);
                        r3 = AE_MOVHALFX4_FROMF16X4(k3);

                        ADD_HX4X2(z0, z1, z0, z1, z2, z3);
                        ADD_HX4X2(r0, r1, r0, r1, r2, r3);
                        ADD_HX4X2(z0, r0, z0, r0, z1, r1);

                        AE_SAHX4X2_IP(z0, r0, alZ, pZ);

                    }
                    //last 6
                    {
                        CONST_HX4X2(Y02, Y03, 0);
                        CONST_HX4X2(Y12, Y13, 0);

                        AE_LHX4X2_IP(z0, z1, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r0, r1, pTread, sizeof(xthalfx8));

                        AE_LHX4X2_IP(z2, z3, pTread, sizeof(xthalfx8));
                        AE_LHX4X2_IP(r2, r3, pTread, sizeof(xthalfx8));

                        AE_LHX4X2_IP(W0123, W456z, pW, sizeof(xthalfx8));
                        MULACNVL_HX4X2(z0, z1, W0123, Y00, Y01);
                        MULACNVH_HX4X2(z2, z3, W0123, Y01, Y02);
                        MULACNVL_HX4X2(z0, z1, W456z, Y01, Y02);
                        MULACNVH_HX4X2(z2, z3, W456z, Y02, Y02);

                        MULACNVL_HX4X2(r0, r1, W0123, Y01, Y02);
                        MULACNVH_HX4X2(r2, r3, W0123, Y02, Y03);
                        MULACNVL_HX4X2(r0, r1, W456z, Y02, Y03);
                        MULACNVH_HX4X2(r2, r3, W456z, Y03, Y03);

                        AE_LHX4X2_IP(W0123, W456z, pW, -(int)sizeof(xthalfx8));
                        MULACNVL_HX4X2(z0, z1, W0123, Y10, Y11);
                        MULACNVH_HX4X2(z2, z3, W0123, Y11, Y12);
                        MULACNVL_HX4X2(z0, z1, W456z, Y11, Y12);
                        MULACNVH_HX4X2(z2, z3, W456z, Y12, Y12);

                        MULACNVL_HX4X2(r0, r1, W0123, Y11, Y12);
                        MULACNVH_HX4X2(r2, r3, W0123, Y12, Y13);
                        MULACNVL_HX4X2(r0, r1, W456z, Y12, Y13);
                        MULACNVH_HX4X2(r2, r3, W456z, Y13, Y13);


                        AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                        AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                        AE_DSEL16X4(t0, t2, t0, t2, sel);
                        AE_DSEL16X4(t1, t3, t1, t3, sel);
                        z0 = AE_MOVHALFX4_FROMF16X4(t0);
                        z1 = AE_MOVHALFX4_FROMF16X4(t2);
                        z2 = AE_MOVHALFX4_FROMF16X4(t1);
                        z3 = AE_MOVHALFX4_FROMF16X4(t3);

                        AE_DSEL16X4(k0, k1, AE_MOVF16X4_FROMHALFX4(r0), AE_MOVF16X4_FROMHALFX4(r1), sel);
                        AE_DSEL16X4(k2, k3, AE_MOVF16X4_FROMHALFX4(r2), AE_MOVF16X4_FROMHALFX4(r3), sel);
                        AE_DSEL16X4(k0, k2, k0, k2, sel);
                        AE_DSEL16X4(k1, k3, k1, k3, sel);
                        r0 = AE_MOVHALFX4_FROMF16X4(k0);
                        r1 = AE_MOVHALFX4_FROMF16X4(k2);
                        r2 = AE_MOVHALFX4_FROMF16X4(k1);
                        r3 = AE_MOVHALFX4_FROMF16X4(k3);

                        ADD_HX4X2(z0, z1, z0, z1, z2, z3);
                        ADD_HX4X2(r0, r1, r0, r1, r2, r3);
                        ADD_HX4X2(z0, r0, z0, r0, z1, r1);

                        AE_SAVHX4X2_XP(z0, r0, alZ, pZ, 6 * sizeof(xthalf));
                    }
                    AE_SA128POS_FP(alZ, pZ);
                }

            }
        }

    }

#   undef M
#   undef N
}

size_t conv2d_11x7_fp16_getScratchSize(int P, int Q)
{
    int M = 11;
    int N = 7;
    return (M * 2 + Q / 4 + (Q+N-1)*4 + 8) * sizeof(xthalfx4);

} // MxN=11x7

#else 

#define MAX(a,b)   ( (a)>(b) ? (a) : (b)  )
#define MIN(a,b)   ( (a)<(b) ? (a) : (b)  )
#define ABS(a)     ( (a)>0   ? (a) : (-a) )

#include "float16.h"

/* simplest reference code */
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
            m0 = XT_MAX(i - P + 1, 0);
            m1 = XT_MIN(i + 1, M);
            n0 = XT_MAX(j - Q + 1, 0);
            n1 = XT_MIN(j + 1, N);
            s = 0;
            for (n = n0; n < n1; n++)
                for (m = m0; m < m1; m++)
                {
                    s = fma_f16(x[m * N + n], y[(i - m) * Q + (j - n)], s);
                }
            z[i * (N + Q - 1) + j] = s;
        }
}/* conv2df() */

size_t conv2d_11x7_fp16_getScratchSize(int M, int N, int P, int Q)
{
    return 0;
}

void conv2d_11x7_fp16(void* pScr, float16_t* z, const float16_t* x, const float16_t* y, int P, int Q)
{
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT_ALIGN(z, 16);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0 && P % 8 == 0);
    NASSERT(Q >= 0 && Q % 8 == 0);
    conv2d_fp16(pScr, z, x, y, 11, 7, P, Q);
}

#endif
