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
size_t conv2d_gen_5x5_fp16_getScratchSize (int P, int Q) 
{ 
    (void)P,(void)Q;
    return 0;
}
DISCARD_FUN(void,conv2d_gen_5x5_fp16,(void* pScr, float16_t *z, const float16_t * x, const float16_t * y, int P, int Q))
#elif HAVE_HPFPU
void conv2d_gen_5x5_fp16(void* pScr, float16_t* z, const float16_t* x, const float16_t* y, int P, int Q)
{
#define M 5
#define N 5
    const xthalfx8* restrict pY;
    const xthalfx8* restrict pY0;
    const xthalfx8* restrict pY1;
    const xthalfx8* restrict pY2;
    const xthalfx8* restrict pY3;
    const xthalfx8* restrict pY4;
    const xthalfx8* restrict pT_read;
    xthalfx8* restrict pW;
    const ae_int16* restrict pX;
    xthalfx8* restrict pT;
    xthalfx8* restrict pZ0;

    xthalfx4 z0,z1,z2,z3,z4,z5,z6,z7;

    xthalfx4 r0, r1, r2, r3;
    ae_int16x4 t0, t1, t2, t3;
    ae_int16x4 k0, k1, k2, k3;

    xthalfx4 S0, S1;
    xthalfx4 Y00, Y01, Y02;
    xthalfx4 Y10, Y11, Y12;
    xthalfx4 Y20, Y21, Y22;
    xthalfx4 Y30, Y31, Y32;
    xthalfx4 Y40, Y41, Y42;

    xthalfx4 W0123_0, W4_0;
    

    ae_valignx2 al0, al1, al2, al3, al4, alZ;

    ae_int16x4 sel = AE_MOVINT16X4_FROMINT64(0x0706050403020100); // SEL7531 + SEL6420

    int i, j, mmin, mmax, xstart, ystart, m;

    NASSERT(x && y && z && pScr);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0);
    NASSERT(Q >= 0);
    if (P <= 0 || Q <= 0) return;

    alZ = AE_ZALIGN128();

    pX = (ae_int16*)x;
    {
        xthalfx8* pWw = (xthalfx8*)(pScr) + 4;
        for (i = 0; i < 5; ++i)
        {
            xthalfx4 w0123, W4, W3, W2, W1, W0;
            ae_int16x4 v_i;
        
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
        
            w0123 = AE_SELH_6543(CONST_HX4(0), W0);
            w0123 = AE_SELH_6543(w0123, W1);
            w0123 = AE_SELH_6543(w0123, W2);
            w0123 = AE_SELH_6543(w0123, W3);
        
            AE_SHX4X2_IP(w0123, W4, pWw, -(int)sizeof(xthalfx8));
        }

    }

    for (i = 0; i < P + M - 1; i += 1)
    {
        mmin = i < M - 1 ? M - 1 - i : 0;
        mmax = i >= P ? M - 1 + P - i : M;
        xstart = i < M - 1 ? (mmin) : 0;
        ystart = i < M - 1 ? 0 : i - (M - 1);

        pW = (xthalfx8*)pScr + xstart;
        pY0 = (xthalfx8*)(y + ystart * Q);
        pZ0 = (xthalfx8*)(z + i * (Q + N - 1));

        if (mmax - mmin == 5)
        {
            Y00 = CONST_HX4(0);
            Y10 = CONST_HX4(0);
            Y20 = CONST_HX4(0);
            Y30 = CONST_HX4(0);
            Y40 = CONST_HX4(0);

            pY1 = (xthalfx8*)((xthalf*)pY0 + Q);
            pY2 = (xthalfx8*)((xthalf*)pY0 + 2*Q);
            pY3 = (xthalfx8*)((xthalf*)pY0 + 3*Q);
            pY4 = (xthalfx8*)((xthalf*)pY0 + 4*Q);

            al0 = AE_LA128_PP(pY0);
            al1 = AE_LA128_PP(pY1);
            al2 = AE_LA128_PP(pY2);
            al3 = AE_LA128_PP(pY3);
            al4 = AE_LA128_PP(pY4);

            for (j = 0; j < (Q&~7); j += 8)
            {

                AE_LAHX4X2_IP(Y01, Y02, al0, pY0);
                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);
                MULCNVH_HX4X2(r0, r1, W0123_0, Y01, Y02);
                MULCNVL_HX4X2(r2, r3, W0123_0, Y01, Y02);
                MUL_HX4X2(S0, S1, W4_0, W4_0, Y01, Y02);
                Y00 = Y02;

                AE_LAHX4X2_IP(Y11, Y12, al1, pY1);
                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123_0, Y10, Y11);
                MULACNVL_HX4X2(z2, z3, W0123_0, Y10, Y11);
                MULACNVH_HX4X2(r0, r1, W0123_0, Y11, Y12);
                MULACNVL_HX4X2(r2, r3, W0123_0, Y11, Y12);
                MADD_HX4X2(S0, S1, W4_0, W4_0, Y11, Y12);
                Y10 = Y12;

                AE_LAHX4X2_IP(Y21, Y22, al2, pY2);
                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123_0, Y20, Y21);
                MULACNVL_HX4X2(z2, z3, W0123_0, Y20, Y21);
                MULACNVH_HX4X2(r0, r1, W0123_0, Y21, Y22);
                MULACNVL_HX4X2(r2, r3, W0123_0, Y21, Y22);
                MADD_HX4X2(S0, S1, W4_0, W4_0, Y21, Y22);
                Y20 = Y22;

                AE_LAHX4X2_IP(Y31, Y32, al3, pY3);
                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123_0, Y30, Y31);
                MULACNVL_HX4X2(z2, z3, W0123_0, Y30, Y31);
                MULACNVH_HX4X2(r0, r1, W0123_0, Y31, Y32);
                MULACNVL_HX4X2(r2, r3, W0123_0, Y31, Y32);
                MADD_HX4X2(S0, S1, W4_0, W4_0, Y31, Y32);
                Y30 = Y32;

                AE_LAHX4X2_IP(Y41, Y42, al4, pY4);
                AE_LHX4X2_IP(W0123_0, W4_0, pW, -4 * (int)sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123_0, Y40, Y41);
                MULACNVL_HX4X2(z2, z3, W0123_0, Y40, Y41);
                MULACNVH_HX4X2(r0, r1, W0123_0, Y41, Y42);
                MULACNVL_HX4X2(r2, r3, W0123_0, Y41, Y42);
                MADD_HX4X2(S0, S1, W4_0, W4_0, Y41, Y42);
                Y40 = Y42;

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
                ADD_HX4X2(S0, S1, z0, r0, S0, S1);

                AE_SAHX4X2_IP(S0, S1, alZ, pZ0);
            }
            // tail
            if (Q&7)
            {
                AE_LAVHX4X2_XP(Y01, Y02, al0, pY0, (Q & 7) * sizeof(xthalf));
                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);
                MULCNVH_HX4X2(r0, r1, W0123_0, Y01, Y02);
                MULCNVL_HX4X2(r2, r3, W0123_0, Y01, Y02);
                MUL_HX4X2(S0, S1, W4_0, W4_0, Y01, Y02);

                AE_LAVHX4X2_XP(Y11, Y12, al1, pY1, (Q & 7) * sizeof(xthalf));
                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123_0, Y10, Y11);
                MULACNVL_HX4X2(z2, z3, W0123_0, Y10, Y11);
                MULACNVH_HX4X2(r0, r1, W0123_0, Y11, Y12);
                MULACNVL_HX4X2(r2, r3, W0123_0, Y11, Y12);
                MADD_HX4X2(S0, S1, W4_0, W4_0, Y11, Y12);
 

                AE_LAVHX4X2_XP(Y21, Y22, al2, pY2, (Q & 7) * sizeof(xthalf));
                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123_0, Y20, Y21);
                MULACNVL_HX4X2(z2, z3, W0123_0, Y20, Y21);
                MULACNVH_HX4X2(r0, r1, W0123_0, Y21, Y22);
                MULACNVL_HX4X2(r2, r3, W0123_0, Y21, Y22);
                MADD_HX4X2(S0, S1, W4_0, W4_0, Y21, Y22);

                AE_LAVHX4X2_XP(Y31, Y32, al3, pY3, (Q & 7) * sizeof(xthalf));
                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123_0, Y30, Y31);
                MULACNVL_HX4X2(z2, z3, W0123_0, Y30, Y31);
                MULACNVH_HX4X2(r0, r1, W0123_0, Y31, Y32);
                MULACNVL_HX4X2(r2, r3, W0123_0, Y31, Y32);
                MADD_HX4X2(S0, S1, W4_0, W4_0, Y31, Y32);

                AE_LAVHX4X2_XP(Y41, Y42, al4, pY4, (Q & 7) * sizeof(xthalf));
                AE_LHX4X2_IP(W0123_0, W4_0, pW, -4 * (int)sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123_0, Y40, Y41);
                MULACNVL_HX4X2(z2, z3, W0123_0, Y40, Y41);
                MULACNVH_HX4X2(r0, r1, W0123_0, Y41, Y42);
                MULACNVL_HX4X2(r2, r3, W0123_0, Y41, Y42);
                MADD_HX4X2(S0, S1, W4_0, W4_0, Y41, Y42);

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
                ADD_HX4X2(S0, S1, z0, r0, S0, S1);

                AE_SAVHX4X2_XP(S0, S1, alZ, pZ0, (Q & 7) * sizeof(xthalf));
                
                for (j=0; j < (Q&7); ++j)
                {
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

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y30), AE_MOVF16X4_FROMHALFX4(Y31), AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    Y30 = AE_MOVHALFX4_FROMF16X4(t0); Y31 = AE_MOVHALFX4_FROMF16X4(t1);
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y31), AE_MOVF16X4_FROMHALFX4(Y32), AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    Y31 = AE_MOVHALFX4_FROMF16X4(t0); Y32 = AE_MOVHALFX4_FROMF16X4(t1);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y40), AE_MOVF16X4_FROMHALFX4(Y41), AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    Y40 = AE_MOVHALFX4_FROMF16X4(t0); Y41 = AE_MOVHALFX4_FROMF16X4(t1);
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y41), AE_MOVF16X4_FROMHALFX4(Y42), AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    Y41 = AE_MOVHALFX4_FROMF16X4(t0); Y42 = AE_MOVHALFX4_FROMF16X4(t1);
                }
            }
            //last N-1 elements
            {
                Y01 = CONST_HX4(0); Y02 = CONST_HX4(0);
                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);
                MULCNVH_HX4X2(r0, r1, W0123_0, Y01, Y02);
                MULCNVL_HX4X2(r2, r3, W0123_0, Y01, Y02);
                MUL_HX4X2(S0, S1, W4_0, W4_0, Y01, Y02);

                Y11 = CONST_HX4(0); Y12 = CONST_HX4(0);
                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123_0, Y10, Y11);
                MULACNVL_HX4X2(z2, z3, W0123_0, Y10, Y11);
                MULACNVH_HX4X2(r0, r1, W0123_0, Y11, Y12);
                MULACNVL_HX4X2(r2, r3, W0123_0, Y11, Y12);
                MADD_HX4X2(S0, S1, W4_0, W4_0, Y11, Y12);

                Y21 = CONST_HX4(0); Y22 = CONST_HX4(0);
                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123_0, Y20, Y21);
                MULACNVL_HX4X2(z2, z3, W0123_0, Y20, Y21);
                MULACNVH_HX4X2(r0, r1, W0123_0, Y21, Y22);
                MULACNVL_HX4X2(r2, r3, W0123_0, Y21, Y22);
                MADD_HX4X2(S0, S1, W4_0, W4_0, Y21, Y22);

                Y31 = CONST_HX4(0); Y32 = CONST_HX4(0);
                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123_0, Y30, Y31);
                MULACNVL_HX4X2(z2, z3, W0123_0, Y30, Y31);
                MULACNVH_HX4X2(r0, r1, W0123_0, Y31, Y32);
                MULACNVL_HX4X2(r2, r3, W0123_0, Y31, Y32);
                MADD_HX4X2(S0, S1, W4_0, W4_0, Y31, Y32);

                Y41 = CONST_HX4(0); Y42 = CONST_HX4(0);
                AE_LHX4X2_IP(W0123_0, W4_0, pW, -4 * (int)sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123_0, Y40, Y41);
                MULACNVL_HX4X2(z2, z3, W0123_0, Y40, Y41);
                MULACNVH_HX4X2(r0, r1, W0123_0, Y41, Y42);
                MULACNVL_HX4X2(r2, r3, W0123_0, Y41, Y42);
                MADD_HX4X2(S0, S1, W4_0, W4_0, Y41, Y42);

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
                ADD_HX4X2(S0, S1, z0, r0, S0, S1);

                AE_SAVHX4X2_XP(z0, z4, alZ, pZ0, 4 * sizeof(xthalf));
            }
            AE_SA128POS_FP(alZ, pZ0);
        }
        else if (mmax - mmin != 1)
        {
            pY = (xthalfx8*)((xthalf*)pY0 + Q);
            // first iteration
            {
                pT = (xthalfx8*)pScr + (M + 1);

                Y00 = CONST_HX4(0);
                al0 = AE_LA128_PP(pY0);

                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                for (j = 0; j < (Q & ~7); j += 8)
                {
                    AE_LAHX4X2_IP(Y01, Y02, al0, pY0);

                    MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                    MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);

                    MULCNVH_HX4X2(z4, z5, W0123_0, Y01, Y02);
                    MULCNVL_HX4X2(z6, z7, W0123_0, Y01, Y02);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z0 = AE_MOVHALFX4_FROMF16X4(t0);
                    z1 = AE_MOVHALFX4_FROMF16X4(t2);
                    z2 = AE_MOVHALFX4_FROMF16X4(t1);
                    z3 = AE_MOVHALFX4_FROMF16X4(t3);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z4), AE_MOVF16X4_FROMHALFX4(z5), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z6), AE_MOVF16X4_FROMHALFX4(z7), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z4 = AE_MOVHALFX4_FROMF16X4(t0);
                    z5 = AE_MOVHALFX4_FROMF16X4(t2);
                    z6 = AE_MOVHALFX4_FROMF16X4(t1);
                    z7 = AE_MOVHALFX4_FROMF16X4(t3);

                    Y00 = Y02;

                    ADD_HX4X2(z0, z4, z0, z4, z1, z5);
                    ADD_HX4X2(z2, z6, z2, z6, z3, z7);
                    ADD_HX4X2(z0, z4, z0, z4, z2, z6);

                    MADDQ_H(z0, z4, Y01, Y02, W4_0);

                    AE_SHX4X2_IP(z0, z4, pT, sizeof(xthalfx8));
                }
                // tail
                if (j < Q)
                {
                    AE_LAHX4X2_IP(Y01, Y02, al0, pY0);

                    MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                    MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);

                    MULCNVH_HX4X2(z4, z5, W0123_0, Y01, Y02);
                    MULCNVL_HX4X2(z6, z7, W0123_0, Y01, Y02);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z0 = AE_MOVHALFX4_FROMF16X4(t0);
                    z1 = AE_MOVHALFX4_FROMF16X4(t2);
                    z2 = AE_MOVHALFX4_FROMF16X4(t1);
                    z3 = AE_MOVHALFX4_FROMF16X4(t3);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z4), AE_MOVF16X4_FROMHALFX4(z5), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z6), AE_MOVF16X4_FROMHALFX4(z7), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z4 = AE_MOVHALFX4_FROMF16X4(t0);
                    z5 = AE_MOVHALFX4_FROMF16X4(t2);
                    z6 = AE_MOVHALFX4_FROMF16X4(t1);
                    z7 = AE_MOVHALFX4_FROMF16X4(t3);

                    ADD_HX4X2(z0, z4, z0, z4, z1, z5);
                    ADD_HX4X2(z2, z6, z2, z6, z3, z7);
                    ADD_HX4X2(z0, z4, z0, z4, z2, z6);

                    MADDQ_H(z0, z4, Y01, Y02, W4_0);

                    AE_SHX4X2_IP(z0, z4, pT, sizeof(xthalfx8));

                    for (; j < Q; ++j)
                    {
                        AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                        Y00 = AE_MOVHALFX4_FROMF16X4(t0); Y01 = AE_MOVHALFX4_FROMF16X4(t1);
                        AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                        Y01 = AE_MOVHALFX4_FROMF16X4(t0); Y02 = AE_MOVHALFX4_FROMF16X4(t1);
                    }
                }
                //last N-1 elements
                {
                    Y01 = CONST_HX4(0); Y02 = CONST_HX4(0);

                    MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                    MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);

                    MULCNVH_HX4X2(z4, z5, W0123_0, Y01, Y02);
                    MULCNVL_HX4X2(z6, z7, W0123_0, Y01, Y02);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z0 = AE_MOVHALFX4_FROMF16X4(t0);
                    z1 = AE_MOVHALFX4_FROMF16X4(t2);
                    z2 = AE_MOVHALFX4_FROMF16X4(t1);
                    z3 = AE_MOVHALFX4_FROMF16X4(t3);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z4), AE_MOVF16X4_FROMHALFX4(z5), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z6), AE_MOVF16X4_FROMHALFX4(z7), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z4 = AE_MOVHALFX4_FROMF16X4(t0);
                    z5 = AE_MOVHALFX4_FROMF16X4(t2);
                    z6 = AE_MOVHALFX4_FROMF16X4(t1);
                    z7 = AE_MOVHALFX4_FROMF16X4(t3);

                    ADD_HX4X2(z0, z4, z0, z4, z1, z5);
                    ADD_HX4X2(z2, z6, z2, z6, z3, z7);
                    ADD_HX4X2(z0, z4, z0, z4, z2, z6);

                    MADDQ_H(z0, z4, Y01, Y02, W4_0);

                    AE_SHX4X2_IP(z0, z4, pT, sizeof(xthalfx8));
                }
            }
            pY0 = pY;
            pY = (xthalfx8*)((xthalf*)pY + Q);

            __Pragma("loop_count min=0, max=2");
            for (m = mmin + 1; m < mmax - 1; ++m)
            {
                pT_read = pT = (xthalfx8*)pScr + (M + 1);

                Y00 = CONST_HX4(0);
                al0 = AE_LA128_PP(pY0);

                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                for (j = 0; j < (Q & ~7); j += 8)
                {
                    AE_LAHX4X2_IP(Y01, Y02, al0, pY0);

                    MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                    MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);

                    MULCNVH_HX4X2(z4, z5, W0123_0, Y01, Y02);
                    MULCNVL_HX4X2(z6, z7, W0123_0, Y01, Y02);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z0 = AE_MOVHALFX4_FROMF16X4(t0);
                    z1 = AE_MOVHALFX4_FROMF16X4(t2);
                    z2 = AE_MOVHALFX4_FROMF16X4(t1);
                    z3 = AE_MOVHALFX4_FROMF16X4(t3);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z4), AE_MOVF16X4_FROMHALFX4(z5), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z6), AE_MOVF16X4_FROMHALFX4(z7), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z4 = AE_MOVHALFX4_FROMF16X4(t0);
                    z5 = AE_MOVHALFX4_FROMF16X4(t2);
                    z6 = AE_MOVHALFX4_FROMF16X4(t1);
                    z7 = AE_MOVHALFX4_FROMF16X4(t3);

                    Y00 = Y02;

                    ADD_HX4X2(z0, z4, z0, z4, z1, z5);
                    ADD_HX4X2(z2, z6, z2, z6, z3, z7);
                    ADD_HX4X2(z0, z4, z0, z4, z2, z6);

                    MADDQ_H(z0, z4, Y01, Y02, W4_0);

                    AE_LHX4X2_IP(S0, S1, pT_read, sizeof(xthalfx8));
                    ADD_HX4X2(z0, z4, z0, z4, S0, S1);

                    AE_SHX4X2_IP(z0, z4, pT, sizeof(xthalfx8));
                }
                // tail
                if (Q&7)
                {
                    AE_LAHX4X2_IP(Y01, Y02, al0, pY0);

                    MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                    MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);

                    MULCNVH_HX4X2(z4, z5, W0123_0, Y01, Y02);
                    MULCNVL_HX4X2(z6, z7, W0123_0, Y01, Y02);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z0 = AE_MOVHALFX4_FROMF16X4(t0);
                    z1 = AE_MOVHALFX4_FROMF16X4(t2);
                    z2 = AE_MOVHALFX4_FROMF16X4(t1);
                    z3 = AE_MOVHALFX4_FROMF16X4(t3);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z4), AE_MOVF16X4_FROMHALFX4(z5), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z6), AE_MOVF16X4_FROMHALFX4(z7), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z4 = AE_MOVHALFX4_FROMF16X4(t0);
                    z5 = AE_MOVHALFX4_FROMF16X4(t2);
                    z6 = AE_MOVHALFX4_FROMF16X4(t1);
                    z7 = AE_MOVHALFX4_FROMF16X4(t3);

                    ADD_HX4X2(z0, z4, z0, z4, z1, z5);
                    ADD_HX4X2(z2, z6, z2, z6, z3, z7);
                    ADD_HX4X2(z0, z4, z0, z4, z2, z6);

                    MADDQ_H(z0, z4, Y01, Y02, W4_0);

                    AE_LHX4X2_IP(S0, S1, pT_read, sizeof(xthalfx8));
                    ADD_HX4X2(z0, z4, z0, z4, S0, S1);

                    AE_SHX4X2_IP(z0, z4, pT, sizeof(xthalfx8));

                    for (j=0; j < (Q&7); ++j)
                    {
                        AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                        Y00 = AE_MOVHALFX4_FROMF16X4(t0); Y01 = AE_MOVHALFX4_FROMF16X4(t1);
                        AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                        Y01 = AE_MOVHALFX4_FROMF16X4(t0); Y02 = AE_MOVHALFX4_FROMF16X4(t1);
                    }
                }
                //last N-1 elements
                {
                    Y01 = CONST_HX4(0); Y02 = CONST_HX4(0);

                    MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                    MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);

                    MULCNVH_HX4X2(z4, z5, W0123_0, Y01, Y02);
                    MULCNVL_HX4X2(z6, z7, W0123_0, Y01, Y02);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z0 = AE_MOVHALFX4_FROMF16X4(t0);
                    z1 = AE_MOVHALFX4_FROMF16X4(t2);
                    z2 = AE_MOVHALFX4_FROMF16X4(t1);
                    z3 = AE_MOVHALFX4_FROMF16X4(t3);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z4), AE_MOVF16X4_FROMHALFX4(z5), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z6), AE_MOVF16X4_FROMHALFX4(z7), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z4 = AE_MOVHALFX4_FROMF16X4(t0);
                    z5 = AE_MOVHALFX4_FROMF16X4(t2);
                    z6 = AE_MOVHALFX4_FROMF16X4(t1);
                    z7 = AE_MOVHALFX4_FROMF16X4(t3);

                    ADD_HX4X2(z0, z4, z0, z4, z1, z5);
                    ADD_HX4X2(z2, z6, z2, z6, z3, z7);
                    ADD_HX4X2(z0, z4, z0, z4, z2, z6);

                    MADDQ_H(z0, z4, Y01, Y02, W4_0);

                    AE_LHX4X2_IP(S0, S1, pT_read, sizeof(xthalfx8));
                    ADD_HX4X2(z0, z4, z0, z4, S0, S1);

                    AE_SHX4X2_IP(z0, z4, pT, sizeof(xthalfx8));
                }

                pY0 = pY;
                pY = (xthalfx8*)((ae_int16*)pY0 + Q);
            }

            //last iteration
            {
                pT_read = pT = (xthalfx8*)pScr + (M+1);

                Y00 = CONST_HX4(0);
                al0 = AE_LA128_PP(pY0);

                AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
                for (j = 0; j < (Q & ~7); j += 8)
                {
                    AE_LAHX4X2_IP(Y01, Y02, al0, pY0);

                    MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                    MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);

                    MULCNVH_HX4X2(z4, z5, W0123_0, Y01, Y02);
                    MULCNVL_HX4X2(z6, z7, W0123_0, Y01, Y02);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z0 = AE_MOVHALFX4_FROMF16X4(t0);
                    z1 = AE_MOVHALFX4_FROMF16X4(t2);
                    z2 = AE_MOVHALFX4_FROMF16X4(t1);
                    z3 = AE_MOVHALFX4_FROMF16X4(t3);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z4), AE_MOVF16X4_FROMHALFX4(z5), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z6), AE_MOVF16X4_FROMHALFX4(z7), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z4 = AE_MOVHALFX4_FROMF16X4(t0);
                    z5 = AE_MOVHALFX4_FROMF16X4(t2);
                    z6 = AE_MOVHALFX4_FROMF16X4(t1);
                    z7 = AE_MOVHALFX4_FROMF16X4(t3);

                    Y00 = Y02;

                    ADD_HX4X2(z0, z4, z0, z4, z1, z5);
                    ADD_HX4X2(z2, z6, z2, z6, z3, z7);
                    ADD_HX4X2(z0, z4, z0, z4, z2, z6);

                    MADDQ_H(z0, z4, Y01, Y02, W4_0);

                    AE_LHX4X2_IP(S0, S1, pT_read, sizeof(xthalfx8));
                    ADD_HX4X2(z0, z4, z0, z4, S0, S1);

                    AE_SAHX4X2_IP(z0, z4, alZ, pZ0);
                }
                // tail
                if (Q&7)
                {
                    AE_LAVHX4X2_XP(Y01, Y02, al0, pY0, (Q & 7) * sizeof(xthalf));

                    MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                    MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);

                    MULCNVH_HX4X2(z4, z5, W0123_0, Y01, Y02);
                    MULCNVL_HX4X2(z6, z7, W0123_0, Y01, Y02);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z0 = AE_MOVHALFX4_FROMF16X4(t0);
                    z1 = AE_MOVHALFX4_FROMF16X4(t2);
                    z2 = AE_MOVHALFX4_FROMF16X4(t1);
                    z3 = AE_MOVHALFX4_FROMF16X4(t3);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z4), AE_MOVF16X4_FROMHALFX4(z5), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z6), AE_MOVF16X4_FROMHALFX4(z7), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z4 = AE_MOVHALFX4_FROMF16X4(t0);
                    z5 = AE_MOVHALFX4_FROMF16X4(t2);
                    z6 = AE_MOVHALFX4_FROMF16X4(t1);
                    z7 = AE_MOVHALFX4_FROMF16X4(t3);

                    ADD_HX4X2(z0, z4, z0, z4, z1, z5);
                    ADD_HX4X2(z2, z6, z2, z6, z3, z7);
                    ADD_HX4X2(z0, z4, z0, z4, z2, z6);

                    MADDQ_H(z0, z4, Y01, Y02, W4_0);

                    AE_LHX4X2_IP(S0, S1, pT_read, sizeof(xthalfx8));
                    ADD_HX4X2(z0, z4, z0, z4, S0, S1);

                    AE_SAVHX4X2_XP(z0, z4, alZ, pZ0, (Q & 7) * sizeof(xthalf));

                    for (j=0; j < (Q&7); ++j)
                    {
                        AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                        Y00 = AE_MOVHALFX4_FROMF16X4(t0); Y01 = AE_MOVHALFX4_FROMF16X4(t1);
                        AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                        Y01 = AE_MOVHALFX4_FROMF16X4(t0); Y02 = AE_MOVHALFX4_FROMF16X4(t1);
                    }
                }
                //last N-1 elements
                {
                    Y01 = CONST_HX4(0); Y02 = CONST_HX4(0);

                    MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                    MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);

                    MULCNVH_HX4X2(z4, z5, W0123_0, Y01, Y02);
                    MULCNVL_HX4X2(z6, z7, W0123_0, Y01, Y02);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z0 = AE_MOVHALFX4_FROMF16X4(t0);
                    z1 = AE_MOVHALFX4_FROMF16X4(t2);
                    z2 = AE_MOVHALFX4_FROMF16X4(t1);
                    z3 = AE_MOVHALFX4_FROMF16X4(t3);

                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z4), AE_MOVF16X4_FROMHALFX4(z5), sel);
                    AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z6), AE_MOVF16X4_FROMHALFX4(z7), sel);
                    AE_DSEL16X4(t0, t2, t0, t2, sel);
                    AE_DSEL16X4(t1, t3, t1, t3, sel);
                    z4 = AE_MOVHALFX4_FROMF16X4(t0);
                    z5 = AE_MOVHALFX4_FROMF16X4(t2);
                    z6 = AE_MOVHALFX4_FROMF16X4(t1);
                    z7 = AE_MOVHALFX4_FROMF16X4(t3);

                    ADD_HX4X2(z0, z4, z0, z4, z1, z5);
                    ADD_HX4X2(z2, z6, z2, z6, z3, z7);
                    ADD_HX4X2(z0, z4, z0, z4, z2, z6);

                    MADDQ_H(z0, z4, Y01, Y02, W4_0);

                    AE_LHX4X2_IP(S0, S1, pT_read, sizeof(xthalfx8));
                    ADD_HX4X2(z0, z4, z0, z4, S0, S1);

                    AE_SAVHX4X2_XP(z0, z4, alZ, pZ0, 4 * sizeof(xthalf));
                }
                AE_SA128POS_FP(alZ, pZ0);
            }
        }
        else
        {
            Y00 = CONST_HX4(0);
            al0 = AE_LA128_PP(pY0);

            AE_LHX4X2_IP(W0123_0, W4_0, pW, sizeof(xthalfx8));
            for (j = 0; j < (Q & ~7); j += 8)
            {
                AE_LAHX4X2_IP(Y01, Y02, al0, pY0);

                MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);

                MULCNVH_HX4X2(z4, z5, W0123_0, Y01, Y02);
                MULCNVL_HX4X2(z6, z7, W0123_0, Y01, Y02);

                AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                AE_DSEL16X4(t0, t2, t0, t2, sel);
                AE_DSEL16X4(t1, t3, t1, t3, sel);
                z0 = AE_MOVHALFX4_FROMF16X4(t0);
                z1 = AE_MOVHALFX4_FROMF16X4(t2);
                z2 = AE_MOVHALFX4_FROMF16X4(t1);
                z3 = AE_MOVHALFX4_FROMF16X4(t3);

                AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z4), AE_MOVF16X4_FROMHALFX4(z5), sel);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z6), AE_MOVF16X4_FROMHALFX4(z7), sel);
                AE_DSEL16X4(t0, t2, t0, t2, sel);
                AE_DSEL16X4(t1, t3, t1, t3, sel);
                z4 = AE_MOVHALFX4_FROMF16X4(t0);
                z5 = AE_MOVHALFX4_FROMF16X4(t2);
                z6 = AE_MOVHALFX4_FROMF16X4(t1);
                z7 = AE_MOVHALFX4_FROMF16X4(t3);

                Y00 = Y02;

                ADD_HX4X2(z0, z4, z0, z4, z1, z5);
                ADD_HX4X2(z2, z6, z2, z6, z3, z7);
                ADD_HX4X2(z0, z4, z0, z4, z2, z6);

                MADDQ_H(z0, z4, Y01, Y02, W4_0);

                AE_SAHX4X2_IP(z0, z4, alZ, pZ0);
            }
            // tail
            if (Q&7)
            {
                AE_LAVHX4X2_XP(Y01, Y02, al0, pY0, (Q&7)*sizeof(xthalf));

                MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);

                MULCNVH_HX4X2(z4, z5, W0123_0, Y01, Y02);
                MULCNVL_HX4X2(z6, z7, W0123_0, Y01, Y02);

                AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                AE_DSEL16X4(t0, t2, t0, t2, sel);
                AE_DSEL16X4(t1, t3, t1, t3, sel);
                z0 = AE_MOVHALFX4_FROMF16X4(t0);
                z1 = AE_MOVHALFX4_FROMF16X4(t2);
                z2 = AE_MOVHALFX4_FROMF16X4(t1);
                z3 = AE_MOVHALFX4_FROMF16X4(t3);

                AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z4), AE_MOVF16X4_FROMHALFX4(z5), sel);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z6), AE_MOVF16X4_FROMHALFX4(z7), sel);
                AE_DSEL16X4(t0, t2, t0, t2, sel);
                AE_DSEL16X4(t1, t3, t1, t3, sel);
                z4 = AE_MOVHALFX4_FROMF16X4(t0);
                z5 = AE_MOVHALFX4_FROMF16X4(t2);
                z6 = AE_MOVHALFX4_FROMF16X4(t1);
                z7 = AE_MOVHALFX4_FROMF16X4(t3);

                ADD_HX4X2(z0, z4, z0, z4, z1, z5);
                ADD_HX4X2(z2, z6, z2, z6, z3, z7);
                ADD_HX4X2(z0, z4, z0, z4, z2, z6);

                MADDQ_H(z0, z4, Y01, Y02, W4_0);

                AE_SAVHX4X2_XP(z0, z4, alZ, pZ0, (Q & 7) * sizeof(xthalf));

                for (j=0; j < (Q&7); ++j)
                {
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y00), AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVINT16X4_FROMINT64(0x0602050104000300));
                    Y00 = AE_MOVHALFX4_FROMF16X4(t0); Y01 = AE_MOVHALFX4_FROMF16X4(t1);
                    AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(Y01), AE_MOVF16X4_FROMHALFX4(Y02), AE_MOVINT16X4_FROMINT64(0x0702060105000300));
                    Y01 = AE_MOVHALFX4_FROMF16X4(t0); Y02 = AE_MOVHALFX4_FROMF16X4(t1);
                }
            }
            //last N-1 elements
            {
                Y01 = CONST_HX4(0); Y02 = CONST_HX4(0);

                MULCNVH_HX4X2(z0, z1, W0123_0, Y00, Y01);
                MULCNVL_HX4X2(z2, z3, W0123_0, Y00, Y01);

                MULCNVH_HX4X2(z4, z5, W0123_0, Y01, Y02);
                MULCNVL_HX4X2(z6, z7, W0123_0, Y01, Y02);

                AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
                AE_DSEL16X4(t0, t2, t0, t2, sel);
                AE_DSEL16X4(t1, t3, t1, t3, sel);
                z0 = AE_MOVHALFX4_FROMF16X4(t0);
                z1 = AE_MOVHALFX4_FROMF16X4(t2);
                z2 = AE_MOVHALFX4_FROMF16X4(t1);
                z3 = AE_MOVHALFX4_FROMF16X4(t3);

                AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z4), AE_MOVF16X4_FROMHALFX4(z5), sel);
                AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z6), AE_MOVF16X4_FROMHALFX4(z7), sel);
                AE_DSEL16X4(t0, t2, t0, t2, sel);
                AE_DSEL16X4(t1, t3, t1, t3, sel);
                z4 = AE_MOVHALFX4_FROMF16X4(t0);
                z5 = AE_MOVHALFX4_FROMF16X4(t2);
                z6 = AE_MOVHALFX4_FROMF16X4(t1);
                z7 = AE_MOVHALFX4_FROMF16X4(t3);

                ADD_HX4X2(z0, z4, z0, z4, z1, z5);
                ADD_HX4X2(z2, z6, z2, z6, z3, z7);
                ADD_HX4X2(z0, z4, z0, z4, z2, z6);

                MADDQ_H(z0, z4, Y01, Y02, W4_0);

                AE_SAVHX4X2_XP(z0, z4, alZ, pZ0, 4 * sizeof(xthalf));

            }
            AE_SA128POS_FP(alZ, pZ0);
        }
    }
#undef M
#undef N
}

size_t conv2d_gen_5x5_fp16_getScratchSize(int P, int Q)
{
    int M = 5;
    return (M + 1) * (sizeof(xthalfx8)) + (Q+16)*sizeof(xthalf);
}

#else 
#include "float16.h"

/* simplest reference code */
void conv2d_gen_fp16(void* pScr, float16_t* z, const float16_t* x, const float16_t* y, int M, int N, int P, int Q)
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
}/* conv2d_genf() */

size_t conv2d_gen_5x5_fp16_getScratchSize(int M, int N, int P, int Q)
{
    return 0;
}

void conv2d_gen_5x5_fp16(void* pScr, float16_t* z, const float16_t* x, const float16_t* y, int P, int Q)
{
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT_ALIGN(z, 16);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0 && P % 8 == 0);
    NASSERT(Q >= 0 && Q % 8 == 0);
    conv2d_gen_fp16(pScr, z, x, y, 5, 5, P, Q);
}
#endif
