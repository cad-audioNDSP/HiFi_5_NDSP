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
size_t conv2d_5x5_fp16_getScratchSize (int P, int Q) 
{ 
    (void)P,(void)Q;
    return 0;
}
DISCARD_FUN(void,conv2d_5x5_fp16,(void* pScr, float16_t *z, const float16_t * x, const float16_t * y, int P, int Q))
#elif HAVE_HPFPU
void conv2d_5x5_fp16(void* pScr, float16_t* z, const float16_t* x, const float16_t* y, int P, int Q)
{
    #   define M 5
    #   define N 5
    int i, j;


    const xthalfx8* restrict pY;
    const xthalfx8* restrict pY0;

    const xthalfx8* pW;
    xthalfx8* restrict pZ;

    xthalfx4 S0, S1, S2, S3;
    xthalfx4 Y00, Y01, Y02;
    xthalfx4 Y10, Y11, Y12;
    xthalfx4 Y20, Y21, Y22;
    xthalfx4 Y30, Y31, Y32;
    xthalfx4 Y40, Y41, Y42;
    xthalfx4 zero;
    xthalfx4 W0, W1, W2, W3, W4;
    xthalfx4 dummy0,dummy1,dummy2,dummy3;

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

    zero = CONST_HX4(0);

    
    //load X and replicate
    {
        xthalfx4 *  pWw = (xthalfx4*)(pScr)+(N+1) * M - 1;
        const xthalfx4* pX;
        pX = (xthalfx4*)(x);
        for (i = 0; i < 5; ++i) {
            ae_int16x4 v_i;

	        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(ae_int16));
	        W0 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
            AE_SHX4XP(W0, pWw, -8);

	        AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(ae_int16));
	        W0 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
            AE_SHX4XP(W0, pWw, -8);

            AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(ae_int16));
	        W0 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
            AE_SHX4XP(W0, pWw, -8);

            AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(ae_int16));
	        W0 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
            AE_SHX4XP(W0, pWw, -8);

            AE_L16_IP(v_i, castxcc(const ae_int16, pX), sizeof(ae_int16));
	        W0 = AE_MOVXTHALFX4_FROMINT16X4(v_i);
            AE_SHX4XP(W0, pWw, -8);

            AE_SHX4XP(zero, pWw, -8);
        }

    }
    {
        xthalfx8* pWw = (xthalfx8*)((xthalfx4*)(pScr)+(N+1) * M);
        pW = (xthalfx8*)pScr;
        for (i = 0; i < 5; ++i)
        {
            xthalfx4 w0123;
            AE_LHX4X2_IP(dummy0, W0, pW, sizeof(xthalfx8));
            AE_LHX4X2_IP(W1, W2, pW, sizeof(xthalfx8));
            AE_LHX4X2_IP(W3, W4, pW, sizeof(xthalfx8));
            
            w0123 = AE_SELH_7362(W0, W1);
            w0123 = AE_SELH_7610(w0123, W2);
            w0123 = AE_MOVHALFX4_FROMF16X4(AE_SEL16X4(AE_MOVF16X4_FROMHALFX4(w0123), 
                                                      AE_MOVF16X4_FROMHALFX4(W3),
                                                      AE_MOVF16X4_FROMINT64(0x0007000600050003)));

            AE_SHX4X2_IP(w0123, W4, pWw, sizeof(xthalfx8));
        }

    }

    alZ = AE_ZALIGN128();

    // Processing of convolution

    //first 4 rows 
    //debuged 
    {
        const xthalfx4* pY0;
        const xthalfx4* pY1;
        const xthalfx4* pY2;
        const xthalfx4* pY3;
        xthalfx4* pZ;
        xthalfx4* pZ1;
        xthalfx4* pZ2;
        xthalfx4* pZ3;


        pW = (xthalfx8*)pScr+3;

        pY0 = (xthalfx4*)(y);
        pY1 = (xthalfx4*)((float16_t*)pY0 + Q);
        pY2 = (xthalfx4*)((float16_t*)pY1 + Q);
        pY3 = (xthalfx4*)((float16_t*)pY2 + Q);

        pZ  = (xthalfx4*)(z + (0) * (N + Q - 1));
        pZ1 = (xthalfx4*)(z + (1) * (N + Q - 1));
        pZ2 = (xthalfx4*)(z + (2) * (N + Q - 1));
        pZ3 = (xthalfx4*)(z + (3) * (N + Q - 1));

        Y00 = zero;

        //preload x[3,43210] for row 3 only
        AE_LHX4X2_IP(dummy0, W0, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W1, W2, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W3, W4, pW, sizeof(xthalfx8));
        // all with X[3,43210] (row #3)
        __Pragma("loop_count min=2, factor=2");
        for (j = 0; j < Q; j += 4)
        {
            AE_LHX4IP(Y01, pY0, sizeof(xthalfx4));
            
            S3 = MUL_HX4(W0,Y00);
            MULACNVH_HX4X2(dummy3, S3, W1, Y00, Y01);
            MULACNVL_HX4X2(S3, dummy3, W2, Y00, Y01);
            MULACNVL_HX4X2(dummy3, S3, W3, Y00, Y01);
            MADD_HX4(S3, W4, Y01);
            AE_SHX4IP(S3, pZ3, sizeof(xthalfx4));
            Y00 = Y01;

        }
        //last 4
        {
            Y01 = zero;
            S3 = MUL_HX4(W0, Y00);
            MULACNVH_HX4X2(dummy3, S3, W1, Y00, Y01);
            MULACNVL_HX4X2(S3, dummy3, W2, Y00, Y01);
            MULACNVL_HX4X2(dummy3, S3, W3, Y00, Y01);
            MADD_HX4(S3, W4, Y01);
            AE_SHX4IP(S3, pZ3, sizeof(xthalfx4));
        }

        pY0 = (xthalfx4*)(y);
        Y00 = zero;
        Y10 = zero;
        pZ3 = (xthalfx4*)(z + (3) * (N + Q - 1));
        //preload x[2,43210] for row 3 and 2
        AE_LHX4X2_IP(dummy0, W0, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W1, W2, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W3, W4, pW, sizeof(xthalfx8));
        // all with X[2,43210] 
        __Pragma("loop_count min=2, factor=2");
        for (j = 0; j < Q; j += 4)
        {
            AE_LHX4IP(Y01, pY0, sizeof(xthalfx4));
            AE_LHX4IP(Y11, pY1, sizeof(xthalfx4));
            S2 = zero;
            S3 = AE_LHX4I((xthalfx4*)pZ3, 0);

            MADD_HX4X2(S2, S3, W0, W0, Y00, Y10);

            MULACNVH_HX4X2(dummy2, S2, W1, Y00, Y01);
            MULACNVH_HX4X2(dummy3, S3, W1, Y10, Y11);

            MULACNVL_HX4X2(S2, dummy2, W2, Y00, Y01);
            MULACNVL_HX4X2(S3, dummy3, W2, Y10, Y11);

            MULACNVL_HX4X2(dummy2, S2, W3, Y00, Y01);
            MULACNVL_HX4X2(dummy3, S3, W3, Y10, Y11);

            MADD_HX4X2(S2, S3, W4, W4, Y01, Y11);

            AE_SHX4IP(S2, pZ2, sizeof(xthalfx4));
            AE_SHX4IP(S3, pZ3, sizeof(xthalfx4));

            Y00 = Y01;
            Y10 = Y11;

        }
        //last 4
        {
            Y01 = zero;
            Y11 = zero;

            S2 = zero;
            S3 = AE_LHX4I(pZ3, 0);

            MADD_HX4X2(S2, S3, W0, W0, Y00, Y10);

            MULACNVH_HX4X2(dummy2, S2, W1, Y00, Y01);
            MULACNVH_HX4X2(dummy3, S3, W1, Y10, Y11);

            MULACNVL_HX4X2(S2, dummy2, W2, Y00, Y01);
            MULACNVL_HX4X2(S3, dummy3, W2, Y10, Y11);

            MULACNVL_HX4X2(dummy2, S2, W3, Y00, Y01);
            MULACNVL_HX4X2(dummy3, S3, W3, Y10, Y11);

            MADD_HX4X2(S2, S3, W4, W4, Y01, Y11);

            AE_SHX4IP(S2, pZ2, sizeof(xthalfx4));
            AE_SHX4IP(S3, pZ3, sizeof(xthalfx4));
        }

        pY0 = (xthalfx4*)(y);
        pY1 = (xthalfx4*)((float16_t*)pY0 + Q);
        Y00 = zero;
        Y10 = zero;
        Y20 = zero;
        pZ2 = (xthalfx4*)(z + (2) * (N + Q - 1));
        pZ3 = (xthalfx4*)(z + (3) * (N + Q - 1));
        //preload x[1,43210] for row 3,2,1
        AE_LHX4X2_IP(dummy0, W0, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W1, W2, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W3, W4, pW, sizeof(xthalfx8));
        // all with X[1,43210] 
        __Pragma("loop_count min=2, factor=2");
        for (j = 0; j < Q; j += 4)
        {
            AE_LHX4IP(Y01, pY0, sizeof(xthalfx4));
            AE_LHX4IP(Y11, pY1, sizeof(xthalfx4));
            AE_LHX4IP(Y21, pY2, sizeof(xthalfx4));

            S2 = AE_LHX4I(pZ2, 0);
            S3 = AE_LHX4I(pZ3, 0);

            S1 = MUL_HX4(W0, Y00);
            MADD_HX4X2(S2, S3, W0, W0, Y10, Y20);

            MULACNVH_HX4X2(dummy1, S1, W1, Y00, Y01);
            MULACNVH_HX4X2(dummy2, S2, W1, Y10, Y11);
            MULACNVH_HX4X2(dummy3, S3, W1, Y20, Y21);

            MULACNVL_HX4X2(S1, dummy1, W2, Y00, Y01);
            MULACNVL_HX4X2(S2, dummy2, W2, Y10, Y11);
            MULACNVL_HX4X2(S3, dummy3, W2, Y20, Y21);

            MULACNVL_HX4X2(dummy1, S1, W3, Y00, Y01);
            MULACNVL_HX4X2(dummy2, S2, W3, Y10, Y11);
            MULACNVL_HX4X2(dummy3, S3, W3, Y20, Y21);

            MADD_HX4(S1, W4, Y01);
            MADD_HX4X2(S2, S3, W4, W4, Y11, Y21);

            AE_SHX4IP(S1, pZ1, sizeof(xthalfx4));
            AE_SHX4IP(S2, pZ2, sizeof(xthalfx4));
            AE_SHX4IP(S3, pZ3, sizeof(xthalfx4));

            Y00 = Y01;
            Y10 = Y11;
            Y20 = Y21;

        }
        //last 4
        {
            Y01 = zero;
            Y11 = zero;
            Y21 = zero;

            S2 = AE_LHX4I(pZ2, 0);
            S3 = AE_LHX4I(pZ3, 0);

            S1 = MUL_HX4( W0, Y00);
            MADD_HX4X2(S2, S3, W0, W0, Y10, Y20);

            MULACNVH_HX4X2(dummy1, S1, W1, Y00, Y01);
            MULACNVH_HX4X2(dummy2, S2, W1, Y10, Y11);
            MULACNVH_HX4X2(dummy3, S3, W1, Y20, Y21);

            MULACNVL_HX4X2(S1, dummy1, W2, Y00, Y01);
            MULACNVL_HX4X2(S2, dummy2, W2, Y10, Y11);
            MULACNVL_HX4X2(S3, dummy3, W2, Y20, Y21);

            MULACNVL_HX4X2(dummy1, S1, W3, Y00, Y01);
            MULACNVL_HX4X2(dummy2, S2, W3, Y10, Y11);
            MULACNVL_HX4X2(dummy3, S3, W3, Y20, Y21);

            MADD_HX4(S1, W4, Y01);
            MADD_HX4X2(S2, S3, W4, W4, Y11, Y21);

            AE_SHX4IP(S1, pZ1, sizeof(xthalfx4));
            AE_SHX4IP(S2, pZ2, sizeof(xthalfx4));
            AE_SHX4IP(S3, pZ3, sizeof(xthalfx4));
        }

        pY0 = (xthalfx4*)(y);
        pY1 = (xthalfx4*)((float16_t*)pY0 + Q);
        pY2 = (xthalfx4*)((float16_t*)pY1 + Q);
        Y00 = zero;
        Y10 = zero;
        Y20 = zero;
        Y30 = zero;
        pZ1 = (xthalfx4*)(z + (1) * (N + Q - 1));
        pZ2 = (xthalfx4*)(z + (2) * (N + Q - 1));
        pZ3 = (xthalfx4*)(z + (3) * (N + Q - 1));
        //preload x[0,43210] for row 3,2,1,0
        AE_LHX4X2_IP(dummy0, W0, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W1, W2, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W3, W4, pW, sizeof(xthalfx8));
        // all with X[0,43210] 
        __Pragma("loop_count min=2, factor=2");
        for (j = 0; j < Q; j += 4)
        {
            AE_LHX4IP(Y01, pY0, sizeof(xthalfx4));
            AE_LHX4IP(Y11, pY1, sizeof(xthalfx4));
            AE_LHX4IP(Y21, pY2, sizeof(xthalfx4));
            AE_LHX4IP(Y31, pY3, sizeof(xthalfx4));

            S0 = zero;
            S1 = AE_LHX4I(pZ1, 0);
            S2 = AE_LHX4I(pZ2, 0);
            S3 = AE_LHX4I(pZ3, 0);

            MADD_HX4X2(S0, S1, W0, W0, Y00, Y10);
            MADD_HX4X2(S2, S3, W0, W0, Y20, Y30);

            MULACNVH_HX4X2(dummy0, S0, W1, Y00, Y01);
            MULACNVH_HX4X2(dummy1, S1, W1, Y10, Y11);
            MULACNVH_HX4X2(dummy2, S2, W1, Y20, Y21);
            MULACNVH_HX4X2(dummy3, S3, W1, Y30, Y31);

            MULACNVL_HX4X2(S0, dummy0, W2, Y00, Y01);
            MULACNVL_HX4X2(S1, dummy1, W2, Y10, Y11);
            MULACNVL_HX4X2(S2, dummy2, W2, Y20, Y21);
            MULACNVL_HX4X2(S3, dummy3, W2, Y30, Y31);

            MULACNVL_HX4X2(dummy0, S0, W3, Y00, Y01);
            MULACNVL_HX4X2(dummy1, S1, W3, Y10, Y11);
            MULACNVL_HX4X2(dummy2, S2, W3, Y20, Y21);
            MULACNVL_HX4X2(dummy3, S3, W3, Y30, Y31);

            MADD_HX4X2(S0, S1, W4, W4, Y01, Y11);
            MADD_HX4X2(S2, S3, W4, W4, Y21, Y31);

            AE_SHX4IP(S0, pZ , sizeof(xthalfx4));
            AE_SHX4IP(S1, pZ1, sizeof(xthalfx4));
            AE_SHX4IP(S2, pZ2, sizeof(xthalfx4));
            AE_SHX4IP(S3, pZ3, sizeof(xthalfx4));

            Y00 = Y01;
            Y10 = Y11;
            Y20 = Y21;
            Y30 = Y31;

        }
        //last 4
        {
            Y01 = zero;
            Y11 = zero;
            Y21 = zero;
            Y31 = zero;

            S0 = zero;
            S1 = AE_LHX4I(pZ1, 0);
            S2 = AE_LHX4I(pZ2, 0);
            S3 = AE_LHX4I(pZ3, 0);

            MADD_HX4X2(S0, S1, W0, W0, Y00, Y10);
            MADD_HX4X2(S2, S3, W0, W0, Y20, Y30);

            MULACNVH_HX4X2(dummy0, S0, W1, Y00, Y01);
            MULACNVH_HX4X2(dummy1, S1, W1, Y10, Y11);
            MULACNVH_HX4X2(dummy2, S2, W1, Y20, Y21);
            MULACNVH_HX4X2(dummy3, S3, W1, Y30, Y31);

            MULACNVL_HX4X2(S0, dummy0, W2, Y00, Y01);
            MULACNVL_HX4X2(S1, dummy1, W2, Y10, Y11);
            MULACNVL_HX4X2(S2, dummy2, W2, Y20, Y21);
            MULACNVL_HX4X2(S3, dummy3, W2, Y30, Y31);

            MULACNVL_HX4X2(dummy0, S0, W3, Y00, Y01);
            MULACNVL_HX4X2(dummy1, S1, W3, Y10, Y11);
            MULACNVL_HX4X2(dummy2, S2, W3, Y20, Y21);
            MULACNVL_HX4X2(dummy3, S3, W3, Y30, Y31);

            MADD_HX4X2(S0, S1, W4, W4, Y01, Y11);
            MADD_HX4X2(S2, S3, W4, W4, Y21, Y31);

            AE_SHX4IP(S0, pZ, sizeof(xthalfx4));
            AE_SHX4IP(S1, pZ1, sizeof(xthalfx4));
            AE_SHX4IP(S2, pZ2, sizeof(xthalfx4));
            AE_SHX4IP(S3, pZ3, sizeof(xthalfx4));
        }
    }
    
    ae_int16x4 sel = AE_MOVINT16X4_FROMINT64(0x0706050403020100); // SEL7531 + SEL6420
    int stride = Q*sizeof(xthalf);
    // main loop
    __Pragma("loop_count min=1");
    for (i = 4; i < P; i += 1)
    {
        xthalfx4 z0, z1, z2, z3;
        xthalfx4 r0, r1, r2, r3;
        ae_int16x4 t0, t1, t2, t3;
        ae_int16x4 k0, k1, k2, k3;
        xthalfx4 W0123;

        pW = (xthalfx8*)((xthalfx4*)(pScr)+(N + 1) * M);
        pY = (xthalfx8*)(y + (i - 4) * Q);

        //setup pointers and init variables
        Y00 = zero;
        Y10 = zero;
        Y20 = zero;
        Y30 = zero;
        Y40 = zero;
        pY0 = pY;
        pZ = (xthalfx8*)(z + (i + 0) * (N + Q - 1));

        __Pragma("loop_count min=1");
        for (j = 0; j < Q; j += 8)
        {

            AE_LHX4X2_XP(Y01, Y02, pY0, stride);
            AE_LHX4X2_IP(W0123, W4, pW, sizeof(xthalfx8));
            MULCNVH_HX4X2(z0, z1, W0123, Y00, Y01);
            MULCNVL_HX4X2(z2, z3, W0123, Y00, Y01);
            MULCNVH_HX4X2(r0, r1, W0123, Y01, Y02);
            MULCNVL_HX4X2(r2, r3, W0123, Y01, Y02);
            MUL_HX4X2(S0, S1, W4, W4, Y01, Y02);
            Y00 = Y02;

            AE_LHX4X2_XP(Y11, Y12, pY0, stride);
            AE_LHX4X2_IP(W0123, W4, pW, sizeof(xthalfx8));
            MULACNVH_HX4X2(z0, z1, W0123, Y10, Y11);
            MULACNVL_HX4X2(z2, z3, W0123, Y10, Y11);
            MULACNVH_HX4X2(r0, r1, W0123, Y11, Y12);
            MULACNVL_HX4X2(r2, r3, W0123, Y11, Y12);
            MADD_HX4X2(S0, S1, W4, W4, Y11, Y12);
            Y10 = Y12;

            AE_LHX4X2_XP(Y21, Y22, pY0, stride);
            AE_LHX4X2_IP(W0123, W4, pW, sizeof(xthalfx8));
            MULACNVH_HX4X2(z0, z1, W0123, Y20, Y21);
            MULACNVL_HX4X2(z2, z3, W0123, Y20, Y21);
            MULACNVH_HX4X2(r0, r1, W0123, Y21, Y22);
            MULACNVL_HX4X2(r2, r3, W0123, Y21, Y22);
            MADD_HX4X2(S0, S1, W4, W4, Y21, Y22);
            Y20 = Y22;

            AE_LHX4X2_XP(Y31, Y32, pY0, stride);
            AE_LHX4X2_IP(W0123, W4, pW, sizeof(xthalfx8));
            MULACNVH_HX4X2(z0, z1, W0123, Y30, Y31);
            MULACNVL_HX4X2(z2, z3, W0123, Y30, Y31);
            MULACNVH_HX4X2(r0, r1, W0123, Y31, Y32);
            MULACNVL_HX4X2(r2, r3, W0123, Y31, Y32);
            MADD_HX4X2(S0, S1, W4, W4, Y31, Y32);
            Y30 = Y32;

            AE_LHX4X2_XP(Y41, Y42, pY0, -4*stride+sizeof(xthalfx8));
            AE_LHX4X2_IP(W0123, W4, pW, -4 * (int)sizeof(xthalfx8));
            MULACNVH_HX4X2(z0, z1, W0123, Y40, Y41);
            MULACNVL_HX4X2(z2, z3, W0123, Y40, Y41);
            MULACNVH_HX4X2(r0, r1, W0123, Y41, Y42);
            MULACNVL_HX4X2(r2, r3, W0123, Y41, Y42);
            MADD_HX4X2(S0, S1, W4, W4, Y41, Y42);
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

            AE_SAHX4X2_IP(S0, S1, alZ, pZ);
        }
        AE_SA128POS_FP(alZ, pZ);
        //last 4
        {
            Y01 = zero;
            Y11 = zero;
            Y21 = zero;
            Y31 = zero;
            Y41 = zero;
            //dupclicate x5 times
            {
                AE_LHX4X2_IP(W0123, W4, pW, sizeof(xthalfx8));
                MULCNVH_HX4X2(z0, z1, W0123, Y00, Y01);
                MULCNVL_HX4X2(z2, z3, W0123, Y00, Y01);

                AE_LHX4X2_IP(W0123, W4, pW, sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123, Y10, Y11);
                MULACNVL_HX4X2(z2, z3, W0123, Y10, Y11);

                AE_LHX4X2_IP(W0123, W4, pW, sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123, Y20, Y21);
                MULACNVL_HX4X2(z2, z3, W0123, Y20, Y21);

                AE_LHX4X2_IP(W0123, W4, pW, sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123, Y30, Y31);
                MULACNVL_HX4X2(z2, z3, W0123, Y30, Y31);

                AE_LHX4X2_XP(W0123, W4, pW, -4*(int)sizeof(xthalfx8));
                MULACNVH_HX4X2(z0, z1, W0123, Y40, Y41);
                MULACNVL_HX4X2(z2, z3, W0123, Y40, Y41);
            }

            AE_DSEL16X4(t0, t1, AE_MOVF16X4_FROMHALFX4(z0), AE_MOVF16X4_FROMHALFX4(z1), sel);
            AE_DSEL16X4(t2, t3, AE_MOVF16X4_FROMHALFX4(z2), AE_MOVF16X4_FROMHALFX4(z3), sel);
            AE_DSEL16X4(t0, t2, t0, t2, sel);
            AE_DSEL16X4(t1, t3, t1, t3, sel);
            z0 = AE_MOVHALFX4_FROMF16X4(t0);
            z1 = AE_MOVHALFX4_FROMF16X4(t2);
            z2 = AE_MOVHALFX4_FROMF16X4(t1);
            z3 = AE_MOVHALFX4_FROMF16X4(t3);

            ADD_HX4X2(z0, z1, z0, z1, z2, z3);
            z0 = ADD_HX4(z0, z1);
            AE_SHX4I(z0, (xthalfx4*)pZ, 0);
        }
    }

    //last 4 rows 
    //debuged 
    {
        const xthalfx4* pY0;
        const xthalfx4* pY1;
        const xthalfx4* pY2;
        const xthalfx4* pY3;
        xthalfx4* pZ;
        xthalfx4* pZ1;
        xthalfx4* pZ2;
        xthalfx4* pZ3;

        pW = (xthalfx8*)pScr;

        pY0 = (xthalfx4*)(y + (P - 4) * Q);
        pY1 = (xthalfx4*)((float16_t*)pY0 + Q);
        pY2 = (xthalfx4*)((float16_t*)pY1 + Q);
        pY3 = (xthalfx4*)((float16_t*)pY2 + Q);

        pZ  = (xthalfx4*)(z + (P+0) * (N + Q - 1));
        pZ1 = (xthalfx4*)(z + (P+1) * (N + Q - 1));
        pZ2 = (xthalfx4*)(z + (P+2) * (N + Q - 1));
        pZ3 = (xthalfx4*)(z + (P+3) * (N + Q - 1));

        Y00 = zero;
        Y10 = zero;
        Y20 = zero;
        Y30 = zero;
        //preload x[4,43210]
        AE_LHX4X2_IP(dummy0, W0, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W1, W2, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W3, W4, pW, sizeof(xthalfx8));

        // all with X[4,43210] 
        __Pragma("loop_count min=2, factor=2");
        for (j = 0; j < Q; j += 4)
        {
            AE_LHX4IP(Y01, pY0, sizeof(xthalfx4));
            AE_LHX4IP(Y11, pY1, sizeof(xthalfx4));
            AE_LHX4IP(Y21, pY2, sizeof(xthalfx4));
            AE_LHX4IP(Y31, pY3, sizeof(xthalfx4));

            MUL_HX4X2(S0, S1, W0, W0, Y00, Y10);
            MUL_HX4X2(S2, S3, W0, W0, Y20, Y30);

            MULACNVH_HX4X2(dummy0, S0, W1, Y00, Y01);
            MULACNVH_HX4X2(dummy1, S1, W1, Y10, Y11);
            MULACNVH_HX4X2(dummy2, S2, W1, Y20, Y21);
            MULACNVH_HX4X2(dummy3, S3, W1, Y30, Y31);

            MULACNVL_HX4X2(S0, dummy0, W2, Y00, Y01);
            MULACNVL_HX4X2(S1, dummy1, W2, Y10, Y11);
            MULACNVL_HX4X2(S2, dummy2, W2, Y20, Y21);
            MULACNVL_HX4X2(S3, dummy3, W2, Y30, Y31);

            MULACNVL_HX4X2(dummy0, S0, W3, Y00, Y01);
            MULACNVL_HX4X2(dummy1, S1, W3, Y10, Y11);
            MULACNVL_HX4X2(dummy2, S2, W3, Y20, Y21);
            MULACNVL_HX4X2(dummy3, S3, W3, Y30, Y31);

            MADD_HX4X2(S0, S1, W4, W4, Y01, Y11);
            MADD_HX4X2(S2, S3, W4, W4, Y21, Y31);

            AE_SHX4IP(S0, pZ, sizeof(xthalfx4));
            AE_SHX4IP(S1, pZ1, sizeof(xthalfx4));
            AE_SHX4IP(S2, pZ2, sizeof(xthalfx4));
            AE_SHX4IP(S3, pZ3, sizeof(xthalfx4));

            Y00 = Y01;
            Y10 = Y11;
            Y20 = Y21;
            Y30 = Y31;

        }
        //last 4
        {
            Y01 = zero;
            Y11 = zero;
            Y21 = zero;
            Y31 = zero;

            MUL_HX4X2(S0, S1, W0, W0, Y00, Y10);
            MUL_HX4X2(S2, S3, W0, W0, Y20, Y30);

            MULACNVH_HX4X2(dummy0, S0, W1, Y00, Y01);
            MULACNVH_HX4X2(dummy1, S1, W1, Y10, Y11);
            MULACNVH_HX4X2(dummy2, S2, W1, Y20, Y21);
            MULACNVH_HX4X2(dummy3, S3, W1, Y30, Y31);

            MULACNVL_HX4X2(S0, dummy0, W2, Y00, Y01);
            MULACNVL_HX4X2(S1, dummy1, W2, Y10, Y11);
            MULACNVL_HX4X2(S2, dummy2, W2, Y20, Y21);
            MULACNVL_HX4X2(S3, dummy3, W2, Y30, Y31);

            MULACNVL_HX4X2(dummy0, S0, W3, Y00, Y01);
            MULACNVL_HX4X2(dummy1, S1, W3, Y10, Y11);
            MULACNVL_HX4X2(dummy2, S2, W3, Y20, Y21);
            MULACNVL_HX4X2(dummy3, S3, W3, Y30, Y31);

            MADD_HX4X2(S0, S1, W4, W4, Y01, Y11);
            MADD_HX4X2(S2, S3, W4, W4, Y21, Y31);

            AE_SHX4IP(S0,pZ, sizeof(xthalfx4));
            AE_SHX4IP(S1,pZ1, sizeof(xthalfx4));
            AE_SHX4IP(S2,pZ2, sizeof(xthalfx4));
            AE_SHX4IP(S3,pZ3, sizeof(xthalfx4));
        }

        pY1 = (xthalfx4*)(y + (P - 3) * Q);
        pY2 = (xthalfx4*)((float16_t*)pY1 + Q);
        pY3 = (xthalfx4*)((float16_t*)pY2 + Q);

        pZ  = (xthalfx4*)(z + (P+0) * (N + Q - 1));
        pZ1 = (xthalfx4*)(z + (P+1) * (N + Q - 1));
        pZ2 = (xthalfx4*)(z + (P+2) * (N + Q - 1));

        Y10 = zero;
        Y20 = zero;
        Y30 = zero;

        //preload x[3,43210] 
        AE_LHX4X2_IP(dummy0, W0, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W1, W2, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W3, W4, pW, sizeof(xthalfx8));
        // all with X[3,43210] 
        __Pragma("loop_count min=2, factor=2");
        for (j = 0; j < Q; j += 4)
        {
            AE_LHX4IP(Y11, pY1, sizeof(xthalfx4));
            AE_LHX4IP(Y21, pY2, sizeof(xthalfx4));
            AE_LHX4IP(Y31, pY3, sizeof(xthalfx4));

            S0 = AE_LHX4I(pZ , 0);
            S1 = AE_LHX4I(pZ1, 0);
            S2 = AE_LHX4I(pZ2, 0);

            MADD_HX4(S0, W0, Y10);
            MADD_HX4X2(S1, S2, W0, W0, Y20, Y30);

            MULACNVH_HX4X2(dummy0, S0, W1, Y10, Y11);
            MULACNVH_HX4X2(dummy1, S1, W1, Y20, Y21);
            MULACNVH_HX4X2(dummy2, S2, W1, Y30, Y31);


            MULACNVL_HX4X2(S0, dummy0, W2, Y10, Y11);
            MULACNVL_HX4X2(S1, dummy1, W2, Y20, Y21);
            MULACNVL_HX4X2(S2, dummy2, W2, Y30, Y31);

            MULACNVL_HX4X2(dummy0, S0, W3, Y10, Y11);
            MULACNVL_HX4X2(dummy1, S1, W3, Y20, Y21);
            MULACNVL_HX4X2(dummy2, S2, W3, Y30, Y31);

            MADD_HX4(S0, W4, Y11);
            MADD_HX4X2(S1, S2, W4, W4, Y21, Y31);

            AE_SHX4IP(S0, pZ , sizeof(xthalfx4));
            AE_SHX4IP(S1, pZ1, sizeof(xthalfx4));
            AE_SHX4IP(S2, pZ2, sizeof(xthalfx4));

            Y10 = Y11;
            Y20 = Y21;
            Y30 = Y31;

        }
        //last 4
        {
            Y11 = zero;
            Y21 = zero;
            Y31 = zero;

            S0 = AE_LHX4I((xthalfx4*)pZ, 0);
            S1 = AE_LHX4I((xthalfx4*)pZ1, 0);
            S2 = AE_LHX4I((xthalfx4*)pZ2, 0);

            MADD_HX4(S0, W0, Y10);
            MADD_HX4X2(S1, S2, W0, W0, Y20, Y30);

            MULACNVH_HX4X2(dummy0, S0, W1, Y10, Y11);
            MULACNVH_HX4X2(dummy1, S1, W1, Y20, Y21);
            MULACNVH_HX4X2(dummy2, S2, W1, Y30, Y31);

            MULACNVL_HX4X2(S0, dummy0, W2, Y10, Y11);
            MULACNVL_HX4X2(S1, dummy1, W2, Y20, Y21);
            MULACNVL_HX4X2(S2, dummy2, W2, Y30, Y31);

            MULACNVL_HX4X2(dummy0, S0, W3, Y10, Y11);
            MULACNVL_HX4X2(dummy1, S1, W3, Y20, Y21);
            MULACNVL_HX4X2(dummy2, S2, W3, Y30, Y31);

            MADD_HX4(S0, W4, Y11);
            MADD_HX4X2(S1, S2, W4, W4, Y21, Y31);

            AE_SHX4IP(S0,pZ, sizeof(xthalfx4));
            AE_SHX4IP(S1,pZ1, sizeof(xthalfx4));
            AE_SHX4IP(S2,pZ2, sizeof(xthalfx4));

        }

        pY2 = (xthalfx4*)(y + (P - 2) * Q);
        pY3 = (xthalfx4*)((float16_t*)pY2 + Q);

        pZ = (xthalfx4*)(z + (P + 0) * (N + Q - 1));
        pZ1 = (xthalfx4*)(z + (P + 1) * (N + Q - 1));

        Y20 = zero;
        Y30 = zero;

        //preload x[2,43210]
        AE_LHX4X2_IP(dummy0, W0, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W1, W2, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W3, W4, pW, sizeof(xthalfx8));
        // all with X[2,43210] 
        __Pragma("loop_count min=2, factor=2");
        for (j = 0; j < Q; j += 4)
        {
            AE_LHX4IP(Y21, pY2, sizeof(xthalfx4));
            AE_LHX4IP(Y31, pY3, sizeof(xthalfx4));

            S0 = AE_LHX4I(pZ, 0);
            S1 = AE_LHX4I(pZ1, 0);

            MADD_HX4X2(S0, S1, W0, W0, Y20, Y30);

            MULACNVH_HX4X2(dummy2, S0, W1, Y20, Y21);
            MULACNVH_HX4X2(dummy3, S1, W1, Y30, Y31);

            MULACNVL_HX4X2(S0, dummy2, W2, Y20, Y21);
            MULACNVL_HX4X2(S1, dummy3, W2, Y30, Y31);

            MULACNVL_HX4X2(dummy2, S0, W3, Y20, Y21);
            MULACNVL_HX4X2(dummy3, S1, W3, Y30, Y31);

            MADD_HX4X2(S0, S1, W4, W4, Y21, Y31);

            AE_SHX4IP(S0, pZ , sizeof(xthalfx4));
            AE_SHX4IP(S1, pZ1, sizeof(xthalfx4));

            Y20 = Y21;
            Y30 = Y31;

        }
        //last 4
        {
            Y21 = zero;
            Y31 = zero;

            S0 = AE_LHX4I(pZ, 0);
            S1 = AE_LHX4I(pZ1, 0);

            MADD_HX4X2(S0, S1, W0, W0, Y20, Y30);

            MULACNVH_HX4X2(dummy2, S0, W1, Y20, Y21);
            MULACNVH_HX4X2(dummy3, S1, W1, Y30, Y31);

            MULACNVL_HX4X2(S0, dummy2, W2, Y20, Y21);
            MULACNVL_HX4X2(S1, dummy3, W2, Y30, Y31);

            MULACNVL_HX4X2(dummy2, S0, W3, Y20, Y21);
            MULACNVL_HX4X2(dummy3, S1, W3, Y30, Y31);

            MADD_HX4X2(S0, S1, W4, W4, Y21, Y31);

            AE_SHX4IP(S0, pZ, sizeof(xthalfx4));
            AE_SHX4IP(S1, pZ1, sizeof(xthalfx4));
        }
        

        pY3 = (xthalfx4*)(y + (P - 1) * Q);

        pZ = (xthalfx4*)(z + (P + 0) * (N + Q - 1));
        Y30 = zero;

        //preload x[1,43210]
        AE_LHX4X2_IP(dummy0, W0, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W1, W2, pW, sizeof(xthalfx8));
        AE_LHX4X2_IP(W3, W4, pW, sizeof(xthalfx8));
        // all with X[1,43210]
        __Pragma("loop_count min=2, factor=2");
        for (j = 0; j < Q; j += 4)
        {
            AE_LHX4IP(Y31, pY3, sizeof(xthalfx4));

            S0 = AE_LHX4I((xthalfx4*)pZ, 0);

            MADD_HX4(S0, W0, Y30);
            MULACNVH_HX4X2(dummy3, S0, W1, Y30, Y31);
            MULACNVL_HX4X2(S0, dummy3, W2, Y30, Y31);
            MULACNVL_HX4X2(dummy3, S0, W3, Y30, Y31);
            MADD_HX4(S0, W4, Y31);
            AE_SHX4IP(S0, pZ, sizeof(xthalfx4));
            Y30 = Y31;

        }
        //last 4
        {
            Y31 = zero;

            S0 = AE_LHX4I(pZ, 0);

            MADD_HX4(S0, W0, Y30);
            MULACNVH_HX4X2(dummy3, S0, W1, Y30, Y31);
            MULACNVL_HX4X2(S0, dummy3, W2, Y30, Y31);
            MULACNVL_HX4X2(dummy3, S0, W3, Y30, Y31);
            MADD_HX4(S0, W4, Y31);
            AE_SHX4IP(S0, pZ, sizeof(xthalfx4));
        }



    }

    #   undef M
    #   undef N
}

size_t conv2d_5x5_fp16_getScratchSize(int P, int Q)
{
    int M = 5;
    int N = 5;
    return (M * (N+3)) * (sizeof(xthalfx4));
}

#else 
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

size_t conv2d_5x5_fp16_getScratchSize(int M, int N, int P, int Q)
{
    return 0;
}

void conv2d_5x5_fp16(void* pScr, float16_t* z, const float16_t* x, const float16_t* y, int P, int Q)
{
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT_ALIGN(z, 16);
    NASSERT_ALIGN(pScr, 16);
    NASSERT(P >= 0 && P % 8 == 0);
    NASSERT(Q >= 0 && Q % 8 == 0);
    conv2d_fp16(pScr, z, x, y, 5, 5, P, Q);
}
#endif
