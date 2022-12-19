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
/* Code optimized for HiFi5 core */
#include "NatureDSP_Signal_matop.h"
#include "NatureDSP_types.h"
#include "common.h"

/*-------------------------------------------------------------------------
  Matrix by Vector Multiply
  These functions compute the expression z = 2^lsh * x * y for the matrices 
  x and vector y. 
  NOTE: lsh factor is not relevant for floating point routines.

  Two versions of functions available: regular version (mtx_vecmpy32x32,  
  mtx_vecmpy16x16, mtx_vecmpy8x8, mtx_vecmpy8x16, mtx_vecmpyf) with arbitrary 
  arguments and faster version (mtx_vecmpy32x32_fast, mtx_vecmpy16x16_fast, 
  mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast,  mtx_vecmpyf_fast) that apply 
  some restrictions

  Precision: 
  32x32 32-bit input, 32-bit output
  16x16 16-bit input, 16-bit output
  8x8   8-bit inputs, 8-bit output
  8x16  8/16-bit inputs, 16-bit output
  f     floating point

  Input:
  x[M*N] input matrix,Q31,Q15 or floating point
  y[N]   input vector,Q31,Q15 or floating point
  M      number of rows in matrix x
  N      number of columns in matrix x
  lsh    additional left shift(applied to the fixed-
         point functions only) 
  Output:
  z[M]   output vector,Q31,Q15 or floating point

  Restriction:
  For regular routines (mtx_vecmpy32x32, mtx_vecmpy16x16, mtx_vecmpy8x8,
  mtx_vecmpy8x16,  mtx_vecmpyf)
  x,y,z should not overlap

  For faster routines  (mtx_vecmpy32x32_fast, mtx_vecmpy16x16_fast, 
  mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast, mtx_vecmpyf_fast)
  x,y,z   should not overlap
  x,y     aligned on 16-byte boundary
  N, M    multiples of 4
  lsh     should be in range:
          -31...31 for mtx_vecmpy32x32, mtx_vecmpy32x32_fast
          -15...15 for mtx_vecmpy16x16, mtx_vecmpy16x16_fast, 
                   mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast   
-------------------------------------------------------------------------*/
void mtx_vecmpy32x32_fast(int32_t* restrict z,
    const int32_t* restrict x,
    const int32_t* restrict y,
    int M, int N, int lsh)
{
    ae_valignx2 az;
    const ae_int32x4 * restrict px0;
    const ae_int32x4 * restrict px1;
    const ae_int32x4 * restrict py;
    ae_int32x4 * restrict pz;
    ae_int32x2 Y0, Y1, Y2,Y3, X0, X1, X2, X3, X4, X5, X6, X7;
    ae_int32x2                Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7;

    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT((N & 3) == 0);
    NASSERT((M & 3) == 0);
    int m, n;
    if (N <= 0 || M <= 0)    /* exceptional situation */
    {
        for (m = 0; m < (M>>2); m++) AE_S32X2X2_IP(0,0,castxcc(ae_int32x4,z),sizeof(ae_int32x4));
        return;
    }

    az = AE_ZALIGN128();
    px0 = (const ae_int32x4 *)x;
    pz = (ae_int32x4 *)z;
    px1=px0+1;

    if (N&4)
    {
        for(m=0; m<(M>>2); m++)
        {
            ae_f64 A0, A1, A2, A3;
            ae_f64 B0, B1, B2, B3;
            py=(const ae_int32x4*)y;
            AE_L32X2X2_IP(Y0, Y1, py, 4 * sizeof(int32_t));

            AE_L32X2X2_X (X2, X3, px0, 1 * N*sizeof(int32_t));
            AE_L32X2X2_X (X4, X5, px0, 2 * N*sizeof(int32_t));
            AE_L32X2X2_X (X6, X7, px0, 3 * N*sizeof(int32_t));
            AE_L32X2X2_IP(X0, X1, px0, 4 * sizeof(int32_t));

            AE_MULZAAF2D32RA_HH_LL(B0, B1, X0, X3, Y0, Y1);
            AE_MULZAAF2D32RA_HH_LL(B2, B3, X4, X7, Y0, Y1);
            AE_MULZAAF2D32RA_HH_LL(A1, A0, X2, X1, Y0, Y1);
            AE_MULZAAF2D32RA_HH_LL(A3, A2, X6, X5, Y0, Y1);
            px1++;
            __Pragma("ymemory(px1)")
            for (n = 0; n < (N >> 3); n++)
            {
                AE_L32X2X2_I (Y2, Y3, py,   sizeof(ae_int32x4));
                AE_L32X2X2_IP(Y0, Y1, py, 2*sizeof(ae_int32x4));

                AE_L32X2X2_X (X2, X3, px0, 1 * N*sizeof(int32_t));
                AE_L32X2X2_X (X4, X5, px0, 2 * N*sizeof(int32_t));
                AE_L32X2X2_X (X6, X7, px0, 3 * N*sizeof(int32_t));
                AE_L32X2X2_IP(X0, X1, px0, 2*sizeof(ae_int32x4));
                AE_L32X2X2_X (Z2, Z3, px1, 1 * N*sizeof(int32_t));
                AE_L32X2X2_X (Z4, Z5, px1, 2 * N*sizeof(int32_t));
                AE_L32X2X2_X (Z6, Z7, px1, 3 * N*sizeof(int32_t));
                AE_L32X2X2_IP(Z0, Z1, px1, 2*sizeof(ae_int32x4));

                AE_MULAAF2D32RA_HH_LL(A3, A2, X6, X5, Y0, Y1);
                AE_MULAAF2D32RA_HH_LL(B3, B2, Z6, Z5, Y2, Y3);
                AE_MULAAF2D32RA_HH_LL(A2, A3, X4, X7, Y0, Y1);
                AE_MULAAF2D32RA_HH_LL(B2, B3, Z4, Z7, Y2, Y3);
                AE_MULAAF2D32RA_HH_LL(A1, A0, X2, X1, Y0, Y1);
                AE_MULAAF2D32RA_HH_LL(B1, B0, Z2, Z1, Y2, Y3);
                AE_MULAAF2D32RA_HH_LL(A0, A1, X0, X3, Y0, Y1);
                AE_MULAAF2D32RA_HH_LL(B0, B1, Z0, Z3, Y2, Y3);
            }
            px0 += 3 * (N >> 2);
            px1 += 3 * (N >> 2);
            Y0 = AE_TRUNCA32X2F64S(B0+A0, B1+A1, 16 + lsh);
            Y1 = AE_TRUNCA32X2F64S(B2+A2, B3+A3, 16 + lsh);
            AE_SA32X2X2_IP(Y0, Y1, az, pz);
        } 
        AE_SA128POS_FP(az, pz);
    }
    else
    {
        NASSERT(N%8==0);
        for(m=0; m<(M>>2); m++)
        {
            ae_f64 A0, A1, A2, A3;
            ae_f64 B0, B1, B2, B3;
            AE_MOVDX2(A0,A1,0,0); AE_MOVDX2(A2,A3,0,0);
            AE_MOVDX2(B0,B1,0,0); AE_MOVDX2(B2,B3,0,0);
            py=(const ae_int32x4*)y;
            __Pragma("loop_count min=1")
            __Pragma("ymemory(px1)")
            for (n = 0; n < (N >> 3); n++)
            {
                AE_L32X2X2_I (Y2, Y3, py,   sizeof(ae_int32x4));
                AE_L32X2X2_IP(Y0, Y1, py, 2*sizeof(ae_int32x4));

                AE_L32X2X2_X (X2, X3, px0, 1 * N*sizeof(int32_t));
                AE_L32X2X2_X (X4, X5, px0, 2 * N*sizeof(int32_t));
                AE_L32X2X2_X (X6, X7, px0, 3 * N*sizeof(int32_t));
                AE_L32X2X2_IP(X0, X1, px0, 2*sizeof(ae_int32x4));
                AE_L32X2X2_X (Z2, Z3, px1, 1 * N*sizeof(int32_t));
                AE_L32X2X2_X (Z4, Z5, px1, 2 * N*sizeof(int32_t));
                AE_L32X2X2_X (Z6, Z7, px1, 3 * N*sizeof(int32_t));
                AE_L32X2X2_IP(Z0, Z1, px1, 2*sizeof(ae_int32x4));

                AE_MULAAF2D32RA_HH_LL(A3, A2, X6, X5, Y0, Y1);
                AE_MULAAF2D32RA_HH_LL(B3, B2, Z6, Z5, Y2, Y3);
                AE_MULAAF2D32RA_HH_LL(A2, A3, X4, X7, Y0, Y1);
                AE_MULAAF2D32RA_HH_LL(B2, B3, Z4, Z7, Y2, Y3);
                AE_MULAAF2D32RA_HH_LL(A1, A0, X2, X1, Y0, Y1);
                AE_MULAAF2D32RA_HH_LL(B1, B0, Z2, Z1, Y2, Y3);
                AE_MULAAF2D32RA_HH_LL(A0, A1, X0, X3, Y0, Y1);
                AE_MULAAF2D32RA_HH_LL(B0, B1, Z0, Z3, Y2, Y3);
            }
            px0 += 3 * (N >> 2);
            px1 += 3 * (N >> 2);
            Y0 = AE_TRUNCA32X2F64S(B0+A0, B1+A1, 16 + lsh);
            Y1 = AE_TRUNCA32X2F64S(B2+A2, B3+A3, 16 + lsh);
            AE_SA32X2X2_IP(Y0, Y1, az, pz);
        } 
        AE_SA128POS_FP(az, pz);
    }
}
