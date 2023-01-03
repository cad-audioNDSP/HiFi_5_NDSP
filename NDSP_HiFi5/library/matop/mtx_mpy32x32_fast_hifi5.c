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
  Matrix Multiply
  These functions compute the expression z = 2^lsh * x * y for the matrices 
  x and y. The columnar dimension of x must match the row dimension of y. 
  The resulting matrix has the same number of rows as x and the same number 
  of columns as y.
  Transposing API allows to interpret input yt as transposed matrix y.

  NOTE: lsh factor is not relevant for floating point routines.

  Functions require scratch memory for storing intermediate data. This 
  scratch memory area should be aligned on 16 byte boundary and its size is 
  calculated by dedicated scratch allocation functions.

  Two versions of functions available: regular version (mtx_mpy[t]32x32, 
  mtx_mpy[t]16x16, mtx_mpy[t]8x16, mtx_mpy[t]8x8, mtx[t]_mpyf) with 
  arbitrary arguments and faster version (mtx_mpy[t]32x32_fast, 
  mtx_mpy[t]16x16_fast, mtx_mpy[t]8x16_fast, mtx_mpy[t]8x8_fast, 
  mtx_mpy[t]f_fast, cntx_mpyt32x32_fast, cntx_mpytf_fast) that apply 
  some restrictions

  Precision:
  32x32 32-bit inputs, 32-bit output
  16x16 16-bit inputs, 16-bit output
  8x8   8-bit inputs, 8-bit output
  8x16  8/16-bit inputs, 16-bit output
  f     floating point

  Input:
  x[M*N]      input matrix x, Q7, Q15, Q31 or floating point
  y[N*P]      input matrix y, Q7, Q15, Q31 or floating point
  yt[P*N]     transposed input matrix y. Q31,Q15, Q7 floating point. (for 
              transposing API only)
  M           number of rows in matrix x and z
  N           number of columns in matrix x and number of rows in matrix y
  P           number of columns in matrices y and z
  lsh         left shift applied to the result (applied to the fixed-
              point functions only) 
  Output:
  z[M*P]      output matrix z, Q7, Q15, Q31 or floating point 
  Scratch:
  pScr        size in bytes defined by corresponding scratch allocation 
              functions

  Restrictions:
  For regular routines mpy[t]32x32, mtx_mpy[t]16x16, mtx_mpy[t]8x8, 
  mtx_mpy[t]8x16, mtx_mpy[t]f):
  pScr    aligned on 16-byte boundary
  x,y,z   should not overlap

  For faster routines (mtx_mpy[t]32x32_fast, mtx_mpy[t]16x16_fast, 
  mtx_mpy[t]8x8_fast, mtx_mpy[t]8x16_fast, 
  mtx_mpy[t]f_fast):
  x,y,z       should not overlap
  x,y,z,pScr  aligned on 16-byte boundary
  M,N,P       multiplies of 4 for mtx_mpy[t]32x32_fast, mtx_mpy[t]16x16_fast, 
              mtx_mpy[t]8x8_fast, mtx_mpy[t]8x16_fast, mtx_mpy[t]f_fast
              multiplies of 32 for cntx_mpyt32x32_fast, cntx_mpytf_fast
  lsh         should be in range:
              -31...31 for mtx_mpy32x32, mtx_mpy32x32_fast, cntx_mpyt32x32_fast, 
                       cntx_mpytf_fast
              -15...15 for mtx_mpy16x16, mtx_mpy16x16_fast, mtx_mpy[t]8x8, 
                       mtx_mpy[t]8x8_fast, mtx_mpy[t]8x16, 
                       mtx_mpy[t]8x16_fast 

-------------------------------------------------------------------------*/
void mtx_mpy32x32_fast(void* pScr,
    int32_t* restrict z,
    const int32_t* restrict x,
    const int32_t* restrict y,
    int M, int N, int P, int lsh)
{
    int m, n, p;
    const ae_int32x4 * restrict px0;
    const ae_int32x4 * restrict px1;
    const ae_int32x4 * restrict py;
    ae_int32x4 * restrict pz;
    ae_f64 B0,B1,B2,B3,B4,B5,B6,B7;
    ae_f64 D0,D1,D2,D3,D4,D5,D6,D7;
    ae_int32x2 Y0,Y1,Y2,Y3,Y4,Y5,Y6,Y7;
    ae_int32x2 X0,X1,X2,X3,X4,X5,X6,X7, C0,C1,C2,C3,C4,C5,C6,C7;
    ae_int16x4 c0, c1, c2, c3, ind;
    static const short ALIGN(32) dsel_ind[4] = { 1797, 1540, 769, 512 };

    NASSERT((N & 3) == 0);
    NASSERT((M & 3) == 0);
    NASSERT((P & 3) == 0);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(z, HIFI_SIMD_WIDTH);

    if ((M <= 0) || (P <= 0)) return;
    if (N <= 0)
    {
        for (m = 0; m < ((M * P)>>2); m++) AE_S32X2X2_IP(AE_ZERO32(),AE_ZERO32(),castxcc(ae_int32x4,z),sizeof(ae_int32x4));
        return;
    }

    ind = AE_L16X4_I((ae_int16x4*)(dsel_ind), 0);

    __Pragma("loop_count min=1");
    for (p = 0; p < (P >> 2); p++)
    {
        pz = (ae_int32x4 *)z;
        px0 = (const ae_int32x4 *)(x);

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            px1 = (const ae_int32x4 *)XT_ADDX4(N,(uintptr_t)px0);
            py = (const ae_int32x4 *)y;

            AE_MOVDX2(B0,B1,AE_ZERO64(),AE_ZERO64());AE_MOVDX2(B2,B3,AE_ZERO64(),AE_ZERO64());
            AE_MOVDX2(B4,B5,AE_ZERO64(),AE_ZERO64());AE_MOVDX2(B6,B7,AE_ZERO64(),AE_ZERO64());
            AE_MOVDX2(D0,D1,AE_ZERO64(),AE_ZERO64());AE_MOVDX2(D2,D3,AE_ZERO64(),AE_ZERO64());
            AE_MOVDX2(D4,D5,AE_ZERO64(),AE_ZERO64());AE_MOVDX2(D6,D7,AE_ZERO64(),AE_ZERO64());

            __Pragma("loop_count min=1");
            for (n = 0; n < (N >> 2); n++)
            {
                /* load data from 'x' */
                AE_L32X2X2_X (X4, X5, px0, 2*N * sizeof(ae_int32));
                AE_L32X2X2_X (X6, X7, px1, 2*N * sizeof(ae_int32));
                AE_L32X2X2_IP(X0, X1, px0, 2 * sizeof(ae_int32x2));
                AE_L32X2X2_IP(X2, X3, px1, 2 * sizeof(ae_int32x2));
                /* load data from 'y' */
                AE_L32X2X2_XP(Y0, Y1, py, P*sizeof(int32_t));
                AE_L32X2X2_XP(Y2, Y3, py, P*sizeof(int32_t));
                AE_L32X2X2_XP(Y4, Y5, py, P*sizeof(int32_t));
                AE_L32X2X2_XP(Y6, Y7, py, P*sizeof(int32_t));

                C0 = AE_SEL32_HH(Y0, Y2);
                C1 = AE_SEL32_LL(Y0, Y2);
                //C[2] = AE_SEL32_HH(Y[1], Y[3]);
                //C[3] = AE_SEL32_LL(Y[1], Y[3]);
                AE_DSEL16X4(c0, c1, AE_MOVINT16X4_FROMINT32X2(Y1), AE_MOVINT16X4_FROMINT32X2(Y3), ind);
                C2 = AE_MOVINT32X2_FROMINT16X4(c0);
                C3 = AE_MOVINT32X2_FROMINT16X4(c1);

                C4 = AE_SEL32_HH(Y4, Y6);
                C5 = AE_SEL32_LL(Y4, Y6);

                //C[6] = AE_SEL32_HH(Y[5], Y[7]);
                //C[7] = AE_SEL32_LL(Y[5], Y[7]);
                AE_DSEL16X4(c2, c3, AE_MOVINT16X4_FROMINT32X2(Y5), AE_MOVINT16X4_FROMINT32X2(Y7), ind);
                C6 = AE_MOVINT32X2_FROMINT16X4(c2);
                C7 = AE_MOVINT32X2_FROMINT16X4(c3);

                /* perform multiplications */
                AE_MULAAF2D32RA_HH_LL(B0, B1, X0, X1, C0, C5);
                AE_MULAAF2D32RA_HH_LL(B0, B1, X1, X0, C4, C1);
                AE_MULAAF2D32RA_HH_LL(B2, B3, X0, X1, C2, C7);
                AE_MULAAF2D32RA_HH_LL(B2, B3, X1, X0, C6, C3);
                AE_MULAAF2D32RA_HH_LL(B4, B5, X2, X3, C0, C5);
                AE_MULAAF2D32RA_HH_LL(B4, B5, X3, X2, C4, C1);
                AE_MULAAF2D32RA_HH_LL(B6, B7, X2, X3, C2, C7);
                AE_MULAAF2D32RA_HH_LL(B6, B7, X3, X2, C6, C3);

                AE_MULAAF2D32RA_HH_LL(D0, D1, X4, X5, C0, C5);
                AE_MULAAF2D32RA_HH_LL(D0, D1, X5, X4, C4, C1);
                AE_MULAAF2D32RA_HH_LL(D2, D3, X4, X5, C2, C7);
                AE_MULAAF2D32RA_HH_LL(D2, D3, X5, X4, C6, C3);
                AE_MULAAF2D32RA_HH_LL(D4, D5, X6, X7, C0, C5);
                AE_MULAAF2D32RA_HH_LL(D4, D5, X7, X6, C4, C1);
                AE_MULAAF2D32RA_HH_LL(D6, D7, X6, X7, C2, C7);
                AE_MULAAF2D32RA_HH_LL(D6, D7, X7, X6, C6, C3);
            }
            /* format values */
            X0 = AE_TRUNCA32X2F64S(B0, B1, 16 + lsh);
            X1 = AE_TRUNCA32X2F64S(B2, B3, 16 + lsh);
            X2 = AE_TRUNCA32X2F64S(B4, B5, 16 + lsh);
            X3 = AE_TRUNCA32X2F64S(B6, B7, 16 + lsh);
            X4 = AE_TRUNCA32X2F64S(D0, D1, 16 + lsh);
            X5 = AE_TRUNCA32X2F64S(D2, D3, 16 + lsh);
            X6 = AE_TRUNCA32X2F64S(D4, D5, 16 + lsh);
            X7 = AE_TRUNCA32X2F64S(D6, D7, 16 + lsh);
            /* save values */
            AE_S32X2X2_XP(X0, X1, pz, P*sizeof(int32_t));
            AE_S32X2X2_XP(X2, X3, pz, P*sizeof(int32_t));
            AE_S32X2X2_XP(X4, X5, pz, P*sizeof(int32_t));
            AE_S32X2X2_XP(X6, X7, pz, P*sizeof(int32_t));
            px0=(const ae_int32x4*)(const ae_int32x4 *)XT_ADDX4(3*N,(uintptr_t)px0);
        }
        y += 4;
        z += 4;
    }
}

size_t mtx_mpy32x32_fast_getScratchSize(int M, int N, int P)
{
    (void)M; (void)N; (void)P;
    return 0;
}
