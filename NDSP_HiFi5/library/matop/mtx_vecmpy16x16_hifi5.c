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
void mtx_vecmpy16x16(int16_t* restrict z,
    const int16_t* restrict x,
    const int16_t* restrict y,
    int M, int N, int lsh)
{
    ae_valignx2 ax0, ax1;
    ae_valignx2 ay;
    const ae_int16x8 * restrict px0;
    const ae_int16x8 * restrict px1;
    const ae_int16   * restrict py;
    ae_int16x4 X0, X1, X2, X3, Y0, Y1;
    ae_int32x2 Z0;
    ae_int64 B0, B1;
    xtbool4 bmask0, bmask1;
    int m, n;
    if (N <= 0 || M <= 0)    /* exceptional situation */
    {
        for (m = 0; m < M; m++) z[m] = 0;
        return;
    }
    py = (const ae_int16 *)y;
    px0 = (const ae_int16x8 *)(x);
    px1 = (const ae_int16x8 *)(x + N);
    if (M >= 2)
    {
        __Pragma("loop_count min=1")
            for (m = 0; m < (M >> 1); m++)
            {
                py = (const ae_int16 *)y;
                px0 = (const ae_int16x8 *)(x);
                px1 = (const ae_int16x8 *)(x + N);
                ax1 = AE_LA128_PP(px1);
                ax0 = AE_LA128_PP(px0);
                B0 = B1 = AE_ZERO64();

                ay = AE_LA128_PP(py);
                AE_LA16X4X2_IP(Y0, Y1, ay, castxcc(ae_int16x8, py));
                for (n = 0; n < (N - 7); n += 8)
                {
                    AE_LA16X4X2_IP(X0, X1, ax0, px0);
                    AE_LA16X4X2_IP(X2, X3, ax1, px1);

                    AE_MULAAAA2Q16(B0, B1, Y0, Y0, X0, X2);
                    AE_MULAAAA2Q16(B0, B1, Y1, Y1, X1, X3);
                    AE_LA16X4X2_IP(Y0, Y1, ay, castxcc(ae_int16x8, py));
                }
                bmask0 = AE_MOVBA4(0xFF00 >> ((N & 7) + 4));
                bmask1 = AE_MOVBA4(0xFF00 >> (N & 7));
                AE_MOVF16X4(Y0, AE_ZERO16(), bmask0);
                AE_MOVF16X4(Y1, AE_ZERO16(), bmask1);

                AE_LA16X4X2_IP(X0, X1, ax0, px0);
                AE_LA16X4X2_IP(X2, X3, ax1, px1);
                AE_MULAAAA2Q16(B0, B1, Y0, Y0, X0, X2);
                AE_MULAAAA2Q16(B0, B1, Y1, Y1, X1, X3);

                Z0 = AE_TRUNCA32X2F64S(B0, B1, lsh + 33);
                Y0 = AE_ROUND16X4F32SASYM(Z0, Z0);
                z[0] = AE_MOVAD16_3(Y0);
                z[1] = AE_MOVAD16_2(Y0);
                z += 2;
                x = (const int16_t  *)XT_ADDX2(2 * N, (int)(uintptr_t)x);
                //x +=4*N;
            }
    }
    // tail: last up to 1 iteration by m
    if (M & 1)
    {
        py = (const ae_int16 *)y;
        px0 = (const ae_int16x8 *)(x);
        ax0 = AE_LA128_PP(px0);
        B0 = B1 = AE_ZERO64();
        ay = AE_LA128_PP(py);
        for (n = 0; n < (N - 7); n += 8)
        {
            AE_LA16X4X2_IP(Y0, Y1, ay, castxcc(ae_int16x8, py));
            AE_LA16X4X2_IP(X0, X1, ax0, px0);
            AE_MULAAAA2Q16(B0, B1, Y0, Y1, X0, X1);
        }
        bmask0 = AE_MOVBA4(0xFF00 >> ((N & 7) + 4));
        bmask1 = AE_MOVBA4(0xFF00 >> (N & 7));

        AE_LA16X4X2_IP(Y0, Y1, ay, castxcc(ae_int16x8, py));
        AE_LA16X4X2_IP(X0, X1, ax0, px0);

        AE_MOVF16X4(Y0, AE_ZERO16(), bmask0);
        AE_MOVF16X4(Y1, AE_ZERO16(), bmask1);

        AE_MULAAAA2Q16(B0, B1, Y0, Y1, X0, X1);

        Z0 = AE_TRUNCA32X2F64S(B0 + B1, B0 + B1, lsh + 33);
        Y0 = AE_ROUND16X4F32SASYM(Z0, Z0);
        z[0] = AE_MOVAD16_3(Y0);
    }
}
