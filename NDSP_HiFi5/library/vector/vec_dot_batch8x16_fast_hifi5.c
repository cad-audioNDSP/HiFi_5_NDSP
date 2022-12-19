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

/* Library API */
#include "NatureDSP_Signal_vector.h"
#include "NatureDSP_types.h"
#include "common.h"

/*-------------------------------------------------------------------------
  Batch Computation of Vector Dot products
  These routines take a set of input vectors and compute their dot product 
  with specific reference data.
  Two versions of routines are available: 
  - regular versions (vec_dot_batch8x8, vec_dot_batch8x16, vec_dot_batch16x16, 
    vec_dot_batchf). They work with arbitratry arguments
  - fast versions (vec_dot_batch8x8_fast, vec_dot_batch8x16_fast, 
    vec_dot_batch16x16_fast, vec_dot_batchf_fast) apply some restrictions.  

  Precision: 
  8x8    8x8-bit data, 16-bit output (fractional multiply Q7xQ7->Q15)
  8x16   8x16-bit data, 16-bit output (fractional multiply Q7xQ15->Q15)
  16x16  16x16-bit data, 32-bit output (fractional multiply Q15xQ15->Q31)
  f      single precision floating point
  fp16   half precision floating point

  Input:
  x[N]     input (reference) data, Q7, Q15 or floating point
  y[M][N]  pointers to M input vectors, Q7, Q15 or floating point
  N        length of vectors
  M        number of vectors
  rsh      right shift for output (for fixed point API only!)
  Output:
  z[M]     dot products between references and M inputs, Q15, Q31 or 
           floating point

  Restrictions:
  Regular versions:
    none
  Faster versions:
    x,y[m] - aligned on 16-byte boundary
    N      - multiple of 8
    M        multiple of 4
-------------------------------------------------------------------------*/

void vec_dot_batch8x16_fast (int16_t   * restrict z, const int8_t * restrict x,const cint16ptr_t *restrict y,int rsh, int N, int M)
{
    int m;
    int n;
    int endN;
    const ae_int16x8* restrict y_0;
    const ae_int16x8* restrict y_1;
    const ae_int16x8* restrict y_2;
    const ae_int16x8* restrict y_3;

    const ae_int8x16* restrict _x;

    ae_int16x4 res16_0;
    ae_int16x4 res16_1;
    ae_int16x4 res;


    ae_int64 acc0, acc4;
    ae_int64 acc1, acc5;
    ae_int64 acc2, acc6;
    ae_int64 acc3, acc7;

    ae_valign uu; // for z

    ae_int8x8 X0_7, X8_15;
    ae_int16x4 Y0_3_0, Y0_3_1, Y0_3_2, Y0_3_3;
    ae_int16x4 Y4_7_0, Y4_7_1, Y4_7_2, Y4_7_3;
    ae_int16x4 Y8_11_0, Y8_11_1, Y8_11_2, Y8_11_3;
    ae_int16x4 Y12_15_0, Y12_15_1, Y12_15_2, Y12_15_3;

    ae_int16x4 * restrict pZ=(ae_int16x4 *) z;


    uu = AE_ZALIGN64();

    NASSERT_ALIGN(x, 16);
    NASSERT(N>=0);
    NASSERT(M>=0);
    NASSERT(N%8==0);
    NASSERT(M%4==0);

    endN = N&(~15);
    rsh = 48 - rsh;

    for (m = 0; m < M; m=m+4) 
    {
        NASSERT_ALIGN(y[m  ], 16);
        NASSERT_ALIGN(y[m+1], 16);
        NASSERT_ALIGN(y[m+2], 16);
        NASSERT_ALIGN(y[m+3], 16);
        acc0 = AE_ZERO64();
        acc1 = AE_ZERO64();
        acc2 = AE_ZERO64();
        acc3 = AE_ZERO64();
        acc4 = AE_ZERO64();
        acc5 = AE_ZERO64();
        acc6 = AE_ZERO64();
        acc7 = AE_ZERO64();
        y_0 = (const ae_int16x8*) y[m  ];
        y_1 = (const ae_int16x8*) y[m+1];
        y_2 = (const ae_int16x8*) y[m+2];
        y_3 = (const ae_int16x8*) y[m+3];
        _x  = (const ae_int8x16*) x;


        for (n=0 ; n < endN; n = n + 16) //5
        {
            AE_L8X8X2_IP (X0_7  , X8_15  ,  _x, 16); 
            AE_L16X4X2_IP (Y0_3_0, Y4_7_0, y_0, 16); 
            AE_L16X4X2_IP (Y0_3_1, Y4_7_1, y_1, 16); 
            AE_L16X4X2_IP (Y0_3_2, Y4_7_2, y_2, 16); 
            AE_L16X4X2_IP (Y0_3_3, Y4_7_3, y_3, 16); 

            AE_L16X4X2_IP (Y8_11_0, Y12_15_0, y_0, 16); 
            AE_L16X4X2_IP (Y8_11_1, Y12_15_1, y_1, 16); 
            AE_L16X4X2_IP (Y8_11_2, Y12_15_2, y_2, 16); 
            AE_L16X4X2_IP (Y8_11_3, Y12_15_3, y_3, 16);  

            AE_MULAAAA2Q16X8(acc0, acc4, Y0_3_0, Y4_7_0, X0_7);
            AE_MULAAAA2Q16X8(acc1, acc5, Y0_3_1, Y4_7_1, X0_7);
            AE_MULAAAA2Q16X8(acc2, acc6, Y0_3_2, Y4_7_2, X0_7);
            AE_MULAAAA2Q16X8(acc3, acc7, Y0_3_3, Y4_7_3, X0_7);

            AE_MULAAAA2Q16X8(acc0, acc4, Y8_11_0, Y12_15_0, X8_15);
            AE_MULAAAA2Q16X8(acc1, acc5, Y8_11_1, Y12_15_1, X8_15);
            AE_MULAAAA2Q16X8(acc2, acc6, Y8_11_2, Y12_15_2, X8_15);
            AE_MULAAAA2Q16X8(acc3, acc7, Y8_11_3, Y12_15_3, X8_15);
        }
        if (N&15) // if N%16!=0 
        {

            X0_7   = AE_L8X8_I ((const ae_int8x8*) _x, 0); 
            AE_L16X4X2_IP (Y0_3_0, Y4_7_0, y_0, 16); 
            AE_L16X4X2_IP (Y0_3_1, Y4_7_1, y_1, 16); 
            AE_L16X4X2_IP (Y0_3_2, Y4_7_2, y_2, 16); 
            AE_L16X4X2_IP (Y0_3_3, Y4_7_3, y_3, 16); 

            AE_MULAAAA2Q16X8(acc0, acc4, Y0_3_0, Y4_7_0, X0_7);
            AE_MULAAAA2Q16X8(acc1, acc5, Y0_3_1, Y4_7_1, X0_7);
            AE_MULAAAA2Q16X8(acc2, acc6, Y0_3_2, Y4_7_2, X0_7);
            AE_MULAAAA2Q16X8(acc3, acc7, Y0_3_3, Y4_7_3, X0_7);
        }

        acc0 = AE_ADD64(acc0, acc4);
        acc1 = AE_ADD64(acc1, acc5);
        acc2 = AE_ADD64(acc2, acc6);
        acc3 = AE_ADD64(acc3, acc7);

        res16_0 = AE_TRUNCA16X4F64S(acc0, acc1, rsh);
        res16_1 = AE_TRUNCA16X4F64S(acc2, acc3, rsh);

        res = AE_SEL16I(res16_0, res16_1,3);

        AE_SA16X4_IP(res, uu, pZ);

    }
    AE_SA64POS_FP(uu,pZ);
}


