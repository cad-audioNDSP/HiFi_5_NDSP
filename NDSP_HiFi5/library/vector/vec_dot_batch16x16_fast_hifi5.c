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

void vec_dot_batch16x16_fast(int32_t * restrict z, const int16_t * restrict x,const cint16ptr_t * restrict y,int rsh, int N, int M)
{
    int m;
    int n;
    const ae_int16x8* restrict y_0;
    const ae_int16x8* restrict y_1;
    const ae_int16x8* restrict y_2;
    const ae_int16x8* restrict y_3;

    const ae_int16x8* restrict _x;

    ae_int32x2 res0;
    ae_int32x2 res1;


    ae_int64 acc0;
    ae_int64 acc1;
    ae_int64 acc2;
    ae_int64 acc3;

    ae_valignx2 uu; // for z

    ae_int16x4 X0123, X4567; 
    ae_int16x4 Y0123_0, Y0123_1, Y0123_2, Y0123_3;
    ae_int16x4 Y4567_0, Y4567_1, Y4567_2, Y4567_3;
    ae_int32x4 * restrict pZ=(ae_int32x4 * )z;

    rsh = 32 - rsh;

    uu = AE_ZALIGN128();

    NASSERT_ALIGN(x, 16);
    NASSERT(N>=0);
    NASSERT(M>=0);
    NASSERT(N%8==0);
    NASSERT(M%4==0);

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
        y_0 = (const ae_int16x8*) y[m  ];
        y_1 = (const ae_int16x8*) y[m+1];
        y_2 = (const ae_int16x8*) y[m+2];
        y_3 = (const ae_int16x8*) y[m+3];
        _x  = (const ae_int16x8*) x;

       
        for (n=0 ; n < N; n = n + 8) 
        {
            AE_L16X4X2_IP (X0123, X4567, _x, 16); 
            AE_L16X4X2_IP (Y0123_0, Y4567_0, y_0, 16); 
            AE_L16X4X2_IP (Y0123_1, Y4567_1, y_1, 16); 
            AE_L16X4X2_IP (Y0123_2, Y4567_2, y_2, 16); 
            AE_L16X4X2_IP (Y0123_3, Y4567_3, y_3, 16); 

            AE_MULAAAA2Q16(acc0, acc1, X0123, X0123, Y0123_0, Y0123_1);
            AE_MULAAAA2Q16(acc2, acc3, X0123, X0123, Y0123_2, Y0123_3);
            AE_MULAAAA2Q16(acc0, acc1, X4567, X4567, Y4567_0, Y4567_1);
            AE_MULAAAA2Q16(acc2, acc3, X4567, X4567, Y4567_2, Y4567_3);

        }

        res0 = AE_TRUNCA32X2F64S(acc0,acc1,rsh);
        res1 = AE_TRUNCA32X2F64S(acc2,acc3,rsh);

        AE_SA32X2X2_IP(res0,res1, uu, pZ);
    }
    AE_SA128POS_FP(uu,pZ);
}


