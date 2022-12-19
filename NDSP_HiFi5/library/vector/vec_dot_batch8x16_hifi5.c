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

#define SMALLER_CODESIZE 1

#if SMALLER_CODESIZE
void vec_dot_batch8x16 (int16_t   *restrict z, const int8_t * restrict x, const cint16ptr_t * restrict y,int rsh, int N, int M)
{
#if 1
  int m, n, sh;

  ae_f64      q0, q1, q2, q3;
  ae_f64      q4, q5, q6, q7;
  ae_int16x4  vai0, vai1;
  ae_int8x8  x0, x1;
  ae_int16x4  y0_0, y0_1, y1_0, y1_1;
  ae_int16x4  y2_0, y2_1, y3_0, y3_1;
  ae_valignx2 ax, ay0, ay1, ay2, ay3;
  ae_valign   az;
  xtbool4     bmask;
  const ae_int8x16 * restrict px;
  const ae_int16x8 * restrict py0;
  const ae_int16x8 * restrict py1;
  const ae_int16x8 * restrict py2;
  const ae_int16x8 * restrict py3;
        ae_int16x4 * restrict pz;

  NASSERT(x);
  NASSERT(y);
  NASSERT(z);

  pz = (ae_int16x4 *)z;
  az = AE_ZALIGN64();

  if (M <= 0) return;

  if (N <= 0)
  {
    vai0 = 0;
    for (m = 0; m < M - 3; m += 4)
    {
      AE_SA16X4_IP(vai0, az, castxcc(ae_int16x4, pz));
    }
    AE_SA64POS_FP(az, castxcc(ae_int16x4, pz));
    for (; m < M; m++)
    {
      z[m] = 0;
    }
    return;
  }
  rsh = 48 - rsh;

  for (m = 0; m < (M >> 2); m++)
  {
    NASSERT(y[4 * m + 0]);
    NASSERT(y[4 * m + 1]);
    NASSERT(y[4 * m + 2]);
    NASSERT(y[4 * m + 3]);

    // Process a batch of 4 vectors;
    px  = (const ae_int8x16 *)x;
    py0 = (const ae_int16x8 *)y[4 * m + 0];
    py1 = (const ae_int16x8 *)y[4 * m + 1];
    py2 = (const ae_int16x8 *)y[4 * m + 2];
    py3 = (const ae_int16x8 *)y[4 * m + 3];
    q0 = AE_ZERO64();
    q1 = AE_ZERO64();
    q2 = AE_ZERO64();
    q3 = AE_ZERO64();
    q4 = AE_ZERO64();
    q5 = AE_ZERO64();
    q6 = AE_ZERO64();
    q7 = AE_ZERO64();
    sh = 0;
    sh = (((uintptr_t)(x)) & 15);
    sh = XT_MAX(0, (XT_MIN(16, N) - sh));
    if (sh > 0)
    {
      sh = sh & 15;

      for (n = 0; n < (sh); n++)
      {
        AE_L8_IP(x0, castxcc(ae_int8, px), sizeof(ae_int8));
        AE_L16_IP(y0_0, castxcc(ae_int16, py0), sizeof(ae_int16));
        AE_L16_IP(y1_0, castxcc(ae_int16, py1), sizeof(ae_int16));
        AE_L16_IP(y2_0, castxcc(ae_int16, py2), sizeof(ae_int16));
        AE_L16_IP(y3_0, castxcc(ae_int16, py3), sizeof(ae_int16));
        bmask = AE_MOVBA4(0xF0 >> 1);
        AE_MOVF16X4(y0_0, AE_ZERO16(), bmask);
        AE_MOVF16X4(y1_0, AE_ZERO16(), bmask);
        AE_MOVF16X4(y2_0, AE_ZERO16(), bmask);
        AE_MOVF16X4(y3_0, AE_ZERO16(), bmask);

        AE_MULAAAA2Q16X8(q0, q1, y0_0, AE_ZERO16(), x0);
        AE_MULAAAA2Q16X8(q2, q3, y1_0, AE_ZERO16(), x0);
        AE_MULAAAA2Q16X8(q4, q5, y2_0, AE_ZERO16(), x0);
        AE_MULAAAA2Q16X8(q6, q7, y3_0, AE_ZERO16(), x0);
      }
    }
    ax  = AE_LA128_PP(px);
    ay0 = AE_LA128_PP(py0);
    ay1 = AE_LA128_PP(py1);
    ay2 = AE_LA128_PP(py2);
    ay3 = AE_LA128_PP(py3);
    for (n = 0; n < ((N - sh) >> 4); n++)
    {
      AE_L8X8X2_IP(x0, x1, px, 2 * sizeof(ae_int8x8));
      AE_LA16X4X2_IP(y0_0, y0_1, ay0, py0);
      AE_LA16X4X2_IP(y1_0, y1_1, ay1, py1);
      AE_LA16X4X2_IP(y2_0, y2_1, ay2, py2);
      AE_LA16X4X2_IP(y3_0, y3_1, ay3, py3);
      AE_MULAAAA2Q16X8(q0, q1, y0_0, y0_1, x0);
      AE_MULAAAA2Q16X8(q2, q3, y1_0, y1_1, x0);
      AE_MULAAAA2Q16X8(q4, q5, y2_0, y2_1, x0);
      AE_MULAAAA2Q16X8(q6, q7, y3_0, y3_1, x0);
      AE_LA16X4X2_IP(y0_0, y0_1, ay0, py0);
      AE_LA16X4X2_IP(y1_0, y1_1, ay1, py1);
      AE_LA16X4X2_IP(y2_0, y2_1, ay2, py2);
      AE_LA16X4X2_IP(y3_0, y3_1, ay3, py3);
      AE_MULAAAA2Q16X8(q0, q1, y0_0, y0_1, x1);
      AE_MULAAAA2Q16X8(q2, q3, y1_0, y1_1, x1);
      AE_MULAAAA2Q16X8(q4, q5, y2_0, y2_1, x1);
      AE_MULAAAA2Q16X8(q6, q7, y3_0, y3_1, x1);
    }

   if (((N - sh) & 15))
   {
     ae_int8x8 tmp;
     sh = ((N - sh) & 15);
     ax = AE_LA128_PP(px);
     AE_LA8X8X2_IP(x0, x1, ax, px);

     ay0 = AE_LA128_PP(py0);
     ay1 = AE_LA128_PP(py1);
     ay2 = AE_LA128_PP(py2);
     ay3 = AE_LA128_PP(py3);
     AE_MOVT8X16_H(x0, tmp, x0, AE_MOVDA8(0), (0xFFFF0000 >> sh));
     AE_MOVT8X16_H(tmp, x1, x1, AE_MOVDA8(0), (0xFFFF0000 >> sh));
     AE_LA16X4X2_IP(y0_0, y0_1, ay0, py0);
     AE_LA16X4X2_IP(y1_0, y1_1, ay1, py1);
     AE_LA16X4X2_IP(y2_0, y2_1, ay2, py2);
     AE_LA16X4X2_IP(y3_0, y3_1, ay3, py3);
     AE_MULAAAA2Q16X8(q0, q1, y0_0, y0_1, x0);
     AE_MULAAAA2Q16X8(q2, q3, y1_0, y1_1, x0);
     AE_MULAAAA2Q16X8(q4, q5, y2_0, y2_1, x0);
     AE_MULAAAA2Q16X8(q6, q7, y3_0, y3_1, x0);
     if (sh>8)
     {
       AE_LA16X4X2_IP(y0_0, y0_1, ay0, py0);
       AE_LA16X4X2_IP(y1_0, y1_1, ay1, py1);
       AE_LA16X4X2_IP(y2_0, y2_1, ay2, py2);
       AE_LA16X4X2_IP(y3_0, y3_1, ay3, py3);
       AE_MULAAAA2Q16X8(q0, q1, y0_0, y0_1, x1);
       AE_MULAAAA2Q16X8(q2, q3, y1_0, y1_1, x1);
       AE_MULAAAA2Q16X8(q4, q5, y2_0, y2_1, x1);
       AE_MULAAAA2Q16X8(q6, q7, y3_0, y3_1, x1);
     }
   }

    q0 = AE_ADD64(q0, q1);
    q1 = AE_ADD64(q2, q3);
    q2 = AE_ADD64(q4, q5);
    q3 = AE_ADD64(q6, q7);
    vai0 = AE_TRUNCA16X4F64S(q0, q0, rsh);
    vai1 = AE_TRUNCA16X4F64S(q1, q1, rsh);
    AE_S16_0_IP(vai0, castxcc(ae_int16, pz), sizeof(ae_int16));
    AE_S16_0_IP(vai1, castxcc(ae_int16, pz), sizeof(ae_int16));
    vai0 = AE_TRUNCA16X4F64S(q2, q2, rsh);
    vai1 = AE_TRUNCA16X4F64S(q3, q3, rsh);
    AE_S16_0_IP(vai0, castxcc(ae_int16, pz), sizeof(ae_int16));
    AE_S16_0_IP(vai1, castxcc(ae_int16, pz), sizeof(ae_int16));
  }
  if (M & 3)
  {
    ae_valign   x_align;
    for (m = (M&(~3)); m < M; m++)
    {
      px = (const ae_int8x16 *)x;
      py0 = (const ae_int16x8 *)y[m];

      x_align = AE_LA64_PP(px);
      ay0 = AE_LA128_PP(py0);
      q0 = q1 = AE_ZERO64();
      for (n = 0; n < (N>>3); n++)
      {
        AE_LA8X8_IP(x0, x_align, castxcc(ae_int8x8, px));
        AE_LA16X4X2_IP(y0_0, y0_1, ay0, py0);
        AE_MULAAAA2Q16X8(q0, q1, y0_0, y0_1, x0);
      }
      // tail
      bmask = AE_MOVBA4(0xF0 >> (N&3));
      if ((N&7)>=4)
      {
        AE_LA8X8_IP(x0, x_align, castxcc(ae_int8x8, px));
        AE_LA16X4X2_IP(y0_0, y0_1, ay0, py0);
        AE_MOVF16X4(y0_1, AE_ZERO16(), bmask);
        AE_MULAAAA2Q16X8(q0, q1, y0_0, y0_1, x0);
      }
      else
      {
        AE_LA8X8_IP(x0, x_align, castxcc(ae_int8x8, px));
        AE_LA16X4X2_IP(y0_0, y0_1, ay0, py0);
        AE_MOVF16X4(y0_0, AE_ZERO16(), bmask);
        AE_MULAAAA2Q16X8(q0, q1, y0_0, AE_ZERO16(), x0);
      }
      q0 = AE_ADD64(q0, q1);
      vai0 = AE_TRUNCA16X4F64S(q0, q0, rsh);
      z[m] = AE_MOVAD16_0(vai0);
    }
  }
#else
    int n;
    int m;
    int Ndiv16, Nmod16;

    const ae_int16x8* restrict y_0;

    const ae_int8x16* restrict _x;

    ae_int16 res16;

    ae_int8x8 zero8;
    ae_int16x4 zero16;

    ae_int64 acc0, acc1;
    ae_int64 acc2, acc3;

    ae_valignx2 uu0,uu3;

    
    ae_int8x8 X0_7, X8_15;
    ae_int16x4 Y0_3, Y4_7, Y8_11, Y12_15;

    ae_int16 * restrict pZ=(ae_int16 *) z;

    Ndiv16 = N>>4;
    Nmod16 = N&15;
    rsh = 48 - rsh;

    zero8  = AE_MOVINT8X8_FROMINT32(0);
    zero16 = AE_MOVINT16X4_FROMINT32(0);

    for (m=0 ; m < M ; m++) 
    {
        acc0=AE_ZERO64();
        acc1=AE_ZERO64();
        acc2=AE_ZERO64();
        acc3=AE_ZERO64();

        y_0 = (const ae_int16x8*) y[m];
        _x  = (const ae_int8x16*) x;
        
        if (Nmod16)
        {
            Y4_7 = zero16;
            __Pragma("no_unroll")
            __Pragma("loop_count min=1,max=15")
            for (n=0; n<Nmod16; n++) 
            {
                AE_L8_IP(X0_7, castxcc(ae_int8,_x), 1);
                AE_L16_IP(Y0_3, castxcc(ae_int16,y_0), 2);
                
                AE_MULAAAA2Q16X8(acc0, acc1, Y0_3, Y4_7, X0_7);
            }
            acc0=AE_SRAI64(acc0,2);
        }

        uu0 = AE_LA128_PP(y_0);
        uu3 = AE_LA128_PP(_x);

        for (n=0 ; n < Ndiv16; n++) // 
        {
            AE_LA8X8X2_IP (X0_7  , X8_15, uu3, _x); 
            AE_LA16X4X2_IP(Y0_3, Y4_7, uu0, y_0); 
            AE_LA16X4X2_IP(Y8_11, Y12_15, uu0, y_0); 

            AE_MULAAAA2Q16X8(acc0, acc1, Y0_3, Y4_7, X0_7);
            AE_MULAAAA2Q16X8(acc2, acc3, Y8_11, Y12_15, X8_15);
        }
       
        acc0 = AE_ADD64(acc0, acc1);
        acc2 = AE_ADD64(acc2, acc3);
        acc0 = AE_ADD64(acc0, acc2);

        res16 = AE_TRUNCA16X4F64S(acc0, acc0, rsh);

        AE_S16_0_IP(res16,pZ,2);
    }
#endif
}
#else
void vec_dot_batch8x16 (int16_t   *z, const int8_t * restrict x, const cint16ptr_t * restrict y,int rsh, int N, int M)
{
    int n;
    int m;
    int Ndiv16, Nmod16;
    int m0,m1,m2;

    const ae_int16x8* restrict y_0;
    const ae_int16x8* restrict y_1;
    const ae_int16x8* restrict y_2;

    const ae_int8x16*restrict _x;


    ae_int16 res16_0,res16_1,res16_2;

    ae_int8x8 zero8;
    ae_int16x4 zero16;

    ae_int64 acc0, acc1;
    ae_int64 acc2, acc3;
    ae_int64 acc4, acc5;

    ae_valignx2 uu0,uu1,uu2,uu3;

    ae_int8x8 X0_7, X8_15;
    ae_int16x4 Y0_3_0, Y0_3_1, Y0_3_2;
    ae_int16x4 Y4_7_0, Y4_7_1, Y4_7_2;
    ae_int16x4 Y8_11_0, Y8_11_1, Y8_11_2;
    ae_int16x4 Y12_15_0, Y12_15_1, Y12_15_2;

    ae_int16 * restrict pZ=(ae_int16 *) z;

    Ndiv16 = N>>4;
    Nmod16 = N&15;
    rsh = 48 - rsh;

    AE_SETCBEGIN0(pZ);
    AE_SETCEND0(pZ+M);

    zero8  = AE_MOVINT8X8_FROMINT32(0);
    zero16 = AE_MOVINT16X4_FROMINT32(0);

    for (m=0 ; m < M ; m+=3) 
    {
        m0 = (m+0) % M;
        m1 = (m+1) % M;
        m2 = (m+2) % M;


        acc0=AE_ZERO64();
        acc1=AE_ZERO64();
        acc2=AE_ZERO64();
        acc3=AE_ZERO64();
        acc4=AE_ZERO64();
        acc5=AE_ZERO64();

        y_0 = (const ae_int16x8*) y[m0];
        y_1 = (const ae_int16x8*) y[m1];
        y_2 = (const ae_int16x8*) y[m2];

        _x  = (const ae_int8x16*) x;
        
        uu0 = AE_LA128_PP(y_0);
        uu1 = AE_LA128_PP(y_1);
        uu2 = AE_LA128_PP(y_2);
        uu3 = AE_LA128_PP(_x);

        for (n=0 ; n < Ndiv16; n++) // 
        {
            AE_LA8X8X2_IP (X0_7  , X8_15, uu3, _x); 
            AE_LA16X4X2_IP(Y0_3_0, Y4_7_0, uu0, y_0); 
            AE_LA16X4X2_IP(Y0_3_1, Y4_7_1, uu1, y_1); 
            AE_LA16X4X2_IP(Y0_3_2, Y4_7_2, uu2, y_2); 

            AE_LA16X4X2_IP(Y8_11_0, Y12_15_0, uu0, y_0); 
            AE_LA16X4X2_IP(Y8_11_1, Y12_15_1, uu1, y_1); 
            AE_LA16X4X2_IP(Y8_11_2, Y12_15_2, uu2, y_2); 

            AE_MULAAAA2Q16X8(acc0, acc1, Y0_3_0, Y4_7_0, X0_7);
            AE_MULAAAA2Q16X8(acc2, acc3, Y0_3_1, Y4_7_1, X0_7);
            AE_MULAAAA2Q16X8(acc4, acc5, Y0_3_2, Y4_7_2, X0_7);

            AE_MULAAAA2Q16X8(acc0, acc1, Y8_11_0, Y12_15_0, X8_15);
            AE_MULAAAA2Q16X8(acc2, acc3, Y8_11_1, Y12_15_1, X8_15);
            AE_MULAAAA2Q16X8(acc4, acc5, Y8_11_2, Y12_15_2, X8_15);
        }

        Y4_7_0 = zero16;
        Y4_7_1 = zero16;
        Y4_7_2 = zero16;

        
        for (n=0 ; n < Nmod16; n++)
        {
            AE_L8_IP(X0_7  ,  castxcc(ae_int8,  _x), 1); 
            AE_L16_IP(Y0_3_0, castxcc(ae_int16,y_0), 2); 
            AE_L16_IP(Y0_3_1, castxcc(ae_int16,y_1), 2); 
            AE_L16_IP(Y0_3_2, castxcc(ae_int16,y_2), 2); 
            
            AE_MOVT8X16_L(X8_15, X0_7  , zero8, X0_7  , 128);
            AE_MOVT16X4(Y0_3_0, zero16, AE_MOVBA4(7));
            AE_MOVT16X4(Y0_3_1, zero16, AE_MOVBA4(7));
            AE_MOVT16X4(Y0_3_2, zero16, AE_MOVBA4(7));

            AE_MULAAAA2Q16X8(acc0, acc1, Y0_3_0, Y4_7_0, X0_7);
            AE_MULAAAA2Q16X8(acc2, acc3, Y0_3_1, Y4_7_1, X0_7);
            AE_MULAAAA2Q16X8(acc4, acc5, Y0_3_2, Y4_7_2, X0_7);
        }

        acc0 = AE_ADD64(acc0, acc1);
        acc2 = AE_ADD64(acc2, acc3);
        acc4 = AE_ADD64(acc4, acc5);

        res16_0 = AE_TRUNCA16X4F64S(acc0, acc0, rsh);
        res16_1 = AE_TRUNCA16X4F64S(acc2, acc2, rsh);
        res16_2 = AE_TRUNCA16X4F64S(acc4, acc4, rsh);
        
        AE_S16_0_XC(res16_0,pZ,2);
        AE_S16_0_XC(res16_1,pZ,2);
        AE_S16_0_XC(res16_2,pZ,2);
    }
}
#endif

