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
#include "vec_alog_table.h"
#include "common.h"
#include "NatureDSP_Signal_math.h"
// code optimized for HiFi5

/*
  NatureDSP Signal Processing Library. Vector Mathematics
   Antilogarithm, natural
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/
/*-------------------------------------------------------------------------
  Antilogarithm
  These routines calculate antilogarithm (base2, natural and base10). 
  Fixed-point functions accept inputs in Q25 and form outputs in Q16.15 
  format and return 0x7FFFFFFF in case of overflow and 0 in case of 
  underflow.

  Precision:
  32x32  32-bit inputs, 32-bit outputs. Accuracy: 8*e-6*y+1LSB
  f      floating point: Accuracy: 2 ULP
  NOTE:
  1.  Although 32 and 24 bit functions provide the similar accuracy, 32-bit
      functions have better input/output resolution (dynamic range).
  2.  Floating point functions are compatible with standard ANSI C routines 
      and set errno and exception flags accordingly.

  Input:
  x[N]  input data,Q25 or floating point 
  N     length of vectors
  Output:
  y[N]  output data,Q16.15 or floating point  

  Restriction:
  x,y should not overlap

  Scalar versions:
  ----------------
  fixed point functions return result in Q16.15

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  x,y - aligned on 16-byte boundary
  N   - multiple of 2
-------------------------------------------------------------------------*/
void vec_antilogn_32x32 (int32_t * restrict y, const int32_t* restrict x, int N)
{
    int n;
          ae_int32x4 * restrict py = (      ae_int32x4 *)y;
    const ae_int32x4 * restrict px = (const ae_int32x4 *)x;
    ae_int32x2 p0,p1,p2,p3,p4;
    ae_valignx2 x_align,y_align;
    p0=AE_L32_I((const ae_int32 *)pow2poly,4*0);
    p1=AE_L32_I((const ae_int32 *)pow2poly,4*1);
    p2=AE_L32_I((const ae_int32 *)pow2poly,4*2);
    p3=AE_L32_I((const ae_int32 *)pow2poly,4*3);
    p4=AE_L32_I((const ae_int32 *)pow2poly,4*4);

    if (N<=0) return;

    x_align = AE_LA128_PP(px);
    y_align = AE_ZALIGN128();
    for (n=0; n<(N>>2); n++)
    {
        ae_int32x2 X0,E0,Y0;
        ae_int32x2 X1,E1,Y1;
        ae_f32x2 t0,t1;
        ae_int64 ah0,al0;
        ae_int64 ah1,al1;

        AE_LA32X2X2_IP(X0,X1, x_align, px);
        AE_MUL32X2S_HH_LL(ah0,al0,1549082005,X0);          
        AE_MUL32X2S_HH_LL(ah1,al1,1549082005,X1);
        X0=AE_TRUNCI32X2F64S(ah0,al0,2); X1=AE_TRUNCI32X2F64S(ah1,al1,2);
        E0=AE_SRAI32(X0,25);             E1=AE_SRAI32(X1,25);
        E0=AE_SUB32(15,E0);              E1=AE_SUB32(15,E1);
        X0=AE_SLLI32(X0,6);              X1=AE_SLLI32(X1,6);
        X0=AE_OR32(X0,0x80000000);       X1=AE_OR32(X1,0x80000000);

        t0=t1=p1; AE_MULAF2P32X4RAS(t0,t1,X0,X1,p0,p0); Y0=t0; Y1=t1;
        t0=t1=p2; AE_MULAF2P32X4RAS(t0,t1,X0,X1,Y0,Y1); Y0=t0; Y1=t1;
        t0=t1=p3; AE_MULAF2P32X4RAS(t0,t1,X0,X1,Y0,Y1); Y0=t0; Y1=t1;
        t0=t1=p4; AE_MULAF2P32X4RAS(t0,t1,X0,X1,Y0,Y1); Y0=t0; Y1=t1;
        Y0=AE_SRAV32RS(Y0,E0);
        Y1=AE_SRAV32RS(Y1,E1);
        AE_SA32X2X2_IP(Y0,Y1, y_align, py);
    }
    AE_SA128POS_FP(y_align, py);
    x += (N&~3);
    y += (N&~3);
    N &= 3;
    if (N>0)
    {
      ae_int32x2 X0, E0, Y0;
      ae_int32x2 X1, E1, Y1;
      ae_f32x2 t0, t1;
      ae_int64 ah0, al0;
      ae_int64 ah1, al1;
      int32_t ALIGN(32) scratch[4];
      ae_int32x4 *pScr;
      pScr = (ae_int32x4*)scratch;
      px = (const ae_int32x4*)x;
      py = (      ae_int32x4*)y;
      AE_S32X2X2_I(0, 0, pScr, 0);
      __Pragma("no_unroll")
      for (n = 0; n<N; n++)
      {
        ae_int32x2 t;
        AE_L32_IP(t, castxcc(ae_int32, px), sizeof(int32_t));
        AE_S32_L_IP(t, castxcc(ae_int32, pScr), sizeof(int32_t));
      }
      pScr = (ae_int32x4*)scratch;
      AE_L32X2X2_I(X0, X1, pScr, 0 * sizeof(ae_int32x4));

      AE_MUL32X2S_HH_LL(ah0, al0, 1549082005, X0);
      AE_MUL32X2S_HH_LL(ah1, al1, 1549082005, X1);
      X0 = AE_TRUNCI32X2F64S(ah0, al0, 2); X1 = AE_TRUNCI32X2F64S(ah1, al1, 2);
      E0 = AE_SRAI32(X0, 25);             E1 = AE_SRAI32(X1, 25);
      E0 = AE_SUB32(15, E0);              E1 = AE_SUB32(15, E1);
      X0 = AE_SLLI32(X0, 6);              X1 = AE_SLLI32(X1, 6);
      X0 = AE_OR32(X0, 0x80000000);       X1 = AE_OR32(X1, 0x80000000);

      t0 = t1 = p1; AE_MULAF2P32X4RAS(t0, t1, X0, X1, p0, p0); Y0 = t0; Y1 = t1;
      t0 = t1 = p2; AE_MULAF2P32X4RAS(t0, t1, X0, X1, Y0, Y1); Y0 = t0; Y1 = t1;
      t0 = t1 = p3; AE_MULAF2P32X4RAS(t0, t1, X0, X1, Y0, Y1); Y0 = t0; Y1 = t1;
      t0 = t1 = p4; AE_MULAF2P32X4RAS(t0, t1, X0, X1, Y0, Y1); Y0 = t0; Y1 = t1;
      Y0 = AE_SRAV32RS(Y0, E0);
      Y1 = AE_SRAV32RS(Y1, E1);
      AE_S32X2X2_I(Y0, Y1, pScr, 0 * sizeof(ae_int32x4));
      __Pragma("no_unroll")
      for (n = 0; n<N; n++)
      {
        ae_int32x2 t;
        AE_L32_IP(t, castxcc(ae_int32, pScr), sizeof(int32_t));
        AE_S32_L_IP(t, castxcc(ae_int32, py), sizeof(int32_t));
      }
    }
} /* vec_antilogn_32x32() */
