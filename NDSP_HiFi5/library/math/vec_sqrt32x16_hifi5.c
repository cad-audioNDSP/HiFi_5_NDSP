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
  NatureDSP Signal Processing Library. Matrix Operations
    Square root 32x16   
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/
#include "NatureDSP_Signal_math.h"
#include "NatureDSP_types.h"
#include "polyrsqrtq23_tbl.h"
#include "common.h"

/*-------------------------------------------------------------------------
  Square Root
  These routines calculate square root.
  NOTE: functions return 0x80000000 on negative argument for 32-bit outputs
  or 0x8000 for 16-bit outputs.
  Two versions of functions available: regular version (vec_sqrt16x16, vec_sqrt32x32, 
  vec_sqrt32x16, vec_sqrt64x32) with arbitrary 
  arguments and faster version (vec_sqrt32x32_fast) that 
  apply some restrictions.

  Precision: 
  16x16  16-bit inputs, 16-bit output. Accuracy: 2LSB
  32x32  32-bit inputs, 32-bit output. Accuracy: (2.6e-7*y+1LSB)
  32x16  32-bit input, 16-bit output.  Accuracy: 2 LSB
  64x32  64-bit inputs, 32-bit output. Accuracy: 2LSB

  Input:
  x[N]  input data, Q15, Q31, Q63 
  N     length of vectors
  Output:
  y[N]  output data, Q15, Q31

  Restriction:
  Regular versions (vec_sqrt16x16, vec_sqrt32x32, vec_sqrt32x16, vec_sqrt64x32):
  x,y - should not overlap

  Faster versions (vec_sqrt32x32_fast):
  x,y - should not overlap
  x,y - aligned on 16-byte boundary
  N   - multiple of 2

  Scalar versions:
  ----------------
  return result, Q15, Q31
-------------------------------------------------------------------------*/
void vec_sqrt32x16      (int16_t* restrict y, const int32_t* restrict x, int N)
{
  int n;
        ae_int16x4 * restrict pY = (      ae_int16x4 *)y;
  const ae_int32x4 * restrict pX = (const ae_int32x4 *)x;
  ae_int32x2 y0, y1, t0, t1;
  ae_int32x2 x0, x1, d0, d1;
  ae_int32x2 c0, c1, c2, c3, c4;
  ae_int16x4 nsa01, h0;
  ae_valignx2 aX;
  ae_valign aY;
  xtbool2 lezero0, lezero1;
  aX = AE_LA128_PP(pX);
  aY = AE_ZALIGN64();

  if (N <= 0) return;

  for (n = 0; n<(N >> 2); n++)
  {
    AE_LA32X2X2_IP(x0, x1, aX, pX);
    /* Normalize x*/
    nsa01 = AE_NSAZ32X4(x0, x1);
    nsa01 = AE_AND16(nsa01, 0xFFFE);
    x0=AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(nsa01)));
    x1=AE_SRAV32RS(x1, AE_SEXT32X2D16_10(AE_NEG16S(nsa01)));

    lezero0 = AE_LT32(x0, 0);
    lezero1 = AE_LT32(x1, 0);

    c0 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 0);
    c1 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 1);
    c2 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 2);
    c3 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 3);
    c4 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 4);

    t0 = t1 = c1; y0 = y1 = c0;
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = c2;
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = c3;
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = c4;
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    /* reiterate rsqrt */
    y0=AE_SLAI32S(y0,11);
    y1=AE_SLAI32S(y1,11);
    AE_MULF2P32X4RAS(d0, d1, y0, y1, y0, y1);

    t0 = t1 = 0x08000000;
    AE_MULSF2P32X4RAS(t0, t1, d0, d1, x0, x1); d0 = t0; d1 = t1;
    AE_MULF2P32X4RAS(d0, d1, d0, d1, y0, y1);
    /* compute sqrt and rescale back */
    AE_MUL2P32X4(d0, d1, d0, d1, 8, 8);
    y0 = AE_ADD32S(y0, d0);
    y1 = AE_ADD32S(y1, d1);
    /* compute sqrt and rescale back */
    AE_MULF2P32X4RAS(y0, y1, x0, x1, y0, y1);
    nsa01 = AE_SRAI16(nsa01, 1);
    nsa01 = AE_SUB16(nsa01, 2);
    y0=AE_SRAV32RS(y0, AE_SEXT32X2D16_32(nsa01));
    y1=AE_SRAV32RS(y1, AE_SEXT32X2D16_10(nsa01));
    AE_MOVT32X2(y0,0x80000000,lezero0);
    AE_MOVT32X2(y1,0x80000000,lezero1);

    h0 = AE_SEL16I(AE_MOVINT16X4_FROMINT32X2(y0), AE_MOVINT16X4_FROMINT32X2(y1), 7);
    AE_SA16X4_IP(h0, aY, pY);
  }
  AE_SA64POS_FP(aY, pY);
  x+=(N&~3);
  y+=(N&~3);
  N&=3;
  if (N>0)
  {
    int32_t ALIGN(32) scratch[4];
    ae_int32x4 *pScr;
    pScr=(      ae_int32x4*)scratch;
    pX  =(const ae_int32x4*)x;
    pY  =(      ae_int16x4*)y;
    AE_S32X2X2_I(0,0,pScr,0);
    __Pragma("no_unroll")
    for(n=0; n<N; n++) 
    {
      ae_int32x2 t;
      AE_L32_IP(t,castxcc(ae_int32,pX),sizeof(int32_t));
      AE_S32_L_IP(t,castxcc(ae_int32,pScr),sizeof(int32_t));
    }
    pScr=(ae_int32x4*)scratch;

    AE_L32X2X2_I (x0,x1,pScr,0*sizeof(ae_int32x4));
    /* Normalize x*/
    nsa01 = AE_NSAZ32X4(x0, x1);
    nsa01 = AE_AND16(nsa01, 0xFFFE);
    x0=AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(nsa01)));
    x1=AE_SRAV32RS(x1, AE_SEXT32X2D16_10(AE_NEG16S(nsa01)));

    lezero0 = AE_LT32(x0, 0);
    lezero1 = AE_LT32(x1, 0);

    c0 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 0);
    c1 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 1);
    c2 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 2);
    c3 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 3);
    c4 = AE_L32_I((const ae_int32 *)polyrsqrtq23, 4 * 4);

    t0 = t1 = c1; y0 = y1 = c0;
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = c2;
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = c3;
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    t0 = t1 = c4;
    AE_MULAF2P32X4RAS(t0, t1, x0, x1, y0, y1); y0 = t0; y1 = t1;
    /* reiterate rsqrt */
    y0 = AE_SLAI32S(y0, 11);
    y1 = AE_SLAI32S(y1, 11);
    AE_MULF2P32X4RAS(d0, d1, y0, y1, y0, y1);

    t0 = t1 = 0x08000000;
    AE_MULSF2P32X4RAS(t0, t1, d0, d1, x0, x1); d0 = t0; d1 = t1;
    AE_MULF2P32X4RAS(d0, d1, d0, d1, y0, y1);
    /* compute sqrt and rescale back */
    AE_MUL2P32X4(d0, d1, d0, d1, 8, 8);
    y0 = AE_ADD32S(y0, d0);
    y1 = AE_ADD32S(y1, d1);
    /* compute sqrt and rescale back */
    AE_MULF2P32X4RAS(y0, y1, x0, x1, y0, y1);
    nsa01 = AE_SRAI16(nsa01, 1);
    nsa01 = AE_SUB16(nsa01, 2);
    y0 = AE_SRAV32RS(y0, AE_SEXT32X2D16_32(nsa01));
    y1 = AE_SRAV32RS(y1, AE_SEXT32X2D16_10(nsa01));
    AE_MOVT32X2(y0, 0x80000000, lezero0);
    AE_MOVT32X2(y1, 0x80000000, lezero1);

    h0 = AE_SEL16I(AE_MOVINT16X4_FROMINT32X2(y0), AE_MOVINT16X4_FROMINT32X2(y1), 7);
    AE_S16X4_I(h0, (ae_int16x4 *) pScr, 0 * sizeof(ae_int16x4));
    __Pragma("no_unroll")
    for (n = 0; n<N; n++)
    {
      ae_int16x4 t;
      AE_L16_IP(t, castxcc(ae_int16, pScr), sizeof(int16_t));
      AE_S16_0_IP(t, castxcc(ae_int16, pY), sizeof(int16_t));
    }
  }
} /* vec_sqrt32x16() */
