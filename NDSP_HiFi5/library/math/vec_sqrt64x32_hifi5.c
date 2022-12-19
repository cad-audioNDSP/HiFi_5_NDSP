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
  NatureDSP Signal Processing Library. Vector matematics
    Square Root
    Optimized code for HiFi5
  IntegrIT, 2006-2019
*/

#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_math.h"
#include "scl_sqrt_table.h"
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
void vec_sqrt64x32 (int32_t* restrict y, const int64_t* restrict x, int N)
{
  int n, nsa0, nsa1, nsa2, nsa3, shft0, shft1, shft2, shft3;
  ae_int32x2 nsa01, nsa23;
  const ae_int64   * restrict pX = (const ae_int64   *)x;
        ae_int32x4 * restrict pY = (      ae_int32x4 *)y;
  ae_int64    v0, v1, v2, v3;
  ae_int32x2  x0, x1, d0, d1;
  ae_int32x2  y0, y1, t0, t1;
  ae_int32x2  p0, p1, p2, p3, p4, p5;
  ae_int32x2  viw, vzw;
  xtbool2     inf0, inf1;
  ae_valignx2   aY;

  NASSERT(x);
  NASSERT(y);
  if (N <= 0) return;

  viw = AE_MOVDA32X2(MIN_INT32, MIN_INT32);
  vzw = AE_MOVDA32X2(0, 0);
  aY = AE_ZALIGN128();

  for (n = 0; n < (N >> 2); n++) 
  {
    AE_L64_IP(v0, pX, sizeof(ae_int64));
    AE_L64_IP(v1, pX, sizeof(ae_int64));
    AE_L64_IP(v2, pX, sizeof(ae_int64));
    AE_L64_IP(v3, pX, sizeof(ae_int64));

    nsa0 = AE_NSA64(v0);
    nsa1 = AE_NSA64(v1);
    nsa2 = AE_NSA64(v2);
    nsa3 = AE_NSA64(v3);

    shft0 = (nsa0 & ~1);
    shft1 = (nsa1 & ~1);
    shft2 = (nsa2 & ~1);
    shft3 = (nsa3 & ~1);
     
    nsa01 = AE_MOVDA32X2(nsa0, nsa1);
    nsa23 = AE_MOVDA32X2(nsa2, nsa3);
    nsa01 = AE_SRAI32(nsa01, 1);
    nsa23 = AE_SRAI32(nsa23, 1);
    v0 = AE_SLAA64(v0, shft0);
    v1 = AE_SLAA64(v1, shft1);
    v2 = AE_SLAA64(v2, shft2);
    v3 = AE_SLAA64(v3, shft3);

    x0 = AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(v0), AE_MOVINT32X2_FROMINT64(v1));
    x1 = AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(v2), AE_MOVINT32X2_FROMINT64(v3));

    p0 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 0);
    p1 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 1);
    p2 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 2);
    p3 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 3);
    p4 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 4);
    p5 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 5);

    t0 = t1 = p0; y0 = y1 = p1;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p2;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p3;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p4;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p5;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    t0 = AE_SLAI32S(t0, 4);
    t1 = AE_SLAI32S(t1, 4);

    y0 = y1 = 0x08000000;
    AE_MULF2P32X4RAS(d0, d1, t0, t1, t0, t1);
    AE_MULSF2P32X4RAS(y0, y1, x0, x1, d0, d1); d0 = y0; d1 = y1;
    AE_MULF2P32X4RAS(d0, d1, d0, d1, t0, t1);
    
    AE_MULA2P32X4(t0, t1, d0, d1, 8, 8);

    /* compute sqrt and rescale back */
    AE_MULF2P32X4RAS(y0, y1, t0, t1, x0, x1);
    AE_MUL2P32X4(y0, y1, y0, y1, 4, 4);
    
    d0 = x0; d1 = x1;
    AE_MULSF2P32X4RAS(d0, d1, y0, y1, y0, y1); 
    d0 = AE_ADD32S(d0, d0);
    d1 = AE_ADD32S(d1, d1);
    AE_MULAF2P32X4RAS(y0, y1, d0, d1, t0, t1);

    inf0 = AE_LT32(x0, 0);
    inf1 = AE_LT32(x1, 0);
       
    y0 = AE_SRAV32RS(y0, nsa01);
    y1 = AE_SRAV32RS(y1, nsa23);

    AE_MOVT32X2(y0, viw, inf0);
    AE_MOVT32X2(y1, viw, inf1);
    AE_SA32X2X2_IP(y0, y1, aY, pY);
  }
  AE_SA128POS_FP(aY, pY);
  x += (N&~3);
  y += (N&~3);
  N &= 3;
  if (N>0)
  {
    int32_t ALIGN(32) scratch[4 * 2];
    ae_int32x4 *pScr;
    pScr = (ae_int32x4*)scratch;
    pX = (const ae_int64*)x;
    pY = (ae_int32x4*)y;
    AE_S32X2X2_I(0, 0, pScr, 0);
    AE_S32X2X2_I(0, 0, pScr, 2 * sizeof(int64_t));
    __Pragma("no_unroll")
    for (n = 0; n<N; n++)
    {
      ae_int64 t;
      AE_L64_IP(t, castxcc(ae_int64, pX), sizeof(int64_t));
      AE_S64_IP(t, castxcc(ae_int64, pScr), sizeof(int64_t));
    }
    pScr = (ae_int32x4*)scratch;

    v0=AE_L64_I ((const ae_int64*) pScr, 0*sizeof(ae_int64));
    v1=AE_L64_I ((const ae_int64*) pScr, 1*sizeof(ae_int64));
    v2=AE_L64_I ((const ae_int64*) pScr, 2*sizeof(ae_int64));
    v3=AE_L64_I ((const ae_int64*) pScr, 3*sizeof(ae_int64));

    nsa0 = AE_NSA64(v0);
    nsa1 = AE_NSA64(v1);
    nsa2 = AE_NSA64(v2);
    nsa3 = AE_NSA64(v3);

    shft0 = (nsa0 & ~1);
    shft1 = (nsa1 & ~1);
    shft2 = (nsa2 & ~1);
    shft3 = (nsa3 & ~1);

    nsa01 = AE_MOVDA32X2(nsa0, nsa1);
    nsa23 = AE_MOVDA32X2(nsa2, nsa3);
    nsa01 = AE_SRAI32(nsa01, 1);
    nsa23 = AE_SRAI32(nsa23, 1);
    v0 = AE_SLAA64(v0, shft0);
    v1 = AE_SLAA64(v1, shft1);
    v2 = AE_SLAA64(v2, shft2);
    v3 = AE_SLAA64(v3, shft3);

    x0 = AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(v0), AE_MOVINT32X2_FROMINT64(v1));
    x1 = AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(v2), AE_MOVINT32X2_FROMINT64(v3));

    p0 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 0);
    p1 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 1);
    p2 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 2);
    p3 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 3);
    p4 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 4);
    p5 = AE_L32_I((const ae_int32 *)sqrt_table, 4 * 5);

    t0 = t1 = p0; y0 = y1 = p1;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p2;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p3;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p4;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    y0 = y1 = p5;
    AE_MULAF2P32X4RAS(y0, y1, x0, x1, t0, t1); t0 = y0; t1 = y1;
    t0 = AE_SLAI32S(t0, 4);
    t1 = AE_SLAI32S(t1, 4);

    y0 = y1 = 0x08000000;
    AE_MULF2P32X4RAS(d0, d1, t0, t1, t0, t1);
    AE_MULSF2P32X4RAS(y0, y1, x0, x1, d0, d1); d0 = y0; d1 = y1;
    AE_MULF2P32X4RAS(d0, d1, d0, d1, t0, t1);

    AE_MULA2P32X4(t0, t1, d0, d1, 8, 8);

    /* compute sqrt and rescale back */
    AE_MULF2P32X4RAS(y0, y1, t0, t1, x0, x1);
    AE_MUL2P32X4(y0, y1, y0, y1, 4, 4);

    d0 = x0; d1 = x1;
    AE_MULSF2P32X4RAS(d0, d1, y0, y1, y0, y1);
    d0 = AE_ADD32S(d0, d0);
    d1 = AE_ADD32S(d1, d1);
    AE_MULAF2P32X4RAS(y0, y1, d0, d1, t0, t1);

    inf0 = AE_LT32(x0, 0);
    inf1 = AE_LT32(x1, 0);

    y0 = AE_SRAV32RS(y0, nsa01);
    y1 = AE_SRAV32RS(y1, nsa23);

    AE_MOVT32X2(y0, viw, inf0);
    AE_MOVT32X2(y1, viw, inf1);

    AE_S32X2X2_I(y0, y1, pScr, 0 * sizeof(ae_int32x4));
    __Pragma("no_unroll")
    for (n = 0; n<N; n++)
    {
      ae_int32x2 t;
      AE_L32_IP(t, castxcc(ae_int32, pScr), sizeof(int32_t));
      AE_S32_L_IP(t, castxcc(ae_int32, pY), sizeof(int32_t));
    }
  }
} /* vec_sqrt64x32() */
