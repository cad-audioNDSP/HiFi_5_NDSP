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
#include "vec_tan32x32_table.h"
#include "common.h"
#include "NatureDSP_Signal_math.h"

/*
  NatureDSP Signal Processing Library. Vector matematics
    Tangent
    Optimized code for HiFi5
  IntegrIT, 2006-2019
*/

/*-------------------------------------------------------------------------
  Tangent 
  Fixed point functions calculate tan(pi*x) for number written in Q31. 
  Floating point functions compute tan(x)
  
  Precision: 
  32x32  32-bit inputs, 32-bit outputs. Accuracy: (1.3e-4*y+1LSB)
                                        if abs(y)<=464873(14.19 in Q15) 
                                        or abs(x)<pi*0.4776
  f      floating point.                Accuracy: 2 ULP

  NOTE:
  1.  Scalar floating point function is compatible with standard ANSI C routines 
      and set errno and exception flags accordingly
  2.  Floating point functions limit the range of allowable input values: [-9099, 9099]
      Whenever the input value does not belong to this range, the result is set to NaN.

  Input:
  x[N]   input data,Q31 or floating point
  N      length of vectors
  Output:
  y[N]   result, Q16.15 or floating point

  Restriction:
  x,y - should not overlap

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  x,z - aligned on 16-byte boundary
  N   - multiple of 2

  Scalar versions:
  ----------------
  Return result, Q16.15 or floating point
-------------------------------------------------------------------------*/

static void mytan32(
        int32_t * scr,
        int32_t * z,
  const int32_t * x,
  int N)
{
  int n;
  const ae_int32x4* restrict pX;
        ae_int32x4* restrict pZ;
  ae_valignx2 aX, aZ;
  const ae_int32x4 * restrict S_rd;
      ae_int32x4 * restrict S_wr;
  ae_int32x2 x0, x1, x2_0, x2_1, r0, r1;
  ae_int32x2 n0, n1, j0, j1, t0, t1;
  ae_int32x2 yt0, yt1, yc0, yc1, sg0, sg1;
  xtbool2 st0, st1;
  ae_int16x4 nsa;

  /* Current block index; overall number of blocks; number of values in the current block */
  int blkLen;
  /* Block size, blkLen <= blkSize */
  const int blkSize = (MAX_ALLOCA_SZ / (3 * (sizeof(int32_t)))&~3);
  /* Allocate a fixed-size scratch area on the stack. */

  NASSERT_ALIGN(scr, 32);
  blkLen = 0;
  for (; N>0; N -= blkLen, x += blkSize, z += blkSize)
  {
    blkLen = XT_MIN(N, blkSize);
    pX = (const ae_int32x4*)x;
    pZ = (      ae_int32x4*)z;
    S_wr = (ae_int32x4*)scr;
    aX = AE_LA128_PP(pX);
    aZ = AE_ZALIGN128();
    for (n = 0; n<(blkLen >> 2); n++)
    {
      AE_LA32X2X2_IP(x0, x1, aX, pX);
      sg0 = x0;
      sg1 = x1;
      x0 = AE_ABS32S(x0);
      x1 = AE_ABS32S(x1);
      /*
      * argument reduction
      */
      n0 = AE_SRAI32(x0, 29);
      n1 = AE_SRAI32(x1, 29);
      n0 = AE_ADD32(n0, 1);
      n1 = AE_ADD32(n1, 1);
      j0 = AE_AND32(n0, -2);
      j1 = AE_AND32(n1, -2);
      AE_MULS2P32X4(x0, x1, j0, j1, 0x20000000, 0x20000000);
      AE_MULF2P32X4RAS(x2_0, x2_1, x0, x1, x0, x1);
      /*Tangent*/
      yt0 = yt1 = AE_L32_I((const ae_int32 *)Pt, 4 * 0);
      t0 = t1 = AE_L32_I((const ae_int32 *)Pt, 4 * 1);
      AE_MULAF2P32X4RAS(t0, t1, x2_0, x2_1, yt0, yt1); yt0 = t0; yt1 = t1;
      t0 = t1 = AE_L32_I((const ae_int32 *)Pt, 4 * 2);
      AE_MULAF2P32X4RAS(t0, t1, x2_0, x2_1, yt0, yt1); yt0 = t0; yt1 = t1;
      t0 = t1 = AE_L32_I((const ae_int32 *)Pt, 4 * 3);
      AE_MULAF2P32X4RAS(t0, t1, x2_0, x2_1, yt0, yt1); yt0 = t0; yt1 = t1;
      t0 = t1 = AE_L32_I((const ae_int32 *)Pt, 4 * 4);
      AE_MULAF2P32X4RAS(t0, t1, x2_0, x2_1, yt0, yt1); yt0 = t0; yt1 = t1;
      t0 = AE_SRAI32(x0, 5);
      t1 = AE_SRAI32(x1, 5);
      AE_MULF2P32X4RAS(yt0, yt1, t0, t1, yt0, yt1); //Q16.15
      /*Cotangent*/
      yc0 = yc1 = AE_L32_I((const ae_int32 *)Qt, 4 * 0);
      t0 = t1 = AE_L32_I((const ae_int32 *)Qt, 4 * 1);
      AE_MULAF2P32X4RAS(t0, t1, x2_0, x2_1, yc0, yc1); yc0 = t0; yc1 = t1;
      t0 = t1 = AE_L32_I((const ae_int32 *)Qt, 4 * 2);
      AE_MULAF2P32X4RAS(t0, t1, x2_0, x2_1, yc0, yc1); yc0 = t0; yc1 = t1;
      t0 = t1 = AE_L32_I((const ae_int32 *)Qt, 4 * 3);
      AE_MULAF2P32X4RAS(t0, t1, x2_0, x2_1, yc0, yc1); yc0 = t0; yc1 = t1;
      t0 = t1 = AE_L32_I((const ae_int32 *)Qt, 4 * 4);
      AE_MULAF2P32X4RAS(t0, t1, x2_0, x2_1, yc0, yc1); yc0 = t0; yc1 = t1;
      t0 = t1 = AE_L32_I((const ae_int32 *)Qt, 4 * 5);
      AE_MULAF2P32X4RAS(t0, t1, x2_0, x2_1, yc0, yc1); yc0 = t0; yc1 = t1;

      AE_S32X2X2_IP(x0, x1, S_wr,  2*sizeof(ae_int32x2));
      AE_S32X2X2_IP(yt0, yt1, S_wr, 2*sizeof(ae_int32x2));
      AE_S32X2X2_IP(yc0, yc1, S_wr, 2*sizeof(ae_int32x2));
    }
    pX = (const ae_int32x4*)x;
    pZ = (ae_int32x4*)z;
    S_rd = (ae_int32x4*)scr;
    aX = AE_LA128_PP(pX);
    aZ = AE_ZALIGN128();
    for (n = 0; n<(blkLen >> 2); n++)
    {
      ae_int32x2 y0, y1, _0x40000000;
      ae_int32x2 e0, e1, xx0, xx1;
      AE_L32X2X2_IP(x0, x1, S_rd, 2*sizeof(ae_int32x2));
      AE_L32X2X2_IP(yt0, yt1, S_rd, 2*sizeof(ae_int32x2));
      AE_L32X2X2_IP(yc0, yc1, S_rd, 2*sizeof(ae_int32x2));
      /* 1/x */
  
      _0x40000000 = AE_MOVDA32(0x40000000);
      nsa = AE_NSA32X4(x0, x1);
      xx0 = AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(nsa)));
      xx1 = AE_SRAV32RS(x1, AE_SEXT32X2D16_10(AE_NEG16S(nsa)));
      xx0 = AE_INT32X2_ABS32S(xx0);
      xx1 = AE_INT32X2_ABS32S(xx1);

      y0 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), xx0);
      y1 = AE_SUB32(AE_MOVDA32((int32_t)0xBAEC0000), xx1);

      t0 = t1 = _0x40000000;
      AE_MULSF2P32X4RAS(t0, t1, xx0, xx1, y0, y1); e0 = t0; e1 = t1;
      e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
      t0 = y0; t1 = y1;
      AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); y0 = t0; y1 = t1;
      t0 = t1 = _0x40000000;
      AE_MULSF2P32X4RAS(t0, t1, xx0, xx1, y0, y1); e0 = t0; e1 = t1;
      e0 = AE_ADD32(e0, e0); e1 = AE_ADD32(e1, e1);
      t0 = y0; t1 = y1;
      AE_MULAF2P32X4RAS(t0, t1, e0, e1, y0, y1); r0 = t0; r1 = t1;

      AE_L32X2X2_X(x0, x1, S_rd, -3*2*(int)sizeof(ae_int32x2));
      yc0=AE_MOVNEG32S_T(yc0,x0);
      yc1=AE_MOVNEG32S_T(yc1,x1);

      AE_MULF2P32X4RAS(yc0, yc1, r0, r1, yc0, yc1); //Q(24-nsa-1)    

      nsa = AE_SUB16(8, nsa);
      yc0 = AE_SRAV32RS(yc0, AE_SEXT32X2D16_32(nsa));
      yc1 = AE_SRAV32RS(yc1, AE_SEXT32X2D16_10(nsa));

      /* adjust sign */
      AE_LA32X2X2_IP(x0, x1, aX, pX);
      sg0 = x0;
      sg1 = x1;
      x0 = AE_ABS32S(x0);
      x1 = AE_ABS32S(x1);
      /*
      * argument reduction
      */

      n0 = AE_ADD32(x0, 0x20000000);
      n1 = AE_ADD32(x1, 0x20000000);
      n0 = AE_SLLI32(n0, 1);
      n1 = AE_SLLI32(n1, 1);
      sg0 = AE_XOR32(sg0, n0);
      sg1 = AE_XOR32(sg1, n1);
      st0 = AE_LT32(n0, 0);
      st1 = AE_LT32(n1, 0);

      AE_MOVT32X2(yt0, yc0, st0);
      AE_MOVT32X2(yt1, yc1, st1);

      t0 = AE_NEG32S(yt0);
      t1 = AE_NEG32S(yt1);
      yt0=AE_MOVNEG32S_T(yt0,sg0);
      yt1=AE_MOVNEG32S_T(yt1,sg1);
      AE_SA32X2X2_IP(yt0, yt1, aZ, pZ);
    } 
    AE_SA128POS_FP(aZ, pZ);
  }
} /* mytan32() */
void vec_tan32x32 ( 
              int32_t* restrict   z, 
              const int32_t *restrict x,
              int N)
{
  const int blkSize = (MAX_ALLOCA_SZ /( (sizeof(int32_t)))&~3);
  /* Allocate a fixed-size scratch area on the stack. */
  int32_t ALIGN(32) scr[blkSize];
  int32_t ALIGN(32) tmpInX[4], tmpOutZ[4];

  int M;
  if ( N<=0 ) return;
  M=N&~3;

  if ( M )
  {
    mytan32(scr, z, x, M);
    z += M;
    x += M;
    N&=3;
  }
  if (N)
  {
    // processing the tail
    int off1,off2;
    ae_int32x2 x0, x1, x2;

    off1=XT_MIN(N-1,1)<<2;
    off2=XT_MIN(N-1,2)<<2;
    x0 = AE_L32_I((const ae_int32*)x, 0);
    x1 = AE_L32_X((const ae_int32*)x, off1);
    x2 = AE_L32_X((const ae_int32*)x, off2);
    AE_S32_L_I(x0, (ae_int32*)tmpInX, 0 * sizeof(int32_t));
    AE_S32_L_I(x1, (ae_int32*)tmpInX, 1 * sizeof(int32_t));
    AE_S32_L_I(x2, (ae_int32*)tmpInX, 2 * sizeof(int32_t));
    

    mytan32(scr, tmpOutZ, tmpInX, 4);

    x0 = AE_L32_I((const ae_int32*)tmpOutZ, 0 * sizeof(int32_t));
    x1 = AE_L32_I((const ae_int32*)tmpOutZ, 1 * sizeof(int32_t));
    x2 = AE_L32_I((const ae_int32*)tmpOutZ,2*sizeof(int32_t));
    AE_S32_L_X(x2,(ae_int32*)z,off2);
    AE_S32_L_X(x1,(ae_int32*)z,off1);
    AE_S32_L_I(x0,(ae_int32*)z,0);

  }
} /* vec_tan32x32() */
