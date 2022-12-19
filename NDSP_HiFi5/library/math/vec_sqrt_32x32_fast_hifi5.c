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
#include "sqrt_table.h"
#include "common.h"
#include "NatureDSP_Signal_math.h"
#include "vec_tail32x32.h"

/*
  NatureDSP Signal Processing Library. Vector Mathematics
   Square Root
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/
#define SMALLER_CODESIZE 1 // 1 - simpler variant, 2 - better codesize, but some performance loss
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
void vec_sqrt32x32_fast (int32_t * restrict y, const int32_t * restrict x, int N)
#if SMALLER_CODESIZE==1
{
  /*
  This variant of algorihm computes almost all in 16-bit domain

  algorithmic note:
  we are computing 2 values: w~-0.5/sqrt(x) and z~sqrt(x)
  via iterations:
  w=w+w*(1-2*w*z)
  z=z+w*(x-z*z);
  first approximation is 
  w=1-0.5*x
  z=0.658*x+0.352
  first 2 iterations might be done in 16-bit precision 
  and only the last for z will be done in 32-bit precision
  to achieve ~ 4e-8 relative accuracy
  */
  int n;
  const ae_int32x4* restrict pX;
        ae_int32x4* restrict pY;
  if (N<=0) return;
  pX  =(const ae_int32x4*)x;
  pY  =(      ae_int32x4*)y;
  __Pragma("concurrent")
  for (n=0; n<(N>>3); n++)
  {
    ae_int32x2 x0,x1,x2,x3;
    ae_int16x4 nsa01,nsa23;
    ae_int16x4 y0,z0,w0,d0;
    ae_int16x4 y1,z1,w1,d1;
    ae_int32x2 t0,t1,t2,t3;
    xtbool2 xnotpos0,xnotpos1,xnotpos2,xnotpos3;
    ae_int32x2 c0,c1,c2,c3;
    
    // read by unaligned load and convert endianess
    AE_L32X2X2_IP(x0, x1, pX, 4 * sizeof(ae_int32));
    AE_L32X2X2_IP(x2, x3, pX, 4 * sizeof(ae_int32));

    xnotpos0=AE_LE32(x0,0);
    xnotpos1=AE_LE32(x1,0);
    xnotpos2=AE_LE32(x2,0);
    xnotpos3=AE_LE32(x3,0);
    c0=MIN_INT32; AE_MOVT32X2(c0,0,AE_EQ32(x0,0));
    c1=MIN_INT32; AE_MOVT32X2(c1,0,AE_EQ32(x1,0));
    c2=MIN_INT32; AE_MOVT32X2(c2,0,AE_EQ32(x2,0));
    c3=MIN_INT32; AE_MOVT32X2(c3,0,AE_EQ32(x3,0));

    nsa01=AE_AND16(AE_MOVDA16(~1),AE_NSA32X4(x0,x1));
    nsa23=AE_AND16(AE_MOVDA16(~1),AE_NSA32X4(x2,x3));
    x0=AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(nsa01)));
    x1=AE_SRAV32RS(x1, AE_SEXT32X2D16_10(AE_NEG16S(nsa01)));
    x2=AE_SRAV32RS(x2, AE_SEXT32X2D16_32(AE_NEG16S(nsa23)));
    x3=AE_SRAV32RS(x3, AE_SEXT32X2D16_10(AE_NEG16S(nsa23)));
    y0=AE_TRUNC16X4F32(x0,x1);
    y1=AE_TRUNC16X4F32(x2,x3);
    // first approximation
    w0=AE_ADD16(AE_MOVDA16(0x8000),AE_SRAI16(y0,1));
    w1=AE_ADD16(AE_MOVDA16(0x8000),AE_SRAI16(y1,1));
    z0=AE_MULFD16X16X4RAS(AE_MOVDA16(21561),AE_MOVDA16(11534),y0,AE_MOVDA16(0x7fff));
    z1=AE_MULFD16X16X4RAS(AE_MOVDA16(21561),AE_MOVDA16(11534),y1,AE_MOVDA16(0x7fff));
    // first iteration
    d0=AE_MULFD16X16X4RAS(w0,AE_MOVDA16(0x4000),z0,AE_MOVDA16(0x7fff));
    d1=AE_MULFD16X16X4RAS(w1,AE_MOVDA16(0x4000),z1,AE_MOVDA16(0x7fff));
    d0=AE_ADD16(d0,d0); 
    d1=AE_ADD16(d1,d1); 
    w0=AE_MULFD16X16X4RAS(w0,w0,AE_MOVDA16(0x7fff),d0);
    w1=AE_MULFD16X16X4RAS(w1,w1,AE_MOVDA16(0x7fff),d1);
    d0=AE_MULFD16X16X4RAS(z0,y0,z0,AE_MOVDA16(MIN_INT16));
    d1=AE_MULFD16X16X4RAS(z1,y1,z1,AE_MOVDA16(MIN_INT16));
    z0=AE_MULFD16X16X4RAS(z0,w0,AE_MOVDA16(0x7fff),d0);
    z1=AE_MULFD16X16X4RAS(z1,w1,AE_MOVDA16(0x7fff),d1);
    // second iteration
    d0=AE_MULFD16X16X4RAS(w0,AE_MOVDA16(0x4000),z0,AE_MOVDA16(0x7fff));
    d1=AE_MULFD16X16X4RAS(w1,AE_MOVDA16(0x4000),z1,AE_MOVDA16(0x7fff));
    d0=AE_ADD16(d0,d0); 
    d1=AE_ADD16(d1,d1); 
    w0=AE_MULFD16X16X4RAS(w0,w0,AE_MOVDA16(0x7fff),d0);
    w1=AE_MULFD16X16X4RAS(w1,w1,AE_MOVDA16(0x7fff),d1);
    d0=AE_MULFD16X16X4RAS(z0,y0,z0,AE_MOVDA16(MIN_INT16));
    d1=AE_MULFD16X16X4RAS(z1,y1,z1,AE_MOVDA16(MIN_INT16));
    z0=AE_MULFD16X16X4RAS(z0,w0,AE_MOVDA16(0x7fff),d0);
    z1=AE_MULFD16X16X4RAS(z1,w1,AE_MOVDA16(0x7fff),d1);
    // third short iteration
    AE_CVTI32X4F16(t0,t1,z0,16);
    AE_CVTI32X4F16(t2,t3,z1,16);
    AE_MULSF16X4SS(x0,x1,z0,z0);
    AE_MULSF16X4SS(x2,x3,z1,z1);
    AE_MULSF2P32X16X4RAS(t0,t1,x0,x1,w0);
    AE_MULSF2P32X16X4RAS(t2,t3,x2,x3,w1);
    // output scaling
    nsa01=AE_SRAI16(nsa01,1);
    nsa23=AE_SRAI16(nsa23,1);
    t0=AE_SRAV32RS(t0, AE_SEXT32X2D16_32(nsa01));
    t1=AE_SRAV32RS(t1, AE_SEXT32X2D16_10(nsa01));
    t2=AE_SRAV32RS(t2, AE_SEXT32X2D16_32(nsa23));
    t3=AE_SRAV32RS(t3, AE_SEXT32X2D16_10(nsa23));
    AE_MOVT32X2(t0,c0,xnotpos0);
    AE_MOVT32X2(t1,c1,xnotpos1);
    AE_MOVT32X2(t2,c2,xnotpos2);
    AE_MOVT32X2(t3,c3,xnotpos3);
    AE_S32X2X2_IP(t0, t1, pY, 4 * sizeof(ae_int32));
    AE_S32X2X2_IP(t2, t3, pY, 4 * sizeof(ae_int32));
  }
  N&=7;
  if (N) vec_tail32x32((int32_t*)pY,(int32_t*)pX,vec_sqrt32x32_fast,N);

} /* vec_sqrt32x32_fast() */
#elif SMALLER_CODESIZE==2 && defined (AE_LAV16X4X2_XP) && defined (AE_SAV16X4X2_XP)
{
  /*
  This variant of algorihm computes almost all in 16-bit domain

  algorithmic note:
  we are computing 2 values: w~-0.5/sqrt(x) and z~sqrt(x)
  via iterations:
  w=w+w*(1-2*w*z)
  z=z+w*(x-z*z);
  first approximation is 
  w=1-0.5*x
  z=0.658*x+0.352
  first 2 iterations might be done in 16-bit precision 
  and only the last for z will be done in 32-bit precision
  to achieve ~ 4e-8 relative accuracy
  */
  static const int16_t ALIGN(16) change1632_tbl[]={6|(2<<8), 7|(3<<8), 4|(0<<8), 5|(1<<8)};
  ae_int16x4 change1632;
  ae_valignx2 aX,aZ;
  int n,nbytesRd,nbytesWr;
  const ae_int32x4* restrict pX;
        ae_int32x4* restrict pY;
  if (N<=0) return;
  change1632=AE_L16X4_I((const ae_int16x4*)change1632_tbl,0);
  pX  =(const ae_int32x4*)x;
  pY  =(      ae_int32x4*)y;
  aX=AE_LA128_PP(pX);
  aZ=AE_ZALIGN128();
  nbytesRd=nbytesWr=N<<2;
  __Pragma("concurrent")
  for (n=0; n<((N+7)>>3); n++)
  {
    ae_int16x4 tmp0,tmp1;
    ae_int32x2 x0,x1,x2,x3;
    ae_int16x4 nsa01,nsa23;
    ae_int16x4 y0,z0,w0,d0;
    ae_int16x4 y1,z1,w1,d1;
    ae_int32x2 t0,t1,t2,t3;
    xtbool2 xnotpos0,xnotpos1,xnotpos2,xnotpos3;
    ae_int32x2 c0,c1,c2,c3;
    
    // read by unaligned load and convert endianess
    AE_LAV16X4X2_XP(tmp0, tmp1,aX,castxcc(ae_int16x8,pX), nbytesRd); nbytesRd-=4 * sizeof(ae_int32);
    AE_DSEL16X4(tmp0,tmp1,tmp1,tmp0,change1632);
    x0=AE_MOVINT32X2_FROMINT16X4(tmp0);
    x1=AE_MOVINT32X2_FROMINT16X4(tmp1);
    AE_LAV16X4X2_XP(tmp0, tmp1,aX,castxcc(ae_int16x8,pX), nbytesRd); nbytesRd-=4 * sizeof(ae_int32);
    AE_DSEL16X4(tmp0,tmp1,tmp1,tmp0,change1632);
    x2=AE_MOVINT32X2_FROMINT16X4(tmp0);
    x3=AE_MOVINT32X2_FROMINT16X4(tmp1);

    xnotpos0=AE_LE32(x0,0);
    xnotpos1=AE_LE32(x1,0);
    xnotpos2=AE_LE32(x2,0);
    xnotpos3=AE_LE32(x3,0);
    c0=MIN_INT32; AE_MOVT32X2(c0,0,AE_EQ32(x0,0));
    c1=MIN_INT32; AE_MOVT32X2(c1,0,AE_EQ32(x1,0));
    c2=MIN_INT32; AE_MOVT32X2(c2,0,AE_EQ32(x2,0));
    c3=MIN_INT32; AE_MOVT32X2(c3,0,AE_EQ32(x3,0));

    nsa01=AE_AND16(AE_MOVDA16(~1),AE_NSA32X4(x0,x1));
    nsa23=AE_AND16(AE_MOVDA16(~1),AE_NSA32X4(x2,x3));
    x0=AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(nsa01)));
    x1=AE_SRAV32RS(x1, AE_SEXT32X2D16_10(AE_NEG16S(nsa01)));
    x2=AE_SRAV32RS(x2, AE_SEXT32X2D16_32(AE_NEG16S(nsa23)));
    x3=AE_SRAV32RS(x3, AE_SEXT32X2D16_10(AE_NEG16S(nsa23)));
    y0=AE_TRUNC16X4F32(x0,x1);
    y1=AE_TRUNC16X4F32(x2,x3);
    // first approximation
    w0=AE_ADD16(AE_MOVDA16(0x8000),AE_SRAI16(y0,1));
    w1=AE_ADD16(AE_MOVDA16(0x8000),AE_SRAI16(y1,1));
    z0=AE_MULFD16X16X4RAS(AE_MOVDA16(21561),AE_MOVDA16(11534),y0,AE_MOVDA16(0x7fff));
    z1=AE_MULFD16X16X4RAS(AE_MOVDA16(21561),AE_MOVDA16(11534),y1,AE_MOVDA16(0x7fff));
    // first iteration
    d0=AE_MULFD16X16X4RAS(w0,AE_MOVDA16(0x4000),z0,AE_MOVDA16(0x7fff));
    d1=AE_MULFD16X16X4RAS(w1,AE_MOVDA16(0x4000),z1,AE_MOVDA16(0x7fff));
    d0=AE_ADD16(d0,d0); 
    d1=AE_ADD16(d1,d1); 
    w0=AE_MULFD16X16X4RAS(w0,w0,AE_MOVDA16(0x7fff),d0);
    w1=AE_MULFD16X16X4RAS(w1,w1,AE_MOVDA16(0x7fff),d1);
    d0=AE_MULFD16X16X4RAS(z0,y0,z0,AE_MOVDA16(MIN_INT16));
    d1=AE_MULFD16X16X4RAS(z1,y1,z1,AE_MOVDA16(MIN_INT16));
    z0=AE_MULFD16X16X4RAS(z0,w0,AE_MOVDA16(0x7fff),d0);
    z1=AE_MULFD16X16X4RAS(z1,w1,AE_MOVDA16(0x7fff),d1);
    // second iteration
    d0=AE_MULFD16X16X4RAS(w0,AE_MOVDA16(0x4000),z0,AE_MOVDA16(0x7fff));
    d1=AE_MULFD16X16X4RAS(w1,AE_MOVDA16(0x4000),z1,AE_MOVDA16(0x7fff));
    d0=AE_ADD16(d0,d0); 
    d1=AE_ADD16(d1,d1); 
    w0=AE_MULFD16X16X4RAS(w0,w0,AE_MOVDA16(0x7fff),d0);
    w1=AE_MULFD16X16X4RAS(w1,w1,AE_MOVDA16(0x7fff),d1);
    d0=AE_MULFD16X16X4RAS(z0,y0,z0,AE_MOVDA16(MIN_INT16));
    d1=AE_MULFD16X16X4RAS(z1,y1,z1,AE_MOVDA16(MIN_INT16));
    z0=AE_MULFD16X16X4RAS(z0,w0,AE_MOVDA16(0x7fff),d0);
    z1=AE_MULFD16X16X4RAS(z1,w1,AE_MOVDA16(0x7fff),d1);
    // third short iteration
    AE_CVTI32X4F16(t0,t1,z0,16);
    AE_CVTI32X4F16(t2,t3,z1,16);
    AE_MULSF16X4SS(x0,x1,z0,z0);
    AE_MULSF16X4SS(x2,x3,z1,z1);
    AE_MULSF2P32X16X4RAS(t0,t1,x0,x1,w0);
    AE_MULSF2P32X16X4RAS(t2,t3,x2,x3,w1);
    // output scaling
    nsa01=AE_SRAI16(nsa01,1);
    nsa23=AE_SRAI16(nsa23,1);
    t0=AE_SRAV32RS(t0, AE_SEXT32X2D16_32(nsa01));
    t1=AE_SRAV32RS(t1, AE_SEXT32X2D16_10(nsa01));
    t2=AE_SRAV32RS(t2, AE_SEXT32X2D16_32(nsa23));
    t3=AE_SRAV32RS(t3, AE_SEXT32X2D16_10(nsa23));
    AE_MOVT32X2(t0,c0,xnotpos0);
    AE_MOVT32X2(t1,c1,xnotpos1);
    AE_MOVT32X2(t2,c2,xnotpos2);
    AE_MOVT32X2(t3,c3,xnotpos3);
    // convert endianess and write by unaligned store
    tmp0=AE_MOVINT16X4_FROMINT32X2(t0);
    tmp1=AE_MOVINT16X4_FROMINT32X2(t1);
    AE_DSEL16X4(tmp0,tmp1,tmp1,tmp0,change1632);
    AE_SAV16X4X2_XP(tmp0, tmp1, aZ,castxcc(ae_int16x8,pY), nbytesWr); nbytesWr-=4 * sizeof(ae_int32);
    tmp0=AE_MOVINT16X4_FROMINT32X2(t2);
    tmp1=AE_MOVINT16X4_FROMINT32X2(t3);
    AE_DSEL16X4(tmp0,tmp1,tmp1,tmp0,change1632);
    AE_SAV16X4X2_XP(tmp0, tmp1, aZ,castxcc(ae_int16x8,pY), nbytesWr); nbytesWr-=4 * sizeof(ae_int32);
  }
  AE_SA128POS_FP(aZ,pY);
} /* vec_sqrt32x32_fast() */
#else
{
  /*
  This variant of algorihm computes almost all in 16-bit domain

  algorithmic note:
  we are computing 2 values: w~-0.5/sqrt(x) and z~sqrt(x)
  via iterations:
  w=w+w*(1-2*w*z)
  z=z+w*(x-z*z);
  first approximation is 
  w=1-0.5*x
  z=0.658*x+0.352
  first 2 iterations might be done in 16-bit precision 
  and only the last for z will be done in 32-bit precision
  to achieve ~ 4e-8 relative accuracy
  */
  int n;
  const ae_int32x4* restrict pX;
        ae_int32x4* restrict pY;
  if (N<=0) return;
  pX  =(const ae_int32x4*)x;
  pY  =(      ae_int32x4*)y;
  for (n=0; n<(N>>2); n++)
  {
    ae_int32x2 x0,x1;
    ae_int16x4 nsa01;
    ae_int16x4 y0,z0,w0,d0;
    ae_int32x2 t0,t1;
    xtbool2 xnotpos0,xnotpos1;
    ae_int32x2 c0,c1;

    AE_L32X2X2_IP(x0, x1, pX, 4 * sizeof(ae_int32));
    xnotpos0=AE_LE32(x0,0);
    xnotpos1=AE_LE32(x1,0);
    c0=MIN_INT32; AE_MOVT32X2(c0,0,AE_EQ32(x0,0));
    c1=MIN_INT32; AE_MOVT32X2(c1,0,AE_EQ32(x1,0));

    nsa01=AE_AND16(AE_MOVDA16(~1),AE_NSA32X4(x0,x1));
    x0=AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(nsa01)));
    x1=AE_SRAV32RS(x1, AE_SEXT32X2D16_10(AE_NEG16S(nsa01)));
    y0=AE_TRUNC16X4F32(x0,x1);
    // first approximation
    w0=AE_ADD16(AE_MOVDA16(0x8000),AE_SRAI16(y0,1));
    z0=AE_MULFD16X16X4RAS(AE_MOVDA16(21561),AE_MOVDA16(11534),y0,AE_MOVDA16(0x7fff));
    // first iteration
    d0=AE_MULFD16X16X4RAS(w0,AE_MOVDA16(0x4000),z0,AE_MOVDA16(0x7fff));
    d0=AE_ADD16(d0,d0); 
    w0=AE_MULFD16X16X4RAS(w0,w0,AE_MOVDA16(0x7fff),d0);
    d0=AE_MULFD16X16X4RAS(z0,y0,z0,AE_MOVDA16(MIN_INT16));
    z0=AE_MULFD16X16X4RAS(z0,w0,AE_MOVDA16(0x7fff),d0);
    // second iteration
    d0=AE_MULFD16X16X4RAS(w0,AE_MOVDA16(0x4000),z0,AE_MOVDA16(0x7fff));
    d0=AE_ADD16(d0,d0); 
    w0=AE_MULFD16X16X4RAS(w0,w0,AE_MOVDA16(0x7fff),d0);
    d0=AE_MULFD16X16X4RAS(z0,y0,z0,AE_MOVDA16(MIN_INT16));
    z0=AE_MULFD16X16X4RAS(z0,w0,AE_MOVDA16(0x7fff),d0);
    // third short iteration
    AE_CVTI32X4F16(t0,t1,z0,16);
    AE_MULSF16X4SS(x0,x1,z0,z0);
    AE_MULSF2P32X16X4RAS(t0,t1,x0,x1,w0);
    // output scaling
    nsa01=AE_SRAI16(nsa01,1);
    t0=AE_SRAV32RS(t0, AE_SEXT32X2D16_32(nsa01));
    t1=AE_SRAV32RS(t1, AE_SEXT32X2D16_10(nsa01));
    AE_MOVT32X2(t0,c0,xnotpos0);
    AE_MOVT32X2(t1,c1,xnotpos1);
    AE_S32X2X2_IP(t0, t1, pY, 4 * sizeof(ae_int32));
  }
  __Pragma("no_reorder");
  x += (N&~3);
  y+=(N&~3);
  if (N&3)
  {
    ae_int32x2 x0,x1;
    ae_int16x4 nsa01;
    ae_int16x4 y0,z0,w0,d0;
    ae_int32x2 t0;
    xtbool2 xnotpos0;
    ae_int32x2 c0;
    pX  =(const ae_int32x4*)x;
    pY  =(      ae_int32x4*)y;
    AE_L32X2_IP(x0, castxcc(ae_int32x2, pX), 2 * sizeof(ae_int32));

    // set masks for outbound conditions
    xnotpos0=AE_LE32(x0,0);
    c0=MIN_INT32; AE_MOVT32X2(c0,0,AE_EQ32(x0,0));
    // normalization
    nsa01=AE_AND16(AE_MOVDA16(~1),AE_NSA32X4(x0,x0));
    x0=AE_SRAV32RS(x0, AE_SEXT32X2D16_32(AE_NEG16S(nsa01)));
    y0=AE_TRUNC16X4F32(x0,x0);
    // first approximation
    w0=AE_ADD16(AE_MOVDA16(0x8000),AE_SRAI16(y0,1));
    z0=AE_MULFD16X16X4RAS(AE_MOVDA16(21561),AE_MOVDA16(11534),y0,AE_MOVDA16(0x7fff));
    // first iteration
    d0=AE_MULFD16X16X4RAS(w0,AE_MOVDA16(0x4000),z0,AE_MOVDA16(0x7fff));
    d0=AE_ADD16(d0,d0); 
    w0=AE_MULFD16X16X4RAS(w0,w0,AE_MOVDA16(0x7fff),d0);
    d0=AE_MULFD16X16X4RAS(z0,y0,z0,AE_MOVDA16(MIN_INT16));
    z0=AE_MULFD16X16X4RAS(z0,w0,AE_MOVDA16(0x7fff),d0);
    // second iteration
    d0=AE_MULFD16X16X4RAS(w0,AE_MOVDA16(0x4000),z0,AE_MOVDA16(0x7fff));
    d0=AE_ADD16(d0,d0); 
    w0=AE_MULFD16X16X4RAS(w0,w0,AE_MOVDA16(0x7fff),d0);
    d0=AE_MULFD16X16X4RAS(z0,y0,z0,AE_MOVDA16(MIN_INT16));
    z0=AE_MULFD16X16X4RAS(z0,w0,AE_MOVDA16(0x7fff),d0);
    // third short iteration
    AE_CVTI32X4F16(t0,t0,z0,16);

    AE_MULSF16X4SS(x0,x1,z0,z0);
    AE_MULSF2P32X16X4RAS(t0,t0,x0,x0,w0);
    // output scaling
    nsa01=AE_SRAI16(nsa01,1);
    t0=AE_SRAV32RS(t0, AE_SEXT32X2D16_32(nsa01));
    // masking
    AE_MOVT32X2(t0,c0,xnotpos0);
    AE_S32X2_IP(t0, castxcc(ae_int32x2, pY), 2 * sizeof(ae_int32));
  }
} /* vec_sqrt32x32_fast() */
#endif
