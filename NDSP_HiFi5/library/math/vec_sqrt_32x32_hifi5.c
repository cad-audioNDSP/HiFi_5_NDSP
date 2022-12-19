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

/*
  NatureDSP Signal Processing Library. Vector Mathematics
   Square Root
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/

#define USE_SCRATCH 0

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
void vec_sqrt32x32 (int32_t * restrict y, const int32_t * restrict x, int N)
#if 0
{
        ae_int32x2     * restrict py  = (      ae_int32x2   *)y;
  const ae_int32x2     * restrict px  = (const ae_int32x2   *)x;
  const ae_int32x2     * restrict tbl = (const ae_int32x2   *)sqrt_table;
  int         n, nsaa, nsab, idxa, idxb;
  ae_int32x2  vcw, vxw, zero, min_int32, vxa, vxb, vxc,  t0, t1, c1, c4;
  ae_f32x2    vxf, vzf;
  xtbool2     inf_ab, is_zero;
  ae_valign    x_align,y_align;
         
  if (N<=0) return;
  zero = AE_ZERO32();
  c1 = AE_MOVI( ~1 );
  c4 = AE_MOVDA32X2( 0xFFFFFFF8,0xFFFFFFF8 );
  min_int32 = AE_MOVDA32X2(MIN_INT32, MIN_INT32);

  /*algorithm*/
  //   f(x) = sqrt(x)  
  //   f(x) = f(xo) + f'(x0)*(x-x0)
  //   x-x0=dx
  //   x in range 0..1
  //   x0 in range 0.25..1
  x_align = AE_LA64_PP(px);
  y_align = AE_ZALIGN64();

  __Pragma("concurrent")
  for (n=0; n<(N>>1); n++)
  {
    /*Load two input numbers */
    AE_LA32X2_IP( vxb, x_align, px);
    is_zero = AE_EQ32(vxb, zero);
    /* Normalize x to x0*/
    nsab = AE_NSAZ32_L(vxb);
    nsaa = AE_NSAZ32_L(AE_SEL32_HH(vxb,vxb));

    t0 = AE_MOVDA32X2( nsaa, nsab );
    t1 = AE_AND32( t0, c1 );

    AE_MOVSARD7(t1);
    vxw = AE_SLAS32S( vxb);
    t0 = AE_SRAI32( t0, 1 );
    inf_ab = AE_LT32(vxw, zero);

    /*get table indexes*/
    vxa = AE_SRAI32(vxw, 21);
    vxa = AE_AND32( vxa, c4 );
    idxa= AE_MOVAD32_H(vxa); //
    idxb= AE_MOVAD32_L(vxa); //
    
    vxa = AE_L32X2_X(tbl, idxa); 
    vxb = AE_L32X2_X(tbl, idxb);

    vcw = AE_SEL32_LL(vxa, vxb);//f'(x0)
    vxc = AE_SEL32_HH(vxa, vxb);// f(x0)
    
    /*Calculate dx for first iteration*/  
    vzf = (vxc);//
    vxf = (vxw);
    AE_MULSFP32X2RAS(vxf, vxc, vxc);// 
    /*Calculate y in first iteration*/
    AE_MULAFP32X2RAS(vzf, vcw, vxf);// vzf = y,Q31
    vxf = (vxw);//
    /*Calculate dx for second iteration*/ 
    AE_MULSFP32X2RAS(vxf, vzf, vzf);//
    /*Calculate y in second iteration*/
    AE_MULAFP32X2RAS(vzf, vcw, vxf);// y, Q31
    {
        ae_int32x2 vxd;
        nsab=AE_MOVAD32_L(t0);
        nsaa=AE_MOVAD32_H(t0);
        vxd = AE_SRAA32(vzf, nsab);
        vxc = AE_SRAA32(vzf, nsaa);
        vxc = AE_SEL32_HL(vxc, vxd);
    }

    AE_MOVT32X2(vxc, min_int32, inf_ab);
    AE_MOVT32X2(vxc, zero, is_zero);
    AE_SA32X2_IP(vxc, y_align, py);
  }
  AE_SA64POS_FP(y_align, py);
  if (N&1)
  {
    /*Load two input numbers */
    vxb=AE_L32_I( (const ae_int32*)px, 0);
    is_zero = AE_EQ32(vxb, zero);
    /* Normalize x to x0*/
    nsab = AE_NSAZ32_L(vxb);
    nsaa = AE_NSAZ32_L(AE_SEL32_HH(vxb,vxb));

    t0 = AE_MOVDA32X2( nsaa, nsab );
    t1 = AE_AND32( t0, c1 );

    AE_MOVSARD7(t1);
    vxw = AE_SLAS32S( vxb);
    t0 = AE_SRAI32( t0, 1 );
    inf_ab = AE_LT32(vxw, zero);

    /*get table indexes*/
    vxa = AE_SRAI32(vxw, 21);
    vxa = AE_AND32( vxa, c4 );
    idxa= AE_MOVAD32_H(vxa); //
    idxb= AE_MOVAD32_L(vxa); //
    
    vxa = AE_L32X2_X(tbl, idxa); 
    vxb = AE_L32X2_X(tbl, idxb);

    vcw = AE_SEL32_LL(vxa, vxb);//f'(x0)
    vxc = AE_SEL32_HH(vxa, vxb);// f(x0)
    
    /*Calculate dx for first iteration*/  
    vzf = (vxc);//
    vxf = (vxw);
    AE_MULSFP32X2RAS(vxf, vxc, vxc);// 
    /*Calculate y in first iteration*/
    AE_MULAFP32X2RAS(vzf, vcw, vxf);// vzf = y,Q31
    vxf = (vxw);//
    /*Calculate dx for second iteration*/ 
    AE_MULSFP32X2RAS(vxf, vzf, vzf);//
    /*Calculate y in second iteration*/
    AE_MULAFP32X2RAS(vzf, vcw, vxf);// y, Q31
    AE_MOVSARD7(t0);
    vxc=AE_SRAS32(vzf);
    {
    ae_int32x2 vxd;
    nsab=AE_MOVAD32_L(t0);
    nsaa=AE_MOVAD32_H(t0);
    vxd = AE_SRAA32(vzf, nsab);
    vxc = AE_SRAA32(vzf, nsaa);
    vxc = AE_SEL32_HL(vxc, vxd);
    }

    AE_MOVT32X2(vxc, min_int32, inf_ab);
    AE_MOVT32X2(vxc, zero, is_zero);
    AE_S32_L_I(vxc, (ae_int32*)py, 0);
  }
}
#elif USE_SCRATCH==1
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
    int n,M;
    int32_t ALIGN(128) scratch[MAX_ALLOCA_SZ/sizeof(int32_t)];
          ae_int32x4* restrict pScr;
    const ae_int32x4* restrict pX;
          ae_int32x4* restrict pY;
    ae_valignx2 aX,aY;
    if (N<=0) return;
    for (M=0;N>=8; N-=M,x+=M,y+=M)
    {
        M=(XT_MIN(MAX_ALLOCA_SZ/sizeof(int32_t),N) & (~7));
        pX  =(const ae_int32x4*)x;
        pScr=(      ae_int32x4*)scratch;
        pY  =(      ae_int32x4*)y;
        aX=AE_LA128_PP(pX);
        for (n=0; n<(M>>3); n++)
        {
            ae_int32x2 x0,x1,x2,x3;
            ae_int16x4 nsa01,nsa23;
            ae_int16x4 y0,y1,z0,z1,w0,w1,d0,d1;
            ae_int32x2 t0,t1,t2,t3;

            AE_LA32X2X2_IP(x0,x1,aX,pX);
            AE_LA32X2X2_IP(x2,x3,aX,pX);

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

            AE_S32X2X2_IP(t0,t1,pScr,sizeof(ae_int32x4));
            AE_S32X2X2_IP(t2,t3,pScr,sizeof(ae_int32x4));
        }
        pX  =(const ae_int32x4*)x;
        pScr=(      ae_int32x4*)scratch;
        pY  =(      ae_int32x4*)y;
        aX=AE_LA128_PP(pX);
        aY=AE_ZALIGN128();
        for (n=0; n<(M>>3); n++)
        {
            ae_int32x2 x0,x1,x2,x3;
            ae_int32x2 t0,t1,t2,t3;
            xtbool2 xnotpos0,xnotpos1,xnotpos2,xnotpos3;
            ae_int32x2 c0,c1,c2,c3;

            AE_LA32X2X2_IP(x0,x1,aX,pX);
            AE_LA32X2X2_IP(x2,x3,aX,pX);
            AE_L32X2X2_IP(t0,t1,pScr,sizeof(ae_int32x4));
            AE_L32X2X2_IP(t2,t3,pScr,sizeof(ae_int32x4));

            xnotpos0=AE_LE32(x0,0);
            xnotpos1=AE_LE32(x1,0);
            xnotpos2=AE_LE32(x2,0);
            xnotpos3=AE_LE32(x3,0);
            c0=MIN_INT32; AE_MOVT32X2(c0,0,AE_EQ32(x0,0));
            c1=MIN_INT32; AE_MOVT32X2(c1,0,AE_EQ32(x1,0));
            c2=MIN_INT32; AE_MOVT32X2(c2,0,AE_EQ32(x2,0));
            c3=MIN_INT32; AE_MOVT32X2(c3,0,AE_EQ32(x3,0));

            AE_MOVT32X2(t0,c0,xnotpos0);
            AE_MOVT32X2(t1,c1,xnotpos1);
            AE_MOVT32X2(t2,c2,xnotpos2);
            AE_MOVT32X2(t3,c3,xnotpos3);
            AE_SA32X2X2_IP(t0,t1,aY,pY);
            AE_SA32X2X2_IP(t2,t3,aY,pY);
        }
        AE_SA128POS_FP(aY,pY);
    }
    if (N>0)
    {
        ae_int32x2 x0,x1,x2,x3;
        ae_int16x4 nsa01,nsa23;
        ae_int16x4 y0,y1,z0,z1,w0,w1,d0,d1;
        ae_int32x2 t0,t1,t2,t3;
            xtbool2 xnotpos0,xnotpos1,xnotpos2,xnotpos3;
            ae_int32x2 c0,c1,c2,c3;
        pScr=(      ae_int32x4*)scratch;
        pX  =(const ae_int32x4*)x;
        pY  =(      ae_int32x4*)y;
        AE_S32X2X2_I(0,0,pScr,0);
        AE_S32X2X2_I(0,0,pScr,sizeof(ae_int32x4));
        __Pragma("no_unroll")
        for(n=0; n<N; n++) 
        {
            ae_int32x2 t;
            AE_L32_IP(t,castxcc(ae_int32,pX),sizeof(int32_t));
            AE_S32_L_IP(t,castxcc(ae_int32,pScr),sizeof(int32_t));
        }
        pScr=(ae_int32x4*)scratch;

        AE_L32X2X2_I (x0,x1,pScr,0*sizeof(ae_int32x4));
        AE_L32X2X2_I (x2,x3,pScr,1*sizeof(ae_int32x4));
        // set masks for outbound conditions
        xnotpos0=AE_LE32(x0,0);
        xnotpos1=AE_LE32(x1,0);
        xnotpos2=AE_LE32(x2,0);
        xnotpos3=AE_LE32(x3,0);
        c0=MIN_INT32; AE_MOVT32X2(c0,0,AE_EQ32(x0,0));
        c1=MIN_INT32; AE_MOVT32X2(c1,0,AE_EQ32(x1,0));
        c2=MIN_INT32; AE_MOVT32X2(c2,0,AE_EQ32(x2,0));
        c3=MIN_INT32; AE_MOVT32X2(c3,0,AE_EQ32(x3,0));
        // normalization
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
        // masking
        AE_MOVT32X2(t0,c0,xnotpos0);
        AE_MOVT32X2(t1,c1,xnotpos1);
        AE_MOVT32X2(t2,c2,xnotpos2);
        AE_MOVT32X2(t3,c3,xnotpos3);
        AE_S32X2X2_I(t0,t1,pScr,0*sizeof(ae_int32x4));
        AE_S32X2X2_I(t2,t3,pScr,1*sizeof(ae_int32x4));
        __Pragma("no_unroll")
        for(n=0; n<N; n++) 
        {
            ae_int32x2 t;
            AE_L32_IP(t,castxcc(ae_int32,pScr),sizeof(int32_t));
            AE_S32_L_IP(t,castxcc(ae_int32,pY),sizeof(int32_t));
        }
    }
}
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
    ae_valignx2 aX,aY;
    if (N<=0) return;
    pX  =(const ae_int32x4*)x;
    pY  =(      ae_int32x4*)y;
    aX=AE_LA128_PP(pX);
    aY=AE_ZALIGN128();
    for (n=0; n<(N>>2); n++)
    {
        ae_int32x2 x0,x1;
        ae_int16x4 nsa01;
        ae_int16x4 y0,z0,w0,d0;
        ae_int32x2 t0,t1;
        xtbool2 xnotpos0,xnotpos1;
        ae_int32x2 c0,c1;

        AE_LA32X2X2_IP(x0,x1,aX,pX);
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
        AE_SA32X2X2_IP(t0,t1,aY,pY);
    }
    AE_SA128POS_FP(aY,pY);
    x+=(N&~3);
    y+=(N&~3);
    N&=3;
    if (N>0)
    {
        int32_t ALIGN(32) scratch[4];
        ae_int32x4 *pScr;
        ae_int32x2 x0,x1;
        ae_int16x4 nsa01;
        ae_int16x4 y0,z0,w0,d0;
        ae_int32x2 t0,t1;
        xtbool2 xnotpos0,xnotpos1;
        ae_int32x2 c0,c1;
        pScr=(      ae_int32x4*)scratch;
        pX  =(const ae_int32x4*)x;
        pY  =(      ae_int32x4*)y;
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
        // set masks for outbound conditions
        xnotpos0=AE_LE32(x0,0);
        xnotpos1=AE_LE32(x1,0);
        c0=MIN_INT32; AE_MOVT32X2(c0,0,AE_EQ32(x0,0));
        c1=MIN_INT32; AE_MOVT32X2(c1,0,AE_EQ32(x1,0));
        // normalization
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
        // masking
        AE_MOVT32X2(t0,c0,xnotpos0);
        AE_MOVT32X2(t1,c1,xnotpos1);
        AE_S32X2X2_I(t0,t1,pScr,0*sizeof(ae_int32x4));
        __Pragma("no_unroll")
        for(n=0; n<N; n++) 
        {
            ae_int32x2 t;
            AE_L32_IP(t,castxcc(ae_int32,pScr),sizeof(int32_t));
            AE_S32_L_IP(t,castxcc(ae_int32,pY),sizeof(int32_t));
        }
    }
} /* vec_sqrt32x32() */
#endif
