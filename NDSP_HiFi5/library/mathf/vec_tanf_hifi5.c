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

/* DSP Library API */
/*    Code optimized for HiFi5 core */
#include "NatureDSP_Signal_math.h"
/* Common helper macros. */
#include "common_fpu.h"
/* Value of 2/pi, 4/pi, etc. */
#include "inv2pif_tbl.h"
/* tan/cotan approximation polynomial coeffs. */
#include "tanf_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"

#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(void,vec_tanf,( float32_t * restrict z, const float32_t * restrict x, int N ))
#elif HAVE_VFPU
#define sz_f32    (int)sizeof(float32_t)

/*
  NatureDSP Signal Processing Library. Vector matematics
    Tangent
    Code optimized for HiFi5 core
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
static void __tanf( float32_t * restrict y, const float32_t * restrict x, int N);
void vec_tanf( float32_t * restrict y, const float32_t * restrict x, int N)
{
    xtfloatx4 * restrict pX;
    xtfloatx4 * restrict pY;
    int n;

    if (N<=0) return;
    if (N&7)
    {
        float32_t ALIGN(32) xbuf[8],ybuf[8];
        pX=(xtfloatx4 *)xbuf;
        pY=(xtfloatx4 *)ybuf;
        AE_SSX2X2_I(XT_CONST_S(0),XT_CONST_S(0),pX,0*sizeof(xtfloatx4));
        AE_SSX2X2_I(XT_CONST_S(0),XT_CONST_S(0),pX,1*sizeof(xtfloatx4));
        AE_SSX2X2_I(XT_CONST_S(0),XT_CONST_S(0),pY,0*sizeof(xtfloatx4));
        AE_SSX2X2_I(XT_CONST_S(0),XT_CONST_S(0),pY,1*sizeof(xtfloatx4));
        for (n=0; n<(N&7); n++) 
        {
            xtfloat t;
            XT_LSIP(t,castxcc(xtfloat,x  ),sizeof(float32_t));
            XT_SSIP(t,castxcc(xtfloat,pX ),sizeof(float32_t));
        }
        pX=(xtfloatx4 *)xbuf;
        __tanf((float32_t*)pY,(float32_t*)pX,8);
        for (n=0; n<(N&7); n++) 
        {
            xtfloat t;
            XT_LSIP(t,castxcc(xtfloat,pY),sizeof(float32_t));
            XT_SSIP(t,castxcc(xtfloat,y ),sizeof(float32_t));
        }
        N&=~7;
    }
    if (N<=0) return;
    __tanf(y,x,N);
}

#if TANF_ALG==0
static void __tanf( float32_t * restrict y, const float32_t * restrict x, int N)
{

    const xtfloatx4 * restrict X;
    xtfloatx4 * restrict Z;
    const xtfloatx4 * restrict S_rd;
    xtfloatx4 * restrict S_wr;
    const xtfloatx4* restrict S_firound_rd;
    const xtfloatx4   * restrict T;
    /* Current block index; overall number of blocks; number of values in the current block */
    ae_valignx2 X_va, Z_va;

    /* Block size, blkLen <= blkSize */
    const int blkSize = MAX_ALLOCA_SZ/(2*sz_f32);
    /* 2/pi splited into 24-bit chunks*/
    xtfloatx2 pi2fc0, pi2fc1, pi2fc2;

    /* init 24bit chunks of pi/2*/
    pi2fc0 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOV32(0x3fc90fdb));
    pi2fc1 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOV32(0xb33bbd2e));
    pi2fc2 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOV32(0xa6f72ced));

    /* Allocate a fixed-size scratch area on the stack. */
    float32_t ALIGN(32) scr[2*blkSize];
    int n,M;
    NASSERT_ALIGN16(scr);
    NASSERT(N%8==0);
    /*
    * Data are processed in blocks of scratch area size. Further, the algorithm
    * implementation is splitted in order to feed the optimizing compiler with a
    * few loops of managable size.
    */


    Z_va = AE_ZALIGN128();
    for (; N>0; N-=M,x+=M,y+=M)
    {
        M = XT_MIN( N, blkSize );

        /*
        *  Loop splited into 3 parts for better optimization.
        */

        X = (xtfloatx4*)(x);
        S_wr = (xtfloatx4*)scr;
        X = (xtfloatx4*)(x);
        X_va = AE_LA128_PP(X);
        /*
        *   Part I.  Argument reduction.
        */
        __Pragma("loop_count factor=2")
        for (n = 0; n < (M >> 2); n++)
        {
            xtfloatx2 p0, p1;
            xtfloatx2 jf0, jf1;

            AE_LASX2X2_IP(p0, p1, X_va, X);

            ABS_SX2X2(p0, p1, p0, p1);
            MULQ_S(jf0, jf1, p0, p1, inv2pif.f);
            jf0 = FIROUND_SX2(jf0);
            jf1 = FIROUND_SX2(jf1);
            AE_SSX2X2_I(jf0, jf1, S_wr, 4 * sz_f32);

            MSUBQ_S(p0, p1, jf0, jf1, pi2fc0);
            MSUBQ_S(p0, p1, jf0, jf1, pi2fc1);
            MSUBQ_S(p0, p1, jf0, jf1, pi2fc2);

            AE_SSX2X2_IP(p0, p1, S_wr, 8 * sz_f32);
        }

        X    = (xtfloatx4*)(x);
        S_rd = (xtfloatx4*)scr;
        S_wr = (xtfloatx4*)scr;
        T = (const xtfloatx4  *)polytanf_tbl;

        /*
        * Part II. Compute tan via minmax polynomial.
        */
        for ( n=0; n<(M>>2); n++ )
        {
            xtfloatx2 p0,p1;
            xtfloatx2 y0,y1;
            xtfloatx2 p2_0, p2_1;
            xtfloatx2 p3_0, p3_1;
            xtfloatx2 p4_0, p4_1;
            xtfloatx2 cf0_0, cf0_1;
            xtfloatx2 cf1_0, cf1_1;
            xtfloatx2 cf2_0, cf2_1;
            xtfloatx2 cf3_0, cf3_1;
            xtfloatx2 cf4_0, cf4_1;
            xtfloatx2 cf5_0, cf5_1;
            xtfloatx2 cf6_0, cf6_1;

            AE_LSX2X2_IP(p0, p1, S_rd, 8*sz_f32);

            AE_LSX2X2_IP(cf0_0, cf0_1, T, 4 * sz_f32);
            AE_LSX2X2_IP(cf1_0, cf1_1, T, 4 * sz_f32);
            AE_LSX2X2_IP(cf2_0, cf2_1, T, 4 * sz_f32);
            AE_LSX2X2_IP(cf3_0, cf3_1, T, 4 * sz_f32);
            AE_LSX2X2_IP(cf4_0, cf4_1, T, 4 * sz_f32);
            AE_LSX2X2_IP(cf5_0, cf5_1, T, 4 * sz_f32);
            AE_LSX2X2_XP(cf6_0, cf6_1, T, -24 * sz_f32);

            MUL_SX2X2(p2_0, p2_1, p0, p1, p0, p1);
            MUL_SX2X2(p3_0, p3_1, p2_0, p2_1, p0, p1);
            MUL_SX2X2(p4_0, p4_1, p2_0, p2_1, p2_0, p2_1);
            
            MADD_SX2X2(cf2_0, cf2_1, p2_0, p2_1, cf1_0, cf1_1); cf1_0 = cf2_0; cf1_1 = cf2_1;
            MADD_SX2X2(cf4_0, cf4_1, p2_0, p2_1, cf3_0, cf3_1); cf2_0 = cf4_0; cf2_1 = cf4_1;
            MADD_SX2X2(cf6_0, cf6_1, p2_0, p2_1, cf5_0, cf5_1); cf3_0 = cf6_0; cf3_1 = cf6_1;
            y0 = cf0_0; y1 = cf0_1;
            MADD_SX2X2(cf1_0, cf1_1, p4_0, p4_1, y0, y1); y0 = cf1_0; y1 = cf1_1;
            MADD_SX2X2(cf2_0, cf2_1, p4_0, p4_1, y0, y1); y0 = cf2_0; y1 = cf2_1;
            MADD_SX2X2(cf3_0, cf3_1, p4_0, p4_1, y0, y1); y0 = cf3_0; y1 = cf3_1;

            MADD_SX2X2(p0, p1, p3_0, p3_1, y0, y1); 

            AE_SSX2X2_IP(p0, p1, S_wr, 8 * sz_f32);
        }

        X = (xtfloatx4*)(x);
        S_rd = (xtfloatx4*)scr;
        S_firound_rd = (xtfloatx4*)scr + 1;
        Z = (xtfloatx4*)(y);
        X_va = AE_LA128_PP(X);

        /* Part III. Sign adjustment, outliers handling and choosing ctg instead tan. */
        __Pragma("loop_count factor=2")
        for (n = 0; n<(M >> 2); n++)
        {
            xtfloatx2 xn0, t0, p0, p1;
            xtfloatx2 xn1, t1;
            xtbool2 b0_dom, b1_dom;
            /* Input value segment number; input and output signs; integer reprentation of output value */
            ae_int32x2 ji0, sx0;
            ae_int32x2 ji1, sx1;

            /* Cosine and sine approximations; output value */
            ae_int32x2 sz0, sz1, s0, s1;
            xtfloatx2 xm0, xm1, z0, z1;

            /* Load tangent approximation from the scratch. */
            AE_LSX2X2_IP(p0, p1, S_rd, +8 * sz_f32);
            AE_LASX2X2_IP(xn0, xn1, X_va, X);
            /* Re-calculate the pi/2-wide segment number. */
            AE_LSX2X2_IP(t0, t1, S_firound_rd, +8 * sz_f32);
            ABS_SX2X2(xm0, xm1, xn0, xn1);
            ji0 = TRUNC_SX2(t0, 0);
            ji1 = TRUNC_SX2(t1, 0);
            /* Compute the sign adjustment term. */

            sx0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xn0);
            sx1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xn1);
            sz0 = AE_XOR32(sx0, AE_SLLI32(ji0, 31));
            sz1 = AE_XOR32(sx1, AE_SLLI32(ji1, 31));
            sz0 = AE_AND32(sz0, 0x80000000);
            sz1 = AE_AND32(sz1, 0x80000000);

            /* Compute the cotangent for odd-numbered segments. */
            z0 = RECIP_SX2(p0);
            z1 = RECIP_SX2(p1);
            //b_cot = AE_LT32(ji, AE_ZERO32());
            MOVF_SX2(z0, p0, AE_MOVBD1X2(ji0));
            MOVF_SX2(z1, p1, AE_MOVBD1X2(ji1));
            /* Adjust the sign. */
            s0 = AE_MOVINT32X2_FROMXTFLOATX2(z0);
            s1 = AE_MOVINT32X2_FROMXTFLOATX2(z1);
            s0 = AE_XOR32(s0, sz0);
            s1 = AE_XOR32(s1, sz1);
            z0 = AE_MOVXTFLOATX2_FROMINT32X2(s0);
            z1 = AE_MOVXTFLOATX2_FROMINT32X2(s1);

            /* Set result to quiet NaN for an out-of-domain input value. */
            b0_dom = OLE_SX2(xm0, tanf_maxval);
            b1_dom = OLE_SX2(xm1, tanf_maxval);
            MOVF_SX2(z0, qNaNf.f, b0_dom);
            MOVF_SX2(z1, qNaNf.f, b1_dom);
            AE_SASX2X2_IP(z0, z1, Z_va, Z);
        }
        AE_SA128POS_FP( Z_va, Z );
    }
} /* vec_tanf() */
#elif TANF_ALG==1 // original code
static void __tanf(float32_t* restrict y, const float32_t* restrict x, int N)
{
    /*
      float32_t x2,y,yt,yc;
      int sx,n,j,st;
      sx=takesignf(x);
      x=sx?-x:x;
      if(x>tanf_maxval) return qNaN.f;
       argument reduction
         process reduces x by integral multiple of pi/4.
         The output is deduced to the sum of two single precision
         floating point values x+dx.

      n = (int)STDLIB_MATH(ceilf)(x*inv4pif.f);
      j = n&~1;

      {
        float32_t dx, t, y = x, jj = (float32_t)j;
        const union ufloat32uint32 c[6] = {
          { 0x3f4a0000 },
          { 0xbb700000 },
          { 0xb6160000 },
          { 0x32080000 },
          { 0x2e060000 },
          { 0xa9b9ee5a } };
        dx = 0.f;
        y -= c[0].f*jj;
        y -= c[1].f*jj;
        y -= c[2].f*jj;
        t = y; y -= c[3].f*jj; t = (t - y); t -= c[3].f*jj; dx = t;
        t = y; y -= c[4].f*jj; t = (t - y); t -= c[4].f*jj; dx = (dx + t);
        t = y; y -= c[5].f*jj; t = (t - y); t -= c[5].f*jj; dx = (dx + t);
        y = (y + dx);
        x = y;
      }

       compute tan via minmax polynomial
      x2 = x*x;
      yt = polytanf_tbl[0].f;
      yt = yt*x2 + polytanf_tbl[1].f;
      yt = yt*x2 + polytanf_tbl[2].f;
      yt = yt*x2 + polytanf_tbl[3].f;
      yt = yt*x2 + polytanf_tbl[4].f;
      yt = yt*x2 + polytanf_tbl[5].f;
      yt = yt*x2 + polytanf_tbl[6].f;
      yt = yt*x2;

      dx is small enough (<3e-8) and wiil be used to modify
      tangent value computed by polynomial using derivation
      tg(x+dx)=tg(x)+dx/cos(x)^2
      for 2 ULP operation it is enough to suppose 1/cos(x)^2 ~ 1
      for 1 ULP operation, it should be computed accurately

      resulting value is decomposed as follows
      tag(x+dx)~(P*x+dx)+x
      The order of summation is crucial!

      yt = (yt*x) + x;
      yc = 1.f / yt;

      adjust sign
      n = (n >> 1) & 1;
      st = sx ^ n;
       select tan/cotan
      y = n ? yc : yt;
       apply the sign
      y = changesignf(y, st);
      return y;
    */
    const xtfloatx4* restrict X;
    const xtfloatx4* restrict X1;
    xtfloatx4* restrict Z;
    const xtfloatx4* restrict S_rd;
    xtfloatx4* restrict S_wr;
    const xtfloat* restrict T;
    /* Current block index; overall number of blocks; number of values in the current block */
    ae_valignx2 X_va, X1_va, Z_va;

    /* Block size, blkLen <= blkSize */
    const int blkSize = MAX_ALLOCA_SZ / 2 * sz_f32;
    /* pi/2 splitted into 7-bit chunks. */
    static const union ufloat32uint32 ALIGN(32) c[6] = {
      { 0x3fca0000 }, { 0xbbf00000 },
      { 0xb6960000 }, { 0x32880000 },
      { 0x2e860000 }, { 0xaa39ee5a }
    };
    /* Allocate a fixed-size scratch area on the stack. */
    float32_t ALIGN(32) scr[blkSize];
    int n, M;
    NASSERT_ALIGN16(scr);
    NASSERT(N % 8 == 0);
    /*
    * Data are processed in blocks of scratch area size. Further, the algorithm
    * implementation is splitted in order to feed the optimizing compiler with a
    * few loops of managable size.
    */

    for (; N > 0; N -= M, x += M, y += M)
    {
        M = XT_MIN(N, blkSize);
        /*
        * Part I, range reduction. Reference C code:
        *
        *   {
        *     float32_t xn, p, dp, t;
        *     float32_t jf;
        *
        *     // pi/2 splitted into 7-bit chunks.
        *     static const union ufloat32uint32 c[6] = {
        *       { 0x3fca0000 }, { 0xbbf00000 },
        *       { 0xb6960000 }, { 0x32880000 },
        *       { 0x2e860000 }, { 0xaa39ee5a }
        *     };
        *
        *     for ( n=0; n<blkLen; n++ )
        *     {
        *       xn = fabsf( x[blkIx*blkSize+n] );
        *
        *       // Determine the pi/2-wide segment the input value belongs to.
        *       jf = roundf( xn*inv2pif.f );
        *
        *       // Calculate the difference between the segment midpoint and input value.
        *       p = xn;
        *       p -= c[0].f*jf;
        *       p -= c[1].f*jf;
        *       p -= c[2].f*jf;
        *       t = p; p -= c[3].f*jf; t = t - p; t -= c[3].f*jf; dp = t;
        *       t = p; p -= c[4].f*jf; t = t - p; t -= c[4].f*jf; dp += t;
        *       t = p; p -= c[5].f*jf; t = t - p; t -= c[5].f*jf; dp += t;
        *       p += dp;
        *
        *       scr[n] = p;
        *     }
        *   }
        */
        X = (xtfloatx4*)(x);
        S_wr = (xtfloatx4*)scr;
        T = (xtfloat*)c;
        X_va = AE_LA128_PP(X);
        __Pragma("loop_count factor=2");
        for (n = 0; n < (M >> 2); n++)
        {
            /* pi/2 splitted into 7-bit chunks. */
            xtfloatx2 c0, c1, c2, c3, c4, c5;
            /* Scalar auxiliary var.  */
            xtfloat cs;
            xtfloatx2 xn0, xn1;
            xtfloatx2 jf0, jf1;
            xtfloatx2 p0, p1, dp0, dp1, t0, t1, r0, r1;

            AE_LASX2X2_IP(xn0, xn1, X_va, X);

            /* Determine the pi/2-wide segment the input value belongs to. */
            ABS_SX2X2(xn0, xn1, xn0, xn1);
            MULQ_S(jf0, jf1, xn0, xn1, inv2pif.f);
            jf0 = XT_FIROUND_SX2(jf0);          jf1 = XT_FIROUND_SX2(jf1);
            /* Calculate the difference between the segment midpoint and input value. */
            /* For this particular loop, XP address update results in a better schedule if compared with IP. */
            XT_LSIP(cs, T, +1 * sz_f32); c0 = cs;
            XT_LSIP(cs, T, +1 * sz_f32); c1 = cs;
            XT_LSIP(cs, T, +1 * sz_f32); c2 = cs;
            XT_LSIP(cs, T, +1 * sz_f32); c3 = cs;
            XT_LSIP(cs, T, +1 * sz_f32); c4 = cs;
            XT_LSXP(cs, T, -5 * sz_f32); c5 = cs;

            p0 = xn0; p1 = xn1;
            MSUBQ_S(p0, p1, jf0, jf1, c0);
            MSUBQ_S(p0, p1, jf0, jf1, c1);
            MSUBQ_S(p0, p1, jf0, jf1, c2);
            MULQ_S(r0, r1, jf0, jf1, c3); t0 = p0; t1 = p1; SUB_SX2X2(p0, p1, p0, p1, r0, r1); SUB_SX2X2(t0, t1, t0, t1, p0, p1); SUB_SX2X2(t0, t1, t0, t1, r0, r1); dp0 = t0; dp1 = t1;
            MULQ_S(r0, r1, jf0, jf1, c4); t0 = p0; t1 = p1; SUB_SX2X2(p0, p1, p0, p1, r0, r1); SUB_SX2X2(t0, t1, t0, t1, p0, p1); SUB_SX2X2(t0, t1, t0, t1, r0, r1); ADD_SX2X2(dp0, dp1, t0, t1, dp0, dp1);
            MULQ_S(r0, r1, jf0, jf1, c5); t0 = p0; t1 = p1; SUB_SX2X2(p0, p1, p0, p1, r0, r1); SUB_SX2X2(t0, t1, t0, t1, p0, p1); SUB_SX2X2(t0, t1, t0, t1, r0, r1); ADD_SX2X2(dp0, dp1, t0, t1, dp0, dp1);
            ADD_SX2X2(p0, p1, p0, p1, dp0, dp1);

            AE_SSX2X2_IP(p0, p1, S_wr, 4 * sz_f32);
        }
        __Pragma("no_reorder");

        /*
        * Part II, tangent approximation via minmax polynomial. Reference C code:
        *
        *   {
        *     float32_t yt, p, p2;
        *
        *     for ( n=0; n<blkLen; n++ )
        *     {
        *       p = scr[n];
        *       p2 = p*p;
        *
        *       yt = polytanf_tbl[0].f;
        *       yt = polytanf_tbl[1].f + yt*p2;
        *       yt = polytanf_tbl[2].f + yt*p2;
        *       yt = polytanf_tbl[3].f + yt*p2;
        *       yt = polytanf_tbl[4].f + yt*p2;
        *       yt = polytanf_tbl[5].f + yt*p2;
        *       yt = polytanf_tbl[6].f + yt*p2;
        *       yt =                     yt*p2;
        *
        *       scr[n] = yt*p + p;
        *     }
        *   }
        */
        X = (xtfloatx4*)(x);
        X1 = (xtfloatx4*)(x);
        S_rd = (xtfloatx4*)scr;
        S_wr = (xtfloatx4*)scr;
        X = (xtfloatx4*)(x);
        X_va = AE_LA128_PP(X);
        X1_va = AE_LA128_PP(X1);
        Z_va = AE_ZALIGN128();
        T = (xtfloat*)&polytanf_tbl[0];
        __Pragma("loop_count factor=2")
            for (n = 0; n < (M >> 2); n++)
            {
                xtfloatx2 p0, p1;
                xtfloatx2 p2_0, p2_1, p3_0, p3_1, p4_0, p4_1, p8_0, p8_1;
                /* Polynomial coefficients for sine and cosine. */
                xtfloatx2 c0, c1, c2, c3, c4, c5, c6;
                xtfloat cs;
                xtfloatx2 c0_0, c1_0, c0_1, c1_1, c0_2, c1_2;
                /*
                * For this loop, XP address updates result in a better schedule if compared with IP.
                */

                AE_LSX2X2_XP(p0, p1, S_rd, +4 * sz_f32);

                XT_LSIP(cs, T, +1 * sz_f32); c0 = cs;
                XT_LSIP(cs, T, +1 * sz_f32); c1 = cs;
                XT_LSIP(cs, T, +1 * sz_f32); c2 = cs;
                XT_LSIP(cs, T, +1 * sz_f32); c3 = cs;
                XT_LSIP(cs, T, +1 * sz_f32); c4 = cs;
                XT_LSIP(cs, T, +1 * sz_f32); c5 = cs;
                XT_LSXP(cs, T, -6 * sz_f32); c6 = cs;

                MUL_SX2X2(p2_0, p2_1, p0, p1, p0, p1);
                MUL_SX2X2(p3_0, p3_1, p2_0, p2_1, p0, p1);
                MUL_SX2X2(p4_0, p4_1, p2_0, p2_1, p2_0, p2_1);
                MUL_SX2X2(p8_0, p8_1, p4_0, p4_1, p4_0, p4_1);

                c0_0 = c1_0 = c2;   MADDQ_S(c0_0, c1_0, p2_0, p2_1, c1);
                c0_1 = c1_1 = c4;   MADDQ_S(c0_1, c1_1, p2_0, p2_1, c3);
                c0_2 = c1_2 = c6;   MADDQ_S(c0_2, c1_2, p2_0, p2_1, c5);

                MADDQ_S(c0_0, c1_0, p4_0, p4_1, c0);
                MADD_SX2X2(c0_2, c1_2, c0_1, c1_1, p4_0, p4_1);
                MADD_SX2X2(c0_2, c1_2, c0_0, c1_0, p8_0, p8_1);
                MADD_SX2X2(p0, p1, p3_0, p3_1, c0_2, c1_2);
                AE_SSX2X2_IP(p0, p1, S_wr, 4 * sz_f32);
            }
        __Pragma("no_reorder");

        /*
        * Part III, estimation of cotangent and finalization of results. Reference C code:
        *
        *   {
        *     float32_t xn, xm, yt, yc, zn;
        *     int ji, sx, sz;
        *
        *     for ( n=0; n<blkLen; n++ )
        *     {
        *       xn = x[blkIx*blkSize+n];
        *       yt = scr[n];
        *
        *       xm = fabsf(xn);
        *       // Determine the pi/2-wide segment the input value belongs to.
        *       ji = ( (int)roundf( xn*inv2pif.f ) & 1 );
        *       sx = takesignf( xn );
        *       sz = sx ^ ji;
        *       yc = 1.f/yt;
        *       zn = ( ji ? yc : yt );
        *       zn = changesignf( zn, sz );
        *
        *       z[blkIx*blkSize+n] = ( xm<=tanf_maxval ? zn : qNaNf.f );
        *     }
        *   }
        */
        X = (xtfloatx4*)(x);
        X1 = (xtfloatx4*)(x);
        S_rd = (xtfloatx4*)scr;
        S_wr = (xtfloatx4*)scr;
        Z = (xtfloatx4*)(y);
        X_va = AE_LA128_PP(X);
        X1_va = AE_LA128_PP(X1);
        Z_va = AE_ZALIGN128();
        __Pragma("loop_count factor=2")
            for (n = 0; n < (M >> 2); n++)
            {
                xtfloatx2 xn0, t0, p0, p1;
                xtfloatx2 xn1, t1;
                xtbool2 b0_dom, b1_dom;
                /* Input value segment number; input and output signs; integer reprentation of output value */
                ae_int32x2 ji0, sx0;
                ae_int32x2 ji1, sx1;

                /* Cosine and sine approximations; output value */
                xtfloatx2 yc0;
                xtfloatx2 yc1;
                ae_int32x2 sz0, sz1, s0, s1;
                xtfloatx2 xm0, xm1, z0, z1;

                /* Load tangent approximation from the scratch. */
                AE_LSX2X2_XP(p0, p1, S_rd, +4 * sz_f32);
                AE_LASX2X2_IP(xn0, xn1, X_va, X);
                /* Re-calculate the pi/2-wide segment number. */
                ABS_SX2X2(xm0, xm1, xn0, xn1);
                MULQ_S(t0, t1, xm0, xm1, inv2pif.f);
                t0 = XT_FIROUND_SX2(t0);
                t1 = XT_FIROUND_SX2(t1);
                ji0 = XT_TRUNC_SX2(t0, 0);
                ji1 = XT_TRUNC_SX2(t1, 0);
                /* Compute the sign adjustment term. */
                sx0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xn0);
                sx1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xn1);
                sz0 = AE_XOR32(sx0, AE_SLLI32(ji0, 31));
                sz1 = AE_XOR32(sx1, AE_SLLI32(ji1, 31));
                sz0 = AE_SRLI32(sz0, 31);
                sz1 = AE_SRLI32(sz1, 31);
                sz0 = AE_SLLI32(sz0, 31);
                sz1 = AE_SLLI32(sz1, 31);

                /* Compute the cotangent for odd-numbered segments. */
#if 1
                yc0 = XT_RECIP_SX2(p0);
                yc1 = XT_RECIP_SX2(p1);
#else
                {
                    xtfloatx2 t20, t21;
                    yc0 = XT_RECIP0_SX2(p0);
                    yc1 = XT_RECIP0_SX2(p1);
                    CONST_SX2X2(t20, t21, 1);
                    MSUB_SX2X2(t20, t21, p0, p1, yc0, yc1);
                    MADD_SX2X2(yc0, yc1, yc0, yc1, t20, t21);
                    CONST_SX2X2(t20, t21, 1);
                    MSUB_SX2X2(t20, t21, p0, p1, yc0, yc1);
                    MADD_SX2X2(yc0, yc1, yc0, yc1, t20, t21);
                }
#endif
                //b_cot = AE_LT32(ji, AE_ZERO32());
                z0 = yc0; XT_MOVF_SX2(z0, p0, AE_MOVBD1X2(ji0));
                z1 = yc1; XT_MOVF_SX2(z1, p1, AE_MOVBD1X2(ji1));
                /* Adjust the sign. */
                s0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(z0);
                s1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(z1);
                s0 = AE_XOR32(s0, sz0);
                s1 = AE_XOR32(s1, sz1);
                z0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(s0);
                z1 = XT_AE_MOVXTFLOATX2_FROMINT32X2(s1);

                /* Set result to quiet NaN for an out-of-domain input value. */
                b0_dom = XT_OLE_SX2(xm0, tanf_maxval);
                b1_dom = XT_OLE_SX2(xm1, tanf_maxval);
                XT_MOVF_SX2(z0, qNaNf.f, b0_dom);
                XT_MOVF_SX2(z1, qNaNf.f, b1_dom);
                AE_SASX2X2_IP(z0, z1, Z_va, Z);
            }
        AE_SA128POS_FP(Z_va, Z);
    }
} /* vec_tanf() */

#if 0
{
    /*
      float32_t x2,y,yt,yc;
      int sx,n,j,st;
      sx=takesignf(x);
      x=sx?-x:x;
      if(x>tanf_maxval) return qNaN.f;
       argument reduction
         process reduces x by integral multiple of pi/4.
         The output is deduced to the sum of two single precision
         floating point values x+dx.

      n = (int)STDLIB_MATH(ceilf)(x*inv4pif.f);
      j = n&~1;

      {
        float32_t dx, t, y = x, jj = (float32_t)j;
        const union ufloat32uint32 c[6] = {
          { 0x3f4a0000 },
          { 0xbb700000 },
          { 0xb6160000 },
          { 0x32080000 },
          { 0x2e060000 },
          { 0xa9b9ee5a } };
        dx = 0.f;
        y -= c[0].f*jj;
        y -= c[1].f*jj;
        y -= c[2].f*jj;
        t = y; y -= c[3].f*jj; t = (t - y); t -= c[3].f*jj; dx = t;
        t = y; y -= c[4].f*jj; t = (t - y); t -= c[4].f*jj; dx = (dx + t);
        t = y; y -= c[5].f*jj; t = (t - y); t -= c[5].f*jj; dx = (dx + t);
        y = (y + dx);
        x = y;
      }

       compute tan via minmax polynomial
      x2 = x*x;
      yt = polytanf_tbl[0].f;
      yt = yt*x2 + polytanf_tbl[1].f;
      yt = yt*x2 + polytanf_tbl[2].f;
      yt = yt*x2 + polytanf_tbl[3].f;
      yt = yt*x2 + polytanf_tbl[4].f;
      yt = yt*x2 + polytanf_tbl[5].f;
      yt = yt*x2 + polytanf_tbl[6].f;
      yt = yt*x2;

      dx is small enough (<3e-8) and wiil be used to modify
      tangent value computed by polynomial using derivation
      tg(x+dx)=tg(x)+dx/cos(x)^2
      for 2 ULP operation it is enough to suppose 1/cos(x)^2 ~ 1
      for 1 ULP operation, it should be computed accurately

      resulting value is decomposed as follows
      tag(x+dx)~(P*x+dx)+x
      The order of summation is crucial!

      yt = (yt*x) + x;
      yc = 1.f / yt;

      adjust sign
      n = (n >> 1) & 1;
      st = sx ^ n;
       select tan/cotan
      y = n ? yc : yt;
       apply the sign
      y = changesignf(y, st);
      return y;
    */
    const xtfloatx4* restrict X;
    const xtfloatx4* restrict X1;
    xtfloatx4* restrict Z;
    const xtfloatx4* restrict S_rd;
    xtfloatx4* restrict S_wr;
    const xtfloat* restrict T;
    /* Current block index; overall number of blocks; number of values in the current block */
    int blkIx, blkNum, blkLen;
    ae_valignx2 X_va, X1_va, Z_va;

    /* Block size, blkLen <= blkSize */
    const int blkSize = MAX_ALLOCA_SZ / 2 * sz_f32;
    /* Allocate a fixed-size scratch area on the stack. */
    float32_t ALIGN(32) scr[blkSize];
    int n;
    NASSERT_ALIGN16(scr);
    /*
    * Data are processed in blocks of scratch area size. Further, the algorithm
    * implementation is splitted in order to feed the optimizing compiler with a
    * few loops of managable size.
    */

    blkNum = (N + blkSize - 1) / blkSize;

    for (blkIx = 0; blkIx < blkNum; blkIx++)
    {
        blkLen = XT_MIN(N - blkIx * blkSize, blkSize);
        /*
        * Part I, range reduction. Reference C code:
        *
        *   {
        *     float32_t xn, p, dp, t;
        *     float32_t jf;
        *
        *     // pi/2 splitted into 7-bit chunks.
        *     static const union ufloat32uint32 c[6] = {
        *       { 0x3fca0000 }, { 0xbbf00000 },
        *       { 0xb6960000 }, { 0x32880000 },
        *       { 0x2e860000 }, { 0xaa39ee5a }
        *     };
        *
        *     for ( n=0; n<blkLen; n++ )
        *     {
        *       xn = fabsf( x[blkIx*blkSize+n] );
        *
        *       // Determine the pi/2-wide segment the input value belongs to.
        *       jf = roundf( xn*inv2pif.f );
        *
        *       // Calculate the difference between the segment midpoint and input value.
        *       p = xn;
        *       p -= c[0].f*jf;
        *       p -= c[1].f*jf;
        *       p -= c[2].f*jf;
        *       t = p; p -= c[3].f*jf; t = t - p; t -= c[3].f*jf; dp = t;
        *       t = p; p -= c[4].f*jf; t = t - p; t -= c[4].f*jf; dp += t;
        *       t = p; p -= c[5].f*jf; t = t - p; t -= c[5].f*jf; dp += t;
        *       p += dp;
        *
        *       scr[n] = p;
        *     }
        *   }
        */
        {
            /* pi/2 splitted into 7-bit chunks. */
            xtfloatx2 c0, c1, c2, c3, c4, c5;
            /* Scalar auxiliary var.  */
            xtfloat cs;

            /* pi/2 splitted into 7-bit chunks. */
            static const union ufloat32uint32 ALIGN(32) c[6] = {
              { 0x3fca0000 }, { 0xbbf00000 },
              { 0xb6960000 }, { 0x32880000 },
              { 0x2e860000 }, { 0xaa39ee5a }
            };

            X = (xtfloatx4*)(x);
            S_wr = (xtfloatx4*)scr;
            T = (xtfloat*)c;

            X_va = AE_LA128_PP(X);

            __Pragma("loop_count min=1");
            for (n = 0; n < ((blkLen + 3) >> 2); n++)
            {
                xtfloatx2 xn0, xn1;
                xtfloatx2 jf0, jf1;
                xtfloatx2 p0, p1, dp0, dp1, t0, t1, r0, r1;

                AE_LASX2X2_IP(xn0, xn1, X_va, X);

                /* Determine the pi/2-wide segment the input value belongs to. */
                ABS_SX2X2(xn0, xn1, xn0, xn1);
                MULQ_S(jf0, jf1, xn0, xn1, inv2pif.f);
                jf0 = XT_FIROUND_SX2(jf0);          jf1 = XT_FIROUND_SX2(jf1);
                /* Calculate the difference between the segment midpoint and input value. */
                /* For this particular loop, XP address update results in a better schedule if compared with IP. */
                XT_LSXP(cs, T, +1 * sz_f32); c0 = cs;
                XT_LSXP(cs, T, +1 * sz_f32); c1 = cs;
                XT_LSXP(cs, T, +1 * sz_f32); c2 = cs;
                XT_LSXP(cs, T, +1 * sz_f32); c3 = cs;
                XT_LSXP(cs, T, +1 * sz_f32); c4 = cs;
                XT_LSXP(cs, T, -5 * sz_f32); c5 = cs;

                p0 = xn0; p1 = xn1;
                MSUBQ_S(p0, p1, jf0, jf1, c0);
                MSUBQ_S(p0, p1, jf0, jf1, c1);
                MSUBQ_S(p0, p1, jf0, jf1, c2);
                MULQ_S(r0, r1, jf0, jf1, c3); t0 = p0; t1 = p1; SUB_SX2X2(p0, p1, p0, p1, r0, r1); SUB_SX2X2(t0, t1, t0, t1, p0, p1); SUB_SX2X2(t0, t1, t0, t1, r0, r1); dp0 = t0; dp1 = t1;
                MULQ_S(r0, r1, jf0, jf1, c4); t0 = p0; t1 = p1; SUB_SX2X2(p0, p1, p0, p1, r0, r1); SUB_SX2X2(t0, t1, t0, t1, p0, p1); SUB_SX2X2(t0, t1, t0, t1, r0, r1); ADD_SX2X2(dp0, dp1, t0, t1, dp0, dp1);
                MULQ_S(r0, r1, jf0, jf1, c5); t0 = p0; t1 = p1; SUB_SX2X2(p0, p1, p0, p1, r0, r1); SUB_SX2X2(t0, t1, t0, t1, p0, p1); SUB_SX2X2(t0, t1, t0, t1, r0, r1); ADD_SX2X2(dp0, dp1, t0, t1, dp0, dp1);
                ADD_SX2X2(p0, p1, p0, p1, dp0, dp1);

                AE_SSX2X2_IP(p0, p1, S_wr, 4 * sz_f32);
            }
        }
        __Pragma("no_reorder");

        /*
        * Part II, tangent approximation via minmax polynomial. Reference C code:
        *
        *   {
        *     float32_t yt, p, p2;
        *
        *     for ( n=0; n<blkLen; n++ )
        *     {
        *       p = scr[n];
        *       p2 = p*p;
        *
        *       yt = polytanf_tbl[0].f;
        *       yt = polytanf_tbl[1].f + yt*p2;
        *       yt = polytanf_tbl[2].f + yt*p2;
        *       yt = polytanf_tbl[3].f + yt*p2;
        *       yt = polytanf_tbl[4].f + yt*p2;
        *       yt = polytanf_tbl[5].f + yt*p2;
        *       yt = polytanf_tbl[6].f + yt*p2;
        *       yt =                     yt*p2;
        *
        *       scr[n] = yt*p + p;
        *     }
        *   }
        */
        {

            xtfloatx2 p0, p1;
            xtfloatx2 p2_0, p2_1, p3_0, p3_1, p4_0, p4_1, p8_0, p8_1;
            /* Polynomial coefficients for sine and cosine. */
            xtfloatx2 c0, c1, c2, c3, c4, c5, c6;
            xtfloat cs;

            X = (xtfloatx4*)(x);
            X1 = (xtfloatx4*)(x);

            S_rd = (xtfloatx4*)scr;
            S_wr = (xtfloatx4*)scr;
            X = (xtfloatx4*)(x);
            X_va = AE_LA128_PP(X);
            X1_va = AE_LA128_PP(X1);
            Z_va = AE_ZALIGN128();
            T = (xtfloat*)&polytanf_tbl[0];
            for (n = 0; n < ((blkLen + 3) >> 2); n++)
            {
                xtfloatx2 c0_0, c1_0, c0_1, c1_1, c0_2, c1_2;
                /*
                * For this loop, XP address updates result in a better schedule if compared with IP.
                */

                AE_LSX2X2_XP(p0, p1, S_rd, +4 * sz_f32);

                XT_LSIP(cs, T, +1 * sz_f32); c0 = cs;
                XT_LSIP(cs, T, +1 * sz_f32); c1 = cs;
                XT_LSIP(cs, T, +1 * sz_f32); c2 = cs;
                XT_LSIP(cs, T, +1 * sz_f32); c3 = cs;
                XT_LSIP(cs, T, +1 * sz_f32); c4 = cs;
                XT_LSIP(cs, T, +1 * sz_f32); c5 = cs;
                XT_LSXP(cs, T, -6 * sz_f32); c6 = cs;

                MUL_SX2X2(p2_0, p2_1, p0, p1, p0, p1);
                MUL_SX2X2(p3_0, p3_1, p2_0, p2_1, p0, p1);
                MUL_SX2X2(p4_0, p4_1, p2_0, p2_1, p2_0, p2_1);
                MUL_SX2X2(p8_0, p8_1, p4_0, p4_1, p4_0, p4_1);

                c0_0 = c1_0 = c2;   MADDQ_S(c0_0, c1_0, p2_0, p2_1, c1);
                c0_1 = c1_1 = c4;   MADDQ_S(c0_1, c1_1, p2_0, p2_1, c3);
                c0_2 = c1_2 = c6;   MADDQ_S(c0_2, c1_2, p2_0, p2_1, c5);

                MADDQ_S(c0_0, c1_0, p4_0, p4_1, c0);
                MADD_SX2X2(c0_2, c1_2, c0_1, c1_1, p4_0, p4_1);
                MADD_SX2X2(c0_2, c1_2, c0_0, c1_0, p8_0, p8_1);
                MADD_SX2X2(p0, p1, p3_0, p3_1, c0_2, c1_2);
                AE_SSX2X2_IP(p0, p1, S_wr, 4 * sz_f32);

            }
        }
        __Pragma("no_reorder");

        /*
        * Part III, estimation of cotangent and finalization of results. Reference C code:
        *
        *   {
        *     float32_t xn, xm, yt, yc, zn;
        *     int ji, sx, sz;
        *
        *     for ( n=0; n<blkLen; n++ )
        *     {
        *       xn = x[blkIx*blkSize+n];
        *       yt = scr[n];
        *
        *       xm = fabsf(xn);
        *       // Determine the pi/2-wide segment the input value belongs to.
        *       ji = ( (int)roundf( xn*inv2pif.f ) & 1 );
        *       sx = takesignf( xn );
        *       sz = sx ^ ji;
        *       yc = 1.f/yt;
        *       zn = ( ji ? yc : yt );
        *       zn = changesignf( zn, sz );
        *
        *       z[blkIx*blkSize+n] = ( xm<=tanf_maxval ? zn : qNaNf.f );
        *     }
        *   }
        */
        {
            xtfloatx2 xn0, t0, p0, p1;
            xtfloatx2 xn1, t1;
            xtbool2 b0_dom, b1_dom;
            /* Input value segment number; input and output signs; integer reprentation of output value */
            ae_int32x2 ji0, sx0;
            ae_int32x2 ji1, sx1;

            /* Cosine and sine approximations; output value */
            xtfloatx2 yc0;
            xtfloatx2 yc1;

            X = (xtfloatx4*)(x);
            X1 = (xtfloatx4*)(x);

            S_rd = (xtfloatx4*)scr;
            S_wr = (xtfloatx4*)scr;
            Z = (xtfloatx4*)(z);
            X_va = AE_LA128_PP(X);
            X1_va = AE_LA128_PP(X1);
            Z_va = AE_ZALIGN128();
            for (n = 0; n < (blkLen >> 2); n++)
            {
                ae_int32x2 sz0, sz1, s0, s1;
                xtfloatx2 xm0, xm1, z0, z1;

                /* Load tangent approximation from the scratch. */
                AE_LSX2X2_XP(p0, p1, S_rd, +4 * sz_f32);
                AE_LASX2X2_IP(xn0, xn1, X_va, X);
                /* Re-calculate the pi/2-wide segment number. */
                ABS_SX2X2(xm0, xm1, xn0, xn1);
                MULQ_S(t0, t1, xm0, xm1, inv2pif.f);
                t0 = XT_FIROUND_SX2(t0);
                t1 = XT_FIROUND_SX2(t1);
                ji0 = XT_TRUNC_SX2(t0, 0);
                ji1 = XT_TRUNC_SX2(t1, 0);
                //  ji0 = AE_SLLI32(ji0, 31);
                //  ji1 = AE_SLLI32(ji1, 31);
                  /* Compute the sign adjustment term. */
                sx0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xn0);
                sx1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(xn1);
                sz0 = AE_XOR32(sx0, AE_SLLI32(ji0, 31));
                sz1 = AE_XOR32(sx1, AE_SLLI32(ji1, 31));
                sz0 = AE_SRLI32(sz0, 31);
                sz1 = AE_SRLI32(sz1, 31);
                sz0 = AE_SLLI32(sz0, 31);
                sz1 = AE_SLLI32(sz1, 31);
                /* Compute the cotangent for odd-numbered segments. */
                yc0 = XT_RECIP_SX2(p0);
                yc1 = XT_RECIP_SX2(p1);
                //b_cot = AE_LT32(ji, AE_ZERO32());
                z0 = yc0; XT_MOVF_SX2(z0, p0, AE_MOVBD1X2(ji0));
                z1 = yc1; XT_MOVF_SX2(z1, p1, AE_MOVBD1X2(ji1));
                /* Adjust the sign. */
                s0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(z0);
                s1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(z1);
                s0 = AE_XOR32(s0, sz0);
                s1 = AE_XOR32(s1, sz1);
                z0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(s0);
                z1 = XT_AE_MOVXTFLOATX2_FROMINT32X2(s1);

                /* Set result to quiet NaN for an out-of-domain input value. */
                b0_dom = XT_OLE_SX2(xm0, tanf_maxval);
                b1_dom = XT_OLE_SX2(xm1, tanf_maxval);
                XT_MOVF_SX2(z0, qNaNf.f, b0_dom);
                XT_MOVF_SX2(z1, qNaNf.f, b1_dom);
                AE_SASX2X2_IP(z0, z1, Z_va, Z);
            }
            AE_SA128POS_FP(Z_va, Z);
            if ((blkLen & 3) > 1)
            {
                ae_valign X_va, Z_va;
                xtfloatx2 xn, xm, yc, t, zn, yt;
                ae_int32x2 ji, sx, sz, s;
                xtbool2  b_cot, b_dom;
                X_va = AE_LA64_PP(castxcc(xtfloatx2, X));
                Z_va = AE_ZALIGN64();

                XT_LASX2IP(xn, X_va, castxcc(xtfloatx2, X));

                /* Re-calculate the pi/2-wide segment number. */
                xm = XT_ABS_SX2(xn);
                t = XT_MUL_SX2(xm, inv2pif.f);
                t = XT_FIROUND_SX2(t);
                ji = XT_TRUNC_SX2(t, 0);
                ji = AE_SLLI32(ji, 31);

                /* Compute the sign adjustment term. */
                sx = XT_AE_MOVINT32X2_FROMXTFLOATX2(xn);
                sz = AE_XOR32(sx, ji);
                sz = AE_SRLI32(sz, 31);
                sz = AE_SLLI32(sz, 31);

                /* Load tangent approximation from the scratch. */
                XT_LSX2IP(yt, castxcc(xtfloatx2, S_rd), +2 * sz_f32);

                /* Compute the cotangent for odd-numbered segments. */
                yc = XT_RECIP_SX2(yt);
                b_cot = AE_LT32(ji, AE_ZERO32());
                zn = yc; XT_MOVF_SX2(zn, yt, b_cot);

                /* Adjust the sign. */
                s = XT_AE_MOVINT32X2_FROMXTFLOATX2(zn);
                s = AE_XOR32(s, sz);
                zn = XT_AE_MOVXTFLOATX2_FROMINT32X2(s);

                /* Set result to quiet NaN for an out-of-domain input value. */
                b_dom = XT_OLE_SX2(xm, tanf_maxval);
                XT_MOVF_SX2(zn, qNaNf.f, b_dom);

                XT_SASX2IP(zn, Z_va, castxcc(xtfloatx2, Z));
                AE_SA64POS_FP(Z_va, castxcc(xtfloatx2, Z));
            }
            /* Process the last input value, if any. */
            if (blkLen & 1)
            {
                xtfloatx2 xn, xm, yc, t, zn, yt;
                ae_int32x2 ji, sx, sz, s;
                xtbool2  b_cot, b_dom;

                xn = XT_LSI((xtfloat*)X, 0);

                /* Re-calculate the pi/2-wide segment number. */
                xm = XT_ABS_SX2(xn);
                t = XT_MUL_SX2(xm, inv2pif.f);
                t = XT_FIROUND_SX2(t);
                ji = XT_TRUNC_SX2(t, 0);
                ji = AE_SLLI32(ji, 31);

                /* Compute the sign adjustment term. */
                sx = XT_AE_MOVINT32X2_FROMXTFLOATX2(xn);
                sz = AE_XOR32(sx, ji);
                sz = AE_SRLI32(sz, 31);
                sz = AE_SLLI32(sz, 31);

                /* Load tangent approximation from the scratch. */
                yt = XT_LSI((xtfloat*)S_rd, 0);

                /* Compute the cotangent for odd-numbered segments. */
                yc = XT_RECIP_SX2(yt);
                b_cot = AE_LT32(ji, AE_ZERO32());
                zn = yt; XT_MOVT_SX2(zn, yc, b_cot);

                /* Adjust the sign. */
                s = XT_AE_MOVINT32X2_FROMXTFLOATX2(zn);
                s = AE_XOR32(s, sz);
                zn = XT_AE_MOVXTFLOATX2_FROMINT32X2(s);

                /* Set result to quiet NaN for an out-of-domain input value. */
                b_dom = XT_OLE_SX2(xm, tanf_maxval);
                XT_MOVF_SX2(zn, qNaNf.f, b_dom);

                XT_SSI(zn, (xtfloat*)Z, 0);
            }
        }
    }
} /* vec_tanf() */
#endif
#else 
#error wrong TANF_ALG
#endif


#elif HAVE_FPU
#define sz_f32    (int)sizeof(float32_t)

/*===========================================================================
  Vector matematics:
  vec_tan             Tangent    
===========================================================================*/

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

  Scalar versions:
  ----------------
  Return result, Q16.15 or floating point
-------------------------------------------------------------------------*/

#if TANF_ALG==0
void vec_tanf(float32_t* restrict z, const float32_t* restrict x, int N)
{
    /*
    * const union ufloat32uint32 pi2fc[3] = {
    *             { 0x3fc90fdb }, { 0xb33bbd2e },
    *             { 0xa6f72ced }
    * };
    * float32_t dx, x2, y, yt, yc;
    * int sx, st, j;
    * sx = takesignf(x);
    * x = sx ? -x : x;
    *
    *    argument reduction
    *    process reduces x by integral multiple of pi/2.
    *    To arrain 1 ULP accuracy, the output is reduced to the sum of
    *    two single precision floating point values x+dx.
    *    For max error of 2 ULP in [-9099,9099] the correction term is unnecessary.
    * 
    * j = roundf(x * inv2pif.f);
    * x = fmaf(-j, pi2fc[0].f, x);
    * x = fmaf(-j, pi2fc[1].f, x);
    * x = fmaf(-j, pi2fc[2].f, x);
    * (void)dx;
    * 
    *   compute tan via minmax polynomial  
    * x2 = x * x;
    * yt = polytanf_tbl[0].f;
    * yt = yt * x2 + polytanf_tbl[1].f;
    * yt = yt * x2 + polytanf_tbl[2].f;
    * yt = yt * x2 + polytanf_tbl[3].f;
    * yt = yt * x2 + polytanf_tbl[4].f;
    * yt = yt * x2 + polytanf_tbl[5].f;
    * yt = yt * x2 + polytanf_tbl[6].f;
    * yt = yt * x2;
    * 
    *    dx is small enough (<3e-8) and wiil be used to modify
    *    tangent value computed by polynomial using derivation
    *    tan(x+dx)=tan(x)+dx/cos(x)^2
    *    Resulting value is decomposed as follows
    *    tag(x+dx)~(P*x+dx/cs2)+x
    *    The order of summation is crucial!
    * 
    * yt = yt * x + x;
    * yc = 1.f / yt;
    * 
    *   adjust sign 
    * st = sx ^ (j & 1);
    *   select tan/cotan 
    * y = (j & 1) ? yc : yt;
    *   apply the sign 
    * y = changesignf(y, st);
    * return y;
    */

    const xtfloat* restrict  X;
    int32_t* restrict  Z;
    const xtfloat* restrict  S_rd;
    xtfloat* restrict  S_wr;
    const xtfloat* restrict  T;

    /* Current block index; overall number of blocks; number of values in the current block */
    int blkIx, blkNum, blkLen;
    /* Block size, blkLen <= blkSize */
    const int blkSize = MAX_ALLOCA_SZ / sz_f32;
    /* Allocate a fixed-size scratch area on the stack. */
    float32_t ALIGN(32) scr[blkSize];

    int n;

    if (N <= 0) return;

    NASSERT_ALIGN8(scr);

    /*
    * Data are processed in blocks of scratch area size. Further, the algorithm
    * implementation is splitted in order to feed the optimizing compiler with a
    * few loops of managable size.
    */

    blkNum = (N + blkSize - 1) / blkSize;

    for (blkIx = 0; blkIx < blkNum; blkIx++)
    {
        blkLen = XT_MIN(N - blkIx * blkSize, blkSize);
        /*
        * Part I, range reduction. Reference C code:
        *
        *   {
        *     float32_t xn, p, dp, t;
        *     int ji;
        *     float32_t jf;
        *
        *     const union ufloat32uint32 pi2fc[3] = {
        *                 { 0x3fc90fdb }, { 0xb33bbd2e },
        *                 { 0xa6f72ced }
        *     };
        *
        *     for ( n=0; n<blkLen; n++ )
        *     {
        *       xn = fabsf( x[blkIx*blkSize+n] );
        *
        *       // Determine the pi/2-wide segment the input value belongs to.
        *       jf = roundf(x * inv2pif.f);
        *       jf = (float32_t)( ji & ~1 );
        *
        *       // Calculate the difference between the segment midpoint and input value.
        *       p = xn;
        *       p -= c[0].f*jf;
        *       p -= c[1].f*jf;
        *       p -= c[2].f*jf;
        *
        *       scr[n] = p;
        *     }
        *   }
        */
        {
            /* Input value; reducted input value; correction term. */
            xtfloat xn, p;
            /* Input value segment number. */
            xtfloat jf;
            /* 2/pi 72-bit splited into 24-bit chunks*/
            xtfloat c0, c1, c2;

            const union ufloat32uint32 c[3] = {
                        { 0x3fc90fdb }, { 0xb33bbd2e },
                        { 0xa6f72ced }
            };

            
            X = (xtfloat*)((uintptr_t)x + blkIx * blkSize * sz_f32);
            S_wr = (xtfloat*)scr;
            T = (xtfloat*)c;

            for (n = 0; n < (blkLen); n++)
            {
                XT_LSIP(xn, X, sz_f32);
                /*
                * Determine the pi/2-wide segment the input value belongs to.
                */
                xn = XT_ABS_S(xn);
                
                jf = XT_MUL_S(xn, inv2pif.f);

                jf = XT_ROUND_S(jf, 0);

                c0 = XT_LSI(T, 0 * sz_f32);
                c1 = XT_LSI(T, 1 * sz_f32);
                c2 = XT_LSI(T, 2 * sz_f32);

                p = xn;
                XT_MSUB_S(p, jf, c0);
                XT_MSUB_S(p, jf, c1);
                XT_MSUB_S(p, jf, c2);

                XT_SSIP(p, S_wr, sz_f32);
            }
        }
        __Pragma("no_reorder");

        /*
        * Part II, tangent approximation via minmax polynomial. Reference C code:
        *
        *   {
        *     float32_t yt, p, p2;
        *
        *     for ( n=0; n<blkLen; n++ )
        *     {
        *       p = scr[n];
        *       p2 = p*p;
        *
        *       yt = polytanf_tbl[0].f;
        *       yt = polytanf_tbl[1].f + yt*p2;
        *       yt = polytanf_tbl[2].f + yt*p2;
        *       yt = polytanf_tbl[3].f + yt*p2;
        *       yt = polytanf_tbl[4].f + yt*p2;
        *       yt = polytanf_tbl[5].f + yt*p2;
        *       yt = polytanf_tbl[6].f + yt*p2;
        *       yt =                     yt*p2;
        *
        *       scr[n] = yt*p + p;
        *     }
        *   }
        */
        {
            /* Reducted input value and its 2nd power; tangent approximation. */
            xtfloat p, p2, yt;
            /* Polynomial coeffs. */
            xtfloat cf0, cf1, cf2, cf3, cf4, cf5, cf6;

            S_rd = (xtfloat*)scr;
            S_wr = (xtfloat*)scr;
            T = (xtfloat*)polytanf_tbl;
            __Pragma("loop_count min=1");
            for (n = 0; n < (blkLen); n++)
            {
                XT_LSIP(p, S_rd, 0 * sz_f32);

                /* Reload polynomial coeffs. */
                cf0 = XT_LSI(T, 0 * sz_f32);
                cf1 = XT_LSI(T, 4 * sz_f32);
                cf2 = XT_LSI(T, 8 * sz_f32);
                cf3 = XT_LSI(T, 12 * sz_f32);
                cf4 = XT_LSI(T, 16 * sz_f32);
                cf5 = XT_LSI(T, 20 * sz_f32);
                cf6 = XT_LSI(T, 24 * sz_f32);

                p2 = XT_MUL_S(p, p);

                yt = cf0;
                XT_MADD_S(cf1, p2, yt); yt = cf1;
                XT_MADD_S(cf2, p2, yt); yt = cf2;
                XT_MADD_S(cf3, p2, yt); yt = cf3;
                XT_MADD_S(cf4, p2, yt); yt = cf4;
                XT_MADD_S(cf5, p2, yt); yt = cf5;
                XT_MADD_S(cf6, p2, yt); yt = cf6;
                yt = XT_MUL_S(p2, yt);
                XT_LSIP(p, S_rd, 1 * sz_f32);
                XT_MADD_S(p, p, yt); yt = p;

                XT_SSIP(yt, S_wr, sz_f32);
            }
        }
        __Pragma("no_reorder");

        /*
        * Part III, estimation of cotangent and finalization of results. Reference C code:
        *
        *   {
        *     float32_t xn, xm, yt, yc, zn;
        *     int ji, sx, sz;
        *
        *     for ( n=0; n<blkLen; n++ )
        *     {
        *       xn = x[blkIx*blkSize+n];
        *       yt = scr[n];
        *
        *       xm = fabsf(xn);
        *       ji = (int)round( xm*inv2pif.f );
        *       ji = ji & 1;
        *       sx = takesignf( xn );
        *       sz = sx ^ ji;
        *       yc = 1.f/yt;
        *       zn = ( ji ? yc : yt );
        *       zn = changesignf( zn, sz );
        *
        *       z[blkIx*blkSize+n] = ( xm<=tanf_maxval ? zn : qNaNf.f );
        *     }
        *   }
        */
        {
            /* Input value and its magnitude; resulting value. */
            xtfloat xn, xm;
            /* Tangent and cotangent */
            xtfloat yt, yc;
            /* Input value segment number. */
            int32_t ji;
            /* Auxiliary floating-point var. */
            xtfloat t;
            /* Input sign; result sign; auxiliary var.  */
            int32_t sx, sz, s;
            /* Cotangent/tangent selection flag; input data validation flag. */
            xtbool b_cot, b_dom;

            X = (xtfloat*)((uintptr_t)x + blkIx * blkSize * sz_f32);
            Z = (int32_t*)((uintptr_t)z + blkIx * blkSize * sz_f32);

            S_rd = (xtfloat*)scr;
            __Pragma("loop_count min=1");
            for (n = 0; n < blkLen; n++)
            {
                XT_LSIP(xn, X, sz_f32);
                /* Load tangent approximation from the scratch. */
                XT_LSIP(yt, S_rd, sz_f32);
                /* Re-calculate the pi/2-wide segment number. */
                xm = XT_ABS_S(xn);
                b_dom = XT_OLE_S(xm, tanf_maxval);
                t = XT_MUL_S(xm, inv2pif.f);
                t = XT_ROUND_S(t,0);
                ji = (int)XT_TRUNC_S(t, 0);
                ji = ji << 31;
                /* Compute the sign adjustment term. */
                sx = XT_RFR(xn);
                sz = sx ^ ji;
                sz = sz & 0x80000000;
                /* Compute the cotangent for odd-numbered segments. */
                yc = XT_RECIP_S(yt);
                b_cot = AE_INT64_LT(AE_MOVINT64_FROMINT32(ji), AE_ZERO64());

                XT_MOVT_S(yt, yc, (b_cot));
                /* Adjust the sign. */
                s = XT_RFR(yt);
                s = XT_XOR(s, sz);

                /* Set result to quiet NaN for an out-of-domain input value. */
                {
                    unsigned int _t = s;
                    XT_MOVF(_t, qNaNf.u, b_dom); s = _t;
                }
                *Z++ = s;
            }
        }
    }
} /* vec_tanf() */
#else
// Taken from Fusion
void vec_tanf( float32_t * restrict z, const float32_t * restrict x, int N )
{
  /*
  * float32_t x2,y,yt,yc;
  * int sx,n,j,st;
  * sx=takesignf(x);
  * x=sx?-x:x;
  * if(x>tanf_maxval) return qNaN.f;
  * argument reduction
  * process reduces x by integral multiple of pi/4.
  * The output is deduced to the sum of two single precision
  * floating point values x+dx.
  * 
  * n = (int)STDLIB_MATH(ceilf)(x*inv4pif.f);
  * j = n&~1;
  * 
  * {
  * float32_t dx, t, y = x, jj = (float32_t)j;
  * const union ufloat32uint32 c[6] = {
  * { 0x3f4a0000 },
  * { 0xbb700000 },
  * { 0xb6160000 },
  * { 0x32080000 },
  * { 0x2e060000 },
  * { 0xa9b9ee5a } };
  * dx = 0.f;
  * y -= c[0].f*jj;
  * y -= c[1].f*jj;
  * y -= c[2].f*jj;
  * t = y; y -= c[3].f*jj; t = (t - y); t -= c[3].f*jj; dx = t;
  * t = y; y -= c[4].f*jj; t = (t - y); t -= c[4].f*jj; dx = (dx + t);
  * t = y; y -= c[5].f*jj; t = (t - y); t -= c[5].f*jj; dx = (dx + t);
  * y = (y + dx);
  * x = y;
  * }
  * 
  * compute tan via minmax polynomial
  * x2 = x*x;
  * yt = polytanf_tbl[0].f;
  * yt = yt*x2 + polytanf_tbl[1].f;
  * yt = yt*x2 + polytanf_tbl[2].f;
  * yt = yt*x2 + polytanf_tbl[3].f;
  * yt = yt*x2 + polytanf_tbl[4].f;
  * yt = yt*x2 + polytanf_tbl[5].f;
  * yt = yt*x2 + polytanf_tbl[6].f;
  * yt = yt*x2;
  * 
  * dx is small enough (<3e-8) and wiil be used to modify
  * tangent value computed by polynomial using derivation
  * tg(x+dx)=tg(x)+dx/cos(x)^2
  * for 2 ULP operation it is enough to suppose 1/cos(x)^2 ~ 1
  * for 1 ULP operation, it should be computed accurately
  * 
  * resulting value is decomposed as follows
  * tag(x+dx)~(P*x+dx)+x
  * The order of summation is crucial!
  * 
  * yt = (yt*x) + x;
  * yc = 1.f / yt;
  * 
  * adjust sign
  * n = (n >> 1) & 1;
  * st = sx ^ n;
  * select tan/cotan
  * y = n ? yc : yt;
  * apply the sign
  * y = changesignf(y, st);
  * return y;
  */

  const xtfloat * restrict  X;
        int32_t * restrict  Z;
  const xtfloat * restrict  S_rd;
        xtfloat * restrict  S_wr;
  const xtfloat * restrict  T;
  const int     * restrict  pT;

  /* Current block index; overall number of blocks; number of values in the current block */
  int blkIx, blkNum, blkLen;
  /* Block size, blkLen <= blkSize */
  const int blkSize = MAX_ALLOCA_SZ / sz_f32;
  /* Allocate a fixed-size scratch area on the stack. */
  float32_t ALIGN(32) scr[blkSize];

  int n;

  if (N <= 0) return;

  NASSERT_ALIGN8(scr);

  /*
  * Data are processed in blocks of scratch area size. Further, the algorithm
  * implementation is splitted in order to feed the optimizing compiler with a
  * few loops of managable size.
  */

  blkNum = (N + blkSize - 1) / blkSize;

  for (blkIx = 0; blkIx<blkNum; blkIx++)
  {
    blkLen = XT_MIN(N - blkIx*blkSize, blkSize);
    /*
    * Part I, range reduction. Reference C code:
    *
    *   {
    *     float32_t xn, p, dp, t;
    *     int ji;
    *     float32_t jf;
    *
    *     static const union ufloat32uint32 c[6] = {
    *       { 0x3f4a0000 }, { 0xbb700000 },
    *       { 0xb6160000 }, { 0x32080000 },
    *       { 0x2e060000 }, { 0xa9b9ee5a }
    *     };
    *
    *     for ( n=0; n<blkLen; n++ )
    *     {
    *       xn = fabsf( x[blkIx*blkSize+n] );
    *
    *       // Determine the pi/2-wide segment the input value belongs to.
    *       ji = (int)ceilf( xn*inv4pif.f );
    *       jf = (float32_t)( ji & ~1 );
    *
    *       // Calculate the difference between the segment midpoint and input value.
    *       p = xn;
    *       p -= c[0].f*jf;
    *       p -= c[1].f*jf;
    *       p -= c[2].f*jf;
    *       t = p; p -= c[3].f*jf; t = t - p; t -= c[3].f*jf; dp = t;
    *       t = p; p -= c[4].f*jf; t = t - p; t -= c[4].f*jf; dp += t;
    *       t = p; p -= c[5].f*jf; t = t - p; t -= c[5].f*jf; dp += t;
    *       p += dp;
    *
    *       scr[n] = p;
    *     }
    *   }
    */
    {
      /* Input value; reducted input value; correction term. */
      xtfloat xn, p, dp;
      /* Auxiliary floating-point vars. */
      xtfloat t, r;
      /* Input value segment number. */
      ae_int32 ji, i0;
      xtfloat jf;
      /* pi/4 splitted into 7-bit chunks. */
      xtfloat c0, c1, c2, c3, c4, c5;

      static const uint32_t ALIGN(32) c[6] = {
        0x3f4a0000 ,  0xbb700000 ,
        0xb6160000 ,  0x32080000 ,
        0x2e060000 ,  0xa9b9ee5a 
      };
      /* 4/pi, ~1 */
      static const uint32_t TAB[2] = { 0x3fa2f983, 0xFFFFFFFE };
      X = (xtfloat*)((uintptr_t)x + blkIx*blkSize*sz_f32);
      S_wr = (xtfloat*)scr;
      T = (xtfloat  *)c;

      pT = (int *)TAB;
      for (n = 0; n<(blkLen); n++)
      {
        XT_LSIP(xn, X, sz_f32);
        /*
        * Determine the pi/2-wide segment the input value belongs to.
        */
        xn = XT_ABS_S(xn);
        XT_LSIP(c0, castxcc(xtfloat, pT), 0*sz_f32);
		t = XT_MUL_S(xn, c0);
        t = XT_FLOAT_S(XT_CEIL_S(t, 0), 0);
        ji = XT_TRUNC_S(t, 0);
        i0 = XT_L32I(pT, 1*sz_f32);
        ji = XT_AND(ji, i0);
        jf = XT_FLOAT_S(ji, 0);

        /*
        * Calculate the difference between the segment midpoint and input value.
        */

        c0 = XT_LSI( T, 0 * sz_f32);
        c1 = XT_LSI( T, 1 * sz_f32);
        c2 = XT_LSI( T, 2 * sz_f32);
        c3 = XT_LSI( T, 3 * sz_f32);
        c4 = XT_LSI( T, 4 * sz_f32);
        c5 = XT_LSI( T, 5 * sz_f32);

        p = xn;
        XT_MSUB_S(p, jf, c0);
        XT_MSUB_S(p, jf, c1);
        XT_MSUB_S(p, jf, c2);

        r = XT_MUL_S(jf, c3); t = p; p = XT_SUB_S(p, r); t = XT_SUB_S(t, p); t = XT_SUB_S(t, r); dp = t;
        r = XT_MUL_S(jf, c4); t = p; p = XT_SUB_S(p, r); t = XT_SUB_S(t, p); t = XT_SUB_S(t, r); dp = XT_ADD_S(t, dp);
        r = XT_MUL_S(jf, c5); t = p; p = XT_SUB_S(p, r); t = XT_SUB_S(t, p); t = XT_SUB_S(t, r); dp = XT_ADD_S(t, dp);

        p = XT_ADD_S(p, dp);

        XT_SSIP(p, S_wr, sz_f32);
      }
    }
    __Pragma("no_reorder");

    /*
    * Part II, tangent approximation via minmax polynomial. Reference C code:
    *
    *   {
    *     float32_t yt, p, p2;
    *
    *     for ( n=0; n<blkLen; n++ )
    *     {
    *       p = scr[n];
    *       p2 = p*p;
    *
    *       yt = polytanf_tbl[0].f;
    *       yt = polytanf_tbl[1].f + yt*p2;
    *       yt = polytanf_tbl[2].f + yt*p2;
    *       yt = polytanf_tbl[3].f + yt*p2;
    *       yt = polytanf_tbl[4].f + yt*p2;
    *       yt = polytanf_tbl[5].f + yt*p2;
    *       yt = polytanf_tbl[6].f + yt*p2;
    *       yt =                     yt*p2;
    *
    *       scr[n] = yt*p + p;
    *     }
    *   }
    */
    {
      /* Reducted input value and its 2nd power; tangent approximation. */
      xtfloat p, p2, yt;
      /* Polynomial coeffs. */
      xtfloat cf0, cf1, cf2, cf3, cf4, cf5, cf6;

      S_rd = (xtfloat*)scr;
      S_wr = (xtfloat*)scr;
      T = (xtfloat  *)polytanf_tbl;
      /* Pre-load polynomial coeffs. */
      cf0 = XT_LSI(T, 0 * sz_f32);
      cf1 = XT_LSI(T, 1 * sz_f32);
      __Pragma("loop_count min=1");
      for (n = 0; n<(blkLen); n++)
      {
		XT_LSIP(p, S_rd, 0*sz_f32);

        /* Reload polynomial coeffs. */
        cf0 = XT_LSI(T, 0 * sz_f32);
        cf1 = XT_LSI(T, 1 * sz_f32);
        cf2 = XT_LSI(T, 2 * sz_f32);
        cf3 = XT_LSI(T, 3 * sz_f32);
        cf4 = XT_LSI(T, 4 * sz_f32);
        cf5 = XT_LSI(T, 5 * sz_f32);
        cf6 = XT_LSI(T, 6 * sz_f32);

        p2 = XT_MUL_S(p, p);

        yt = cf0;
        XT_MADD_S(cf1, p2, yt); yt = cf1;
        XT_MADD_S(cf2, p2, yt); yt = cf2;
        XT_MADD_S(cf3, p2, yt); yt = cf3;
        XT_MADD_S(cf4, p2, yt); yt = cf4;
        XT_MADD_S(cf5, p2, yt); yt = cf5;
        XT_MADD_S(cf6, p2, yt); yt = cf6;
        yt = XT_MUL_S(p2, yt);
        XT_LSIP(p, S_rd, 1*sz_f32);
        XT_MADD_S(p, p, yt); yt = p;
       
        XT_SSIP(yt, S_wr, sz_f32);
      }
    }
    __Pragma("no_reorder");

    /*
    * Part III, estimation of cotangent and finalization of results. Reference C code:
    *
    *   {
    *     float32_t xn, xm, yt, yc, zn;
    *     int ji, sx, sz;
    *
    *     for ( n=0; n<blkLen; n++ )
    *     {
    *       xn = x[blkIx*blkSize+n];
    *       yt = scr[n];
    *
    *       xm = fabsf(xn);
    *       ji = (int)ceilf( xm*inv4pif.f );
    *       ji = (ji>>1) & 1;
    *       sx = takesignf( xn );
    *       sz = sx ^ ji;
    *       yc = 1.f/yt;
    *       zn = ( ji ? yc : yt );
    *       zn = changesignf( zn, sz );
    *
    *       z[blkIx*blkSize+n] = ( xm<=tanf_maxval ? zn : qNaNf.f );
    *     }
    *   }
    */
    {
      /* Input value and its magnitude; resulting value. */
      xtfloat xn, xm;
      /* Tangent and cotangent */
      xtfloat yt, yc;
      /* Input value segment number. */
      int32_t ji;
      /* Auxiliary floating-point var. */
      xtfloat t;
      /* Input sign; result sign; auxiliary var.  */
      int32_t sx, sz, s;
      /* Cotangent/tangent selection flag; input data validation flag. */
      xtbool b_cot, b_dom;

      X = (xtfloat*)((uintptr_t)x + blkIx*blkSize*sz_f32);
      Z = (int32_t*)((uintptr_t)z + blkIx*blkSize*sz_f32);

      S_rd = (xtfloat*)scr;
      __Pragma("loop_count min=1");
      for (n = 0; n<blkLen; n++)
      {
        XT_LSIP(xn, X, sz_f32);
        /* Load tangent approximation from the scratch. */
        XT_LSIP(yt, S_rd, sz_f32);
        /* Re-calculate the pi/2-wide segment number. */
        xm = XT_ABS_S(xn);
        b_dom = XT_OLE_S(xm, tanf_maxval);
        t = XT_MUL_S(xm, inv4pif.f);
        t = XT_FLOAT_S(XT_CEIL_S(t, 0), 0);
        ji = (int)XT_TRUNC_S(t, 0);
        ji = ji << 30;
        /* Compute the sign adjustment term. */
        sx = XT_RFR(xn);
        sz = sx ^ ji;
        sz = sz & 0x80000000;
        /* Compute the cotangent for odd-numbered segments. */
        yc = XT_RECIP_S(yt);
        b_cot = AE_INT64_LT(AE_MOVINT64_FROMINT32(ji), AE_ZERO64());

        XT_MOVT_S(yt, yc, (b_cot));
        /* Adjust the sign. */
        s = XT_RFR(yt);
        s = XT_XOR(s, sz);

		/* Set result to quiet NaN for an out-of-domain input value. */
        {
          unsigned int _t = s;
          XT_MOVF(_t, qNaNf.u, b_dom); s = _t;
        }
        *Z++=s;
      }
    }
  }
} /* vec_tanf() */
#endif

#endif
