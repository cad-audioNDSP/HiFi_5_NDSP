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
#include "NatureDSP_Signal_math.h"
/* Common helper macros. */
#include "common_fpu.h"
/* Tables */
#include "pif_tbl.h"
#include "atanf_tbl.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"

#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(void,vec_atanf,( float32_t * restrict z, const float32_t * restrict x, int N ))
#elif HAVE_VFPU
#define sz_f32    (int)sizeof(float32_t)

/*
  NatureDSP Signal Processing Library. Vector Mathematics
  Arctangent
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/


/*-------------------------------------------------------------------------
  Arctangent 
  Functions calculate arctangent of number. Fixed point functions 
  scale output to pi so it is always in range -0x20000000 : 0x20000000 
  which corresponds to the real phases +pi/4 and represent input and output 
  in Q31
  NOTE:
  1.  Scalar floating point function is compatible with standard ANSI C
      routines and sets errno and exception flags accordingly

  Accuracy:
  24 bit version: 74000 (3.4e-5) 
  32 bit version: 42    (2.0e-8)
  floating point: 2 ULP

  Precision: 
  32x32  32-bit inputs, 32-bit output.
  f      floating point

  Input:
  x[N]   input data, Q31 or floating point
  N      length of vectors
  Output:
  z[N]   result, Q31 or floating point

  Restriction:
  x,z should not overlap

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  x,z - aligned on 16-byte boundary
  N   - multiple of 2

  Scalar versions:
  ----------------
  return result, Q31 or floating point
-------------------------------------------------------------------------*/
static void __atanf (    float32_t * restrict y, const float32_t * restrict x, int N );

void vec_atanf (    float32_t * restrict y, 
              const float32_t * restrict x, 
              int N )
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
        __atanf((float32_t*)pY,(float32_t*)pX,8);
        for (n=0; n<(N&7); n++) 
        {
            xtfloat t;
            XT_LSIP(t,castxcc(xtfloat,pY),sizeof(float32_t));
            XT_SSIP(t,castxcc(xtfloat,y ),sizeof(float32_t));
        }
        N&=~7;
    }
    if (N<=0) return;
    __atanf(y,x,N);
} /* vec_atanf() */

static void __atanf (    float32_t * restrict y, const float32_t * restrict x, int N )
{
  /*
    float32_t y;
    int sx,big;
    const union ufloat32uint32* p;
     range reduction 
    sx = x<0;
    x = sx ? -x : x;
    big = x>1.0f;
    if (big) x = 1.0f / x;
    p = (x<0.5f) ? atanftbl1 : atanftbl2;
     approximate atan(x)/x-1 
    y = p[0].f;
    y = x*y + p[1].f;
    y = x*y + p[2].f;
    y = x*y + p[3].f;
    y = x*y + p[4].f;
    y = x*y + p[5].f;
    y = x*y + p[6].f;
    y = x*y + p[7].f;
     convert result to true atan(x) 
    y = x*y + x;

    if (big) y = pi2f.f - y;
    y = sx ? -y : y; apply original sign 
    return y;
  */

  const xtfloatx2 *          X;
  const xtfloatx4 *          X0;
        xtfloatx4 * restrict Z0;

  const xtfloatx4 *          S0_rd;
        xtfloatx2 * restrict S_wr;
  const xtfloat   *          T1;
  const xtfloat   *          T2;

  ae_valign X_va, Z_va;
  ae_valignx2 X0_va, Z0_va;
  /* Block size, blkLen <= blkSize */
  const int blkSize = MAX_ALLOCA_SZ/2*sz_f32;
  /* Allocate a fixed-size scratch area on the stack. */
  float32_t ALIGN(32) scr[blkSize];
  int n,M;

  NASSERT(N%8==0);
  NASSERT_ALIGN16( scr );

  /*
   * Data are processed in blocks of scratch area size. Further, the algorithm
   * implementation is splitted in order to feed the optimizing compiler with a
   * few loops of managable size.
   */
    for (; N>0; N-=M, x+=M,y+=M)
    {
        M = XT_MIN( N, blkSize );

        /*
         * Part I, range reduction. Reference C code:
         *
         *   {
         *     float32_t x0, y0;
         *   
         *     for ( n=0; n<blkLen; n++ )
         *     {
         *       x0 = fabsf( x[blkIx*blkSize+n] );
         *       y0 = ( x0>1.f ? 1.f/x0 : x0 );
         *       scr[n] = y0;
         *     }
         *   }
         */

        X0   = (xtfloatx4*)(x );
        S_wr = (xtfloatx2*)scr;
        X0_va = AE_LA128_PP(X0);
        __Pragma( "loop_count factor=2" )
        for ( n=0; n<(N>>2); n++ )
        {
            /* Input value; reducted value. */
            xtfloatx2 x0, y0;
            xtfloatx2 x1, y1;
            /* Is greater than one; is a +/-infinity  */
            xtbool2 b_gt10, b_inf0;
            xtbool2 b_gt11, b_inf1;
            AE_LASX2X2_IP( x0,x1, X0_va, X0 );
            ABS_SX2X2( x0,x1,x0,x1 );
            b_inf0 = XT_OEQ_SX2( plusInff.f, x0 );
            b_inf1 = XT_OEQ_SX2( plusInff.f, x1 );
            b_gt10 = XT_OLT_SX2(XT_CONST_S(1), x0 );
            b_gt11 = XT_OLT_SX2(XT_CONST_S(1), x1 );
            /* y <- 1.f/x */
            y0 = XT_RECIP_SX2( x0 );
            y1 = XT_RECIP_SX2( x1 );
            /* Fast reciprocal refinement produces NaN for an infinity on input! */
            XT_MOVT_SX2( y0, XT_CONST_S(0), b_inf0 );
            XT_MOVT_SX2( y1, XT_CONST_S(0), b_inf1 );
            /* Select reciprocal for x>1.f */
            XT_MOVF_SX2( y0, x0, b_gt10);
            XT_MOVF_SX2( y1, x1, b_gt11);
            AE_SSX2X2_IP(y0,y1,castxcc(xtfloatx4,S_wr),sizeof(xtfloatx4));
        }
        __Pragma( "no_reorder" );

    /*
     * Part II, polynomial approximation. Reference C code:
     *
     *   {
     *     const union ufloat32uint32 * ptbl;
     *     float32_t x0, y0, z0;
     *   
     *     for ( n=0; n<blkLen; n++ )
     *     {
     *       x0 = x[blkIx*blkSize+n];
     *       y0 = scr[n];
     *       
     *       ptbl = ( y0<0.5f ? atanftbl1 : atanftbl2 );
     *   
     *       // Approximate atan(x)/x-1
     *       z0 = ptbl[0].f;
     *       z0 = ptbl[1].f + y0*z0;
     *       z0 = ptbl[2].f + y0*z0;
     *       z0 = ptbl[3].f + y0*z0;
     *       z0 = ptbl[4].f + y0*z0;
     *       z0 = ptbl[5].f + y0*z0;
     *       z0 = ptbl[6].f + y0*z0;
     *       z0 = ptbl[7].f + y0*z0;
     *       z0 =        y0 + y0*z0;
     *   
     *       if ( fabsf(x0)>1.f ) z0 = pi2f.f - z0;
     *   
     *       // Restore the input sign.
     *       z0 = setsignf( z0, takesignf(x0) );
     *   
     *       z[blkIx*blkSize+n] = z0;
     *     }
     *   }
     */
        X    = (xtfloatx2*)x;
        X0   = (xtfloatx4*)x;
        Z0   = (xtfloatx4*)y;

        S0_rd = (xtfloatx4*)scr;
        T1   = (xtfloat  *)atanftbl1;
        T2   = (xtfloat  *)atanftbl2;

        X_va = XT_LASX2PP( X ); 
        Z_va = AE_ZALIGN64();
        X0_va = AE_LA128_PP(X0);
        Z0_va = AE_ZALIGN128();
        /* Process input data packed in quadruples. */
        for ( n=0; n<(M>>2); n++ )
        {
            /* Input value; reducted input value; output value. */
            xtfloatx2 x0, x1, y0, y1, z0, z1;
            /* Polynomial coeffs for 0.f<=y<0.5f (#1) and 0.5f<=y<=1.f (#2). */
            xtfloatx2 cf1_0, cf1_1, cf1_2, cf1_3, cf1_4, cf1_5, cf1_6, cf1_7;
            xtfloatx2 cf2_0, cf2_1, cf2_2, cf2_3, cf2_4, cf2_5, cf2_6, cf2_7;
            /* Temporary scalar for auto-incrementing coeff loads. */
            xtfloat cf;
            /* Input sign; integer representation of output value. */
            ae_int32x2  z0_i, z1_i, sx0, sx1;
            /* Is greater than one; is less than 0.5f */
            xtbool2 b0_gt1, b0_lt05;
            xtbool2 b1_gt1, b1_lt05;
            xtfloatx2 f0, f1;//, n0, n1;
            xtfloatx2 af1_0,af1_1,af1_2,af1_3;
            xtfloatx2 af2_0,af2_1,af2_2,af2_3;
            xtfloatx2 bf1_0,bf1_1,bf1_2,bf1_3;
            xtfloatx2 bf2_0,bf2_1,bf2_2,bf2_3;
            xtfloatx2 a0,a1,b0,b1;
            xtfloatx2 t0,t1;           
            /* Load coeff sets for both polynomial variants */

            XT_LSIP( cf, T1, +1*sz_f32 ); cf1_0 = cf;
            XT_LSIP( cf, T1, +1*sz_f32 ); cf1_1 = cf;
            XT_LSIP( cf, T1, +1*sz_f32 ); cf1_2 = cf;
            XT_LSIP( cf, T1, +1*sz_f32 ); cf1_3 = cf;
            XT_LSIP( cf, T1, +1*sz_f32 ); cf1_4 = cf;
            XT_LSIP( cf, T1, +1*sz_f32 ); cf1_5 = cf;
            XT_LSIP( cf, T1, +1*sz_f32 ); cf1_6 = cf;
            XT_LSXP( cf, T1, -7*sz_f32 ); cf1_7 = cf;

            XT_LSIP( cf, T2, +1*sz_f32 ); cf2_0 = cf;
            XT_LSIP( cf, T2, +1*sz_f32 ); cf2_1 = cf;
            XT_LSIP( cf, T2, +1*sz_f32 ); cf2_2 = cf;
            XT_LSIP( cf, T2, +1*sz_f32 ); cf2_3 = cf;
            XT_LSIP( cf, T2, +1*sz_f32 ); cf2_4 = cf;
            XT_LSIP( cf, T2, +1*sz_f32 ); cf2_5 = cf;
            XT_LSIP( cf, T2, +1*sz_f32 ); cf2_6 = cf;
            XT_LSXP( cf, T2, -7*sz_f32 ); cf2_7 = cf;
            /*-----------------------------------------------------------------------------*
            *                     Process input values #0 and #1; #2 and #3                          */

            /* Load reduced input values from the scratch area. For this loop, XP address
            * update results in a better schedule if compared with IP. */

            AE_LSX2X2_IP(y0, y1, S0_rd, 4 * sz_f32);
            b0_lt05 = XT_OLT_SX2(y0, XT_CONST_S(3));
            b1_lt05 = XT_OLT_SX2(y1, XT_CONST_S(3));

            /*
            * Compute the approximation to z(y) = atan(y)/y-1. Here we use a combination
            * of Estrin's rule and Horner's method of polynomial evaluation to shorten the
            * dependency path at the cost of additional multiplication.
            */

            /* compute 2 different polynomials and select one of them */
            af1_0=cf1_1; af2_0=cf2_1; MADDQ_S(af1_0, af2_0, cf1_0, cf2_0, y0); 
            bf1_0=cf1_1; bf2_0=cf2_1; MADDQ_S(bf1_0, bf2_0, cf1_0, cf2_0, y1); 
            af1_1=cf1_3; af2_1=cf2_3; MADDQ_S(af1_1, af2_1, cf1_2, cf2_2, y0); 
            bf1_1=cf1_3; bf2_1=cf2_3; MADDQ_S(bf1_1, bf2_1, cf1_2, cf2_2, y1); 
            af1_2=cf1_5; af2_2=cf2_5; MADDQ_S(af1_2, af2_2, cf1_4, cf2_4, y0); 
            bf1_2=cf1_5; bf2_2=cf2_5; MADDQ_S(bf1_2, bf2_2, cf1_4, cf2_4, y1); 
            af1_3=cf1_7; af2_3=cf2_7; MADDQ_S(af1_3, af2_3, cf1_6, cf2_6, y0); 
            bf1_3=cf1_7; bf2_3=cf2_7; MADDQ_S(bf1_3, bf2_3, cf1_6, cf2_6, y1); 

            MUL_SX2X2(x0, x1, y0, y1, y0, y1); 
            a0 = af1_1; a1 = af2_1; MADDQ_S(a0, a1, af1_0, af2_0,x0); 
            b0 = bf1_1; b1 = bf2_1; MADDQ_S(b0, b1, bf1_0, bf2_0,x1); 
            t0 = af1_2; t1 = af2_2; MADDQ_S(t0, t1, a0, a1, x0); a0 = t0; a1 = t1;
            t0 = bf1_2; t1 = bf2_2; MADDQ_S(t0, t1, b0, b1, x1); b0 = t0; b1 = t1;
            t0 = af1_3; t1 = af2_3; MADDQ_S(t0, t1, a0, a1, x0); a0 = t0; a1 = t1;
            t0 = bf1_3; t1 = bf2_3; MADDQ_S(t0, t1, b0, b1, x1); b0 = t0; b1 = t1;
            t0=y0; t1=y0;           MADDQ_S(t0, t1, a0, a1, y0); a0 = t0; a1 = t1;
            t0=y1; t1=y1;           MADDQ_S(t0, t1, b0, b1, y1); b0 = t0; b1 = t1;

            f0 = a0; XT_MOVF_SX2(f0, a1, b0_lt05);
            f1 = b0; XT_MOVF_SX2(f1, b1, b1_lt05);

            /* Load original input values (not the reduced ones). */
            AE_LASX2X2_IP(x0, x1, X0_va, X0);
         
            /* Extract the sign bit and take absolute value. */
            sx0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(x0);
            sx1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(x1);
            sx0 = AE_SRLI32(sx0, 31);
            sx0 = AE_SLLI32(sx0, 31);
            sx1 = AE_SRLI32(sx1, 31);
            sx1 = AE_SLLI32(sx1, 31);
            ABS_SX2X2(x0, x1, x0, x1);
            /* Account for the range reduction. */
            b0_gt1 = XT_OLT_SX2((xtfloatx2)1.0f, x0);
            b1_gt1 = XT_OLT_SX2((xtfloatx2)1.0f, x1);
            SUB_SX2X2(z0, z1, pi2f.f, pi2f.f, f0, f1);
            XT_MOVF_SX2(z0, f0, b0_gt1);
            XT_MOVF_SX2(z1, f1, b1_gt1);

            /* Propagate the input sign. */
            z0_i = XT_AE_MOVINT32X2_FROMXTFLOATX2(z0);
            z1_i = XT_AE_MOVINT32X2_FROMXTFLOATX2(z1);
            z0_i = AE_OR32(z0_i, sx0);
            z1_i = AE_OR32(z1_i, sx1);
            z0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(z0_i);
            z1 = XT_AE_MOVXTFLOATX2_FROMINT32X2(z1_i);
            /* Save resulting values. */
            AE_SASX2X2_IP(z0, z1, Z0_va, Z0);
        }
        AE_SA128POS_FP( Z0_va, Z0);
    }
} /* vec_atanf() */
#elif HAVE_FPU
#define sz_f32    (int)sizeof(float32_t)

/*===========================================================================
  Vector matematics:
  vec_tan          Arctangent        
===========================================================================*/

/*-------------------------------------------------------------------------
  Arctangent 
  Functions calculate arctangent of number. Fixed point functions 
  scale output to pi so it is always in range -0x20000000 : 0x20000000 
  which corresponds to the real phases +pi/4 and represent input and output 
  in Q31
  NOTE:
  1.  Scalar floating point function is compatible with standard ANSI C
      routines and sets errno and exception flags accordingly

  Accuracy:
  32 bit version: 42    (2.0e-8)
  floating point: 2 ULP

  Precision: 
  32x32  32-bit inputs, 32-bit output.
  f      floating point
 
  Input:
  x[N]   input data, Q31 or floating point
  N      length of vectors
  Output:
  z[N]   result, Q31 or floating point

  Restriction:
  x,z should not overlap

  Scalar versions:
  ----------------
  return result, Q31 or floating point
-------------------------------------------------------------------------*/

// Taken from Fusion
void vec_atanf( float32_t * restrict z, const float32_t * restrict x, int N )
{
  /*
  * float32_t y;
  * int sx,big;
  * const union ufloat32uint32* p;
  * range reduction
  * sx = x<0;
  * x = sx ? -x : x;
  * big = x>1.0f;
  * if (big) x = 1.0f / x;
  * p = (x<0.5f) ? atanftbl1 : atanftbl2;
  * approximate atan(x)/x-1
  * y = p[0].f;
  * y = x*y + p[1].f;
  * y = x*y + p[2].f;
  * y = x*y + p[3].f;
  * y = x*y + p[4].f;
  * y = x*y + p[5].f;
  * y = x*y + p[6].f;
  * y = x*y + p[7].f;
  * convert result to true atan(x)
  * y = x*y + x;
  * 
  * if (big) y = pi2f.f - y;
  * y = sx ? -y : y; apply original sign
  * return y;
  */
  const xtfloat *          X;
        int32_t * restrict Z;
  const xtfloat *          S_rd;
        xtfloat * restrict S_wr;
  const xtfloat *          POLY_TBL1;
  const xtfloat *          POLY_TBL2;
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
  POLY_TBL1 = (xtfloat*)atanftbl1;
  POLY_TBL2 = (xtfloat*)atanftbl2;
  for (blkIx = 0; blkIx<blkNum; blkIx++)
  {
    blkLen = XT_MIN(N - blkIx*blkSize, blkSize);
    /*
    * Part I, range reduction. Reference C code:
    *
    *   {
    *     float32_t x0, y0;
    *
    *     for ( n=0; n<blkLen; n++ )
    *     {
    *       x0 = fabsf( x[blkIx*blkSize+n] );
    *       y0 = ( x0>1.f ? 1.f/x0 : x0 );
    *       scr[n] = y0;
    *     }
    *   }
    */

    {
      /* Input value; reducted value. */
      xtfloat x0, y0;
      /* Is greater than one; is a +/-infinity  */
      xtbool b_gt1, b_inf;

      X = (xtfloat*)((uintptr_t)x + blkIx*blkSize*sz_f32);
      S_wr = (xtfloat*)scr;
      __Pragma("loop_count min=1");
      for (n = 0; n<(blkLen ); n++)
      {
       XT_LSIP(x0, X, sz_f32);
  
       x0 = XT_ABS_S(x0);
       b_inf = XT_OEQ_S(plusInff.f, x0);
       b_gt1 = XT_OLT_S(XT_CONST_S(1), x0);
  
       /* y <- 1.f/x */
       y0 = XT_RECIP_S(x0);
  
       /* Fast reciprocal refinement produces NaN for an infinity on input! */
       XT_MOVT_S(y0, XT_CONST_S(0), b_inf);
       /* Select reciprocal for x>1.f */
       XT_MOVF_S(y0, x0, b_gt1);
  
       XT_SSIP(y0, S_wr, sz_f32);
      }
    }
    __Pragma("no_reorder");

    /*
    * Part II, polynomial approximation. Reference C code:
    *
    *   {
    *     const union ufloat32uint32 * ptbl;
    *     float32_t x0, y0, z0;
    *
    *     for ( n=0; n<blkLen; n++ )
    *     {
    *       x0 = x[blkIx*blkSize+n];
    *       y0 = scr[n];
    *
    *       ptbl = ( y0<0.5f ? atanftbl1 : atanftbl2 );
    *
    *       // Approximate atan(x)/x-1
    *       z0 = ptbl[0].f;
    *       z0 = ptbl[1].f + y0*z0;
    *       z0 = ptbl[2].f + y0*z0;
    *       z0 = ptbl[3].f + y0*z0;
    *       z0 = ptbl[4].f + y0*z0;
    *       z0 = ptbl[5].f + y0*z0;
    *       z0 = ptbl[6].f + y0*z0;
    *       z0 = ptbl[7].f + y0*z0;
    *       z0 =        y0 + y0*z0;
    *
    *       if ( fabsf(x0)>1.f ) z0 = pi2f.f - z0;
    *
    *       // Restore the input sign.
    *       z0 = setsignf( z0, takesignf(x0) );
    *
    *       z[blkIx*blkSize+n] = z0;
    *     }
    *   }
    */
    {
      /* Input value; reducted input value; output value. */
      xtfloat x0, y0, z0, z1;
      /* Polynomial coeffs for 0.f<=y<0.5f (#1) and 0.5f<=y<=1.f (#2). */
      xtfloat cf1_0, cf1_1, cf1_2, cf1_3, cf1_4, cf1_5, cf1_6, cf1_7;
      xtfloat cf2_0, cf2_1, cf2_2, cf2_3, cf2_4, cf2_5, cf2_6, cf2_7;
      /* Selected polynomial coeffs. */
      xtfloat cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7;
      /* Input sign; integer representation of output value. */
      int32_t sx, z0_i;
      /* Is greater than one; is less than 0.5f */
      xtbool b_gt1, b_lt05;

      X = (xtfloat*)((uintptr_t)x + blkIx*blkSize*sz_f32);
      Z = (int32_t*)((uintptr_t)z + blkIx*blkSize*sz_f32);
      S_rd = (xtfloat*)scr;

      /* Pre-load first input value */
      XT_LSIP(x0, X, sz_f32);
      for (n = 0; n<blkLen ; n++)
      {
        /*Reload polynomial coeff set #2. */

        cf2_0 = XT_LSI(POLY_TBL2, 0 * sz_f32);
        cf2_1 = XT_LSI(POLY_TBL2, 1 * sz_f32);
        cf2_2 = XT_LSI(POLY_TBL2, 2 * sz_f32);
        cf2_3 = XT_LSI(POLY_TBL2, 3 * sz_f32);
        cf2_4 = XT_LSI(POLY_TBL2, 4 * sz_f32);
        cf2_5 = XT_LSI(POLY_TBL2, 5 * sz_f32);
        cf2_6 = XT_LSI(POLY_TBL2, 6 * sz_f32);
        cf2_7 = XT_LSI(POLY_TBL2, 7 * sz_f32);

        /* Extract the sign bit and take absolute value. */
        sx = XT_RFR(x0);
        sx = sx & 0x80000000;
        x0 = XT_ABS_S(x0);
      
        XT_LSIP(y0, S_rd, sz_f32);
      
        b_lt05 = XT_OLT_S(y0, XT_CONST_S(3));
      
        /* Reload coeff set #1 on each iteration. */

        cf1_0 = XT_LSI(POLY_TBL1, 0 * sz_f32);
        cf1_1 = XT_LSI(POLY_TBL1, 1 * sz_f32);
        cf1_2 = XT_LSI(POLY_TBL1, 2 * sz_f32);
        cf1_3 = XT_LSI(POLY_TBL1, 3 * sz_f32);
        cf1_4 = XT_LSI(POLY_TBL1, 4 * sz_f32);
        cf1_5 = XT_LSI(POLY_TBL1, 5 * sz_f32);
        cf1_6 = XT_LSI(POLY_TBL1, 6 * sz_f32);
        cf1_7 = XT_LSI(POLY_TBL1, 7 * sz_f32);
        /* Select coeffs from sets #1, #2 by reducted value's magnitude. */
        {
          xtfloat p0, p1;
          p0 = cf1_0;
          p1 = cf2_0;
          XT_MOVF_S(p0, p1, b_lt05); cf0 = p0;
          p0 = cf1_1;
          p1 = cf2_1;
          XT_MOVF_S(p0, p1, b_lt05); cf1 = p0;
          p0 = cf1_2;
          p1 = cf2_2;
          XT_MOVF_S(p0, p1, b_lt05); cf2 = p0;
          p0 = cf1_3;
          p1 = cf2_3;
          XT_MOVF_S(p0, p1, b_lt05); cf3 = p0;
          p0 = cf1_4;
          p1 = cf2_4;
          XT_MOVF_S(p0, p1, b_lt05); cf4 = p0;
          p0 = cf1_5;
          p1 = cf2_5;
          XT_MOVF_S(p0, p1, b_lt05); cf5 = p0;
          p0 = cf1_6;
          p1 = cf2_6;
          XT_MOVF_S(p0, p1, b_lt05); cf6 = p0;
          p0 = cf1_7;
          p1 = cf2_7;
          XT_MOVF_S(p0, p1, b_lt05); cf7 = p0;
        }
        /* Compute the approximation to z(y) = tan(y)/y-1. We use
        * Horner's method for better pipelining of a few iterations. */
        z0 = cf0;
        XT_MADD_S(cf1, y0, z0); z0 = cf1;
        XT_MADD_S(cf2, y0, z0); z0 = cf2;
        XT_MADD_S(cf3, y0, z0); z0 = cf3;
        XT_MADD_S(cf4, y0, z0); z0 = cf4;
        XT_MADD_S(cf5, y0, z0); z0 = cf5;
        XT_MADD_S(cf6, y0, z0); z0 = cf6;
        XT_MADD_S(cf7, y0, z0); z0 = cf7;
        XT_MADD_S( y0, y0, z0); z0 = y0;
      
        /* Account for the range reduction. */
        b_gt1 = XT_OLT_S(XT_CONST_S(1), x0);
        z1 = XT_SUB_S(pi2f.f, z0);
        XT_MOVT_S(z0, z1, b_gt1);
      
        /* Propagate the input sign. */
        z0_i = XT_RFR(z0);
        z0_i = z0_i | sx;
        /* Load next input value */
        XT_LSIP(x0, X, sz_f32);
        /* Save output value */
		*Z++=z0_i;
      }
    }
  }
} /* vec_atanf() */
#endif
