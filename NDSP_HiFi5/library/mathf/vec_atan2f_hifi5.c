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

#include <float.h>

/* DSP Library API */
#include "NatureDSP_Signal_math.h"
/* Common helper macros. */
#include "common_fpu.h"
/* Tables */
#include "pif_tbl.h"
#include "atanf_tbl.h"
/* +/-Infinity, single precision */
#include "inff_tbl.h"
/* sNaN/qNaN, single precision. */
#include "nanf_tbl.h"

#if !HAVE_VFPU && !HAVE_FPU
DISCARD_FUN(void,vec_atan2f,( float32_t * z, const float32_t * y, const float32_t * x,  int N ))
#elif HAVE_VFPU
#define sz_f32    (int)sizeof(float32_t)

/*
  NatureDSP Signal Processing Library. Vector Mathematics
  Full-Quadrant Arc Tangent
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/*-------------------------------------------------------------------------
  Full-Quadrant Arc Tangent
  The functions compute the arc tangent of the ratios y[N]/x[N] and store the
  result to output vector z[N]. 
  Floating point functions output is in radians. Fixed point functions
  scale its output by pi.

  NOTE:
  1. Scalar floating point function is compatible with standard ANSI C routines and set 
     errno and exception flags accordingly
  2. Scalar floating point function assigns EDOM to errno whenever y==0 and x==0.

  Accuracy:
  24 bit version: 768 (3.57e-7)
  floating point: 2 ULP

  Special cases:
       y    |   x   |  result   |  extra conditions    
    --------|-------|-----------|---------------------
     +/-0   | -0    | +/-pi     |
     +/-0   | +0    | +/-0      |
     +/-0   |  x    | +/-pi     | x<0
     +/-0   |  x    | +/-0      | x>0
     y      | +/-0  | -pi/2     | y<0
     y      | +/-0  |  pi/2     | y>0
     +/-y   | -inf  | +/-pi     | finite y>0
     +/-y   | +inf  | +/-0      | finite y>0
     +/-inf | x     | +/-pi/2   | finite x
     +/-inf | -inf  | +/-3*pi/4 | 
     +/-inf | +inf  | +/-pi/4   |

  Input:
    y[N]  vector of numerator values, Q31 or floating point
    x[N]  vector of denominator values, Q31 or floating point
    N     length of vectors
  Output:
    z[N]  results, Q31 or floating point

---------------------------------------------------------------------------*/
static void __atan2f( float32_t * z, const float32_t * y, const float32_t * x,  int N );
void vec_atan2f( float32_t * restrict z, const float32_t * restrict y, const float32_t * restrict x,  int N )
{
    xtfloatx4 * restrict pX;
    xtfloatx4 * restrict pY;
    xtfloatx4 * restrict pZ;
    int n;
    if (N<=0) return;
    if (N&7)
    {
        float32_t ALIGN(32) xbuf[8],zbuf[8],ybuf[8];
        pX=(xtfloatx4 *)xbuf;
        pY=(xtfloatx4 *)ybuf;
        pZ=(xtfloatx4 *)zbuf;
        AE_SSX2X2_I(XT_CONST_S(0),XT_CONST_S(0),pX,0*sizeof(xtfloatx4));
        AE_SSX2X2_I(XT_CONST_S(0),XT_CONST_S(0),pX,1*sizeof(xtfloatx4));
        AE_SSX2X2_I(XT_CONST_S(0),XT_CONST_S(0),pY,0*sizeof(xtfloatx4));
        AE_SSX2X2_I(XT_CONST_S(0),XT_CONST_S(0),pY,1*sizeof(xtfloatx4));
        AE_SSX2X2_I(XT_CONST_S(0),XT_CONST_S(0),pZ,0*sizeof(xtfloatx4));
        AE_SSX2X2_I(XT_CONST_S(0),XT_CONST_S(0),pZ,1*sizeof(xtfloatx4));
        for (n=0; n<(N&7); n++) 
        {
            xtfloat t;
            XT_LSIP(t,castxcc(xtfloat,x  ),sizeof(float32_t));
            XT_SSIP(t,castxcc(xtfloat,pX ),sizeof(float32_t));
            XT_LSIP(t,castxcc(xtfloat,y  ),sizeof(float32_t));
            XT_SSIP(t,castxcc(xtfloat,pY ),sizeof(float32_t));
        }
        pX=(xtfloatx4 *)xbuf;
        pY=(xtfloatx4 *)ybuf;
        __atan2f((float32_t*)pZ,(float32_t*)pY,(float32_t*)pX,8);
        for (n=0; n<(N&7); n++) 
        {
            xtfloat t;
            XT_LSIP(t,castxcc(xtfloat,pZ),sizeof(float32_t));
            XT_SSIP(t,castxcc(xtfloat,z ),sizeof(float32_t));
        }
        N&=~7;
    }
    if (N<=0) return;
    __atan2f(z,y,x,N);
}

static void __atan2f( float32_t * z, const float32_t * y, const float32_t * x,  int N )
{
  /*
    const union ufloat32uint32* p;
    int sx,sy,big;
    sx=takesignf(x);
    sy=takesignf(y);
    x=fabs(x);
    y=fabs(y);
    if(x==0.f && y==0.f)
    {
      // The actual result depends on input signs.
      x = 1.f;
      y = 0.f;
    }

    big=x>y;
    if(big)
    {
        x=y/x;
    }
    else
    {
      // compare x==y is necessary to support (+/-Inf, +/-Inf) cases 
      x = (x == y) ? 1.0f : x / y;
    }
    p = (x<0.5f) ? atanftbl1 : atanftbl2;
    // approximate atan(x)/x-1 
    y = p[0].f;
    y = x*y + p[1].f;
    y = x*y + p[2].f;
    y = x*y + p[3].f;
    y = x*y + p[4].f;
    y = x*y + p[5].f;
    y = x*y + p[6].f;
    y = x*y + p[7].f;
    // convert result to true atan(x) 
    y = x*y + x;

    if (!big) y = pi2f.f - y;
    if (sx)   y = pif.f - y;
    if (sy)   y = -y;
    return   y;
  */
  const xtfloat   *          T1;
  const xtfloat   *          T2;
  const xtfloatx4  *          X;
  const xtfloatx4  *          Y;
        xtfloatx4  * restrict Z;
  const xtfloatx4 *          S_rd;
        xtfloatx4 * restrict S_wr;

  ae_valignx2 X_va, Y_va, Z_va;


  /* Current block index; overall number of blocks; number of values in the current block */
  /* Block size, blkLen <= blkSize */
  const int blkSize = MAX_ALLOCA_SZ / 2 * sz_f32;
  /* Allocate a fixed-size scratch area on the stack. */
  float32_t ALIGN(32) scr[blkSize];
  int n,M;

  NASSERT(N%8==0);
  NASSERT_ALIGN16(scr);

  /*
  * Data are processed in blocks of scratch area size. Further, the algorithm
  * implementation is splitted in order to feed the optimizing compiler with a
  * few loops of managable size.
  */

  for (; N>0; N-=M, x+=M,y+=M,z+=M)
  {
    M = XT_MIN( N , blkSize );

    /*
    * Part I, reduction to [0,pi/4]. Reference C code:
    *
    *   {
    *     float32_t x0, y0, p0;
    *   
    *     for ( n=0; n<blkLen; n++ )
    *     {
    *       x0 = fabsf( x[blkIx*blkSize+n] );
    *       y0 = fabsf( y[blkIx*blkSize+n] );
    *   
    *       // The actual result depends on input signs.
    *       if ( x0==0.f && y0==0.f ) { x0 = 1.f; y0 = 0.f; };
    *   
    *       if ( x0>y0 ) p0 = y0/x0;
    *       // Special case of x==y is necessary to support (+/-Inf, +/-Inf) cases.
    *       else p0 = ( x0==y0 ? 1.f : x0/y0 );
    *   
    *       scr[n] = p0;
    *     }
    *   }
    */
        X    = (xtfloatx4*)(x);
        Y    = (xtfloatx4*)(y);
        S_wr = (xtfloatx4*)scr;
        X_va = AE_LA128_PP(X);
        Y_va = AE_LA128_PP(Y);
        __Pragma( "loop_count min=1" );
        /*21*/
        __Pragma("loop_count factor=2")
        for ( n=0; n<(M>>2); n++ )
        {
            /* Input values */
            xtfloatx2 x0, y0, x1, y1;
            /* Numerator; denominator; reciprocal; quotient */
            xtfloatx2 num0, den0, rcp0, quo0;
            xtfloatx2 num1, den1, rcp1, quo1;
            /* Scaling factor; error term */
            xtfloatx2 scl0, scl1, eps0, eps1;
            /* Is NaN; Inf/Inf; x/Inf; 0/0; x and y are subnormal */
            xtbool2 b0_nan, b0_num_inf, b0_den_inf, b0_eqz, b0_subn;
            xtbool2 b1_nan, b1_num_inf, b1_den_inf, b1_eqz, b1_subn;
            AE_LASX2X2_IP(x0, x1, X_va, X);
            AE_LASX2X2_IP(y0, y1, Y_va, Y);
            /* Replicate NaNs in both x and y to ensure NaN propagation. */
            b0_nan = XT_UN_SX2(x0, y0);
            b1_nan = XT_UN_SX2(x1, y1);
            ABS_SX2X2(x0, y0, x0, y0);
            ABS_SX2X2(x1, y1, x1, y1);
            XT_MOVT_SX2(x0, qNaNf.f, b0_nan);
            XT_MOVT_SX2(y0, qNaNf.f, b0_nan);
            XT_MOVT_SX2(x1, qNaNf.f, b1_nan);
            XT_MOVT_SX2(y1, qNaNf.f, b1_nan);
            /* num <= den */
            num0 = XT_MIN_SX2(x0, y0);
            num1 = XT_MIN_SX2(x1, y1);
            den0 = XT_MAX_SX2(y0, x0);
            den1 = XT_MAX_SX2(y1, x1);
            /* Scale up numerator and denominator if BOTH are subnormal. */
            b0_subn = XT_OLT_SX2(num0, FLT_MIN);
            b1_subn = XT_OLT_SX2(num1, FLT_MIN);
            scl0 = scl1 = (xtfloatx2)8388608.f; 
            XT_MOVF_SX2(scl0, XT_CONST_S(1), b0_subn);
            XT_MOVF_SX2(scl1, XT_CONST_S(1), b1_subn);
            MUL_SX2X2(num0, num1, num0, num1, scl0, scl1);
            MUL_SX2X2(den0, den1, den0, den1, scl0, scl1);
            /* Classify numerator and denominator. */
            b0_num_inf = XT_OEQ_SX2(num0, plusInff.f);           /* Inf/Inf */
            b0_den_inf = XT_OEQ_SX2(den0, plusInff.f);           /* x/Inf   */
            b0_eqz = XT_OEQ_SX2(den0, XT_CONST_S(0)); /* 0/0     */
            /* Classify numerator and denominator. */
            b1_num_inf = XT_OEQ_SX2(num1, plusInff.f);           /* Inf/Inf */
            b1_den_inf = XT_OEQ_SX2(den1, plusInff.f);           /* x/Inf   */
            b1_eqz = XT_OEQ_SX2(den1, XT_CONST_S(0)); /* 0/0     */
            /* Initial appromimation for 1/den. */
            rcp0 = XT_RECIP0_SX2(den0);
            rcp1 = XT_RECIP0_SX2(den1);
            eps0 = eps1 = XT_CONST_S(1);
            MSUB_SX2X2(eps0, eps1, rcp0, rcp1, den0, den1);
            MADD_SX2X2(rcp0, rcp1, rcp0, rcp1, eps0, eps1);
            /* Approximation for the quotient num/den. */
            MUL_SX2X2(quo0, quo1, num0, num1, rcp0, rcp1);
            /* Refine the quotient by a modified Newton-Raphson iteration. */
            eps0 = num0; eps1 = num1;
            MSUB_SX2X2(eps0, eps1, quo0, quo1, den0, den1);
            MADD_SX2X2(quo0, quo1, rcp0, rcp1, eps0, eps1);
            /* Force conventional results for special cases. */
            XT_MOVT_SX2(quo0, XT_CONST_S(0), b0_den_inf); /* x/Inf -> 0   */
            XT_MOVT_SX2(quo1, XT_CONST_S(0), b1_den_inf); /* x/Inf -> 0   */
            XT_MOVT_SX2(quo1, XT_CONST_S(1), b1_num_inf); /* Inf/Inf -> 1 */
            XT_MOVT_SX2(quo0, XT_CONST_S(1), b0_num_inf); /* Inf/Inf -> 1 */
            XT_MOVT_SX2(quo0, XT_CONST_S(0), b0_eqz); /* 0/0 -> 0     */
            XT_MOVT_SX2(quo1, XT_CONST_S(0), b1_eqz); /* 0/0 -> 0     */
            AE_SSX2X2_IP(quo0, quo1, S_wr, +4 * sz_f32);
        }
        __Pragma("no_reorder");
     /*
     * Part II, polynomial approximation and full quadrant restoration.
     * Reference C code:
     *
     *   {
     *     const union ufloat32uint32 * ptbl;
     *     float32_t x0, y0, z0, p0;
     *     int sx, sy;
     *   
     *     for ( n=0; n<blkLen; n++ )
     *     {
     *       x0 = x[blkIx*blkSize+n];
     *       y0 = y[blkIx*blkSize+n];
     *       p0 = scr[n];
     *   
     *       sx = takesignf( x0 ); x0 = fabsf( x0 );
     *       sy = takesignf( y0 ); y0 = fabsf( y0 );
     *   
     *       ptbl = ( p0<0.5f ? atanftbl1 : atanftbl2 );
     *    
     *       // Approximate atan(p)/p-1
     *       z0 = ptbl[0].f;
     *       z0 = ptbl[1].f + p0*z0;
     *       z0 = ptbl[2].f + p0*z0;
     *       z0 = ptbl[3].f + p0*z0;
     *       z0 = ptbl[4].f + p0*z0;
     *       z0 = ptbl[5].f + p0*z0;
     *       z0 = ptbl[6].f + p0*z0;
     *       z0 = ptbl[7].f + p0*z0;
     *       z0 =        p0 + p0*z0;
     *   
     *       if ( x0<y0 ) z0 = pi2f.f - z0;
     *       if ( sx    ) z0 = pif.f - z0;
     *       if ( sy    ) z0 = -z0;
     *   
     *       z[blkIx*blkSize+n] = z0;
     *     }
     *   }
     */
        T1   = (xtfloat  *)atanftbl1;
        T2   = (xtfloat  *)atanftbl2;
        X = (xtfloatx4*)(x);
        Y = (xtfloatx4*)(y);
        Z = (xtfloatx4*)(z);
        S_rd = (xtfloatx4*)scr;
        X_va = AE_LA128_PP(X);
        Y_va = AE_LA128_PP(Y);
        Z_va = AE_ZALIGN128();
        __Pragma("loop_count factor=2")
        for ( n=0; n<(M>>2); n++ )
        {
            /* Input values; output value; reducted input value and its 2nd power. */
            xtfloatx2 x0, y0, x1, y1, z0, z1, p0, p1;
            xtfloat cf;
            xtfloatx2/* n0, n1,*/ f0, f1;
            /* Temporary; input value signs */
            ae_int32x2 t0, sx0, sy0, t1, sx1, sy1;
            /* Polynomial coeffs for 0.f<=p<0.5f (#1) and 0.5f<=p<=1.f (#2). */
            xtfloatx2 cf1_0, cf1_1, cf1_2, cf1_3, cf1_4, cf1_5, cf1_6, cf1_7;
            xtfloatx2 cf2_0, cf2_1, cf2_2, cf2_3, cf2_4, cf2_5, cf2_6, cf2_7;
            /* Selected polynomial coeffs. */
            //xtfloatx2 cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7;
            /* x less than y; x is negative; p is less than 0.5f. */
            xtbool2 b0_xlty, b0_sx, b0_lt05;
            xtbool2 b1_xlty, b1_sx, b1_lt05;
            xtfloatx2 af1_0,af1_1,af1_2,af1_3;
            xtfloatx2 af2_0,af2_1,af2_2,af2_3;
            xtfloatx2 bf1_0,bf1_1,bf1_2,bf1_3;
            xtfloatx2 bf2_0,bf2_1,bf2_2,bf2_3;
            xtfloatx2 a0,a1,b0,b1;

            AE_LASX2X2_IP(x0, x1, X_va, X);
            AE_LASX2X2_IP(y0, y1, Y_va, Y); 

            /* Keep sign of x as a boolean. */
            sx0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(x0);
            b0_sx = AE_LT32(sx0, AE_ZERO32());
            sx1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(x1);
            b1_sx = AE_LT32(sx1, AE_ZERO32());
            /* Keep y sign as a binary value. */
            sy0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(y0);
            //sy0 = AE_AND32(sy0, 0x80000000);
            sy0 = AE_SRLI32(sy0, 31);
            sy0 = AE_SLLI32(sy0, 31);
            sy1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(y1);
            //sy1 = AE_AND32(sy1, 0x80000000);
            sy1 = AE_SRLI32(sy1, 31);
            sy1 = AE_SLLI32(sy1, 31);
            ABS_SX2X2(x0, x1, x0, x1);
            ABS_SX2X2(y0, y1, y0, y1);
            b0_xlty = XT_OLT_SX2(x0, y0);
            b1_xlty = XT_OLT_SX2(x1, y1);
            //
            AE_LSX2X2_IP(p0, p1, S_rd, +4 * sz_f32);
            b0_lt05 = XT_OLT_SX2(p0, (xtfloatx2)0.5f);
            b1_lt05 = XT_OLT_SX2(p1, (xtfloatx2)0.5f);

            /* Reload coeff sets on each iteration. */
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
            /*
            * Compute the approximation to z(p) = atan(p)/p-1. Here we use a combination
            * of Estrin's rule and Horner's method of polynomial evaluation to shorten the
            * dependency path at the cost of additional multiplication.
            */
#if 0
            n0 = cf1_0; XT_MOVF_SX2(n0, cf2_0, b0_lt05);
            n1 = cf1_0; XT_MOVF_SX2(n1, cf2_0, b1_lt05);
            f0 = cf1_1; XT_MOVF_SX2(f0, cf2_1, b0_lt05);
            f1 = cf1_1; XT_MOVF_SX2(f1, cf2_1, b1_lt05);
            MADD_SX2X2(f0, f1, n0, n1, p0, p1); cf1_0 = f0; cf2_0 = f1;
            n0 = cf1_2; XT_MOVF_SX2(n0, cf2_2, b0_lt05);
            n1 = cf1_2; XT_MOVF_SX2(n1, cf2_2, b1_lt05);
            f0 = cf1_3; XT_MOVF_SX2(f0, cf2_3, b0_lt05);
            f1 = cf1_3; XT_MOVF_SX2(f1, cf2_3, b1_lt05);
            MADD_SX2X2(f0, f1, n0, n1, p0, p1); cf1_1 = f0; cf2_1 = f1;
            n0 = cf1_4; XT_MOVF_SX2(n0, cf2_4, b0_lt05);
            n1 = cf1_4; XT_MOVF_SX2(n1, cf2_4, b1_lt05);
            f0 = cf1_5; XT_MOVF_SX2(f0, cf2_5, b0_lt05);
            f1 = cf1_5; XT_MOVF_SX2(f1, cf2_5, b1_lt05);
            MADD_SX2X2(f0, f1, n0, n1, p0, p1); cf1_2 = f0; cf2_2 = f1;
            n0 = cf1_6; XT_MOVF_SX2(n0, cf2_6, b0_lt05);
            n1 = cf1_6; XT_MOVF_SX2(n1, cf2_6, b1_lt05);
            f0 = cf1_7; XT_MOVF_SX2(f0, cf2_7, b0_lt05);
            f1 = cf1_7; XT_MOVF_SX2(f1, cf2_7, b1_lt05);
            MADD_SX2X2(f0, f1, n0, n1, p0, p1); cf1_3 = f0; cf2_3 = f1;

            MUL_SX2X2(x0, x1, p0, p1, p0, p1);
            n0 = cf1_0; n1 = cf2_0;
            f0 = cf1_1; f1 = cf2_1;
            MADD_SX2X2(f0, f1, n0, n1, x0, x1); n0 = f0; n1 = f1;
            f0 = cf1_2; f1 = cf2_2;
            MADD_SX2X2(f0, f1, n0, n1, x0, x1); n0 = f0; n1 = f1;
            f0 = cf1_3; f1 = cf2_3;
            MADD_SX2X2(f0, f1, n0, n1, x0, x1); n0 = f0; n1 = f1;

            MADD_SX2X2(p0, p1, n0, n1, p0, p1); f0 = p0; f1 = p1;
#else
            /*
            * Compute the approximation to z(y) = atan(y)/y-1. Here we use a combination
            * of Estrin's rule and Horner's method of polynomial evaluation to shorten the
            * dependency path at the cost of additional multiplication.
            */
            /* compute 2 different polynomials and select one of them */
            af1_0=cf1_1; af2_0=cf2_1; MADDQ_S(af1_0, af2_0, cf1_0, cf2_0, p0); 
            bf1_0=cf1_1; bf2_0=cf2_1; MADDQ_S(bf1_0, bf2_0, cf1_0, cf2_0, p1); 
            af1_1=cf1_3; af2_1=cf2_3; MADDQ_S(af1_1, af2_1, cf1_2, cf2_2, p0); 
            bf1_1=cf1_3; bf2_1=cf2_3; MADDQ_S(bf1_1, bf2_1, cf1_2, cf2_2, p1); 
            af1_2=cf1_5; af2_2=cf2_5; MADDQ_S(af1_2, af2_2, cf1_4, cf2_4, p0); 
            bf1_2=cf1_5; bf2_2=cf2_5; MADDQ_S(bf1_2, bf2_2, cf1_4, cf2_4, p1); 
            af1_3=cf1_7; af2_3=cf2_7; MADDQ_S(af1_3, af2_3, cf1_6, cf2_6, p0); 
            bf1_3=cf1_7; bf2_3=cf2_7; MADDQ_S(bf1_3, bf2_3, cf1_6, cf2_6, p1); 

            MUL_SX2X2(x0, x1, p0, p1, p0, p1); 
            a0 = af1_1; a1 = af2_1; MADDQ_S(a0, a1, af1_0, af2_0,x0); 
            b0 = bf1_1; b1 = bf2_1; MADDQ_S(b0, b1, bf1_0, bf2_0,x1); 
            f0 = af1_2; f1 = af2_2; MADDQ_S(f0, f1, a0, a1, x0); a0 = f0; a1 = f1;
            f0 = bf1_2; f1 = bf2_2; MADDQ_S(f0, f1, b0, b1, x1); b0 = f0; b1 = f1;
            f0 = af1_3; f1 = af2_3; MADDQ_S(f0, f1, a0, a1, x0); a0 = f0; a1 = f1;
            f0 = bf1_3; f1 = bf2_3; MADDQ_S(f0, f1, b0, b1, x1); b0 = f0; b1 = f1;
            f0=p0; f1=p0;           MADDQ_S(f0, f1, a0, a1, p0); a0 = f0; a1 = f1;
            f0=p1; f1=p1;           MADDQ_S(f0, f1, b0, b1, p1); b0 = f0; b1 = f1;

            f0 = a0; XT_MOVF_SX2(f0, a1, b0_lt05);
            f1 = b0; XT_MOVF_SX2(f1, b1, b1_lt05);
#endif
            SUB_SX2X2(z0, z1, pi2f.f, pi2f.f, f0, f1);
            XT_MOVT_SX2(f0, z0, b0_xlty);
            XT_MOVT_SX2(f1, z1, b1_xlty);

            SUB_SX2X2(z0, z1, pif.f, pif.f, f0, f1);
            XT_MOVT_SX2(f0, z0, b0_sx);
            XT_MOVT_SX2(f1, z1, b1_sx);
        
            t0 = XT_AE_MOVINT32X2_FROMXTFLOATX2(f0);
            t1 = XT_AE_MOVINT32X2_FROMXTFLOATX2(f1);

            t0 = AE_XOR32(t0, sy0);
            t1 = AE_XOR32(t1, sy1);
            z0 = XT_AE_MOVXTFLOATX2_FROMINT32X2(t0);
            z1 = XT_AE_MOVXTFLOATX2_FROMINT32X2(t1);
            AE_SASX2X2_IP(z0, z1, Z_va, Z);
        }
        AE_SA128POS_FP( Z_va, Z); 
    }
}
#elif HAVE_FPU
#define sz_f32    (int)sizeof(float32_t)

/*===========================================================================
  Scalar matematics:
  scl_atan2          full quadrant Arctangent        
===========================================================================*/
/*-------------------------------------------------------------------------
Floating-Point Full-Quadrant Arc Tangent
The functions compute the full quadrant arc tangent of the ratio y/x. 
Floating point functions output is in radians. Fixed point functions 
scale its output by pi.

NOTE:
1. Scalar function is compatible with standard ANSI C routines and set 
   errno and exception flags accordingly
2. Scalar function assigns EDOM to errno whenever y==0 and x==0.

Special cases:
     y    |   x   |  result   |  extra conditions    
  --------|-------|-----------|---------------------
   +/-0   | -0    | +/-pi     |
   +/-0   | +0    | +/-0      |
   +/-0   |  x    | +/-pi     | x<0
   +/-0   |  x    | +/-0      | x>0
   y      | +/-0  | -pi/2     | y<0
   y      | +/-0  |  pi/2     | y>0
   +/-y   | -inf  | +/-pi     | finite y>0
   +/-y   | +inf  | +/-0      | finite y>0
   +/-inf | x     | +/-pi/2   | finite x
   +/-inf | -inf  | +/-3*pi/4 | 
   +/-inf | +inf  | +/-pi/4   |

Input:
  y[N]  input data, Q15 or floating point
  x[N]  input data, Q15 or floating point
  N     length of vectors
Output:
  z[N]  result, Q15 or floating point
  
Restrictions:
x, y, z should not overlap
---------------------------------------------------------------------------*/

// Taken from Fusion
void vec_atan2f( float32_t * z, const float32_t * y, const float32_t * x, int N )
{
  /*
  * const union ufloat32uint32* p;
  * int sx,sy,big;
  * sx=takesignf(x);
  * sy=takesignf(y);
  * x=fabs(x);
  * y=fabs(y);
  * if(x==0.f && y==0.f)
  * {
  * // The actual result depends on input signs.
  * x = 1.f;
  * y = 0.f;
  * }
  * 
  * big=x>y;
  * if(big)
  * {
  * x=y/x;
  * }
  * else
  * {
  * // compare x==y is necessary to support (+/-Inf, +/-Inf) cases
  * x = (x == y) ? 1.0f : x / y;
  * }
  * p = (x<0.5f) ? atanftbl1 : atanftbl2;
  * // approximate atan(x)/x-1
  * y = p[0].f;
  * y = x*y + p[1].f;
  * y = x*y + p[2].f;
  * y = x*y + p[3].f;
  * y = x*y + p[4].f;
  * y = x*y + p[5].f;
  * y = x*y + p[6].f;
  * y = x*y + p[7].f;
  * // convert result to true atan(x)
  * y = x*y + x;
  * 
  * if (!big) y = pi2f.f - y;
  * if (sx)   y = pif.f - y;
  * if (sy)   y = -y;
  * return   y;
  */
  const xtfloat * restrict X;
  const xtfloat * restrict Y;
        int32_t * restrict Z;
  const xtfloat * restrict S_rd;
        xtfloat * restrict S_wr;
  const xtfloat * restrict POLY_TBL1;
  const xtfloat * restrict POLY_TBL2;

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
    * Part I, reduction to [0,pi/4]. Reference C code:
    *
    *   {
    *     float32_t x0, y0, p0;
    *
    *     for ( n=0; n<blkLen; n++ )
    *     {
    *       y0 = fabsf( y[blkIx*blkSize+n] );
    *       x0 = fabsf( x[blkIx*blkSize+n] );
    *
    *       // The actual result depends on input signs.
    *       if ( x0==0.f && y0==0.f ) { x0 = 1.f; y0 = 0.f; };
    *
    *       if ( x0>y0 ) p0 = y0/x0;
    *       // Special case of x==y is necessary to support (+/-Inf, +/-Inf) cases.
    *       else p0 = ( x0==y0 ? 1.f : x0/y0 );
    *
    *       scr[n] = p0;
    *     }
    *   }
    */

    {
      /* Input values */
      xtfloat x0, y0, i0;
      /* Numerator; denominator; reciprocal; quotient */
      xtfloat num, den, rcp, quo;
      /* Auxiliary vars */
      xtfloat s, eps;
      /* Is NaN; Inf/Inf; x/Inf; 0/0; x and y are subnormal */
      xtbool b_nan, b_num_inf, b_den_inf, b_eqz, b_subn;
      const xtfloat * pT;

      X = (xtfloat*)((uintptr_t)x + blkIx*blkSize*sz_f32);
      Y = (xtfloat*)((uintptr_t)y + blkIx*blkSize*sz_f32);
      S_wr = (xtfloat*)scr;

      static const uint32_t TAB[4] = { 0x7fc00000, 0x00800000,
        0x4b000000, 0x7f800000
      };
      pT = (xtfloat *)TAB;
      __Pragma("loop_count min=1");
      for (n = 0; n<blkLen; n++)
      {
        XT_LSIP(x0, X, sz_f32);
        XT_LSIP(y0, Y, sz_f32);

        /* Reproduce NaN in both x and y to ensure NaN propagation. */
        b_nan = XT_UN_S(x0, y0);
        i0 = pT[0];
        
        XT_MOVT_S(x0, i0, b_nan);

        x0 = XT_ABS_S(x0);
        y0 = XT_ABS_S(y0);

        /* num <= den */
        num = XT_MIN_S(y0, x0);
        den = XT_MAX_S(y0, x0);

        /* Classify numerator and denominator. */
        i0 = pT[1];
        b_subn = XT_OLT_S(num, i0);
        
        /* Scale up numerator and denominator if BOTH are subnormal. */
        i0 = pT[2];
        s = XT_MUL_S(num, i0); XT_MOVT_S(num, s, b_subn);
        s = XT_MUL_S(den, i0); XT_MOVT_S(den, s, b_subn);

        /* Initial appromimation for 1/den. */
        rcp = XT_RECIP0_S(den);
        /* Newton-Raphson iteration for 1/den. */
        eps = XT_CONST_S(1);
        XT_MSUB_S(eps, rcp, den);
        XT_MADD_S(rcp, rcp, eps);
        /* Approximation for the quotient num/den. */
        quo = XT_MUL_S(num, rcp);
        /* Refine the quotient by a modified Newton-Raphson iteration. */
        eps = num;
        XT_MSUB_S(eps, quo, den);
        XT_MADD_S(quo, rcp, eps);

        i0 = pT[3];
        b_num_inf = XT_OEQ_S(num, i0); /* Inf/Inf! */
        b_den_inf = XT_OEQ_S(den, i0);
        b_eqz = XT_OEQ_S(den, XT_CONST_S(0)); /* 0/0! */
        b_eqz = XT_ORB(b_eqz, b_den_inf);

        XT_MOVT_S(quo, XT_CONST_S(0), b_eqz);     /* 0/0 -> 0 or x/Inf -> 0*/
        XT_MOVT_S(quo, XT_CONST_S(1), b_num_inf); /* Inf/Inf -> 1 */

        XT_SSIP(quo, S_wr, sz_f32);
      }
    }
    __Pragma("no_reorder");

    /*
    * Part II, polynomial approximation and full quadrant restoration.
    * Reference C code:
    *
    *   {
    *     const union ufloat32uint32 * ptbl;
    *     float32_t x0, y0, z0, p0;
    *     int sx, sy;
    *
    *     for ( n=0; n<blkLen; n++ )
    *     {
    *       y0 = y[blkIx*blkSize+n];
    *       x0 = x[blkIx*blkSize+n];
    *       p0 = scr[n];
    *
    *       sy = takesignf( y0 ); y0 = fabsf( y0 );
    *       sx = takesignf( x0 ); x0 = fabsf( x0 );
    *
    *       ptbl = ( p0<0.5f ? atanftbl1 : atanftbl2 );
    *
    *       // Approximate atan(p)/p-1
    *       z0 = ptbl[0].f;
    *       z0 = ptbl[1].f + p0*z0;
    *       z0 = ptbl[2].f + p0*z0;
    *       z0 = ptbl[3].f + p0*z0;
    *       z0 = ptbl[4].f + p0*z0;
    *       z0 = ptbl[5].f + p0*z0;
    *       z0 = ptbl[6].f + p0*z0;
    *       z0 = ptbl[7].f + p0*z0;
    *       z0 =        p0 + p0*z0;
    *
    *       if ( x0<y0 ) z0 = pi2f.f - z0;
    *       if ( sx    ) z0 = pif.f - z0;
    *       if ( sy    ) z0 = -z0;
    *
    *       z[blkIx*blkSize+n] = z0;
    *     }
    *   }
    */
    {
      const xtfloat   *          pT;
      /* Input values; output value; reducted input value*/
      xtfloat x0, y0, z0, z1, p0;
      /* Temporary; input values' sign */
      int32_t sx, sy;
      /* Polynomial coeffs for 0.f<=p<0.5f (#1) and 0.5f<=p<=1.f (#2). */
      xtfloat cf1_0, cf1_1, cf1_2, cf1_3, cf1_4, cf1_5, cf1_6, cf1_7;
      xtfloat cf2_0, cf2_1, cf2_2, cf2_3, cf2_4, cf2_5, cf2_6, cf2_7;
      /* Selected polynomial coeffs. */
      xtfloat cf0, cf1, cf2, cf3, cf4, cf5, cf6, cf7;
      /* x less than y; x is negative; num/den is less than 0.5f. */
      xtbool b_xlty, b_sx, b_lt05;

      X = (xtfloat*)((uintptr_t)x + blkIx*blkSize*sz_f32);
      Y = (xtfloat*)((uintptr_t)y + blkIx*blkSize*sz_f32);
      Z = (int32_t*)((uintptr_t)z + blkIx*blkSize*sz_f32);

      S_rd = (xtfloat*)scr;
      /* pi/2, pi */
      static const uint32_t TAB[2] = { 0x3fc90fdb, 0x40490fdb
      };
      pT = (xtfloat  *)TAB;
      __Pragma("loop_count min=1");
      for (n = 0; n<blkLen; n++)
      {
        xtfloat i0;
        XT_LSIP(x0, X, 0*sz_f32);
        XT_LSIP(y0, Y, 0*sz_f32);

        x0 = XT_ABS_S(x0);
        y0 = XT_ABS_S(y0);
        b_xlty = XT_OLT_S(x0, y0);

        XT_LSIP(p0, S_rd, sz_f32);

        b_lt05 = XT_OLT_S(p0, XT_CONST_S(3));

        /*Reload polynomial coeff sets. */

        cf1_0 = XT_LSI(POLY_TBL1, 0 * sz_f32);
        cf2_0 = XT_LSI(POLY_TBL2, 0 * sz_f32);
        cf1_1 = XT_LSI(POLY_TBL1, 1 * sz_f32);
        cf2_1 = XT_LSI(POLY_TBL2, 1 * sz_f32);
        cf1_2 = XT_LSI(POLY_TBL1, 2 * sz_f32);
        cf2_2 = XT_LSI(POLY_TBL2, 2 * sz_f32);
        cf1_3 = XT_LSI(POLY_TBL1, 3 * sz_f32);
        cf2_3 = XT_LSI(POLY_TBL2, 3 * sz_f32);
        cf1_4 = XT_LSI(POLY_TBL1, 4 * sz_f32);
        cf2_4 = XT_LSI(POLY_TBL2, 4 * sz_f32);
        cf1_5 = XT_LSI(POLY_TBL1, 5 * sz_f32);
        cf2_5 = XT_LSI(POLY_TBL2, 5 * sz_f32);
        cf1_6 = XT_LSI(POLY_TBL1, 6 * sz_f32);
        cf2_6 = XT_LSI(POLY_TBL2, 6 * sz_f32);
        cf1_7 = XT_LSI(POLY_TBL1, 7 * sz_f32);
        cf2_7 = XT_LSI(POLY_TBL2, 7 * sz_f32);

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

        /* Compute the approximation to z(p) = tan(p)/p-1. We use
        * Horner's method for better pipelining of a few iterations. */
        z0 = cf0;
        XT_MADD_S(cf1, p0, z0); z0 = cf1;
        XT_MADD_S(cf2, p0, z0); z0 = cf2;
        XT_MADD_S(cf3, p0, z0); z0 = cf3;
        XT_MADD_S(cf4, p0, z0); z0 = cf4;
        XT_MADD_S(cf5, p0, z0); z0 = cf5;
        XT_MADD_S(cf6, p0, z0); z0 = cf6;
        XT_MADD_S(cf7, p0, z0); z0 = cf7;

        XT_MADD_S( p0, p0, z0); z0 = p0;

        /* Keep signs of x and y. */
        sx = (int32_t)((int *)X)[0]; X++;
        sy = (int32_t)((int *)Y)[0]; Y++;

		sy = sy & 0x80000000;

        b_sx = AE_INT64_LT(AE_MOVINT64_FROMINT32(sx), AE_ZERO64());

        /* if ( x0<y0 ) z0 = pi2f.f - z0; */
        i0 = XT_LSI(pT, 0*sz_f32);
        z1 = XT_SUB_S(i0, z0); XT_MOVT_S(z0, z1, b_xlty);
        /* if ( sx ) z0 = pif.f - z0; */
        i0 = XT_LSI(pT, 1*sz_f32);
        z1 = XT_SUB_S(i0, z0); XT_MOVT_S(z0, z1, b_sx);
        /* if ( sy ) z0 = -z0; */
        sx = XT_RFR(z0);
        sx = sx ^ sy;
        *Z++ = sx;
      }
    }
  }
} /* vec_atan2f() */
#endif
