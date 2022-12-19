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
 * Half precision floating-point arithmetic and manipulation functions.
 */

#ifndef __FLOAT16_H
#define __FLOAT16_H

/* Portable data types. */
#include "NatureDSP_types.h"
#include <fenv.h>

#ifdef __cplusplus
extern "C"
{
#endif

/* Half precision floating-point limits and characteristics. */
#define FLT16_DIG         3                  /* # of decimal digits of precision */
#define FLT16_EPSILON     (float16_t)0x1400  /* 2^-5 == 9.766e-4, diff. between 1 and the least representable number greater than 1 */
#define FLT16_MANT_DIG    11                 /* # of bits in mantissa */
#define FLT16_MAX         (float16_t)0x7bff  /* 65504, max finite value */
#define FLT16_MAX_10_EXP  4                  /* Max decimal exponent */
#define FLT16_MAX_EXP     16                 /* Max binary exponent */
#define FLT16_MIN         (float16_t)0x0400  /* 6.1035e-5, min positive normalized value */
#define FLT16_MIN_10_EXP  (-4)               /* Min decimal exponent */
#define FLT16_MIN_EXP     (-13)              /* Min binary exponent */

/* Half precision floating-point special numbers */
#define FLT16_SNAN        (float16_t)0x7c01  /* Signalling NaN */
#define FLT16_QNAN        (float16_t)0x7e00  /* Quiet NaN */
#define FLT16_MINUSINF    (float16_t)0xfc00  /* -Inf */
#define FLT16_PLUSINF     (float16_t)0x7c00  /* +Inf */

/* Conversions */
float16_t conv_f32_to_f16( float32_t x );
float16_t conv_f64_to_f16( float64_t x );
float32_t conv_f16_to_f32( float16_t x );
float64_t conv_f16_to_f64( float16_t x );
float16_t conv_i16_to_f16( int16_t x, int p ); /* x*2^-p, bi-directional */

/* Elementary arithmetics */
float16_t add_f16( float16_t x, float16_t y );
float16_t sub_f16( float16_t x, float16_t y );
float16_t mul_f16( float16_t x, float16_t y );
float16_t div_f16( float16_t x, float16_t y );

/* Fused multiply-add: return x*y+z. */
float16_t fma_f16( float16_t x, float16_t y, float16_t z );
/* Fused multiply-subtract: return z-x*y. */
float16_t fms_f16( float16_t x, float16_t y, float16_t z );
/* Return absolute value of half precision floating-point number. */
float16_t fabs_f16( float16_t x );
/* Multiply the first argument by two, raised to the power of the second argument. */
float16_t ldexp_f16( float16_t x, int n );
/* Round the argument to the nearest integer value in half precision floating-
 * point format, rounding halfway cases away from zero. */
float16_t round_f16( float16_t x );
/* Round the argument to the nearest integer value, rounding halfway cases 
 * away from zero. For infinite input argument, the function returns either
 * INT_MIN (for negative infinity on input) or INT_MAX (for positive infinity).
 * If the input argument is Not-a-Number (NaN), the returned value is not defined. */
int iround_f16( float16_t x );
/* Round the argument to integer value, nearest to but no larger in magnitude 
 * than the argument. For infinite input argument, the function returns either
 * INT_MIN (for negative infinity on input) or INT_MAX (for positive infinity).
 * If the input argument is Not-a-Number (NaN), the returned value is not defined. */
int itrunc_f16( float16_t x );
/* Break a floating-point number x into a normalized fraction f, 0.5<=f<1, and an
 * integral power of 2 e, such that x == f*2^e. For a numeric input x, return the
 * fraction f and set output argument *pe to the integral power e. For a non-numeric
 * input x, just return x. */
float16_t frexp_f16(float16_t x, int * pe);
/* Negate the argument. */
inline_ float16_t neg_f16( float16_t x ) { return (float16_t)((int)x^(1<<15)); }
/* Take the absolute value of the argument. */
inline_ float16_t abs_f16( float16_t x ) { return (float16_t)((int)x&~(1<<15)); }

/* Return nonzero value if argument is +/-infinity */
inline_ int isinf_f16( float16_t x )
{
  return ( 0x7c00 == ( (uint16_t)x & 0x7fff ) );
}

/* Return nonzero value if argument is not a Not-a-Number (NaN) as well as not +/-infinity. */
inline_ int isfinite_f16( float16_t x )
{
  return ( 0x7c00 != ( (uint16_t)x & 0x7c00 ) );
}

/* Return nonzero value if argument is a signaling Not-a-Number (sNaN). */
inline_ int is_snan_f16( float16_t x )
{
  return ( ( 0x7c00 == ( (uint16_t)x & 0x7c00 ) ) &&
           (      0 == ( (uint16_t)x & 0x0200 ) ) &&
           (      0 != ( (uint16_t)x & 0x01ff ) ) );
}

/* Return nonzero value if argument is a signaling or quiet Not-a-Number (NaN). */
inline_ int isnan_f16( float16_t x )
{
  if (is_snan_f16(x)) feraiseexcept(FE_INVALID);
  return ( ( 0x7c00 == ( (uint16_t)x & 0x7c00 ) ) &&
           (      0 != ( (uint16_t)x & 0x03ff ) ) );
}

/* Return value of the sign bit. */
inline_ int takesign_f16( float16_t x )
{
  return ( 0 != ( (uint16_t)x & 0x8000 ) );
}

/* Set the sign bit. */
inline_ float16_t setsign_f16( float16_t x, int sgn )
{
  return (float16_t)( ( (uint16_t)x & 0x7fff ) | ( (0!=sgn) << 15 ) );
}

/* Toggle the sign bit if sgn argument is nonzero. */
inline_ float16_t changesign_f16( float16_t x, int sgn )
{
  return (float16_t)( (uint16_t)x ^ ( (0!=sgn) << 15 ) );
}

/* Return nonzero value if arguments are unordered: a and/or b is NaN */
inline_ int un_f16 ( float16_t a, float16_t b )
{
  return ( isnan_f16(a) || isnan_f16(b) );
}

/* All comparison functions return zero for unordered arguments, but do not
 * raise any floating-point exceptions. */
int lt_f16  ( float16_t a, float16_t b ); /* a<b  */
int lte_f16 ( float16_t a, float16_t b ); /* a<=b */
int gt_f16  ( float16_t a, float16_t b ); /* a>b  */
int gte_f16 ( float16_t a, float16_t b ); /* a>=b */
int eq_f16  ( float16_t a, float16_t b ); /* a==b */
int oneq_f16( float16_t a, float16_t b ); /* a!=b, o stands for "ordered" */

#ifdef __cplusplus
};
#endif

#endif /* __FLOAT16_H */
