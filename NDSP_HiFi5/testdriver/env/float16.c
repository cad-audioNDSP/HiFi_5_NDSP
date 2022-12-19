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

/* Partable data types. */
#include "types.h"
/* Half precision floating-point arithmetic and manipulation functions. */
#include "float16.h"
#include <fenv.h>
#include <math.h>
#include <limits.h>

#define MIN(a,b)  ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)  ( (a)>(b) ? (a) : (b) )


static int16_t _S_exp0_l (int32_t x)
{
  int16_t z=0;
  if ( x==0 )  return 0x1F;
  while ( (int32_t)(x^(x<<1))>0 ) //MSB != MSB-1
  {
    x<<=1;
    z++;
  }
  return z;
}

/* Convert single precision floating-point data to half precision format. */
float16_t conv_f32_to_f16( float32_t x )
{
  union ufloat32uint32 _x;
  int32_t s, e, f;
  int32_t tb, rbit, tail;

  _x.f = x;

  s = (int32_t)(_x.u>>31);
  e = ( (int32_t)(_x.u>>23) & 255L ) - 127L;
  
  /* Restore the implicit leading significand bit. */
  f = (int32_t)( _x.u & ( (1L<<23) - 1 ) );

  if ( 128==e && 0!=f ) /* NaN */
  {
    /* Return a quiet (!) NaN with original payload bits. */
    return (float16_t)( (s<<15) | 0x7e00 | (f>>(23-10)) ) ; 
  }
  else if ( e>15 ) /* Too large magnitude */
  {
    return (float16_t)( (s<<15) | 0x7c00 );
  }
  else if ( e<-25 ) /* Too small magnitude */
  {
    return (float16_t)(s<<15);
  }
  else if ( -14>e ) /* Subnormal half precision value */
  {
    /* Restore the implicit leading significand bit. */ 
    f += (1L<<23);
    /* 10 significand trailing bits */
    tb = ( f >> (23-10-e-14) ); 
    /* Rounding bit value */
    rbit = ( f >> (23-10-1-e-14) ) & 1L; 
    /* Are there any nonzero bits left? */
    tail = ( 0 != ( f & ( (1L<<(23-10-1-e-14)) - 1 ) ) );
    e = -15;
  }
  else /* Normal half precision value */
  {
    /* 10 significand trailing bits */
    tb = ( f >> (23-10) ); 
    /* Rounding bit value */
    rbit = ( f >> (23-10-1) ) & 1L; 
    /* Are there any nonzero bits left? */
    tail = ( 0 != ( f & ( (1L<<(23-10-1)) - 1 ) ) );
  }

  if ( 0!=rbit ) 
  {
    /* Round to nearest. In case of a tie round to even mantissa. */
    tb += ( 0==tail ? (tb&1) : 1 );
    /* Look if the exponent shall be adjusted for the rounded mantissa. */
    if ( tb & (1L<<10) )
    { 
      /* Correct the exponent and mantissa, and check if the rounded value overflows. */
      tb = ( ++e<=15 ? ( (tb+1024) >> 1 ) : 0 );
    }
  }

  return (float16_t)( (s<<15) | ( (e+15) << 10 ) | ( tb & 0x3ff ) );

} /* conv_f32_to_f16() */

/* Convert double precision floating-point data to half precision format. */
float16_t conv_f64_to_f16( float64_t x )
{
  union ufloat64uint64 _x;
  int32_t s, e;
  uint64_t f;
  int32_t tb, rbit, tail;

  _x.f = x;

  s = (int32_t)(_x.u>>63);
  e = ( (int32_t)(_x.u>>52) & 2047L ) - 1023L;
  f = ( _x.u & ( (1ULL<<52) - 1 ) );

  if ( 1024==e && 0!=f ) /* NaN */
  {
    /* Return a quiet (!) NaN with original payload bits. */
    return (float16_t)( (s<<15) | 0x7e00 | (int32_t)(f>>(52-10)) ) ; 
  }
  else if ( e>15 ) /* Too large magnitude */
  {
    return (float16_t)( (s<<15) | 0x7c00 );
  }
  else if ( e<-25 ) /* Too small magnitude */
  {
    return (float16_t)(s<<15);
  }
  else if ( -14>e ) /* Subnormal half precision value */
  {
    /* Restore the implicit leading significand bit. */ 
    f += (1ULL<<52);
    /* 10 significand trailing bits */
    tb = (int32_t)(f>>(52-10-e-14)); 
    /* Rounding bit value */
    rbit = (int32_t)( (f>>(52-10-1-e-14)) & 1L ); 
    /* Are there any nonzero bits left? */
    tail = ( 0 != ( f & ( (1ULL<<(52-10-1-e-14)) - 1 ) ) );
    e = -15;
  }
  else /* Normal half precision value */
  {
    /* 10 significand trailing bits */
    tb = (int32_t)(f>>(52-10)); 
    /* Rounding bit value */
    rbit = (int32_t)( (f>>(52-10-1)) & 1L ); 
    /* Are there any nonzero bits left? */
    tail = ( 0 != ( f & ( (1ULL<<(52-10-1)) - 1 ) ) );
  }

  if ( 0!=rbit ) 
  {
    /* Round to nearest. In case of a tie round to even mantissa. */
    tb += ( 0==tail ? (tb&1) : 1 );
    /* Look if the exponent shall be adjusted for the rounded mantissa. */
    if ( tb & (1L<<10) )
    { 
      /* Correct the exponent and mantissa, and check if the rounded value overflows. */
      tb = ( ++e<=15 ? ( (tb+1024) >> 1 ) : 0 );
    }
  }

  return (float16_t)( (s<<15) | ( (e+15) << 10 ) | ( tb & 0x3ff ) );

} /* conv_f64_to_f16() */

/* Convert half precision floating-point data to single precision format. */
float32_t conv_f16_to_f32( float16_t x )
{
  int32_t s,e,f;
  union ufloat32uint32 y = {0};

  s = ( (uint16_t)x >> 15 );
  e = ( (uint16_t)x >> 10 ) & 31;
  f = ( (uint16_t)x >>  0 ) & 1023;

  if ( 31==e )
  {
    if ( 0==f ) /* +/-Infinity */
    {
      y.u = ( (uint32_t)s << 31 ) | ( 0xffUL << 23 );
    }
    else /* NaN */
    {
      /* Return a quiet (!) NaN with original payload bits. */
      y.u = ( (uint32_t)s << 31 ) | ( 0x1ffUL << 22 ) | ( ( (uint32_t)f & 0x1ffUL ) << (23-10) ); 
    }
  }
  else if ( 0==e && 0!=f ) /* Subnormal */
  {
    int nsa = _S_exp0_l( f );
    /* Clear the leading bit of the significand. */
    f = (f<<(nsa-7)) ^ (1L<<23);
    /* Q(23-e) <- [ Q(10+14) + nsa ] + (23-30) */
    y.u = ( (uint32_t)s << 31 ) | ( (uint32_t)(6-nsa+127) << 23 ) | (uint32_t)f;
  }
  else if ( 0==e ) /* +/-0 */
  {
    y.u = ( (uint32_t)s << 31 );
  }
  else /* Normal finite value. */
  {
    y.u = ( (uint32_t)s << 31 ) | ( (uint32_t)(e-15+127) << 23 ) | ( (uint32_t)f << (23-10) );
  }

  return (y.f);

} /* conv_f16_to_f32() */

/* Convert half precision floating-point data to double precision format. */
float64_t conv_f16_to_f64( float16_t x )
{
  int32_t s,e,f;
  union ufloat64uint64 y = {0};

  s = ( (uint16_t)x >> 15 );
  e = ( (uint16_t)x >> 10 ) & 31;
  f = ( (uint16_t)x >>  0 ) & 1023;

  if ( 31==e )
  {
    if ( 0==f ) /* +/-Infinity */
    {
      y.u = ( (uint64_t)s << 63 ) | ( 0x7ffULL << 52 );
    }
    else /* NaN */
    {
      /* Return a quiet (!) NaN with original payload bits. */
      y.u = ( (uint64_t)s << 63 ) | ( 0xfffULL << 51 ) | ( ( (uint64_t)f & 0x1ffULL ) << (52-10) ); 
    }
  }
  else if ( 0==e && 0!=f ) /* Subnormal */
  {
    int nsa = _S_exp0_l( f );
    /* Clear the leading bit of the significand. */
    f = (f<<nsa) & ~(1L<<30);
    /* Q(52-e) <- [ Q(10+14) + nsa ] + (52-30) */
    y.u = ( (uint64_t)s << 63 ) | ( (uint64_t)(6-nsa+1023) << 52 ) | ( (uint64_t)f << (52-30) );
  }
  else if ( 0==e ) /* +/-0 */
  {
    y.u = ( (uint64_t)s << 63 );
  }
  else /* Normal finite value. */
  {
    y.u = ( (uint64_t)s << 63 ) | ( (uint64_t)(e-15+1023) << 52 ) | ( (uint64_t)f << (52-10) );
  }

  return (y.f);

} /* conv_f16_to_f64() */

/* Scale a 16-bit signed integer value and round it to the nearest half precision
 * floating-point value, with halfway cases rounded to even mantissa. Scaling factor
 * equals 2 to negative power of the input argument e (bi-directional). */
float16_t conv_i16_to_f16( int16_t x, int p )
{
    int32_t s, e, f, nsa, sft, mask, rt;
    /* Remember the sign of the input integer. */
    s = x<0;
    /* Normalize and implicitly scale the input integer: 
     * Q(nsa+p) = 2^-p*Q0 + nsa + p */
    nsa = _S_exp0_l(x); f = (int32_t)x<<nsa;
    /* Take the absolute value. */
    if (s) f = -f;
    /* Check if the integer is a power of two and negative. */
    if (f<0) { f = (uint32_t)f>>1; nsa--; }
    /* Compute the floating-point exponent, given that the leading 1 
     * resides is in 30th bit position of the absolute value. Also,
     * take care of a zero integer on input. */
    e = f>0 ? 30-nsa-p : -15;
    /* Compute the right shift amount such that when applied to the 31-bit significand:
     *  A) For a normal number, the leading one moves into position 10, and
     *     the trailing bits occupy positions 0..9
     *  B) For a subnormal number, the leading one also moves into the appropriate
     *     position between 0 and 9.
     *  C) For a tiny non-zero number, bits 0..10 are flushed to zero. */
    sft = 20 + MIN(12, MAX(-14, e)-e); /* 20<=sft<=32 */
    /* Rounding term and round-off mask. */
    rt = 1L<<(sft-1); mask = (rt<<1)-1;
    /* Round to nearest representable number, ties are rounded to even mantissa. */
    if ((f&mask)==rt) {
        /* A tie, round to even ULP. */
        f = (f+rt) & ~(rt<<1);
    } else {
        /* Not a tie, round to nearest. */
        f += rt;
    }
    /* Check if the number has been rounded to a power of 2. Note that the nested
     * conditionals have an intersection for e==-15: there exist numbers that are
     * subnormal before the rounding, and become normal after the rounding, for 
     * example +/-2047*2^-25. */
    if (f<0) { 
        /* Move the leading one back to position 30. */
        f = (uint32_t)f>>1;
        if (e<-14) {
            /* If the number is subnormal before the rounding, then the leading
             * one should appear one position to the left of what we have initially
             * supposed. */
            sft--;
        }
        if (e>-16) {
            /* If the number is normal after the rounding, adjust the exponent to
             * match the increased magnitude. */
            e++;
        }
    }
    if (e<16) {
        /* For a rounded value of magnitude less or equal to the maximim representable
         * number, form the half precision significand and cut off the leading one, 
         * unless the value is subnormal. */
        f = (f>>sft)&1023;
    } else {
        /* For a huge value, zero out the significand. */
        f = 0;
    }
    /* Bias the exponent and clamp it to the conventional range. */
    e = MIN(31, MAX(0, e+15));
    /* Form the half precision representation. */
    return (float16_t)((s<<15)|(e<<10)|f);
} /* conv_i16_to_f16() */

/* Sum of two half precision floating-point values */
float16_t add_f16( float16_t x, float16_t y )
{
  return ( conv_f64_to_f16( conv_f16_to_f64(x) +
                            conv_f16_to_f64(y) ) );
}

/* Sum of two half precision floating-point values */
float16_t sub_f16( float16_t x, float16_t y )
{
  return ( conv_f64_to_f16( conv_f16_to_f64(x) -
                            conv_f16_to_f64(y) ) );
}

/* Multiply two half precision floating-point values */
float16_t mul_f16( float16_t x, float16_t y )
{
  return ( conv_f64_to_f16( conv_f16_to_f64(x) *
                            conv_f16_to_f64(y) ) );
}

/* Multiply two half precision floating-point values */
float16_t div_f16( float16_t x, float16_t y )
{
  return ( conv_f64_to_f16( conv_f16_to_f64(x) /
                            conv_f16_to_f64(y) ) );
}

/* Fused multiply-add: return x*y+z. */
float16_t fma_f16( float16_t x, float16_t y, float16_t z )
{
  float64_t xd, yd, zd;

  xd = conv_f16_to_f64(x);
  yd = conv_f16_to_f64(y);
  zd = conv_f16_to_f64(z);

  return ( conv_f64_to_f16( xd*yd+zd ) );
}

/* Fused multiply-subtract: return z-x*y. */
float16_t fms_f16( float16_t x, float16_t y, float16_t z )
{
  float64_t xd, yd, zd;

  xd = conv_f16_to_f64(x);
  yd = conv_f16_to_f64(y);
  zd = conv_f16_to_f64(z);

  return ( conv_f64_to_f16( zd-xd*yd ) );
}

/* Return absolute value of half precision floating-point number. */
float16_t fabs_f16( float16_t x )
{
  return ( isnan_f16(x) ? x : (float16_t)( (uint16_t)x & 0x7fff ) );
}

/* Multiply the first argument by two, raised to the power of the second argument. */
float16_t ldexp_f16( float16_t x, int n )
{
  int n0, n1, n2;
  float16_t y, s0, s1, s2;

  /* Sensible range of exponent adjustment is [-41,40], as shown below:
   *  - max finite -> zero: pow2(it_half((2-2^-10)*2^15),-41) == 0
   *  - min subnormal -> Infinity: pow2(it_half(2^-24),40) == Infinity */
  n = MAX( -41, MIN( 40, n ) );

  /* Multiplication by a power of 2 is the most convenient way to modify
   * the exponent. Maximum representable power of two for the half precision
   * is 15, thus we have to split n into 3 components to cover the full range
   * of exponent adjustment: [-41,40]. The value is partitioned in such a way
   * that a subnormal result may appear only once in three products, excepting
   * trivial multiplucations by 1. This guarantees that the subnormal result is
   * rounded only once. */
  n2 = MAX( -14, MIN( 15, n    ) ); s2 = (float16_t)( (uint16_t)(n2+15)<<10 );
  n1 = MAX( -14, MIN( 15, n-n2 ) ); s1 = (float16_t)( (uint16_t)(n1+15)<<10 );
  n0 = n-n2-n1;                     s0 = (float16_t)( (uint16_t)(n0+15)<<10 );

  /* 1. The order of multiplications is crucial.
   * 2. We rely on the * operator to raise the FE_OVERFLOW floating-point
   *    exception whenever result overflows. */
  y = mul_f16( mul_f16( mul_f16(x,s0), s1 ), s2 );

  return (y);

} /* ldexp_f16() */

/* Round the argument to the nearest integer value in half precision floating-
 * point format, rounding halfway cases away from zero. */
float16_t round_f16( float16_t x )
{
  int32_t e;
  e = ((uint16_t)x>>10) & 31;
  if (e>24)
  {
    return x;
  }
  else if (e<14)
  {
    return (float16_t)(x & (1<<15));
  }
  else
  {
    int32_t h = 1<<(24-e);
    return (float16_t)((x+h) & ~(MIN(1024, 2*h)-1));
  }
} /* round_f16() */

/* Round the argument to the nearest integer value, rounding halfway cases 
 * away from zero. For infinite input argument, the function returns either
 * INT_MIN (for negative infinity on input) or INT_MAX (for positive infinity).
 * If the input argument is Not-a-Number (NaN), the returned value is not defined. */
int iround_f16( float16_t x )
{
  int32_t s, e, f, t=0;
  s = ((uint16_t)x>>15) & 1;
  e = ((uint16_t)x>>10) & 31;
  f = (uint16_t)x & 1023;
  if (e>30)
  {
    return s ? INT_MIN : INT_MAX;
  }
  else if (e>24)
  {
    t = (1024+f)<<(e-25);
  }
  else if (e>13)
  {
    int32_t h = 1<<(24-e);
    t = (1024+f+h)>>(25-e);
  }
  return s ? -t : t;
} /* iround_f16() */

/* Round the argument to integer value, nearest to but no larger in magnitude 
 * than the argument. For infinite input argument, the function returns either
 * INT_MIN (for negative infinity on input) or INT_MAX (for positive infinity).
 * If the input argument is Not-A-Number (NaN), the returned value is not defined. */
int itrunc_f16( float16_t x )
{
  int32_t s, e, f, t=0;
  s = ((uint16_t)x>>15) & 1;
  e = ((uint16_t)x>>10) & 31;
  f = (uint16_t)x & 1023;
  if (e>30)
  {
    return s ? INT_MIN : INT_MAX;
  }
  else if (e>24)
  {
    t = (1024+f)<<(e-25);
  }
  else if (e>13)
  {
    t = (1024+f)>>(25-e);
  }
  return s ? -t : t;
} /* itrunc_f16() */

/* Break a floating-point number x into a normalized fraction f, 0.5<=f<1, and an
 * integral power of 2 e, such that x == f*2^e. For a numeric input x, return the
 * fraction f and set output argument *pe to the integral power e. For a non-numeric
 * input x, just return x. */
float16_t frexp_f16(float16_t x, int * pe)
{
  int32_t s, e, f, t;
  NASSERT(pe);
  s = ((uint16_t)x>>15) & 1;
  e = ((uint16_t)x>>10) & 31;
  f = (uint16_t)x & 1023;
  if (e==31) { /* Non-numeric value on input (+/-Inf, NaN) */
      return x;
  } else if (e>0) { /* Normal number on input */
      *pe = e-14;
      return (float16_t)((s<<15)|(14<<10)|f);
  } else if (f>0) { /* Subnormal number on input */
      t = 30-_S_exp0_l(f);
      *pe = t-23;
      return (float16_t)((s<<15)|(14<<10)|((f<<(10-t))&1023));
  } else { /* Zero on input */
      *pe = 0;
      return (float16_t)(s<<15);
  }
} /* frexp_f16() */

typedef enum { _lt_f16, _lte_f16, _gt_f16, 
               _gte_f16, _eq_f16, _oneq_f16 } f16_order_t;

static int compare_f16( float16_t a, float16_t b, f16_order_t ord )
{
  int r;
  int32_t ai, bi;
  if ( un_f16(a,b) ) return (0);
  ai = ( (uint16_t)a & 0x7fff ); if ( (uint16_t)a & 0x8000 ) ai = -ai;
  bi = ( (uint16_t)b & 0x7fff ); if ( (uint16_t)b & 0x8000 ) bi = -bi;
  r = ( ( ord == _lt_f16 ) ? (ai<bi ) : ( ord == _lte_f16  ) ? (ai<=bi) :
        ( ord == _gt_f16 ) ? (ai>bi ) : ( ord == _gte_f16  ) ? (ai>=bi) :
        ( ord == _eq_f16 ) ? (ai==bi) : ( ord == _oneq_f16 ) ? (ai!=bi) : -1 );
  NASSERT( r>=0 );
  return (r);
}

int lt_f16  ( float16_t a, float16_t b ) { return compare_f16( a, b, _lt_f16   ); } /* a<b  */
int lte_f16 ( float16_t a, float16_t b ) { return compare_f16( a, b, _lte_f16  ); } /* a<b  */
int gt_f16  ( float16_t a, float16_t b ) { return compare_f16( a, b, _gt_f16   ); } /* a<b  */
int gte_f16 ( float16_t a, float16_t b ) { return compare_f16( a, b, _gte_f16  ); } /* a<b  */
int eq_f16  ( float16_t a, float16_t b ) { return compare_f16( a, b, _eq_f16   ); } /* a<b  */
int oneq_f16( float16_t a, float16_t b ) { return compare_f16( a, b, _oneq_f16 ); } /* a!=b */
