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
  RMS calculation routines for 16-bit and 32-bit data.

  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "types.h"
/* Fixed point arithmetic. */
#include "NatureDSP_Math.h"
/* RMS routines declaration. */
#include "rms.h"

#define MIN(a,b)  ( (a)<(b) ? (a) : (b) )

/* value for conversion from ln to log2 */
#define LOG2_CONST (0x38aa)   // 1-1/ln(2)

static int16_t _S_log2_l( int32_t accA );

/*-------------------------------------------------------------------------
Calculate relative Root-Mean-Square (RMS) power of a real 16-bit or 32-bit
signal with a full-scale sine wave used as a reference. Result is a Q8.7
value measured in dB.

NOTE
  Functions may be also used to obtain the relative RMS power of a complex
  waveform against the full-scale complex wave.

Parameters:
  Input:
    x[N]          Real input signal
  Returned value:
                  Relative RMS power in dB, Q8.7
-------------------------------------------------------------------------*/

int16_t rms16( const int16_t * x, int N )
{
  int32_t acc0, acc1;
  int16_t y, z, rms;

  int exp, n, m;

  int16_t N_exp = S_exp0_l( N );

  ASSERT( x && N > 0 );

  for ( exp=15, n=0; n<N; n++ )
  {
    exp = MIN( exp, S_exp0_s( x[n] ) );
  }

  for ( acc0=0, n=0; n<N/4; n++ )
  {
    for ( acc1=0, m=0; m<4; m++ )
    {
      // Q(13+exp) <- Q15 + exp - 2 w/ rounding
      y = ( S_sature_l( (int32_t)( x[4*n+m]<<exp ) + 1L ) >> 2 );
      // Q(26+2*exp) <- Q(13+exp)*Q(13+exp)
      acc1 += L_mul_ss( y, y );
    }

    // Q(2*exp+N_exp) <- Q(26+2*exp) - (30-N_exp) + 4
    acc0 += L_shl_l( acc1, N_exp - 26 );
  }

  for ( acc1=0, m=0; m<(N&3); m++ )
  {
    // Q(13+exp) <- Q15 + exp - 2 w/ rounding
    y = ( S_sature_l( (int32_t)( x[4*n+m]<<exp ) + 1L ) >> 2 );
    // Q(26+2*exp) <- Q(13+exp)*Q(13+exp)
    acc1 += L_mul_ss( y, y );
  }

  // Q(2*exp+N_exp) <- Q(26+2*exp) - (30-N_exp) + 4
  acc0 += L_shl_l( acc1, N_exp - 26 );

  //
  // rms = 10*log10(acc0/N/0.5) - 10/log2(10)*( 2*exp + N_exp )
  //

  if ( acc0 > 0 )
  {
    // z = log10(acc0/N/0.5) = ( log2(acc0) - log2(N) + 1 )*1/log2(10);
    // Q4.11 <- Q6.9*Q17 - 15 w/o rounding
    z = (int16_t)( (int32_t)( _S_log2_l( acc0 ) - 
                              _S_log2_l( N    ) + (1<<9) )*39457 >> 15 );
  }
  else
  {
    z = (int16_t)0x8000;
  }

  // Q8.7 <- sat16( ( Q0*Q4.11 - Q2.11*Q0 ) - 4 w/ rounding )
  rms = S_sature_l( ( L_mul_ss( 10, z ) - 
                      L_mul_ss( 6165, 2*exp+N_exp ) + (1<<3) ) >> 4 );

  return (rms);

} // rms16()

int16_t rms32( const int32_t * x, int N )
{
  int32_t acc0, acc1;
  int16_t y, z, rms;

  int exp, n, m;

  int16_t N_exp = S_exp0_l( N );

  ASSERT( x && N > 0 );

  for ( exp=30, n=0; n<N; n++ )
  {
    exp = MIN( exp, S_exp0_l( x[n] ) );
  }

  for ( acc0=0, n=0; n<N/4; n++ )
  {
    for ( acc1=0, m=0; m<4; m++ )
    {
      // Q(13+exp) <- Q31 + exp - 18 w/ rounding
      y = S_extract_l( L_add_ll( x[4*n+m]<<exp, 1L<<17 ) >> 2 );
      // Q(26+2*exp) <- Q(13+exp)*Q(13+exp)
      acc1 += L_mul_ss( y, y );
    }

    // Q(2*exp+N_exp) <- Q(26+2*exp) - (30-N_exp) + 4
    acc0 += L_shl_l( acc1, N_exp - 26 );
  }

  for ( acc1=0, m=0; m<(N&3); m++ )
  {
    // Q(13+exp) <- Q31 + exp - 18 w/ rounding
    y = S_extract_l( L_add_ll( x[4*n+m]<<exp, 1L<<17 ) >> 2 );
    // Q(26+2*exp) <- Q(13+exp)*Q(13+exp)
    acc1 += L_mul_ss( y, y );
  }

  // Q(2*exp+N_exp) <- Q(26+2*exp) - (30-N_exp) + 4
  acc0 += L_shl_l( acc1, N_exp - 26 );

  //
  // rms = 10*log10(acc0/N/0.5) - 10/log2(10)*( 2*exp + N_exp )
  //

  if ( acc0 > 0 )
  {
    // z = log10(acc0/N/0.5) = ( log2(acc0) - log2(N) + 1 )*1/log2(10);
    // Q4.11 <- Q6.9*Q17 - 15 w/o rounding
    z = (int16_t)( (int32_t)( _S_log2_l( acc0 ) - 
                              _S_log2_l( N    ) + (1<<9) )*39457 >> 15 );
  }
  else
  {
    z = (int16_t)0x8000;
  }

  // Q8.7 <- sat16( ( Q0*Q4.11 - Q2.11*Q0 ) - 4 w/ rounding )
  rms = S_sature_l( ( L_mul_ss( 10, z ) - 
                      L_mul_ss( 6165, 2*exp+N_exp ) + (1<<3) ) >> 4 );

  return (rms);

} // rms32()

/*
    Algorithm: ln(1+x) ~ sum(a * T'(x))
    where T'(x) -- shifted Chebychev polynomials of first kind
    a(0) =  0.37645281291919543163
    a(1) =  0.34314575050761980479
    a(2) = -0.02943725152285941438
    a(3) =  0.00336708925556438925
    a(4) = -0.00043327588861004446

    MATLAB code:
    function y=myLog(x)
    x=x-1;
    a=[0.37645281291919543163 0.34314575050761980479 -0.02943725152285941438 0.00336708925556438925 -0.00043327588861004446];
    y = 2*x-1;
    T0 = 1;
    T1 = y;
    T2 = 2*y.*T1 - T0;
    T3 = 2*y.*T2 - T1;
    T4 = 2*y.*T3 - T2;
    y=a(1)*T0 + a(2)*T1 + a(3)*T2 + a(4)*T3 + a(5)*T4;;

    Accuracy: 1.1e-3

    References:
    Y.L Luk, Mathematical functions and their approximations, Academic Press, 1975
*/
//--------------------------------------------------
//    Extracted from NatureDSP_Math library
//    Description:
//    calculates estimation of log2(x)
//
//    x - Q0 argument
//
//    returns 0x8000 on negative of zero argument
//    or number in Q9 format
//--------------------------------------------------
int16_t _S_log2_l (int32_t accA)
{
  int16_t T1,T0,w16Res,w16Exp,y;
  int32_t accB;

  if (accA<=0) return (int16_t)0x8000;

  w16Exp = S_exp_l (accA);
  w16Res = 30 - w16Exp;
  accA ^= ((int32_t)1) << w16Res;
  accA <<= (w16Exp+1);//ln(1 + x) = x - x*x/2 + x*x*x/4

//  accA <<= (w16Exp);
  accA -= 0x40000000UL;
  y = (int16_t)(accA >> 15);        // 2*x-1

  T0 = 0x7FFF;
  T1 = y;

//  accB = 808426260L + L_mpy_ss(T1,11244);
  accB = L_mac_ss(808426260L,T1,11244);
  accA  = L_mpy_ss(y,(int16_t)(T1>>1))>>14;
  accA -= T0;
  T0 = S_sature_l(accA);
  accB += L_mpy_ss(T0, -30867)>>5;

  accA  = L_mpy_ss(y,(int16_t)(T0>>1))>>14;
  accA -= T1;
  T1 = S_sature_l(accA);
  accA = L_mpy_ss(y, T1);
  accB += L_mpy_ss(T1,28245)>>8;

  accA  = L_mpy_ss(y,(int16_t)(T1>>1))>>14;
  accA -= T0;
  T0 = S_sature_l(accA);
  accB += L_mpy_ss(T0, -3635)>>8;

    //ln(1+x) -> log2(1+x)
  accA  =  L_mpy_ss(S_extract_l(accB), LOG2_CONST)>>1;
  accA +=  accB>>1;
  accA += (1L<<20);
  accA >>= 21;

  accA += w16Res << 9;

  return (int16_t)accA; //|sign-1bit|integer part-6bits|fractional part-9bits|
}
