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
    NatureDSP Signal Processing Library. IIR part
    Biquad Real Block IIR, 32x16-bit, Direct Form II
    C code optimized for HiFi4
    IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_iir.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Filter instance structure. */
#include "bqriir32x16_df2_common.h"
#define sz_i16 sizeof(int16_t)
#define sz_i32 sizeof(int32_t)

/*-------------------------------------------------------------------------
  Bi-quad Real Block IIR
  Computes a real IIR filter (cascaded IIR direct form I or II using 5 
  coefficients per bi-quad + gain term). Real input data are stored
  in vector x. The filter output result is stored in vector r. The filter 
  calculates N output samples using SOS and G matrices.
  NOTE:
  1. Bi-quad coefficients may be derived from standard SOS and G matrices 
  generated by MATLAB. However, typically biquad stages have big peaks 
  in their step response which may cause undesirable overflows at the 
  intermediate outputs. To avoid that the additional scale factors 
  coef_g[M] may be applied. These per-section scale factors may require 
  some tuning to find a compromise between quantization noise and possible 
  overflows. Output of the last section is directed to an additional 
  multiplier, with the gain factor being a power of two, either negative 
  or non-negative. It is specified through the total gain shift amount 
  parameter gain of each filter initialization function.
  2. 16x16 filters may suffer more from accumulation of the roundoff errors,
  so filters should be properly designed to match noise requirements
  3. Due to the performance reasons, IIR biquad filters may introduce 
  additional algorithmic delay of several sampless. Amount of that delay
  might be requested by the  xxx_groupDelay API. For sensitive applications
  all the filters have delayless implementations (with  _nd  suffix in the name).
  Formally, the xxx_groupDelay APIs is also implemented for that kind of filters,
  but return zero.
  
  Precision: 
  16x16         16-bit data, 16-bit coefficients, 16-bit intermediate 
                stage outputs (DF1, DF1 stereo, DF II form)
  32x16         32-bit data, 16-bit coefficients, 32-bit intermediate 
                (DF1, DF1 stereo, DF II form) stage outputs
  32x32         32-bit data, 32-bit coefficients, 32-bit intermediate 
                (DF I, DF1 stereo,  DF II form) stage outputs 
  f             floating point (DF I, DF1 stereo, DF II and DF IIt)

  Input:
  N             length of input sample block
  M             number of bi-quad sections
  S             1 for mono, 2 for stereo API
  s[]           scratch memory area (for fixed-point functions only). 
                Minimum number of bytes depends on selected filter structure 
                and precision:
                  BQRIIR16X16_DF1_SCRATCH_SIZE( N,  M ), or
                  BQRIIR16X16_DF2_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X16_DF1_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X16_DF2_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X32_DF1_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X32_DF2_SCRATCH_SIZE( N,  M ),
                  STEREO_BQRIIR16X16_DF1_SCRATCH_SIZE( N, M ) , or
                  STEREO_BQRIIR32X32_DF1_SCRATCH_SIZE( N, M ) , or
                  STEREO_BQRIIRF_DF1_SCRATCH_SIZE    ( N, M )
                 If a particular macro returns zero, then the corresponding
                 IIR doesn't require a scratch area and parameter s may 
                 hold zero.
  coef_sos[M*5]  filter coefficients stored in blocks of 5 numbers: 
                 b0 b1 b2 a1 a2. 
                 For fixed-point funcions, fixed point format of filter 
                 coefficients is Q1.14 for 16x16 and 32x16, or Q1.30 for 32x32 
  coef_sosl[M*5] filter coefficients for the left channel (stereo filters only)
  coef_sosr[M*5] filter coefficients for the left channel (stereo filters only)
  coef_g[M]      scale factor for each section, Q15 (for fixed-point 
                 functions only)
  coef_gl[M]     scale factor for the left channel (stereo filters only)
  coef_gr[M]     scale factor for the right channel (stereo filters only)
  gain           total gain shift amount applied to output signal of the
                 last section, -48..15
  gainl          total gain shift amount  for the left channel (stereo filters 
                 only)
  gainr          total gain shift amount for the right channel (stereo filters 
                 only)

  x[N*S]         input samples, Q15, Q31 or floating point. Stereo samples 
                 go in interleaved order (left, right)
  Output:
  r[N*S]         output data, Q15, Q31 or floating point. Stereo samples go 
                 in interleaved order (left, right) 

  Restriction:
  x,r,s,coef_g,coef_sos should not overlap
  N   - must be a multiple of 2
  s[] - whenever supplied must be aligned on an 16-bytes boundary
-------------------------------------------------------------------------*/
void bqriir32x16_df2( bqriir32x16_df2_handle_t _bqriir,
                      void    * restrict       s,
                      int32_t * restrict       r,
                const int32_t *                x, int N )
{
#if (ALG32X16_DF2_ND == 1)
    return bqriir32x16_df2_nd(_bqriir, s, r, x, N);
#else
  bqriir32x16_df2_ptr_t bqriir = (bqriir32x16_df2_ptr_t)_bqriir;

  const ae_int16x8 * restrict coef;
        ae_int32x4 * restrict state;
  const ae_int32x2 * restrict px;
        ae_int32x2 * restrict pr;
  ae_valign alx, alr;

  int M;
  int m, n;

  NASSERT( bqriir && bqriir->magic == MAGIC && r && x );
  if (N<=0) return;
  NASSERT( N%2==0 );

  M = bqriir->M;
  alr = AE_ZALIGN64();

  coef  = (const ae_int16x8*)bqriir->coef;
  state = (      ae_int32x4*)bqriir->state;

  px = (const ae_int32x2*)x;
  pr = (      ae_int32x2*)r;

  //
  // Perform data block processing for each pair of successive sections. Use the
  // output array r[N] for temporal storage of inter-section signal.
  //

  for ( m=0; m<(M/4); m++ )
  {
    ae_int16x4 cf0_a2a1Zg, cf0_a2a1gZ, cf0_b2b1b0Z, cf0_Zb2b1b0;
    ae_int16x4 cf1_a2a1Zg, cf1_a2a1gZ, cf1_b2b1b0Z, cf1_Zb2b1b0;
    ae_int16x4 cf2_a2a1Zg, cf2_a2a1gZ, cf2_b2b1b0Z, cf2_Zb2b1b0;
    ae_int16x4 cf3_a2a1Zg, cf3_a2a1gZ, cf3_b2b1b0Z, cf3_Zb2b1b0;

    ae_int32x2 x0, x1, z01, u01, v01, y01;

    ae_int64 fb0, ff0; // fbX == Feedback, ffX == Feedforward
    ae_int64 fb1, ff1; // X - sample number
    ae_int32x2 d2d1_0, dxx_0, d2d1_1, dxx_1;
    ae_int32x2 d2d1_2, dxx_2, d2d1_3, dxx_3;

    //
    // Load coefficients for the m-th pair of sections.
    //

    // Q14
    AE_L16X4X2_IP(cf0_a2a1Zg, cf0_b2b1b0Z, coef, sz_i16*4*2);
    AE_L16X4X2_IP(cf1_a2a1Zg, cf1_b2b1b0Z, coef, sz_i16*4*2);
    cf0_a2a1gZ = AE_SEL16_7632(cf0_a2a1Zg, AE_SHORTSWAP(cf0_a2a1Zg));
    cf1_a2a1gZ = AE_SEL16_7632(cf1_a2a1Zg, AE_SHORTSWAP(cf1_a2a1Zg));
    cf0_Zb2b1b0 = AE_SEL16_4321(cf0_b2b1b0Z, cf0_b2b1b0Z);
    cf1_Zb2b1b0 = AE_SEL16_4321(cf1_b2b1b0Z, cf1_b2b1b0Z);
    // Q14
    AE_L16X4X2_IP(cf2_a2a1Zg, cf2_b2b1b0Z, coef, sz_i16*4*2);
    AE_L16X4X2_IP(cf3_a2a1Zg, cf3_b2b1b0Z, coef, sz_i16*4*2);
    cf2_a2a1gZ = AE_SEL16_7632(cf2_a2a1Zg, AE_SHORTSWAP(cf2_a2a1Zg));
    cf3_a2a1gZ = AE_SEL16_7632(cf3_a2a1Zg, AE_SHORTSWAP(cf3_a2a1Zg));
    cf2_Zb2b1b0 = AE_SEL16_4321(cf2_b2b1b0Z, cf2_b2b1b0Z);
    cf3_Zb2b1b0 = AE_SEL16_4321(cf3_b2b1b0Z, cf3_b2b1b0Z);

    //
    // Load sections' state.
    //

    // Q31
    AE_L32X2X2_I(d2d1_0, d2d1_1, state, 0*sz_i32*2*2);
    AE_L32X2X2_I(d2d1_2, d2d1_3, state, 1*sz_i32*2*2);

    //
    // Pass N/2 sample pairs through 2 sections.
    //

    alx = AE_LA64_PP(px);
    __Pragma( "loop_count min=1" );
    for ( n=0; n<N/2; n++ )
    {
      // Q31
      AE_L32_IP(x0, castxcc(ae_int32,px), sizeof(int32_t));
      AE_L32_IP(x1, castxcc(ae_int32,px), sizeof(int32_t));

      //------------------------------------------------------------------
      // Section 0; samples 0 and 1.
      dxx_0 = d2d1_0;
      // Q45 <- Q31*Q14
      fb0 = AE_MULZAAAAQ32X16( d2d1_0, x0, cf0_a2a1Zg );
      // Q31 <- Q45 + 2 - 16 w/ saturation
      AE_PKSR32(d2d1_0, fb0, 2); // update delay elements
      fb1 = AE_MULZAAAAQ32X16( d2d1_0, x1, cf0_a2a1Zg );
      AE_PKSR32(d2d1_0, fb1, 2); // update delay elements

      // Q45 <- Q31*Q14
      AE_MULZAAAA2Q32X16(ff0, ff1, dxx_0, d2d1_0, cf0_b2b1b0Z, cf0_Zb2b1b0);
      z01 = AE_TRUNCA32X2F64S(ff0, ff1, 18);

      //------------------------------------------------------------------
      // Section 1; samples 0 and 1.
      dxx_1 = d2d1_1;
      // Q45 <- Q31*Q14
      fb0 = AE_MULZAAAAQ32X16( d2d1_1, z01, cf1_a2a1gZ );
      // Q31 <- Q45 + 2 - 16 w/ saturation
      AE_PKSR32(d2d1_1, fb0, 2); // update delay elements
      fb1 = AE_MULZAAAAQ32X16( d2d1_1, z01, cf1_a2a1Zg );
      AE_PKSR32(d2d1_1, fb1, 2); // update delay elements

      // Q45 <- Q31*Q14
      AE_MULZAAAA2Q32X16(ff0, ff1, dxx_1, d2d1_1, cf1_b2b1b0Z, cf1_Zb2b1b0);
      u01 = AE_TRUNCA32X2F64S(ff0, ff1, 18);

      //------------------------------------------------------------------
      // Section 2; samples 0 and 1.
      dxx_2 = d2d1_2;
      // Q45 <- Q31*Q14
      fb0 = AE_MULZAAAAQ32X16( d2d1_2, u01, cf2_a2a1gZ );
      // Q31 <- Q45 + 2 - 16 w/ saturation
      AE_PKSR32(d2d1_2, fb0, 2); // update delay elements
      fb1 = AE_MULZAAAAQ32X16( d2d1_2, u01, cf2_a2a1Zg );
      AE_PKSR32(d2d1_2, fb1, 2); // update delay elements

      // Q45 <- Q31*Q14
      AE_MULZAAAA2Q32X16(ff0, ff1, dxx_2, d2d1_2, cf2_b2b1b0Z, cf2_Zb2b1b0);
      v01 = AE_TRUNCA32X2F64S(ff0, ff1, 18);

      //------------------------------------------------------------------
      // Section 3; samples 0 and 1.
      dxx_3 = d2d1_3;
      // Q45 <- Q31*Q14
      fb0 = AE_MULZAAAAQ32X16( d2d1_3, v01, cf3_a2a1gZ );
      // Q31 <- Q45 + 2 - 16 w/ saturation
      AE_PKSR32(d2d1_3, fb0, 2); // update delay elements
      fb1 = AE_MULZAAAAQ32X16( d2d1_3, v01, cf3_a2a1Zg );
      AE_PKSR32(d2d1_3, fb1, 2); // update delay elements

      // Q45 <- Q31*Q14
      AE_MULZAAAA2Q32X16(ff0, ff1, dxx_3, d2d1_3, cf3_b2b1b0Z, cf3_Zb2b1b0);
      y01 = AE_TRUNCA32X2F64S(ff0, ff1, 18);

      // Q31
      AE_SA32X2_IP( y01, alr, pr );
    }
    AE_SA64POS_FP(alr, pr);

    //
    // Save sections' state.
    //

    // Q31
    AE_S32X2X2_IP(d2d1_0, d2d1_1, state, sz_i32*2*2);
    AE_S32X2X2_IP(d2d1_2, d2d1_3, state, sz_i32*2*2);

    // Second to last pair of sections are fed with output signal of the
    // previous pair.
    px = pr = (ae_int32x2*)(r);
  }

  //
  // Process remaining sections (last 0..3)
  //
  M = M & 3;

  if (M == 3)
  {
    ae_int16x4 cf0_a2a1Zg, cf0_a2a1gZ, cf0_b2b1b0Z, cf0_Zb2b1b0;
    ae_int16x4 cf1_a2a1Zg, cf1_a2a1gZ, cf1_b2b1b0Z, cf1_Zb2b1b0;
    ae_int16x4 cf2_a2a1Zg, cf2_a2a1gZ, cf2_b2b1b0Z, cf2_Zb2b1b0;
    ae_int32x2 x01, z01, y01;
    ae_int64 fb00, fb01, ff00, ff01; // fbXY == Feedback, ffXY == Feedforward
    ae_int64 fb10, fb11, ff10, ff11; // X - section number, Y - sample number
    ae_int64 fb20, fb21, ff20, ff21;
    ae_int32x2 d2d1_0, dxx_0, d2d1_1, dxx_1, d2d1_2, dxx_2;

    //
    // Load coefficients for the m-th pair of sections.
    //

    // Q14
    AE_L16X4X2_IP(cf0_a2a1Zg, cf0_b2b1b0Z, coef, sz_i16*4*2);
    AE_L16X4X2_IP(cf1_a2a1Zg, cf1_b2b1b0Z, coef, sz_i16*4*2);
    AE_L16X4X2_IP(cf2_a2a1Zg, cf2_b2b1b0Z, coef, sz_i16*4*2);

    cf0_a2a1gZ = AE_SEL16_7632(cf0_a2a1Zg, AE_SHORTSWAP(cf0_a2a1Zg));
    cf1_a2a1gZ = AE_SEL16_7632(cf1_a2a1Zg, AE_SHORTSWAP(cf1_a2a1Zg));
    cf2_a2a1gZ = AE_SEL16_7632(cf2_a2a1Zg, AE_SHORTSWAP(cf2_a2a1Zg));
    cf0_Zb2b1b0 = AE_SEL16_4321(cf0_b2b1b0Z, cf0_b2b1b0Z);
    cf1_Zb2b1b0 = AE_SEL16_4321(cf1_b2b1b0Z, cf1_b2b1b0Z);
    cf2_Zb2b1b0 = AE_SEL16_4321(cf2_b2b1b0Z, cf2_b2b1b0Z);
    // Set final scale
    WUR_AE_SAR(bqriir->gain+2);

    //
    // Load sections' state.
    //

    // Q31
    AE_L32X2X2_I(d2d1_0, d2d1_1, state, 0);
    d2d1_2 = AE_L32X2_I((const ae_int32x2 *)state, 4*sz_i32);

    //
    // Pass N/2 sample pairs through 2 sections.
    //

    alx = AE_LA64_PP(px);
    for ( n=0; n<(N>>1); n++ )
    {
      // Q31
      AE_LA32X2_IP( x01, alx, px );

      //------------------------------------------------------------------
      // Section 0; samples 0 and 1.

      dxx_0 = d2d1_0;
      // Q45 <- Q31*Q14
      fb00 = AE_MULZAAAAQ32X16( dxx_0, x01, cf0_a2a1gZ );
      // Q31 <- Q45 + 2 - 16 w/ saturation
      AE_PKSR32(dxx_0, fb00, 2);
      fb01 = AE_MULZAAAAQ32X16( dxx_0, x01, cf0_a2a1Zg );
      AE_PKSR32(dxx_0, fb01, 2);

      // Q45 <- Q31*Q14
      AE_MULZAAAA2Q32X16(ff00, ff01, d2d1_0, dxx_0, cf0_b2b1b0Z, cf0_Zb2b1b0);
      //ff00 = AE_MULZAAAAQ32X16( d2d1_0, dxx_0, cf0_b2b1b0Z );
      //ff01 = AE_MULZAAAAQ32X16( d2d1_0, dxx_0, cf0_Zb2b1b0 );
      // Update delay elements
      d2d1_0 = dxx_0;
      AE_PKSR32(z01, ff00, 2);
      AE_PKSR32(z01, ff01, 2);

      //------------------------------------------------------------------
      // Section 1; samples 0 and 1.

      dxx_1 = d2d1_1;
      // Q45 <- Q31*Q14
      fb10 = AE_MULZAAAAQ32X16( dxx_1, z01, cf1_a2a1gZ );
      // Q31 <- Q45 + 2 - 16 w/ saturation
      AE_PKSR32(dxx_1, fb10, 2);
      fb11 = AE_MULZAAAAQ32X16( dxx_1, z01, cf1_a2a1Zg );
      AE_PKSR32(dxx_1, fb11, 2);

      // Q45 <- Q31*Q14
      AE_MULZAAAA2Q32X16(ff10, ff11, d2d1_1, dxx_1, cf1_b2b1b0Z, cf1_Zb2b1b0);
      // Update delay elements
      d2d1_1 = dxx_1;
      AE_PKSR32(y01, ff10, 2);
      AE_PKSR32(y01, ff11, 2);

      //------------------------------------------------------------------
      // Section 2; samples 0 and 1.

      dxx_2 = d2d1_2;
      // Q45 <- Q31*Q14
      fb20 = AE_MULZAAAAQ32X16(dxx_2, y01, cf2_a2a1gZ);
      // Q31 <- Q45 + 2 - 16 w/ saturation
      AE_PKSR32(dxx_2, fb20, 2);
      fb21 = AE_MULZAAAAQ32X16(dxx_2, y01, cf2_a2a1Zg);
      AE_PKSR32(dxx_2, fb21, 2);

      // Q45 <- Q31*Q14
      AE_MULZAAAA2Q32X16(ff20, ff21, d2d1_2, dxx_2, cf2_b2b1b0Z, cf2_Zb2b1b0);
      // Update delay elements
      d2d1_2 = dxx_2;

      x01 = AE_ROUND32X2F48SASYM(AE_SLAS64S(ff20), AE_SLAS64S(ff21));
      // Q31
      AE_SA32X2_IP(x01, alr, pr);
    }
    AE_SA64POS_FP(alr, pr);

    //
    // Save sections' state.
    //

    // Q31
    AE_S32X2X2_IP(d2d1_0, d2d1_1, state, sz_i32*2*2);
    AE_S32X2_IP(d2d1_2, castxcc(ae_int32x2,state), sz_i32*2);
  }
  else if (M == 2)
  {
    ae_int16x4 cf0_a2a1Zg, cf0_a2a1gZ, cf0_b2b1b0Z, cf0_Zb2b1b0;
    ae_int16x4 cf1_a2a1Zg, cf1_a2a1gZ, cf1_b2b1b0Z, cf1_Zb2b1b0;

    ae_int32x2 x01, z01, y01;

    ae_int64 fb00, fb01, ff00, ff01; // fbXY == Feedback, ffXY == Feedforward
    ae_int64 fb10, fb11, ff10, ff11; // X - section number, Y - sample number
    ae_int32x2 d2d1_0, dxx_0, d2d1_1, dxx_1;

    //
    // Load coefficients for the m-th pair of sections.
    //

    // Q14
    AE_L16X4X2_IP(cf0_a2a1Zg, cf0_b2b1b0Z, coef, sz_i16*4*2);
    AE_L16X4X2_IP(cf1_a2a1Zg, cf1_b2b1b0Z, coef, sz_i16*4*2);

    cf0_a2a1gZ = AE_SEL16_7632(cf0_a2a1Zg, AE_SHORTSWAP(cf0_a2a1Zg));
    cf1_a2a1gZ = AE_SEL16_7632(cf1_a2a1Zg, AE_SHORTSWAP(cf1_a2a1Zg));
    cf0_Zb2b1b0 = AE_SEL16_4321(cf0_b2b1b0Z, cf0_b2b1b0Z);
    cf1_Zb2b1b0 = AE_SEL16_4321(cf1_b2b1b0Z, cf1_b2b1b0Z);
    // Set final scale
    WUR_AE_SAR(bqriir->gain+2);

    //
    // Load sections' state.
    //

    // Q31
    AE_L32X2X2_I(d2d1_0, d2d1_1, state, 0*sz_i32*2*2);

    //
    // Pass N/2 sample pairs through 2 sections.
    //

    alx = AE_LA64_PP(px);
    for ( n=0; n<(N>>1); n++ )
    {
      // Q31
      AE_LA32X2_IP( x01, alx, px );

      //------------------------------------------------------------------
      // Section 0; samples 0 and 1.

      dxx_0 = d2d1_0;
      // Q45 <- Q31*Q14
      fb00 = AE_MULZAAAAQ32X16( dxx_0, x01, cf0_a2a1gZ );
      // Q31 <- Q45 + 2 - 16 w/ saturation
      AE_PKSR32(dxx_0, fb00, 2);
      fb01 = AE_MULZAAAAQ32X16( dxx_0, x01, cf0_a2a1Zg );
      AE_PKSR32(dxx_0, fb01, 2);

      // Q45 <- Q31*Q14
      AE_MULZAAAA2Q32X16(ff00, ff01, d2d1_0, dxx_0, cf0_b2b1b0Z, cf0_Zb2b1b0);
      // Update delay elements
      d2d1_0 = dxx_0;
      AE_PKSR32(z01, ff00, 2);
      AE_PKSR32(z01, ff01, 2);

      //------------------------------------------------------------------
      // Section 1; samples 0 and 1.

      dxx_1 = d2d1_1;
      // Q45 <- Q31*Q14
      fb10 = AE_MULZAAAAQ32X16( dxx_1, z01, cf1_a2a1gZ );
      // Q31 <- Q45 + 2 - 16 w/ saturation
      AE_PKSR32(dxx_1, fb10, 2);
      fb11 = AE_MULZAAAAQ32X16( dxx_1, z01, cf1_a2a1Zg );
      AE_PKSR32(dxx_1, fb11, 2);

      // Q45 <- Q31*Q14
      AE_MULZAAAA2Q32X16(ff10, ff11, d2d1_1, dxx_1, cf1_b2b1b0Z, cf1_Zb2b1b0);
      // Update delay elements
      d2d1_1 = dxx_1;

      y01 = AE_ROUND32X2F48SASYM(AE_SLAS64S(ff10), AE_SLAS64S(ff11));
      // Q31
      AE_SA32X2_IP( y01, alr, pr );
    }
    AE_SA64POS_FP(alr, pr);

    //
    // Save sections' state.
    //

    // Q31
    AE_S32X2X2_IP(d2d1_0, d2d1_1, state, sz_i32*2*2);
  }
  // Process the last section if the number of biquads is odd.
  else if (M == 1)
  {
    ae_int16x4 cf_a2a1Zg, cf_a2a1gZ, cf_b2b1Zb0, cf_b2b1b0Z;

    ae_int32x2 x01, z01;
    ae_f64     fb0, fb1, ff0, ff1;
    ae_f32x2   st0, t0;

    //
    // Load coefficients for the last section.
    //

    // Q14
    AE_L16X4X2_IP(cf_a2a1Zg, cf_b2b1b0Z, coef, sz_i16*4*2);
    cf_a2a1gZ  = AE_SEL16_7632(cf_a2a1Zg, AE_SHORTSWAP(cf_a2a1Zg));
    cf_b2b1Zb0 = AE_SEL16_7632(cf_b2b1b0Z, AE_SHORTSWAP(cf_b2b1b0Z));
    // Set final scale
    WUR_AE_SAR(bqriir->gain+1);

    //
    // Load last section's state.
    //

    // Q31
    st0 = AE_L32X2_I((const ae_int32x2 *)state, 0);

    //
    // Pass N/2 sample pairs through the section.
    //

    alx = AE_LA64_PP(px);
    __Pragma( "loop_count min=1" );
    for ( n=0; n<(N>>1); n++ )
    {
      // Q31
      AE_LA32X2_IP( x01, alx, px );

      // Q46 <- Q31*Q14 + 1
      fb0 = AE_MULZAAAAFQ32X16( st0, x01, cf_a2a1gZ );
      // Q31 <- Q46 + 17 - 32 w/ saturation
      t0 = AE_TRUNCA32F64S_L( st0, fb0, 17 );
      // Q46 <- Q31*Q14 + 1
      ff0 = AE_MULZAAAAFQ32X16( st0, t0, cf_b2b1Zb0 );
      // Q46 <- Q31*Q14 + 1
      fb1 = AE_MULZAAAAFQ32X16( t0, x01, cf_a2a1Zg );
      // Q31 <- Q46 + 17 - 32 w/ saturation
      st0 = AE_TRUNCA32F64S_L( t0, fb1, 17 );
      // Q46 <- Q31*Q14 + 1
      ff1 = AE_MULZAAAAFQ32X16( t0, st0, cf_b2b1Zb0 );

      // Q31 <- Q46 + gain + 1 - 16 w/ rounding and saturation
      z01 = AE_ROUND32X2F48SASYM(AE_SLAS64S(ff0), AE_SLAS64S(ff1));

      // Q31
      AE_SA32X2_IP( z01, alr, pr );
    }
    AE_SA64POS_FP(alr, pr);

    //
    // Save last section's state.
    //

    // Q31
    AE_S32X2_I(st0, (ae_int32x2 *)state, 0);
  }
  // Scale output data when M is a multiple of 4, if required.
  else if ( bqriir->gain != 0 )
  {
    ae_int32x2 t0, t1;
    ae_valignx2 al2x, al2r;

    WUR_AE_SAR( bqriir->gain );
    px = (const ae_int32x2*)r;
    pr = (      ae_int32x2*)r;

    al2r = AE_ZALIGN128();
    al2x = AE_LA128_PP(px);
    __Pragma("ymemory(pr)");
    for ( n=0; n<(N>>2); n++ )
    {
      AE_LA32X2X2_IP(t0, t1, al2x, castxcc(ae_int32x4,px));
      t0 = AE_SLAS32S( t0 );
      t1 = AE_SLAS32S( t1 );
      AE_SA32X2X2_IP(t0, t1, al2r, castxcc(ae_int32x4,pr));
    }
    AE_SA128POS_FP(al2r, pr);

    if ( N & 2 )
    {
      alx = AE_LA64_PP(px);
      alr = AE_ZALIGN64();
      AE_LA32X2_IP(t0, alx, px);
      t0 = AE_SLAS32S(t0);
      AE_SA32X2_IP(t0, alr, pr);
      AE_SA64POS_FP(alr, pr);
    }
  }
#endif
} // bqriir32x16_df2_process()
