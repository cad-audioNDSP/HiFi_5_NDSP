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
    Biquad Real Block IIR, 32x32-bit, Direct Form II
    C code optimized for HiFi4
    IntegrIT, 2006-2018
*/



/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_iir.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Filter instance structure. */
#include "bqriir32x32_df2_nd_common.h"

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
void bqriir32x32_df2_nd( bqriir32x32_df2_nd_handle_t _bqriir,
                      void    * restrict       s,
                      int32_t * restrict       r,
                const int32_t *                x, int N )
{
  bqriir32x32_df2_nd_ptr_t bqriir = (bqriir32x32_df2_nd_ptr_t)_bqriir;

  const ae_int32x2*          coef_gsos;
        ae_int32x2* restrict state;
  const ae_int32x2* restrict px;
        ae_int32x2* restrict pr;
  ae_valign alx, alr;

  int M;
  int m, n;

  NASSERT( bqriir && bqriir->magic == MAGIC && r && x );
  if (N<=0) return;
  NASSERT( N%2==0 );

  M = bqriir->M;
  alr = AE_ZALIGN64();
  coef_gsos = (const ae_int32x2*)bqriir->coef_gsos;
  state     = (      ae_int32x2*)bqriir->state;

  px = (const ae_int32x2*)x;
  pr = (      ae_int32x2*)r;

  //
  // Perform data block processing for each 4 of successive sections. Use the
  // output array r[N] for temporal storage of inter-section signal.
  //

  for ( m=0; m<((M-1)>>2); m++ )
  {
    ae_int32x2 cf0_gb2, cf0_0b0, cf0_b2b1, cf0_a2a1;
    ae_int32x2 cf1_gb2, cf1_0b0, cf1_b2b1, cf1_a2a1;
    ae_int32x2 cf2_gb2, cf2_0b0, cf2_b2b1, cf2_a2a1;
    ae_int32x2 cf3_gb2, cf3_0b0, cf3_b2b1, cf3_a2a1;

    ae_int32x2 st0, st1, st2, st3;
    ae_int32x2 x01, y01, z01, v01;

    ae_int64 fb00, fb01, ff00, ff01; // fbXY == Feedback, ffXY == Feedforward
    ae_int64 fb10, fb11, ff10, ff11; // X - section number, Y - sample number
    ae_int64 fb20, fb21, ff20, ff21;
    ae_int64 fb30, fb31, ff30, ff31;

    //
    // Load coefficients for the m-th pair of sections.
    //

    AE_L32X2_IP( cf0_gb2 , coef_gsos, +2*4 );
    AE_L32X2_IP( cf0_0b0 , coef_gsos, +2*4 );
    AE_L32X2_IP( cf0_b2b1, coef_gsos, +2*4 );
    AE_L32X2_IP( cf0_a2a1, coef_gsos, +2*4 );
                                       
    AE_L32X2_IP( cf1_gb2 , coef_gsos, +2*4 );
    AE_L32X2_IP( cf1_0b0 , coef_gsos, +2*4 );
    AE_L32X2_IP( cf1_b2b1, coef_gsos, +2*4 );
    AE_L32X2_IP( cf1_a2a1, coef_gsos, +2*4 );

    AE_L32X2_IP( cf2_gb2 , coef_gsos, +2*4 );
    AE_L32X2_IP( cf2_0b0 , coef_gsos, +2*4 );
    AE_L32X2_IP( cf2_b2b1, coef_gsos, +2*4 );
    AE_L32X2_IP( cf2_a2a1, coef_gsos, +2*4 );

    AE_L32X2_IP( cf3_gb2 , coef_gsos, +2*4 );
    AE_L32X2_IP( cf3_0b0 , coef_gsos, +2*4 );
    AE_L32X2_IP( cf3_b2b1, coef_gsos, +2*4 );
    AE_L32X2_IP( cf3_a2a1, coef_gsos, +2*4 );



    //
    // Load sections' state.
    //

    st0 = AE_L32X2_I( state, 0*2*4 );
    st1 = AE_L32X2_I( state, 1*2*4 );
    st2 = AE_L32X2_I( state, 2*2*4 );
    st3 = AE_L32X2_I( state, 3*2*4 );
    ae_int32x2 cf0_a2a1_neg = AE_NEG32(cf0_a2a1);
    ae_int32x2 cf1_a2a1_neg = AE_NEG32(cf1_a2a1);
    ae_int32x2 cf2_a2a1_neg = AE_NEG32(cf2_a2a1);
    WUR_AE_SAR(2);

    //
    // Pass N/2 sample pairs through 4 sections.
    //

    alx = AE_LA64_PP(px);
    __Pragma( "loop_count min=1" );
    for ( n=0; n<(N>>1); n++ )
    {
      AE_LA32X2_IP( x01, alx, px );

      //
      // Section 0, samples 0 and 1.
      //

      //Q62 <- Q30*Q31 + 1
      AE_MULZAAF2D32S_HH_LL(fb00, ff00, cf0_a2a1_neg, cf0_b2b1, st0, st0);
      AE_MULAF32S_HH(fb00, cf0_gb2, x01 );
      
      // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
      AE_PKSRF32( st0, fb00, 1 );

      //Q62
      AE_MULAF32S_LL( ff00, cf0_0b0, st0 );
      //Q62
      AE_MULZAAF2D32S_HH_LL(fb01, ff01, cf0_a2a1_neg, cf0_b2b1, st0, st0);
      AE_MULAF32S_HL(fb01, cf0_gb2, x01 );
      // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
      AE_PKSRF32( st0, fb01, 1 );
      //Q62
      AE_MULAF32S_LL( ff01, cf0_0b0, st0 );

      // Q31 <- Q62 + 1 - 32 w/ saturation
      y01 = AE_TRUNCA32X2F64S( ff00, ff01, 1 );

      //
      // Section 1, samples 0 and 1.
      //

      //Q62 <- Q30*Q31 + 1
      AE_MULZAAF2D32S_HH_LL(fb10, ff10, cf1_a2a1_neg, cf1_b2b1, st1, st1);
      AE_MULAF32S_HH(fb10, cf1_gb2, y01 );

      // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
      AE_PKSRF32( st1, fb10, 1 );
      //Q62
      AE_MULAF32S_LL( ff10, cf1_0b0, st1 );
      //Q62
      AE_MULZAAF2D32S_HH_LL(fb11, ff11, cf1_a2a1_neg, cf1_b2b1, st1, st1);
      AE_MULAF32S_HL(fb11, cf1_gb2, y01 );
      // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
      AE_PKSRF32( st1, fb11, 1 );
      //Q62
      AE_MULAF32S_LL( ff11, cf1_0b0, st1 );

      // Q31 <- Q62 + 1 - 32 w/ saturation
      z01 = AE_TRUNCA32X2F64S( ff10, ff11, 1 );

      //
      // Section 2, samples 0 and 1.
      //

      // Q62 <- Q30*Q31 + 1    
      AE_MULZAAF2D32S_HH_LL(fb20, ff20, cf2_a2a1_neg, cf2_b2b1, st2, st2);
      AE_MULAF32S_HH(fb20, cf2_gb2, z01 );

      // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
      AE_PKSRF32( st2, fb20, 1 );

      //Q62
      AE_MULAF32S_LL( ff20, cf2_0b0, st2 );
      //Q62
      AE_MULZAAF2D32S_HH_LL(fb21, ff21, cf2_a2a1_neg, cf2_b2b1, st2, st2);
      AE_MULAF32S_HL(fb21, cf2_gb2, z01 );
      // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
      AE_PKSRF32( st2, fb21, 1 );
      //Q62
      AE_MULAF32S_LL( ff21, cf2_0b0, st2 );

      // Q31 <- Q62 + 1 - 32 w/ saturation
      v01 = AE_TRUNCA32X2F64S( ff20, ff21, 1 );

      //
      // Section 3, samples 0 and 1.
      //
      //Q61 <- Q30*Q31
      fb30 = AE_MUL32_HH( cf3_gb2, v01 );
      AE_MULSSD32_HH_LL( fb30, cf3_a2a1, st3 );
      //Q61
      ff30 = AE_MULZAAD32_HH_LL( cf3_b2b1, st3 );

      // Q31 <- Q61 + 2 - 32 w/ rounding, saturation
      AE_PKSRF32( st3, fb30, 2 );
      //Q61
      AE_MULA32_LL( ff30, cf3_0b0, st3 );
      //Q61
      fb31 = AE_MUL32_HL( cf3_gb2, v01 );
      AE_MULSSD32_HH_LL( fb31, cf3_a2a1, st3 );
      //Q61
      ff31 = AE_MULZAAD32_HH_LL( cf3_b2b1, st3 );
      // Q31 <- Q61 + 2 - 32 w/ rounding, saturation
      AE_PKSRF32( st3, fb31, 2 );
      //Q61
      AE_MULA32_LL( ff31, cf3_0b0, st3 );

      // Q63 <- Q61 + 2 /w saturation
      ff30 = AE_SLAS64S(ff30);
      ff31 = AE_SLAS64S(ff31);

      // Q31 <- Q63 - 32 w/ rounding
      x01 = AE_ROUND32X2F64SASYM(ff30, ff31);
      AE_SA32X2_IP( x01, alr, pr );
    }
    AE_SA64POS_FP(alr, pr);

    //
    // Save sections' state.
    //

    AE_S32X2_IP( st0, state, +2*4 );
    AE_S32X2_IP( st1, state, +2*4 );
    AE_S32X2_IP( st2, state, +2*4 );
    AE_S32X2_IP( st3, state, +2*4 );

    // Second to last pair of sections are fed with output signal of the
    // previous pair.
    px = pr = (ae_int32x2*)( (uintptr_t)pr - N*4 );
  }
  //
  // Process last 1-4 sections
  //
  M = M - ((M-1)&(~3));
  if (M == 4)
  {
      ae_int32x2 cf0_gb2, cf0_0b0, cf0_b2b1, cf0_a2a1;
      ae_int32x2 cf1_gb2, cf1_0b0, cf1_b2b1, cf1_a2a1;
      ae_int32x2 cf2_gb2, cf2_0b0, cf2_b2b1, cf2_a2a1;
      ae_int32x2 cf3_gb2, cf3_0b0, cf3_b2b1, cf3_a2a1;

      ae_int32x2 st0, st1, st2, st3;
      ae_int32x2 x01, y01, z01, v01;

      ae_int64 fb00, fb01, ff00, ff01; // fbXY == Feedback, ffXY == Feedforward
      ae_int64 fb10, fb11, ff10, ff11; // X - section number, Y - sample number
      ae_int64 fb20, fb21, ff20, ff21;
      ae_int64 fb30, fb31, ff30, ff31;

      //
      // Load coefficients for the m-th pair of sections.
      //

      AE_L32X2_IP( cf0_gb2 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf0_0b0 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf0_b2b1, coef_gsos, +2*4 );
      AE_L32X2_IP( cf0_a2a1, coef_gsos, +2*4 );

      AE_L32X2_IP( cf1_gb2 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf1_0b0 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf1_b2b1, coef_gsos, +2*4 );
      AE_L32X2_IP( cf1_a2a1, coef_gsos, +2*4 );

      AE_L32X2_IP( cf2_gb2 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf2_0b0 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf2_b2b1, coef_gsos, +2*4 );
      AE_L32X2_IP( cf2_a2a1, coef_gsos, +2*4 );

      AE_L32X2_IP( cf3_gb2 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf3_0b0 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf3_b2b1, coef_gsos, +2*4 );
      AE_L32X2_IP( cf3_a2a1, coef_gsos, +2*4 );

      ae_int32x2 cf0_a2a1_neg = AE_NEG32(cf0_a2a1);
      ae_int32x2 cf1_a2a1_neg = AE_NEG32(cf1_a2a1);
      ae_int32x2 cf2_a2a1_neg = AE_NEG32(cf2_a2a1);

      //
      // Load sections' state.
      //

      st0 = AE_L32X2_I( state, 0*2*4 );
      st1 = AE_L32X2_I( state, 1*2*4 );
      st2 = AE_L32X2_I( state, 2*2*4 );
      st3 = AE_L32X2_I( state, 3*2*4 );

      //
      // Pass N/2 sample pairs through 4 sections.
      //
      WUR_AE_SAR( bqriir->gain+2 );

      alx = AE_LA64_PP(px);
      __Pragma( "loop_count min=1" );
      for ( n=0; n<(N>>1); n++ )
      {
          AE_LA32X2_IP( x01, alx, px );

          //
          // Section 0, samples 0 and 1.
          //

          //Q62 <- Q30*Q31 + 1
          AE_MULZAAF2D32S_HH_LL(fb00, ff00, cf0_a2a1_neg, cf0_b2b1, st0, st0);
          AE_MULAF32S_HH(fb00, cf0_gb2, x01 );

          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st0, fb00, 1 );

          //Q62
          AE_MULAF32S_LL( ff00, cf0_0b0, st0 );
          //Q62
          AE_MULZAAF2D32S_HH_LL(fb01, ff01, cf0_a2a1_neg, cf0_b2b1, st0, st0);
          AE_MULAF32S_HL(fb01, cf0_gb2, x01 );
          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st0, fb01, 1 );
          //Q62
          AE_MULAF32S_LL( ff01, cf0_0b0, st0 );

          // Q31 <- Q62 + 1 - 32 w/ saturation
          y01 = AE_TRUNCA32X2F64S( ff00, ff01, 1 );

          //
          // Section 1, samples 0 and 1.
          //

          //Q62 <- Q30*Q31 + 1
          AE_MULZAAF2D32S_HH_LL(fb10, ff10, cf1_a2a1_neg, cf1_b2b1, st1, st1);
          AE_MULAF32S_HH(fb10, cf1_gb2, y01 );

          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st1, fb10, 1 );
          //Q62
          AE_MULAF32S_LL( ff10, cf1_0b0, st1 );
          //Q62
          AE_MULZAAF2D32S_HH_LL(fb11, ff11, cf1_a2a1_neg, cf1_b2b1, st1, st1);
          AE_MULAF32S_HL(fb11, cf1_gb2, y01 );
          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st1, fb11, 1 );
          //Q62
          AE_MULAF32S_LL( ff11, cf1_0b0, st1 );

          // Q31 <- Q62 + 1 - 32 w/ saturation
          z01 = AE_TRUNCA32X2F64S( ff10, ff11, 1 );


          //
          // Section 2, samples 0 and 1.
          //

          // Q62 <- Q30*Q31 + 1    
          AE_MULZAAF2D32S_HH_LL(fb20, ff20, cf2_a2a1_neg, cf2_b2b1, st2, st2);
          AE_MULAF32S_HH(fb20, cf2_gb2, z01 );

          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st2, fb20, 1 );

          //Q62
          AE_MULAF32S_LL( ff20, cf2_0b0, st2 );
          //Q62
          AE_MULZAAF2D32S_HH_LL(fb21, ff21, cf2_a2a1_neg, cf2_b2b1, st2, st2);
          AE_MULAF32S_HL(fb21, cf2_gb2, z01 );
          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st2, fb21, 1 );
          //Q62
          AE_MULAF32S_LL( ff21, cf2_0b0, st2 );

          // Q31 <- Q62 + 1 - 32 w/ saturation
          v01 = AE_TRUNCA32X2F64S( ff20, ff21, 1 );

          //
          // Section 3, samples 0 and 1.
          //

          fb30 = AE_MUL32_HH( cf3_gb2, v01 );
          AE_MULSSD32_HH_LL( fb30, cf3_a2a1, st3 );

          ff30 = AE_MULZAAD32_HH_LL( cf3_b2b1, st3 );

          // Q31 <- Q61 + 2 - 32 w/ rounding, saturation
          AE_PKSRF32( st3, fb30, 2 );

          AE_MULA32_LL( ff30, cf3_0b0, st3 );

          fb31 = AE_MUL32_HL( cf3_gb2, v01 );
          AE_MULSSD32_HH_LL( fb31, cf3_a2a1, st3 );

          ff31 = AE_MULZAAD32_HH_LL( cf3_b2b1, st3 );

          AE_PKSRF32( st3, fb31, 2 );

          AE_MULA32_LL( ff31, cf3_0b0, st3 );

          // Q63 <- Q61 + 2 /w saturation
          ff30 = AE_SLAS64S(ff30);
          ff31 = AE_SLAS64S(ff31);

          // Q31 <- Q63 - 32 w/ rounding
          x01 = AE_ROUND32X2F64SASYM(ff30, ff31);
          AE_SA32X2_IP( x01, alr, pr );
      }
      AE_SA64POS_FP(alr, pr);

      //
      // Save sections' state.
      //

      AE_S32X2_IP( st0, state, +2*4 );
      AE_S32X2_IP( st1, state, +2*4 );
      AE_S32X2_IP( st2, state, +2*4 );
      AE_S32X2_IP( st3, state, +2*4 );

  }
  else if (M == 3)
  {
      ae_int32x2 cf0_gb2, cf0_0b0, cf0_b2b1, cf0_a2a1;
      ae_int32x2 cf1_gb2, cf1_0b0, cf1_b2b1, cf1_a2a1;
      ae_int32x2 cf2_gb2, cf2_0b0, cf2_b2b1, cf2_a2a1;

      ae_int32x2 st0, st1, st2;
      ae_int32x2 x01, y01, z01;

      ae_int64 fb00, fb01, ff00, ff01; // fbXY == Feedback, ffXY == Feedforward
      ae_int64 fb10, fb11, ff10, ff11; // X - section number, Y - sample number
      ae_int64 fb20, fb21, ff20, ff21;

      //
      // Load coefficients for the m-th pair of sections.
      //

      AE_L32X2_IP( cf0_gb2 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf0_0b0 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf0_b2b1, coef_gsos, +2*4 );
      AE_L32X2_IP( cf0_a2a1, coef_gsos, +2*4 );

      AE_L32X2_IP( cf1_gb2 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf1_0b0 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf1_b2b1, coef_gsos, +2*4 );
      AE_L32X2_IP( cf1_a2a1, coef_gsos, +2*4 );

      AE_L32X2_IP( cf2_gb2 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf2_0b0 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf2_b2b1, coef_gsos, +2*4 );
      AE_L32X2_IP( cf2_a2a1, coef_gsos, +2*4 );


      //
      // Load sections' state.
      //

      st0 = AE_L32X2_I( state, 0*2*4 );
      st1 = AE_L32X2_I( state, 1*2*4 );
      st2 = AE_L32X2_I( state, 2*2*4 );

      ae_int32x2 cf0_a2a1_neg = AE_NEG32(cf0_a2a1);
      ae_int32x2 cf1_a2a1_neg = AE_NEG32(cf1_a2a1);
      ae_int32x2 cf2_a2a1_neg = AE_NEG32(cf2_a2a1);

      WUR_AE_SAR( bqriir->gain+1 );

      //
      // Pass N/2 sample pairs through 3 sections.
      //

      alx = AE_LA64_PP(px);
      __Pragma( "loop_count min=1" );
      for ( n=0; n<(N>>1); n++ )
      {
          AE_LA32X2_IP( x01, alx, px );

          //
          // Section 0, samples 0 and 1.
          //

          //Q62 <- Q30*Q31 + 1
          AE_MULZAAF2D32S_HH_LL(fb00, ff00, cf0_a2a1_neg, cf0_b2b1, st0, st0);
          AE_MULAF32S_HH(fb00, cf0_gb2, x01 );

          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st0, fb00, 1 );

          //Q62
          AE_MULAF32S_LL( ff00, cf0_0b0, st0 );
          //Q62
          AE_MULZAAF2D32S_HH_LL(fb01, ff01, cf0_a2a1_neg, cf0_b2b1, st0, st0);
          AE_MULAF32S_HL(fb01, cf0_gb2, x01 );
          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st0, fb01, 1 );
          //Q62
          AE_MULAF32S_LL( ff01, cf0_0b0, st0 );

          // Q31 <- Q62 + 1 - 32 w/ saturation
          y01 = AE_TRUNCA32X2F64S( ff00, ff01, 1 );

          //
          // Section 1, samples 0 and 1.
          //

          //Q62 <- Q30*Q31 + 1
          AE_MULZAAF2D32S_HH_LL(fb10, ff10, cf1_a2a1_neg, cf1_b2b1, st1, st1);
          AE_MULAF32S_HH(fb10, cf1_gb2, y01 );

          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st1, fb10, 1 );
          //Q62
          AE_MULAF32S_LL( ff10, cf1_0b0, st1 );
          //Q62
          AE_MULZAAF2D32S_HH_LL(fb11, ff11, cf1_a2a1_neg, cf1_b2b1, st1, st1);
          AE_MULAF32S_HL(fb11, cf1_gb2, y01 );
          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st1, fb11, 1 );
          //Q62
          AE_MULAF32S_LL( ff11, cf1_0b0, st1 );

          // Q31 <- Q62 + 1 - 32 w/ saturation
          z01 = AE_TRUNCA32X2F64S( ff10, ff11, 1 );

          //
          // Section 2, samples 0 and 1.
          //

          // Q62 <- Q30*Q31 + 1    
          AE_MULZAAF2D32S_HH_LL(fb20, ff20, cf2_a2a1_neg, cf2_b2b1, st2, st2);
          AE_MULAF32S_HH(fb20, cf2_gb2, z01 );

          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st2, fb20, 1 );

          //Q62
          AE_MULAF32S_LL( ff20, cf2_0b0, st2 );
          //Q62
          AE_MULZAAF2D32S_HH_LL(fb21, ff21, cf2_a2a1_neg, cf2_b2b1, st2, st2);
          AE_MULAF32S_HL(fb21, cf2_gb2, z01 );
          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st2, fb21, 1 );
          //Q62
          AE_MULAF32S_LL( ff21, cf2_0b0, st2 );

          // Q63 <- Q62 + 1 /w saturation
          ff20 = AE_SLAS64S(ff20);
          ff21 = AE_SLAS64S(ff21);

          // Q31 <- Q63 - 32 w/ rounding
          x01 = AE_ROUND32X2F64SASYM(ff20, ff21);
          AE_SA32X2_IP( x01, alr, pr );
      }
      AE_SA64POS_FP(alr, pr);

      //
      // Save sections' state.
      //

      AE_S32X2_IP( st0, state, +2*4 );
      AE_S32X2_IP( st1, state, +2*4 );
      AE_S32X2_IP( st2, state, +2*4 );

  }
  else if (M == 2)
  {
      ae_int32x2 cf0_gb2, cf0_0b0, cf0_b2b1, cf0_a2a1;
      ae_int32x2 cf1_gb2, cf1_0b0, cf1_b2b1, cf1_a2a1;

      ae_int32x2 st0, st1;
      ae_int32x2 x01, y01;

      ae_int64 fb00, fb01, ff00, ff01; // fbXY == Feedback, ffXY == Feedforward
      ae_int64 fb10, fb11, ff10, ff11; // X - section number, Y - sample number

      //
      // Load coefficients for the m-th pair of sections.
      //

      AE_L32X2_IP( cf0_gb2 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf0_0b0 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf0_b2b1, coef_gsos, +2*4 );
      AE_L32X2_IP( cf0_a2a1, coef_gsos, +2*4 );

      AE_L32X2_IP( cf1_gb2 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf1_0b0 , coef_gsos, +2*4 );
      AE_L32X2_IP( cf1_b2b1, coef_gsos, +2*4 );
      AE_L32X2_IP( cf1_a2a1, coef_gsos, +2*4 );

      //
      // Load sections' state.
      //

      st0 = AE_L32X2_I( state, 0*2*4 );
      st1 = AE_L32X2_I( state, 1*2*4 );

      ae_int32x2 cf0_a2a1_neg = AE_NEG32(cf0_a2a1);
      ae_int32x2 cf1_a2a1_neg = AE_NEG32(cf1_a2a1);

      WUR_AE_SAR( bqriir->gain+1 );

      //
      // Pass N/2 sample pairs through 2 sections.
      //

      alx = AE_LA64_PP(px);
      __Pragma( "loop_count min=1" );
      for ( n=0; n<(N>>1); n++ )
      {
          AE_LA32X2_IP( x01, alx, px );

          //
          // Section 0, samples 0 and 1.
          //

          //Q62 <- Q30*Q31 + 1
          AE_MULZAAF2D32S_HH_LL(fb00, ff00, cf0_a2a1_neg, cf0_b2b1, st0, st0);
          AE_MULAF32S_HH(fb00, cf0_gb2, x01 );

          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st0, fb00, 1 );

          //Q62
          AE_MULAF32S_LL( ff00, cf0_0b0, st0 );
          //Q62
          AE_MULZAAF2D32S_HH_LL(fb01, ff01, cf0_a2a1_neg, cf0_b2b1, st0, st0);
          AE_MULAF32S_HL(fb01, cf0_gb2, x01 );
          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st0, fb01, 1 );
          //Q62
          AE_MULAF32S_LL( ff01, cf0_0b0, st0 );

          // Q31 <- Q62 + 1 - 32 w/ saturation
          y01 = AE_TRUNCA32X2F64S( ff00, ff01, 1 );

          //
          // Section 2, samples 0 and 1.
          //

          //Q62 <- Q30*Q31 + 1
          AE_MULZAAF2D32S_HH_LL(fb10, ff10, cf1_a2a1_neg, cf1_b2b1, st1, st1);
          AE_MULAF32S_HH(fb10, cf1_gb2, y01 );

          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st1, fb10, 1 );
          //Q62
          AE_MULAF32S_LL( ff10, cf1_0b0, st1 );
          //Q62
          AE_MULZAAF2D32S_HH_LL(fb11, ff11, cf1_a2a1_neg, cf1_b2b1, st1, st1);
          AE_MULAF32S_HL(fb11, cf1_gb2, y01 );
          // Q31 <- Q62 + 1 - 32 w/ rounding, saturation
          AE_PKSRF32( st1, fb11, 1 );
          //Q62
          AE_MULAF32S_LL( ff11, cf1_0b0, st1 );

          // Q63 <- Q62 + 1 /w saturation
          ff10 = AE_SLAS64S(ff10);
          ff11 = AE_SLAS64S(ff11);

          // Q31 <- Q63 - 32 w/ rounding
          x01 = AE_ROUND32X2F64SASYM(ff10, ff11);
          AE_SA32X2_IP( x01, alr, pr );

      }
      AE_SA64POS_FP(alr, pr);

      //
      // Save sections' state.
      //

      AE_S32X2_IP( st0, state, +2*4 );
      AE_S32X2_IP( st1, state, +2*4 );

  }
  else if ( M == 1 )
  {
    ae_int32x2 cf_gb2, cf_0b0, cf_b2b1, cf_a2a1;

    ae_int32x2 st0, x01;

    WUR_AE_SAR( bqriir->gain+2 );

    ae_int64 fb0, ff0; // fbY == Feedback, ffY == Feedforward
    ae_int64 fb1, ff1; // Y - sample number

    //
    // Load coefficients for the last section.
    //

    cf_gb2  = AE_L32X2_I( coef_gsos, 0*2*4 );
    cf_0b0  = AE_L32X2_I( coef_gsos, 1*2*4 );
    cf_b2b1 = AE_L32X2_I( coef_gsos, 2*2*4 );
    cf_a2a1 = AE_L32X2_I( coef_gsos, 3*2*4 );

    //
    // Load last section's state.
    //

    st0 = AE_L32X2_I( state, 0*2*4 );

    //
    // Pass N/2 sample pairs through the section.
    //

    alx = AE_LA64_PP(px);
    for ( n=0; n<(N>>1); n++ )
    {
      AE_LA32X2_IP( x01, alx, px );

      // Q61 <- Q30*Q31
      fb0 = AE_MUL32_HH( cf_gb2, x01 );
      AE_MULSSD32_HH_LL( fb0, cf_a2a1, st0 );

      ff0 = AE_MULZAAD32_HH_LL( cf_b2b1, st0 );

      // Q31 <- Q61 + 2 - 32 w/ rounding, saturation
      AE_PKSRF32( st0, fb0, 2 );

      AE_MULA32_LL( ff0, cf_0b0, st0 );

      fb1 = AE_MUL32_HL( cf_gb2, x01 );
      AE_MULSSD32_HH_LL( fb1, cf_a2a1, st0 );

      ff1 = AE_MULZAAD32_HH_LL( cf_b2b1, st0 );

      AE_PKSRF32( st0, fb1, 2 );

      AE_MULA32_LL( ff1, cf_0b0, st0 );

      // Q63 <- Q61 + 2 /w saturation
      ff0 = AE_SLAS64S(ff0);
      ff1 = AE_SLAS64S(ff1);

      // Q31 <- Q63 - 32 w/ rounding
      x01 = AE_ROUND32X2F64SASYM(ff0, ff1);
      AE_SA32X2_IP( x01, alr, pr );
    }
    AE_SA64POS_FP(alr, pr);

    //
    // Save last section's state.
    //

    AE_S32X2_I( st0, state, 0 );
  } 
} // bqriir32x32_df2_nd_process()
