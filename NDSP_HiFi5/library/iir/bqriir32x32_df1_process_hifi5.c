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
    Biquad Real Block IIR, 32x32-bit, Direct Form I
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
#include "bqriir32x32_df1_common.h"
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
void bqriir32x32_df1( bqriir32x32_df1_handle_t _bqriir,
                      void    * restrict       s,
                      int32_t * restrict       r,
                const int32_t *                x , int N)
{
  bqriir32x32_df1_ptr_t bqriir = (bqriir32x32_df1_ptr_t)_bqriir;

  const ae_int32x4* restrict coef_gsos;
        ae_int32x4* restrict state;
  const ae_int32x2* restrict px;
        ae_int32x2* restrict pr;
  ae_valign alx, alr;
  int sh=bqriir->gain+1;

  int M;
  int m, n;

  NASSERT( bqriir && bqriir->magic == MAGIC && r && x );
  if (N<=0) return;
  NASSERT( N%2==0 );

  M = bqriir->M;
  alr = AE_ZALIGN64();
  coef_gsos = (const ae_int32x4*)bqriir->coef_gsos;
  state     = (      ae_int32x4*)bqriir->state;

  px = (const ae_int32x2*)x;
  pr = (      ae_int32x2*)r;

  //
  // Perform data block processing for each of 4 successive sections. Use the
  // output array r[N] for temporal storage of inter-section signal.
  //

  for ( m=0; m<(M/4); m++ )
  {
    ae_int32x2 cf0_b21, cf0_b10, cf0_a21;
    ae_int32x2 cf1_b21, cf1_b10, cf1_a21;
    ae_int32x2 cf2_b21, cf2_b10, cf2_a21;
    ae_int32x2 cf3_b21, cf3_b10, cf3_a21;
    ae_int32x2 st_0ff, st_0fb, st_1ff, st_1fb, st_2ff, st_2fb, st_3ff, st_3fb;
    ae_int32x2 x01;
    // accXY: X - section id, Y - sample id
    ae_int64 acc00, acc01, acc10, acc11, acc20, acc21, acc30, acc31;

    //
    // Load coefficients of 4 sections.
    //

    AE_L32X2X2_IP(cf0_b21, cf0_a21, coef_gsos, 2*sz_i32*2);
    AE_L32X2X2_IP(cf0_b10, cf1_b21, coef_gsos, 2*sz_i32*2);
    AE_L32X2X2_IP(cf1_a21, cf1_b10, coef_gsos, 2*sz_i32*2);
    AE_L32X2X2_IP(cf2_b21, cf2_a21, coef_gsos, 2*sz_i32*2);
    AE_L32X2X2_IP(cf2_b10, cf3_b21, coef_gsos, 2*sz_i32*2);
    AE_L32X2X2_IP(cf3_a21, cf3_b10, coef_gsos, 2*sz_i32*2);

    //
    // Load sections' state.
    //

    AE_L32X2X2_I(st_0ff, st_0fb, state, 0*2*sz_i32*2 );
    AE_L32X2X2_I(st_1ff, st_1fb, state, 1*2*sz_i32*2 );
    AE_L32X2X2_I(st_2ff, st_2fb, state, 2*2*sz_i32*2 );
    AE_L32X2X2_I(st_3ff, st_3fb, state, 3*2*sz_i32*2 );

    //
    // Pass N/2 sample pairs through 4 sections.
    //

    alx = AE_LA64_PP(px);
    for ( n=0; n<(N>>1); n++ )
    {
        // Load 2 input samples
        AE_LA32X2_IP(x01, alx, px);

        acc00 = AE_MULF32S_LH(cf0_b10, x01);
        acc10 = AE_MULF32S_LL(cf1_b10, st_0fb);
        acc20 = AE_MULF32S_LL(cf2_b10, st_1fb);
        acc30 = AE_MULF32S_LL(cf3_b10, st_2fb);
        acc01 = AE_MULF32S_HL(cf0_b21, st_0ff);
        acc11 = AE_MULF32S_HL(cf1_b21, st_1ff);
        acc21 = AE_MULF32S_HL(cf2_b21, st_2ff);
        acc31 = AE_MULF32S_HL(cf3_b21, st_3ff);
        AE_MULAAF2D32S_HH_LL(acc00, acc10, cf0_b21, cf1_b21, st_0ff, st_1ff);
        AE_MULSSF2D32S_HH_LL(acc00, acc10, cf0_a21, cf1_a21, st_0fb, st_1fb);
        AE_MULAAF2D32S_HH_LL(acc20, acc30, cf2_b21, cf3_b21, st_2ff, st_3ff);
        AE_MULSSF2D32S_HH_LL(acc20, acc30, cf2_a21, cf3_a21, st_2fb, st_3fb);
        AE_PKSRF32(st_0fb, acc00, 1);
        AE_PKSRF32(st_1fb, acc10, 1);
        AE_PKSRF32(st_2fb, acc20, 1);
        AE_PKSRF32(st_3fb, acc30, 1);

        AE_MOVD32X4(st_2ff,st_3ff,st_1fb,st_2fb);;
        AE_MOVD32X4(st_0ff,st_1ff,x01,st_0fb);;
        AE_MULAAF2D32S_HH_LL(acc01, acc11, cf0_b10, cf1_b10, x01, st_0fb);
        AE_MULSSF2D32S_HH_LL(acc01, acc11, cf0_a21, cf1_a21, st_0fb, st_1fb);
        AE_MULAAF2D32S_HH_LL(acc21, acc31, cf2_b10, cf3_b10, st_1fb, st_2fb);
        AE_MULSSF2D32S_HH_LL(acc21, acc31, cf2_a21, cf3_a21, st_2fb, st_3fb);
        AE_PKSRF32(st_0fb, acc01, 1);
        AE_PKSRF32(st_1fb, acc11, 1);
        AE_PKSRF32(st_2fb, acc21, 1);
        AE_PKSRF32(st_3fb, acc31, 1);

        AE_SA32X2_IP(st_3fb, alr, pr);
    }
    AE_SA64POS_FP(alr, pr);

    //
    // Save sections' state.
    //

    AE_S32X2X2_IP(st_0ff, st_0fb, state, 2*sz_i32*2);
    AE_S32X2X2_IP(st_1ff, st_1fb, state, 2*sz_i32*2);
    AE_S32X2X2_IP(st_2ff, st_2fb, state, 2*sz_i32*2);
    AE_S32X2X2_IP(st_3ff, st_3fb, state, 2*sz_i32*2);

    // Second to last pair of sections are fed with output signal of the
    // previous pair.
    px = pr = (ae_int32x2*)r;
  }

  //
  // Process remaining sections (last 0..3)
  //
  M = M & 3;
  WUR_AE_SAR( bqriir->gain+1 );

  if (M == 3)
  {
    ae_int32x2 cf0_b21, cf0_b10, cf0_a21;
    ae_int32x2 cf1_b21, cf1_b10, cf1_a21;
    ae_int32x2 cf2_b21, cf2_b10, cf2_a21;
    ae_int32x2 st_0ff, st_0fb_1ff, st_1fb_2ff, st_2fb, stx, sty;
    ae_int32x2 x01, y01;
    // accXY: X - section id, Y - sample id
    ae_int64 acc00, acc01, acc10, acc11, acc20, acc21;

    //
    // Load coefficients.
    //

    AE_L32X2X2_IP(cf0_b21, cf0_a21, coef_gsos, 2*sz_i32*2);
    AE_L32X2X2_IP(cf0_b10, cf1_b21, coef_gsos, 2*sz_i32*2);
    AE_L32X2X2_IP(cf1_a21, cf1_b10, coef_gsos, 2*sz_i32*2);
    AE_L32X2X2_IP(cf2_b21, cf2_a21, coef_gsos, 2*sz_i32*2);
    cf2_b10 = AE_L32X2_I((const ae_int32x2 *)coef_gsos, 0);

    //
    // Load sections' state.
    //

    AE_L32X2X2_I(st_0ff    , st_0fb_1ff, state, 0*4*sz_i32);
    AE_L32X2X2_I(st_1fb_2ff, st_2fb    , state, 1*4*sz_i32);

    //
    // Pass N/2 sample pairs through 2 sections.
    //

    alx = AE_LA64_PP(px);
    for ( n=0; n<(N>>1); n++ )
    {
      AE_LA32X2_IP(x01, alx, px);

      stx = st_0fb_1ff;

      acc00 = AE_MULZAAFD32S_HH_LL(cf0_b21, st_0ff);
      AE_MULAF32S_LH(acc00, cf0_b10, x01);
      AE_MULSSFD32S_HH_LL(acc00, cf0_a21, st_0fb_1ff);
      AE_PKSRF32(st_0fb_1ff, acc00, 1);

      acc01 = AE_MULF32S_HL(cf0_b21, st_0ff);
      AE_MULAAFD32S_HH_LL(acc01, cf0_b10, x01);
      AE_MULSSFD32S_HH_LL(acc01, cf0_a21, st_0fb_1ff);
      AE_PKSRF32(st_0fb_1ff, acc01, 1);

      st_0ff = x01;
      sty = st_1fb_2ff;

      acc10 = AE_MULZAAFD32S_HH_LL(cf1_b21, stx);
      AE_MULAF32S_LH(acc10, cf1_b10, st_0fb_1ff);
      AE_MULSSFD32S_HH_LL(acc10, cf1_a21, st_1fb_2ff);
      AE_PKSRF32(st_1fb_2ff, acc10, 1);

      acc11 = AE_MULF32S_HL(cf1_b21, stx);
      AE_MULAAFD32S_HH_LL(acc11, cf1_b10, st_0fb_1ff);
      AE_MULSSFD32S_HH_LL(acc11, cf1_a21, st_1fb_2ff);
      AE_PKSRF32(st_1fb_2ff, acc11, 1);

      acc20 = AE_MULZAAFD32S_HH_LL(cf2_b21, sty);
      AE_MULAF32S_LH(acc20, cf2_b10, st_1fb_2ff);
      AE_MULSSFD32S_HH_LL(acc20, cf2_a21, st_2fb);
      AE_PKSRF32(st_2fb, acc20, 1);

      acc21 = AE_MULF32S_HL(cf2_b21, sty);
      AE_MULAAFD32S_HH_LL(acc21, cf2_b10, st_1fb_2ff);
      AE_MULSSFD32S_HH_LL(acc21, cf2_a21, st_2fb);
      AE_PKSRF32(st_2fb, acc21, 1);

      //y01 = AE_ROUND32X2F64SASYM(AE_SLAS64S(acc20), AE_SLAS64S(acc21));
      y01 = AE_TRUNCA32X2F64S(acc20,acc21,sh);
      AE_SA32X2_IP(y01, alr, pr);
    }
    AE_SA64POS_FP(alr, pr);

    //
    // Save sections' state.
    //

    AE_S32X2X2_IP(st_0ff    , st_0fb_1ff, state, 4*sz_i32);
    AE_S32X2X2_IP(st_1fb_2ff, st_2fb    , state, 4*sz_i32);
  }
  else if (M == 2) // process 2 sections
  {
    ae_int32x2 cf0_b21, cf0_b10, cf0_a21;
    ae_int32x2 cf1_b21, cf1_b10, cf1_a21;
    ae_int32x2 st_0ff, st_0fb_1ff, st_1fb, stx;
    ae_int32x2 x01, y01;
    // accXY: X - section id, Y - sample id
    ae_int64 acc00, acc01, acc10, acc11;

    //
    // Load coefficients.
    //

    AE_L32X2X2_IP(cf0_b21, cf0_a21, coef_gsos, 2*sz_i32*2);
    AE_L32X2X2_IP(cf0_b10, cf1_b21, coef_gsos, 2*sz_i32*2);
    AE_L32X2X2_IP(cf1_a21, cf1_b10, coef_gsos, 2*sz_i32*2);

    //
    // Load sections' state.
    //

    st_0ff     = AE_L32X2_I((const ae_int32x2 *)state, 0*2*sz_i32);
    st_0fb_1ff = AE_L32X2_I((const ae_int32x2 *)state, 1*2*sz_i32);
    st_1fb     = AE_L32X2_I((const ae_int32x2 *)state, 2*2*sz_i32);

    //
    // Pass N/2 sample pairs through 2 sections.
    //

    alx = AE_LA64_PP(px);
    for ( n=0; n<(N>>1); n++ )
    {
      AE_LA32X2_IP(x01, alx, px);

      stx = st_0fb_1ff;

      acc00 = AE_MULZAAFD32S_HH_LL(cf0_b21, st_0ff);
      AE_MULAF32S_LH(acc00, cf0_b10, x01);
      AE_MULSSFD32S_HH_LL(acc00, cf0_a21, st_0fb_1ff);
      AE_PKSRF32(st_0fb_1ff, acc00, 1);

      acc01 = AE_MULF32S_HL(cf0_b21, st_0ff);
      AE_MULAAFD32S_HH_LL(acc01, cf0_b10, x01);
      AE_MULSSFD32S_HH_LL(acc01, cf0_a21, st_0fb_1ff);
      AE_PKSRF32(st_0fb_1ff, acc01, 1);

      st_0ff = x01;

      acc10 = AE_MULZAAFD32S_HH_LL(cf1_b21, stx);
      AE_MULAF32S_LH(acc10, cf1_b10, st_0fb_1ff);
      AE_MULSSFD32S_HH_LL(acc10, cf1_a21, st_1fb);
      AE_PKSRF32(st_1fb, acc10, 1);

      acc11 = AE_MULF32S_HL(cf1_b21, stx);
      AE_MULAAFD32S_HH_LL(acc11, cf1_b10, st_0fb_1ff);
      AE_MULSSFD32S_HH_LL(acc11, cf1_a21, st_1fb);
      AE_PKSRF32(st_1fb, acc11, 1);

     // y01 = AE_ROUND32X2F64SASYM(AE_SLAS64S(acc10), AE_SLAS64S(acc11));
      y01 = AE_TRUNCA32X2F64S(acc10,acc11,sh);
      AE_SA32X2_IP(y01, alr, pr);
    }
    AE_SA64POS_FP(alr, pr);

    //
    // Save sections' state.
    //

    AE_S32X2_IP( st_0ff    , castxcc(ae_int32x2,state), 2*sz_i32);
    AE_S32X2_IP( st_0fb_1ff, castxcc(ae_int32x2,state), 2*sz_i32);
    AE_S32X2_IP( st_1fb    , castxcc(ae_int32x2,state), 2*sz_i32);
  }
  else if ( M == 1 )
  {
    ae_int32x2 cf_b21, cf_b10, cf_a21;
    ae_int32x2 st0, st1, st2;
    ae_int32x2 x01, y01;
    // accY: Y - sample id
    ae_int64 acc0, acc1;


    //
    // Load coefficients for the last section.
    //

    cf_b21 = AE_L32X2_I((const ae_int32x2 *)coef_gsos, 0*2*sz_i32);
    cf_a21 = AE_L32X2_I((const ae_int32x2 *)coef_gsos, 1*2*sz_i32);
    cf_b10 = AE_L32X2_I((const ae_int32x2 *)coef_gsos, 2*2*sz_i32);

    //
    // Load last section's state.
    //

    st0 = AE_L32X2_I((const ae_int32x2 *)state, 0*2*sz_i32);
    st1 = AE_L32X2_I((const ae_int32x2 *)state, 1*2*sz_i32);
    st2 = AE_L32X2_I((const ae_int32x2 *)state, 2*2*sz_i32);

    //
    // Pass N/2 sample pairs through the last section.
    //

    alx = AE_LA64_PP(px);
    __Pragma( "loop_count min=1" );
    for ( n=0; n<(N>>1); n++ )
    {
      AE_LA32X2_IP( x01, alx, px );

      acc0 = AE_MULF32S_LH( cf_b10, x01 );

      AE_MULAAFD32S_HH_LL( acc0, cf_b21, st0 );
      AE_MULSSFD32S_HH_LL( acc0, cf_a21, st1 );

      AE_PKSRF32( st1, acc0, 1 );

      acc1 = AE_MULF32S_HL( cf_b21, st0 );

      AE_MULAAFD32S_HH_LL( acc1, cf_b10, x01 );
      AE_MULSSFD32S_HH_LL( acc1, cf_a21, st1 );

      AE_PKSRF32( st1, acc1, 1 );

      //y01 = AE_ROUND32X2F64SASYM(AE_SLAS64S(acc0), AE_SLAS64S(acc1));
      y01 = AE_TRUNCA32X2F64S(acc0,acc1,sh);
      AE_SA32X2_IP( y01, alr, pr );

      st0 = x01;
    }
    AE_SA64POS_FP(alr, pr);

    //
    // Save section's state.
    //

    AE_S32X2_I(st0, (ae_int32x2 *)state, 0*2*sz_i32);
    AE_S32X2_I(st1, (ae_int32x2 *)state, 1*2*sz_i32);
    AE_S32X2_I(st2, (ae_int32x2 *)state, 2*2*sz_i32);
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
} // bqriir32x32_df1_process()
