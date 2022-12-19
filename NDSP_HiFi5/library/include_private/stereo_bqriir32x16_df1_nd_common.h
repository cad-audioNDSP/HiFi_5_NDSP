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
#ifndef STEREO_BQRIIR32X16_DF1_ND_COMMON_H__
#define STEREO_BQRIIR32X16_DF1_ND_COMMON_H__

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"

/* Instance pointer validation number. */
#define STEREO_BQRIIR32X16_DF1_ND_MAGIC 0xd936b63a

/* Filter instance structure. */
typedef struct tag_stereo_bqriir32x16_df1_nd_t
{
  uint32_t       magic; // Instance pointer validation number
  int            M;     // Number of sections
  int16_t        gainl; // Total gain shift amount for the last biquad, left channel
  int16_t        gainr; // Total gain shift amount for the last biquad, right channel
  const int16_t *coef;  // Num/den coefs (Q14) and gain (Q15) for each biquad, both channels;
                        // coefficients are stored in the interleaved manner
  int32_t       *state; // 4 state elements per section, Q31, for both channels;
                        // state elements are stored in the interleaved manner
} stereo_bqriir32x16_df1_nd_t, *stereo_bqriir32x16_df1_nd_ptr_t;

#endif /* STEREO_BQRIIR32X16_DF1_ND_COMMON_H__ */
