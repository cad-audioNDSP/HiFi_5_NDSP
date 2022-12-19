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
  NatureDSP Signal Processing Library. FIR part
    Polyphase FIR filters
    Internal declarations
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"

#ifndef __PPFIR_INTERNAL_H
#define __PPFIR_INTERNAL_H

/*
  Additional padding is used to push the delay lines to an another memory bank
  so this helps to avoid memory conflicts when concurrently loading
  successive delay elements
*/
#define DELAY_PAD_I32  (4)

/*
  Additional padding is used to push the delay line allocation to an another memory bank
  so this helps to avoid memory conflicts when concurrently loading
  filter coefficients and delay elements 
*/
#define ADD_PAD_CI32   (2)

/*
  Additional padding is used to push the delay line allocation to an another memory bank
  so this helps to avoid memory conflicts when concurrently loading
  filter coefficients and delay elements 
*/
#define ADD_PAD_F32    (4)

/*
 * Polyphase FIR filter kernel.
 * This function performs polyphase filtering for a block of M samples
 * specified via the input argument x[M] (type-1 polyphase decomposition).
 * This is accomplished in two steps:
 *   1. Insert M new samples into the current position (input argument p)
 *      of the circular delay line (input/output argument d) of size N*M. Current
 *      position p is updated by M sample positions in a round-robin fashion,
 *      and will be returned to the caller for subsequent calls of this function.
 *   2. Compute output samples of M subfilters, each of N taps (input argument 
 *      h[N*M]), and store results to the output argument y[M].  
 * Domain:
 *   rpp      Real data, real coefficients
 *   cpp      Complex data, real coefficients
 * Data type and precision
 *   pp32x32  32-bit fixed-point data, 32-bit fixed-point coefficients
 *   ppf      Floating point
 * Input:
 *   M        Number of subfilters (or the number of phases)
 *   N        Subfilter length
 *   lsh      (fixed-point variants only) Bi-directional saturating left shift amount
 *            to be applied to input samples prior to storing them into the delay line.
 *            For a negative (i.e. right) shift amount, results are asymmetrically
 *            rounded.
 *   x[M]     A block of input samples, assuming type-1 polyphase decomposition.
 *            That is, x[0] is the oldest sample in a block, and x[M-1] is the 
 *            newest one.
 *   h[N*M]   Subfilters' coefficients organized into an N-by-M matrix, where N is
 *            the number of subfilter taps, and M is the number of subfilters. That
 *            is, N coefficients of the m-th subfilter are stored to the m-th column:
 *            h[(0..N-1)*M+m], where coefficient h[0*M+m] corresponds to the oldest
 *            input sample in subfilter's delay line.
 * Input/Output:
 *   d[N*M]   Subfilters' delay lines organized into an N-by-M matrix, analogously to
 *            subfilter coefficients storage h[N*M]. That is, N-sample delay line of
 *            the m-th subfilter occupies the m-th column: d[(0..N-1)*M+m]. 
 *   p        Current position in the delay line, that is a pointer to one of N rows
 *            of the N-by-M matrix d[N*M]. In a single call of this function, this
 *            row gets filled with M input samples x[M], and the current position 
 *            switches to the next row. After the last row of the delay line is 
 *            filled with new data, the current position wraps around to the 0th row.
 *            The function returns the updated current position to be used for 
 *            subsequent calls.
 * Output:
 *   y[M]     Output samples of M subfilters. Output sample of 0th subfilter is stored
 *            to y[0], output sample of 1st subfilter - to y[1], and so forth, up to
 *            M-1st subfilter's output sample stored to y[M-1].
 * Restrictions:
 *   d,p,h,x  Must not overlap, and must be 16-byte aligned.
 *   N        4..24, multiple of 2
 *   M        32...640, multiple of 32
 */

int32_t * rppfir32x32(
                int32_t * restrict y,
                int32_t * restrict d,
                int32_t * restrict p,
          const int32_t * restrict h,
          const int32_t * restrict x,
          int M, int N, int lsh);
complex_fract32 * cppfir32x32(
                complex_fract32 * restrict y,
                complex_fract32 * restrict d,
                complex_fract32 * restrict p,
          const int32_t         * restrict h,
          const complex_fract32 * restrict x,
          int M, int N, int lsh);
float32_t * rppfirf(
                float32_t * restrict y,
                float32_t * restrict d,
                float32_t * restrict p,
          const float32_t * restrict h,
          const float32_t * restrict x,
          int M, int N);
complex_float * cppfirf(
                complex_float * restrict y,
                complex_float * restrict d,
                complex_float * restrict p,
          const float32_t     * restrict h,
          const complex_float * restrict x,
          int M, int N);

#endif /* __PPFIR_INTERNAL_ */
