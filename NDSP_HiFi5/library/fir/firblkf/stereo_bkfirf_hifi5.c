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
    Real block FIR filter, floating point
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "common_fpu.h"
/* Filter instance structure. */
#include "stereo_bkfirf_common.h"

/*-------------------------------------------------------------------------
  Real FIR filter.
  Computes a real FIR filter (direct-form) using IR stored in vector h. The 
  real data input is stored in vector x. The filter output result is stored 
  in vector y. The filter calculates N output samples using M coefficients 
  and requires last M-1 samples in the delay line.
  NOTE: 
  1. User application is not responsible for management of delay lines
  2. User has an option to set IR externally or copy from original location 
     (i.e. from the slower constant memory). In the first case, user is 
     responsible for right alignment, ordering and zero padding of filter 
     coefficients - usually array is composed from zeroes (left padding), 
     reverted IR and right zero padding.


  Precision: 
  16x16    16-bit data, 16-bit coefficients, 16-bit outputs. Ordinary variant 
           and stereo
  32x16    32-bit data, 16-bit coefficients, 32-bit outputs
  32x32    32-bit data, 32-bit coefficients, 32-bit outputs. Ordinary variant 
           and stereo
  32x32ep  the same as above but using 72-bit accumulator for intermediate 
           computations
  f        floating point. Ordinary variant and stereo

  Input:
  x[N*S]   input samples, Q31, Q15, floating point
  h[M]     filter coefficients in normal order, Q31, Q15, floating point
  hl[M]    for stereo filters: filter coefficients for left channel
  hr[M]    for stereo filters: filter coefficients for right channel
  N        length of sample block, should be a multiple of 4
  M        length of filter, should be a multiple of 4
  extIR    if zero, IR is copied from original location, otherwise not
           but user should keep alignment, order of coefficients 
           and zero padding requirements shown below
  S        1 for ordinary (single channel) filters, 2 - for stereo variant
  
  Output:
  y[N*S]   output samples, Q31, Q15, floating point

  Alignment, ordering and zero padding for external IR  (extIR!=0)
  ------------------------+----------+--------------+--------------+----------------
  Function                |Alignment,|Left zero     |   Coefficient| Right zero 
                          | bytes    |padding, bytes|   order      | padding, bytes
  ------------------------+----------+--------------+--------------+----------------
  bkfir16x16_init         |    16    |      2       |  inverted    |  6
  bkfir32x16_init         |    16    |      2       |  inverted    |  6
  bkfir32x32_init         |    16    |      4       |  inverted    |  12
  bkfir32x32ep_init       |    16    |      4       |  inverted    |  12
  bkfirf_init             |    16    |      4       |  inverted    |  12
  stereo_bkfir16x16_init  |    16    |      2       |  inverted    |  6
  stereo_bkfir32x32_init  |    16    |      4       |  inverted    |  12
  stereo_bkfirf_init      |    16    |      4       |  inverted    |  12
  ------------------------+----------+--------------+--------------+----------------

  Restrictions:
  x, y     should not be overlapping
  x, h     aligned on a 16-bytes boundary
  N, M     multiples of 4 
-------------------------------------------------------------------------*/

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(size_t,stereo_bkfirf_alloc,( int M , int extIR))
DISCARD_FUN(stereo_bkfirf_handle_t,stereo_bkfirf_init,( void * objmem, int M, int extIR, const float32_t * hl, const float32_t * hr ))
#else

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
      ( (size_t)(size) + (align) - 1 )

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
      (void*)( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

#define sz_f32      sizeof(float32_t)

/* Calculate the memory block size for an FIR filter with given attributes. */
size_t stereo_bkfirf_alloc( int M , int extIR)
{
#if !(HAVE_VFPU) //for scalar code
    NASSERT(M > 0 && M % 4 == 0);
    return (ALIGNED_SIZE(sizeof(stereo_bkfirf_t), 4)
            + 2 * bkfirf_alloc(M, extIR)); // Sub structs
#else
    NASSERT(M > 0 && M % 4 == 0);
    return (ALIGNED_SIZE(sizeof(stereo_bkfirf_t), 4)
        + // Delay line
        2*(ALIGNED_SIZE((M + 8)*sz_f32, 16))
        + // Filter coefficients
        2*(extIR ? 0 : ALIGNED_SIZE((M + 8)*sz_f32, 16)));
#endif
} /* stereo_bkfirf_alloc() */

/* Initialize the filter structure. The delay line is zeroed. */
stereo_bkfirf_handle_t stereo_bkfirf_init( void * objmem, int M, int extIR, const float32_t * hl, const float32_t * hr )
{
#if !(HAVE_VFPU) //scalar version
    stereo_bkfirf_ptr_t   stereo_bkfir;
    void                * ptr;
	size_t szbkfir;
	szbkfir = bkfirf_alloc(M,extIR);

    NASSERT(objmem && hl && hr);
    NASSERT_ALIGN(hl, 16);
    NASSERT_ALIGN(hr, 16);
    NASSERT(M > 0 && M % 4 == 0);

    // Partition the memory block
	ptr = objmem;
    stereo_bkfir = (stereo_bkfirf_ptr_t)ALIGNED_ADDR(ptr, 4);
	stereo_bkfir->magic = STEREO_BKFIRF_MAGIC;
	ptr = stereo_bkfir + 1;
    stereo_bkfir->bkfir_left_mem = ptr;
	ptr = (void*)(((uintptr_t)ptr)+szbkfir);
    stereo_bkfir->bkfir_right_mem = ptr;
	stereo_bkfir->bkfir_left  = bkfirf_init(stereo_bkfir->bkfir_left_mem ,M,extIR,hl);
	stereo_bkfir->bkfir_right = bkfirf_init(stereo_bkfir->bkfir_right_mem,M,extIR,hr);
	return stereo_bkfir;
#else
    stereo_bkfirf_t  * bkfir;
    void      * ptr;
    float32_t * delLine;
    float32_t * coefl;
    float32_t * coefr;

    int m;

    NASSERT(objmem && hl && hr);
    NASSERT_ALIGN(hl, 16);
    NASSERT_ALIGN(hr, 16);
    NASSERT(M > 0 && M % 4 == 0);

    //
    // Partition the memory block
    //
    ptr     = objmem;
    bkfir   = (stereo_bkfirf_ptr_t)ALIGNED_ADDR(ptr, 4);
    ptr     = bkfir + 1;
    delLine = (float32_t *)ALIGNED_ADDR(ptr, 16);
    ptr     = delLine + 2*(M+8);
    if (extIR)
    {
        coefl = (float32_t *)hl;
        coefr = (float32_t *)hr;
    }
    else
    {
        coefl = (float32_t *)ALIGNED_ADDR(ptr, 16);
        ptr  = coefl + M + 4;
        coefr = (float32_t *)ALIGNED_ADDR(ptr, 16);
        ptr  = coefr + M + 4;
    }
    ASSERT((int8_t*)ptr - (int8_t*)objmem <= (int)stereo_bkfirf_alloc(M, extIR));

    //
    // Copy the filter coefficients and zero the delay line. Original impulse
    // response is padded with zeros: three zeros go before the first tap
    // (corresponds to the newest sample), one zero follows the last tap,
    // which matches the oldest sample. After that the order of filter
    // coefficients is reverted.
    //
    if (extIR == 0)
    {
        coefl[0] = 0;
        coefr[0] = 0;

        for (m = 1; m < M + 1; m++)
        {
            coefl[m] = hl[M - m];
            coefr[m] = hr[M - m];
        }

        for (; m < M + 4; m++)
        {
            coefl[m] = 0;
            coefr[m] = 0;
        }
    }

    for (m = 0; m < 2*(M+8); m++)
    {
        delLine[m] = 0;
    }

    //
    // Initialize the filter instance.
    //
    bkfir->magic     = STEREO_BKFIRF_MAGIC;
    bkfir->M         = M;
    bkfir->coefLeft  = coefl;
    bkfir->coefRight = coefr;
    //if (extIR == 1) 
    //{
    //    bkfir->coefRight = coefr;
    //}
    bkfir->delayLine = delLine;
    bkfir->delayLen  = 2*(M +8);
    bkfir->delayPos  = delLine;
    return (bkfir);

#endif
} /* stereo_bkfirf_init() */
#endif
