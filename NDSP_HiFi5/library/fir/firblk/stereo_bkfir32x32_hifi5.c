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
    Real block FIR filter, 32x32-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Filter instance structure. */
#include "stereo_bkfir32x32_common.h"

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

/* Reserve memory for alignment. */
#define ALIGNED_SIZE( size, align ) \
      ( (size_t)(size) + (align) - 1 )

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
      (void*)( ( (uintptr_t)(addr) + ( (align) - 1 ) ) & ~( (align) - 1 ) )

#if 0
/* Allocation routine for filters. Returns: size of memory in bytes to be allocated */
size_t stereo_bkfir32x32_alloc( int M, int extIR )
{
    NASSERT(M > 0 && M % 4 == 0);
    return (ALIGNED_SIZE(sizeof(stereo_bkfir32x32_t), 4)
            + 2 * bkfir32x32_alloc(M, extIR)); // Sub structs
} /* stereo_bkfir32x32_alloc() */

/* Initialization for filters. Returns: handle to the object */
stereo_bkfir32x32_handle_t stereo_bkfir32x32_init( void * objmem, int M, int extIR, const int32_t * hl, const int32_t * hr )
{
    stereo_bkfir32x32_ptr_t   stereo_bkfir;
    void                    * ptr;
	size_t szbkfir;
	szbkfir = bkfir32x32_alloc(M,extIR);

    NASSERT(objmem && hl && hr);
    NASSERT_ALIGN(hl, 16);
    NASSERT_ALIGN(hr, 16);
    NASSERT(M > 0 && M % 4 == 0);

    // Partition the memory block
	ptr = objmem;
    stereo_bkfir = (stereo_bkfir32x32_ptr_t)ALIGNED_ADDR(ptr, 4);
	stereo_bkfir->magic = STEREO_BKFIR32X32_MAGIC;
	ptr = stereo_bkfir + 1;
    stereo_bkfir->bkfir_left_mem = ptr;
	ptr = (void*)(((uintptr_t)ptr)+szbkfir);
    stereo_bkfir->bkfir_right_mem = ptr;
	stereo_bkfir->bkfir_left  = bkfir32x32_init(stereo_bkfir->bkfir_left_mem ,M,extIR,hl);
	stereo_bkfir->bkfir_right = bkfir32x32_init(stereo_bkfir->bkfir_right_mem,M,extIR,hr);
	return stereo_bkfir;
} /* stereo_bkfir32x32_init() */
#else

#define sz_i32   sizeof(int32_t)

/* Allocation routine for filters. Returns: size of memory in bytes to be allocated */
size_t stereo_bkfir32x32_alloc( int M, int extIR )
{
  NASSERT(M > 0 && M % 4 == 0);

  return ( ALIGNED_SIZE( sizeof( stereo_bkfir32x32_t ), 4 )
           + // Delay line
           ALIGNED_SIZE( 2*( M + 12 )*sz_i32, 16 )
           + // Filter coefficients
           (extIR?0:ALIGNED_SIZE( 2*( M + 4 )*sz_i32, 16 ) ));
} /* stereo_bkfir32x32_alloc() */

/* Initialization for filters. Returns: handle to the object */
stereo_bkfir32x32_handle_t stereo_bkfir32x32_init( void * objmem, int M, int extIR, const int32_t * hl, const int32_t * hr )
{
  stereo_bkfir32x32_t * stereo_bkfir;
  void *                ptr;
  int32_t *             delLineLeft;
  int32_t *             delLineRight;
  int32_t *             coefLeft;
  int32_t *             coefRight;

  int m;

  NASSERT(objmem && hl && hr);
  NASSERT_ALIGN(hl, 16);
  NASSERT_ALIGN(hr, 16);
  NASSERT(M > 0 && M % 4 == 0);

  //
  // Partition the memory block
  //

  ptr          = objmem;
  stereo_bkfir = (stereo_bkfir32x32_ptr_t)ALIGNED_ADDR( ptr, 4 );
  ptr          = stereo_bkfir + 1;
  delLineLeft  = (int32_t *)ALIGNED_ADDR( ptr, 16 );
  ptr          = delLineLeft + M + 12;
  delLineRight = (int32_t *)ALIGNED_ADDR( ptr, 16 );
  ptr          = delLineRight + M + 12;
  if(extIR)
  {
    coefLeft  = (int32_t *) hl;
    coefRight = (int32_t *) hr;
  }
  else
  {
    coefLeft  = (int32_t *)ALIGNED_ADDR( ptr, 16 );
    ptr       = coefLeft + M + 4;
    coefRight = (int32_t *)ALIGNED_ADDR( ptr, 16 );
    ptr       = coefRight + M + 4;
  }

  ASSERT( (int8_t*)ptr - (int8_t*)objmem <= (int)stereo_bkfir32x32_alloc( M, extIR ) );

  //
  // Copy the filter coefficients and zero the delay line. Original impulse
  // response is padded with zeros: three zeros go before the first tap
  // (corresponds to the newest sample), one zero follows the last tap,
  // which matches the oldest sample. After that the order of filter
  // coefficients is reverted.
  //
  if (!extIR)
  {
    coefLeft [0] = 0;
    coefRight[0] = 0;

    for ( m=1; m<M+1; m++ )
    {
      coefLeft [m] = hl[M-m];
      coefRight[m] = hr[M-m];
    }

    for ( ; m<M+4; m++ )
    {
      coefLeft [m] = 0;
      coefRight[m] = 0;
    }
  }

  for ( m=0; m<M+12; m++ )
  {
    delLineLeft [m] = 0;
    delLineRight[m] = 0;
  }

  //
  // Initialize the filter instance.
  //

  //stereo_bkfir->magic          = MAGIC;
  stereo_bkfir->M              = M;
  stereo_bkfir->coefLeft       = coefLeft;
  stereo_bkfir->coefRight      = coefRight;
  stereo_bkfir->delayLineLeft  = delLineLeft;
  stereo_bkfir->delayLineRight = delLineRight;
  stereo_bkfir->delayLen       = M + 12;
  stereo_bkfir->delayPosLeft   = delLineLeft;
  stereo_bkfir->delayPosRight  = delLineRight;

  return ((stereo_bkfir32x32_handle_t)stereo_bkfir);
} /* stereo_bkfir32x32_init() */
#endif
