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
 * Test engine implementation
 * Collection of test case loading functions
 */

#include <stdio.h>
#include <string.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"
/* Test engine API. */
#include "testeng.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

/* Allocate vectors and load the data set: vector/scalar X (in), vector/scalar Z (out). */
static int te_loadFxn_vXvZ_helper( tTestEngContext * context, int fmtX, int nElemX, int fmtZ, int nElemZ );
/* Allocate vectors and load the data set:
 * vector/scalar X (in), vector/scalar Y (in), vector/scalar Z (out) */
static int te_loadFxn_vXvYvZ_helper( tTestEngContext * context, 
                                     int fmtX, int nElemX, 
                                     int fmtY, int nElemY, 
                                     int fmtZ, int nElemZ );
/* Allocate vectors and load the data set:
 * vector/scalar X (in), vector/scalar Y (in), vector/scalar U (in), vector/scalar Z (out) */
static int te_loadFxn_vXvYvUvZ_helper( tTestEngContext * context, int fmtX, int nElemX, 
                                       int fmtY, int nElemY, int fmtU, int nElemU, int fmtZ, int nElemZ );
/* Allocate vectors and load the data set:
 * vector/scalar X (in), vector/scalar Y (in), vector/scalar Z (out), vector/scalar W (out) */
static int te_loadFxn_vXvYvZvW_helper( tTestEngContext * context,
                                       int fmtX, int nElemX, int fmtY, int nElemY,
                                       int fmtZ, int nElemZ, int fmtW, int nElemW );
/* Allocate vectors and load the data set:
 * vector/scalar X (in), vector/scalar Y (in), vector/scalar Z (out), vector/scalar W (out) */
static int te_loadFxn_vXvZvW_helper( tTestEngContext * context, int fmtX, int nElemX, 
                                     int fmtZ, int nElemZ, int fmtW, int nElemW );
/* Allocate vectors and load the data set:
 * vector/scalar X (in), vector/scalar Y (in), vector/scalar U (in), vector/scalar Z (out), vector/scalar W (out) */
static int te_loadFxn_vXvYvUvZvW_helper( tTestEngContext * context,
                                         int fmtX, int nElemX, int fmtY, int nElemY, int fmtU, int nElemU,
                                         int fmtZ, int nElemZ, int fmtW, int nElemW );

/* Allocate vectors and load the data set:
 * vector X (in), vector Z (out) */
int te_loadFxn_vXvZ( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmt;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;

  fmt = context->desc->fmt;

  return ( te_loadFxn_vXvZ_helper( context, fmt, nElem, fmt, nElem ) );

} /* te_loadFxn_vXvZ() */

/* Allocate vectors and load the data set:
 * vector X (in), vector Z (out, single precision floating-point) */
int te_loadFxn_vXvZf( tTestEngContext * context )
{
  int nElem, fmt;

  ASSERT( context && context->seqFile );

  nElem = MAX( 0, context->args.dim[0] )*
          MAX( 0, context->args.dim[1] )*
          MAX( 0, context->args.dim[2] )*
          MAX( 0, context->args.dim[3] );
  fmt = context->desc->fmt;

  return ( te_loadFxn_vXvZ_helper( context, fmt, nElem, FMT_FLOAT32, nElem ) );

} /* te_loadFxn_vXvZf() */

/* Allocate vectors and load the data set:
 * vector X (in), vector Z (out, double precision floating-point) */
int te_loadFxn_vXvZd( tTestEngContext * context )
{
  int nElem, fmt;

  ASSERT( context && context->seqFile );

  nElem = MAX( 0, context->args.dim[0] )*
          MAX( 0, context->args.dim[1] )*
          MAX( 0, context->args.dim[2] )*
          MAX( 0, context->args.dim[3] );
  fmt = context->desc->fmt;

  return ( te_loadFxn_vXvZ_helper( context, fmt, nElem, FMT_FLOAT64, nElem ) );

} /* te_loadFxn_vXvZd() */

/* Allocate vectors and load the data set:
 * vector X (in), vector Z (out) */
int te_loadFxn_vXvZh( tTestEngContext * context )
{
  int nElem, fmt;

  ASSERT( context && context->seqFile );

  nElem = MAX( 0, context->args.dim[0] )*
          MAX( 0, context->args.dim[1] )*
          MAX( 0, context->args.dim[2] )*
          MAX( 0, context->args.dim[3] );
  fmt = context->desc->fmt;

  return ( te_loadFxn_vXvZ_helper( context, fmt, nElem, FMT_FLOAT16, nElem ) );

} /* te_loadFxn_vXvZh() */

/* Allocate vectors and load the data set:
 * vector int64 X (in), vector Z (out) */
int te_loadFxn_vX64vZ( tTestEngContext * context )
{
  int nElem, fmt;

  ASSERT( context && context->seqFile );

  nElem = MAX( 0, context->args.dim[0] )*
          MAX( 0, context->args.dim[1] )*
          MAX( 0, context->args.dim[2] )*
          MAX( 0, context->args.dim[3] );
  fmt = context->desc->fmt;

  return ( te_loadFxn_vXvZ_helper( context, FMT_INT64, nElem, fmt, nElem ) );

} /* te_loadFxn_vX64vZ() */

/* Allocate vectors and load the data set:
 * * vector complex X (in), vector Z (out) */
int te_loadFxn_vXcvZ( tTestEngContext * context )
{
  int nElem, fmt;

  ASSERT( context && context->seqFile );

  nElem = MAX( 0, context->args.dim[0] )*
          MAX( 0, context->args.dim[1] )*
          MAX( 0, context->args.dim[2] )*
          MAX( 0, context->args.dim[3] );
  fmt = context->desc->fmt;

  return ( te_loadFxn_vXvZ_helper( context, (fmt & ~FMT_DOMAIN_MASK) | FMT_CPLX, nElem, fmt, nElem ) );

} /* te_loadFxn_vXcvZ() */

/* Allocate vectors and load the data set:
 * vector X (in), scalar Z (out) */
int te_loadFxn_vXsZ( tTestEngContext * context )
{
  int M, N, L;
  int nElemX, nElemZ, fmt;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemX = M*N*L;
  nElemZ = L;

  fmt = context->desc->fmt;

  return ( te_loadFxn_vXvZ_helper( context, fmt, nElemX, fmt, nElemZ ) );

} /* te_loadFxn_vXsZ() */

/* Allocate vectors and load the data set:
 * vector X (in), scalar int16 Z (out) */
int te_loadFxn_vXsZ16( tTestEngContext * context )
{
  int M, N, L;
  int nElemX, nElemZ, fmtX, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  fmtX = context->desc->fmt;
  fmtZ = FMT_INT16;

  nElemX = M*N*L;
  nElemZ = L;

  return ( te_loadFxn_vXvZ_helper( context, fmtX, nElemX, fmtZ, nElemZ ) );

} /* te_loadFxn_vXsZ16() */

/* Allocate vectors and load the data set:
 * vector X (in), scalar int32 Z (out) */
int te_loadFxn_vXsZ32( tTestEngContext * context )
{
  int M, N, L;
  int nElemX, nElemZ, fmtX, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  fmtX = context->desc->fmt;
  fmtZ = FMT_INT32;

  nElemX = M*N*L;
  nElemZ = L;

  return ( te_loadFxn_vXvZ_helper( context, fmtX, nElemX, fmtZ, nElemZ ) );

} /* te_loadFxn_vXsZ32() */

/* Allocate vectors and load the data set:
 * vector int32 X (in), scalar Z (out) */
int te_loadFxn_vX32sZ( tTestEngContext * context )
{
  int M, N, L;
  int nElemX, nElemZ, fmtX, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  fmtX = FMT_INT32;
  fmtZ = context->desc->fmt;

  nElemX = M*N*L;
  nElemZ = L;

  return ( te_loadFxn_vXvZ_helper( context, fmtX, nElemX, fmtZ, nElemZ ) );

} /* te_loadFxn_vX32sZ() */

/* Allocate vectors and load the data set:
 * 32-bit vector X (in), 16-bit vector Z (out) */
int te_loadFxn_vX32vZ16( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmtX, fmtZ;

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;

  ASSERT( context );

  fmtX = ( ( context->desc->fmt & FMT_DOMAIN_MASK ) | FMT_FRACT32 );
  fmtZ = ( ( context->desc->fmt & FMT_DOMAIN_MASK ) | FMT_FRACT16 );

  return ( te_loadFxn_vXvZ_helper( context, fmtX, nElem, fmtZ, nElem ) );

} /* te_loadFxn_vX32vZ16() */

/* Allocate vectors and load the data set:
 * 16-bit vector X (in), 32-bit vector Z (out) */
int te_loadFxn_vX16vZ32( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmtX, fmtZ;

  ASSERT( context );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;

  fmtX = ( ( context->desc->fmt & FMT_DOMAIN_MASK ) | FMT_FRACT16 );
  fmtZ = ( ( context->desc->fmt & FMT_DOMAIN_MASK ) | FMT_FRACT32 );

  return ( te_loadFxn_vXvZ_helper( context, fmtX, nElem, fmtZ, nElem ) );

} /* te_loadFxn_vX16vZ32() */

/* Allocate vectors and load the data set:
 * real/complex vector X (in), real vector Z (out) */
int te_loadFxn_vXvZr( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmtX, fmtZ;

  ASSERT( context );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;

  fmtX = context->desc->fmt;
  fmtZ = ( ( context->desc->fmt & FMT_DTYPE_MASK ) | FMT_REAL );

  return ( te_loadFxn_vXvZ_helper( context, fmtX, nElem, fmtZ, nElem ) );

} /* te_loadFxn_vXvZr() */

/* Allocate vectors and load the data set:
 * real/complex vector X (in), complex vector Z (out) */
int te_loadFxn_vXvZc( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmtX, fmtZ;

  ASSERT( context );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;

  fmtX = context->desc->fmt;
  fmtZ = ( ( context->desc->fmt & FMT_DTYPE_MASK ) | FMT_CPLX );

  return ( te_loadFxn_vXvZ_helper( context, fmtX, nElem, fmtZ, nElem ) );

} /* te_loadFxn_vXvZc() */

/* Allocate vectors and load the data set:
 * vector X (in), vector Y (in), vector Z (out) */
int te_loadFxn_vXvYvZ( tTestEngContext * context )
{
  int M, N, L;
  int fmt, nElem;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;
  fmt = context->desc->fmt;

  return ( te_loadFxn_vXvYvZ_helper( context, fmt, nElem, fmt, nElem, fmt, nElem ) );

} /* te_loadFxn_vXvYvZ() */

/* Allocate vectors and load the data set:
 * vector X (in), vector Y (in), scalar Z (out) */
int te_loadFxn_vXvYsZ( tTestEngContext * context )
{
  int M, N, L;
  int nElemXY, nElemZ, fmt;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemXY = M*N*L;
  nElemZ = L;
  fmt = context->desc->fmt;

  return ( te_loadFxn_vXvYvZ_helper( context, fmt, nElemXY, fmt, nElemXY, fmt, nElemZ ) );

} /* te_loadFxn_vXvYsZ() */

/* Allocate vectors and load the data set:
 * vector X (in), scalar Y (in), vector Z (out) */
int te_loadFxn_vXsYvZ( tTestEngContext * context )
{
  int M, N, L;
  int fmt, nElemXZ, nElemY;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemXZ = M*N*L;
  nElemY  = L;

  fmt = context->desc->fmt;

  return ( te_loadFxn_vXvYvZ_helper( context, fmt, nElemXZ, fmt, nElemY, fmt, nElemXZ ) );

} /* te_loadFxn_vXsYvZ() */

/* Allocate vectors and load the data set:
 * vector X (in), scalar int32 Y (in), vector Z (out) */
int te_loadFxn_vXsY32vZ( tTestEngContext * context )
{
  int M, N, L;
  int fmtXZ, fmtY, nElemXZ, nElemY;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemXZ = M*N*L;
  nElemY  = L;

  fmtXZ = context->desc->fmt;
  fmtY  = (context->desc->fmt & FMT_DOMAIN_MASK) | FMT_INT32;

  return ( te_loadFxn_vXvYvZ_helper( context, fmtXZ, nElemXZ, fmtY, nElemY, fmtXZ, nElemXZ ) );

} /* te_loadFxn_vXsY32vZ() */

/* Allocate vectors and load the data set:
 * vector int32 X (in), scalar int32 Y (in), vector Z (out) */
int te_loadFxn_vX32sY32vZ( tTestEngContext * context )
{
  int M, N, L;
  int fmtXY, fmtZ, nElemXZ, nElemY;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemXZ = M*N*L;
  nElemY  = L;

  fmtXY = (context->desc->fmt & FMT_DOMAIN_MASK) | FMT_INT32;
  fmtZ  = context->desc->fmt;

  return ( te_loadFxn_vXvYvZ_helper( context, fmtXY, nElemXZ, fmtXY, nElemY, fmtZ, nElemXZ ) );

} /* te_loadFxn_vX32sY32vZ() */

/* Allocate vectors and load the data set:
 * vector X (in), scalar int32 Y (in), vector int32 Z (out) */
int te_loadFxn_vXsY32vZ32( tTestEngContext * context )
{
  int M, N, L;
  int fmtX, fmtYZ, nElemXZ, nElemY;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemXZ = M*N*L;
  nElemY  = L;

  fmtX  = context->desc->fmt;
  fmtYZ = (context->desc->fmt & FMT_DOMAIN_MASK) | FMT_INT32;

  return ( te_loadFxn_vXvYvZ_helper( context, fmtX, nElemXZ, fmtYZ, nElemY, fmtYZ, nElemXZ ) );

} /* te_loadFxn_vXsY32vZ32() */

/* Allocate vectors and load the data set:
 * vector X (in), scalar Y (in), vector Z (out) */
int te_loadFxn_vXsYrvZ(tTestEngContext * context)
{
    int M, N, L;
    int fmtXZ, fmtY, nElemXZ, nElemY;

    ASSERT(context && context->seqFile);

    N = MAX(0, context->args.dim[0]);
    M = MAX(0, context->args.dim[1]);
    L = MAX(0, context->args.dim[2]);

    nElemXZ = M*N*L;
    nElemY = L;

    fmtXZ = 
    fmtY = context->desc->fmt;

    return (te_loadFxn_vXvYvZ_helper(context, fmtXZ, nElemXZ, fmtY, nElemY, fmtXZ, nElemXZ));

} /* te_loadFxn_vXsYrvZ() */

/* Allocate vectors and load the data set:
 * vector X (in), real scalar Y (in), scalar int16 Z (out) */
int te_loadFxn_vXsYr_sZ16( tTestEngContext * context )  
{
  int M, N, L;
  int nElemX, nElemYZ, fmtX, fmtY, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemX  = M*N*L;
  nElemYZ = L;

  fmtX = context->desc->fmt;
  fmtY = FMT_INT32;
  fmtZ = ( ( context->desc->fmt & FMT_DOMAIN_MASK ) | FMT_FRACT16 );

  return ( te_loadFxn_vXvYvZ_helper( context, fmtX, nElemX, fmtY, nElemYZ, fmtZ, nElemYZ ) );

} /* te_loadFxn_vXsYr_sZ16() */

/* Allocate vectors and load the data set:
 * vector X (in), real scalar Y (in), scalar int32 Z (out) */
int te_loadFxn_vXsYr_sZ32( tTestEngContext * context )  
{
  int M, N, L;
  int nElemX, nElemYZ, fmtX, fmtY, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemX  = M*N*L;
  nElemYZ = L;

  fmtX = context->desc->fmt;
  fmtY = FMT_INT32;
  fmtZ = ( ( context->desc->fmt & FMT_DOMAIN_MASK ) | FMT_FRACT32 );

  return ( te_loadFxn_vXvYvZ_helper( context, fmtX, nElemX, fmtY, nElemYZ, fmtZ, nElemYZ ) );

} /* te_loadFxn_vXsYr_sZ32() */

/* Allocate vectors and load the data set:
 * real/complex vector X (in), real scalar Y (in), real vector Z (out) */
int te_loadFxn_vXsYr_vZr( tTestEngContext * context )
{
  int M, N, L;
  int nElemXZ, nElemY, fmtX, fmtY, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemXZ = M*N*L;
  nElemY  = L;

  fmtX = context->desc->fmt;
  fmtY = FMT_INT32;
  fmtZ = ( ( context->desc->fmt & FMT_DTYPE_MASK ) | FMT_REAL );

  return ( te_loadFxn_vXvYvZ_helper( context, fmtX, nElemXZ, fmtY, nElemY, fmtZ, nElemXZ ) );

} /* te_loadFxn_vXsYr_vZr() */

/* Allocate vectors and load the data set:
 * real/complex vector X (in), real scalar Y (in), complex vector Z (out) */
int te_loadFxn_vXsYr_vZc( tTestEngContext * context )
{
  int M, N, L;
  int nElemXZ, nElemY, fmtX, fmtY, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemXZ = M*N*L;
  nElemY  = L;

  fmtX = ( ( context->desc->fmt & FMT_DTYPE_MASK ) | FMT_REAL );
  fmtY = FMT_INT32;
  fmtZ = ( ( context->desc->fmt & FMT_DTYPE_MASK ) | FMT_CPLX );

  return ( te_loadFxn_vXvYvZ_helper( context, fmtX, nElemXZ, fmtY, nElemY, fmtZ, nElemXZ ) );

} /* te_loadFxn_vXsYr_vZc() */

/* Allocate vectors and load the data set:
 * vector X (in), vector Y (in), float32 vector Z (out) */
int te_loadFxn_vXvYvZf( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmtXY, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;

  fmtXY = context->desc->fmt;
  fmtZ  = ( ( context->desc->fmt & FMT_DOMAIN_MASK ) | FMT_FLOAT32 );

  return ( te_loadFxn_vXvYvZ_helper( context, fmtXY, nElem, fmtXY, nElem, fmtZ, nElem ) );

} /* te_loadFxn_vXvYvZf() */

/* Allocate vectors and load the data set:
 * vector X (in), vector Y (in), float64 vector Z (out) */
int te_loadFxn_vXvYvZd( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmtXY, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;

  fmtXY = context->desc->fmt;
  fmtZ  = ( ( context->desc->fmt & FMT_DOMAIN_MASK ) | FMT_FLOAT64 );

  return ( te_loadFxn_vXvYvZ_helper( context, fmtXY, nElem, fmtXY, nElem, fmtZ, nElem ) );

} /* te_loadFxn_vXvYvZf() */

/* Allocate vectors and load the data set:
 * vector int64 X (in), vector Y (in), vector Z (out) */
int te_loadFxn_vX64vYvZ( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmtX, fmtYZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;

  fmtX  = ( ( context->desc->fmt & FMT_DOMAIN_MASK ) | FMT_INT64 );
  fmtYZ = context->desc->fmt;

  return ( te_loadFxn_vXvYvZ_helper( context, fmtX, nElem, fmtYZ, nElem, fmtYZ, nElem ) );

} /* te_loadFxn_vX64vYvZ() */

/* Allocate vectors and load the data set:
 * vector X (in), scalar int32 Y (in), scalar int64 Z (out) */
int te_loadFxn_vXsY32sZ64( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmtX, fmtY, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;

  fmtX = context->desc->fmt;
  fmtY = (context->desc->fmt & FMT_DOMAIN_MASK) | FMT_INT32;
  fmtZ = (context->desc->fmt & FMT_DOMAIN_MASK) | FMT_INT64;

  return ( te_loadFxn_vXvYvZ_helper( context, fmtX, nElem, fmtY, 1, fmtZ, 1 ) );

} /* te_loadFxn_vXsY32sZ64() */

/* Allocate vectors and load the data set:
 * vector X (in), vector Y (in), scalar int32 Z (out) */
int te_loadFxn_vXvYsZ32( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmtXY, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;

  fmtXY = context->desc->fmt;
  fmtZ  = (context->desc->fmt & FMT_DOMAIN_MASK) | FMT_INT32;

  return ( te_loadFxn_vXvYvZ_helper( context, fmtXY, nElem, fmtXY, nElem, fmtZ, 1 ) );

} /* te_loadFxn_vXvYsZ32() */

/* Allocate vectors and load the data set:
 * vector X (in), vector Y (in), scalar int64 Z (out) */
int te_loadFxn_vXvYsZ64( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmtXY, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;

  fmtXY = context->desc->fmt;
  fmtZ  = (context->desc->fmt & FMT_DOMAIN_MASK) | FMT_INT64;

  return ( te_loadFxn_vXvYvZ_helper( context, fmtXY, nElem, fmtXY, nElem, fmtZ, 1 ) );

} /* te_loadFxn_vXvYsZ64() */

/* Allocate vectors and load the data set:
 * vector X (in), vector int16 Y (in), scalar int64 Z (out) */
int te_loadFxn_vXvY16sZ64( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmtX, fmtY, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;

  fmtX = context->desc->fmt;
  fmtY = (context->desc->fmt & FMT_DOMAIN_MASK) | FMT_INT16;
  fmtZ = (context->desc->fmt & FMT_DOMAIN_MASK) | FMT_INT64;

  return ( te_loadFxn_vXvYvZ_helper( context, fmtX, nElem, fmtY, nElem, fmtZ, 1 ) );

} /* te_loadFxn_vXvY16sZ64() */

/* Allocate vectors and load the data set:
 * vector X (in), vector int32 Y (in), scalar int64 Z (out) */
int te_loadFxn_vXvY32sZ64( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmtX, fmtY, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;

  fmtX = context->desc->fmt;
  fmtY = (context->desc->fmt & FMT_DOMAIN_MASK) | FMT_INT32;
  fmtZ = (context->desc->fmt & FMT_DOMAIN_MASK) | FMT_INT64;

  return ( te_loadFxn_vXvYvZ_helper( context, fmtX, nElem, fmtY, nElem, fmtZ, 1 ) );

} /* te_loadFxn_vXvY32sZ64() */

/* Allocate vectors and load the data set:
 * vector X (in), vector Z (out), vector W (out) */
int te_loadFxn_vXvZvW( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmt;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;
  fmt = context->desc->fmt;

  return ( te_loadFxn_vXvZvW_helper( context, fmt, nElem, fmt, nElem, fmt, nElem ) );

} /* te_loadFxn_vXvZvW() */

/* Allocate vectors and load the data set:
 * vector X (in), vector Y (in), vector Z (out), vector W (out) */
int te_loadFxn_vXvYvZvW( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmt;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;
  fmt = context->desc->fmt;

  return ( te_loadFxn_vXvYvZvW_helper( context, fmt, nElem, fmt, nElem, fmt, nElem, fmt, nElem ) );

} /* te_loadFxn_vXvYvZvW() */

/* Allocate vectors and load the data set:
 * scalar X (in), scalar Y (in), vector Z (out), scalar W (out) */
int te_loadFxn_sXsYvZsW( tTestEngContext * context )
{
  int M, N, L;
  int nElemXYW, nElemZ, fmtXYW, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemXYW = L;
  nElemZ = M*N*L;

  fmtXYW = FMT_INT32;
  fmtZ = context->desc->fmt;

  return ( te_loadFxn_vXvYvZvW_helper( context,
                                       fmtXYW, nElemXYW, fmtXYW, nElemXYW, 
                                       fmtZ, nElemZ, fmtXYW, nElemXYW ) );

} /* te_loadFxn_sXsYvZsW() */

/* Allocate vectors and load the data set:
 * vector X (in), vector Y (in), real scalar U (in), scalar int32 Z (out) */
int te_loadFxn_vXvYsUr_sZ32( tTestEngContext * context )
{
  int M, N, L;
  int nElemXY, nElemUZ, fmtXY, fmtU, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemXY = M*N*L;
  nElemUZ = L;

  fmtXY = context->desc->fmt;
  fmtU  = FMT_INT32;
  fmtZ = ((context->desc->fmt & FMT_DOMAIN_MASK) | FMT_FRACT32 );

  return ( te_loadFxn_vXvYvUvZ_helper( context, fmtXY, nElemXY, fmtXY, 
                                      nElemXY, fmtU, nElemUZ, fmtZ, nElemUZ ) );

} /* te_loadFxn_vXvYsUr_sZ32() */

/* Allocate vectors and load the data set:
 * vector X (in), vector Y (in), real scalar U (in), complex vector Z (out) */
int te_loadFxn_vXvYsUr_vZc( tTestEngContext * context )
{
  int M, N, L;
  int nElemXYZ, nElemU, fmtXY, fmtU, fmtZ;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemXYZ = M*N*L;
  nElemU   = L;

  fmtXY = context->desc->fmt;
  fmtU  = FMT_INT32;
  fmtZ  = ( ( context->desc->fmt & FMT_DTYPE_MASK ) | FMT_CPLX );

  return ( te_loadFxn_vXvYvUvZ_helper( context, fmtXY, nElemXYZ, fmtXY, nElemXYZ,
                                       fmtU, nElemU, fmtZ, nElemXYZ ) );

} /* te_loadFxn_vXvYsUr_vZc() */

/* Allocate vectors and load the data set:
 * vector X (in), real scalars Y,U (in), vector Z (out) */
int te_loadFxn_vXsYrsUr_vZ( tTestEngContext * context )
{
  int M, N, L;
  int nElemXZ, nElemYU, fmtXZ, fmtYU;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemXZ = M*N*L;
  nElemYU = L;

  fmtXZ = context->desc->fmt;
  fmtYU = FMT_INT32;

  return ( te_loadFxn_vXvYvUvZ_helper( context, fmtXZ, nElemXZ, fmtYU, nElemYU,
                                       fmtYU, nElemYU, fmtXZ, nElemXZ ) );

} /* te_loadFxn_vXsYrsUr_vZ() */

/* Allocate vectors and load the data set:
 * vector X (in), real scalars Y.U (in), vector Z (out), real scale W (out) */
int te_loadFxn_vXsYrsUr_vZsWr( tTestEngContext * context )
{
  int M, N, L;
  int fmtXZ, fmtYUW, nElemXZ, nElemYUW;

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElemXZ  = M*N*L;
  nElemYUW = L;

  fmtXZ  = context->desc->fmt;
  switch(context->desc->fmt & FMT_DTYPE_MASK)
  {
  case FMT_FLOAT32: fmtYUW = FMT_FLOAT32; break;
  case FMT_FLOAT16: fmtYUW = FMT_FLOAT16; break;
  default:          fmtYUW = FMT_INT32;   break;
  }

  return ( te_loadFxn_vXvYvUvZvW_helper( context, fmtXZ, nElemXZ, fmtYUW, nElemYUW, 
                                         fmtYUW, nElemYUW, fmtXZ, nElemXZ, fmtYUW, nElemYUW ) );

} /* te_loadFxn_vXsYrsUr_vZsWr() */

/* Allocate vectors and load the data set: vector/scalar X (in), vector/scalar Z (out). */
int te_loadFxn_vXvZ_helper( tTestEngContext * context, int fmtX, int nElemX, int fmtZ, int nElemZ )
{
  int res = 0;

  ASSERT( context && context->seqFile );

  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  /* Allocate input data vector memory. */
  if ( !vecAlloc( &context->dataSet.X, nElemX,
                   context->desc->isAligned, fmtX, 0 ) )
  {
    printf( "te_loadFxn_vXvZ_helper(): failed to allocate vector X; "
            "fmtX = 0x%02x, nElemX = %d\n", (unsigned)fmtX, nElemX );
  }
  /* Allocate output/reference vectors memory. */
  else if ( 3 != vecsAlloc( context->desc->isAligned, fmtZ,
                           &context->dataSet.Z  , nElemZ,
                           &context->dataSet.Zlo, nElemZ,
                           &context->dataSet.Zhi, nElemZ, 0 ) )
  {
    printf( "te_loadFxn_vXvZ_helper(): failed to allocate vectors Z/Zlo/Zhi; "
            "fmtZ = 0x%02x, nElemZ = %d\n", (unsigned)fmtZ, nElemZ );
  }
  /* Load vectors data from the SEQ-file. */
  else if ( !seqFileReadVecs( context->seqFile,
                             &context->dataSet.X,
                             &context->dataSet.Zlo,
                             &context->dataSet.Zhi, 0 ) )
  {
    printf( "te_loadFxn_vXvZ_helper(): failed to read vectors data; "
            "fmtX = 0x%02x, nElemX = %d, fmtZ = 0x%02x, nElemZ = %d\n",
            (unsigned)fmtX, nElemX, (unsigned)fmtZ, nElemZ );
  }
  else
  {
    res = 1;
  }

  /* Free vectors data if failed. */
  if ( !res )
  {
    if ( context->dataSet.X  .szBulk ) vecFree( &context->dataSet.X   );
    if ( context->dataSet.Z  .szBulk ) vecFree( &context->dataSet.Z   );
    if ( context->dataSet.Zlo.szBulk ) vecFree( &context->dataSet.Zlo );
    if ( context->dataSet.Zhi.szBulk ) vecFree( &context->dataSet.Zhi );
  }

  return (res);

} /* te_loadFxn_vXvZ_helper() */

/* Allocate vectors and load the data set:
 * vector/scalar X (in), vector/scalar Y (in), vector/scalar Z (out) */
int te_loadFxn_vXvYvZ_helper( tTestEngContext * context, 
                              int fmtX, int nElemX, 
                              int fmtY, int nElemY, 
                              int fmtZ, int nElemZ )  
{
  int res = 0;

  ASSERT( context && context->seqFile );

  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  /* Allocate input data vector X. */
  if ( !vecAlloc( &context->dataSet.X, nElemX, context->desc->isAligned, fmtX, 0 ) )
  {
    printf( "te_loadFxn_vXvYvZ_helper(): failed to allocate vector X; "
            "fmtX = 0x%02x, nElemX = %d\n", (unsigned)fmtX, nElemX );
  }
  /* Allocate input data vector Y. */
  else if ( !vecAlloc( &context->dataSet.Y, nElemY, context->desc->isAligned, fmtY, 0 ) )
  {
    printf( "te_loadFxn_vXvYvZ_helper(): failed to allocate vector Y; "
            "fmtY = 0x%02x, nElemY = %d\n", (unsigned)fmtY, nElemY );
  }
  /* Allocate output/reference data vectors memory. */
  else if ( 3 != vecsAlloc( context->desc->isAligned, fmtZ,
                           &context->dataSet.Z  , nElemZ,
                           &context->dataSet.Zlo, nElemZ,
                           &context->dataSet.Zhi, nElemZ, 0 ) )
  {
    printf( "te_loadFxn_vXvYvZ_helper(): failed to allocate vectors Z/Zlo/Zhi; "
            "fmtZ = 0x%02x, nElemZ = %d\n", (unsigned)fmtZ, nElemZ );
  }
  /* Load vectors data from the SEQ-file. */
  else if ( !seqFileReadVecs( context->seqFile,
                             &context->dataSet.X,
                             &context->dataSet.Y,
                             &context->dataSet.Zlo,
                             &context->dataSet.Zhi, 0 ) )
  {
    printf( "te_loadFxn_vXvYvZ_helper(): failed to read vectors data; "
            "fmtX = 0x%02x, nElemX = %d, fmtY = 0x%02x, nElemY = %d, fmtZ = 0x%02x, nElemZ = %d\n",
            (unsigned)fmtX, nElemX, (unsigned)fmtY, nElemY, (unsigned)fmtZ, nElemZ );
  }
  else
  {
    res = 1;
  }

  /* Free vectors data if failed. */
  if ( !res )
  {
    if ( context->dataSet.X  .szBulk ) vecFree( &context->dataSet.X   );
    if ( context->dataSet.Y  .szBulk ) vecFree( &context->dataSet.Y   );
    if ( context->dataSet.Z  .szBulk ) vecFree( &context->dataSet.Z   );
    if ( context->dataSet.Zlo.szBulk ) vecFree( &context->dataSet.Zlo );
    if ( context->dataSet.Zhi.szBulk ) vecFree( &context->dataSet.Zhi );
  }

  return (res);

} /* te_loadFxn_vXvYvZ_helper() */

/* Allocate vectors and load the data set:
 * vector/scalar X (in), vector/scalar Y (in), vector/scalar U (in), vector/scalar Z (out) */
int te_loadFxn_vXvYvUvZ_helper( tTestEngContext * context, int fmtX, int nElemX, 
                                int fmtY, int nElemY, int fmtU, int nElemU, int fmtZ, int nElemZ )
{
  int res = 0;

  ASSERT( context && context->seqFile );

  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  /* Allocate input data vector X. */
  if ( !vecAlloc( &context->dataSet.X, nElemX, context->desc->isAligned, fmtX, 0 ) )
  {
    printf( "te_loadFxn_vXvYvUvZ_helper(): failed to allocate vector X; "
            "fmtX = 0x%02x, nElemX = %d\n", (unsigned)fmtX, nElemX );
  }
  /* Allocate input data vector Y. */
  else if ( !vecAlloc( &context->dataSet.Y, nElemY, context->desc->isAligned, fmtY, 0 ) )
  {
    printf( "te_loadFxn_vXvYvUvZ_helper(): failed to allocate vector Y; "
            "fmtY = 0x%02x, nElemY = %d\n", (unsigned)fmtY, nElemY );
  }
  /* Allocate input data vector U. */
  else if ( !vecAlloc( &context->dataSet.U, nElemU, context->desc->isAligned, fmtU, 0 ) )
  {
    printf( "te_loadFxn_vXvYvUvZ_helper(): failed to allocate vector U; "
            "fmtU = 0x%02x, nElemU = %d\n", (unsigned)fmtU, nElemU );
  }
  /* Allocate output/reference data vectors Z/Zlo/Zhi memory. */
  else if ( 3 != vecsAlloc( context->desc->isAligned, fmtZ,
                           &context->dataSet.Z  , nElemZ,
                           &context->dataSet.Zlo, nElemZ,
                           &context->dataSet.Zhi, nElemZ, 0 ) )
  {
    printf( "te_loadFxn_vXvYvUvZ_helper(): failed to allocate vectors Z/Zlo/Zhi; "
            "fmtZ = 0x%02x, nElemZ = %d\n", (unsigned)fmtZ, nElemZ );
  }
  /* Load vectors data from the SEQ-file. */
  else if ( !seqFileReadVecs( context->seqFile,
                             &context->dataSet.X,
                             &context->dataSet.Y,
                             &context->dataSet.U,
                             &context->dataSet.Zlo,
                             &context->dataSet.Zhi, 0 ) )
  {
    printf( "te_loadFxn_vXvYvUvZ_helper(): failed to read vectors data; "
            "fmtX = 0x%02x, nElemX = %d, fmtY = 0x%02x, nElemY = %d, fmtU = 0x%02x, nElemU = %d, "
            "fmtZ = 0x%02x, nElemZ = %d\n",
            (unsigned)fmtX, nElemX, (unsigned)fmtY, nElemY, 
            (unsigned)fmtU, nElemU, (unsigned)fmtZ, nElemZ );
  }
  else
  {
    res = 1;
  }

  /* Free vectors data if failed. */
  if ( !res )
  {
    if ( context->dataSet.X  .szBulk ) vecFree( &context->dataSet.X   );
    if ( context->dataSet.Y  .szBulk ) vecFree( &context->dataSet.Y   );
    if ( context->dataSet.U  .szBulk ) vecFree( &context->dataSet.U   );
    if ( context->dataSet.Z  .szBulk ) vecFree( &context->dataSet.Z   );
    if ( context->dataSet.Zlo.szBulk ) vecFree( &context->dataSet.Zlo );
    if ( context->dataSet.Zhi.szBulk ) vecFree( &context->dataSet.Zhi );
  }

  return (res);

} /* te_loadFxn_vXvYvUvZ_helper() */

/* Allocate vectors and load the data set:
 * vector/scalar X (in), vector/scalar Y (in), vector/scalar Z (out), vector/scalar W (out) */
int te_loadFxn_vXvYvZvW_helper( tTestEngContext * context,
                                int fmtX, int nElemX, int fmtY, int nElemY,
                                int fmtZ, int nElemZ, int fmtW, int nElemW )
{
  int res = 0;

  ASSERT( context && context->seqFile );

  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  /* Allocate input data vector X. */
  if ( !vecAlloc( &context->dataSet.X, nElemX, context->desc->isAligned, fmtX, 0 ) )
  {
    printf( "te_loadFxn_vXvYvZvW_helper(): failed to allocate vector X; "
            "fmtX = 0x%02x, nElemX = %d\n", (unsigned)fmtX, nElemX );
  }
  /* Allocate input data vector Y. */
  else if ( !vecAlloc( &context->dataSet.Y, nElemY, context->desc->isAligned, fmtY, 0 ) )
  {
    printf( "te_loadFxn_vXvYvZvW_helper(): failed to allocate vector Y; "
            "fmtY = 0x%02x, nElemY = %d\n", (unsigned)fmtY, nElemY );
  }
  /* Allocate output/reference data vectors Z/Zlo/Zhi memory. */
  else if ( 3 != vecsAlloc( context->desc->isAligned, fmtZ,
                           &context->dataSet.Z  , nElemZ,
                           &context->dataSet.Zlo, nElemZ,
                           &context->dataSet.Zhi, nElemZ, 0 ) )
  {
    printf( "te_loadFxn_vXvYvZvW_helper(): failed to allocate vectors Z/Zlo/Zhi; "
            "fmtZ = 0x%02x, nElemZ = %d\n", (unsigned)fmtZ, nElemZ );
  }
  /* Allocate output/reference data vectors W/Wlo/Whi memory. */
  else if ( 3 != vecsAlloc( context->desc->isAligned, fmtW,
                           &context->dataSet.W  , nElemW,
                           &context->dataSet.Wlo, nElemW,
                           &context->dataSet.Whi, nElemW, 0 ) )
  {
    printf( "te_loadFxn_vXvYvZvW_helper(): failed to allocate vectors W/Wlo/Whi; "
            "fmtW = 0x%02x, nElemW = %d\n", (unsigned)fmtW, nElemW );
  }
  /* Load vectors data from the SEQ-file. */
  else if ( !seqFileReadVecs( context->seqFile,
                             &context->dataSet.X,
                             &context->dataSet.Y,
                             &context->dataSet.Zlo,
                             &context->dataSet.Zhi,
                             &context->dataSet.Wlo,
                             &context->dataSet.Whi, 0 ) )
  {
    printf( "te_loadFxn_vXvYvZvW_helper(): failed to read vectors data; "
            "fmtX = 0x%02x, nElemX = %d, fmtY = 0x%02x, nElemY = %d, "
            "fmtZ = 0x%02x, nElemZ = %d, fmtW = 0x%02x, nElemW = %d\n",
            (unsigned)fmtX, nElemX, (unsigned)fmtY, nElemY, 
            (unsigned)fmtZ, nElemZ, (unsigned)fmtW, nElemW );
  }
  else
  {
    res = 1;
  }

  /* Free vectors data if failed. */
  if ( !res )
  {
    if ( context->dataSet.X  .szBulk ) vecFree( &context->dataSet.X   );
    if ( context->dataSet.Y  .szBulk ) vecFree( &context->dataSet.Y   );
    if ( context->dataSet.Z  .szBulk ) vecFree( &context->dataSet.Z   );
    if ( context->dataSet.Zlo.szBulk ) vecFree( &context->dataSet.Zlo );
    if ( context->dataSet.Zhi.szBulk ) vecFree( &context->dataSet.Zhi );
    if ( context->dataSet.W  .szBulk ) vecFree( &context->dataSet.W   );
    if ( context->dataSet.Wlo.szBulk ) vecFree( &context->dataSet.Wlo );
    if ( context->dataSet.Whi.szBulk ) vecFree( &context->dataSet.Whi );
  }

  return (res);

} /* te_loadFxn_vXvYvZvW_helper() */

/* Allocate vectors and load the data set:
 * vector/scalar X (in), vector/scalar Y (in), vector/scalar Z (out), vector/scalar W (out) */
int te_loadFxn_vXvZvW_helper( tTestEngContext * context, int fmtX, int nElemX, 
                              int fmtZ, int nElemZ, int fmtW, int nElemW )
{
  int res = 0;

  ASSERT( context && context->seqFile );

  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  /* Allocate input data vector X. */
  if ( !vecAlloc( &context->dataSet.X, nElemX, context->desc->isAligned, fmtX, 0 ) )
  {
    printf( "te_loadFxn_vXvZvW_helper(): failed to allocate vector X; "
            "fmtX = 0x%02x, nElemX = %d\n", (unsigned)fmtX, nElemX );
  }
  /* Allocate output/reference data vectors Z/Zlo/Zhi memory. */
  else if ( 3 != vecsAlloc( context->desc->isAligned, fmtZ,
                           &context->dataSet.Z  , nElemZ,
                           &context->dataSet.Zlo, nElemZ,
                           &context->dataSet.Zhi, nElemZ, 0 ) )
  {
    printf( "te_loadFxn_vXvZvW_helper(): failed to allocate vectors Z/Zlo/Zhi; "
            "fmtZ = 0x%02x, nElemZ = %d\n", (unsigned)fmtZ, nElemZ );
  }
  /* Allocate output/reference data vectors W/Wlo/Whi memory. */
  else if ( 3 != vecsAlloc( context->desc->isAligned, fmtW,
                           &context->dataSet.W  , nElemW,
                           &context->dataSet.Wlo, nElemW,
                           &context->dataSet.Whi, nElemW, 0 ) )
  {
    printf( "te_loadFxn_vXvZvW_helper(): failed to allocate vectors W/Wlo/Whi; "
            "fmtW = 0x%02x, nElemW = %d\n", (unsigned)fmtW, nElemW );
  }
  /* Load vectors data from the SEQ-file. */
  else if ( !seqFileReadVecs( context->seqFile,
                             &context->dataSet.X,
                             &context->dataSet.Zlo,
                             &context->dataSet.Zhi,
                             &context->dataSet.Wlo,
                             &context->dataSet.Whi, 0 ) )
  {
    printf( "te_loadFxn_vXvZvW_helper(): failed to read vectors data; "
            "fmtX = 0x%02x, nElemX = %d, fmtZ = 0x%02x, nElemZ = %d, "
            "fmtW = 0x%02x, nElemW = %d\n", (unsigned)fmtX, nElemX,
            (unsigned)fmtZ, nElemZ, (unsigned)fmtW, nElemW );
  }
  else
  {
    res = 1;
  }

  /* Free vectors data if failed. */
  if ( !res )
  {
    if ( context->dataSet.X  .szBulk ) vecFree( &context->dataSet.X   );
    if ( context->dataSet.Z  .szBulk ) vecFree( &context->dataSet.Z   );
    if ( context->dataSet.Zlo.szBulk ) vecFree( &context->dataSet.Zlo );
    if ( context->dataSet.Zhi.szBulk ) vecFree( &context->dataSet.Zhi );
    if ( context->dataSet.W  .szBulk ) vecFree( &context->dataSet.W   );
    if ( context->dataSet.Wlo.szBulk ) vecFree( &context->dataSet.Wlo );
    if ( context->dataSet.Whi.szBulk ) vecFree( &context->dataSet.Whi );
  }

  return (res);

} /* te_loadFxn_vXvZvW_helper() */

/* Allocate vectors and load the data set:
 * vector/scalar X (in), vector/scalar Y (in), vector/scalar U (in), vector/scalar Z (out), vector/scalar W (out) */
int te_loadFxn_vXvYvUvZvW_helper( tTestEngContext * context,
                                  int fmtX, int nElemX, int fmtY, int nElemY, int fmtU, int nElemU,
                                  int fmtZ, int nElemZ, int fmtW, int nElemW )
{
  int res = 0;

  ASSERT( context && context->seqFile );

  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  /* Allocate input data vector X. */
  if ( !vecAlloc( &context->dataSet.X, nElemX, context->desc->isAligned, fmtX, 0 ) )
  {
    printf( "te_loadFxn_vXvYvUvZvW_helper(): failed to allocate vector X; "
            "fmtX = 0x%02x, nElemX = %d\n", (unsigned)fmtX, nElemX );
  }
  /* Allocate input data vector Y. */
  else if ( !vecAlloc( &context->dataSet.Y, nElemY, context->desc->isAligned, fmtY, 0 ) )
  {
    printf( "te_loadFxn_vXvYvUvZvW_helper(): failed to allocate vector Y; "
            "fmtY = 0x%02x, nElemY = %d\n", (unsigned)fmtY, nElemY );
  }
  /* Allocate input data vector U. */
  else if ( !vecAlloc( &context->dataSet.U, nElemU, context->desc->isAligned, fmtU, 0 ) )
  {
    printf( "te_loadFxn_vXvYvUvZvW_helper(): failed to allocate vector U; "
            "fmtU = 0x%02x, nElemU = %d\n", (unsigned)fmtU, nElemU );
  }
  /* Allocate output/reference data vectors Z/Zlo/Zhi memory. */
  else if ( 3 != vecsAlloc( context->desc->isAligned, fmtZ,
                           &context->dataSet.Z  , nElemZ,
                           &context->dataSet.Zlo, nElemZ,
                           &context->dataSet.Zhi, nElemZ, 0 ) )
  {
    printf( "te_loadFxn_vXvYvUvZvW_helper(): failed to allocate vectors Z/Zlo/Zhi; "
            "fmtZ = 0x%02x, nElemZ = %d\n", (unsigned)fmtZ, nElemZ );
  }
  /* Allocate output/reference data vectors W/Wlo/Whi memory. */
  else if ( 3 != vecsAlloc( context->desc->isAligned, fmtW,
                           &context->dataSet.W  , nElemW,
                           &context->dataSet.Wlo, nElemW,
                           &context->dataSet.Whi, nElemW, 0 ) )
  {
    printf( "te_loadFxn_vXvYvUvZvW_helper(): failed to allocate vectors W/Wlo/Whi; "
            "fmtW = 0x%02x, nElemW = %d\n", (unsigned)fmtW, nElemW );
  }
  /* Load vectors data from the SEQ-file. */
  else if ( !seqFileReadVecs( context->seqFile,
                             &context->dataSet.X,
                             &context->dataSet.Y,
                             &context->dataSet.U,
                             &context->dataSet.Zlo,
                             &context->dataSet.Zhi,
                             &context->dataSet.Wlo,
                             &context->dataSet.Whi, 0 ) )
  {
    printf( "te_loadFxn_vXvYvUvZvW_helper(): failed to read vectors data; "
            "fmtX = 0x%02x, nElemX = %d, fmtY = 0x%02x, nElemY = %d, fmtU = 0x%02x, nElemU = %d, "
            "fmtZ = 0x%02x, nElemZ = %d, fmtW = 0x%02x, nElemW = %d\n",
            (unsigned)fmtX, nElemX, (unsigned)fmtY, nElemY, (unsigned)fmtU, nElemU, 
            (unsigned)fmtZ, nElemZ, (unsigned)fmtW, nElemW );
  }
  else
  {
    res = 1;
  }

  /* Free vectors data if failed. */
  if ( !res )
  {
    if ( context->dataSet.X  .szBulk ) vecFree( &context->dataSet.X   );
    if ( context->dataSet.Y  .szBulk ) vecFree( &context->dataSet.Y   );
    if ( context->dataSet.U  .szBulk ) vecFree( &context->dataSet.U   );
    if ( context->dataSet.Z  .szBulk ) vecFree( &context->dataSet.Z   );
    if ( context->dataSet.Zlo.szBulk ) vecFree( &context->dataSet.Zlo );
    if ( context->dataSet.Zhi.szBulk ) vecFree( &context->dataSet.Zhi );
    if ( context->dataSet.W  .szBulk ) vecFree( &context->dataSet.W   );
    if ( context->dataSet.Wlo.szBulk ) vecFree( &context->dataSet.Wlo );
    if ( context->dataSet.Whi.szBulk ) vecFree( &context->dataSet.Whi );
  }

  return (res);

} /* te_loadFxn_vXvYvUvZvW_helper() */

/* Allocate vectors and load the data set:
 * vector X (in), vector Y (in), vector Z (out, double precision) */
int te_loadFxn_vXvYvZdp( tTestEngContext * context )
{
  int M, N, L;
  int fmt, nElem;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );

  nElem = M*N*L;
  fmt = context->desc->fmt;

  return ( te_loadFxn_vXvYvZ_helper( context, fmt, nElem, fmt, nElem, FMT_REAL|FMT_FLOAT64, nElem ) );

} 
