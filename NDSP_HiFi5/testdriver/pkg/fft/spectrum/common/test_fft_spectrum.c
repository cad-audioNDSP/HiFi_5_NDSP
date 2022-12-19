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
 * Test procedures for FFT magnitude functions
 */

#include "test_fft_spectrum.h"

/* Perform all tests for FFT magnitude API functions. */
int main_spectrum( int phaseNum, int isFull, int isVerbose, int breakOnError, const TestDef_t *tbl )
{
  int ix, res;

  for ( ix=0, res=1; tbl[ix].seqFilename && ( res || !breakOnError ); ix++ )
  {
    if ( phaseNum == 0 || phaseNum == tbl[ix].phaseNum )
    {
      res &= ( 0 != TestEngRun( tbl[ix].target,
                                &tbl[ix].desc,
                                tbl[ix].seqFilename,
                                isFull, isVerbose, breakOnError,0 ) );
    }
  }

  return (res);

} /* main_fft_spectrum() */

/* Test Engine method: allocate vectors and load the data set: */
int loadFxn_fft_spectrum( tTestEngContext * context )
{
  int N, logN, nElem, res = 0;
  int isFract, baseFmt,outFmt;
  int mode = 0, bexp = -1;

  NASSERT( context );

  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  baseFmt = ( context->desc->fmt & ~FMT_CPLX );
  isFract = ( baseFmt == FMT_FRACT16 || baseFmt == FMT_FRACT32 );

  /* For a fractional variant, read the block_exponent argument. */
  if ( isFract && 1 != seqFileScanf( context->seqFile, "%d", &bexp ) )
  {
    printf( "loadFxn_fft_spectrum(): failed to read the block_exponent argument\n" );
  }
  /* Read the mode argument */
  else if ( 1 != seqFileScanf( context->seqFile, "%d", &mode ) )
  {
    printf( "loadFxn_fft_spectrum(): failed to read the mode argument\n" );
  }
  else res = 1;

  if ( !res ) return (0); else res = 0;
  
  N = context->args.dim[0]; logN = 30 - S_exp0_l( N );
  /* Allocate empty data vectors whenever:
   *   A) FFT size is out of domain
   *   B) FFT size is not a power of two 
   *   C) mode argument value is invalid
   *   D) block_exponent argument is invalid */
  if ( N>=2 && !(N&(N-1)) && ( mode == 0 || mode == 1 ) && bexp <= logN )
  {
    nElem = ( mode ? N/2+1 : N );
  }
  else
  {
    nElem = 0;
  }
  outFmt=(baseFmt==FMT_FRACT16) ?  FMT_FRACT32| FMT_REAL: baseFmt | FMT_REAL;

  /* Allocate input data vector memory. */
  if ( !vecAlloc( &context->dataSet.X, nElem, TE_ALIGN_NO, baseFmt | FMT_CPLX, 0 ) )
  {
    printf( "loadFxn_fft_spectrum(): failed to allocate vector X; "
            "fmt = 0x%02x, nElem = %d\n", (unsigned)(baseFmt|FMT_CPLX), nElem );
  }
  /* Allocate output/reference data vectors memory. */
  else if ( 3 != vecsAlloc( TE_ALIGN_NO, outFmt,
                            &context->dataSet.Z  , nElem,
                            &context->dataSet.Zlo, nElem,
                            &context->dataSet.Zhi, nElem, 0 ) )
  {
    printf( "loadFxn_fft_spectrum(): failed to allocate vectors Z/Zlo/Zhi; "
            "fmt = 0x%02x, nElem = %d\n", (unsigned)(baseFmt|FMT_REAL), nElem );
  }
  /* Load input/reference data vectors from the SEQ-file. */
  else if ( !seqFileReadVecs( context->seqFile,
                              &context->dataSet.X,
                              &context->dataSet.Zlo,
                              &context->dataSet.Zhi, 0 ) )
  {
    printf( "loadFxn_fft_spectrum(): failed to input/reference vectors data; "
            "fmt = 0x%02x, nElem = %d\n", (unsigned)baseFmt, nElem );
  }
  else res = 1;

  if ( !res )
  {
    if ( context->dataSet.X  .szBulk > 0 ) vecFree( &context->dataSet.X   );
    if ( context->dataSet.Z  .szBulk > 0 ) vecFree( &context->dataSet.Z   );
    if ( context->dataSet.Zlo.szBulk > 0 ) vecFree( &context->dataSet.Zlo );
    if ( context->dataSet.Zhi.szBulk > 0 ) vecFree( &context->dataSet.Zhi );
  }
  else
  {
    context->args.dim[fft_mag_aux_bexp] = bexp;
    context->args.dim[fft_mag_aux_mode] = mode;
  }

  return (res);

} /* loadFxn_fft_spectrum() */

/* Test Engine method: apply the function under test to test case data set. */
void processFxn_fft_spectrum( tTestEngContext * context )
{
  typedef void tFxn_fract( void * y,const void * x, 
                           int N, int block_exponent, int mode );
  typedef void tFxn_float(  void * y,const void * x,
                           int N, int mode );

  tTestEngTarget fxn;
  int baseFmt, isFract;
  int N, bexp, mode;
  const void * X;
  void * Z;

  NASSERT( context );

  baseFmt = ( context->desc->fmt & ~FMT_CPLX );
  isFract = ( baseFmt == FMT_FRACT16 || baseFmt == FMT_FRACT32 );

  N    = context->args.dim[0];
  bexp = context->args.dim[fft_mag_aux_bexp];
  mode = context->args.dim[fft_mag_aux_mode];
  fxn  = context->target.fut;

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  /* logging */
  te_vReportStd(context);

  if ( isFract )
  {
    ( *(tFxn_fract*)fxn )(Z, X,  N, bexp, mode );
  }
  else
  {
    ( *(tFxn_float*)fxn )(Z, X,  N, mode );
  }

} /* processFxn_fft_spectrum() */

