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
 * Test procedures for real DCT functions.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API. */
#include LIBRARY_HEADER(fft)
/* Test engine API. */
#include "testeng.h"
/* Test engine extension for DCT. */
#include "../common/testeng_dct.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"

/* Suppress Visual C warnings on +/-HUGE_VALF macros. */
#ifdef COMPILER_MSVC
#pragma warning(disable:4056)
#pragma warning(disable:4756)
#endif

/* Test data file processing for a standalone test. Applies the target DCT function
 * to test data, compares the DCT output with reference data and estimates the SINAD. */
static void fileProc_stnd( tTestEngContext * context )
{
  tTestEngContext_dct           * context_dct;
  tTestEngFrameProcFxn_fft_stnd * frameProcFxn;

  tVec inVec, outVec, refVec; /* In/out/reference vectors related to test data files. */
  tVec xVec, yVec;            /* In/out vectors for the target DCT routine. */
  FILE *fIn = 0, *fRef = 0;
  int res = 0;
  int lenIn=0, lenOut=0;
  int fmtX=0, fmtY=0;
  int n, N;

  NASSERT( context && context->desc );
  NASSERT( context->target.fut && context->target.handle );

  context_dct = (tTestEngContext_dct*)context->target.handle;

  /* DCT size */
  N = context->args.dim[0];

  /* Select length and format (real/complex) for all vectors. */
  switch ( context->args.caseType )
  {
  case TE_FFT_TESTCASE_FORWARD:
    lenIn  = N;
    lenOut = N;
    if (context->desc->extraParam & TE_MDCT) lenIn = 2*N;
    break;
  case TE_FFT_TESTCASE_INVERSE:
    lenIn  = N;
    lenOut = N;
    if (context->desc->extraParam & TE_MDCT) lenOut = 2*N;
    break;
  default:
    ASSERT( 0 );
  }
  fmtX   = FMT_REAL;
  fmtY   = FMT_REAL;

  /* Preset an invalid SINAD for the test result. */
  *vecGetElem_fl32( &context->dataSet.Z, 0 ) = -HUGE_VALF;

  if ( context->isVerbose )
  {
    if ( context->desc->extraParam & TE_FFT_OPT_SCALE_METH )
    {
      printf( "scale_mtd %d ", context_dct->scale_method );
    }

    printf( "%-40s  ", context_dct->fRefName );
  }

  memset( &inVec , 0, sizeof(inVec ) );
  memset( &outVec, 0, sizeof(outVec) );
  memset( &refVec, 0, sizeof(refVec) );
    {
        char * fname;
        const char *dir;
		dir = getVectorsDir(context->isFull);
        fname=(char*)malloc(strlen(context_dct->fInName)+strlen(dir)+2+5);
        sprintf(fname,"%s/%s",dir,context_dct->fInName);
        fIn = fopen(fname, "rb");
        free(fname);
        fname=(char*)malloc(strlen(context_dct->fRefName)+strlen(dir)+2+5);
        sprintf(fname,"%s/%s",dir,context_dct->fRefName);
        fRef = fopen(fname, "rb");
        free(fname);
    }

  /* Allocate vectors for in/out/reference data. */
  if ( !vecAlloc( &inVec , lenIn , TE_ALIGN_YES, FMT_FRACT16|fmtX, 0 ) ||
       !vecAlloc( &outVec, lenOut, TE_ALIGN_YES, FMT_FLOAT64|fmtY, 0 ) ||
       !vecAlloc( &refVec, lenOut, TE_ALIGN_YES, FMT_FLOAT64|fmtY, 0 ) )
  {
    printf( "fileProc_stnd(): failed to allocate in/out/ref vectors, N=%d\n", N );
  }
  /* Open input data file. */
  else if (fIn==NULL )
  {
    printf( "fileProc_stnd(): failed to open %s for reading\n", context_dct->fInName );
  }
  /* Open reference data file. */
  else if (fRef == NULL )
  {
    printf( "fileProc_stnd(): failed to open %s for reading\n", context_dct->fRefName );
  }
  else res = 1;

  /*
   * Process the input file frame-by-frame.
   */

  if ( res )
  {
    int baseFmt = ( context->desc->fmt & FMT_DTYPE_MASK );
    int isAligned = context->desc->isAligned;

    fract16   * pin  = (fract16   *)vecGetElem( &inVec , 0 );
    float64_t * pout = (float64_t *)vecGetElem( &outVec, 0 );
    float64_t * pref = (float64_t *)vecGetElem( &refVec, 0 );

    float32_t sinadAvg, sinadMin = HUGE_VALF;
    float64_t errSum = 0, refSum = 0;

    int efbMin = 32;

    memset( &xVec, 0, sizeof(xVec) );
    memset( &yVec, 0, sizeof(yVec) );

    context_dct->frameCnt = 0;

    frameProcFxn = ( (tTestEngFrameProcFxn_fft_stnd*)context->target.fut );

    /* Read 16-bit complex samples from the input file. */
    while ( ( n = fread( pin, inVec.szElem, ( lenIn & ~1 ), fIn ) ) > 0 )
    {
      /* Augment the (last) incomplete frame with zeros. */
      memset( (uint8_t*)pin + n*inVec.szElem, 0, ( ( lenIn & ~1 ) - n )*inVec.szElem );
      /* Zero the output frame. */
      memset( pout, 0, lenOut*outVec.szElem );

      /* Allocate in/out buffers for the target DCT routine. */
      if ( ( !xVec.szBulk && !vecAlloc( &xVec, lenIn , isAligned, baseFmt | fmtX, 0 ) ) ||
           ( !yVec.szBulk && !vecAlloc( &yVec, lenOut, isAligned, baseFmt | fmtY, 0 ) ) )
      {
        printf( "fileProc_stnd(): failed to allocate xVec/yVec, "
                "frameCnt=%d, N=%d\n", context_dct->frameCnt, N );
        res = 0; break;
      }

      /* Use a proprietary frame processing function to perform the target DCT. */
      if ( !( res = frameProcFxn( context, pin, pout, &xVec, &yVec ) ) ) break;

      /* When in unaligned mode, in/out vectors should be periodically reallocated to try various
       * address offsets. */
      if ( !( ++context_dct->frameCnt % XY_VEC_REALLOC_PERIOD ) && !isAligned )
      {
        if ( !( res = vecsFree( &xVec, &yVec, 0 ) ) )
        {
          printf( "fileProc_stnd(): vecsFree() failed for xVec/yVec, "
                  "frameCnt=%d, N=%d\n", context_dct->frameCnt, N );
          break;
        }

        memset( &xVec, 0, sizeof(xVec) );
        memset( &yVec, 0, sizeof(yVec) );
      }

      /* Read double precision reference samples */
      if ( (int)fread( pref, refVec.szElem, lenOut, fRef ) < lenOut )
      {
        printf( "fileProc_stnd(): failed to read reference data, "
                "frameCnt=%d, N=%d\n", context_dct->frameCnt, lenOut );
        res = 0; break;
      }

      /* Estimate the SINAD ratio for the current frame, and update error/reference
       * power sums for the whole file SINAD. */
      {
        float64_t refMax = 0, errMax = 0;
        float64_t p, q;
        int m, M;

        M = ( ( fmtY & FMT_CPLX ) ? 2*lenOut : lenOut );

        for ( p=q=0, m=0; m<M; m++ )
        {
          float64_t err = pout[m] - pref[m];

          /* |signal|^2 */
          p += pref[m]*pref[m];
          /* |noise+distortion|^2 */
          q += err*err;

          refMax = MAX( refMax, fabs( pref[m] ) );
          errMax = MAX( errMax, fabs( err     ) );
        }

        refSum += p;
        errSum += q;

        if ( p>0 )
        {
          sinadMin = MIN( sinadMin, (float32_t)(p/q) );
          efbMin   = MAX( 0, MIN( efbMin, L_sub_ll( (int)logb(refMax), (int)logb(errMax) ) ) );
        }
      }
    }

    /*
     * Finalize the min SINAD estimation and print a summary.
     */

    if ( res )
    {
      sinadMin = 10.f*STDLIB_MATH(log10f)( sinadMin );

      /* Set the test result for test verification. */
      *vecGetElem_fl32( &context->dataSet.Z, 0 ) = sinadMin;

      if ( context->isVerbose )
      {
        sinadAvg = ( refSum > 0 ? 10.f*STDLIB_MATH(log10f)( (float32_t)(refSum/errSum) ) : -HUGE_VALF );

        if ( ( baseFmt & FMT_DTYPE_MASK ) == FMT_FRACT16 || ( baseFmt & FMT_DTYPE_MASK ) == FMT_FRACT32 )
        {
          printf( "Error-Free Bits %2d  SINAD min %5.1f dB  avg %5.1f dB  ",
                  efbMin, sinadMin, sinadAvg );
        }
        else
        {
          printf( "SINAD min %5.1f dB  avg %5.1f dB  ", sinadMin, sinadAvg );
        }
      }
    }
  }

  /*
   * Close files and free vectors.
   */

  if ( fIn  ) fclose( fIn  );
  if ( fRef ) fclose( fRef );

  if ( inVec .szBulk > 0 ) vecFree( &inVec  );
  if ( outVec.szBulk > 0 ) vecFree( &outVec );
  if ( refVec.szBulk > 0 ) vecFree( &refVec );
  if ( xVec  .szBulk > 0 ) vecFree( &xVec   );
  if ( yVec  .szBulk > 0 ) vecFree( &yVec   );

} /* fileProc_stnd() */

#define sz_fr16c    sizeof(complex_fract16)
#define sz_fl64c    sizeof(complex_double)

/* Initializer for a test description structure. */
#define TEST_DESC_STND( fmt, extraParam, align,fwddct,invdct )   {{ (fmt),(extraParam),NULL,TE_DIM_NUM_1,(align), \
                                                     &te_create_dct, &te_destroy_dct, \
                                                     &te_load_dct, &fileProc_stnd },(tFrwTransFxn)fwddct,(tInvTransFxn)invdct}
#define TEST_DESC_DCT2D( fmt, extraParam, align,fwddct,invdct )   {{ (fmt),(extraParam),NULL,TE_DIM_NUM_2,(align), \
                                                     &te_create_dct, &te_destroy_dct, \
                                                     &te_loadFxn_dct2d, &te_processFxn_dct2d },(tFrwTransFxn)fwddct,(tInvTransFxn)invdct}

static const struct tagTestDef
{
  int                phaseNum;
  const char *       seqFilename;
  tTestEngTarget     target;
  tTestEngDesc_dct   desc;

} testTbl[] = 
{
  { 1, "dct1/dct32_16x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT16, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_16x16 , NULL   ) },
  { 1, "dct1/dct64_16x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT16, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_16x16 , NULL   ) },
#if 0 // HiFi3/3z API
  { 1, "dct1/dct32_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_24x24 , NULL   ) },
  { 1, "dct1/dct64_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_24x24 , NULL   ) },
#endif
  { 1, "dct1/dct32_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_32x16 , NULL   ) },
  { 1, "dct1/dct64_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_32x16 , NULL   ) },
  { 1, "dct1/dct32_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_32x32 , NULL   ) },
  { 1, "dct1/dct64_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct_32x32 , NULL   ) },
#if 0 // HiFi3/3z API
  { 1, "dct1/dct4_32_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct1/dct4_64_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct1/dct4_128_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct1/dct4_256_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct1/dct4_512_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct1/dct4_32_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct1/dct4_64_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct1/dct4_128_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct1/dct4_256_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
  { 1, "dct1/dct4_512_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_24x24, NULL   ) },
#endif
  { 1, "dct1/dct4_32_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x16, NULL   ) },
  { 1, "dct1/dct4_64_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x16, NULL   ) },
  { 1, "dct1/dct4_128_32x16.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x16, NULL   ) },
  { 1, "dct1/dct4_256_32x16.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x16, NULL   ) },
  { 1, "dct1/dct4_512_32x16.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x16, NULL   ) },

  { 1, "dct1/dct4_32_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x32, NULL   ) },
  { 1, "dct1/dct4_64_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x32, NULL   ) },
  { 1, "dct1/dct4_128_32x32.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x32, NULL   ) },
  { 1, "dct1/dct4_256_32x32.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x32, NULL   ) },
  { 1, "dct1/dct4_512_32x32.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, dct4_32x32, NULL   ) },

#if 0 // HiFi3/3z API
  { 1, "dct1/mdct32_24x24.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_24x24, NULL    ) },
  { 1, "dct1/mdct64_24x24.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_24x24, NULL    ) },
  { 1, "dct1/mdct128_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_24x24, NULL    ) },
  { 1, "dct1/mdct256_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_24x24, NULL    ) },
  { 1, "dct1/mdct512_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_24x24, NULL    ) },
#endif
  { 1, "dct1/mdct32_32x16.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x16, NULL    ) },
  { 1, "dct1/mdct64_32x16.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x16, NULL    ) },
  { 1, "dct1/mdct128_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x16, NULL    ) },
  { 1, "dct1/mdct256_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x16, NULL    ) },
  { 1, "dct1/mdct512_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x16, NULL    ) },

  { 1, "dct1/mdct32_32x32.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x32, NULL    ) },
  { 1, "dct1/mdct64_32x32.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x32, NULL    ) },
  { 1, "dct1/mdct128_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x32, NULL    ) },
  { 1, "dct1/mdct256_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x32, NULL    ) },
  { 1, "dct1/mdct512_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, mdct_32x32, NULL    ) },
#if 0 // HiFi3/3z API
  { 1, "dct1/imdct32_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_24x24   ) },
  { 1, "dct1/imdct64_24x24.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_24x24   ) },
  { 1, "dct1/imdct128_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_24x24   ) },
  { 1, "dct1/imdct256_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_24x24   ) },
  { 1, "dct1/imdct512_24x24.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_24x24   ) },
#endif
  { 1, "dct1/imdct32_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x16   ) },
  { 1, "dct1/imdct64_32x16.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x16   ) },
  { 1, "dct1/imdct128_32x16.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x16   ) },
  { 1, "dct1/imdct256_32x16.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x16   ) },
  { 1, "dct1/imdct512_32x16.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x16   ) },

  { 1, "dct1/imdct32_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x32   ) },
  { 1, "dct1/imdct64_32x32.seq" , (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x32   ) },
  { 1, "dct1/imdct128_32x32.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x32   ) },
  { 1, "dct1/imdct256_32x32.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x32   ) },
  { 1, "dct1/imdct512_32x32.seq", (tTestEngTarget)te_frameProc_stnd_scl_dct , TEST_DESC_STND( FMT_FRACT32, ATTR_FAST|TE_FFT_OPT_SCALE_METH|TE_MDCT, TE_ALIGN_YES, NULL, imdct_32x32   ) },

  { 1, "dct1/dct2d_8x16.seq" , (tTestEngTarget)dct2d_8x16 , TEST_DESC_DCT2D( FMT_UINT8|FMT_REAL, TE_FFT_OPT_SCALE_METH|TE_DCT2D , TE_ALIGN_YES, dct2d_8x16, NULL ) },
  { 1, "dct1/idct2d_16x8.seq", (tTestEngTarget)idct2d_16x8, TEST_DESC_DCT2D( FMT_INT16|FMT_REAL, TE_FFT_OPT_SCALE_METH|TE_IDCT2D, TE_ALIGN_YES, NULL, idct2d_16x8) },

  { 0 } /* End of table */
};
 
/* Perform all tests for dct API functions. */
int func_dct1(int isFull, int isVerbose, int breakOnError)
{
  int ix, res;

  for ( ix=0,res=1; testTbl[ix].seqFilename && ( res || !breakOnError ); ix++ )
  {
    res &= ( 0 != TestEngRun( (tTestEngTarget )testTbl[ix].target,
                              &testTbl[ix].desc.desc,
                              testTbl[ix].seqFilename,
                              isFull, isVerbose, breakOnError,0 ) );
  }

  return (res);
}
