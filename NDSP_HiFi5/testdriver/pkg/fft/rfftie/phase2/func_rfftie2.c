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
* Test procedures for complex FFT functions.
*/
#include "../common/test_rfftie.h"

/* Test data file processing for a standalone test. Applies the target FFT function
 * to test data, compares the FFT output with reference data and estimates the SINAD. */
static void fileProc_stnd_rfftie( tTestEngContext * context )
{
  tTestEngContext_fft           * context_fft;
  tTestEngFrameProcFxn_fft_stnd * frameProcFxn;

  tVec inVec, outVec, refVec; /* In/out/reference vectors related to test data files. */
  tVec xVec, yVec;            /* In/out vectors for the target FFT routine. */
  FILE *fIn = 0, *fRef = 0;
  int res = 0;
  int lenIn=0, lenOut=0;
  int lenFread = 0;  /* Number of samples readed from file */
  int fmtX=0, fmtY=0;
  int n, N;

  NASSERT( context && context->desc );
  NASSERT( context->target.fut && context->target.handle );

  context_fft = (tTestEngContext_fft*)context->target.handle;

  /* FFT size */
  N = context->args.N_DIM_;

  /* Select length and format (real/complex) for all vectors. */
  switch ( context->args.caseType )
  {
  case TE_FFT_TESTCASE_FORWARD:
    lenIn  = N;
    lenFread = N; 
    lenOut = ( ( context->desc->extraParam & TE_FFT_FULL_SPECTRUM ) ? N : N/2+1 );
    fmtX   = FMT_REAL;
    fmtY   = FMT_CPLX;
    break;
  case TE_FFT_TESTCASE_INVERSE:
    lenIn  = N/2+1;
    lenFread = lenIn - 1; 
    lenOut = N;
    fmtX   = FMT_CPLX;
    fmtY   = FMT_REAL;
    break;
  default:
    ASSERT( 0 );
  }

  /* Preset an invalid SINAD for the test result. */
  *vecGetElem_fl32( &context->dataSet.Z, 0 ) = -HUGE_VALF;

  if ( context->isVerbose )
  {
    if ( context->desc->extraParam & TE_FFT_OPT_SCALE_METH )
    {
      printf( "scale_mtd %d ", context_fft->scale_method );
    }

    printf( "%-40s  ", context_fft->fRefName );
  }

  memset( &inVec , 0, sizeof(inVec ) );
  memset( &outVec, 0, sizeof(outVec) );
  memset( &refVec, 0, sizeof(refVec) );
  memset( &xVec, 0, sizeof(xVec) );
  memset( &yVec, 0, sizeof(yVec) );
    {
        char * fname;
        const char* dir;
		dir = getVectorsDir(context->isFull);
        fname=(char*)malloc(strlen(context_fft->fInName)+strlen(dir)+2+6);
        sprintf(fname,"%s/%s",dir,context_fft->fInName);
        fIn = fopen(fname, "rb");
        free(fname);
        fname=(char*)malloc(strlen(context_fft->fRefName)+strlen(dir)+2+6);
        sprintf(fname,"%s/%s",dir,context_fft->fRefName);
        fRef = fopen(fname, "rb");
        free(fname);
    }

  /* Allocate vectors for in/out/reference data. */
  if ( !vecAlloc( &inVec , lenIn , TE_ALIGN_YES, FMT_FRACT16|fmtX, 0 ) ||
       !vecAlloc( &outVec, lenOut, TE_ALIGN_YES, FMT_FLOAT64|fmtY, 0 ) ||
       !vecAlloc( &refVec, N     , TE_ALIGN_YES, FMT_FLOAT64|fmtY, 0 ) )
  {
    printf( "fileProc_stnd_rfftie(): failed to allocate in/out/ref vectors, N=%d\n", N );
  }
  /* Open input data file. */
  else if (fIn==NULL )
  {
    printf( "fileProc_stnd_rfftie(): failed to open %s for reading\n", context_fft->fInName );
  }
  /* Open reference data file. */
  else if (fRef == NULL )
  {
    printf( "fileProc_stnd_rfftie(): failed to open %s for reading\n", context_fft->fRefName );
  }
  else res = 1;

  /*
   * Process the input file frame-by-frame.
   */

  if ( res )
  {
    int baseFmt = ( context->desc->fmt & FMT_DTYPE_MASK );
    int isAligned = context->desc->isAligned;

    int16_t   * pin  = (int16_t   *)vecGetElem( &inVec , 0 );
    float64_t * pout = (float64_t *)vecGetElem( &outVec, 0 );
    float64_t * pref = (float64_t *)vecGetElem( &refVec, 0 );

    float32_t sinadAvg, sinadMin = HUGE_VALF;
    float64_t errSum = 0, refSum = 0;

    int efbMin = 32;

    memset( &xVec, 0, sizeof(xVec) );
    memset( &yVec, 0, sizeof(yVec) );

    context_fft->frameCnt = 0;

    frameProcFxn = ( (tTestEngFrameProcFxn_fft_stnd*)context->target.fut );

    /* Read 16-bit complex samples from the input file. */
    while ((n = fread(pin, inVec.szElem, lenFread, fIn)) > 0)
    {
      /* Augment the (last) incomplete frame with zeros. */
        memset((uint8_t*)pin + n*inVec.szElem, 0, (lenFread - n)*inVec.szElem);
      /* Zero the output frame. */
      memset( pout, 0, lenOut*outVec.szElem );

      /* For inverse FFT test, unpack the Nyquist frequency bin. */
      if ( context->args.caseType == TE_FFT_TESTCASE_INVERSE )
      {

        pin[2*(lenIn-1)+0] = pin[1];
        pin[2*(lenIn-1)+1] = pin[1] = 0;
      }

      /* Allocate in/out buffers for the target FFT routine. */
      if ( ( !xVec.szBulk && !vecAlloc( &xVec, lenIn , isAligned, baseFmt | fmtX, 0 ) ) ||
           ( !yVec.szBulk && !vecAlloc( &yVec, lenOut, isAligned, baseFmt | fmtY, 0 ) ) )
      {
        printf( "fileProc_stnd_rfftie(): failed to allocate xVec/yVec, "
                "frameCnt=%d, N=%d\n", context_fft->frameCnt, N );
        res = 0; break;
      }
	  context_fft->dataSize = (size_t)lenIn+(size_t)lenOut;

      /* Use a proprietary frame processing function to perform the target FFT. */
      if ( !( res = frameProcFxn( context, pin, pout, &xVec, &yVec ) ) ) break;

      /* When in unaligned mode, in/out vectors should be periodically reallocated to try various
       * address offsets. */
      if ( !( ++context_fft->frameCnt % XY_VEC_REALLOC_PERIOD ) && !isAligned )
      {
        if ( !( res = vecsFree( &xVec, &yVec, 0 ) ) )
        {
          printf( "fileProc_stnd_rfftie(): vecsFree() failed for xVec/yVec, "
                  "frameCnt=%d, N=%d\n", context_fft->frameCnt, N );
          break;
        }

        memset( &xVec, 0, sizeof(xVec) );
        memset( &yVec, 0, sizeof(yVec) );
      }

      /* Read double precision reference samples. Reference frame always contains
       * exactly N samples (real or complex)! */
      if ( (int)fread( pref, refVec.szElem, N, fRef ) < N )
      {
        printf( "fileProc_stnd_rfftie(): failed to read reference data, "
                "frameCnt=%d, N=%d\n", context_fft->frameCnt, N );
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
      sinadMin = 10.f*log10f( sinadMin );

      /* Set the test result for test verification. */
      *vecGetElem_fl32( &context->dataSet.Z, 0 ) = sinadMin;

      if ( context->isVerbose )
      {
        sinadAvg = ( refSum > 0 ? 10.f*log10f( (float32_t)(refSum/errSum) ) : -HUGE_VALF );

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

} /* fileProc_stnd_rfftie() */

/* Initializer for a test description structure. */
#define TEST_DESC_STND( fmt, extraParam, align,fwdfft,invfft )   {{ (fmt),(extraParam),NULL, TE_ARGNUM_1,(align), \
                                                     &te_create_fft, &te_destroy_fft, \
                                                     &te_load_fft, &fileProc_stnd_rfftie },(tFrwTransFxn)fwdfft,(tInvTransFxn)invfft}

static TestDef_t testTbl_rfftie2[] =
{
#if 1
    { 2, "rfftie2/fft8_realf_ie.seq"   , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "rfftie2/fft16_realf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "rfftie2/fft32_realf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "rfftie2/fft64_realf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "rfftie2/fft128_realf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "rfftie2/fft256_realf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "rfftie2/fft512_realf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "rfftie2/fft1024_realf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "rfftie2/fft2048_realf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
    { 2, "rfftie2/fft4096_realf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, fft_realf_ie, NULL) },
              
    { 2, "rfftie2/ifft8_realf_ie.seq"   , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "rfftie2/ifft16_realf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "rfftie2/ifft32_realf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "rfftie2/ifft64_realf_ie.seq"  , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "rfftie2/ifft128_realf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "rfftie2/ifft256_realf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "rfftie2/ifft512_realf_ie.seq" , (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "rfftie2/ifft1024_realf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "rfftie2/ifft2048_realf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
    { 2, "rfftie2/ifft4096_realf_ie.seq", (tTestEngTarget)te_frameProc_stnd_fft, TEST_DESC_STND(FMT_FLOAT32, ATTR_FAST, TE_ALIGN_YES, NULL, &ifft_realf_ie) },
#endif
    { 0 } /* End of table */
}; // testTbl_rfftie

int func_rfftie2(int isFull, int isVerbose, int breakOnError)
{
    return main_rfftie(2, isFull, isVerbose, breakOnError, testTbl_rfftie2);
}
