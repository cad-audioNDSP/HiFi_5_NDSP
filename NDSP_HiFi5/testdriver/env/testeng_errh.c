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
 * Test engine implementation.
 * Error Handling (ERRH) verification support functions.
 */

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <fenv.h>

/* Portable data types */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"
/* Utility functions. */
#include "utils.h"
/* Test engine API. */
#include "testeng.h"

/* Initialize the ERRH support; is called internally by the Test Engine before 
 * it reads a SEQ-file for the first time. */
void te_errh_init( tTestEngContext * context )
{
  int isEnabled = 0, extMask;
  char strBuf[16];

  ASSERT( context && context->seqFile );
  memset(strBuf,0,sizeof(strBuf));
  if ( 1 == seqFileScanf( context->seqFile, "ERRH_ENABLE: %15s", strBuf ) )
  {
    isEnabled = !strcmp( strBuf, "ON" ); 
  }

  extMask = TE_ERRH_EXTENDED_TEST_DISABLE | TE_ERRH_EXTENDED_TEST_ENABLE;
  if ( isEnabled && ( context->desc->extraParam & extMask ) == TE_ERRH_EXTENDED_TEST_ENABLE )
  {
    printf( "[Extended FP exceptions test] " );
    fflush( stdout );
  }

  memset( &context->errh, 0, sizeof(context->errh) );

  context->errh.isEnabled = isEnabled;
  context->errh.isPassed  = 0;
  context->errh.entranceExcepts = 0;
  context->errh.exceptEnable    = 0;
} /* te_errh_init() */

/* Read the SEQ-file to load reference error states for a test case.
 * Return zero if failed. */
int te_errh_loadRef( tTestEngContext * context )
{
  char sname[16];
  int snum, rc = 1;

  ASSERT( context && context->seqFile );

  memset( &context->errh.refStates, 0, sizeof( context->errh.refStates ) );

  while ( 2 == seqFileScanf( context->seqFile, "ERRH_%15s %d", sname, &snum ) )
  {
    tVec * svec = 0;
    rc = 0;

         if ( !strcmp( sname, "EDOM"    ) ) svec = &context->errh.refStates.edom;
    else if ( !strcmp( sname, "ERANGE"  ) ) svec = &context->errh.refStates.erange;
    else if ( !strcmp( sname, "FE_INV"  ) ) svec = &context->errh.refStates.fe_inv;
    else if ( !strcmp( sname, "FE_DIVZ" ) ) svec = &context->errh.refStates.fe_divz;
    else if ( !strcmp( sname, "FE_OVFL" ) ) svec = &context->errh.refStates.fe_ovfl;

    if ( !svec )
    {
      printf( "te_errh_loadRef(): unknown error state %s\n", sname );
    }
    else if ( svec->szBulk > 0 )
    {
      printf( "te_errh_loadRef(): multiple definitions of error state %s\n", sname );
    }
    else if ( !vecAlloc( svec, snum, TE_ALIGN_YES, FMT_INT32, 0 ) )
    {
      printf( "te_errh_loadRef(): failed to allocate vector of %d "
              "indices for error state %s\n", snum, sname );
    }
    else if ( !seqFileReadVec( context->seqFile, svec ) )
    {
      printf( "te_errh_loadRef(): failed to load %d indices "
              "for error state %s\n", snum, sname );
    }
    else rc = 1;

    if ( !rc ) break;
  }

  if ( !rc )
  {
    if ( context->errh.refStates.edom   .szBulk > 0 ) vecFree( &context->errh.refStates.edom    );
    if ( context->errh.refStates.erange .szBulk > 0 ) vecFree( &context->errh.refStates.erange  );
    if ( context->errh.refStates.fe_inv .szBulk > 0 ) vecFree( &context->errh.refStates.fe_inv  );
    if ( context->errh.refStates.fe_divz.szBulk > 0 ) vecFree( &context->errh.refStates.fe_divz );
    if ( context->errh.refStates.fe_ovfl.szBulk > 0 ) vecFree( &context->errh.refStates.fe_ovfl );

    memset( &context->errh.refStates, 0, sizeof( context->errh.refStates ) );
  }
  /* Reset the ERRH validation result flag for a test case. We have to reset it here,
   * because te_errh_resetStates()/te_errh_verifyStates() may be called multiple times,
   * (for a scalar function under test), while the validation result must cover the
   * whole test case. */
  else context->errh.isPassed = 1;

  return (rc);

} /* te_errh_loadRef() */

/* Free reference error state vectors. */
void te_errh_freeRef( tTestEngContext * context )
{
  ASSERT( context );

  if ( context->errh.refStates.edom   .szBulk > 0 ) vecFree( &context->errh.refStates.edom    );
  if ( context->errh.refStates.erange .szBulk > 0 ) vecFree( &context->errh.refStates.erange  );
  if ( context->errh.refStates.fe_inv .szBulk > 0 ) vecFree( &context->errh.refStates.fe_inv  );
  if ( context->errh.refStates.fe_divz.szBulk > 0 ) vecFree( &context->errh.refStates.fe_divz );
  if ( context->errh.refStates.fe_ovfl.szBulk > 0 ) vecFree( &context->errh.refStates.fe_ovfl );

  memset( &context->errh.refStates, 0, sizeof( context->errh.refStates ) );
} /* te_errh_freeRef() */

/* Prepare to call the function under test: reset the errno and clear
 * FP exception flags. */
void te_errh_resetStates( tTestEngContext * context )
{
  int mask, excepts;

  ASSERT( context );
 // (void)excepts;
  excepts = 0;
  mask = ( FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW | FE_INEXACT );
  (void)mask;

  if ( context->errh.isEnabled && context->errh.isPassed )
  {
    errno = 0;
    if (GET_ISA_OPT(HAVE_FP))
    {
      /* Reset all exception flags to zero. */
      if ( 0 != FECLEAREXCEPT(FE_ALL_EXCEPT) )
      {
        printf( "te_errh_resetStates(): FECLEAREXCEPT() failed\n" );
      }

      /* For an extended error handling test, randomize the exception flags. When the
       * FUT returns we will check that none of asserted flags is removed. */
      if ( (context->desc->extraParam & TE_ERRH_EXTENDED_TEST_ENABLE) !=0)
      {
        if ( 0 != FERAISEEXCEPT(excepts = ( Rand() & mask )) )
        {
          printf( "te_errh_resetStates(): FERAISEEXCEPT() failed\n" );
        }

        context->errh.entranceExcepts = excepts;
      }
      else
      {
        context->errh.entranceExcepts = 0;
      }
    }
  }

} /* te_errh_resetStates() */

/* Look if a reference state is asserted for the given index. For a vector
 * function under test idx must hold a negative value. */
static int getRefState( const tVec * svec, int idx )
{
  int n, sflag;

  ASSERT( svec );

  if ( idx >= 0 )
  {
    for ( n=sflag=0; n<svec->nElem && !sflag; n++ )
    {
      sflag = ( idx == *vecGetElem_i32( svec, n ) );
    }
  }
  else
  {
    sflag = ( svec->nElem > 0 );
  }

  return (sflag);

} /* getRefState() */

/* Compare actual error state with reference, and print a report on an inconsistency. */
static int checkState( int sactual, int sref, const char * sname, int idx, int isVerbose )
{
  int smatch = ( !sactual == !sref );

  ASSERT( sname );

  if ( !smatch && isVerbose )
  {
    printf( "Error Handling verification: %s %s for idx=%d\n",
            ( sref ? "missing" : "unexpected" ), sname, idx );
  }

  return (smatch);

} /* checkState() */

/* Sample errno and FP exception flags and check them against the reference
 * error states. For a scalar function under test the second argument must hold
 * the current index value for in/out test data vectors. For a vector function,
 * assign it a negative value. */
int te_errh_verifyStates( tTestEngContext * context, int idx )
{
  struct {
    int edom, erange, fe_inv, fe_divz, fe_ovfl, fe_unfl, fe_inex;
  } sref, sactual;

  int excepts, rc;

  ASSERT( context );
  (void)excepts;
  if ( !context->errh.isEnabled )
  {
    return ( context->errh.isPassed = 1 );
  }
  else if ( !context->errh.isPassed )
  {
    return (0);
  }

  /*
   * Figure out expected (reference) error states for index idx.
   */ 

  memset( &sref, 0, sizeof(sref) );

  sref.edom    = getRefState( &context->errh.refStates.edom   , idx );
  sref.erange  = getRefState( &context->errh.refStates.erange , idx );
  /* For vector functions, EDOM takes precedence over ERANGE. */
  if ( idx<0 ) sref.erange &= !sref.edom;

  sref.fe_inv  = getRefState( &context->errh.refStates.fe_inv , idx );
  sref.fe_divz = getRefState( &context->errh.refStates.fe_divz, idx );
  sref.fe_ovfl = getRefState( &context->errh.refStates.fe_ovfl, idx );

  /*
   * If any FP exception flags were non-zero before the FUT call, then
   * they must have been preserved.
   */

  if ( 0 != context->errh.entranceExcepts )
  {
    sref.fe_inv  |= ( 0 != ( context->errh.entranceExcepts & FE_INVALID   ) );
    sref.fe_divz |= ( 0 != ( context->errh.entranceExcepts & FE_DIVBYZERO ) );
    sref.fe_ovfl |= ( 0 != ( context->errh.entranceExcepts & FE_OVERFLOW  ) );
    sref.fe_unfl |= ( 0 != ( context->errh.entranceExcepts & FE_UNDERFLOW ) );
    sref.fe_inex |= ( 0 != ( context->errh.entranceExcepts & FE_INEXACT   ) );
  }

  /* FE_OVERFLOW and FE_UNDERFLOW entail the FE_INEXACT! */
  sref.fe_inex |= ( sref.fe_ovfl | sref.fe_unfl );

  /*
   * Determine actual error states.
   */

  memset( &sactual, 0, sizeof(sactual) );

  sactual.edom    = ( errno == EDOM   );
  sactual.erange  = ( errno == ERANGE );

  if (GET_ISA_OPT(HAVE_FP))
  {
    excepts = FETESTEXCEPT( FE_ALL_EXCEPT );

    sactual.fe_inv  = ( 0 != ( excepts & FE_INVALID   ) );
    sactual.fe_divz = ( 0 != ( excepts & FE_DIVBYZERO ) );
    sactual.fe_ovfl = ( 0 != ( excepts & FE_OVERFLOW  ) );
    sactual.fe_unfl = ( 0 != ( excepts & FE_UNDERFLOW ) );
    sactual.fe_inex = ( 0 != ( excepts & FE_INEXACT   ) );

    /*
     * Compare actual error states with reference.
     */

    rc  = ( 0 != checkState( sactual.edom   , sref.edom   , "EDOM"   , idx, context->isVerbose ) );
    rc &= ( 0 != checkState( sactual.erange , sref.erange , "ERANGE" , idx, context->isVerbose ) );
    rc &= ( 0 != checkState( sactual.fe_inv , sref.fe_inv , "FE_INVALID"  , idx, context->isVerbose ) );
    rc &= ( 0 != checkState( sactual.fe_divz, sref.fe_divz, "FE_DIVBYZERO", idx, context->isVerbose ) );
    rc &= ( 0 != checkState( sactual.fe_ovfl, sref.fe_ovfl, "FE_OVERFLOW" , idx, context->isVerbose ) );

    if ( 0 == ( context->desc->extraParam & TE_ERRH_IGNORE_FE_UNDERFLOW ) )
    {
      rc &= ( 0 != checkState( sactual.fe_unfl, sref.fe_unfl, "FE_UNDERFLOW", idx, context->isVerbose ) );
    }
    if ( 0 == ( context->desc->extraParam & TE_ERRH_IGNORE_FE_INEXACT ) )
    {
      rc &= ( 0 != checkState( sactual.fe_inex, sref.fe_inex, "FE_INEXACT"  , idx, context->isVerbose ) );
    }
  }
  else
  {
    rc = 1;
  }

  return ( context->errh.isPassed &= rc );

} /* te_errh_verifyStates() */
