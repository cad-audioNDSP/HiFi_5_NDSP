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
 */

#include <stdio.h>
#include <string.h>

/* Portable data types */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"
/* Test engine API. */
#include "testeng.h"
/* Verification report. */
#include "vreport.h"
/*  Test case type textual representation */
#include "testcase.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

/* Validate the test result. In the case of a failure optionally prints
 * a detailed report and returns zero. */
static int validateTestResult( tTestEngContext * context, int isVerbose ); 
static int validateBitexactness(tTestEngContext * context, const tVec * Zvec, const tVec * Wvec, int isVerbose );
/* Free vectors allocated for the test case. Return zero if guard memory areas check fails. */
int te_freeVectors( tTestEngContext * context );

/* Test executive function. Performs the specified test on a brief or full version
 * of the designated SEQ-file. Return the test result (non-zero if passed). */
int TestEngRun( tTestEngTarget targetFxn, const tTestEngDesc * desc,
                const char * seqName, int isFull, int isVerbose, int breakOnError, int testBitexactness )
{
  tTestEngContext context;
  tSeqFile sf;

  const char * seqDir;
  char fnameBuf[256];
  int n = 0, res = 1, testResult = 1, tcaseResult=0;

  const int dimNum = (int)sizeof(context.args.dim)/sizeof(context.args.dim[0]);
  const int caseTypeLim = (int)sizeof(testCaseTypeStr)/sizeof(testCaseTypeStr[0]) - 1;

  ASSERT( desc && seqName );

  /*
   * Select either full or brief variant of the SEQ-file.
   */

  /* Select either full or brief variant of the SEQ-file. */
  seqDir = getVectorsDir(isFull);
  ASSERT( strlen(seqDir)+strlen(seqName) < sizeof(fnameBuf) );
  sprintf( fnameBuf, "%s/%s", seqDir, seqName );
  
  printf( "%s... ", fnameBuf );
  fflush( stdout );

  /*
   * Check if the target function is included in the current build.
   */

  if ( !IS_PRESENT(targetFxn) )
  {
    printf( "NOT TESTED\n" );
    return (1);
  }

  /*
   * Open the SEQ-file.
   */

  res = ( 0 != ( sf = seqFileOpen( fnameBuf ) ) );
  if (testBitexactness)
  {
      printf( "[bitexactness vs scalar function] " );
      fflush( stdout );
  }
  if ( isVerbose ) printf( "\n" );

  /*
   * Prepare the test case context structure.
   */

  memset( &context, 0, sizeof(context) );
  context.desc       = desc;
  context.seqFile    = sf;
  context.isFull     = isFull;
  context.isVerbose  = isVerbose;
  context.doRetake   = 0;
  context.target.fut = targetFxn;

  ASSERT( strlen(seqDir) < sizeof(context.seqDir) );
  strcpy( context.seqDir, seqDir );

  if ( res )
  {
    ASSERT( TE_DIM_NUM_1 <= desc->dimNum && desc->dimNum <= TE_DIM_NUM_5 );
    ASSERT( desc->dimNum <= dimNum );

    for ( n=0; n<dimNum; n++ )
    {
      context.args.dim[n] = 1;
    }

    /* Initialize the Error Handling verification support. */
    te_errh_init( &context );
  }

  /*
   * Create the target algorithm instance (optional).
   */

  if ( desc->targetCreateFxn && res )
  {
    res = desc->targetCreateFxn( &context );
  }

  /*
   * Perform test cases until:
   *   A) SEQ-file is exhausted or
   *   B) test engine fails to perform a test case (e.g. because of a memory allocation error)
   *   C) some test fails, and breakOnError is true or isVerbose is false (there is no sense
   *      to continue testing of once failed function whenever isVerbose is false).
   */

  if ( !isVerbose ) breakOnError = 1;

  while ( sf && res>0 && ( testResult || !breakOnError ) )
  {
    /* Read the test case number and semantic type. */
    if ( 2 != seqFileScanf( sf, "%d %d",
                            &context.args.caseNum,
                            &context.args.caseType ) ) 
    {
      /* If the SEQ-file is not exhausted, but we cannot read the test case
       * number and type then we clearly face a format error. */
      if ( !( res = feof(sf->file) ) ) printf( "bad SEQ-file format (1)\n" );
      break;
    }

    if ( isVerbose )
    {
      n = MIN( context.args.caseType, caseTypeLim );
      printf( "#%-3d <%s> ", context.args.caseNum, testCaseTypeStr[n] );
    }

    /* Read dimensional parameters. */
    for ( n=0; n<desc->dimNum; n++ )
    {
      if ( 1 != seqFileScanf( sf, "%d", &context.args.dim[n] ) ) break;
    }

    /* Verify that dimensional parameters have been read successfully. */
    if ( !( res = ( n == desc->dimNum ) ) )
    {
      printf( "bad SEQ-file format (2)\n" );
      break;
    }

    /* Keep iterating the load/process/validate sequence until the test case 
     * is done or it fails. */
    do
    {
      context.doRetake = 0;
      /* Allocate data vectors and load test set data from the SEQ-file. */
      res = desc->testCaseLoadFxn( &context );
      if (res<0) tcaseResult=1;
      if (res<=0) break;

      /* Load reference error states, if any. */
      if ( !( res = te_errh_loadRef( &context ) ) ) break;

      /* Apply the target function. */
      desc->testCaseProcessFxn( &context );

      /* Validate the test case result. */
      tcaseResult = ( 0 != validateTestResult( &context, isVerbose ) );
      if (testBitexactness && IS_PRESENT(context.target.fut))
      {
          const tTestEngVecSclTbl *tbl;
          te_fun_ptr_t vec=context.target.fut;
          /* clone Z and W arrays */
          tVec cloneZ,cloneW;
          if (0==vecClone(&cloneZ, &context.dataSet.Z) || 0==vecClone(&cloneW, &context.dataSet.W))
          {
              printf( "Failed to clone Z,W vectors\n" );
              tcaseResult = 0;
          }
          else
          {
              /* find scalar test case processing function from the map */
              for (tbl=(const tTestEngVecSclTbl *)desc->extraPtr; tbl->testCaseProcessFxn; tbl++)
              {
                  if (tbl->vec==vec) break;
              }
              if (tbl->vec==vec) 
              {
                  te_fun_ptr_t oldfxn= context.target.fut;
                  context.target.fut=tbl->scl;
                  tbl->testCaseProcessFxn( &context );
                  context.target.fut=oldfxn;
                  tcaseResult &= ( 0 != validateTestResult( &context, isVerbose ) );
                  tcaseResult &= ( 0 != validateBitexactness( &context, &cloneZ, &cloneW, isVerbose));
              }
              else
              {
                  ASSERT(0);  /* mapping scalar->vector is not found */
              }
          }
          /* free arrays */
          if (cloneZ.szBulk>0) vecFree(&cloneZ);
          if (cloneW.szBulk>0) vecFree(&cloneW);
      }
      vReportSetResult(tcaseResult ? evr_OK : evr_FAIL, evr_NOTTESTED,evr_NOTTESTED);
      /* Attach the result of Error Handling validation, if applicable. */
      if ( context.errh.isEnabled ) 
      {
          tcaseResult &= ( 0 != context.errh.isPassed );
          vReportSetResult(evr_NOTTESTED, context.errh.isPassed ? evr_OK : evr_FAIL, evr_NOTTESTED);
      }

    } while ( context.doRetake && tcaseResult );

    /* Update the overall test result. */
    testResult &= tcaseResult;

    /* Free reference error states. */
    te_errh_freeRef( &context );

    /* Free test vectors. */
    if ( !te_freeVectors( &context ) ) res = 0;

    /* Print the test case result. */
    if ( isVerbose )
    {
      if (res == -1) printf("NOT TESTED\n");
      else  { //      printf((res && testResult) ? "OK\n" : "FAIL\n");
      if ( res && tcaseResult ) printf( "OK\n"   );
      else if ( !breakOnError ) printf( "FAIL\n" );
      }
    }
  }

  /* Destroy the algorithm instance. */
  if ( desc->targetDestroyFxn )
  {
    desc->targetDestroyFxn( &context );
  }

  /* Close the SEQ-file. */
  if ( sf ) seqFileClose( sf );

  /* Print the overall test result. */
  if (res==-1) printf( "NOT TESTED\n" );
  else         printf( ( res && testResult ) ? "OK\n" : "FAIL\n" );
  fflush(stdout);
  vReportFinish();
  return ( res && testResult );

} /* TestEngRun() */

/* Validate the test result. In the case of a failure optionally prints
 * a detailed report and returns zero. */
static int validateTestResult( tTestEngContext * context, int isVerbose )
{
  int failIx, ix, resZ = 1, resW = 1;
  char buf[128];
  int n;

  const tVec *X, *Y, *U, *V, *Z, *W;
  const tVec *Zlo, *Zhi, *Wlo, *Whi;

  ASSERT( context );

  X   = &context->dataSet.X;
  Y   = &context->dataSet.Y;
  U   = &context->dataSet.U;
  V   = &context->dataSet.V;
  Z   = &context->dataSet.Z;
  W   = &context->dataSet.W;
  Zlo = &context->dataSet.Zlo;
  Zhi = &context->dataSet.Zhi;
  Wlo = &context->dataSet.Wlo;
  Whi = &context->dataSet.Whi;

  if ( Z->nElem > 0 ) resZ = ( 0 != vecCheckRange( Z, Zlo, Zhi, &failIx ) );
  if ( W->nElem > 0 ) resW = ( 0 != vecCheckRange( W, Wlo, Whi, &failIx ) );

  if ( isVerbose && ( !resZ || !resW ) )
  {
    const char * legend[] = { "s Z,W", " Z", " W", "" };

    printf( "\n##### start of error report #####\n" );

    printf( "Dimensions: " );

    for ( n=0; n<context->desc->dimNum; n++ )
    {
      printf( "%d ", context->args.dim[n] );
    }

    printf( "\n" );
    if ((context->desc->extraParam & OVLP)!=0)
    {
        printf( "Test with overlapped input/output arguments (OVLP_xxx)\n" );
    }

    printf( "Test result validation failed for vector%s at index %d\n", legend[(resZ<<1)|resW], failIx );
    printf( "Input: " );

    if ( X->nElem > 0 )    
    {
      ix = MIN( X->nElem-1, failIx );
      vecPrintElem( buf, X, ix );
      printf( "x[%d] = %s; ", ix, buf );
    }

    if ( Y->nElem > 0 )    
    {
      ix = MIN( Y->nElem-1, failIx );
      vecPrintElem( buf, Y, ix );
      printf( "y[%d] = %s; ", ix, buf );
    }

    if ( U->nElem > 0 )    
    {
      ix = MIN( U->nElem-1, failIx );
      vecPrintElem( buf, U, ix );
      printf( "u[%d] = %s; ", ix, buf );
    }

    if ( V->nElem > 0 )    
    {
      ix = MIN( V->nElem-1, failIx );
      vecPrintElem( buf, V, ix );
      printf( "v[%d] = %s; ", ix, buf );
    }

    if ( Z->nElem > 0 && Zlo->nElem > 0 && Zhi->nElem > 0 )
    {
      printf( "\nResult range: " );

      ix = MIN( Zlo->nElem-1, failIx );
      vecPrintElem( buf, Zlo, ix );
      printf( "zlo[%d] = %s; ", ix, buf );

      ix = MIN( Zhi->nElem-1, failIx );
      vecPrintElem( buf, Zhi, ix );
      printf( "zhi[%d] = %s; ", ix, buf );

      printf( "\nActual result: " );

      ix = MIN( Z->nElem-1, failIx );
      vecPrintElem( buf, Z, ix );
      printf( "z[%d] = %s; ", ix, buf );
    }

    if ( W->nElem > 0 && Wlo->nElem > 0 && Whi->nElem > 0 )
    {
      printf( "\nResult range: " );

      ix = MIN( Wlo->nElem-1, failIx );
      vecPrintElem( buf, Wlo, ix );
      printf( "wlo[%d] = %s; ", ix, buf );

      ix = MIN( Whi->nElem-1, failIx );
      vecPrintElem( buf, Whi, ix );
      printf( "whi[%d] = %s; ", ix, buf );

      printf( "\nActual result: " );

      ix = MIN( W->nElem-1, failIx );
      vecPrintElem( buf, W, ix );
      printf( "w[%d] = %s; ", ix, buf );
    }

    printf( "\n##### end of error report #####\n" );
  }

  return ( resZ && resW );

} /* validateTestResult() */

/* Validate bitexactness. In the case of a failure optionally prints
 * a detailed report and returns zero. */
static int validateBitexactness(tTestEngContext * context, const tVec * Zvec, const tVec * Wvec, int isVerbose )
{
  int failIx, ix, resZ = 1, resW = 1;
  const tVec *Z, *W;
  char buf[128];
  Z   = &context->dataSet.Z;
  W   = &context->dataSet.W;

  if ( Zvec->nElem > 0 ) resZ = ( 0 != vecCheckRange( Z, Zvec, Zvec, &failIx));
  if ( Wvec->nElem > 0 ) resW = ( 0 != vecCheckRange( W, Wvec, Wvec, &failIx));

  if ( isVerbose && ( !resZ || !resW ) )
  {
    const char * legend[] = { "s Z,W", " Z", " W", "" };
    printf( "\n##### start of error report #####\n" );
    printf( "Bitexacness test result validation failed for vector%s at index %d\n", legend[(resZ<<1)|resW], failIx );

    if ( Z->nElem > 0 )
    {
      printf( "\nVector function: " );
      ix = MIN( Z->nElem-1, failIx );
      vecPrintElem( buf, Zvec, ix );
      printf( "z[%d] = %s; ", ix, buf );
      printf( "\nScalar function: " );
      ix = MIN( Z->nElem-1, failIx );
      vecPrintElem( buf, Z, ix );
      printf( "z[%d] = %s; ", ix, buf );
    }

    if ( W->nElem > 0 )
    {
      printf( "\nVector function: " );
      ix = MIN( W->nElem-1, failIx );
      vecPrintElem( buf, Wvec, ix );
      printf( "w[%d] = %s; ", ix, buf );
      printf( "\nScalar function: " );
      ix = MIN( W->nElem-1, failIx );
      vecPrintElem( buf, W, ix );
      printf( "w[%d] = %s; ", ix, buf );
    }

    printf( "\n##### end of error report #####\n" );
  }

  return ( resZ && resW );

} /* validateBitexactness() */

/* Free vectors allocated for the test case. Return zero if guard memory areas check fails. */
int te_freeVectors( tTestEngContext * context )
{
  const int auxNum = (int)sizeof(context->dataSet.auxVec)/sizeof(context->dataSet.auxVec[0]);
  int n, res;

  ASSERT( context );

  res = ( ( !context->dataSet.X  .szBulk || 0 != vecFree( &context->dataSet.X   ) ) &
          ( !context->dataSet.Y  .szBulk || 0 != vecFree( &context->dataSet.Y   ) ) &
          ( !context->dataSet.U  .szBulk || 0 != vecFree( &context->dataSet.U   ) ) &
          ( !context->dataSet.V  .szBulk || 0 != vecFree( &context->dataSet.V   ) ) &
          ( !context->dataSet.Z  .szBulk || 0 != vecFree( &context->dataSet.Z   ) ) &
          ( !context->dataSet.Zlo.szBulk || 0 != vecFree( &context->dataSet.Zlo ) ) &
          ( !context->dataSet.Zhi.szBulk || 0 != vecFree( &context->dataSet.Zhi ) ) &
          ( !context->dataSet.W  .szBulk || 0 != vecFree( &context->dataSet.W   ) ) &
          ( !context->dataSet.Wlo.szBulk || 0 != vecFree( &context->dataSet.Wlo ) ) &
          ( !context->dataSet.Whi.szBulk || 0 != vecFree( &context->dataSet.Whi ) ) );

  for ( n=0; n<auxNum; n++ )
  {
    if ( context->dataSet.auxVec[n].szBulk )
    {
      res &= ( 0 != vecFree( &context->dataSet.auxVec[n] ) );
    }
  }
  vReportSetResult(evr_NOTTESTED,evr_NOTTESTED,res ? evr_OK:evr_FAIL);
  return (res);

} /* te_freeVectors() */

/* Test executive function. Look through the test definition table, find an approptiate test 
 * description, and perform the specified test on a brief or full version of the designated
 * SEQ-file. Return the test result (non-zero if passed). 
 *
 * Test definition table (first argument) ties together test target functions and test
 * descriptions. It must contain entries that match the following structure:
 *   struct 
 *   {
 *     tTestEngTarget * funcList[maxFuncNum];
 *     tTestEngDesc     testDesc;
 *   };
 * 
 * Test definition table size and function list size must be specified through tblSize
 * and maxFuncNum arguments, respectively. Unused elements of funcList must be NULL.
 */
int te_Exec( const void * pTestDefTbl, size_t tblSize, size_t maxFuncNum,
             tTestEngTarget  targetFxn, const char * seqName,
             int isFull, int isVerbose, int breakOnError )
{
  tTestEngTarget *  pFxns;
  const tTestEngDesc * pDesc;
  int tblIx, funcIx;
  for ( tblIx=0; tblIx<(int)tblSize; tblIx++ )
  {
    pFxns = (tTestEngTarget*)pTestDefTbl;
    pDesc = (tTestEngDesc*)( pFxns + maxFuncNum );

    for ( funcIx=0; funcIx<(int)maxFuncNum; funcIx++ )
    {
      if ( targetFxn == pFxns[funcIx] )
      {
        return ( TestEngRun( targetFxn, pDesc, seqName, 
                             isFull, isVerbose, breakOnError, 0 ) );
      }
    }

    pTestDefTbl = pDesc + 1;
  }
  
  ASSERT( !"Test not defined" );
  return (0);

} /* te_Exec() */

/* function returns total data size allocated by input/output data for given test case */
size_t te_vGetDataSize(const tTestEngContext * context)
{
    int n;
    size_t sz=0;
    if (context->dataSet.X.szBulk) sz+=context->dataSet.X.nElem*context->dataSet.X.szElem;
    if (context->dataSet.Y.szBulk) sz+=context->dataSet.Y.nElem*context->dataSet.Y.szElem;
    if (context->dataSet.U.szBulk) sz+=context->dataSet.U.nElem*context->dataSet.U.szElem;
    if (context->dataSet.V.szBulk) sz+=context->dataSet.V.nElem*context->dataSet.V.szElem;
    if (context->dataSet.Z.szBulk) sz+=context->dataSet.Z.nElem*context->dataSet.Z.szElem;
    if (context->dataSet.W.szBulk) sz+=context->dataSet.W.nElem*context->dataSet.W.szElem;
    for (n=0; n<(int)(sizeof(context->dataSet.auxVec)/sizeof(context->dataSet.auxVec[0])); n++) 
    {
        if (context->dataSet.auxVec[n].szBulk) sz+=context->dataSet.auxVec[n].nElem*context->dataSet.auxVec[n].szElem;
    }
    return sz;
}

/* standard verfication reporting function 
   should be called if FUT is written to context->target.fut
   and there is no specific variants of calling function 
*/
void te_vReportStd(tTestEngContext * context)
{
    vReportAdd((tReportFUT*)&context->target.fut,1,NULL,context->seqFile->filename,context->args.caseType,te_vGetDataSize(context));
}
/*--------------------------------------------------------------
    function returns vectors directory
--------------------------------------------------------------*/
const char* getVectorsDir(int isFull)
{
    switch(isFull)
    {
    case 0:  return BRIEF_VECTOR_DIR PACKAGE_SUFFIX;
	  case 1:  return FULL_VECTOR_DIR PACKAGE_SUFFIX;
	  case 2:  return SANITY_VECTOR_DIR;
	  default: return SANITY_VECTOR_DIR;
    }
}
