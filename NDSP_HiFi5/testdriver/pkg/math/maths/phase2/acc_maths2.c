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
 * accuracy testing for scalar math
 */
#include <string.h>
#include <errno.h>
#include <fenv.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* Library API */
#include LIBRARY_HEADER(math)
/* Test Engine API. */
#include "testeng.h"
#include "fpstat.h" // test statistics
#include "testcase.h" // test case types
/* Test Engine add-on for Vector Mathematics functions */
#include "../../../vector/common/testeng_vector.h"
#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )
#include "float16.h"

static int  createFn(tTestEngContext* context);
static int  destroyFn(tTestEngContext* context);
static int  te_loadFxn_vXvZacc( tTestEngContext * context );
static int  te_loadFxn_vXvYvZacc( tTestEngContext * context );
static void te_math_processFxn_s_vXvZ_acc(tTestEngContext * context);
static void te_math_processFxn_s_vXvYvZ_acc(tTestEngContext * context);

static void te_math_processFxn_s_atan2_acc(tTestEngContext * context);

static void te_math_processFxn_s_vXvZ_exc(tTestEngContext * context);
static void te_math_processFxn_s_vXvYvZ_exc(tTestEngContext * context);

static void te_math_processFxn_s_atan2_exc(tTestEngContext * context);

/* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }
/* Initializer for a test description structure. */
#define MAX_FUNC_NUM   16
  /* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }

/* vec API test definitions. */
typedef struct tagTestDef
{
tTestEngTarget   funcList[MAX_FUNC_NUM];
tTestEngDesc     testDesc;
}
tTestDef;

#define TEST_DESC( fmt, argNum, align, loadFxn, procFxn ) { (fmt),0,NULL,(argNum),(align),createFn,destroyFn,(loadFxn),(procFxn) }
static const tTestDef testDefTblAccuracy[] =
{
  /*  Stage 2  */
  {
    FUNC_LIST(    (tTestEngTarget)&scl_log2f  ,
                  (tTestEngTarget)&scl_lognf  ,
                  (tTestEngTarget)&scl_log10f ,
                  (tTestEngTarget)&scl_antilog2f,
                  (tTestEngTarget)&scl_antilognf,
                  (tTestEngTarget)&scl_antilog10f),
    TEST_DESC(FMT_REAL | FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZacc, &te_math_processFxn_s_vXvZ_acc) },
  {
    FUNC_LIST(    (tTestEngTarget)&scl_sinef,
                  (tTestEngTarget)&scl_cosinef ,
                  (tTestEngTarget)&scl_tanf ,
                  (tTestEngTarget)&scl_atanf,
                  (tTestEngTarget)&scl_tanhf,
                  (tTestEngTarget)&scl_sigmoidf),
    TEST_DESC(FMT_REAL | FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZacc, &te_math_processFxn_s_vXvZ_acc) },
  {
    FUNC_LIST(    (tTestEngTarget)&scl_powf,
                  (tTestEngTarget)&scl_reluf),
    TEST_DESC(FMT_REAL | FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZacc, &te_math_processFxn_s_vXvYvZ_acc) },
  {
    FUNC_LIST((tTestEngTarget)&scl_atan2f),
    TEST_DESC(FMT_REAL | FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZacc, &te_math_processFxn_s_atan2_acc) },
  {
    FUNC_LIST(NULL), TEST_DESC(0, 0, 0, NULL, NULL) } /* End of table */
};
#undef TEST_DESC

#define TEST_DESC( fmt, argNum, align, loadFxn, procFxn ) { (fmt),0,NULL,(argNum),(align),NULL,NULL,(loadFxn),(procFxn) }
static const tTestDef testDefTblExc[] =
{
  /*  Stage 2  */
  {
    FUNC_LIST(    (tTestEngTarget)&scl_log2f  ,
                  (tTestEngTarget)&scl_lognf  ,
                  (tTestEngTarget)&scl_log10f ,
                  (tTestEngTarget)&scl_antilog2f,
                  (tTestEngTarget)&scl_antilognf,
                  (tTestEngTarget)&scl_antilog10f),
    TEST_DESC(FMT_REAL | FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_math_processFxn_s_vXvZ_exc) },
  {
    FUNC_LIST(    (tTestEngTarget)&scl_sinef,
                  (tTestEngTarget)&scl_cosinef ,
                  (tTestEngTarget)&scl_tanf ,
                  (tTestEngTarget)&scl_atanf,
                  (tTestEngTarget)&scl_tanhf,
                  (tTestEngTarget)&scl_sigmoidf),
    TEST_DESC(FMT_REAL | FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_math_processFxn_s_vXvZ_exc) },
  {
    FUNC_LIST(    (tTestEngTarget)&scl_powf,
                  (tTestEngTarget)&scl_reluf),
    TEST_DESC(FMT_REAL | FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_math_processFxn_s_vXvYvZ_exc) },
  {
    FUNC_LIST((tTestEngTarget)&scl_atan2f),
    TEST_DESC(FMT_REAL | FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_math_processFxn_s_atan2_exc) },
  {
		FUNC_LIST(NULL), TEST_DESC(0, 0, 0, NULL, NULL) } /* End of table */
};
#undef TEST_DESC

/* Test executive function. Performs the specified test on a brief or full version
* of the designated SEQ-file. Return the test result (non-zero if passed). */
static int testExec(const tTestDef* pTestDef, size_t szTestDef,
  tTestEngTarget targetFxn, const char * seqName,
  int errhExtendedTest, int isFull, int isVerbose, int breakOnError);

// prinf SP value
static void printfSP(float32_t* x)
{
  uint32_t* u=(uint32_t*)x;
  while(1)
  {
      if (u[0]==0)          { printf("positive zero"); break;}
      if (u[0]==0x80000000) { printf("negative zero"); break;}
      if (x[0]==1.f)        { printf("positive one"); break;}
      if (x[0]==-1.f)       { printf("negative one"); break;}
      if (u[0]==0x7f800000) { printf("positive infinity"); break;}
      if (u[0]==0xff800000) { printf("negative infinity"); break;}
      if (u[0]==0x00800000) { printf("+smallest_norm");break;}
      if (u[0]==0x80800000) { printf("-smallest_norm");break;}
      if (u[0]==0x7f7fffff) { printf("+largest_norm");break;}
      if (u[0]==0xff7ffffF) { printf("-largest_norm");break;}
      if (u[0]==0x00000001) { printf("+smallest_denorm");break;}
      if (u[0]==0x80000001) { printf("-smallest_denorm");break;}
      if (u[0]==0x007fffff) { printf("+largest_denorm");break;}
      if (u[0]==0x807fffff) { printf("-largest_denorm");break;}
      if (x[0]==x[0])       { printf("%e",x[0]);break;}
      if ((u[0]&0x80400000)==0x80400000)  { printf("-qNaN");break;}
      if ((u[0]&0x00400000)==0x00400000)  { printf("+qNaN");break;}
      if (u[0]&0x80000000)  { printf("-sNaN");break;}
                              printf("+sNaN");break;
  }
  printf ("(0x%08x)",u[0]);
}

// prinf HP value
static void printfHP(float16_t* x)
{
  uint16_t* u=(uint16_t*)x;
  while(1)
  {
      if (u[0]==0)          { printf("positive zero"); break;}
      if (u[0]==0x8000)     { printf("negative zero"); break;}
      if (u[0]==0x3c00)     { printf("positive one"); break;}
      if (u[0]==0xbc00)     { printf("negative one"); break;}
      if (u[0]==0x7c00)     { printf("positive infinity"); break;}
      if (u[0]==0xfc00)     { printf("negative infinity"); break;}
      if (u[0]==0x0400)     { printf("+smallest_norm");break;}
      if (u[0]==0x8400)     { printf("-smallest_norm");break;}
      if (u[0]==0x7bff)     { printf("+largest_norm");break;}
      if (u[0]==0xfbff)     { printf("-largest_norm");break;}
      if (u[0]==0x0001)     { printf("+smallest_denorm");break;}
      if (u[0]==0x8001)     { printf("-smallest_denorm");break;}
      if (u[0]==0x03ff)     { printf("+largest_denorm");break;}
      if (u[0]==0x83ff)     { printf("-largest_denorm");break;}
      if (eq_f16(x[0],x[0])){ printf("%e",conv_f16_to_f32(x[0]));break;}
      if ((u[0]&0x8200)==0x8200)  { printf("-qNaN");break;}
      if ((u[0]&0x0200)==0x0200)  { printf("+qNaN");break;}
      if (u[0]&0x8000)     { printf("-sNaN");break;}
                              printf("+sNaN");break;
  }
  printf ("(0x%04x)",u[0]);
}

#if 0
// prinf int value
static void printfInt32(int32_t* x)
{
    uint32_t* u=(uint32_t*)x;
    while(1)
    {
        if (u[0]==0)        { printf("positive zero"); break;}
        if (x[0]==1)        { printf("positive one"); break;}
        if (x[0]==-1)       { printf("negative one"); break;}
        printf("%d",x[0]);break;
    }
    printf ("(0x%08x)",u[0]);
}
static void printfInt16(int16_t* x)
{
    uint16_t* u=(uint16_t*)x;
    while(1)
    {
        if (u[0]==0)        { printf("positive zero"); break;}
        if (x[0]==1)        { printf("positive one"); break;}
        if (x[0]==-1)       { printf("negative one"); break;}
        printf("%d",x[0]);break;
    }
    printf ("(0x%04x)",(int)u[0]);
}
#endif
// prinf errno value
static void printfErrno()
{
  switch(errno)
  {
  case EDOM: printf("EDOM"); break;
  case ERANGE: printf("ERANGE"); break;
  }
}

static void printfExceptions()
{
	int excepts;
	if (GET_ISA_OPT(HAVE_FP))
	{
		excepts = FETESTEXCEPT(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
		if (excepts & FE_INVALID)   { excepts &= ~FE_INVALID  ; printf("FE_INVALID%s", excepts ? "," : ""); }
		if (excepts & FE_DIVBYZERO) { excepts &= ~FE_DIVBYZERO; printf("FE_DIVBYZERO%s", excepts ? "," : ""); }
		if (excepts & FE_OVERFLOW)  { excepts &= ~FE_OVERFLOW ; printf("FE_OVERFLOW%s", excepts ? "," : ""); }
		if (excepts & FE_UNDERFLOW) { excepts &= ~FE_UNDERFLOW; printf("FE_UNDERFLOW%s", excepts ? "," : ""); }
		if (excepts & FE_INEXACT)   { excepts &= ~FE_INEXACT  ; printf("FE_INEXACT%s", excepts ? "," : ""); }
	}
}

/*
    print input/output,exception and errno
*/
static void printExceptionsErrno(void* Z, void* X, size_t sz)
{
  switch(sz)
  {
  case 4:
      printf("\t");
      printfSP((float32_t*)X);
      printf("\t");
      printfSP((float32_t*)Z);
      printf("\t");
      printfErrno();
      printf("\t");
      printfExceptions();
      printf("\n");
      break;
  case 2: 
      printf("\t");
      printfHP((float16_t*)X);
      printf("\t");
      printfHP((float16_t*)Z);
      printf("\t");
      printfErrno();
      printf("\t");
      printfExceptions();
      printf("\n");
      break;
  default: 
      ASSERT(0); // not implemented yet
  }
}

static void printExceptionsErrno2(void* Z, void* X,void* Y, size_t sz)
{
  switch(sz)
  {
  case 2:
      printf("\t");
      printfHP((float16_t*)X); printf(", "); printfHP((float16_t*)Y);
      printf("\t");
      printfHP((float16_t*)Z);
      printf("\t");
      printfErrno();
      printf("\t");
      printfExceptions();
      printf("\n");
      break;
  case 4:
      printf("\t");
      printfSP((float32_t*)X); printf(", "); printfSP((float32_t*)Y);
      printf("\t");
      printfSP((float32_t*)Z);
      printf("\t");
      printfErrno();
      printf("\t");
      printfExceptions();
      printf("\n");
      break;
  default: ASSERT(0); // not implemented yet
  }
}

#if 0
static void printExceptionsErrnoLdexp(void* Z, void* X, size_t szX, void* Y,size_t szY)
{
    switch((szX<<4)+szY)
    {
    case (4<<4)+4:
        printf("\t");
        printfSP((float32_t*)X); printf(", "); printfInt32((int32_t*)Y);
        printf("\t");
        printfSP((float32_t*)Z);
        printf("\t");
        printfErrno();
        printf("\t");
        printfExceptions();
        printf("\n");
        break;
    case (2<<4)+2:
        printf("\t");
        printfHP((float16_t*)X); printf(", "); printfInt16((int16_t*)Y);
        printf("\t");
        printfHP((float16_t*)Z);
        printf("\t");
        printfErrno();
        printf("\t");
        printfExceptions();
        printf("\n");
        break;
    default: ASSERT(0); // not implemented yet
    }
}
#endif
static tFPStat fpstat;

/* default create/destroy functions */
static int createFn(tTestEngContext* context)
{
    FPStatCreate(&fpstat);
    return 1;
}

/* default create/destroy functions */
static int destroyFn(tTestEngContext* context)
{
    fprintf(stdout,"\n");
    FPStatPrint(&fpstat,stdout);
    FPStatClose(&fpstat);
    return 1;
}


/* Allocate vectors and load the data set:
 * vector/scalar X (in), vector/scalar Y (in), vector/scalar Z (out)
   vector V is used for accuracy tests
 */
static int te_loadFxn_vXvYvZacc_helper( tTestEngContext * context, 
                              int fmtX, int nElemX, 
                              int fmtY, int nElemY, 
                              int fmtZ, int nElemZ,
                              int fmtV, int nElemV)  
{
  int res = 0;

  ASSERT( context && context->seqFile );

  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  /* Allocate input data vector X. */
  if ( !vecAlloc( &context->dataSet.X, nElemX, context->desc->isAligned, fmtX, 0 ) )
  {
    printf( "te_loadFxn_vXvYvZacc_helper(): failed to allocate vector X; "
            "fmtX = 0x%02x, nElemX = %d\n", (unsigned)fmtX, nElemX );
  }
  /* Allocate input data vector Y. */
  else if ( !vecAlloc( &context->dataSet.Y, nElemY, context->desc->isAligned, fmtY, 0 ) )
  {
    printf( "te_loadFxn_vXvYvZacc_helper(): failed to allocate vector Y; "
            "fmtY = 0x%02x, nElemY = %d\n", (unsigned)fmtY, nElemY );
  }
  /* Allocate input data vector V. */
  else if ( !vecAlloc( &context->dataSet.V, nElemV, context->desc->isAligned, fmtV, 0 ) )
  {
    printf( "te_loadFxn_vXvYvZacc_helper(): failed to allocate vector Z; "
            "fmtZ = 0x%02x, nElemZ = %d\n", (unsigned)fmtZ, nElemZ );
  }
  /* Allocate output/reference data vectors memory. */
  else if ( 3 != vecsAlloc( context->desc->isAligned, fmtZ,
                           &context->dataSet.Z  , nElemZ,
                           &context->dataSet.Zlo, nElemZ,
                           &context->dataSet.Zhi, nElemZ, 0 ) )
  {
    printf( "te_loadFxn_vXvYvZacc_helper(): failed to allocate vectors Zlo/Zhi; "
            "fmtZ = 0x%02x, nElemZ = %d\n", (unsigned)fmtZ, nElemZ );
  }
  /* Load vectors data from the SEQ-file. */
  else if ( !seqFileReadVecs( context->seqFile,
                             &context->dataSet.X,
                             &context->dataSet.Y,
                             &context->dataSet.Zlo,
                             &context->dataSet.Zhi, 
                             &context->dataSet.V,
                             0 ) )
  {
    printf( "te_loadFxn_vXvYvZacc_helper(): failed to read vectors data; "
            "fmtX = 0x%02x, nElemX = %d, fmtY = 0x%02x, nElemY = %d, fmtZ = 0x%02x, nElemZ = %d\n",
            (unsigned)fmtX, nElemX, (unsigned)fmtY, nElemY, (unsigned)fmtZ, nElemZ );
  }
  else
  {
    res = 1;
  }

  /* Free vectors data if failed. */
  if ( !res ) te_freeVectors(context);
  return (res);
} /* te_loadFxn_vXvYvZacc_helper() */


/* Allocate vectors and load the data set:
 * vector X (in), vector Z (out) */
static int te_loadFxn_vXvZacc( tTestEngContext * context )
{
  int M, N, L;
  int nElem, fmt;
  ASSERT( context && context->seqFile );
  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );
  nElem = M*N*L;
  fmt = context->desc->fmt;
  return te_loadFxn_vXvYvZacc_helper( context, fmt, nElem, fmt, 0, fmt, nElem, (fmt&~FMT_DTYPE_MASK)|FMT_FLOAT64,nElem );
}

/* Allocate vectors and load the data set:
 * vector X (in), vector Y (in), vector Z (out) */
static int te_loadFxn_vXvYvZacc( tTestEngContext * context )
{
  int M, N, L;
  int fmt, nElem;
  ASSERT( context && context->seqFile );
  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );
  nElem = M*N*L;
  fmt = context->desc->fmt;
  return te_loadFxn_vXvYvZacc_helper( context, fmt, nElem, fmt, nElem, fmt, nElem, (fmt&~FMT_DTYPE_MASK)|FMT_FLOAT64,nElem );
} /* te_loadFxn_vXvYvZacc() */

/* Apply a function to the test case data set:
 * scalar functions with single argument, e.g. cos() 
   
   function also update accuracy statistics
 */
static void te_math_processFxn_s_vXvZ_acc(tTestEngContext * context)
{
  typedef float16_t tFxn_fl16( float16_t );
  typedef float32_t tFxn_fl32( float32_t );

  #define CALL_FXN( typeFxn, Fxn, typeXZ, X, Z, n )           \
      { te_errh_resetStates( context );                       \
        *(typeXZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeXZ*)(X) ); \
        te_errh_verifyStates( context, n );                   \
        (X) = (typeXZ*)(X) + 1; (Z) = (typeXZ*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Z;
  int n, N;

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );
  fxn = context->target.fut;
  N = context->args.dim[0];

  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FLOAT16: for ( n=0; n<N; n++ ) 
      CALL_FXN( tFxn_fl16, fxn, float16_t, X, Z, n ); 
      FPStatAdd_fp16(&fpstat, testCaseTypeStr[context->args.caseType], vecGetElem_fl64( &context->dataSet.V, 0 ), vecGetElem_fl16( &context->dataSet.Z, 0 ), N);
      break;
  case FMT_FLOAT32: 
      for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, float32_t, X, Z, n ); 
      FPStatAdd(&fpstat, testCaseTypeStr[context->args.caseType], vecGetElem_fl64( &context->dataSet.V, 0 ), vecGetElem_fl32( &context->dataSet.Z, 0 ), N);
      break;
  default: ASSERT( 0 );
  }
  #undef CALL_FXN
} /* te_math_processFxn_s_vXvZ_acc() */

/* Apply a function to the test case data set:
 * scalar atan2(y,x) */
static void te_math_processFxn_s_atan2_acc(tTestEngContext * context)
{
  typedef float32_t tFxn_fl32( float32_t, float32_t );
  typedef float16_t tFxn_fl16( float16_t, float16_t );

  #define CALL_FXN( typeFxn, Fxn, typeXYZ, Y, X, Z, n )                            \
      { te_errh_resetStates( context );                                            \
        *(typeXYZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeXYZ*)(Y), *(typeXYZ*)(X) );    \
        te_errh_verifyStates( context, n );                                        \
        (Y) = (typeXYZ*)(Y) + 1; (X) = (typeXYZ*)(X) + 1; (Z) = (typeXYZ*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Y, *Z;
  int n, N;

  te_vReportStd(context);
  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  fxn = context->target.fut;
  N = context->args.dim[0];

  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FLOAT16: 
      for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl16, fxn, float16_t, Y, X, Z, n ); 
      FPStatAdd_fp16(&fpstat, testCaseTypeStr[context->args.caseType], vecGetElem_fl64( &context->dataSet.V, 0 ), vecGetElem_fl16( &context->dataSet.Z, 0 ), N);
      break;
  case FMT_FLOAT32: 
      for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, float32_t, Y, X, Z, n ); 
      FPStatAdd(&fpstat, testCaseTypeStr[context->args.caseType], vecGetElem_fl64( &context->dataSet.V, 0 ), vecGetElem_fl32( &context->dataSet.Z, 0 ), N);
      break;
  default: ASSERT( 0 );
  }
  #undef CALL_FXN
} /* te_math_processFxn_s_atan2_acc() */

/* Apply a function to the test case data set:
 * scalar functions with two arguments, e.g. powf() */
static void te_math_processFxn_s_vXvYvZ_acc(tTestEngContext * context)
{
  typedef float32_t tFxn_fl32( float32_t, float32_t );
  typedef float16_t tFxn_fl16( float16_t, float16_t );

  #define CALL_FXN( typeFxn, Fxn, typeXYZ, X, Y, Z, n )                         \
      { te_errh_resetStates( context );                                         \
        *(typeXYZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeXYZ*)(X), *(typeXYZ*)(Y) ); \
        te_errh_verifyStates( context, n );                                     \
        (X) = (typeXYZ*)(X) + 1; (Y) = (typeXYZ*)(Y) + 1; (Z) = (typeXYZ*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Y, *Z;
  int n, N;

  te_vReportStd(context);
  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  fxn = context->target.fut;
  N = context->args.dim[0];

  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FLOAT16: 
      for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl16, fxn, float16_t, X, Y, Z, n ); 
      FPStatAdd_fp16(&fpstat, testCaseTypeStr[context->args.caseType], vecGetElem_fl64( &context->dataSet.V, 0 ), vecGetElem_fl16( &context->dataSet.Z, 0 ), N);
      break;
  case FMT_FLOAT32: 
      for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, float32_t, X, Y, Z, n ); 
      FPStatAdd(&fpstat, testCaseTypeStr[context->args.caseType], vecGetElem_fl64( &context->dataSet.V, 0 ), vecGetElem_fl32( &context->dataSet.Z, 0 ), N);
      break;
  default: ASSERT( 0 );
  }
  #undef CALL_FXN
} /* te_math_processFxn_s_vXvYvZ_acc() */

/* Apply a function to the test case data set:
 * scalar functions with single argument, e.g. cos() 
   
   function also prints exception outputs
 */
static void te_math_processFxn_s_vXvZ_exc(tTestEngContext * context)
{
  typedef float16_t tFxn_fl16( float16_t );
  typedef float32_t tFxn_fl32( float32_t );

  #define CALL_FXN( typeFxn, Fxn, typeXZ, X, Z, n )           \
      { FECLEAREXCEPT(FE_ALL_EXCEPT); \
        te_errh_resetStates( context );                       \
        *(typeXZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeXZ*)(X) ); \
        te_errh_verifyStates( context, n );                   \
        printExceptionsErrno((typeXZ*)(Z),(typeXZ*)(X),sizeof(typeXZ));      \
        (X) = (typeXZ*)(X) + 1; (Z) = (typeXZ*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Z;
  int n, N;

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );
  fxn = context->target.fut;
  N = context->args.dim[0];
  printf("\n");
  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FLOAT16: 
      for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl16, fxn, float16_t, X, Z, n ); 
      break;
  case FMT_FLOAT32: 
      for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, float32_t, X, Z, n ); 
      break;
  default: ASSERT( 0 );
  }
  #undef CALL_FXN
}
static void te_math_processFxn_s_atan2_exc(tTestEngContext * context)
{
  typedef float16_t tFxn_fl16( float16_t, float16_t );
  typedef float32_t tFxn_fl32( float32_t, float32_t );

  #define CALL_FXN( typeFxn, Fxn, typeXYZ, Y, X, Z, n )                            \
      { FECLEAREXCEPT(FE_ALL_EXCEPT); \
        te_errh_resetStates( context );                                            \
        *(typeXYZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeXYZ*)(Y), *(typeXYZ*)(X) );    \
        te_errh_verifyStates( context, n );                                        \
        printExceptionsErrno2((typeXYZ*)(Z),(typeXYZ*)(X),(typeXYZ*)(Y),sizeof(typeXYZ));      \
        (Y) = (typeXYZ*)(Y) + 1; (X) = (typeXYZ*)(X) + 1; (Z) = (typeXYZ*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Y, *Z;
  int n, N;

  te_vReportStd(context);
  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  fxn = context->target.fut;
  N = context->args.dim[0];
  printf("\n");
  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FLOAT16: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl16, fxn, float16_t, Y, X, Z, n ); break;
  case FMT_FLOAT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, float32_t, Y, X, Z, n ); break;
  default: ASSERT( 0 );
  }
  #undef CALL_FXN
}

/* Apply a function to the test case data set:
* scalar functions with two arguments, e.g. powf() */
static void te_math_processFxn_s_vXvYvZ_exc(tTestEngContext * context)
{
	typedef float16_t tFxn_fl16(float16_t, float16_t);
	typedef float32_t tFxn_fl32(float32_t, float32_t);

#define CALL_FXN( typeFxn, Fxn, typeXYZ, X, Y, Z, n )                         \
      { FECLEAREXCEPT(FE_ALL_EXCEPT); \
        te_errh_resetStates( context );                                         \
        *(typeXYZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeXYZ*)(X), *(typeXYZ*)(Y) ); \
        te_errh_verifyStates( context, n );                                     \
        printExceptionsErrno2((typeXYZ*)(Z),(typeXYZ*)(X),(typeXYZ*)(Y),sizeof(typeXYZ));      \
        (X) = (typeXYZ*)(X) + 1; (Y) = (typeXYZ*)(Y) + 1; (Z) = (typeXYZ*)(Z) + 1; }

	tTestEngTarget  fxn;
	void *X, *Y, *Z;
	int n, N;

	te_vReportStd(context);
	X = vecGetElem(&context->dataSet.X, 0);
	Y = vecGetElem(&context->dataSet.Y, 0);
	Z = vecGetElem(&context->dataSet.Z, 0);

	fxn = context->target.fut;
	N = context->args.dim[0];
	printf("\n");
	switch (context->desc->fmt & FMT_DTYPE_MASK)
	{
	case FMT_FLOAT16: for (n = 0; n<N; n++) CALL_FXN(tFxn_fl16, fxn, float16_t, X, Y, Z, n); break;
	case FMT_FLOAT32: for (n = 0; n<N; n++) CALL_FXN(tFxn_fl32, fxn, float32_t, X, Y, Z, n); break;
	default: ASSERT(0);
	}
#undef CALL_FXN
} /* te_math_processFxn_s_vXvYvZ_exc() */

/* Vector Mathematics API test definitions. */

/* Perform all accuracy tests for Vector Mathematics API functions. */
int acc_maths2(int isFull, int isVerbose, int breakOnError, int optAccuracy, int optException)
{
  int res = 1;

#define DO_TEST_ACC(fxn, seqFile, extraFlags)                                                                   \
  if ( res || !breakOnError ) res &= ( 0 != testExec( testDefTblAccuracy,sizeof(testDefTblAccuracy)/sizeof(testDefTblAccuracy[0]), \
                                                      (tTestEngTarget)(fxn), "math2/" seqFile,                   \
                                                      (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE,             \
                                                      isFull, isVerbose, breakOnError ) )

#define DO_TEST_EXC(fxn, seqFile, extraFlags)                                                                   \
  if ( res || !breakOnError ) res &= ( 0 != testExec( testDefTblExc,sizeof(testDefTblExc)/sizeof(testDefTblExc[0]), \
                                                      (tTestEngTarget)(fxn), "math2/" seqFile,                   \
                                                      (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE,             \
                                                      isFull, isVerbose, breakOnError ) )

  // tests for accuracy
  if (optAccuracy)
  {
    /*
    * Stage 2
    */
    DO_TEST_ACC( &scl_log2f         , "accuracy_scl_log2f.seq"           , 0 );
    DO_TEST_ACC( &scl_lognf         , "accuracy_scl_lognf.seq"           , 0 );
    DO_TEST_ACC( &scl_log10f        , "accuracy_scl_log10f.seq"          , 0 );
    DO_TEST_ACC( &scl_antilog2f     , "accuracy_scl_antilog2f.seq"       , 0 );
    DO_TEST_ACC( &scl_antilog10f    , "accuracy_scl_antilog10f.seq"      , 0 );
    DO_TEST_ACC( &scl_antilognf     , "accuracy_scl_antilognf.seq"       , 0 );
    DO_TEST_ACC( &scl_powf          , "accuracy_scl_powf.seq"            , 0 );
    DO_TEST_ACC( &scl_sinef         , "accuracy_scl_sinef.seq"           , 0 );
    DO_TEST_ACC( &scl_cosinef       , "accuracy_scl_cosinef.seq"         , 0 );
    DO_TEST_ACC( &scl_tanf          , "accuracy_scl_tanf.seq"            , 0 );
    DO_TEST_ACC( &scl_atan2f        , "accuracy_scl_atan2f.seq"          , 0 );
    DO_TEST_ACC( &scl_atanf         , "accuracy_scl_atanf.seq"           , 0 );
    DO_TEST_ACC( &scl_tanhf         , "accuracy_scl_tanhf.seq"           , 0 );
  }
  // tests for exceptions
  if (optException)
  {
    /*
    * Stage 2
    */
    DO_TEST_EXC( &scl_log2f         , "exceptions_scl_log2f.seq"           , 0);
    DO_TEST_EXC( &scl_lognf         , "exceptions_scl_lognf.seq"           , 0);
    DO_TEST_EXC( &scl_log10f        , "exceptions_scl_log10f.seq"          , 0);
    DO_TEST_EXC( &scl_antilog2f     , "exceptions_scl_antilog2f.seq"       , 0);
    DO_TEST_EXC( &scl_antilog10f    , "exceptions_scl_antilog10f.seq"      , 0);
    DO_TEST_EXC( &scl_antilognf     , "exceptions_scl_antilognf.seq"       , 0);
    DO_TEST_EXC( &scl_powf          , "exceptions_scl_powf.seq"            , 0);
    DO_TEST_EXC( &scl_sinef         , "exceptions_scl_sinef.seq"           , 0);
    DO_TEST_EXC( &scl_cosinef       , "exceptions_scl_cosinef.seq"         , 0);
    DO_TEST_EXC( &scl_tanf          , "exceptions_scl_tanf.seq"            , 0);
    DO_TEST_EXC( &scl_atan2f        , "exceptions_scl_atan2f.seq"          , 0);
    DO_TEST_EXC( &scl_atanf         , "exceptions_scl_atanf.seq"           , 0);
    DO_TEST_EXC( &scl_tanhf         , "exceptions_scl_tanhf.seq"           , 0);
  }
  return (res);
}

/* Test executive function. Performs the specified test on a brief or full version
* of the designated SEQ-file. Return the test result (non-zero if passed). */
static int testExec(const tTestDef* pTestDef, size_t szTestDef,
                    tTestEngTarget   targetFxn, const char * seqName,
                    int errhExtendedTest, int isFull, int isVerbose, int breakOnError )
{
  int tblIx, funcIx;
  for (tblIx = 0; tblIx<(int)szTestDef; tblIx++)
  {
    for (funcIx = 0; funcIx<MAX_FUNC_NUM; funcIx++)
    {
      if (targetFxn == pTestDef[tblIx].funcList[funcIx])
      {
        tTestEngDesc testDesc = pTestDef[tblIx].testDesc;
        testDesc.extraParam = (uint32_t)errhExtendedTest;

        return (TestEngRun(targetFxn, &testDesc,
          seqName, isFull,
          isVerbose, breakOnError, 0));
      }
    }
  }
  ASSERT(!"Test not defined");
  return (0);
}
