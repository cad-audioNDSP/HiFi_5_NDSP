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
 * Test engine API
 */

#ifndef __TESTENG_H
#define __TESTENG_H

/* Portable data types */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"
#include "vreport.h"

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------------------------*
 *                Forward declarations of Test Engine structures                *
 *------------------------------------------------------------------------------*/

/* Test definition: target function + test attributes + test means. */
typedef struct tagTestEngDesc tTestEngDesc;
/* Test engine operation context. */
typedef struct tagTestEngContext tTestEngContext;

/*------------------------------------------------------------------------------*
 *  "Virtual" functions used by the Test Engine throughout the test procedure   *
 *------------------------------------------------------------------------------*/

/* Function pointer. */
typedef void (*te_fun_ptr_t)();

/* Pointer to a function(s) under test. */
typedef te_fun_ptr_t tTestEngTarget;
/* Create a target algorithm instance and set tTestEngContext::target fields.
 * Return zero if failed. */
typedef int tTestEngTargetCreateFxn( tTestEngContext * context );
/* Destroy the target algorithm instance and free memory block(s) allocated
 * for the target object. Return zero whenever any violation of memory area
 * bounds is detected. */
typedef int tTestEngTargetDestroyFxn( tTestEngContext * context );

/* Allocate in/out vectors for the next test case, and load the data set
 * from the SEQ-file. Return zero if failed. */
typedef int  tTestEngLoadFxn( tTestEngContext * context );
/* Apply the target function to the test case data set. */
typedef void tTestEngProcessFxn( tTestEngContext * context );

/*------------------------------------------------------------------------------*
 *                        Test engine operation context                         *
 *------------------------------------------------------------------------------*/

struct tagTestEngContext
{
  const tTestEngDesc * desc;      /* Test definition                             */
  char                 seqDir[192];/* Path to SEQ-file location                  */
  tSeqFile             seqFile;   /* SEQ-file with test data                     */
  int                  isFull;    /* Full test if nonzero, otherwise brief test  */
  int                  isVerbose; /* Print additional info on test progress      */
  int                  doRetake;  /* If appears non-zero after                   *
                                   * testCaseProcessFxn(), the test engine       *
                                   * replays the load/process/validate sequence  *
                                   * within the scope of the current test case.  *
                                   * This allows to break a large test case      *
                                   * into a few subcases.                        */

  struct {                        /* Test target spec                            */
    tTestEngTarget     fut;       /* Function Under Test                         */
    void *             handle;    /* User-defined handle (optional)              */
  } target;

  struct {                        /* Test case arguments:                        */
    int caseNum;                  /*  - test case number                         */
    int caseType;                 /*  - semantic type (see testCaseTypeStr)      */
    int dim[5];                   /*  - dimensional arguments                    */
  } args;

  struct {                        /* Test case data set:                         */
    tVec X, Y;                    /*  - input vectors/scalars                    */
    tVec U, V;                    /*  - additional input/scratch vectors/scalars */
    tVec Z, W;                    /*  - output vectors/scalars                   */
    tVec Zlo, Zhi, Wlo, Whi;      /*  - output range boundaries                  */
    tVec auxVec[8];               /*  - auxiliary vectors for exceptional        *
                                   *    cases when there is not enough named     *
                                   *    vectors to hold all data entities        */
  } dataSet;

  struct {                        /* Error handling verfication support:         */
    int isEnabled;                /*  - verification enable flag                 */
    int isPassed;                 /*  - verification result for a test case      *
                                   *    (this may involve multiple invocations   *
                                   *    of a scalar FUT, or a single call of a   *
                                   *    vector FUT).                             */
    struct {                      /*  - reference error states:                  */
      tVec edom;                  /*    - errno=EDOM assertion indices           */
      tVec erange;                /*    - errno=ERANGE asseertion indices        */
      tVec fe_inv;                /*    - FE_INVALID exception raise indices     */
      tVec fe_divz;               /*    - FE_DIVBYZERO exception raise indices   */
      tVec fe_ovfl;               /*    - FE_OVERFLOW exception raise indices    */
    } refStates;                    /*    All indices refer to in/out data vector *
                                     *    positions.                              */
    int entranceExcepts;            /*  - exception flags set just beofore a FUT  *
                                     *    is called (for preserved state check)   */
    int exceptEnable;               /*  - exception enable controls (FUT must not *
                                     *    distort the exc. enable controls)       */
  } errh;
};

/*------------------------------------------------------------------------------*
 *      User-supplied test description and test engine executive function       *
 *------------------------------------------------------------------------------*/

/* Number of dimensional arguments of a test case (M,N,etc.) */
#define TE_DIM_NUM_1   1
#define TE_DIM_NUM_2   2
#define TE_DIM_NUM_3   3
#define TE_DIM_NUM_4   4
#define TE_DIM_NUM_5   5

/* Data vectors alignment options */
#define TE_ALIGN_YES   1
#define TE_ALIGN_NO    0

/* Extra parameter for tests with error handling verification: whether or not
 * the Test Engine verifies that a FUT does not clear any FP exception flags. */
#define TE_ERRH_EXTENDED_TEST_ENABLE   1
#define TE_ERRH_EXTENDED_TEST_DISABLE  0
/* Do not check particular floating-point exception flags. */
#define TE_ERRH_IGNORE_FE_INEXACT      2
#define TE_ERRH_IGNORE_FE_UNDERFLOW    4

#define TE_VECSCL_BITEXACTNESS  0x10    /* flag to check bitexactness 
                                           For checking the bitexacness, extraPtr should point to 
                                           the table containing the mapping between vector functions 
                                           and their scalar counterparts (tTestEngVecSclTbl)
                                           table should be ended with NULL pointers
                                        */
typedef struct tagTestEngVecSclTbl
{
    tTestEngProcessFxn*   testCaseProcessFxn; /* testcase processing function for scalar version */
    te_fun_ptr_t          scl;                /* scalar function                                 */
    te_fun_ptr_t          vec;                /* vector function                                 */
}
tTestEngVecSclTbl;

/* flags in tagTestEngDesc::extraParam controlling overlapping of input/outputs */
#define OVLP    0xFF000000
#define OVLP_XZ 0x01000000    /* overlap X and Z */
#define OVLP_YZ 0x02000000    /* overlap Y and Z */
#define OVLP_XW 0x04000000    /* overlap X and W */
#define OVLP_YW 0x08000000    /* overlap Y and W */

/* Test definition: target function + test attributes + test means. All-zero entry marks the
 * end of descriptions table. */
struct tagTestEngDesc
{
  int                        fmt;                /* Data format, a combination of FMT_x symbols (see vectool.h)   */
  uint32_t                   extraParam;         /* Extra constant parameter of function, usually zero.           *
                                                  * some standard values are defined in OVLP_xxx.                 *
                                                  * 3 LSBs are used for error handling verification control.      */
  const void *               extraPtr;           /* Auxiliary pointer to be used by extensions, e.g. TestEngFir.  */
  int                        dimNum;             /* The number of dimensional arguments of a test case (M,N,etc.) */
  int                        isAligned;          /* Use aligned data buffers.                                     */
  tTestEngTargetCreateFxn  * targetCreateFxn;    /* Target algorithm instance creation (optional)                 */
  tTestEngTargetDestroyFxn * targetDestroyFxn;   /* Target algorithm instance deletion (optional)                 */
  tTestEngLoadFxn          * testCaseLoadFxn;    /* Prepares the data set for a test case.                        */
  tTestEngProcessFxn       * testCaseProcessFxn; /* Applies the target function to the data set.                  */
};

/* Test executive function. Performs the specified test on a brief or full version
 * of the designated SEQ-file. Return the test result (non-zero if passed). */
int TestEngRun( tTestEngTarget targetFxn, const tTestEngDesc * desc,
                const char * seqName, int isFull, int isVerbose, int breakOnError, int testBitexactness );

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
             int isFull, int isVerbose, int breakOnError );

/* function returns total data size allocated by input/output data for given test case */
size_t te_vGetDataSize(const tTestEngContext * context);

/* standard verfication reporting function 
   should be called if FUT is written to context->target.fut
   and there is no specific variants of calling function 
*/
void te_vReportStd(tTestEngContext * context);

/* helper functions for load/process callbacks */
int te_freeVectors( tTestEngContext * context );

/*------------------------------------------------------------------------------*
 *            Repertory of common test case load/process functions              *
 *------------------------------------------------------------------------------*/

/* Allocate vectors and load the data set: */
tTestEngLoadFxn    te_loadFxn_vXvZ;            /* vector X (in), vector Z (out)                                            */
tTestEngLoadFxn    te_loadFxn_vXvZf;           /* vector X (in), vector Z (out, single precision floating-point)           */
tTestEngLoadFxn    te_loadFxn_vXvZd;           /* vector X (in), vector Z (out, double precision floating-point)           */
tTestEngLoadFxn    te_loadFxn_vXvZh;           /* vector X (in), vector Z (out, half precision floating-point)             */
tTestEngLoadFxn    te_loadFxn_vX64vZ;          /* vector int64 X (in), vector Z (out)                                      */
tTestEngLoadFxn    te_loadFxn_vXcvZ;           /* vector complex X (in), vector Z (out)                                    */
tTestEngLoadFxn    te_loadFxn_vXsZ;            /* vector X (in), scalar Z (out)                                            */
tTestEngLoadFxn    te_loadFxn_vXsZ16;          /* vector X (in), scalar int16 Z (out)                                      */
tTestEngLoadFxn    te_loadFxn_vXsZ32;          /* vector X (in), scalar int32 Z (out)                                      */
tTestEngLoadFxn    te_loadFxn_vX32sZ;          /* vector int32 X (in), scalar Z (out)                                      */
tTestEngLoadFxn    te_loadFxn_vX32vZ16;        /* 32-bit vector X (in), 16-bit vector Z (out)                              */
tTestEngLoadFxn    te_loadFxn_vX16vZ32;        /* 16-bit vector X (in), 32-bit vector Z (out)                              */
tTestEngLoadFxn    te_loadFxn_vXvZr;           /* real/complex vector X (in), real vector Z (out)                          */
tTestEngLoadFxn    te_loadFxn_vXvZc;           /* real/complex vector X (in), complex vector Z (out)                       */
tTestEngLoadFxn    te_loadFxn_vXvYvZ;          /* vector X (in), vector Y (in), vector Z (out)                             */
tTestEngLoadFxn    te_loadFxn_vXvYsZ;          /* vector X (in), vector Y (in), scalar Z (out)                             */
tTestEngLoadFxn    te_loadFxn_vXsYvZ;          /* vector X (in), scalar Y (in), vector Z (out)                             */
tTestEngLoadFxn    te_loadFxn_vXsY32vZ;        /* vector X (in), scalar int32 Y (in), vector Z (out)                       */
tTestEngLoadFxn    te_loadFxn_vX32sY32vZ;      /* vector int32 X (in), scalar int32 Y (in), vector Z (out)                 */
tTestEngLoadFxn    te_loadFxn_vXsY32vZ32;      /* vector X (in), scalar int32 Y (in), vector int32 Z (out)                 */
tTestEngLoadFxn    te_loadFxn_vXsYrvZ;         /* vector X (in), real scalar Y (in), vector Z (out)                        */
tTestEngLoadFxn    te_loadFxn_vXsYr_sZ16;      /* vector X (in), real scalar Y (in), scalar int16 Z (out)                  */
tTestEngLoadFxn    te_loadFxn_vXsYr_sZ32;      /* vector X (in), real scalar Y (in), scalar int32 Z (out)                  */
tTestEngLoadFxn    te_loadFxn_vXsYr_vZr;       /* real/complex vector X (in), real scalar Y (in), real vector Z (out)      */
tTestEngLoadFxn    te_loadFxn_vXsYr_vZc;       /* real/complex vector X (in), real scalar Y (in), complex vector Z (out)   */
tTestEngLoadFxn    te_loadFxn_vXvYvZf;         /* vector X (in), vector Y (in), float32 vector Z (out)                     */
tTestEngLoadFxn    te_loadFxn_vXvYvZd;         /* vector X (in), vector Y (in), float64 vector Z (out)                     */
tTestEngLoadFxn    te_loadFxn_vX64vYvZ;        /* vector int64 X (in), vector Y (in), vector Z (out)                       */
tTestEngLoadFxn    te_loadFxn_vXsY32sZ64;      /* vector X (in), scalar int32 Y (in), scalar int64 Z (out)                 */
tTestEngLoadFxn    te_loadFxn_vXvYsZ32;        /* vector X (in), vector Y (in), scalar int32 Z (out)                       */
tTestEngLoadFxn    te_loadFxn_vXvYsZ64;        /* vector X (in), vector Y (in), scalar int64 Z (out)                       */
tTestEngLoadFxn    te_loadFxn_vXvY16sZ64;      /* vector X (in), vector int16 Y (in), scalar int64 Z (out)                 */
tTestEngLoadFxn    te_loadFxn_vXvY32sZ64;      /* vector X (in), vector int32 Y (in), scalar int64 Z (out)                 */
tTestEngLoadFxn    te_loadFxn_vXvZvW;          /* vector X (in), vector Z (out), vector W (out)                            */
tTestEngLoadFxn    te_loadFxn_vXvYvZvW;        /* vector X (in), vector Y (in), vector Z (out), vector W (out)             */
tTestEngLoadFxn    te_loadFxn_sXsYvZsW;        /* scalar X (in), scalar Y (in), vector Z (out), scalar W (out)             */
tTestEngLoadFxn    te_loadFxn_vXvYsUr_sZ32;    /* vector X (in), vector Y (in), real scalar U (in), scalar int32 Z (out)   */
tTestEngLoadFxn    te_loadFxn_vXvYsUr_vZc;     /* vector X (in), vector Y (in), real scalar U (in), complex vector Z (out) */
tTestEngLoadFxn    te_loadFxn_vXsYrsUr_vZ;     /* vector X (in), real scalars Y,U (in), vector Z (out)                     */
tTestEngLoadFxn    te_loadFxn_vXsYrsUr_vZsWr;  /* vector X (in), real scalars Y.U (in), vector Z (out), real scale W (out) */
tTestEngLoadFxn    te_loadFxn_vXvYvZdp;        /* Allocate vectors and load the data set: * vector X (in), vector Y (in), vector Z (out, double precision) */

/* Apply the target function to the test case data set, non-streaming vector variants. */
tTestEngProcessFxn te_processFxn_vXsYvZ;       /* vector X (in), scalar Y (in), vector Z (out)                             */
tTestEngProcessFxn te_processFxn_vZvXsY;       /* vector Z (out), vector X (in), scalar Y (in)                             */
tTestEngProcessFxn te_processFxn_vZvXsY32;     /* vector Z (out), vector X (in), scalar int32 Y (in)                       */
tTestEngProcessFxn te_processFxn_vXvYvZ;       /* vector X (in), vector Y (in), vector Z (out)                             */
tTestEngProcessFxn te_processFxn_vZvXvY;       /* vector Z (out), vector X (in), vector Y (in)                             */
tTestEngProcessFxn te_processFxn_vZvWvX;       /* vector Z (out), vector W (out), vector X (in)                            */
tTestEngProcessFxn te_processFxn_vZvWvXvY;     /* vector Z (out), vector W (out), vector X (in), vector Y (in)             */
tTestEngProcessFxn te_processFxn_vYvXvZ;       /* vector Y (in), vector X (in), vector Z (out)                             */
tTestEngProcessFxn te_processFxn_vXvYsZ;       /* vector X (in), vector Y (in), scalar Z (out)                             */
tTestEngProcessFxn te_processFxn_vZvYvX;       /* vector Z (out), vector Y (in), vector X (in)                             */
tTestEngProcessFxn te_processFxn_vXvZ;         /* vector X (in), vector Z (out)                                            */
tTestEngProcessFxn te_processFxn_vZvX;         /* vector Z (out), vector X (in)                                            */
tTestEngProcessFxn te_processFxn_vXsZ;         /* vector X (in), scalar Z (out)                                            */
tTestEngProcessFxn te_processFxn_vXsZ32;       /* vector X (in), scalar int32 Z (out)                                      */
tTestEngProcessFxn te_processFxn_vXvZvW;       /* vector X (in), vector Z (out), vector W (out)                            */
tTestEngProcessFxn te_processFxn_sZ32_vXvYsUr; /* scalar int32 Z (out), vector X (in), vector Y (in), real scalar U (in)   */
tTestEngProcessFxn te_processFxn_vZ_vXvYsUr;   /* vector Z (out), vector X (in), vector Y (in), real scalar U (in)         */
tTestEngProcessFxn te_processFxn_sZ16_vXsYr;   /* scalar int16 Z (out), vector X (in), real scalar Y (in)                  */
tTestEngProcessFxn te_processFxn_sZ32_vXsYr;   /* scalar int32 Z (out), vector X (in), real scalar Y (in)                  */
tTestEngProcessFxn te_processFxn_vZ_vXsY;      /* vector Z (out), vector X (in), scalar Y (in)                             */
tTestEngProcessFxn te_processFxn_sZi_vX;       /* scalar int Z (out), vector X (in)                                        */
tTestEngProcessFxn te_processFxn_vXsY32sZ64;   /* vector X (in), scalar int32 Y (in), scalar int64 Z (out)                 */
tTestEngProcessFxn te_processFxn_sZ32vXvY;     /* scalar int32 Z (out), vector X (in), vector Y (in)                       */
tTestEngProcessFxn te_processFxn_sZ64vXvY;     /* scalar int64 Z (out), vector X (in), vector Y (in)                       */

/* Apply the target function to the test case data set, scalar target functions. */
tTestEngProcessFxn te_processFxn_scl_vZvX;     /* vector Z (out), vector X (in)                                            */

/* Apply the target function to the test case data set, streaming vector variants. */
tTestEngProcessFxn te_processFxn_s_vXvYvZ;     /* vector X (in), vector Y (in), vector Z (out)                             */
tTestEngProcessFxn te_processFxn_s_vXvZ;       /* vector X (in), vector Z (out)                                            */

/*------------------------------------------------------------------------------*
 *             Error Handling (ERRH) verification support functions             *
 *------------------------------------------------------------------------------*/

/* Initialize the ERRH support; is called internally by the Test Engine before 
 * it reads a SEQ-file for the first time. */
void te_errh_init( tTestEngContext * context );
/* Read the SEQ-file to load reference error states for a test case.
 * Return zero if failed. */
int te_errh_loadRef( tTestEngContext * context );
/* Free reference error state vectors. */
void te_errh_freeRef( tTestEngContext * context );
/* Prepare to call the function under test: reset the errno and clear
 * FP exception flags. */
void te_errh_resetStates( tTestEngContext * context );
/* Sample errno and FP exception flags and check them against the reference
 * error states. For a scalar function under test the second argument must hold
 * the current index value for in/out test data vectors. For a vector function,
 * assign it a negative value. */
int te_errh_verifyStates( tTestEngContext * context, int idx );

/*--------------------------------------------------------------
    function returns vectors directory
--------------------------------------------------------------*/
const char* getVectorsDir(int isFull);
#ifdef __cplusplus
};
#endif

#endif /* __TESTENG_H */
