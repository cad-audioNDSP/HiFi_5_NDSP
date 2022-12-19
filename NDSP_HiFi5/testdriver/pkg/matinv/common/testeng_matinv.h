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
 * Test-engine add-on for matrix inversion categories 
 */
#ifndef TESTENG_MTXINV_H__
#define TESTENG_MTXINV_H__

/* Cross-platform data type definitions. */
#include "types.h"
/* Test engine API. */
#include "testeng.h"

/* extraParam flags */
#define MTXINV_PLAIN  0 /* non-streaming data */
#define MTXINV_STREAM 1 /* streaming data     */

int te_loadFxn_mtxinv(tTestEngContext * context);
int te_loadFxn_mtxgjelim(tTestEngContext * context);
void te_processFxn_matinv( tTestEngContext * context );
void te_processFxn_gjelim( tTestEngContext * context );

/* API definition structure. */
typedef struct 
{
    int    isFloatingPoint;     // 1 for floating point API, 0 for fixed point
    size_t (*getScratchSize)(); 
} tMtxInvApi;

typedef struct 
{
  const tTestEngDesc *pFunDescr;
  tTestEngTarget      fxns;
  int                 isFull;   /* 1 - brief & full, 0 - full only */
  const char*         seqFile;
}
tTestMtxinv;

/* Perform all tests for matrix inversion API functions. */
int te_ExecMtxInv(const tTestMtxinv* pTbl, int isFull, int isVerbose, int breakOnError );

#endif /*TESTENG_MTXINV_H__*/

