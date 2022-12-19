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
 * Test procedures for Kalman Filter functions
 */
#include "../common/test_iirkal.h"

static const t_kalmanupd1_api kalmanupd1f_api      = { &kalmanupd1f_getScratchSize      };

/* Test descriptions */
static const struct {
    tTestEngTarget funcList[MAX_FUNC_NUM];
    tTestEngDesc   testDesc;
} testDefTbl[] = {
    { FUNC_LIST( (tTestEngTarget)&kalmanupd1f     ), TEST_DESC( FMT_REAL|FMT_FLOAT32, 0, &kalmanupd1f_api     , TE_DIM_NUM_1, TE_ALIGN_YES,  &te_loadFxn_kalmanupd1, &te_processFxn_kalmanupd1) },
    { FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, 0, 0, NULL, NULL ) } /* End of table */
};

/* Perform all tests for Kalman Filter API functions. */
int func_iirkal2(int isFull, int isVerbose, int breakOnError)
{
    int res = 1;
    #define DO_TEST(fxn, seqFile)                                                                    \
        if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                               \
                                                           sizeof(testDefTbl)/sizeof(testDefTbl[0]), \
                                                           MAX_FUNC_NUM,                             \
                                                           (tTestEngTarget)(fxn), "iirkal2/" seqFile,    \
                                                           isFull, isVerbose, breakOnError ) )
    DO_TEST(&kalmanupd1f, "kalmanupd1f.seq");
    return res;
}
