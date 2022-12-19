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
 * Test procedures for FIR
 */
#include "../common/test_firint.h"

static tFirOldDescr api_firinterpf          ={{(tFirOldFxnAlloc*)firinterpf_alloc,    (tFirOldFxnInit*)firinterpf_init,    (tFirOldFxnProcess*)firinterpf_process    }};
static const tTestEngDesc descr_firinterpf  = { FMT_REAL | FMT_FLOAT32, TE_FIR_UP|TE_FIR_OLDDECIMATOR             ,NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };

static const tTbl tests[] =
{
  { 2, &descr_firinterpf, (tTestEngTarget)&api_firinterpf,1,"firint2/firinterpf_up2x.seq"},
  { 2, &descr_firinterpf, (tTestEngTarget)&api_firinterpf,1,"firint2/firinterpf_up3x.seq"},
  { 2, &descr_firinterpf, (tTestEngTarget)&api_firinterpf,1,"firint2/firinterpf_up6x.seq"},
  { 2, &descr_firinterpf, (tTestEngTarget)&api_firinterpf,0,"firint2/firinterpf_up4x.seq"},
  { 2, &descr_firinterpf, (tTestEngTarget)&api_firinterpf,0,"firint2/firinterpf_up5x.seq"},
};

/* Perform all tests for FIR API functions. */
int func_firint2(int isFull, int isVerbose, int breakOnError)
{
    return test_firint(2, isFull, isVerbose, breakOnError, tests, (int)(sizeof(tests)/sizeof(tests[0])));
}
