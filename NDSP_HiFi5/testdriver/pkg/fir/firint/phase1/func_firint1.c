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

#if 0// for HiFi3/3z
static tFirOldDescr api_firinterp24x24  ={{(tFirOldFxnAlloc*)firinterp24x24_alloc,  (tFirOldFxnInit*)firinterp24x24_init,  (tFirOldFxnProcess*)firinterp24x24_process  }};
static const tTestEngDesc descr_firinterp24x24  = { FMT_REAL | FMT_FRACT32, TE_FIR_UP|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
#endif
static tFirOldDescr api_firinterp16x16  ={{(tFirOldFxnAlloc*)firinterp16x16_alloc,  (tFirOldFxnInit*)firinterp16x16_init,  (tFirOldFxnProcess*)firinterp16x16_process  }};
static tFirOldDescr api_firinterp32x16  ={{(tFirOldFxnAlloc*)firinterp32x16_alloc,  (tFirOldFxnInit*)firinterp32x16_init,  (tFirOldFxnProcess*)firinterp32x16_process  }};
static tFirOldDescr api_firinterp32x32  ={{(tFirOldFxnAlloc*)firinterp32x32_alloc,  (tFirOldFxnInit*)firinterp32x32_init,  (tFirOldFxnProcess*)firinterp32x32_process  }};
static tFirOldDescr api_firinterp32x32ep={{(tFirOldFxnAlloc*)firinterp32x32ep_alloc,(tFirOldFxnInit*)firinterp32x32ep_init,(tFirOldFxnProcess*)firinterp32x32ep_process}};

static const tTestEngDesc descr_firinterp16x16  = { FMT_REAL | FMT_FRACT16, TE_FIR_UP|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_firinterp32x16  = { FMT_REAL | FMT_FRACT32, TE_FIR_UP|TE_FIR_OLDDECIMATOR|TE_FIR_FILTER_32X16, NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_firinterp32x32  = { FMT_REAL | FMT_FRACT32, TE_FIR_UP|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_firinterp32x32ep= { FMT_REAL | FMT_FRACT32, TE_FIR_UP|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };

static const tTbl tests[] =
{
  { 1, &descr_firinterp16x16, (tTestEngTarget)&api_firinterp16x16,1,"firint1/firinterp16x16_up2x.seq" },
  { 1, &descr_firinterp16x16, (tTestEngTarget)&api_firinterp16x16,1,"firint1/firinterp16x16_up3x.seq" },
  { 1, &descr_firinterp16x16, (tTestEngTarget)&api_firinterp16x16,1,"firint1/firinterp16x16_up6x.seq" },
  { 1, &descr_firinterp16x16, (tTestEngTarget)&api_firinterp16x16,0,"firint1/firinterp16x16_up4x.seq" },
  { 1, &descr_firinterp16x16, (tTestEngTarget)&api_firinterp16x16,0,"firint1/firinterp16x16_up5x.seq" },

  { 1, &descr_firinterp32x16, (tTestEngTarget)&api_firinterp32x16,1,"firint1/firinterp32x16_up2x.seq" },
  { 1, &descr_firinterp32x16, (tTestEngTarget)&api_firinterp32x16,1,"firint1/firinterp32x16_up3x.seq" },
  { 1, &descr_firinterp32x16, (tTestEngTarget)&api_firinterp32x16,1,"firint1/firinterp32x16_up6x.seq" },
  { 1, &descr_firinterp32x16, (tTestEngTarget)&api_firinterp32x16,0,"firint1/firinterp32x16_up4x.seq" },
  { 1, &descr_firinterp32x16, (tTestEngTarget)&api_firinterp32x16,0,"firint1/firinterp32x16_up5x.seq" },
#if 0// for HiFi3/3z
  { 1, &descr_firinterp24x24, (tTestEngTarget)&api_firinterp24x24,1,"firinterp24x24_up2x.seq" },
  { 1, &descr_firinterp24x24, (tTestEngTarget)&api_firinterp24x24,1,"firinterp24x24_up3x.seq" },
  { 1, &descr_firinterp24x24, (tTestEngTarget)&api_firinterp24x24,1,"firinterp24x24_up6x.seq" },
  { 1, &descr_firinterp24x24, (tTestEngTarget)&api_firinterp24x24,0,"firinterp24x24_up4x.seq" },
  { 1, &descr_firinterp24x24, (tTestEngTarget)&api_firinterp24x24,0,"firinterp24x24_up5x.seq" },
#endif
  { 1, &descr_firinterp32x32, (tTestEngTarget)&api_firinterp32x32,1,"firint1/firinterp32x32_up2x.seq" },
  { 1, &descr_firinterp32x32, (tTestEngTarget)&api_firinterp32x32,1,"firint1/firinterp32x32_up3x.seq" },
  { 1, &descr_firinterp32x32, (tTestEngTarget)&api_firinterp32x32,1,"firint1/firinterp32x32_up6x.seq" },
  { 1, &descr_firinterp32x32, (tTestEngTarget)&api_firinterp32x32,0,"firint1/firinterp32x32_up4x.seq" },
  { 1, &descr_firinterp32x32, (tTestEngTarget)&api_firinterp32x32,0,"firint1/firinterp32x32_up5x.seq" },

  { 1, &descr_firinterp32x32ep, (tTestEngTarget)&api_firinterp32x32ep,1,"firint1/firinterp32x32ep_up2x.seq" },
  { 1, &descr_firinterp32x32ep, (tTestEngTarget)&api_firinterp32x32ep,1,"firint1/firinterp32x32ep_up3x.seq" },
  { 1, &descr_firinterp32x32ep, (tTestEngTarget)&api_firinterp32x32ep,1,"firint1/firinterp32x32ep_up6x.seq" },
  { 1, &descr_firinterp32x32ep, (tTestEngTarget)&api_firinterp32x32ep,0,"firint1/firinterp32x32ep_up4x.seq" },
  { 1, &descr_firinterp32x32ep, (tTestEngTarget)&api_firinterp32x32ep,0,"firint1/firinterp32x32ep_up5x.seq" },
};

/* Perform all tests for FIR API functions. */
int func_firint1(int isFull, int isVerbose, int breakOnError)
{
    return test_firint(1, isFull, isVerbose, breakOnError, tests, (int)(sizeof(tests)/sizeof(tests[0])));
}
