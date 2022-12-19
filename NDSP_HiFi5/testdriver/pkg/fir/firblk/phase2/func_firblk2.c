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
#include "../common/test_firblk.h"

static tFirOldDescr     api_bkfirf        ={{(tFirOldFxnAlloc*)bkfirf_alloc,        (tFirOldFxnInit*)bkfirf_init,        (tFirOldFxnProcess*)bkfirf_process        ,(tFirOldextIrFxnAllocExtIr*)bkfirf_allocExtIr ,(tFirOldextIrFxnCopyExtIr*)bkfirf_copyExtIr }};
static tFirOldDescr     api_bkfiraf       ={{(tFirOldFxnAlloc*)bkfiraf_alloc,       (tFirOldFxnInit*)bkfiraf_init,       (tFirOldFxnProcess*)bkfiraf_process       }};
static tFirOldDescr     api_cxfirf        ={{(tFirOldFxnAlloc*)cxfirf_alloc,        (tFirOldFxnInit*)cxfirf_init,        (tFirOldFxnProcess*)cxfirf_process        ,(tFirOldextIrFxnAllocExtIr*)cxfirf_allocExtIr ,(tFirOldextIrFxnCopyExtIr*)cxfirf_copyExtIr }};
static tFirStereoDescr  api_stereo_bkfirf ={{(tFirStereoFxnAlloc*)stereo_bkfirf_alloc,        (tFirStereoFxnInit*)stereo_bkfirf_init,        (tFirStereoFxnProcess*)stereo_bkfirf_process        ,(tFirStereoextIrFxnAllocExtIr*)stereo_bkfirf_allocExtIr ,(tFirStereoextIrFxnCopyExtIr*)stereo_bkfirf_copyExtIr }};
static const tTestEngDesc descr_bkfirf          = { FMT_REAL | FMT_FLOAT32, TE_FIR_FIR|TE_FIR_OLDEXTIR  , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_bkfiraf         = { FMT_REAL | FMT_FLOAT32, TE_FIR_FIR|TE_FIR_OLDGENERIC, NULL, TE_DIM_NUM_1, TE_ALIGN_NO , te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_cxfirf          = { FMT_CPLX | FMT_FLOAT32, TE_FIR_FIR|TE_FIR_OLDEXTIR  , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_stereo_bkfirf          = { FMT_REAL | FMT_FLOAT32, TE_FIR_FIR|TE_FIR_OLDEXTIR             , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_stereo, te_destroy_fir_stereo, &te_loadFxn_fir_stereo, &te_processFxn_fir_stereo };

static const tTbl tests[] =
{
  /* Block Real FIR Filter */
  { 2, &descr_bkfirf, (tTestEngTarget)&api_bkfirf, 1,"firblk2/bkfirf_lpf1.seq"},
  { 2, &descr_bkfirf, (tTestEngTarget)&api_bkfirf, 1,"firblk2/bkfirf_bpf1.seq"},
  { 2, &descr_bkfirf, (tTestEngTarget)&api_bkfirf, 1,"firblk2/bkfirf_8taps.seq"},
  { 2, &descr_bkfirf, (tTestEngTarget)&api_bkfirf, 0,"firblk2/bkfirf_lpf2.seq"},
  { 2, &descr_bkfirf, (tTestEngTarget)&api_bkfirf, 0,"firblk2/bkfirf_hpf1.seq"},
  { 2, &descr_bkfirf, (tTestEngTarget)&api_bkfirf, 0,"firblk2/bkfirf_bpf2.seq"},
  { 2, &descr_bkfirf, (tTestEngTarget)&api_bkfirf, 0,"firblk2/bkfirf_4taps.seq"},
 
  { 2, &descr_bkfiraf, (tTestEngTarget)&api_bkfiraf, 1,"firblk2/bkfiraf_lpf1.seq"},
  { 2, &descr_bkfiraf, (tTestEngTarget)&api_bkfiraf, 1,"firblk2/bkfiraf_bpf1.seq"},
  { 2, &descr_bkfiraf, (tTestEngTarget)&api_bkfiraf, 1,"firblk2/bkfiraf_3taps.seq"},
  { 2, &descr_bkfiraf, (tTestEngTarget)&api_bkfiraf, 0,"firblk2/bkfiraf_lpf2.seq"},
  { 2, &descr_bkfiraf, (tTestEngTarget)&api_bkfiraf, 0,"firblk2/bkfiraf_hpf1.seq"},
  { 2, &descr_bkfiraf, (tTestEngTarget)&api_bkfiraf, 0,"firblk2/bkfiraf_bpf2.seq"},
  { 2, &descr_bkfiraf, (tTestEngTarget)&api_bkfiraf, 0,"firblk2/bkfiraf_4taps.seq"},
  { 2, &descr_bkfiraf, (tTestEngTarget)&api_bkfiraf, 0,"firblk2/bkfiraf_2taps.seq"},
  { 2, &descr_bkfiraf, (tTestEngTarget)&api_bkfiraf, 0,"firblk2/bkfiraf_1tap.seq"},
 
  { 2, &descr_stereo_bkfirf, (tTestEngTarget)&api_stereo_bkfirf, 1,"firblk2/stereo_bkfirf_lpf1.seq"},
  { 2, &descr_stereo_bkfirf, (tTestEngTarget)&api_stereo_bkfirf, 1,"firblk2/stereo_bkfirf_bpf1.seq"},
  { 2, &descr_stereo_bkfirf, (tTestEngTarget)&api_stereo_bkfirf, 1,"firblk2/stereo_bkfirf_tap8.seq"},
  { 2, &descr_stereo_bkfirf, (tTestEngTarget)&api_stereo_bkfirf, 0,"firblk2/stereo_bkfirf_lpf2.seq"},
  { 2, &descr_stereo_bkfirf, (tTestEngTarget)&api_stereo_bkfirf, 0,"firblk2/stereo_bkfirf_hpf1.seq"},
  { 2, &descr_stereo_bkfirf, (tTestEngTarget)&api_stereo_bkfirf, 0,"firblk2/stereo_bkfirf_bpf2.seq"},
  { 2, &descr_stereo_bkfirf, (tTestEngTarget)&api_stereo_bkfirf, 0,"firblk2/stereo_bkfirf_tap4.seq"},

  /* Block Complex FIR Filter */
  { 2, &descr_cxfirf, (tTestEngTarget)&api_cxfirf,1,"firblk2/cxfirf_lpf1.seq"},
  { 2, &descr_cxfirf, (tTestEngTarget)&api_cxfirf,1,"firblk2/cxfirf_bpf1.seq"},
  { 2, &descr_cxfirf, (tTestEngTarget)&api_cxfirf,1,"firblk2/cxfirf_lpf2.seq"},
  { 2, &descr_cxfirf, (tTestEngTarget)&api_cxfirf,0,"firblk2/cxfirf_hpf1.seq"},
  { 2, &descr_cxfirf, (tTestEngTarget)&api_cxfirf,0,"firblk2/cxfirf_bpf2.seq"},
  { 2, &descr_cxfirf, (tTestEngTarget)&api_cxfirf,0,"firblk2/cxfirf_4taps.seq"},
};

/* Perform all tests for FIR API functions. */
int func_firblk2(int isFull, int isVerbose, int breakOnError)
{
    return test_firblk(2, isFull, isVerbose, breakOnError, tests, (int)(sizeof(tests)/sizeof(tests[0])));
}
