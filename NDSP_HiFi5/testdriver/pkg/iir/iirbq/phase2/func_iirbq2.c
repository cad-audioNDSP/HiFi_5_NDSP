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
 * Test procedures for IIR
 */
#include "../common/test_iirbq.h"

static const tIirDescr api_bqriirf_df1 = { (tIirFxnAlloc*)bqriirf_df1_alloc, (tIirFxnInit*)bqriirf_df1_init, NULL, (tIirFxnProcess*)bqriirf_df1, (tIirFxnDelay*)bqriirf_df1_groupDelay };
static const tIirDescr api_bqriirf_df2 = { (tIirFxnAlloc*)bqriirf_df2_alloc, (tIirFxnInit*)bqriirf_df2_init, NULL, (tIirFxnProcess*)bqriirf_df2, (tIirFxnDelay*)bqriirf_df2_groupDelay };
static const tIirDescr api_bqriirf_df2t = { (tIirFxnAlloc*)bqriirf_df2t_alloc, (tIirFxnInit*)bqriirf_df2t_init, NULL, (tIirFxnProcess*)bqriirf_df2t, (tIirFxnDelay*)bqriirf_df2t_groupDelay };
static const tIirDescr api_bqciirf_df1 = { (tIirFxnAlloc*)bqciirf_df1_alloc, (tIirFxnInit*)bqciirf_df1_init, NULL, (tIirFxnProcess*)bqciirf_df1, (tIirFxnDelay*)bqciirf_df1_groupDelay };
static const tIirStereoDescr api_stereo_bqriirf_df1 = { (tIirStereoFxnAlloc*)stereo_bqriirf_df1_alloc, (tIirStereoFxnInit*)stereo_bqriirf_df1_init, NULL, (tIirStereoFxnProcess*)stereo_bqriirf_df1, (tIirStereoFxnDelay*)stereo_bqriirf_df1_groupDelay };


static const tTestEngDesc descr_bqriirf_df1 = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF1 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriirf_df2 = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF2 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriirf_df2t = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF2T | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqciirf_df1 = { FMT_CPLX | FMT_FLOAT32, TE_IIR_DF1 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };

static const tTestEngDesc descr_stereo_bqriirf_df1 = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF1 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_stereo, te_destroy_iir_stereo, &te_loadFxn_iir_stereo, &te_processFxn_iir_stereo };

static const tTbl tests[] =
{
  /*
  * Stage 2
  */
  { 2, &descr_bqriirf_df1, (tTestEngTarget)&api_bqriirf_df1, 1, "iirbq2/bqriirf_df1_lpf1.seq" },
  { 2, &descr_bqriirf_df1, (tTestEngTarget)&api_bqriirf_df1, 1, "iirbq2/bqriirf_df1_bpf1.seq" },
  { 2, &descr_bqriirf_df1, (tTestEngTarget)&api_bqriirf_df1, 0, "iirbq2/bqriirf_df1_bpf2.seq" },
  { 2, &descr_bqriirf_df1, (tTestEngTarget)&api_bqriirf_df1, 0, "iirbq2/bqriirf_df1_bpf3.seq" },
  { 2, &descr_bqriirf_df1, (tTestEngTarget)&api_bqriirf_df1, 0, "iirbq2/bqriirf_df1_bsf1.seq" },
  { 2, &descr_bqriirf_df1, (tTestEngTarget)&api_bqriirf_df1, 0, "iirbq2/bqriirf_df1_hpf1.seq" },
                                                                
  { 2, &descr_bqriirf_df2, (tTestEngTarget)&api_bqriirf_df2, 1, "iirbq2/bqriirf_df2_lpf1.seq" },
  { 2, &descr_bqriirf_df2, (tTestEngTarget)&api_bqriirf_df2, 1, "iirbq2/bqriirf_df2_bpf1.seq" },
  { 2, &descr_bqriirf_df2, (tTestEngTarget)&api_bqriirf_df2, 0, "iirbq2/bqriirf_df2_bpf2.seq" },
  { 2, &descr_bqriirf_df2, (tTestEngTarget)&api_bqriirf_df2, 0, "iirbq2/bqriirf_df2_bpf3.seq" },
  { 2, &descr_bqriirf_df2, (tTestEngTarget)&api_bqriirf_df2, 0, "iirbq2/bqriirf_df2_bsf1.seq" },
  { 2, &descr_bqriirf_df2, (tTestEngTarget)&api_bqriirf_df2, 0, "iirbq2/bqriirf_df2_hpf1.seq" },

  { 2, &descr_bqriirf_df2t, (tTestEngTarget)&api_bqriirf_df2t, 1, "iirbq2/bqriirf_df2t_lpf1.seq" },
  { 2, &descr_bqriirf_df2t, (tTestEngTarget)&api_bqriirf_df2t, 1, "iirbq2/bqriirf_df2t_bpf1.seq" },
  { 2, &descr_bqriirf_df2t, (tTestEngTarget)&api_bqriirf_df2t, 0, "iirbq2/bqriirf_df2t_bpf2.seq" },
  { 2, &descr_bqriirf_df2t, (tTestEngTarget)&api_bqriirf_df2t, 0, "iirbq2/bqriirf_df2t_bpf3.seq" },
  { 2, &descr_bqriirf_df2t, (tTestEngTarget)&api_bqriirf_df2t, 0, "iirbq2/bqriirf_df2t_bsf1.seq" },
  { 2, &descr_bqriirf_df2t, (tTestEngTarget)&api_bqriirf_df2t, 0, "iirbq2/bqriirf_df2t_hpf1.seq" },

  { 2, &descr_bqciirf_df1, (tTestEngTarget)&api_bqciirf_df1, 1, "iirbq2/bqciirf_df1_lpf1.seq" },
  { 2, &descr_bqciirf_df1, (tTestEngTarget)&api_bqciirf_df1, 1, "iirbq2/bqciirf_df1_bpf1.seq" },
  { 2, &descr_bqciirf_df1, (tTestEngTarget)&api_bqciirf_df1, 0, "iirbq2/bqciirf_df1_bpf2.seq" },
  { 2, &descr_bqciirf_df1, (tTestEngTarget)&api_bqciirf_df1, 0, "iirbq2/bqciirf_df1_bsf1.seq" },
  { 2, &descr_bqciirf_df1, (tTestEngTarget)&api_bqciirf_df1, 0, "iirbq2/bqciirf_df1_hpf1.seq" },

  { 2, &descr_stereo_bqriirf_df1, (tTestEngTarget)&api_stereo_bqriirf_df1, 1, "iirbq2/stereo_bqriirf_df1_lpf1.seq" },
  { 2, &descr_stereo_bqriirf_df1, (tTestEngTarget)&api_stereo_bqriirf_df1, 1, "iirbq2/stereo_bqriirf_df1_bpf1.seq" },
  { 2, &descr_stereo_bqriirf_df1, (tTestEngTarget)&api_stereo_bqriirf_df1, 0, "iirbq2/stereo_bqriirf_df1_bpf2.seq" },
  { 2, &descr_stereo_bqriirf_df1, (tTestEngTarget)&api_stereo_bqriirf_df1, 0, "iirbq2/stereo_bqriirf_df1_bpf3.seq" },
  { 2, &descr_stereo_bqriirf_df1, (tTestEngTarget)&api_stereo_bqriirf_df1, 0, "iirbq2/stereo_bqriirf_df1_bsf1.seq" },
  { 2, &descr_stereo_bqriirf_df1, (tTestEngTarget)&api_stereo_bqriirf_df1, 0, "iirbq2/stereo_bqriirf_df1_hpf1.seq" },
};

int func_iirbq2(int isFull, int isVerbose, int breakOnError)
{
    return main_iirbq(2, isFull, isVerbose, breakOnError, tests, (int)(sizeof(tests)/sizeof(tests[0])));
}
