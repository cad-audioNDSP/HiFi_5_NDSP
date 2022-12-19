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

#if 0 //HiFi3/3z API
static size_t bqriir24x24_df1_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR24X24_DF1_SCRATCH_SIZE( N, M ) ; }
static size_t bqriir24x24_df2_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR24X24_DF2_SCRATCH_SIZE( N, M ) ; }
static const tIirOldDescr api_bqriir24x24_df1={(tIirOldFxnAlloc*)bqriir24x24_df1_alloc,(tIirOldFxnInit*)bqriir24x24_df1_init,bqriir24x24_df1_getScratch,(tIirOldFxnProcess*)bqriir24x24_df1};
static const tIirOldDescr api_bqriir24x24_df2={(tIirOldFxnAlloc*)bqriir24x24_df2_alloc,(tIirOldFxnInit*)bqriir24x24_df2_init,bqriir24x24_df2_getScratch,(tIirOldFxnProcess*)bqriir24x24_df2};
static const tTestEngDesc descr_bqriir24x24_df1 = { FMT_REAL | FMT_FRACT32, TE_IIR_DF1,                NULL,TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_old, te_destroy_iir_old, &te_loadFxn_iir_old, &te_processFxn_iir_old };
static const tTestEngDesc descr_bqriir24x24_df2 = { FMT_REAL | FMT_FRACT32, TE_IIR_DF2,                NULL,TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_old, te_destroy_iir_old, &te_loadFxn_iir_old, &te_processFxn_iir_old };
#endif
static size_t bqriir16x16_df1_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR16X16_DF1_SCRATCH_SIZE( N, M ) ; }
static size_t bqriir16x16_df2_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR16X16_DF2_SCRATCH_SIZE( N, M ) ; }

static size_t bqriir32x16_df1_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR32X16_DF1_SCRATCH_SIZE( N, M ) ; }
static size_t bqriir32x16_df2_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR32X16_DF2_SCRATCH_SIZE( N, M ) ; }
static size_t bqriir32x32_df1_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR32X32_DF1_SCRATCH_SIZE( N, M ) ; }
static size_t bqriir32x32_df2_getScratch(int N, int M) { (void)N; (void)M; return BQRIIR32X32_DF2_SCRATCH_SIZE( N, M ) ; }
static size_t stereo_bqriir16x16_df1_getScratch(int N, int M) { (void)N; (void)M; return STEREO_BQRIIR16X16_DF1_SCRATCH_SIZE( N, M ) ; }
static size_t stereo_bqriir32x16_df1_getScratch(int N, int M) { (void)N; (void)M; return STEREO_BQRIIR32X16_DF1_SCRATCH_SIZE( N, M ) ; }
static size_t stereo_bqriir32x32_df1_getScratch(int N, int M) { (void)N; (void)M; return STEREO_BQRIIR32X32_DF1_SCRATCH_SIZE( N, M ) ; }

static const tIirDescr api_bqriir16x16_df1 = { (tIirFxnAlloc*)bqriir16x16_df1_alloc, (tIirFxnInit*)bqriir16x16_df1_init, bqriir16x16_df1_getScratch, (tIirFxnProcess*)bqriir16x16_df1, (tIirFxnDelay*)bqriir16x16_df1_groupDelay };
static const tIirDescr api_bqriir16x16_df2 = { (tIirFxnAlloc*)bqriir16x16_df2_alloc, (tIirFxnInit*)bqriir16x16_df2_init, bqriir16x16_df2_getScratch, (tIirFxnProcess*)bqriir16x16_df2, (tIirFxnDelay*)bqriir16x16_df2_groupDelay };
static const tIirDescr api_bqriir32x16_df1 = { (tIirFxnAlloc*)bqriir32x16_df1_alloc, (tIirFxnInit*)bqriir32x16_df1_init, bqriir32x16_df1_getScratch, (tIirFxnProcess*)bqriir32x16_df1, (tIirFxnDelay*)bqriir32x16_df1_groupDelay };
static const tIirDescr api_bqriir32x16_df2 = { (tIirFxnAlloc*)bqriir32x16_df2_alloc, (tIirFxnInit*)bqriir32x16_df2_init, bqriir32x16_df2_getScratch, (tIirFxnProcess*)bqriir32x16_df2, (tIirFxnDelay*)bqriir32x16_df2_groupDelay };
static const tIirDescr api_bqriir32x32_df1 = { (tIirFxnAlloc*)bqriir32x32_df1_alloc, (tIirFxnInit*)bqriir32x32_df1_init, bqriir32x32_df1_getScratch, (tIirFxnProcess*)bqriir32x32_df1, (tIirFxnDelay*)bqriir32x32_df1_groupDelay };
static const tIirDescr api_bqriir32x32_df2 = { (tIirFxnAlloc*)bqriir32x32_df2_alloc, (tIirFxnInit*)bqriir32x32_df2_init, bqriir32x32_df2_getScratch, (tIirFxnProcess*)bqriir32x32_df2, (tIirFxnDelay*)bqriir32x32_df2_groupDelay };
static const tIirStereoDescr api_stereo_bqriir16x16_df1 = { (tIirStereoFxnAlloc*)stereo_bqriir16x16_df1_alloc, (tIirStereoFxnInit*)stereo_bqriir16x16_df1_init, stereo_bqriir16x16_df1_getScratch, (tIirStereoFxnProcess*)stereo_bqriir16x16_df1, (tIirStereoFxnDelay*)stereo_bqriir16x16_df1_groupDelay };
static const tIirStereoDescr api_stereo_bqriir32x16_df1 = { (tIirStereoFxnAlloc*)stereo_bqriir32x16_df1_alloc, (tIirStereoFxnInit*)stereo_bqriir32x16_df1_init, stereo_bqriir32x16_df1_getScratch, (tIirStereoFxnProcess*)stereo_bqriir32x16_df1, (tIirStereoFxnDelay*)stereo_bqriir32x16_df1_groupDelay };
static const tIirStereoDescr api_stereo_bqriir32x32_df1 = { (tIirStereoFxnAlloc*)stereo_bqriir32x32_df1_alloc, (tIirStereoFxnInit*)stereo_bqriir32x32_df1_init, stereo_bqriir32x32_df1_getScratch, (tIirStereoFxnProcess*)stereo_bqriir32x32_df1, (tIirStereoFxnDelay*)stereo_bqriir32x32_df1_groupDelay };

static const tTestEngDesc descr_bqriir16x16_df1 = { FMT_REAL | FMT_FRACT16, TE_IIR_DF1 | TE_IIR_16X16, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir16x16_df2 = { FMT_REAL | FMT_FRACT16, TE_IIR_DF2 | TE_IIR_16X16, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir32x16_df1 = { FMT_REAL | FMT_FRACT16, TE_IIR_DF1, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir32x16_df2 = { FMT_REAL | FMT_FRACT16, TE_IIR_DF2, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir32x32_df1 = { FMT_REAL | FMT_FRACT32, TE_IIR_DF1, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriir32x32_df2 = { FMT_REAL | FMT_FRACT32, TE_IIR_DF2, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_stereo_bqriir16x16_df1 = { FMT_REAL | FMT_FRACT16, TE_IIR_DF1 | TE_IIR_16X16, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_stereo, te_destroy_iir_stereo, &te_loadFxn_iir_stereo, &te_processFxn_iir_stereo };
static const tTestEngDesc descr_stereo_bqriir32x16_df1 = { FMT_REAL | FMT_FRACT16, TE_IIR_DF1, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_stereo, te_destroy_iir_stereo, &te_loadFxn_iir_stereo, &te_processFxn_iir_stereo };
static const tTestEngDesc descr_stereo_bqriir32x32_df1 = { FMT_REAL | FMT_FRACT32, TE_IIR_DF1, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_stereo, te_destroy_iir_stereo, &te_loadFxn_iir_stereo, &te_processFxn_iir_stereo };

static const tTbl tests[] =
{
  //tests with delay imitation 
  { 1, &descr_bqriir16x16_df1    , (tTestEngTarget)&api_bqriir16x16_df1    , 1, "iirbq1/bqriir16x16_df1_lpf1.seq" },
  { 1, &descr_bqriir16x16_df1    , (tTestEngTarget)&api_bqriir16x16_df1    , 1, "iirbq1/bqriir16x16_df1_bpf1.seq" },
  { 1, &descr_bqriir16x16_df1    , (tTestEngTarget)&api_bqriir16x16_df1    , 0, "iirbq1/bqriir16x16_df1_bpf2.seq" },
  { 1, &descr_bqriir16x16_df1    , (tTestEngTarget)&api_bqriir16x16_df1    , 0, "iirbq1/bqriir16x16_df1_bsf1.seq" },
  { 1, &descr_bqriir16x16_df1    , (tTestEngTarget)&api_bqriir16x16_df1    , 0, "iirbq1/bqriir16x16_df1_hpf1.seq" },

  { 1, &descr_bqriir16x16_df2    , (tTestEngTarget)&api_bqriir16x16_df2    , 1, "iirbq1/bqriir16x16_df2_lpf1.seq" },
  { 1, &descr_bqriir16x16_df2    , (tTestEngTarget)&api_bqriir16x16_df2    , 1, "iirbq1/bqriir16x16_df2_bpf1.seq" },
  { 1, &descr_bqriir16x16_df2    , (tTestEngTarget)&api_bqriir16x16_df2    , 0, "iirbq1/bqriir16x16_df2_bpf2.seq" },
  { 1, &descr_bqriir16x16_df2    , (tTestEngTarget)&api_bqriir16x16_df2    , 0, "iirbq1/bqriir16x16_df2_bsf1.seq" },
  { 1, &descr_bqriir16x16_df2    , (tTestEngTarget)&api_bqriir16x16_df2    , 0, "iirbq1/bqriir16x16_df2_hpf1.seq" },

  { 1, &descr_bqriir32x16_df1    , (tTestEngTarget)&api_bqriir32x16_df1    , 1, "iirbq1/bqriir32x16_df1_lpf1.seq" },
  { 1, &descr_bqriir32x16_df1    , (tTestEngTarget)&api_bqriir32x16_df1    , 1, "iirbq1/bqriir32x16_df1_bpf1.seq" },
  { 1, &descr_bqriir32x16_df1    , (tTestEngTarget)&api_bqriir32x16_df1    , 0, "iirbq1/bqriir32x16_df1_bpf2.seq" },
  { 1, &descr_bqriir32x16_df1    , (tTestEngTarget)&api_bqriir32x16_df1    , 0, "iirbq1/bqriir32x16_df1_bsf1.seq" },
  { 1, &descr_bqriir32x16_df1    , (tTestEngTarget)&api_bqriir32x16_df1    , 0, "iirbq1/bqriir32x16_df1_hpf1.seq" },
  { 1, &descr_bqriir32x16_df1    , (tTestEngTarget)&api_bqriir32x16_df1    , 0, "iirbq1/bqriir32x16_df1_hpf2.seq" },

  { 1, &descr_bqriir32x16_df2    , (tTestEngTarget)&api_bqriir32x16_df2    , 1, "iirbq1/bqriir32x16_df2_lpf1.seq" },
  { 1, &descr_bqriir32x16_df2    , (tTestEngTarget)&api_bqriir32x16_df2    , 1, "iirbq1/bqriir32x16_df2_bpf1.seq" },
  { 1, &descr_bqriir32x16_df2    , (tTestEngTarget)&api_bqriir32x16_df2    , 0, "iirbq1/bqriir32x16_df2_bpf2.seq" },
  { 1, &descr_bqriir32x16_df2    , (tTestEngTarget)&api_bqriir32x16_df2    , 0, "iirbq1/bqriir32x16_df2_bsf1.seq" },
  { 1, &descr_bqriir32x16_df2    , (tTestEngTarget)&api_bqriir32x16_df2    , 0, "iirbq1/bqriir32x16_df2_hpf1.seq" },
  { 1, &descr_bqriir32x16_df2    , (tTestEngTarget)&api_bqriir32x16_df2    , 0, "iirbq1/bqriir32x16_df2_hpf2.seq" },
#if 0 //HiFi3/3z API
  { 1, &descr_bqriir24x24_df1    , (tTestEngTarget)&api_bqriir24x24_df1    , 1, "iirbq1/bqriir24x24_df1_lpf1.seq" },
  { 1, &descr_bqriir24x24_df1    , (tTestEngTarget)&api_bqriir24x24_df1    , 1, "iirbq1/bqriir24x24_df1_bpf1.seq" },
  { 1, &descr_bqriir24x24_df1    , (tTestEngTarget)&api_bqriir24x24_df1    , 0, "iirbq1/bqriir24x24_df1_bpf2.seq" },
  { 1, &descr_bqriir24x24_df1    , (tTestEngTarget)&api_bqriir24x24_df1    , 0, "iirbq1/bqriir24x24_df1_bsf1.seq" },
  { 1, &descr_bqriir24x24_df1    , (tTestEngTarget)&api_bqriir24x24_df1    , 0, "iirbq1/bqriir24x24_df1_hpf1.seq" },

  { 1, &descr_bqriir24x24_df2    , (tTestEngTarget)&api_bqriir24x24_df2    , 1, "iirbq1/bqriir24x24_df2_lpf1.seq" },
  { 1, &descr_bqriir24x24_df2    , (tTestEngTarget)&api_bqriir24x24_df2    , 1, "iirbq1/bqriir24x24_df2_bpf1.seq" },
  { 1, &descr_bqriir24x24_df2    , (tTestEngTarget)&api_bqriir24x24_df2    , 0, "iirbq1/bqriir24x24_df2_bpf2.seq" },
  { 1, &descr_bqriir24x24_df2    , (tTestEngTarget)&api_bqriir24x24_df2    , 0, "iirbq1/bqriir24x24_df2_bsf1.seq" },
  { 1, &descr_bqriir24x24_df2    , (tTestEngTarget)&api_bqriir24x24_df2    , 0, "iirbq1/bqriir24x24_df2_hpf1.seq" },
#endif
  { 1, &descr_bqriir32x32_df1    , (tTestEngTarget)&api_bqriir32x32_df1    , 1, "iirbq1/bqriir32x32_df1_lpf1.seq" },
  { 1, &descr_bqriir32x32_df1    , (tTestEngTarget)&api_bqriir32x32_df1    , 1, "iirbq1/bqriir32x32_df1_bpf1.seq" },
  { 1, &descr_bqriir32x32_df1    , (tTestEngTarget)&api_bqriir32x32_df1    , 0, "iirbq1/bqriir32x32_df1_bpf2.seq" },
  { 1, &descr_bqriir32x32_df1    , (tTestEngTarget)&api_bqriir32x32_df1    , 0, "iirbq1/bqriir32x32_df1_bsf1.seq" },
  { 1, &descr_bqriir32x32_df1    , (tTestEngTarget)&api_bqriir32x32_df1    , 0, "iirbq1/bqriir32x32_df1_bsf2.seq" },
  { 1, &descr_bqriir32x32_df1    , (tTestEngTarget)&api_bqriir32x32_df1    , 0, "iirbq1/bqriir32x32_df1_hpf1.seq" },

  { 1, &descr_bqriir32x32_df2    , (tTestEngTarget)&api_bqriir32x32_df2    , 1, "iirbq1/bqriir32x32_df2_lpf1.seq" },
  { 1, &descr_bqriir32x32_df2    , (tTestEngTarget)&api_bqriir32x32_df2    , 1, "iirbq1/bqriir32x32_df2_bpf1.seq" },
  { 1, &descr_bqriir32x32_df2    , (tTestEngTarget)&api_bqriir32x32_df2    , 0, "iirbq1/bqriir32x32_df2_bpf2.seq" },
  { 1, &descr_bqriir32x32_df2    , (tTestEngTarget)&api_bqriir32x32_df2    , 0, "iirbq1/bqriir32x32_df2_bsf1.seq" },
  { 1, &descr_bqriir32x32_df2    , (tTestEngTarget)&api_bqriir32x32_df2    , 0, "iirbq1/bqriir32x32_df2_bsf2.seq" },
  { 1, &descr_bqriir32x32_df2    , (tTestEngTarget)&api_bqriir32x32_df2    , 0, "iirbq1/bqriir32x32_df2_hpf1.seq" },

  { 1, &descr_stereo_bqriir16x16_df1    , (tTestEngTarget)&api_stereo_bqriir16x16_df1    , 1, "iirbq1/stereo_bqriir16x16_df1_lpf1.seq" },
  { 1, &descr_stereo_bqriir16x16_df1    , (tTestEngTarget)&api_stereo_bqriir16x16_df1    , 1, "iirbq1/stereo_bqriir16x16_df1_bpf1.seq" },
  { 1, &descr_stereo_bqriir16x16_df1    , (tTestEngTarget)&api_stereo_bqriir16x16_df1    , 0, "iirbq1/stereo_bqriir16x16_df1_bpf2.seq" },
  { 1, &descr_stereo_bqriir16x16_df1    , (tTestEngTarget)&api_stereo_bqriir16x16_df1    , 0, "iirbq1/stereo_bqriir16x16_df1_bsf1.seq" },
  { 1, &descr_stereo_bqriir16x16_df1    , (tTestEngTarget)&api_stereo_bqriir16x16_df1    , 0, "iirbq1/stereo_bqriir16x16_df1_hpf1.seq" },

  { 1, &descr_stereo_bqriir32x16_df1    , (tTestEngTarget)&api_stereo_bqriir32x16_df1    , 1, "iirbq1/stereo_bqriir32x16_df1_lpf1.seq" },
  { 1, &descr_stereo_bqriir32x16_df1    , (tTestEngTarget)&api_stereo_bqriir32x16_df1    , 1, "iirbq1/stereo_bqriir32x16_df1_bpf1.seq" },
  { 1, &descr_stereo_bqriir32x16_df1    , (tTestEngTarget)&api_stereo_bqriir32x16_df1    , 0, "iirbq1/stereo_bqriir32x16_df1_bpf2.seq" },
  { 1, &descr_stereo_bqriir32x16_df1    , (tTestEngTarget)&api_stereo_bqriir32x16_df1    , 0, "iirbq1/stereo_bqriir32x16_df1_bsf1.seq" },
  { 1, &descr_stereo_bqriir32x16_df1    , (tTestEngTarget)&api_stereo_bqriir32x16_df1    , 0, "iirbq1/stereo_bqriir32x16_df1_hpf1.seq" },

  { 1, &descr_stereo_bqriir32x32_df1    , (tTestEngTarget)&api_stereo_bqriir32x32_df1    , 1, "iirbq1/stereo_bqriir32x32_df1_lpf1.seq" },
  { 1, &descr_stereo_bqriir32x32_df1    , (tTestEngTarget)&api_stereo_bqriir32x32_df1    , 1, "iirbq1/stereo_bqriir32x32_df1_bpf1.seq" },
  { 1, &descr_stereo_bqriir32x32_df1    , (tTestEngTarget)&api_stereo_bqriir32x32_df1    , 0, "iirbq1/stereo_bqriir32x32_df1_bpf2.seq" },
  { 1, &descr_stereo_bqriir32x32_df1    , (tTestEngTarget)&api_stereo_bqriir32x32_df1    , 0, "iirbq1/stereo_bqriir32x32_df1_bsf1.seq" },
  { 1, &descr_stereo_bqriir32x32_df1    , (tTestEngTarget)&api_stereo_bqriir32x32_df1    , 0, "iirbq1/stereo_bqriir32x32_df1_hpf1.seq" },
};

int func_iirbq1(int isFull, int isVerbose, int breakOnError)
{
    return main_iirbq(1, isFull, isVerbose, breakOnError, tests, (int)(sizeof(tests)/sizeof(tests[0])));
}
