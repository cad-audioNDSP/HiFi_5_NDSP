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
#include "../../iirbq/common/test_iirbq.h"

static const tIirDescr api_bqriirf_df1_nd = { (tIirFxnAlloc*)bqriirf_df1_nd_alloc, (tIirFxnInit*)bqriirf_df1_nd_init, NULL, (tIirFxnProcess*)bqriirf_df1_nd, (tIirFxnDelay*)bqriirf_df1_nd_groupDelay };
static const tIirDescr api_bqriirf_df2_nd = { (tIirFxnAlloc*)bqriirf_df2_nd_alloc, (tIirFxnInit*)bqriirf_df2_nd_init, NULL, (tIirFxnProcess*)bqriirf_df2_nd, (tIirFxnDelay*)bqriirf_df2_nd_groupDelay };
static const tIirDescr api_bqriirf_df2t_nd = { (tIirFxnAlloc*)bqriirf_df2t_nd_alloc, (tIirFxnInit*)bqriirf_df2t_nd_init, NULL, (tIirFxnProcess*)bqriirf_df2t_nd, (tIirFxnDelay*)bqriirf_df2t_nd_groupDelay };
static const tIirDescr api_bqciirf_df1_nd = { (tIirFxnAlloc*)bqciirf_df1_nd_alloc, (tIirFxnInit*)bqciirf_df1_nd_init, NULL, (tIirFxnProcess*)bqciirf_df1_nd, (tIirFxnDelay*)bqciirf_df1_nd_groupDelay };
static const tIirStereoDescr api_stereo_bqriirf_df1_nd = { (tIirStereoFxnAlloc*)stereo_bqriirf_df1_nd_alloc, (tIirStereoFxnInit*)stereo_bqriirf_df1_nd_init, NULL, (tIirStereoFxnProcess*)stereo_bqriirf_df1_nd, (tIirStereoFxnDelay*)stereo_bqriirf_df1_nd_groupDelay };

static const tTestEngDesc descr_bqriirf_df1_nd = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF1 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriirf_df2_nd = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF2 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqriirf_df2t_nd = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF2T | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_bqciirf_df1_nd = { FMT_CPLX | FMT_FLOAT32, TE_IIR_DF1 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir, te_destroy_iir, &te_loadFxn_iir, &te_processFxn_iir };
static const tTestEngDesc descr_stereo_bqriirf_df1_nd = { FMT_REAL | FMT_FLOAT32, TE_IIR_DF1 | TE_IIR_FLOAT, NULL, TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_stereo, te_destroy_iir_stereo, &te_loadFxn_iir_stereo, &te_processFxn_iir_stereo };

static const tTbl tests_nd[] =
{
    //no delay functions
    { 2, &descr_bqriirf_df1_nd, (tTestEngTarget)&api_bqriirf_df1_nd, 1, "iirbqnd2/bqriirf_df1_nd_lpf1.seq" },
    { 2, &descr_bqriirf_df1_nd, (tTestEngTarget)&api_bqriirf_df1_nd, 1, "iirbqnd2/bqriirf_df1_nd_bpf1.seq" },
    { 2, &descr_bqriirf_df1_nd, (tTestEngTarget)&api_bqriirf_df1_nd, 0, "iirbqnd2/bqriirf_df1_nd_bpf2.seq" },
    { 2, &descr_bqriirf_df1_nd, (tTestEngTarget)&api_bqriirf_df1_nd, 0, "iirbqnd2/bqriirf_df1_nd_bpf3.seq" },
    { 2, &descr_bqriirf_df1_nd, (tTestEngTarget)&api_bqriirf_df1_nd, 0, "iirbqnd2/bqriirf_df1_nd_bsf1.seq" },
    { 2, &descr_bqriirf_df1_nd, (tTestEngTarget)&api_bqriirf_df1_nd, 0, "iirbqnd2/bqriirf_df1_nd_hpf1.seq" },

    { 2, &descr_bqriirf_df2_nd, (tTestEngTarget)&api_bqriirf_df2_nd, 1, "iirbqnd2/bqriirf_df2_nd_lpf1.seq" },
    { 2, &descr_bqriirf_df2_nd, (tTestEngTarget)&api_bqriirf_df2_nd, 1, "iirbqnd2/bqriirf_df2_nd_bpf1.seq" },
    { 2, &descr_bqriirf_df2_nd, (tTestEngTarget)&api_bqriirf_df2_nd, 0, "iirbqnd2/bqriirf_df2_nd_bpf2.seq" },
    { 2, &descr_bqriirf_df2_nd, (tTestEngTarget)&api_bqriirf_df2_nd, 0, "iirbqnd2/bqriirf_df2_nd_bpf3.seq" },
    { 2, &descr_bqriirf_df2_nd, (tTestEngTarget)&api_bqriirf_df2_nd, 0, "iirbqnd2/bqriirf_df2_nd_bsf1.seq" },
    { 2, &descr_bqriirf_df2_nd, (tTestEngTarget)&api_bqriirf_df2_nd, 0, "iirbqnd2/bqriirf_df2_nd_hpf1.seq" },

    { 2, &descr_bqriirf_df2t_nd, (tTestEngTarget)&api_bqriirf_df2t_nd, 1, "iirbqnd2/bqriirf_df2t_nd_lpf1.seq" },
    { 2, &descr_bqriirf_df2t_nd, (tTestEngTarget)&api_bqriirf_df2t_nd, 1, "iirbqnd2/bqriirf_df2t_nd_bpf1.seq" },
    { 2, &descr_bqriirf_df2t_nd, (tTestEngTarget)&api_bqriirf_df2t_nd, 0, "iirbqnd2/bqriirf_df2t_nd_bpf2.seq" },
    { 2, &descr_bqriirf_df2t_nd, (tTestEngTarget)&api_bqriirf_df2t_nd, 0, "iirbqnd2/bqriirf_df2t_nd_bpf3.seq" },
    { 2, &descr_bqriirf_df2t_nd, (tTestEngTarget)&api_bqriirf_df2t_nd, 0, "iirbqnd2/bqriirf_df2t_nd_bsf1.seq" },
    { 2, &descr_bqriirf_df2t_nd, (tTestEngTarget)&api_bqriirf_df2t_nd, 0, "iirbqnd2/bqriirf_df2t_nd_hpf1.seq" },

    { 2, &descr_bqciirf_df1_nd, (tTestEngTarget)&api_bqciirf_df1_nd, 1, "iirbqnd2/bqciirf_df1_nd_lpf1.seq" },
    { 2, &descr_bqciirf_df1_nd, (tTestEngTarget)&api_bqciirf_df1_nd, 1, "iirbqnd2/bqciirf_df1_nd_bpf1.seq" },
    { 2, &descr_bqciirf_df1_nd, (tTestEngTarget)&api_bqciirf_df1_nd, 0, "iirbqnd2/bqciirf_df1_nd_bpf2.seq" },
    { 2, &descr_bqciirf_df1_nd, (tTestEngTarget)&api_bqciirf_df1_nd, 0, "iirbqnd2/bqciirf_df1_nd_bsf1.seq" },
    { 2, &descr_bqciirf_df1_nd, (tTestEngTarget)&api_bqciirf_df1_nd, 0, "iirbqnd2/bqciirf_df1_nd_hpf1.seq" },

    { 2, &descr_stereo_bqriirf_df1_nd, (tTestEngTarget)&api_stereo_bqriirf_df1_nd, 1, "iirbqnd2/stereo_bqriirf_df1_nd_lpf1.seq" },
    { 2, &descr_stereo_bqriirf_df1_nd, (tTestEngTarget)&api_stereo_bqriirf_df1_nd, 1, "iirbqnd2/stereo_bqriirf_df1_nd_bpf1.seq" },
    { 2, &descr_stereo_bqriirf_df1_nd, (tTestEngTarget)&api_stereo_bqriirf_df1_nd, 0, "iirbqnd2/stereo_bqriirf_df1_nd_bpf2.seq" },
    { 2, &descr_stereo_bqriirf_df1_nd, (tTestEngTarget)&api_stereo_bqriirf_df1_nd, 0, "iirbqnd2/stereo_bqriirf_df1_nd_bpf3.seq" },
    { 2, &descr_stereo_bqriirf_df1_nd, (tTestEngTarget)&api_stereo_bqriirf_df1_nd, 0, "iirbqnd2/stereo_bqriirf_df1_nd_bsf1.seq" },
    { 2, &descr_stereo_bqriirf_df1_nd, (tTestEngTarget)&api_stereo_bqriirf_df1_nd, 0, "iirbqnd2/stereo_bqriirf_df1_nd_hpf1.seq" },
};

int func_iirbqnd2(int isFull, int isVerbose, int breakOnError)
{
    return main_iirbq(2, isFull, isVerbose, breakOnError, tests_nd, (int)(sizeof(tests_nd)/sizeof(tests_nd[0])));
}
