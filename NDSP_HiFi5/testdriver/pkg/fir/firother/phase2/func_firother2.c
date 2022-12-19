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
#include "../common/test_firother.h"

const te_gccphat_api_t gccphatf_api     = { &gccphatf_getScratchSize      };

/* API test definitions. */
static const struct 
{
  tTestEngTarget   funcList[MAX_FUNC_NUM];
  tTestEngDesc     testDesc;
}
testDefTbl[] =
{
    {  FUNC_LIST( (tTestEngTarget)&fir_convolf),        TEST_DESC_CONVOLVE ( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorrf),         TEST_DESC_AUTOCORR ( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorrf),         TEST_DESC_CROSSCORR( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_xcorrf),       TEST_DESC_CROSSCORR( FMT_CPLX|FMT_FLOAT32, 0, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_convolaf),       TEST_DESC_CONVOLVE ( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_acorraf),        TEST_DESC_AUTOCORR ( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_xcorraf),        TEST_DESC_CROSSCORR( FMT_REAL|FMT_FLOAT32, 0, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_xcorraf),      TEST_DESC_CROSSCORR( FMT_CPLX|FMT_FLOAT32, 0, TE_ALIGN_NO  ) },
    {  FUNC_LIST( (tTestEngTarget)&gccphatf),           TEST_DESC_GCCPHAT( FMT_REAL|FMT_FLOAT32, 0, &gccphatf_api, TE_ALIGN_YES ) },
    {  FUNC_LIST( (tTestEngTarget)&fir_blmsf),          { FMT_REAL|FMT_FLOAT32, 0,NULL,TE_DIM_NUM_2, TE_ALIGN_YES, &te_create_lms,&te_destroy_lms, &te_loadFxn_lms, &te_processFxn_lms } },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_blmsf),        { FMT_CPLX|FMT_FLOAT32, 0,NULL,TE_DIM_NUM_2, TE_ALIGN_YES, &te_create_lms,&te_destroy_lms, &te_loadFxn_lms, &te_processFxn_lms } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms16x16_convergence),      { FMT_REAL|FMT_FRACT16, 0                 ,(void *)fir_blms16x16  ,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms16x32_convergence),      { FMT_REAL|FMT_FRACT16, TE_FIR_OTHER_32X16,(void *)fir_blms16x32  ,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms32x32_convergence),      { FMT_REAL|FMT_FRACT32, 0,                 (void *)fir_blms32x32  ,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blms32x32ep_convergence),    { FMT_REAL|FMT_FRACT32, 0,                 (void *)fir_blms32x32ep,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_blms32x32_convergence),    { FMT_CPLX|FMT_FRACT32, 0,                 (void *)cxfir_blms32x32,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( (tTestEngTarget)&fir_blmsf_convergence),          { FMT_REAL|FMT_FLOAT32, 0,                 (void *)fir_blmsf      ,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( (tTestEngTarget)&cxfir_blmsf_convergence),        { FMT_CPLX|FMT_FLOAT32, 0,                 (void *)cxfir_blmsf    ,TE_DIM_NUM_3, TE_ALIGN_YES, te_create_convergence,NULL, &te_loadFxn_lmsconv, &te_processFxn_lmsconv } },
    {  FUNC_LIST( NULL ),              TEST_DESC(  0, 0, 0, 0, NULL, NULL ) } /* End of table */
};

static const tTestEngDesc descr_cppaf      = { FMT_CPLX | FMT_FLOAT32, TE_FIR_PP|TE_FIR_POLYPHASE           , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_cppsf      = { FMT_CPLX | FMT_FLOAT32, TE_FIR_PP|TE_FIR_POLYPHASE           , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_rppaf      = { FMT_REAL | FMT_FLOAT32, TE_FIR_PP|TE_FIR_POLYPHASE           , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_rppsf      = { FMT_REAL | FMT_FLOAT32, TE_FIR_PP|TE_FIR_POLYPHASE           , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };

static tFirOldDescr api_cppaf        ={{(tFirOldFxnAlloc*)cppafirf_alloc,        (tFirOldFxnInit*)cppafirf_init,        (tFirOldFxnProcess*)cppafirf_process        ,NULL,NULL }};
static tFirOldDescr api_cppsf        ={{(tFirOldFxnAlloc*)cppsfirf_alloc,        (tFirOldFxnInit*)cppsfirf_init,        (tFirOldFxnProcess*)cppsfirf_process        ,NULL,NULL }};
static tFirOldDescr api_rppaf        ={{(tFirOldFxnAlloc*)rppafirf_alloc,        (tFirOldFxnInit*)rppafirf_init,        (tFirOldFxnProcess*)rppafirf_process        ,NULL,NULL }};
static tFirOldDescr api_rppsf        ={{(tFirOldFxnAlloc*)rppsfirf_alloc,        (tFirOldFxnInit*)rppsfirf_init,        (tFirOldFxnProcess*)rppsfirf_process        ,NULL,NULL }};

static const tTbl tests[] =
{
  { 2, &descr_rppaf, (tTestEngTarget)&api_rppaf, 1,"firother2/rppaf_filter_M32N24.seq"  },
  { 2, &descr_rppaf, (tTestEngTarget)&api_rppaf, 1,"firother2/rppaf_filter_M64N22.seq"  },
  { 2, &descr_rppaf, (tTestEngTarget)&api_rppaf, 1,"firother2/rppaf_filter_M96N20.seq"  },
  { 2, &descr_rppaf, (tTestEngTarget)&api_rppaf, 0,"firother2/rppaf_filter_M128N18.seq" },
  { 2, &descr_rppaf, (tTestEngTarget)&api_rppaf, 0,"firother2/rppaf_filter_M160N16.seq" },
  { 2, &descr_rppaf, (tTestEngTarget)&api_rppaf, 0,"firother2/rppaf_filter_M192N14.seq" },
  { 2, &descr_rppaf, (tTestEngTarget)&api_rppaf, 0,"firother2/rppaf_filter_M224N12.seq" },
  { 2, &descr_rppaf, (tTestEngTarget)&api_rppaf, 0,"firother2/rppaf_filter_M256N10.seq" },
  { 2, &descr_rppaf, (tTestEngTarget)&api_rppaf, 0,"firother2/rppaf_filter_M384N8.seq"  },
  { 2, &descr_rppaf, (tTestEngTarget)&api_rppaf, 0,"firother2/rppaf_filter_M512N6.seq"  },
  { 2, &descr_rppaf, (tTestEngTarget)&api_rppaf, 0,"firother2/rppaf_filter_M640N4.seq"  },

  { 2, &descr_rppsf, (tTestEngTarget)&api_rppsf, 1,"firother2/rppsf_filter_M32N24.seq"  },
  { 2, &descr_rppsf, (tTestEngTarget)&api_rppsf, 1,"firother2/rppsf_filter_M64N22.seq"  },
  { 2, &descr_rppsf, (tTestEngTarget)&api_rppsf, 1,"firother2/rppsf_filter_M96N20.seq"  },
  { 2, &descr_rppsf, (tTestEngTarget)&api_rppsf, 0,"firother2/rppsf_filter_M128N18.seq" },
  { 2, &descr_rppsf, (tTestEngTarget)&api_rppsf, 0,"firother2/rppsf_filter_M160N16.seq" },
  { 2, &descr_rppsf, (tTestEngTarget)&api_rppsf, 0,"firother2/rppsf_filter_M192N14.seq" },
  { 2, &descr_rppsf, (tTestEngTarget)&api_rppsf, 0,"firother2/rppsf_filter_M224N12.seq" },
  { 2, &descr_rppsf, (tTestEngTarget)&api_rppsf, 0,"firother2/rppsf_filter_M256N10.seq" },
  { 2, &descr_rppsf, (tTestEngTarget)&api_rppsf, 0,"firother2/rppsf_filter_M384N8.seq"  },
  { 2, &descr_rppsf, (tTestEngTarget)&api_rppsf, 0,"firother2/rppsf_filter_M512N6.seq"  },
  { 2, &descr_rppsf, (tTestEngTarget)&api_rppsf, 0,"firother2/rppsf_filter_M640N4.seq"  },

  { 2, &descr_cppaf, (tTestEngTarget)&api_cppaf, 1,"firother2/cppaf_filter_M32N24.seq"  },
  { 2, &descr_cppaf, (tTestEngTarget)&api_cppaf, 1,"firother2/cppaf_filter_M64N22.seq"  },
  { 2, &descr_cppaf, (tTestEngTarget)&api_cppaf, 1,"firother2/cppaf_filter_M96N20.seq"  },
  { 2, &descr_cppaf, (tTestEngTarget)&api_cppaf, 0,"firother2/cppaf_filter_M128N18.seq" },
  { 2, &descr_cppaf, (tTestEngTarget)&api_cppaf, 0,"firother2/cppaf_filter_M160N16.seq" },
  { 2, &descr_cppaf, (tTestEngTarget)&api_cppaf, 0,"firother2/cppaf_filter_M192N14.seq" },
  { 2, &descr_cppaf, (tTestEngTarget)&api_cppaf, 0,"firother2/cppaf_filter_M224N12.seq" },
  { 2, &descr_cppaf, (tTestEngTarget)&api_cppaf, 0,"firother2/cppaf_filter_M256N10.seq" },
  { 2, &descr_cppaf, (tTestEngTarget)&api_cppaf, 0,"firother2/cppaf_filter_M384N8.seq"  },
  { 2, &descr_cppaf, (tTestEngTarget)&api_cppaf, 0,"firother2/cppaf_filter_M512N6.seq"  },
  { 2, &descr_cppaf, (tTestEngTarget)&api_cppaf, 0,"firother2/cppaf_filter_M640N4.seq"  },

  { 2, &descr_cppsf, (tTestEngTarget)&api_cppsf, 1,"firother2/cppsf_filter_M32N24.seq"  },
  { 2, &descr_cppsf, (tTestEngTarget)&api_cppsf, 1,"firother2/cppsf_filter_M64N22.seq"  },
  { 2, &descr_cppsf, (tTestEngTarget)&api_cppsf, 1,"firother2/cppsf_filter_M96N20.seq"  },
  { 2, &descr_cppsf, (tTestEngTarget)&api_cppsf, 0,"firother2/cppsf_filter_M128N18.seq" },
  { 2, &descr_cppsf, (tTestEngTarget)&api_cppsf, 0,"firother2/cppsf_filter_M160N16.seq" },
  { 2, &descr_cppsf, (tTestEngTarget)&api_cppsf, 0,"firother2/cppsf_filter_M192N14.seq" },
  { 2, &descr_cppsf, (tTestEngTarget)&api_cppsf, 0,"firother2/cppsf_filter_M224N12.seq" },
  { 2, &descr_cppsf, (tTestEngTarget)&api_cppsf, 0,"firother2/cppsf_filter_M256N10.seq" },
  { 2, &descr_cppsf, (tTestEngTarget)&api_cppsf, 0,"firother2/cppsf_filter_M384N8.seq"  },
  { 2, &descr_cppsf, (tTestEngTarget)&api_cppsf, 0,"firother2/cppsf_filter_M512N6.seq"  },
  { 2, &descr_cppsf, (tTestEngTarget)&api_cppsf, 0,"firother2/cppsf_filter_M640N4.seq"  },
};

/* Perform all tests for FIR API functions. */
int func_firother2(int isFull, int isVerbose, int breakOnError)
{
    int res = 1;

    #define DO_TEST(fxn, seqFile)                                                                    \
        if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                               \
                                                           sizeof(testDefTbl)/sizeof(testDefTbl[0]), \
                                                           MAX_FUNC_NUM,                             \
                                                           (tTestEngTarget)(fxn), "firother2/" seqFile,    \
                                                           isFull, isVerbose, breakOnError ) )
    
    res &= test_firother(2, isFull, isVerbose, breakOnError, tests, (int)(sizeof(tests)/sizeof(tests[0])));

    DO_TEST( &fir_acorrf       , "fir_acorrf.seq"       );
    DO_TEST( &fir_xcorrf       , "fir_xcorrf.seq"       );
    DO_TEST( &cxfir_xcorrf     , "cxfir_xcorrf.seq"     );
    DO_TEST( &fir_convolf      , "fir_convolf.seq"      );
    DO_TEST( &fir_acorraf      , "fir_acorraf.seq"      );
    DO_TEST( &fir_xcorraf      , "fir_xcorraf.seq"      );
    DO_TEST( &cxfir_xcorraf    , "cxfir_xcorraf.seq"    );
    DO_TEST( &fir_convolaf     , "fir_convolaf.seq"     );
    DO_TEST( &gccphatf         , "gccphatf.seq"         );
    DO_TEST( &fir_blmsf        , "fir_blmsf.seq"        );
    DO_TEST( &cxfir_blmsf      , "cxfir_blmsf.seq"      );
    DO_TEST( &cxfir_blmsf_convergence      , "cxfir_blmsf_convergence.seq"      );
    DO_TEST( &fir_blmsf_convergence        , "fir_blmsf_convergence.seq"        );

    return (res);
}
