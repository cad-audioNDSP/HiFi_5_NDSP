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
 * Test procedures for matrix functions
 */

#include <string.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API. */
#include LIBRARY_HEADER(matop)
/* Test engine API. */
#include "../common/testeng_matop.h"

#define MAX_FUNC_NUM   10
/* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }
/* Initializer for a test description structure. */
#define TEST_DESC( fmt, extra, api, dimNum, align, loadFxn, procFxn ) { (fmt),extra,api,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }

static const tMatmulApi_MNP apimtx_mpyf             = { &mtx_mpyf_getScratchSize             };
static const tMatmulApi_MNP apimtx_mpyf_fast        = { &mtx_mpyf_fast_getScratchSize        };
static const tMatmulApi_MNP apimtx_mpytf            = { &mtx_mpytf_getScratchSize            };
static const tMatmulApi_MNP apimtx_mpytf_fast       = { &mtx_mpytf_fast_getScratchSize       };
static const tMatmulApi_MNP apicmtx_mpytf_fast      = { &cmtx_mpytf_fast_getScratchSize      };
static const tMatmulApi_N   apicmtx_lrmpyf_fast     = { &cmtx_lrmpyf_fast_getScratchSize     };

/* vec API test definitions. */
static const struct 
{
tTestEngTarget   funcList[MAX_FUNC_NUM];
tTestEngDesc     testDesc;
}
testDefTbl[] =
{
    { FUNC_LIST( (tTestEngTarget)&mtx_mpyf ),          TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,&apimtx_mpyf      ,TE_DIM_NUM_4, TE_ALIGN_NO , &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpytf ),         TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,&apimtx_mpytf     ,TE_DIM_NUM_4, TE_ALIGN_NO , &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpyf ),       TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,NULL              ,TE_DIM_NUM_4, TE_ALIGN_NO , &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpyf_fast),      TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,&apimtx_mpyf_fast ,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_mpytf_fast),     TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,&apimtx_mpytf_fast,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_mmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_vecmpyf_fast ),  TEST_DESC( FMT_REAL|FMT_FLOAT32, MTX_PLAIN,NULL              ,TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_mmlt, &te_processFxn_vecmmlt ) },
    { FUNC_LIST( (tTestEngTarget)&mtx_transposef     ),TEST_DESC( FMT_REAL|FMT_FLOAT32, 0,NULL      ,TE_DIM_NUM_3, TE_ALIGN_NO,  &te_loadFxn_vXvZ, &te_processFxn_mtx_transpose) },
    { FUNC_LIST( (tTestEngTarget)&mtx_transposef_fast),TEST_DESC( FMT_REAL|FMT_FLOAT32, 0,NULL      ,TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_vXvZ, &te_processFxn_mtx_transpose) },

    { FUNC_LIST( (tTestEngTarget)&cmtx_mpytf_fast    ), TEST_DESC( FMT_CPLX|FMT_FLOAT32, 0, &apicmtx_mpytf_fast , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_mmlt   , &te_processFxn_mmlt    ) },
    { FUNC_LIST( (tTestEngTarget)&cmtx_vecmpytf_fast ), TEST_DESC( FMT_CPLX|FMT_FLOAT32, 0, NULL                , TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vecmpyt, &te_processFxn_vecmpyt ) },
    { FUNC_LIST( (tTestEngTarget)&cmtx_lrmpyf_fast   ), TEST_DESC( FMT_CPLX|FMT_FLOAT32, 0, &apicmtx_lrmpyf_fast, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_lrmpy  , &te_processFxn_lrmpy   ) },
    
    { FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, 0, 0, NULL, NULL ) } /* End of table */
};

/* Perform all tests for mat API functions. */
int func_matop2(int isFull, int isVerbose, int breakOnError)
{
  int res = 1;

  #define DO_TEST(fxn, seqFile)                                                                  \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                               \
                                                       sizeof(testDefTbl)/sizeof(testDefTbl[0]), \
                                                       MAX_FUNC_NUM,                             \
                                                       (tTestEngTarget)(fxn), "mtx2/" seqFile,    \
                                                       isFull, isVerbose, breakOnError ) )

    DO_TEST( &mtx_mpyf           , "mtx_mpyf.seq"           );
    DO_TEST( &mtx_mpyf_fast      , "mtx_mpyf_fast.seq"      );
    DO_TEST( &mtx_mpytf          , "mtx_mpytf.seq"          );
    DO_TEST( &mtx_mpytf_fast     , "mtx_mpytf_fast.seq"     );
    DO_TEST( &mtx_vecmpyf        , "mtx_vecmpyf.seq"        );
    DO_TEST( &mtx_vecmpyf_fast   , "mtx_vecmpyf_fast.seq"   );
    DO_TEST( &mtx_transposef     , "mtx_transposef.seq"     );
    DO_TEST( &mtx_transposef_fast, "mtx_transposef_fast.seq");
    DO_TEST( &cmtx_mpytf_fast    , "cmtx_mpytf_fast.seq"    );
    DO_TEST( &cmtx_lrmpyf_fast   , "cmtx_lrmpyf_fast.seq"   );
    DO_TEST( &cmtx_vecmpytf_fast , "cmtx_vecmpytf_fast.seq" );

    return (res);
}
