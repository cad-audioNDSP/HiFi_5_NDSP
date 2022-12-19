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
 * Test procedures for Cholesky decomposition (block data)
 */

#include "packages.h"
#include "../common/test_chol.h"

static const struct 
{
  tTestEngTarget funcList[MAX_FUNC_NUM];
  tTestEngDesc   testDesc;
}
testDefTbl[] =
{
    // 32-bit fixed point
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp4x4_32x32       ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol4x4_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp6x6_32x32       ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol6x6_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp8x8_32x32       ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol8x8_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp10x10_32x32     ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol10x10_32x32_Api, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },
    { FUNC_LIST((tTestEngTarget)&matcholdecomp4x4_32x32        ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol4x4_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },
    { FUNC_LIST((tTestEngTarget)&matcholdecomp6x6_32x32        ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol6x6_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },
    { FUNC_LIST((tTestEngTarget)&matcholdecomp8x8_32x32        ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol8x8_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },
    { FUNC_LIST((tTestEngTarget)&matcholdecomp10x10_32x32      ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol10x10_32x32_Api, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_chol,     &te_processFxn_chol) },

    { FUNC_LIST( (tTestEngTarget)&cmatcholbkwsubst4x4_32x32    ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol4x4_32x32_Api  , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholbkwsubst6x6_32x32    ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol6x6_32x32_Api  , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholbkwsubst8x8_32x32    ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol8x8_32x32_Api  , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholbkwsubst10x10_32x32  ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol10x10_32x32_Api, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },
    { FUNC_LIST( (tTestEngTarget)&matcholbkwsubst4x4_32x32     ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol4x4_32x32_Api   , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },
    { FUNC_LIST( (tTestEngTarget)&matcholbkwsubst6x6_32x32     ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol6x6_32x32_Api   , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },
    { FUNC_LIST( (tTestEngTarget)&matcholbkwsubst8x8_32x32     ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol8x8_32x32_Api   , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },
    { FUNC_LIST( (tTestEngTarget)&matcholbkwsubst10x10_32x32   ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol10x10_32x32_Api , TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkw, &te_processFxn_cholbkw ) },

    { FUNC_LIST( (tTestEngTarget)&cmatcholfwdsubst4x4_32x32    ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol4x4_32x32_Api  , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholfwdsubst6x6_32x32    ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol6x6_32x32_Api  , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholfwdsubst8x8_32x32    ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol8x8_32x32_Api  , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholfwdsubst10x10_32x32  ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &cmatchol10x10_32x32_Api, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },
    { FUNC_LIST( (tTestEngTarget)&matcholfwdsubst4x4_32x32     ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol4x4_32x32_Api   , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },
    { FUNC_LIST( (tTestEngTarget)&matcholfwdsubst6x6_32x32     ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol6x6_32x32_Api   , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },
    { FUNC_LIST( (tTestEngTarget)&matcholfwdsubst8x8_32x32     ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol8x8_32x32_Api   , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },
    { FUNC_LIST( (tTestEngTarget)&matcholfwdsubst10x10_32x32   ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED, &matchol10x10_32x32_Api , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwd, &te_processFxn_cholfwd ) },

    { FUNC_LIST( (tTestEngTarget)&cmatcholmmsesolver4x4_32x32  ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol4x4_32x32_Api  , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholmmsesolver6x6_32x32  ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol6x6_32x32_Api  , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholmmsesolver8x8_32x32  ),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol8x8_32x32_Api  , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },
    { FUNC_LIST( (tTestEngTarget)&cmatcholmmsesolver10x10_32x32),TEST_DESC( FMT_CPLX|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol10x10_32x32_Api, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },
    { FUNC_LIST( (tTestEngTarget)&matcholmmsesolver4x4_32x32   ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol4x4_32x32_Api   , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },
    { FUNC_LIST( (tTestEngTarget)&matcholmmsesolver6x6_32x32   ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol6x6_32x32_Api   , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },
    { FUNC_LIST( (tTestEngTarget)&matcholmmsesolver8x8_32x32   ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol8x8_32x32_Api   , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },
    { FUNC_LIST( (tTestEngTarget)&matcholmmsesolver10x10_32x32 ),TEST_DESC( FMT_REAL|FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol10x10_32x32_Api , TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmse, &te_processFxn_cholmmse ) },

    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess4x4_32x32   ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol4x4_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess6x6_32x32   ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol6x6_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess8x8_32x32   ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol8x8_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess10x10_32x32 ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol10x10_32x32_Api, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess4x4_32x32    ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol4x4_32x32_Api,    TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess6x6_32x32    ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol6x6_32x32_Api,    TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess8x8_32x32    ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol8x8_32x32_Api,    TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess10x10_32x32  ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol10x10_32x32_Api,  TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf,     &te_processFxn_cholpreprocessf) },

    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv4x4_32x32    ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol4x4_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv6x6_32x32    ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol6x6_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv8x8_32x32    ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol8x8_32x32_Api,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv10x10_32x32  ), TEST_DESC(FMT_CPLX | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &cmatchol10x10_32x32_Api, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv4x4_32x32     ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol4x4_32x32_Api,    TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv6x6_32x32     ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol6x6_32x32_Api,    TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv8x8_32x32     ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol8x8_32x32_Api,    TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv10x10_32x32   ), TEST_DESC(FMT_REAL | FMT_FRACT32, CHOLAPI_FIX_SZ|CHOLAPI_PACKED|CHOLAPI_USESIGMA2, &matchol10x10_32x32_Api,  TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },

    { FUNC_LIST(NULL), TEST_DESC(0,0,0,0,0, NULL, NULL) } /* End of table */
};

int func_chol1(int isFull, int isVerbose, int breakOnError)
{
  int res = 1;

  #define DO_TEST(fxn, seqFile)                                                                  \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                               \
                                                       sizeof(testDefTbl)/sizeof(testDefTbl[0]), \
                                                       MAX_FUNC_NUM,                             \
                                                       (tTestEngTarget)(fxn), "chol1/" seqFile,   \
                                                       isFull, isVerbose, breakOnError ) )

  DO_TEST( &cmatcholdecomp4x4_32x32       , "cmatcholdecomp4x4_32x32_cond20dynRange29qRAm3_0.seq"     );
  DO_TEST( &cmatcholdecomp4x4_32x32       , "cmatcholdecomp4x4_32x32_cond200dynRange29qRAm3_0.seq"    );
  DO_TEST( &cmatcholdecomp4x4_32x32       , "cmatcholdecomp4x4_32x32_cond2000dynRange29qRAm3_0.seq"   );
  DO_TEST( &cmatcholdecomp6x6_32x32       , "cmatcholdecomp6x6_32x32_cond20dynRange29qRAm3_0.seq"     );
  DO_TEST( &cmatcholdecomp6x6_32x32       , "cmatcholdecomp6x6_32x32_cond200dynRange29qRAm3_0.seq"    );
  DO_TEST( &cmatcholdecomp6x6_32x32       , "cmatcholdecomp6x6_32x32_cond2000dynRange29qRAm3_0.seq"   );
  DO_TEST( &cmatcholdecomp8x8_32x32       , "cmatcholdecomp8x8_32x32_cond20dynRange29qRAm3_0.seq"     );
  DO_TEST( &cmatcholdecomp8x8_32x32       , "cmatcholdecomp8x8_32x32_cond200dynRange29qRAm3_0.seq"    );
  DO_TEST( &cmatcholdecomp8x8_32x32       , "cmatcholdecomp8x8_32x32_cond2000dynRange29qRAm3_0.seq"   );
  DO_TEST( &cmatcholdecomp10x10_32x32     , "cmatcholdecomp10x10_32x32_cond20dynRange29qRAm3_0.seq"   );
  DO_TEST( &cmatcholdecomp10x10_32x32     , "cmatcholdecomp10x10_32x32_cond200dynRange29qRAm3_0.seq"  );
  DO_TEST( &cmatcholdecomp10x10_32x32     , "cmatcholdecomp10x10_32x32_cond2000dynRange29qRAm3_0.seq" );
  DO_TEST( &cmatcholfwdsubst4x4_32x32     , "cmatcholfwdsubst4x4_32x32_cond20dynRange29qRAm3_0.seq"   );
  DO_TEST( &cmatcholfwdsubst6x6_32x32     , "cmatcholfwdsubst6x6_32x32_cond20dynRange29qRAm3_0.seq"   );
  DO_TEST( &cmatcholfwdsubst8x8_32x32     , "cmatcholfwdsubst8x8_32x32_cond20dynRange29qRAm3_0.seq"   );
  DO_TEST( &cmatcholfwdsubst10x10_32x32   , "cmatcholfwdsubst10x10_32x32_cond20dynRange29qRAm3_0.seq" );
  DO_TEST( &cmatcholbkwsubst4x4_32x32     , "cmatcholbkwsubst4x4_32x32_cond20dynRange29qRAm3_0.seq"   );
  DO_TEST( &cmatcholbkwsubst6x6_32x32     , "cmatcholbkwsubst6x6_32x32_cond20dynRange29qRAm3_0.seq"   );
  DO_TEST( &cmatcholbkwsubst8x8_32x32     , "cmatcholbkwsubst8x8_32x32_cond20dynRange29qRAm3_0.seq"   );
  DO_TEST( &cmatcholbkwsubst10x10_32x32   , "cmatcholbkwsubst10x10_32x32_cond20dynRange29qRAm3_0.seq" );
  DO_TEST( &cmatcholmmsesolver4x4_32x32   , "cmatcholmmsesolver4x4_32x32_cond2dynRange29qRAm3_0.seq"  );
  DO_TEST( &cmatcholmmsesolver4x4_32x32   , "cmatcholmmsesolver4x4_32x32_cond20dynRange29qRAm3_0.seq" );
  DO_TEST( &cmatcholmmsesolver6x6_32x32   , "cmatcholmmsesolver6x6_32x32_cond2dynRange29qRAm3_0.seq"  );
  DO_TEST( &cmatcholmmsesolver6x6_32x32   , "cmatcholmmsesolver6x6_32x32_cond20dynRange29qRAm3_0.seq" );
  DO_TEST( &cmatcholmmsesolver8x8_32x32   , "cmatcholmmsesolver8x8_32x32_cond2dynRange29qRAm3_0.seq"  );
  DO_TEST( &cmatcholmmsesolver8x8_32x32   , "cmatcholmmsesolver8x8_32x32_cond20dynRange29qRAm3_0.seq" );
  DO_TEST( &cmatcholmmsesolver10x10_32x32 , "cmatcholmmsesolver10x10_32x32_cond2dynRange29qRAm3_0.seq");
  DO_TEST( &cmatcholmmsesolver10x10_32x32 , "cmatcholmmsesolver10x10_32x32_cond20dynRange29qRAm3_0.seq");
  DO_TEST( &cmatcholpreprocess4x4_32x32   , "cmatcholpreprocess4x4_32x32_cond20dynRange29qRAm3_0.seq" );
  DO_TEST( &cmatcholpreprocess6x6_32x32   , "cmatcholpreprocess6x6_32x32_cond20dynRange29qRAm3_0.seq" );
  DO_TEST( &cmatcholpreprocess8x8_32x32   , "cmatcholpreprocess8x8_32x32_cond20dynRange29qRAm3_0.seq" );
  DO_TEST( &cmatcholpreprocess10x10_32x32 , "cmatcholpreprocess10x10_32x32_cond20dynRange29qRAm3_0.seq");
  DO_TEST( &cmatcholpseudoinv4x4_32x32    , "cmatcholpseudoinv4x4_32x32_cond2.seq"                    );
  DO_TEST( &cmatcholpseudoinv6x6_32x32    , "cmatcholpseudoinv6x6_32x32_cond2.seq"                    );
  DO_TEST( &cmatcholpseudoinv8x8_32x32    , "cmatcholpseudoinv8x8_32x32_cond2.seq"                    );
  DO_TEST( &cmatcholpseudoinv10x10_32x32  , "cmatcholpseudoinv10x10_32x32_cond2.seq"                  );
  DO_TEST( &cmatcholpseudoinv4x4_32x32    , "cmatcholpseudoinv4x4_32x32_cond20.seq"                   );
  DO_TEST( &cmatcholpseudoinv6x6_32x32    , "cmatcholpseudoinv6x6_32x32_cond20.seq"                   );
  DO_TEST( &cmatcholpseudoinv8x8_32x32    , "cmatcholpseudoinv8x8_32x32_cond20.seq"                   );
  DO_TEST( &cmatcholpseudoinv10x10_32x32  , "cmatcholpseudoinv10x10_32x32_cond20.seq"                 );
  
  DO_TEST( &matcholdecomp4x4_32x32        , "matcholdecomp4x4_32x32_cond20dynRange29qRAm3_0.seq"      );
  DO_TEST( &matcholdecomp4x4_32x32        , "matcholdecomp4x4_32x32_cond200dynRange29qRAm3_0.seq"     );
  DO_TEST( &matcholdecomp4x4_32x32        , "matcholdecomp4x4_32x32_cond2000dynRange29qRAm3_0.seq"    );
  DO_TEST( &matcholdecomp6x6_32x32        , "matcholdecomp6x6_32x32_cond20dynRange29qRAm3_0.seq"      );
  DO_TEST( &matcholdecomp6x6_32x32        , "matcholdecomp6x6_32x32_cond200dynRange29qRAm3_0.seq"     );
  DO_TEST( &matcholdecomp6x6_32x32        , "matcholdecomp6x6_32x32_cond2000dynRange29qRAm3_0.seq"    );
  DO_TEST( &matcholdecomp8x8_32x32        , "matcholdecomp8x8_32x32_cond20dynRange29qRAm3_0.seq"      );
  DO_TEST( &matcholdecomp8x8_32x32        , "matcholdecomp8x8_32x32_cond200dynRange29qRAm3_0.seq"     );
  DO_TEST( &matcholdecomp8x8_32x32        , "matcholdecomp8x8_32x32_cond2000dynRange29qRAm3_0.seq"    );
  DO_TEST( &matcholdecomp10x10_32x32      , "matcholdecomp10x10_32x32_cond20dynRange29qRAm3_0.seq"    );
  DO_TEST( &matcholdecomp10x10_32x32      , "matcholdecomp10x10_32x32_cond200dynRange29qRAm3_0.seq"   );
  DO_TEST( &matcholdecomp10x10_32x32      , "matcholdecomp10x10_32x32_cond2000dynRange29qRAm3_0.seq"  );
  DO_TEST( &matcholfwdsubst4x4_32x32      , "matcholfwdsubst4x4_32x32_cond20dynRange29qRAm3_0.seq"    );
  DO_TEST( &matcholfwdsubst6x6_32x32      , "matcholfwdsubst6x6_32x32_cond20dynRange29qRAm3_0.seq"    );
  DO_TEST( &matcholfwdsubst8x8_32x32      , "matcholfwdsubst8x8_32x32_cond20dynRange29qRAm3_0.seq"    );
  DO_TEST( &matcholfwdsubst10x10_32x32    , "matcholfwdsubst10x10_32x32_cond20dynRange29qRAm3_0.seq"  );
  DO_TEST( &matcholbkwsubst4x4_32x32      , "matcholbkwsubst4x4_32x32_cond20dynRange29qRAm3_0.seq"    );
  DO_TEST( &matcholbkwsubst6x6_32x32      , "matcholbkwsubst6x6_32x32_cond20dynRange29qRAm3_0.seq"    );
  DO_TEST( &matcholbkwsubst8x8_32x32      , "matcholbkwsubst8x8_32x32_cond20dynRange29qRAm3_0.seq"    );
  DO_TEST( &matcholbkwsubst10x10_32x32    , "matcholbkwsubst10x10_32x32_cond20dynRange29qRAm3_0.seq"  );
  DO_TEST( &matcholmmsesolver4x4_32x32    , "matcholmmsesolver4x4_32x32_cond2dynRange29qRAm3_0.seq"   );
  DO_TEST( &matcholmmsesolver4x4_32x32    , "matcholmmsesolver4x4_32x32_cond20dynRange29qRAm3_0.seq"  );
  DO_TEST( &matcholmmsesolver6x6_32x32    , "matcholmmsesolver6x6_32x32_cond2dynRange29qRAm3_0.seq"   );
  DO_TEST( &matcholmmsesolver6x6_32x32    , "matcholmmsesolver6x6_32x32_cond20dynRange29qRAm3_0.seq"  );
  DO_TEST( &matcholmmsesolver8x8_32x32    , "matcholmmsesolver8x8_32x32_cond2dynRange29qRAm3_0.seq"   );
  DO_TEST( &matcholmmsesolver8x8_32x32    , "matcholmmsesolver8x8_32x32_cond20dynRange29qRAm3_0.seq"  );
  DO_TEST( &matcholmmsesolver10x10_32x32  , "matcholmmsesolver10x10_32x32_cond2dynRange29qRAm3_0.seq" );
  DO_TEST( &matcholmmsesolver10x10_32x32  , "matcholmmsesolver10x10_32x32_cond20dynRange29qRAm3_0.seq");
  DO_TEST( &matcholpreprocess4x4_32x32    , "matcholpreprocess4x4_32x32_cond20dynRange29qRAm3_0.seq"  );
  DO_TEST( &matcholpreprocess6x6_32x32    , "matcholpreprocess6x6_32x32_cond20dynRange29qRAm3_0.seq"  );
  DO_TEST( &matcholpreprocess8x8_32x32    , "matcholpreprocess8x8_32x32_cond20dynRange29qRAm3_0.seq"  );
  DO_TEST( &matcholpreprocess10x10_32x32  , "matcholpreprocess10x10_32x32_cond20dynRange29qRAm3_0.seq");
  DO_TEST( &matcholpseudoinv4x4_32x32     , "matcholpseudoinv4x4_32x32_cond2.seq"                     );
  DO_TEST( &matcholpseudoinv6x6_32x32     , "matcholpseudoinv6x6_32x32_cond2.seq"                     );
  DO_TEST( &matcholpseudoinv8x8_32x32     , "matcholpseudoinv8x8_32x32_cond2.seq"                     );
  DO_TEST( &matcholpseudoinv10x10_32x32   , "matcholpseudoinv10x10_32x32_cond2.seq"                   );
  DO_TEST( &matcholpseudoinv4x4_32x32     , "matcholpseudoinv4x4_32x32_cond20.seq"                    );
  DO_TEST( &matcholpseudoinv6x6_32x32     , "matcholpseudoinv6x6_32x32_cond20.seq"                    );
  DO_TEST( &matcholpseudoinv8x8_32x32     , "matcholpseudoinv8x8_32x32_cond20.seq"                    );
  DO_TEST( &matcholpseudoinv10x10_32x32   , "matcholpseudoinv10x10_32x32_cond20.seq"                  );

  return res;
}
