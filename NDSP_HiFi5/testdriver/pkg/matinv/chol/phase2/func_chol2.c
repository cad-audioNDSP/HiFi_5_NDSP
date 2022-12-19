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

#include "../common/test_chol.h"
#include "packages.h"

static const struct 
{
  tTestEngTarget funcList[MAX_FUNC_NUM];
  tTestEngDesc   testDesc;
}
testDefTbl[] =
{
    // floating point
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp4x4f           ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol4x4x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp6x6f           ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol6x6x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp8x8f           ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol8x8x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholdecomp10x10f         ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholfwdsubst4x4f         ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol4x4x1nfApi,   TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholfwdsubst6x6f         ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol6x6x1nfApi,   TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholfwdsubst8x8f         ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol8x8x1nfApi,   TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholfwdsubst10x10f       ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol10x10x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholbkwsubst4x4f         ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol4x4x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholbkwsubst6x6f         ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol6x6x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholbkwsubst8x8f         ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol8x8x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholbkwsubst10x10f       ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholmmsesolver4x4f       ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol4x4x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&cmatcholmmsesolver6x6f       ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol6x6x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&cmatcholmmsesolver8x8f       ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol8x8x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&cmatcholmmsesolver10x10f     ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol10x10x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess4x4f       ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol4x4x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess6x6f       ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol6x6x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess8x8f       ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol8x8x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpreprocess10x10f     ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv4x4f        ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol4x4x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv6x6f        ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol6x6x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv8x8f        ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol8x8x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&cmatcholpseudoinv10x10f      ), TEST_DESC(FMT_CPLX | FMT_FLOAT32, 0, &chol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
	
    { FUNC_LIST((tTestEngTarget)&matcholdecomp4x4f            ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol4x4x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&matcholdecomp6x6f            ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol6x6x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&matcholdecomp8x8f            ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol8x8x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&matcholdecomp10x10f          ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholf,     &te_processFxn_cholf) },
    { FUNC_LIST((tTestEngTarget)&matcholfwdsubst4x4f          ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol4x4x1nfApi,   TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&matcholfwdsubst6x6f          ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol6x6x1nfApi,   TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&matcholfwdsubst8x8f          ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol8x8x1nfApi,   TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&matcholfwdsubst10x10f        ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol10x10x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholfwdf, &te_processFxn_cholfwdf) },
    { FUNC_LIST((tTestEngTarget)&matcholbkwsubst4x4f          ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol4x4x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&matcholbkwsubst6x6f          ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol6x6x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&matcholbkwsubst8x8f          ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol8x8x1nfApi,   TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&matcholbkwsubst10x10f        ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholbkwf, &te_processFxn_cholbkwf) },
    { FUNC_LIST((tTestEngTarget)&matcholmmsesolver4x4f        ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol4x4x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&matcholmmsesolver6x6f        ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol6x6x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&matcholmmsesolver8x8f        ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol8x8x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&matcholmmsesolver10x10f      ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol10x10x1nfApi, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_cholmmsef, &te_processFxn_cholmmsef) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess4x4f        ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol4x4x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess6x6f        ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol6x6x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess8x8f        ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol8x8x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpreprocess10x10f      ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpreprocessf, &te_processFxn_cholpreprocessf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv4x4f         ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol4x4x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv6x6f         ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol6x6x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv8x8f         ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol8x8x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },
    { FUNC_LIST((tTestEngTarget)&matcholpseudoinv10x10f       ), TEST_DESC(FMT_REAL | FMT_FLOAT32, 0, &rchol10x10x1nfApi, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_cholpinvf, &te_processFxn_cholpinvf) },

    { FUNC_LIST(NULL), TEST_DESC(0,0,0,0,0, NULL, NULL) } /* End of table */
};

int func_chol2(int isFull, int isVerbose, int breakOnError)
{
  int res = 1;
  
  #define DO_TEST(fxn, seqFile)                                                                  \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                               \
                                                       sizeof(testDefTbl)/sizeof(testDefTbl[0]), \
                                                       MAX_FUNC_NUM,                             \
                                                       (tTestEngTarget)(fxn), "chol2/" seqFile,   \
                                                       isFull, isVerbose, breakOnError ) )

  DO_TEST( &cmatcholdecomp4x4f      , "cmatcholdecomp4x4f_cond2.seq"          );
  DO_TEST( &cmatcholdecomp4x4f      , "cmatcholdecomp4x4f_cond10.seq"         );
  DO_TEST( &cmatcholdecomp6x6f      , "cmatcholdecomp6x6f_cond2.seq"          );
  DO_TEST( &cmatcholdecomp6x6f      , "cmatcholdecomp6x6f_cond10.seq"         );
  DO_TEST( &cmatcholdecomp8x8f      , "cmatcholdecomp8x8f_cond2.seq"          );
  DO_TEST( &cmatcholdecomp8x8f      , "cmatcholdecomp8x8f_cond10.seq"         );
  DO_TEST( &cmatcholdecomp10x10f    , "cmatcholdecomp10x10f_cond2.seq"        );
  DO_TEST( &cmatcholdecomp10x10f    , "cmatcholdecomp10x10f_cond10.seq"       );

  DO_TEST( &cmatcholfwdsubst4x4f    , "cmatcholfwdsubst4x4x1f_cond2.seq"      );
  DO_TEST( &cmatcholfwdsubst4x4f    , "cmatcholfwdsubst4x4x1f_cond10.seq"     );
  DO_TEST( &cmatcholfwdsubst6x6f    , "cmatcholfwdsubst6x6x1f_cond2.seq"      );
  DO_TEST( &cmatcholfwdsubst6x6f    , "cmatcholfwdsubst6x6x1f_cond10.seq"     );
  DO_TEST( &cmatcholfwdsubst8x8f    , "cmatcholfwdsubst8x8x1f_cond2.seq"      );
  DO_TEST( &cmatcholfwdsubst8x8f    , "cmatcholfwdsubst8x8x1f_cond10.seq"     );
  DO_TEST( &cmatcholfwdsubst10x10f  , "cmatcholfwdsubst10x10x1f_cond2.seq"    );
  DO_TEST( &cmatcholfwdsubst10x10f  , "cmatcholfwdsubst10x10x1f_cond10.seq"   );

  DO_TEST( &cmatcholbkwsubst4x4f    , "cmatcholbkwsubst4x1f_cond2.seq"        );
  DO_TEST( &cmatcholbkwsubst4x4f    , "cmatcholbkwsubst4x1f_cond10.seq"       );
  DO_TEST( &cmatcholbkwsubst6x6f    , "cmatcholbkwsubst6x1f_cond2.seq"        );
  DO_TEST( &cmatcholbkwsubst6x6f    , "cmatcholbkwsubst6x1f_cond10.seq"       );
  DO_TEST( &cmatcholbkwsubst8x8f    , "cmatcholbkwsubst8x1f_cond2.seq"        );
  DO_TEST( &cmatcholbkwsubst8x8f    , "cmatcholbkwsubst8x1f_cond10.seq"       );
  DO_TEST( &cmatcholbkwsubst10x10f  , "cmatcholbkwsubst10x1f_cond2.seq"       );
  DO_TEST( &cmatcholbkwsubst10x10f  , "cmatcholbkwsubst10x1f_cond10.seq"      );

  DO_TEST( &cmatcholmmsesolver4x4f  , "cmatcholmmsesolver4x4x1f_cond2.seq"    );
  DO_TEST( &cmatcholmmsesolver4x4f  , "cmatcholmmsesolver4x4x1f_cond10.seq"   );
  DO_TEST( &cmatcholmmsesolver6x6f  , "cmatcholmmsesolver6x6x1f_cond2.seq"    );
  DO_TEST( &cmatcholmmsesolver6x6f  , "cmatcholmmsesolver6x6x1f_cond10.seq"   );
  DO_TEST( &cmatcholmmsesolver8x8f  , "cmatcholmmsesolver8x8x1f_cond2.seq"    );
  DO_TEST( &cmatcholmmsesolver8x8f  , "cmatcholmmsesolver8x8x1f_cond10.seq"   );
  DO_TEST( &cmatcholmmsesolver10x10f, "cmatcholmmsesolver10x10x1f_cond2.seq"  );
  DO_TEST( &cmatcholmmsesolver10x10f, "cmatcholmmsesolver10x10x1f_cond10.seq" );

  DO_TEST( &cmatcholpreprocess4x4f  , "cmatcholpreprocess4x4f_cond2.seq"      );
  DO_TEST( &cmatcholpreprocess4x4f  , "cmatcholpreprocess4x4f_cond10.seq"     );
  DO_TEST( &cmatcholpreprocess6x6f  , "cmatcholpreprocess6x6f_cond2.seq"      );
  DO_TEST( &cmatcholpreprocess6x6f  , "cmatcholpreprocess6x6f_cond10.seq"     );
  DO_TEST( &cmatcholpreprocess8x8f  , "cmatcholpreprocess8x8f_cond2.seq"      );
  DO_TEST( &cmatcholpreprocess8x8f  , "cmatcholpreprocess8x8f_cond10.seq"     );
  DO_TEST( &cmatcholpreprocess10x10f, "cmatcholpreprocess10x10f_cond2.seq"    );
  DO_TEST( &cmatcholpreprocess10x10f, "cmatcholpreprocess10x10f_cond10.seq"   );

  DO_TEST( &cmatcholpseudoinv4x4f   , "cmatcholpseudoinv4x4f_cond2.seq"       );
  DO_TEST( &cmatcholpseudoinv4x4f   , "cmatcholpseudoinv4x4f_cond10.seq"      );
  DO_TEST( &cmatcholpseudoinv6x6f   , "cmatcholpseudoinv6x6f_cond2.seq"       );
  DO_TEST( &cmatcholpseudoinv6x6f   , "cmatcholpseudoinv6x6f_cond10.seq"      );
  DO_TEST( &cmatcholpseudoinv8x8f   , "cmatcholpseudoinv8x8f_cond2.seq"       );
  DO_TEST( &cmatcholpseudoinv8x8f   , "cmatcholpseudoinv8x8f_cond10.seq"      );
  DO_TEST( &cmatcholpseudoinv10x10f , "cmatcholpseudoinv10x10f_cond2.seq"     );
  DO_TEST( &cmatcholpseudoinv10x10f , "cmatcholpseudoinv10x10f_cond10.seq"    );

  DO_TEST( &matcholdecomp4x4f       , "matcholdecomp4x4f_cond2.seq"           );
  DO_TEST( &matcholdecomp4x4f       , "matcholdecomp4x4f_cond10.seq"          );
  DO_TEST( &matcholdecomp6x6f       , "matcholdecomp6x6f_cond2.seq"           );
  DO_TEST( &matcholdecomp6x6f       , "matcholdecomp6x6f_cond10.seq"          );
  DO_TEST( &matcholdecomp8x8f       , "matcholdecomp8x8f_cond2.seq"           );
  DO_TEST( &matcholdecomp8x8f       , "matcholdecomp8x8f_cond10.seq"          );
  DO_TEST( &matcholdecomp10x10f     , "matcholdecomp10x10f_cond2.seq"         );
  DO_TEST( &matcholdecomp10x10f     , "matcholdecomp10x10f_cond10.seq"        );

  DO_TEST( &matcholfwdsubst4x4f     , "matcholfwdsubst4x4x1f_cond2.seq"       );
  DO_TEST( &matcholfwdsubst4x4f     , "matcholfwdsubst4x4x1f_cond10.seq"      );
  DO_TEST( &matcholfwdsubst6x6f     , "matcholfwdsubst6x6x1f_cond2.seq"       );
  DO_TEST( &matcholfwdsubst6x6f     , "matcholfwdsubst6x6x1f_cond10.seq"      );
  DO_TEST( &matcholfwdsubst8x8f     , "matcholfwdsubst8x8x1f_cond2.seq"       );
  DO_TEST( &matcholfwdsubst8x8f     , "matcholfwdsubst8x8x1f_cond10.seq"      );
  DO_TEST( &matcholfwdsubst10x10f   , "matcholfwdsubst10x10x1f_cond2.seq"     );
  DO_TEST( &matcholfwdsubst10x10f   , "matcholfwdsubst10x10x1f_cond10.seq"    );

  DO_TEST( &matcholbkwsubst4x4f     , "matcholbkwsubst4x1f_cond2.seq"         );
  DO_TEST( &matcholbkwsubst4x4f     , "matcholbkwsubst4x1f_cond10.seq"        );
  DO_TEST( &matcholbkwsubst6x6f     , "matcholbkwsubst6x1f_cond2.seq"         );
  DO_TEST( &matcholbkwsubst6x6f     , "matcholbkwsubst6x1f_cond10.seq"        );
  DO_TEST( &matcholbkwsubst8x8f     , "matcholbkwsubst8x1f_cond2.seq"         );
  DO_TEST( &matcholbkwsubst8x8f     , "matcholbkwsubst8x1f_cond10.seq"        );
  DO_TEST( &matcholbkwsubst10x10f   , "matcholbkwsubst10x1f_cond2.seq"        );
  DO_TEST( &matcholbkwsubst10x10f   , "matcholbkwsubst10x1f_cond10.seq"       );

  DO_TEST( &matcholmmsesolver4x4f  , "matcholmmsesolver4x4x1f_cond2.seq"     );
  DO_TEST( &matcholmmsesolver4x4f  , "matcholmmsesolver4x4x1f_cond10.seq"    );
  DO_TEST( &matcholmmsesolver6x6f  , "matcholmmsesolver6x6x1f_cond2.seq"     );
  DO_TEST( &matcholmmsesolver6x6f  , "matcholmmsesolver6x6x1f_cond10.seq"    );
  DO_TEST( &matcholmmsesolver8x8f  , "matcholmmsesolver8x8x1f_cond2.seq"     );
  DO_TEST( &matcholmmsesolver8x8f  , "matcholmmsesolver8x8x1f_cond10.seq"    );
  DO_TEST( &matcholmmsesolver10x10f, "matcholmmsesolver10x10x1f_cond2.seq"   );
  DO_TEST( &matcholmmsesolver10x10f, "matcholmmsesolver10x10x1f_cond10.seq"  );

  DO_TEST( &matcholpreprocess4x4f  , "matcholpreprocess4x4f_cond2.seq"       );
  DO_TEST( &matcholpreprocess4x4f  , "matcholpreprocess4x4f_cond10.seq"      );
  DO_TEST( &matcholpreprocess6x6f  , "matcholpreprocess6x6f_cond2.seq"       );
  DO_TEST( &matcholpreprocess6x6f  , "matcholpreprocess6x6f_cond10.seq"      );
  DO_TEST( &matcholpreprocess8x8f  , "matcholpreprocess8x8f_cond2.seq"       );
  DO_TEST( &matcholpreprocess8x8f  , "matcholpreprocess8x8f_cond10.seq"      );
  DO_TEST( &matcholpreprocess10x10f, "matcholpreprocess10x10f_cond2.seq"     );
  DO_TEST( &matcholpreprocess10x10f, "matcholpreprocess10x10f_cond10.seq"    );

  DO_TEST( &matcholpseudoinv4x4f   , "matcholpseudoinv4x4f_cond2.seq"        );
  DO_TEST( &matcholpseudoinv4x4f   , "matcholpseudoinv4x4f_cond10.seq"       );
  DO_TEST( &matcholpseudoinv6x6f   , "matcholpseudoinv6x6f_cond2.seq"        );
  DO_TEST( &matcholpseudoinv6x6f   , "matcholpseudoinv6x6f_cond10.seq"       );
  DO_TEST( &matcholpseudoinv8x8f   , "matcholpseudoinv8x8f_cond2.seq"        );
  DO_TEST( &matcholpseudoinv8x8f   , "matcholpseudoinv8x8f_cond10.seq"       );
  DO_TEST( &matcholpseudoinv10x10f , "matcholpseudoinv10x10f_cond2.seq"      );
  DO_TEST( &matcholpseudoinv10x10f , "matcholpseudoinv10x10f_cond10.seq"     );

  return res;
}
