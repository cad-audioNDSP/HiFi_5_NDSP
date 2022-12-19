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
 * Test procedures for FIR filters
 */

#include <string.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* Library API */
#include LIBRARY_HEADER(fir)
/* Test Engine expansion for FIR filters. */
#include "../common/test_firdec.h"

/* ------------------------------------------------------------------------ */
/* Copyright (c) 2020 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
/*          Copyright (C) 2009-2020 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
 * Test procedures for FIR
 */
#include "../common/test_firdec.h"

static tFirOldDescr api_firdec16x16     ={{(tFirOldFxnAlloc*)firdec16x16_alloc,     (tFirOldFxnInit*)firdec16x16_init,     (tFirOldFxnProcess*)firdec16x16_process     }};
static tFirOldDescr api_firdec32x16     ={{(tFirOldFxnAlloc*)firdec32x16_alloc,     (tFirOldFxnInit*)firdec32x16_init,     (tFirOldFxnProcess*)firdec32x16_process     }};
static tFirOldDescr api_firdec32x32     ={{(tFirOldFxnAlloc*)firdec32x32_alloc,     (tFirOldFxnInit*)firdec32x32_init,     (tFirOldFxnProcess*)firdec32x32_process     }};
static tFirOldDescr api_firdec32x32ep   ={{(tFirOldFxnAlloc*)firdec32x32ep_alloc,   (tFirOldFxnInit*)firdec32x32ep_init,   (tFirOldFxnProcess*)firdec32x32ep_process   }};

static const tTestEngDesc descr_firdec16x16     = { FMT_REAL|FMT_FRACT16, TE_FIR_DN|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_firdec32x16     = { FMT_REAL|FMT_FRACT32, TE_FIR_DN|TE_FIR_OLDDECIMATOR|TE_FIR_FILTER_32X16, NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_firdec32x32     = { FMT_REAL|FMT_FRACT32, TE_FIR_DN|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };
static const tTestEngDesc descr_firdec32x32ep   = { FMT_REAL|FMT_FRACT32, TE_FIR_DN|TE_FIR_OLDDECIMATOR                    , NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };

static const tTbl tests[] =
{
  { 1, &descr_firdec16x16, (tTestEngTarget)&api_firdec16x16,1,"firdec1/firdec16x16_dn2x.seq" },
  { 1, &descr_firdec16x16, (tTestEngTarget)&api_firdec16x16,1,"firdec1/firdec16x16_dn3x.seq" },
  { 1, &descr_firdec16x16, (tTestEngTarget)&api_firdec16x16,1,"firdec1/firdec16x16_dn11x.seq"},
  { 1, &descr_firdec16x16, (tTestEngTarget)&api_firdec16x16,1,"firdec1/firdec16x16_dn6x.seq" },
  { 1, &descr_firdec16x16, (tTestEngTarget)&api_firdec16x16,0,"firdec1/firdec16x16_dn4x.seq" },
  { 1, &descr_firdec16x16, (tTestEngTarget)&api_firdec16x16,0,"firdec1/firdec16x16_dn5x.seq" },
  { 1, &descr_firdec16x16, (tTestEngTarget)&api_firdec16x16,0,"firdec1/firdec16x16_dn23x.seq"},

  { 1, &descr_firdec32x16, (tTestEngTarget)&api_firdec32x16,1,"firdec1/firdec32x16_dn2x.seq" },
  { 1, &descr_firdec32x16, (tTestEngTarget)&api_firdec32x16,1,"firdec1/firdec32x16_dn3x.seq" },
  { 1, &descr_firdec32x16, (tTestEngTarget)&api_firdec32x16,1,"firdec1/firdec32x16_dn11x.seq"},
  { 1, &descr_firdec32x16, (tTestEngTarget)&api_firdec32x16,1,"firdec1/firdec32x16_dn6x.seq" },
  { 1, &descr_firdec32x16, (tTestEngTarget)&api_firdec32x16,0,"firdec1/firdec32x16_dn4x.seq" },
  { 1, &descr_firdec32x16, (tTestEngTarget)&api_firdec32x16,0,"firdec1/firdec32x16_dn5x.seq" },
  { 1, &descr_firdec32x16, (tTestEngTarget)&api_firdec32x16,0,"firdec1/firdec32x16_dn23x.seq"},
#if 0// for HiFi3/3z
  { 1, &descr_firdec24x24, (tTestEngTarget)&api_firdec24x24,1,"firdec24x24_dn2x.seq" },
  { 1, &descr_firdec24x24, (tTestEngTarget)&api_firdec24x24,1,"firdec24x24_dn3x.seq" },
  { 1, &descr_firdec24x24, (tTestEngTarget)&api_firdec24x24,1,"firdec24x24_dn11x.seq"},
  { 1, &descr_firdec24x24, (tTestEngTarget)&api_firdec24x24,1,"firdec24x24_dn6x.seq" },
  { 1, &descr_firdec24x24, (tTestEngTarget)&api_firdec24x24,0,"firdec24x24_dn4x.seq" },
  { 1, &descr_firdec24x24, (tTestEngTarget)&api_firdec24x24,0,"firdec24x24_dn5x.seq" },
  { 1, &descr_firdec24x24, (tTestEngTarget)&api_firdec24x24,0,"firdec24x24_dn23x.seq"},
#endif
  { 1, &descr_firdec32x32, (tTestEngTarget)&api_firdec32x32,1,"firdec1/firdec32x32_dn2x.seq" },
  { 1, &descr_firdec32x32, (tTestEngTarget)&api_firdec32x32,1,"firdec1/firdec32x32_dn3x.seq" },
  { 1, &descr_firdec32x32, (tTestEngTarget)&api_firdec32x32,1,"firdec1/firdec32x32_dn11x.seq"},
  { 1, &descr_firdec32x32, (tTestEngTarget)&api_firdec32x32,1,"firdec1/firdec32x32_dn6x.seq" },
  { 1, &descr_firdec32x32, (tTestEngTarget)&api_firdec32x32,0,"firdec1/firdec32x32_dn4x.seq" },
  { 1, &descr_firdec32x32, (tTestEngTarget)&api_firdec32x32,0,"firdec1/firdec32x32_dn5x.seq" },
  { 1, &descr_firdec32x32, (tTestEngTarget)&api_firdec32x32,0,"firdec1/firdec32x32_dn23x.seq"},

  { 1, &descr_firdec32x32ep, (tTestEngTarget)&api_firdec32x32ep,1,"firdec1/firdec32x32ep_dn2x.seq" },
  { 1, &descr_firdec32x32ep, (tTestEngTarget)&api_firdec32x32ep,1,"firdec1/firdec32x32ep_dn3x.seq" },
  { 1, &descr_firdec32x32ep, (tTestEngTarget)&api_firdec32x32ep,1,"firdec1/firdec32x32ep_dn11x.seq"},
  { 1, &descr_firdec32x32ep, (tTestEngTarget)&api_firdec32x32ep,1,"firdec1/firdec32x32ep_dn6x.seq" },
  { 1, &descr_firdec32x32ep, (tTestEngTarget)&api_firdec32x32ep,0,"firdec1/firdec32x32ep_dn4x.seq" },
  { 1, &descr_firdec32x32ep, (tTestEngTarget)&api_firdec32x32ep,0,"firdec1/firdec32x32ep_dn5x.seq" },
  { 1, &descr_firdec32x32ep, (tTestEngTarget)&api_firdec32x32ep,0,"firdec1/firdec32x32ep_dn23x.seq"},
};

/* Perform all tests for FIR API functions. */
int func_firdec1(int isFull, int isVerbose, int breakOnError)
{
    return test_firdec(1, isFull, isVerbose, breakOnError, tests, (int)(sizeof(tests)/sizeof(tests[0])));
}
