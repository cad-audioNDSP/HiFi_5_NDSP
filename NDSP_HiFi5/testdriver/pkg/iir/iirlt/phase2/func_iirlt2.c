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
#include "../common/test_iirlt.h"

static const tIirLatDescr api_latrf    ={(tIirLatFxnAlloc*)latrf_alloc,    (tIirLatFxnInit*)latrf_init,    (tIirLatFxnProcess*)latrf_process    };
static const tTestEngDesc descr_latrf         = { FMT_REAL | FMT_FLOAT32, 0, NULL,TE_DIM_NUM_1, TE_ALIGN_NO, te_create_iir_lat, te_destroy_iir_lat, &te_loadFxn_iir_lat, &te_processFxn_iir_lat };

static const tTbl tests[] =
{
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 1, "iirlt2/latrf_lpf1.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 0, "iirlt2/latrf_lpf2.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 0, "iirlt2/latrf_lpf3.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 0, "iirlt2/latrf_lpf4.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 1, "iirlt2/latrf_lpf5.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 0, "iirlt2/latrf_lpf6.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 0, "iirlt2/latrf_lpf7.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 0, "iirlt2/latrf_lpf8.seq"           },
  { 2, &descr_latrf              , (tTestEngTarget)&api_latrf              , 1, "iirlt2/latrf_lpf9.seq"           },
};

int func_iirlt2(int isFull, int isVerbose, int breakOnError)
{
    return main_iirlt(2, isFull, isVerbose, breakOnError, tests, (int)(sizeof(tests)/sizeof(tests[0])));
}
