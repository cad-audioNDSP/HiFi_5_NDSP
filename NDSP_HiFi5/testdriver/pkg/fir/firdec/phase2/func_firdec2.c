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

static tFirOldDescr api_firdecf         ={{(tFirOldFxnAlloc*)firdecf_alloc,       (tFirOldFxnInit*)firdecf_init,       (tFirOldFxnProcess*)firdecf_process       }};
static const tTestEngDesc descr_firdecf = { FMT_REAL|FMT_FLOAT32, TE_FIR_DN|TE_FIR_OLDDECIMATOR             ,NULL, TE_DIM_NUM_1, TE_ALIGN_YES, te_create_fir_old, te_destroy_fir_old, &te_loadFxn_fir_old, &te_processFxn_fir_old };

static const tTbl tests[] =
{
  { 2, &descr_firdecf, (tTestEngTarget)&api_firdecf,1,"firdec2/firdecf_dn2x.seq"},
  { 2, &descr_firdecf, (tTestEngTarget)&api_firdecf,1,"firdec2/firdecf_dn3x.seq"},
  { 2, &descr_firdecf, (tTestEngTarget)&api_firdecf,1,"firdec2/firdecf_dn11x.seq"},
  { 2, &descr_firdecf, (tTestEngTarget)&api_firdecf,1,"firdec2/firdecf_dn6x.seq"},
  { 2, &descr_firdecf, (tTestEngTarget)&api_firdecf,0,"firdec2/firdecf_dn4x.seq"},
  { 2, &descr_firdecf, (tTestEngTarget)&api_firdecf,0,"firdec2/firdecf_dn5x.seq"},
  { 2, &descr_firdecf, (tTestEngTarget)&api_firdecf,0,"firdec2/firdecf_dn23x.seq"},
};

/* Perform all tests for FIR API functions. */
int func_firdec2(int isFull, int isVerbose, int breakOnError)
{
    return test_firdec(2, isFull, isVerbose, breakOnError, tests, (int)(sizeof(tests)/sizeof(tests[0])));
}
