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
 * Test-engine add-on for log mel filterbank APIs.
 */

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* DSP Library API. */
#include LIBRARY_HEADER(audio)
/* Test engine API. */
#include "testeng.h"

#ifndef __TESTENG_LOGMEL_H
#define __TESTENG_LOGMEL_H

typedef void * te_logmel_handle_t;

typedef size_t             te_logmel_alloc_fxn_t( const logmel_params_t * params );
typedef te_logmel_handle_t te_logmel_init_fxn_t( void * objmem, const logmel_params_t * params );
typedef void               te_logmel_process_fxn_t( te_logmel_handle_t handle, void * pScr, void * logFbe, const void * spectra, int scaleExp );
typedef size_t             te_logmel_getScratchSize_fxn_t( const logmel_params_t * params );

typedef struct te_logmel_api_tag {
    te_logmel_alloc_fxn_t          * alloc;
    te_logmel_init_fxn_t           * init;
    te_logmel_process_fxn_t        * process;
    te_logmel_getScratchSize_fxn_t * getScratchSize;
} te_logmel_api_t;

/* Prepare to run test cases. */
int te_createFxn_logmel( tTestEngContext * context );
/* Release resources after the last test case, */
int te_destroyFxn_logmel( tTestEngContext * context );
/* Allocate vectors and load the data set for a test case. */
int te_loadFxn_logmel( tTestEngContext * context );
/* Apply the target function to the test case data set. */
void te_processFxn_logmel( tTestEngContext * context );

#endif /* __TESTENG_LOGMEL_H */

