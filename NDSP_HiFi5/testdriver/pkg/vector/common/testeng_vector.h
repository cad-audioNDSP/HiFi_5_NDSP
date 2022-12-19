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
* Test-engine add-on for arithmetic and logic functions on data vectors
*/

#ifndef TESTENG_VECTOR_H__
#define TESTENG_VECTOR_H__

/* Cross-platform data type definitions. */
#include "types.h"
/* Test engine API. */
#include "testeng.h"

/* Apply the target function to the test case data set:
* vector X (in), scalar F0 (in), scalar F1 (in), scalar Y (in), vector Z (out) */
void te_vector_processFxn_vXsF0sF1sYvZ(tTestEngContext * context);
/* Allocate vectors and load the data set:
* vector X (in), scalar F0 (in), scalar F1 (in), scalar Y (in), vector Z (out) */
int te_vector_loadFxn_vXsF0sF1sYvZ(tTestEngContext * context);

void processFxn_scl_vXvZ32( tTestEngContext * context );

#endif

