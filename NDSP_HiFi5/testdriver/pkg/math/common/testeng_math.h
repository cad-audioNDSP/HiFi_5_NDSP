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
* Test-engine add-on for vector mathematics
*/

#ifndef TESTENG_MATH_H__
#define TESTENG_MATH_H__

/* Cross-platform data type definitions. */
#include "types.h"
/* Test engine API. */
#include "testeng.h"

/* processing function for vec_recip16x16 */
void te_math_processFxn_recip16x16(tTestEngContext * context);
/* processing function for vec_divide16x16 */
void te_math_processFxn_divide16x16(tTestEngContext * context);
/* processing function for vec_divide32x32 */
void te_math_processFxn_divide32x32(tTestEngContext * context);
/* processing function for vec_recip32x32 vec_recip24x24 */
void te_math_processFxn_recip32x32(tTestEngContext * context);
/* load/processing function for vec_recip64x64/vec_divide64x64 */
int  te_math_loadFxn_recip64x64   (tTestEngContext * context);
void te_math_processFxn_vec_recip64x64(tTestEngContext * context);
void te_math_processFxn_scl_recip64x64(tTestEngContext * context);
int  te_math_loadFxn_divide64x64   (tTestEngContext * context);
void te_math_processFxn_vec_divide64x64(tTestEngContext * context);
void te_math_processFxn_scl_divide64x64(tTestEngContext * context);

/* Apply a function to the test case data set:
* scalar functions with single argument, e.g. cos() */
void te_math_processFxn_scl_vXvZ(tTestEngContext * context);
void te_math_processFxn_scl_vXvZ16(tTestEngContext * context);
void te_math_processFxn_scl_vXvZ32(tTestEngContext * context);
void te_math_processFxn_scl_vX32sY32vZ(tTestEngContext * context);
void te_math_processFxn_scl_vXsY32vZ32(tTestEngContext * context);
void te_math_processFxn_scl_vXcvZ(tTestEngContext * context);
/* Apply a function to the test case data set:
* scalar atan2(y,x) */
void te_math_processFxn_scl_atan2(tTestEngContext * context);
void te_math_processFxn_scl_dividef(tTestEngContext * context);
/* processing function for scl_recip16x16 */
void te_math_processFxn_scl_recip16x16(tTestEngContext * context);
/* processing function for scl_recip32x32 scl_recip24x24 */
void te_math_processFxn_scl_recip32x32(tTestEngContext * context);
/* processing function for scl_divide32x32 scl_divide24x24 */
void te_math_processFxn_scl_divide32x32(tTestEngContext * context);
/* processing function for scl_divide16x16 */
void te_math_processFxn_scl_divide16x16(tTestEngContext * context);
/* processing function for scl_divide64x32 */
void te_math_processFxn_scl_divide64x32(tTestEngContext * context);
void te_math_processFxn_pow32x32( tTestEngContext * context );
void te_math_processFxn_scl_vXvYvZ(tTestEngContext * context);
#endif

