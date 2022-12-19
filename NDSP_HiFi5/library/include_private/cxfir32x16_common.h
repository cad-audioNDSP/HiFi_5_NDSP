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
#ifndef CXFIR32X16_COMMON_H_
#define CXFIR32X16_COMMON_H_

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"

/* Instance pointer validation number. */
#define CXFIR32X16_MAGIC     0x34423c06

/* Filter instance structure. */
typedef struct tag_cxfir32x16_t
{
    uint32_t                magic;     // Instance pointer validation number
    int                     M;         // Number of complex filter coefficients
    const complex_fract16 * coef;      // M complex coefficients
    complex_fract32       * delayLine; // Delay line for complex samples
    int                     delayLen;  // Delay line length, in complex samples
    complex_fract32       * delayPos;  // Delay line slot to be filled next
} cxfir32x16_t, *cxfir32x16_ptr_t;

#endif
