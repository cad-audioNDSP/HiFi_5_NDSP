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
#ifndef BKFIRF_H__
#define BKFIRF_H__

/* Portable data types. */
#include "NatureDSP_types.h"

#define BKFIRF_MAGIC     0xc47d93fa

/* Filter instance structure. */
typedef struct tag_bkfirf_t
{
    uint32_t          magic;     // Instance pointer validation number
    int               M;         // Number of filter coefficients
    const float32_t * coef;      // M filter coefficients
    float32_t       * delayLine; // Delay line for samples
    int               delayLen;  // Delay line length, in samples
    float32_t       * delayPos;  // Delay line slot to be filled next
} bkfirf_t, *bkfirf_ptr_t;

#endif /* BKFIRF_H__ */
