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
#ifndef STEREO_BKFIRF_COMMON_H__
#define STEREO_BKFIRF_COMMON_H__

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"

#define STEREO_BKFIRF_MAGIC    0x830a4f21

/* Filter instance structure. */
typedef struct tag_stereo_bkfirf_t
{
#if !(HAVE_VFPU) //For scalar code
	uint32_t magic; // Should be STEREO_BKFIRF_MAGIC
	bkfirf_handle_t bkfir_left;
	bkfirf_handle_t bkfir_right;
	void * bkfir_left_mem; // Memory allocated for filters
	void * bkfir_right_mem;
#else
    uint32_t          magic;          // Instance pointer validation number
    int               M;              // Number of filter coefficients
    const float32_t * coefLeft;       // M filter coefficients
    const float32_t * coefRight;      // M filter coefficients
    float32_t       * delayLine;  // Delay line for samples
    int               delayLen;       // Delay line length, in samples
    float32_t       * delayPos;   // Delay line slot to be filled next
#endif
} stereo_bkfirf_t, *stereo_bkfirf_ptr_t;

#endif /* STEREO_BKFIRF_COMMON_H__ */
