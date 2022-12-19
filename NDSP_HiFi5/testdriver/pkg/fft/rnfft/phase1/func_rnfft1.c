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
* Test procedures for complex FFT functions.
*/
#include "../../rfft/common/test_rfft.h"

static TestDef_t testTbl_rnfft[] =
{
#if 1
    // Mixed radix rfft
    { 1, "rnfft1/fft160_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rnfft1/fft192_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rnfft1/fft240_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rnfft1/fft320_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rnfft1/fft384_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rnfft1/fft480_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },

    { 1, "rnfft1/fft160_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rnfft1/fft192_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rnfft1/fft240_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rnfft1/fft320_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rnfft1/fft384_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rnfft1/fft480_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },

    { 1, "rnfft1/ifft160_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rnfft1/ifft192_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rnfft1/ifft240_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rnfft1/ifft320_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rnfft1/ifft384_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rnfft1/ifft480_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },

    { 1, "rnfft1/ifft160_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rnfft1/ifft192_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rnfft1/ifft240_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rnfft1/ifft320_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rnfft1/ifft384_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rnfft1/ifft480_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
#endif
#if 1
    // Mixed radix rfft
    { 1, "rnfft1/fft160_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rnfft1/fft192_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rnfft1/fft240_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rnfft1/fft320_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rnfft1/fft384_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rnfft1/fft480_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },

    { 1, "rnfft1/fft160_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rnfft1/fft192_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rnfft1/fft240_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rnfft1/fft320_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rnfft1/fft384_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rnfft1/fft480_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },

    { 1, "rnfft1/ifft160_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rnfft1/ifft192_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rnfft1/ifft240_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rnfft1/ifft320_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rnfft1/ifft384_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rnfft1/ifft480_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },

    { 1, "rnfft1/ifft160_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rnfft1/ifft192_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rnfft1/ifft240_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rnfft1/ifft320_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rnfft1/ifft384_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rnfft1/ifft480_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
#endif
#if 1

    // Mixed radix rfft
    { 1, "rnfft1/fft12_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft24_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft30_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft36_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft48_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft60_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft72_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft90_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft96_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft108_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft120_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft144_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft160_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft180_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft192_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft216_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft240_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft288_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft300_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft320_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft324_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft360_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft384_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft432_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft480_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft540_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft576_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft720_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft768_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft960_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft1152_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft1440_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft1536_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft1920_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },

    { 1, "rnfft1/fft12_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft24_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft30_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft36_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft48_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft60_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft72_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft90_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft96_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft108_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft120_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft144_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft160_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft180_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft192_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft216_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft240_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft288_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft300_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft320_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft324_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft360_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft384_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft432_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft480_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft540_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft576_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft720_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft768_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft960_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft1152_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft1440_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft1536_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rnfft1/fft1920_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },

    { 1, "rnfft1/ifft12_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft24_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft30_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft36_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft48_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft60_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft72_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft90_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft96_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft108_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft120_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft144_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft160_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft180_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft192_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft216_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft240_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft288_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft300_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft320_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft324_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft360_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft384_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft432_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft480_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft540_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft576_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft720_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft768_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft960_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft1152_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft1440_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft1536_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft1920_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },

    { 1, "rnfft1/ifft12_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft24_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft30_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft36_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft48_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft60_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft72_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft90_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft96_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft108_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft120_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft144_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft160_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft180_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft192_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft216_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft240_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft288_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft300_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft320_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft324_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft360_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft384_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft432_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft480_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft540_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft576_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft720_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft768_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft960_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft1152_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft1440_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft1536_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rnfft1/ifft1920_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
#endif
    { 0 } /* End of table */
}; //testTbl_rnfft

/* Perform all tests for fft_realMxN, ifft_realMxN API functions. */
int main_rnfft( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
  int ix, res;
  for ( ix=0,res=1; testTbl_rnfft[ix].seqFilename && ( res || !breakOnError ); ix++ )
  {
    if ( phaseNum == 0 || phaseNum == testTbl_rnfft[ix].phaseNum )
    {
        tTestEngTarget target = testTbl_rnfft[ix].target;
        /* Make sure that all functions is present */
        if (!IS_PRESENT(testTbl_rnfft[ix].desc.frwTransFxn) &&
            !IS_PRESENT(testTbl_rnfft[ix].desc.invTransFxn))
        {
            target = (tTestEngTarget)testTbl_rnfft[ix].desc.frwTransFxn;
        }

        res &= ( 0 != TestEngRun(  target,
                                   &testTbl_rnfft[ix].desc.desc,
                                   testTbl_rnfft[ix].seqFilename,
                                   isFull, isVerbose, breakOnError,0 ) );
    }
  }
  return (res);
} /* main_rnfft() */

int func_rnfft1(int isFull, int isVerbose, int breakOnError)
{
    return main_rnfft(1, isFull, isVerbose, breakOnError);
}
