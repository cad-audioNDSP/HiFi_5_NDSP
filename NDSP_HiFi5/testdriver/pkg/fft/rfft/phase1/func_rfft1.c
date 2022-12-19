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
#include "../common/test_rfft.h"

static TestDef_t testTbl_rfft[] =
{
#if 1
    { 1, "rfft1/fft32_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft64_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft128_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft256_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft512_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft1024_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft2048_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft4096_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft8192_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },

    { 1, "rfft1/fft32_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft64_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft128_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft256_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft512_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft1024_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft2048_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft4096_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },
    { 1, "rfft1/fft8192_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real16x16, NULL) },

    { 1, "rfft1/ifft32_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft64_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft128_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft256_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft512_real16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft1024_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft2048_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft4096_real16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft8192_real16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },

    { 1, "rfft1/ifft32_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft64_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft128_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft256_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft512_real16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft1024_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft2048_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft4096_real16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
    { 1, "rfft1/ifft8192_real16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real16x16) },
#endif
#if 0 //HiFi3/3z API
    { 1, "rfft1/fft_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
    { 1, "rfft1/fft_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
    { 1, "rfft1/fft_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
    { 1, "rfft1/fft_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
    { 1, "rfft1/ifft_real24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
    { 1, "rfft1/ifft_real24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
    { 1, "rfft1/ifft_real24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
    { 1, "rfft1/ifft_real24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_UNPACKED24), TE_ALIGN_YES, fft_real24x24, ifft_real24x24) },
#endif 
#if 1
    { 1, "rfft1/fft32_real32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft64_real32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft128_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft256_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft512_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft1024_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft2048_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft4096_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft8192_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },

    { 1, "rfft1/fft32_real32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft64_real32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft128_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft256_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft512_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft1024_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft2048_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft4096_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },
    { 1, "rfft1/fft8192_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, fft_real32x16, NULL) },

    { 1, "rfft1/ifft32_real32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft64_real32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft128_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft256_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft512_real32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft1024_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft2048_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft4096_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft8192_real32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },

    { 1, "rfft1/ifft32_real32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft64_real32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft128_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft256_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft512_real32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft1024_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft2048_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft4096_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
    { 1, "rfft1/ifft8192_real32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_32X16), TE_ALIGN_YES, NULL, ifft_real32x16) },
#endif

#if 1
    { 1, "rfft1/fft32_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft64_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft128_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft256_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft512_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft1024_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft2048_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft4096_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft8192_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },

    { 1, "rfft1/fft32_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft64_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft128_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft256_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft512_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft1024_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft2048_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft4096_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },
    { 1, "rfft1/fft8192_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,    TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, fft_real32x32, NULL) },

    { 1, "rfft1/ifft32_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft64_real32x32s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft128_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft256_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft512_real32x32s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft1024_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft2048_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft4096_real32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },

    { 1, "rfft1/ifft32_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft64_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft128_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft256_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft512_real32x32s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft1024_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft2048_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft4096_real32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
    { 1, "rfft1/ifft8192_real32x32s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft,   TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH), TE_ALIGN_YES, NULL, ifft_real32x32) },
#endif 
    { 0 } /* End of table */
}; //testTbl_rfft

/* Perform all tests for fft_realMxN, ifft_realMxN API functions. */
int main_rfft( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
  int ix, res;
  for ( ix=0,res=1; testTbl_rfft[ix].seqFilename && ( res || !breakOnError ); ix++ )
  {
    if ( phaseNum == 0 || phaseNum == testTbl_rfft[ix].phaseNum )
    {
        tTestEngTarget target = testTbl_rfft[ix].target;
        /* Make sure that all functions is present */
        if (!IS_PRESENT(testTbl_rfft[ix].desc.frwTransFxn) &&
            !IS_PRESENT(testTbl_rfft[ix].desc.invTransFxn))
        {
            target = (tTestEngTarget)testTbl_rfft[ix].desc.frwTransFxn;
        }

        res &= ( 0 != TestEngRun(  target,
                                   &testTbl_rfft[ix].desc.desc,
                                   testTbl_rfft[ix].seqFilename,
                                   isFull, isVerbose, breakOnError,0 ) );
    }
  }
  return (res);
} /* main_rfft() */

int func_rfft1(int isFull, int isVerbose, int breakOnError)
{
    return main_rfft(1, isFull, isVerbose, breakOnError);
}
