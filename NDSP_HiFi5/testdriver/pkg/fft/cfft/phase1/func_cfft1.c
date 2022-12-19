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
#include "../common/test_cfft.h"

TestDef_t testTbl_cfft[] =
{
#if 1
    { 1, "cfft1/fft16_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft32_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft64_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft128_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft256_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft512_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft1024_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft2048_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft4096_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
  
    { 1, "cfft1/fft16_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft32_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft64_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft128_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft256_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft512_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft1024_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft2048_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/fft4096_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },

    { 1, "cfft1/ifft16_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft32_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft64_cplx16x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft128_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft256_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft512_cplx16x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft1024_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft2048_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft4096_cplx16x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },

    { 1, "cfft1/ifft16_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft32_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft64_cplx16x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft128_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft256_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft512_cplx16x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft1024_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft2048_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
    { 1, "cfft1/ifft4096_cplx16x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT16, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX), TE_ALIGN_YES, &fft_cplx16x16, &ifft_cplx16x16) },
#endif
#if 0 //HiFi3/3z API
    { 1, "cfft1/fft_cplx24x24s0.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
    { 1, "cfft1/fft_cplx24x24s1.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
    { 1, "cfft1/fft_cplx24x24s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
    { 1, "cfft1/fft_cplx24x24s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
    { 1, "cfft1/ifft_cplx24x24s0.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
    { 1, "cfft1/ifft_cplx24x24s1.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
    { 1, "cfft1/ifft_cplx24x24s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
    { 1, "cfft1/ifft_cplx24x24s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_UNPACKED24), TE_ALIGN_YES, &fft_cplx24x24, &ifft_cplx24x24) },
#endif

#if 1
    { 1, "cfft1/fft16_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft32_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft64_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft128_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft256_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft512_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft1024_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft2048_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft4096_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft16_cplx32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft32_cplx32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft64_cplx32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft128_cplx32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft256_cplx32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft512_cplx32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft1024_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft2048_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },
    { 1, "cfft1/fft4096_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) },

    { 1, "cfft1/ifft16_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft32_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft64_cplx32x16s3.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft128_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft256_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft512_cplx32x16s3.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft1024_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft2048_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft4096_cplx32x16s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft16_cplx32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft32_cplx32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft64_cplx32x16s2.seq"  , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft128_cplx32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft256_cplx32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft512_cplx32x16s2.seq" , (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft1024_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft2048_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
    { 1, "cfft1/ifft4096_cplx32x16s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, (ATTR_FAST | TE_FFT_OPT_SCALE_METH | TE_FFT_CPLX | TE_FFT_32X16), TE_ALIGN_YES, &fft_cplx32x16, &ifft_cplx32x16) }, 
#endif

#if 1
    { 1, "cfft1/fft16_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft32_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft64_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft128_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft256_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft512_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft1024_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft2048_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft4096_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },

    { 1, "cfft1/fft16_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft32_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft64_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft128_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft256_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft512_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft1024_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft2048_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },
    { 1, "cfft1/fft4096_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, &fft_cplx32x32, NULL) },

    { 1, "cfft1/ifft16_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft32_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft64_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft128_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft256_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft512_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft1024_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft2048_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft4096_cplx32x32s3.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },

    { 1, "cfft1/ifft16_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft32_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft64_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft128_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft256_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft512_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft1024_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft2048_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
    { 1, "cfft1/ifft4096_cplx32x32s2.seq", (tTestEngTarget)te_frameProc_stnd_scl_fft, TEST_DESC_STND(FMT_FRACT32, TE_FFT_CPLX | TE_FFT_OPT_REUSE_INBUF | TE_FFT_OPT_SCALE_METH, TE_ALIGN_YES, NULL, &ifft_cplx32x32) },
#endif
    { 0 } /* End of table */
}; // testTbl_cfft

/* Perform all tests for fft_cplx, ifft_cplx API functions. */
int main_cfft(int phaseNum, int isFull, int isVerbose, int breakOnError)
{
    int ix, res;
    for (ix = 0, res = 1; testTbl_cfft[ix].seqFilename && (res || !breakOnError); ix++)
    {
        if (phaseNum == 0 || phaseNum == testTbl_cfft[ix].phaseNum)
        {
            tTestEngTarget target = testTbl_cfft[ix].target;

            if (!IS_PRESENT(testTbl_cfft[ix].desc.frwTransFxn) &&
                !IS_PRESENT(testTbl_cfft[ix].desc.invTransFxn))
            {
                target = (tTestEngTarget)testTbl_cfft[ix].desc.frwTransFxn;
            }

            res &= (0 != TestEngRun(target,
                &testTbl_cfft[ix].desc.desc,
                testTbl_cfft[ix].seqFilename,
                isFull, isVerbose, breakOnError,0));
        }
    }
    return (res);

} /* main_cfft() */

int func_cfft1(int isFull, int isVerbose, int breakOnError)
{
    return main_cfft(1, isFull, isVerbose, breakOnError);
}
