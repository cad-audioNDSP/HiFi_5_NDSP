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
* Test module for testing cycle performance (DCT)
*/
#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(fft)
#include "mips.h"

#define PROFILE_DCT4_32X16(isFull, isVerbose, N)  PROFILE_SIMPLE(isFull, isVerbose, dct4_32x16 , (mips.out0.i32, mips.inp0.i32, dct4_16_##N, 3), fout, "N="#N", scalingOpt=3"     , prf_cycle);
#define PROFILE_DCT4_32X32(isFull, isVerbose, N)  PROFILE_SIMPLE(isFull, isVerbose, dct4_32x32 , (mips.out0.i32, mips.inp0.i32, dct4_32_##N, 3), fout, "N="#N", scalingOpt=3"     , prf_cycle);
#define PROFILE_DCT4_24X24(isFull, isVerbose, N)  PROFILE_SIMPLE(isFull, isVerbose, dct4_24x24 , (mips.out0.i32, mips.inp0.i32, dct4_32_##N, 3), fout, "N="#N", scalingOpt=3"     , prf_cycle);
#define PROFILE_MDCT_32X16(isFull, isVerbose, N)  PROFILE_SIMPLE(isFull, isVerbose, mdct_32x16 , (mips.out0.i32, mips.inp0.i32, mdct_16_##N, 3), fout, "N="#N", scalingOpt=3"     , prf_cycle);
#define PROFILE_MDCT_32X32(isFull, isVerbose, N)  PROFILE_SIMPLE(isFull, isVerbose, mdct_32x32 , (mips.out0.i32, mips.inp0.i32, mdct_32_##N, 3), fout, "N="#N", scalingOpt=3"     , prf_cycle);
#define PROFILE_MDCT_24X24(isFull, isVerbose, N)  PROFILE_SIMPLE(isFull, isVerbose, mdct_24x24 , (mips.out0.i32, mips.inp0.i32, mdct_32_##N, 3), fout, "N="#N", scalingOpt=3"     , prf_cycle);
#define PROFILE_IMDCT_32X16(isFull, isVerbose, N) PROFILE_SIMPLE(isFull, isVerbose, imdct_32x16, (mips.out0.i32, mips.inp0.i32, mdct_16_##N, 3), fout, "N="#N", scalingOpt=3"     , prf_cycle);
#define PROFILE_IMDCT_32X32(isFull, isVerbose, N) PROFILE_SIMPLE(isFull, isVerbose, imdct_32x32, (mips.out0.i32, mips.inp0.i32, mdct_32_##N, 3), fout, "N="#N", scalingOpt=3"     , prf_cycle);
#define PROFILE_IMDCT_24X24(isFull, isVerbose, N) PROFILE_SIMPLE(isFull, isVerbose, imdct_24x24, (mips.out0.i32, mips.inp0.i32, mdct_32_##N, 3), fout, "N="#N", scalingOpt=3"     , prf_cycle);
#define PROFILE_DCT2D_8X16(isFull, isVerbose, Larg)  PROFILE_NORMALIZED(isFull, isVerbose, dct2d_8x16 , (mips.out0.i16, mips.inp0.u8, dct2d_16_8, Larg, 0), fout, "N=8, L="#Larg", scalingOpt=0", prf_cyclesblk,Larg);
#define PROFILE_IDCT2D_16X8(isFull, isVerbose, Larg) PROFILE_NORMALIZED(isFull, isVerbose, idct2d_16x8, (mips.out0.u8, mips.inp0.i16, idct2d_16_8, Larg, 0), fout, "N=8, L="#Larg", scalingOpt=0", prf_cyclesblk,Larg);

void mips_dct1(int isFull, int isVerbose, FILE * fout)
{
#if 0 // HiFi3/3z API
    PROFILE_SIMPLE(1     , isVerbose, dct_24x24, (mips.out0.i32, mips.inp0.i32, dct2_32_32, 3), fout, "N=32, scalingOpt=3", prf_cycle);
    PROFILE_SIMPLE(isFull, isVerbose, dct_24x24, (mips.out0.i32, mips.inp0.i32, dct2_32_64, 3), fout, "N=64, scalingOpt=3", prf_cycle);
#endif
    PROFILE_SIMPLE(1     , isVerbose, dct_32x16, (mips.out0.i32, mips.inp0.i32, dct2_16_32, 3), fout, "N=32, scalingOpt=3", prf_cycle);
    PROFILE_SIMPLE(isFull, isVerbose, dct_32x16, (mips.out0.i32, mips.inp0.i32, dct2_16_64, 3), fout, "N=64, scalingOpt=3", prf_cycle);
    PROFILE_SIMPLE(1     , isVerbose, dct_32x32, (mips.out0.i32, mips.inp0.i32, dct2_32_32, 3), fout, "N=32, scalingOpt=3", prf_cycle);
    PROFILE_SIMPLE(isFull, isVerbose, dct_32x32, (mips.out0.i32, mips.inp0.i32, dct2_32_64, 3), fout, "N=64, scalingOpt=3", prf_cycle);
    PROFILE_SIMPLE(1     , isVerbose, dct_16x16, (mips.out0.i16, mips.inp0.i16, dct2_16_32, 3), fout, "N=32, scalingOpt=3", prf_cycle);
    PROFILE_SIMPLE(isFull, isVerbose, dct_16x16, (mips.out0.i16, mips.inp0.i16, dct2_16_64, 3), fout, "N=64, scalingOpt=3", prf_cycle);

    PROFILE_DCT4_32X16(1     , isVerbose, 32 );
    PROFILE_DCT4_32X16(isFull, isVerbose, 64 );
    PROFILE_DCT4_32X16(isFull, isVerbose, 128);
    PROFILE_DCT4_32X16(isFull, isVerbose, 256);
    PROFILE_DCT4_32X16(isFull, isVerbose, 512);
    PROFILE_DCT4_32X32(1     , isVerbose, 32 );
    PROFILE_DCT4_32X32(isFull, isVerbose, 64 );
    PROFILE_DCT4_32X32(isFull, isVerbose, 128);
    PROFILE_DCT4_32X32(isFull, isVerbose, 256);
    PROFILE_DCT4_32X32(isFull, isVerbose, 512);
#if 0 // HiFi3/3z API
    PROFILE_DCT4_24X24(1     , isVerbose, 32 );
    PROFILE_DCT4_24X24(isFull, isVerbose, 64 );
    PROFILE_DCT4_24X24(isFull, isVerbose, 128);
    PROFILE_DCT4_24X24(isFull, isVerbose, 256);
    PROFILE_DCT4_24X24(isFull, isVerbose, 512);
#endif
    PROFILE_MDCT_32X16(1     , isVerbose, 32 );
    PROFILE_MDCT_32X16(isFull, isVerbose, 64 );
    PROFILE_MDCT_32X16(isFull, isVerbose, 128);
    PROFILE_MDCT_32X16(isFull, isVerbose, 256);
    PROFILE_MDCT_32X16(isFull, isVerbose, 512);
    PROFILE_MDCT_32X32(1     , isVerbose, 32 );
    PROFILE_MDCT_32X32(isFull, isVerbose, 64 );
    PROFILE_MDCT_32X32(isFull, isVerbose, 128);
    PROFILE_MDCT_32X32(isFull, isVerbose, 256);
    PROFILE_MDCT_32X32(isFull, isVerbose, 512);
#if 0 // HiFi3/3z API
    PROFILE_MDCT_24X24(1     , isVerbose, 32 );
    PROFILE_MDCT_24X24(isFull, isVerbose, 64 );
    PROFILE_MDCT_24X24(isFull, isVerbose, 128);
    PROFILE_MDCT_24X24(isFull, isVerbose, 256);
    PROFILE_MDCT_24X24(isFull, isVerbose, 512);
#endif
    PROFILE_IMDCT_32X16(1     , isVerbose, 32 );
    PROFILE_IMDCT_32X16(isFull, isVerbose, 64 );
    PROFILE_IMDCT_32X16(isFull, isVerbose, 128);
    PROFILE_IMDCT_32X16(isFull, isVerbose, 256);
    PROFILE_IMDCT_32X16(isFull, isVerbose, 512);
    PROFILE_IMDCT_32X32(1     , isVerbose, 32 );
    PROFILE_IMDCT_32X32(isFull, isVerbose, 64 );
    PROFILE_IMDCT_32X32(isFull, isVerbose, 128);
    PROFILE_IMDCT_32X32(isFull, isVerbose, 256);
    PROFILE_IMDCT_32X32(isFull, isVerbose, 512);
#if 0 // HiFi3/3z API
    PROFILE_IMDCT_24X24(1     , isVerbose, 32 );
    PROFILE_IMDCT_24X24(isFull, isVerbose, 64 );
    PROFILE_IMDCT_24X24(isFull, isVerbose, 128);
    PROFILE_IMDCT_24X24(isFull, isVerbose, 256);
    PROFILE_IMDCT_24X24(isFull, isVerbose, 512);
#endif
    PROFILE_DCT2D_8X16(isFull, isVerbose, 1   );
    PROFILE_DCT2D_8X16(isFull, isVerbose, 32  );
    PROFILE_DCT2D_8X16(1     , isVerbose, 1024);

    PROFILE_IDCT2D_16X8(isFull, isVerbose, 1   );
    PROFILE_IDCT2D_16X8(isFull, isVerbose, 32  );
    PROFILE_IDCT2D_16X8(1     , isVerbose, 1024);
}
