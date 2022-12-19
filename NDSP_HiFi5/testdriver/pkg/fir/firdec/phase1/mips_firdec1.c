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
* Test module for testing cycle performance (Decimation)
*/

#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(fir)
#include "mips.h"

#define FIRDEC16X16_PROFILE(cond,verb,D,N,M)   OBJ_PROFILE_INVERTED(cond,verb,firdec16x16  ,(D,M),(objinstance_memory, D,M, mips.inp2.i16),(firdec16x16  ,mips.out2.i16, mips.inp2.i16,N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M);
#define FIRDEC32X32_PROFILE(cond,verb,D,N,M)   OBJ_PROFILE_INVERTED(cond,verb,firdec32x32  ,(D,M),(objinstance_memory, D,M, mips.inp2.i32),(firdec32x32  ,mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M);
#define FIRDEC32X32EP_PROFILE(cond,verb,D,N,M) OBJ_PROFILE_INVERTED(cond,verb,firdec32x32ep,(D,M),(objinstance_memory, D,M, mips.inp2.i32),(firdec32x32ep,mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M);
#define FIRDEC32X16_PROFILE(cond,verb,D,N,M)   OBJ_PROFILE_INVERTED(cond,verb,firdec32x16  ,(D,M),(objinstance_memory, D,M, mips.inp2.i16),(firdec32x16  ,mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M);
#define FIRDEC24X24_PROFILE(cond,verb,D,N,M)   OBJ_PROFILE_INVERTED(cond,verb,firdec24x24  ,(D,M),(objinstance_memory, D,M, mips.inp2.i32),(firdec24x24  ,mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M);

void mips_firdec16x16(int isFull, int isVerbose, FILE * fout)
{
    FIRDEC16X16_PROFILE(isFull, isVerbose, 2, 1024,   2)
    FIRDEC16X16_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 2, 1024, 261)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 2,   80, 256)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 3, 1024,   2)
    FIRDEC16X16_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 3, 1024, 261)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 4, 1024,   2)
    FIRDEC16X16_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 4, 1024, 261)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRDEC16X16_PROFILE(isFull, isVerbose, 7, 1024, 260)
}

void mips_firdec32x16(int isFull, int isVerbose, FILE * fout)
{
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 1024, 2)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 1024, 4 )
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 1024, 8 )
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 1024, 16)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 1024, 32)
    FIRDEC32X16_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2, 1024, 261)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 2,   80, 256)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 3, 1024,   2)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 3, 1024,  4 )
    FIRDEC32X16_PROFILE(isFull, isVerbose, 3, 1024,  8 )
    FIRDEC32X16_PROFILE(isFull, isVerbose, 3, 1024,  16)
    FIRDEC32X16_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 3, 1024, 261)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 4, 1024,   2)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 4, 1024,   4)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 4, 1024,   8)
    FIRDEC32X16_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 4, 1024, 261)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 6, 1024, 256)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 6, 1024, 260)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRDEC32X16_PROFILE(isFull, isVerbose, 7, 1024, 260)
}

void mips_firdec24x24(int isFull, int isVerbose, FILE * fout)
{
#if 0// for HiFi3/3z
    FIRDEC24X24_PROFILE(isFull, isVerbose, 2, 1024,   2)
    FIRDEC24X24_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 2, 1024, 261)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 3, 1024,   2)
    FIRDEC24X24_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 3, 1024, 261)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 4, 1024,   2)
    FIRDEC24X24_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 4, 1024, 261)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 7, 1024, 260)
    FIRDEC24X24_PROFILE(isFull, isVerbose, 2,   80, 256)
#endif
}

void mips_firdec32x32(int isFull, int isVerbose, FILE * fout)
{
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2, 1024,   2)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2, 1024,   4)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2, 1024,   8)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2, 1024,  16)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2, 1024,  32)
    FIRDEC32X32_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2, 1024, 261)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 3, 1024,   2)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 3, 1024,  4 )
    FIRDEC32X32_PROFILE(isFull, isVerbose, 3, 1024,  8 )
    FIRDEC32X32_PROFILE(isFull, isVerbose, 3, 1024,  16)
    FIRDEC32X32_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 3, 1024, 261)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 4, 1024,   2)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 4, 1024,   4)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 4, 1024,   8)
    FIRDEC32X32_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 4, 1024, 261)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 6, 1024, 256)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 6, 1024, 260)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 7, 1024, 260)
    FIRDEC32X32_PROFILE(isFull, isVerbose, 2,   80, 256)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2,   80, 256)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 1024,   2)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 1024,  4 )
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 1024,  8 )
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 1024,  16)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 1024,  32)
    FIRDEC32X32EP_PROFILE(     1, isVerbose, 2, 1024, 256)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 1024, 260)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 2, 1024, 261)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 3, 1024,   2)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 3, 1024,  4 )
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 3, 1024,  8 )
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 3, 1024,  16)
    FIRDEC32X32EP_PROFILE(     1, isVerbose, 3, 1024, 256)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 3, 1024, 260)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 3, 1024, 261)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 4, 1024,   2)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 4, 1024,   4)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 4, 1024,   8)
    FIRDEC32X32EP_PROFILE(     1, isVerbose, 4, 1024, 256)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 4, 1024, 260)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 4, 1024, 261)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 5, 1024, 256)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 5, 1024, 260)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 6, 1024, 256)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 6, 1024, 260)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 7, 1024, 256)
    FIRDEC32X32EP_PROFILE(isFull, isVerbose, 7, 1024, 260)
}

void mips_firdec1(int isFull, int isVerbose, FILE * fout)
{
  mips_firdec16x16(isFull, isVerbose, fout);
  mips_firdec32x16(isFull, isVerbose, fout);
  mips_firdec24x24(isFull, isVerbose, fout);
  mips_firdec32x32(isFull, isVerbose, fout);
}
