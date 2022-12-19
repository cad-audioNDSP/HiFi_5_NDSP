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
* Test module for testing cycle performance (Filtering)
*/

#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(fir)
#include "mips.h"

#define BKFIRF_PROFILE(cond,verb,N,M)        OBJ_PROFILE_INVERTED(cond,verb,bkfirf,      (M),(objinstance_memory, M, 0, mips.inp2.f32),(bkfirf      , mips.out2.f32, mips.inp2.f32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M);
#define BKFIRAF_PROFILE(cond,verb,N,M)       OBJ_PROFILE_INVERTED(cond,verb,bkfiraf    ,  (M),(objinstance_memory, M, mips.inp2.f32),(bkfiraf    ,  mips.out2.f32, mips.inp2.f32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M);
#define CXFIRF_PROFILE(cond,verb,N,M)        OBJ_PROFILE_INVERTED(cond,verb,cxfirf,      (M),(objinstance_memory, M, 0, mips.inp2.cf32),(cxfirf      , mips.out2.cf32, mips.inp2.cf32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M*4);
#define STEREO_BKFIRF_PROFILE(cond,verb,N,M)        OBJ_PROFILE_INVERTED(cond,verb,stereo_bkfirf,      (M),(objinstance_memory, M, 0, mips.inp2.f32, mips.inp2.f32),(stereo_bkfirf      , mips.out2.f32, mips.inp2.f32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M*2);

void mips_firblk2(int isFull, int isVerbose, FILE * fout)
{
  BKFIRAF_PROFILE(      1, isVerbose, 512,  32);
  BKFIRAF_PROFILE( isFull, isVerbose, 1024,  32);
  BKFIRAF_PROFILE( isFull, isVerbose, 1024, 256);
  BKFIRAF_PROFILE( isFull, isVerbose, 1024, 512);
  BKFIRF_PROFILE(      1, isVerbose, 512,  32);
  BKFIRF_PROFILE( isFull, isVerbose, 1024,  32);
  BKFIRF_PROFILE( isFull, isVerbose, 1024, 256);
  BKFIRF_PROFILE( isFull, isVerbose, 1024, 512);
  STEREO_BKFIRF_PROFILE(      1, isVerbose, 512,  32);
  STEREO_BKFIRF_PROFILE( isFull, isVerbose, 1024,  32);
  STEREO_BKFIRF_PROFILE( isFull, isVerbose, 1024, 256);
  STEREO_BKFIRF_PROFILE( isFull, isVerbose, 1024, 512);
  CXFIRF_PROFILE(      1, isVerbose, 512,32);
  CXFIRF_PROFILE( isFull, isVerbose, 512,256);
}
