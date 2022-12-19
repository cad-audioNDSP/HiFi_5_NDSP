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

#define FIRDECF_PROFILE(cond,verb,D,N,M)       OBJ_PROFILE_INVERTED(cond,verb,firdecf      ,(D,M),(objinstance_memory, D,M, mips.inp2.f32),(firdecf      ,mips.out2.f32, mips.inp2.f32,N),fout,"N: " #N "; M: " #M "; D: " #D ,prf_maccycle,N*M);

void mips_firdec2(int isFull, int isVerbose, FILE * fout)
{
  FIRDECF_PROFILE(     1, isVerbose, 2,  1024, 256);
  FIRDECF_PROFILE(isFull, isVerbose, 2,  1024, 512);
  FIRDECF_PROFILE(     1, isVerbose, 3,  1024, 256);
  FIRDECF_PROFILE(isFull, isVerbose, 3,  1024, 512);
  FIRDECF_PROFILE(     1, isVerbose, 4,  1024, 256);
  FIRDECF_PROFILE(isFull, isVerbose, 4,  1024, 512);
  FIRDECF_PROFILE(isFull, isVerbose, 8,  1024, 256);
  FIRDECF_PROFILE(isFull, isVerbose, 8,  1024, 512);
  FIRDECF_PROFILE(isFull, isVerbose, 11, 1024, 256 );
  FIRDECF_PROFILE(isFull, isVerbose, 11, 1024, 512 );
  FIRDECF_PROFILE(isFull, isVerbose, 23, 1024, 256 );
  FIRDECF_PROFILE(isFull, isVerbose, 23, 1024, 512 );
}
