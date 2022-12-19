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
* Test module for testing cycle performance (Scalar Math)
*/

#include "mips.h"
#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(math)

void mips_maths2(int isFull, int isVerbose, FILE * fout)
{
  PROFILE_SIMPLE(isFull, isVerbose, scl_int2float,(9621325  ,1),fout,""     ,prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_float2int,(9621325.f,1),fout,""     ,prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_sinef     ,(1.2f),      fout,""     ,prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_cosinef   ,(1.2f),      fout,""     ,prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_tanf      ,(0.4f),      fout,"x=0.4",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_tanf      ,(1.2f),      fout,"x=1.2",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_log2f     ,(1.2f),      fout,""     ,prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_log10f    ,(1.2f),      fout,""     ,prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_lognf     ,(1.2f),      fout,""     ,prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_antilog2f ,(1.2f),      fout,""     ,prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_antilog10f,(1.2f),      fout,""     ,prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_antilognf ,(1.2f),      fout,""     ,prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_powf      ,(1.f, 1.f)    , fout, "x=1 y=1"        , prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_powf      ,(1.25f, 0.75f), fout, "x=1.25 y=0.75"  , prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_atanf     ,(0.7f),      fout,"x=0.7",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_atanf     ,(1.3f),      fout,"x=1.3",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_atan2f,    (1.2f,2.f),  fout,""     ,prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_sigmoidf,  (1.2f),      fout,"",     prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_tanhf,     (1.2f),      fout,"",     prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_reluf,     (1.2f,2.f),  fout,""     ,prf_cycle);
}
