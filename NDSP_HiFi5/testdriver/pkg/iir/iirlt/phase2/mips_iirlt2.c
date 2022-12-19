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
* Test module for testing cycle performance (Lattice Filters)
*/
#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(iir)
#include "mips.h"

/* IIR performance measurement tests */

//----latr--- 
static const float32_t reflf[9] =
{
  -0.60830f,0.64465f,-0.63682f,0.58399f,-0.44456f,0.22212f,-0.60830f,0.64465f,-0.63682f
};
static const float32_t scalef = 0.23847f;

#define OBJ_PROFILE_NORMALIZED_LAT(_cond,_verb,_objname, _a_arg, _i_arg, _p_arg, _file, _info_,_fmt, _norm)  \
{                                                                               \
   _objname##_handle_t handle=NULL;                                             \
  int isPresent;                                                                \
  isPresent =IS_PRESENT(_objname##_alloc);                                      \
  isPresent|=IS_PRESENT(_objname##_init);                                       \
  isPresent |= IS_PRESENT(_objname##_process);                                  \
  if (isPresent ) handle = _objname##_init _i_arg;                              \
  if (handle == NULL) isPresent = 0;                                            \
  PROFILE_NORMALIZED(_cond,_verb,_objname##_process, _p_arg, _file, _info_, _fmt, _norm)    \
}

/* IIR performance measurement tests */
#define PROFILE_LATRF(cond,verb,N,M)     OBJ_PROFILE_NORMALIZED_LAT(cond,verb,latrf,(M),(objinstance_memory,M,reflf,scalef),(handle, mips.out0.f32, mips.inp0.f32, N),fout,"N="#N", M="#M,prf_cyclessampleM,N*M);

void mips_iirlt2(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_LATRF(isFull, isVerbose, 256, 1)
    PROFILE_LATRF(isFull, isVerbose, 256, 2)
    PROFILE_LATRF(isFull, isVerbose, 256, 3)
    PROFILE_LATRF(isFull, isVerbose, 256, 4)
    PROFILE_LATRF(isFull, isVerbose, 256, 5)
    PROFILE_LATRF(isFull, isVerbose, 256, 6)
    PROFILE_LATRF(isFull, isVerbose, 256, 7)
    PROFILE_LATRF(     1, isVerbose, 256, 8)
    PROFILE_LATRF(isFull, isVerbose, 256, 9)
    PROFILE_LATRF(isFull, isVerbose, 80, 6)
}
