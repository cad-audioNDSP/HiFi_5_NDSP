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
static const int16_t scale_q15 = 7815;                
static const int32_t scale_q31 = 512117760;
static const int16_t refl_q15[9] = 
{
  -19932,21123,-20868,19136,-14564,7276,-19932,21123,-20868
};
static const int32_t refl_q31[9] =
{
  -1306313216,1384367872,-1367562496,1254099968,-954675200,476997120,-1306313216,1384367872,-1367562496
};

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

#define PROFILE_LATR16X16(cond,verb,N,M) OBJ_PROFILE_NORMALIZED_LAT(cond,verb,latr16x16,(M),(objinstance_memory,M,refl_q15,scale_q15),(handle, mips.out0.i16, mips.inp0.i16, N),fout,"N="#N", M="#M,prf_cyclessampleM,N*M);
#define PROFILE_LATR32X16(cond,verb,N,M) OBJ_PROFILE_NORMALIZED_LAT(cond,verb,latr32x16,(M),(objinstance_memory,M,refl_q15,scale_q15),(handle, mips.out0.i32, mips.inp0.i32, N),fout,"N="#N", M="#M,prf_cyclessampleM,N*M);
#define PROFILE_LATR24x24(cond,verb,N,M) OBJ_PROFILE_NORMALIZED_LAT(cond,verb,latr24x24,(M),(objinstance_memory,M,refl_q31,scale_q31),(handle, mips.out0.i32, mips.inp0.i32, N),fout,"N="#N", M="#M,prf_cyclessampleM,N*M);
#define PROFILE_LATR32X32(cond,verb,N,M) OBJ_PROFILE_NORMALIZED_LAT(cond,verb,latr32x32,(M),(objinstance_memory,M,refl_q31,scale_q31),(handle, mips.out0.i32, mips.inp0.i32, N),fout,"N="#N", M="#M,prf_cyclessampleM,N*M);

void mips_iirlt1(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_LATR16X16(isFull, isVerbose, 256, 1)
    PROFILE_LATR16X16(isFull, isVerbose, 256, 2)
    PROFILE_LATR16X16(isFull, isVerbose, 256, 3)
    PROFILE_LATR16X16(isFull, isVerbose, 256, 4)
    PROFILE_LATR16X16(isFull, isVerbose, 256, 5)
    PROFILE_LATR16X16(isFull, isVerbose, 256, 6)
    PROFILE_LATR16X16(isFull, isVerbose, 256, 7)
    PROFILE_LATR16X16(     1, isVerbose, 256, 8)
    PROFILE_LATR16X16(isFull, isVerbose, 256, 9)
    PROFILE_LATR16X16(isFull, isVerbose, 80, 6)

    PROFILE_LATR32X16(isFull, isVerbose, 256, 1)
    PROFILE_LATR32X16(isFull, isVerbose, 256, 2)
    PROFILE_LATR32X16(isFull, isVerbose, 256, 3)
    PROFILE_LATR32X16(isFull, isVerbose, 256, 4)
    PROFILE_LATR32X16(isFull, isVerbose, 256, 5)
    PROFILE_LATR32X16(isFull, isVerbose, 256, 6)
    PROFILE_LATR32X16(isFull, isVerbose, 256, 7)
    PROFILE_LATR32X16(     1, isVerbose, 256, 8)
    PROFILE_LATR32X16(isFull, isVerbose, 256, 9)
    PROFILE_LATR32X16(isFull, isVerbose, 80, 6)
#if 0 //HiFi3/3z API
    PROFILE_LATR24x24(isFull, isVerbose, 256, 1)
    PROFILE_LATR24x24(isFull, isVerbose, 256, 2)
    PROFILE_LATR24x24(isFull, isVerbose, 256, 3)
    PROFILE_LATR24x24(isFull, isVerbose, 256, 4)
    PROFILE_LATR24x24(isFull, isVerbose, 256, 5)
    PROFILE_LATR24x24(isFull, isVerbose, 256, 6)
    PROFILE_LATR24x24(isFull, isVerbose, 256, 7)
    PROFILE_LATR24x24(     1, isVerbose, 256, 8)
    PROFILE_LATR24x24(isFull, isVerbose, 256, 9)
    PROFILE_LATR24x24(isFull, isVerbose, 80, 6)
#endif
    PROFILE_LATR32X32(isFull, isVerbose, 256, 1)
    PROFILE_LATR32X32(isFull, isVerbose, 256, 2)
    PROFILE_LATR32X32(isFull, isVerbose, 256, 3)
    PROFILE_LATR32X32(isFull, isVerbose, 256, 4)
    PROFILE_LATR32X32(isFull, isVerbose, 256, 5)
    PROFILE_LATR32X32(isFull, isVerbose, 256, 6)
    PROFILE_LATR32X32(isFull, isVerbose, 256, 7)
    PROFILE_LATR32X32(     1, isVerbose, 256, 8)
    PROFILE_LATR32X32(isFull, isVerbose, 256, 9)
}
