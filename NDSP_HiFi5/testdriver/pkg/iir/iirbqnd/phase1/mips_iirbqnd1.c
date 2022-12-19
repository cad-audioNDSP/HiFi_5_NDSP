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
* Test module for testing cycle performance (Biquad Filters)
*/
#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(iir)
#include "mips.h"

/* IIR performance measurement tests */
//----bqriir---
static const int16_t coef_sos_16[5*8] =
{ // b0      b1      b2      a1     a2
   16384,  31383,  16384,  12682, 14622,
   16384, -31383,  16384, -12682, 14622,
   16384,  32215,  16384,   7625, 11691,
   16384, -32215,  16384,  -7625, 11691,
   16384,      0, -16384,      0, 10377,
   16384,  31383,  16384,  12682, 14622,
   16384, -31383,  16384, -12682, 14622,
   16384,  32215,  16384,   7625, 11691,
};
static const int32_t coef_sos_32[5*8] =
{ //  b0          b1         b2         a1         a2
  1073741824, 2056704919, 1073741824, 831104644,958261518,
  1073741824,-2056704919, 1073741824,-831104644,958261518,
  1073741824, 2111239901, 1073741824, 499713750,766176384,
  1073741824,-2111239901, 1073741824,-499713750,766176384,
  1073741824,          0,-1073741824,         0,680063938,
  1073741824, 2056704919, 1073741824, 831104644,958261518,
  1073741824,-2056704919, 1073741824,-831104644,958261518,
  1073741824, 2111239901, 1073741824, 499713750,766176384,
};

static const int16_t coef_g[8] =
{
  2460,19682,9107,9107,22170,2460,19682,9107
};

#define OBJ_PROFILE_NORMALIZED_IIR(_cond,_verb,_objname, _a_arg, _i_arg, _p_arg, _file, _info_,_fmt, _norm)  \
{                                                                               \
   _objname##_handle_t handle=NULL;                                             \
  int isPresent;                                                                \
  isPresent =IS_PRESENT(_objname##_alloc);                                      \
  isPresent|=IS_PRESENT(_objname##_init);                                       \
  isPresent |= IS_PRESENT(_objname);                                            \
  if (isPresent )     handle = _objname##_init _i_arg;                          \
  if (handle == NULL) isPresent = 0;                                            \
  PROFILE_NORMALIZED(_cond,_verb,_objname, _p_arg, _file, _info_, _fmt, _norm)  \
}

/* IIR performance measurement tests */
#define PROFILE_BQRIIR16X16_DF1_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir16x16_df1_nd,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i16, mips.inp1.i16, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR16X16_DF2_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir16x16_df2_nd,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i16, mips.inp1.i16, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR32X16_DF1_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir32x16_df1_nd,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR32X16_DF2_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir32x16_df2_nd,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR32X32_DF1_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir32x32_df1_nd,(M),(objinstance_memory,M,coef_sos_32, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIR32X32_DF2_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriir32x32_df2_nd,(M),(objinstance_memory,M,coef_sos_32, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_STEREO_BQRIIR16X16_DF1_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,stereo_bqriir16x16_df1_nd,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i16, mips.inp1.i16, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_STEREO_BQRIIR32X16_DF1_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,stereo_bqriir32x16_df1_nd,(M),(objinstance_memory,M,coef_sos_16, coef_g, gain,coef_sos_16, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);
#define PROFILE_STEREO_BQRIIR32X32_DF1_ND(cond,verb,N,M,gain) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,stereo_bqriir32x32_df1_nd,(M),(objinstance_memory,M,coef_sos_32, coef_g, gain,coef_sos_32, coef_g, gain),(handle, mips.scratch0.u8, mips.out0.i32, mips.inp1.i32, N),fout,"N=" #N ", M=" #M ", gain=" #gain ,prf_cyclesbqd,N*M);

void mips_iirbqnd1(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 1, 0)
    PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 2, 1)
    PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 3, 0)
    PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 4, 1)
    PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 5, 0)
    PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 6, 1)
    PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 7, 0)
    PROFILE_BQRIIR16X16_DF1_ND(1, isVerbose, 256, 8, 1)
    PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 80, 5, 0)
    PROFILE_BQRIIR16X16_DF1_ND(isFull, isVerbose, 80, 5, 1)

    PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 256, 1, 0)
    PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 256, 2, 1)
    PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 256, 3, 0)
    PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 256, 4, 1)
    PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 256, 5, 0)
    PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 256, 6, 1)
    PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 256, 7, 0)
    PROFILE_BQRIIR16X16_DF2_ND(1, isVerbose, 256, 8, 1)
    PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 80, 5, 0)
    PROFILE_BQRIIR16X16_DF2_ND(isFull, isVerbose, 80, 5, 1)

    PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 1, 0)
    PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 2, 1)
    PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 3, 0)
    PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 4, 1)
    PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 5, 0)
    PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 6, 1)
    PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 7, 0)
    PROFILE_BQRIIR32X16_DF1_ND(1, isVerbose, 256, 8, 1)
    PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 80, 5, 0)
    PROFILE_BQRIIR32X16_DF1_ND(isFull, isVerbose, 80, 5, 1)

    PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 256, 1, 0)
    PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 256, 2, 1)
    PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 256, 3, 0)
    PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 256, 4, 1)
    PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 256, 5, 0)
    PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 256, 6, 1)
    PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 256, 7, 0)
    PROFILE_BQRIIR32X16_DF2_ND(1, isVerbose, 256, 8, 1)
    PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 80, 5, 0)
    PROFILE_BQRIIR32X16_DF2_ND(isFull, isVerbose, 80, 5, 1)

    PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 1, 0)
    PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 2, 1)
    PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 3, 0)
    PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 4, 1)
    PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 5, 0)
    PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 6, 1)
    PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 7, 0)
    PROFILE_BQRIIR32X32_DF1_ND(1, isVerbose, 256, 8, 1)
    PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 80, 5, 0)
    PROFILE_BQRIIR32X32_DF1_ND(isFull, isVerbose, 80, 5, 1)

    PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 256, 1, 0)
    PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 256, 2, 1)
    PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 256, 3, 0)
    PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 256, 4, 1)
    PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 256, 5, 0)
    PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 256, 6, 1)
    PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 256, 7, 0)
    PROFILE_BQRIIR32X32_DF2_ND(1, isVerbose, 256, 8, 1)
    PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 80, 5, 0)
    PROFILE_BQRIIR32X32_DF2_ND(isFull, isVerbose, 80, 5, 1)

    PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 1, 0)
    PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 2, 1)
    PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 3, 0)
    PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 4, 1)
    PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 5, 0)
    PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 6, 1)
    PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 256, 7, 0)
    PROFILE_STEREO_BQRIIR16X16_DF1_ND(1, isVerbose, 256, 8, 1)
    PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 80, 5, 0)
    PROFILE_STEREO_BQRIIR16X16_DF1_ND(isFull, isVerbose, 80, 5, 1)

    PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 1, 0)
    PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 2, 1)
    PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 3, 0)
    PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 4, 1)
    PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 5, 0)
    PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 6, 1)
    PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 256, 7, 0)
    PROFILE_STEREO_BQRIIR32X16_DF1_ND(1, isVerbose, 256, 8, 1)
    PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 80, 5, 0)
    PROFILE_STEREO_BQRIIR32X16_DF1_ND(isFull, isVerbose, 80, 5, 1)

    PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 1, 0)
    PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 2, 1)
    PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 3, 0)
    PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 4, 1)
    PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 5, 0)
    PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 6, 1)
    PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 256, 7, 0)
    PROFILE_STEREO_BQRIIR32X32_DF1_ND(1, isVerbose, 256, 8, 1)
    PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 80, 5, 0)
    PROFILE_STEREO_BQRIIR32X32_DF1_ND(isFull, isVerbose, 80, 5, 1)
}
