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
static const float32_t coef_sos_f[5*16] =
{ //  b0          b1         b2         a1         a2
1.0000000f, 1.9154557f, 1.0000000f, 0.7740265f, 0.8924506f,
1.0000000f,-1.9154557f, 1.0000000f,-0.7740265f, 0.8924506f,
1.0000000f, 1.9662454f, 1.0000000f, 0.4653947f, 0.7135574f,
1.0000000f,-1.9662454f, 1.0000000f,-0.4653947f, 0.7135574f,
1.0000000f, 0.0000000f,-1.0000000f, 0.0000000f, 0.6333589f,
1.0000000f, 1.9154557f, 1.0000000f, 0.7740265f, 0.8924506f,
1.0000000f,-1.9154557f, 1.0000000f,-0.7740265f, 0.8924506f,
1.0000000f, 1.9662454f, 1.0000000f, 0.4653947f, 0.7135574f,
1.0000000f, 1.9154557f, 1.0000000f, 0.7740265f, 0.8924506f,
1.0000000f,-1.9154557f, 1.0000000f,-0.7740265f, 0.8924506f,
1.0000000f, 1.9662454f, 1.0000000f, 0.4653947f, 0.7135574f,
1.0000000f,-1.9662454f, 1.0000000f,-0.4653947f, 0.7135574f,
1.0000000f, 0.0000000f,-1.0000000f, 0.0000000f, 0.6333589f,
1.0000000f, 1.9154557f, 1.0000000f, 0.7740265f, 0.8924506f,
1.0000000f,-1.9154557f, 1.0000000f,-0.7740265f, 0.8924506f,
1.0000000f, 1.9662454f, 1.0000000f, 0.4653947f, 0.7135574f
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
#define PROFILE_BQRIIRF_DF1(cond,verb,N,M)  OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriirf_df1 ,(M),(objinstance_memory,M,coef_sos_f,1),(handle, mips.out0.f32, mips.inp1.f32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIRF_DF2(cond,verb,N,M)  OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriirf_df2 ,(M),(objinstance_memory,M,coef_sos_f,1),(handle, mips.out0.f32, mips.inp1.f32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);
#define PROFILE_BQRIIRF_DF2T(cond,verb,N,M) OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqriirf_df2t,(M),(objinstance_memory,M,coef_sos_f,1),(handle, mips.out0.f32, mips.inp1.f32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);
#define PROFILE_BQCIIRF_DF1(cond,verb,N,M)  OBJ_PROFILE_NORMALIZED_IIR(cond,verb,bqciirf_df1 ,(M),(objinstance_memory,M,coef_sos_f,1),(handle, mips.out0.cf32, mips.inp1.cf32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);
#define PROFILE_STEREO_BQRIIRF_DF1(cond,verb,N,M)  OBJ_PROFILE_NORMALIZED_IIR(cond,verb,stereo_bqriirf_df1 ,(M),(objinstance_memory,M,coef_sos_f,1,coef_sos_f,1),(handle, mips.out0.f32, mips.inp1.f32, N),fout,"N=" #N ", M=" #M ,prf_cyclesbqd,N*M);

void mips_iirbq2(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_BQRIIRF_DF1(isFull, isVerbose, 512, 1) 
    PROFILE_BQRIIRF_DF1(isFull, isVerbose, 512, 2) 
    PROFILE_BQRIIRF_DF1(isFull, isVerbose, 512, 3) 
    PROFILE_BQRIIRF_DF1(isFull, isVerbose, 512, 4) 
    PROFILE_BQRIIRF_DF1(isFull, isVerbose, 512, 8) 
    PROFILE_BQRIIRF_DF1(isFull, isVerbose, 512,12) 
    PROFILE_BQRIIRF_DF1(     1, isVerbose, 512,16) 

    PROFILE_BQRIIRF_DF2(isFull, isVerbose, 512, 1) 
    PROFILE_BQRIIRF_DF2(isFull, isVerbose, 512, 2) 
    PROFILE_BQRIIRF_DF2(isFull, isVerbose, 512, 3) 
    PROFILE_BQRIIRF_DF2(isFull, isVerbose, 512, 4) 
    PROFILE_BQRIIRF_DF2(isFull, isVerbose, 512, 8) 
    PROFILE_BQRIIRF_DF2(isFull, isVerbose, 512,12) 
    PROFILE_BQRIIRF_DF2(     1, isVerbose, 512,16) 

    PROFILE_BQRIIRF_DF2T(isFull, isVerbose, 512, 1) 
    PROFILE_BQRIIRF_DF2T(isFull, isVerbose, 512, 2) 
    PROFILE_BQRIIRF_DF2T(isFull, isVerbose, 512, 3) 
    PROFILE_BQRIIRF_DF2T(isFull, isVerbose, 512, 4) 
    PROFILE_BQRIIRF_DF2T(isFull, isVerbose, 512, 8) 
    PROFILE_BQRIIRF_DF2T(isFull, isVerbose, 512,12) 
    PROFILE_BQRIIRF_DF2T(     1, isVerbose, 512,16) 

    PROFILE_BQCIIRF_DF1(isFull, isVerbose, 512, 1) 
    PROFILE_BQCIIRF_DF1(isFull, isVerbose, 512, 2) 
    PROFILE_BQCIIRF_DF1(isFull, isVerbose, 512, 3) 
    PROFILE_BQCIIRF_DF1(isFull, isVerbose, 512, 4) 
    PROFILE_BQCIIRF_DF1(isFull, isVerbose, 512, 8) 
    PROFILE_BQCIIRF_DF1(isFull, isVerbose, 512,12) 
    PROFILE_BQCIIRF_DF1(     1, isVerbose, 512,16) 

    PROFILE_STEREO_BQRIIRF_DF1(isFull, isVerbose, 512, 1) 
    PROFILE_STEREO_BQRIIRF_DF1(isFull, isVerbose, 512, 2) 
    PROFILE_STEREO_BQRIIRF_DF1(isFull, isVerbose, 512, 3) 
    PROFILE_STEREO_BQRIIRF_DF1(isFull, isVerbose, 512, 4) 
    PROFILE_STEREO_BQRIIRF_DF1(isFull, isVerbose, 512, 8) 
    PROFILE_STEREO_BQRIIRF_DF1(isFull, isVerbose, 512,12) 
    PROFILE_STEREO_BQRIIRF_DF1(     1, isVerbose, 512,16) 
}
