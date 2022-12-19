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
 * constants and definitions for testing cycle performance of library
 */

#ifndef MIPS_H__
#define MIPS_H__
#include "types.h"
#include "config.h"
#include "profiler.h"
#include "utils.h"
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include "malloc16.h"
/*#include "hwoptions.h"*/

/* Target function I/O/Scratch buffer size, in 32-byte units */
#define PROFILE_DATA_INP_SIZE   8192//32768
#define PROFILE_DATA_OUT_SIZE   9192//65536
#define PROFILE_SCRATCH_SIZE    24576

typedef union tag_profiler_dataInp
{
  uint8_t           u8[PROFILE_DATA_INP_SIZE*32];
  int8_t            i8[PROFILE_DATA_INP_SIZE*32];
  int16_t          i16[PROFILE_DATA_INP_SIZE*16];
  float16_t        f16[PROFILE_DATA_INP_SIZE*16];
  int32_t          i32[PROFILE_DATA_INP_SIZE* 8];
  float32_t        f32[PROFILE_DATA_INP_SIZE* 8];
  int64_t          i64[PROFILE_DATA_INP_SIZE* 4];
  float64_t        f64[PROFILE_DATA_INP_SIZE* 4];
  complex_fract16 ci16[PROFILE_DATA_INP_SIZE* 8];
  complex_float16 cf16[PROFILE_DATA_INP_SIZE* 8];
  complex_fract32 ci32[PROFILE_DATA_INP_SIZE* 4];
  complex_fract64 ci64[PROFILE_DATA_INP_SIZE* 2];
  complex_float   cf32[PROFILE_DATA_INP_SIZE* 4];
  complex_double  cf64[PROFILE_DATA_INP_SIZE* 2];
}
tProfiler_dataInp;

typedef union tag_profiler_dataOut
{
  uint8_t           u8[PROFILE_DATA_OUT_SIZE*32];
  int8_t            i8[PROFILE_DATA_OUT_SIZE*32];
  int16_t          i16[PROFILE_DATA_OUT_SIZE*16];
  float16_t        f16[PROFILE_DATA_OUT_SIZE*16];
  int32_t          i32[PROFILE_DATA_OUT_SIZE* 8];
  float32_t        f32[PROFILE_DATA_OUT_SIZE* 8];
  int64_t          i64[PROFILE_DATA_OUT_SIZE* 4];
  float64_t        f64[PROFILE_DATA_OUT_SIZE* 4];
  complex_fract16 ci16[PROFILE_DATA_OUT_SIZE* 8];
  complex_float16 cf16[PROFILE_DATA_OUT_SIZE* 8];
  complex_fract32 ci32[PROFILE_DATA_OUT_SIZE* 4];
  complex_fract64 ci64[PROFILE_DATA_OUT_SIZE* 2];
  complex_float   cf32[PROFILE_DATA_OUT_SIZE* 4];
  complex_double  cf64[PROFILE_DATA_OUT_SIZE* 2];
}
tProfiler_dataOut;

typedef union tag_profiler_scratch
{
  uint8_t           u8[PROFILE_SCRATCH_SIZE*32];
  int8_t            i8[PROFILE_SCRATCH_SIZE*32];
  int16_t          i16[PROFILE_SCRATCH_SIZE*16];
  float16_t        f16[PROFILE_SCRATCH_SIZE*16];
  int32_t          i32[PROFILE_SCRATCH_SIZE* 8];
  float32_t        f32[PROFILE_SCRATCH_SIZE* 8];
  int64_t          i64[PROFILE_SCRATCH_SIZE* 4];
  float64_t        f64[PROFILE_SCRATCH_SIZE* 4];
  complex_fract16 ci16[PROFILE_SCRATCH_SIZE* 8];
  complex_float16 cf16[PROFILE_SCRATCH_SIZE* 8];
  complex_fract32 ci32[PROFILE_SCRATCH_SIZE* 4];
  complex_fract64 ci64[PROFILE_SCRATCH_SIZE* 2];
  complex_float   cf32[PROFILE_SCRATCH_SIZE* 4];
  complex_double  cf64[PROFILE_SCRATCH_SIZE* 2];
}
tProfiler_scratch;

#define OBJINSTANCE_SIZE  ( 10*4096)

/* external data (all aligned) */
extern uint8_t objinstance_memory[OBJINSTANCE_SIZE];/* PLACE_IN_DRAM1; */
#if 0
extern tProfiler_scratch scratch0;
extern tProfiler_dataInp inp0,inp1,inp2;
extern tProfiler_dataOut out0,out1,out2;
#else
typedef struct
{
    tProfiler_dataInp inp0,inp1,inp2;
    tProfiler_dataOut out0,out1,out2;
    tProfiler_scratch scratch0;
}
tProfilerData;
extern tProfilerData mips;
#endif

#define PROFILE_SIMPLE(_cond, _verbose, _f_name, _f_arg, _file, _info, _fmt) \
{                                                                 \
  if (_cond)                                                      \
  {                                                               \
      if (_verbose)                                               \
      {                                                           \
           extern const char ANNOTATE_FUN_REF(_f_name)[];         \
           perf_info(_file, "%-12s\t", #_f_name);                 \
           perf_info(_file,"%s\t", ANNOTATE_FUN_REF(_f_name));    \
           perf_info(_file, "%-16s\t",_info);                     \
      } else                                                      \
      {                                                           \
           perf_info(_file, "%-12s\t%-16s\t", #_f_name,_info);    \
      }                                                           \
      if (IS_PRESENT(_f_name))                                    \
      {                                                           \
           PROFILE(_f_name _f_arg);                               \
           perf_info(_file, _fmt, profiler_clk);                  \
           perf_info(_file, "\n");                                \
      }                                                           \
      else                                                        \
      {                                                           \
           perf_info(_file, "NOT TESTED\n");                      \
      }                                                           \
  }                                                               \
}

#define PROFILE_INVERTED(_cond, _verbose, _f_name, _f_arg, _file, _info, _fmt, _norm) \
{                                                                 \
  if (_cond)                                                      \
  {                                                               \
      float32_t __norm;                                           \
      if (_verbose)                                               \
      {                                                           \
           extern const char ANNOTATE_FUN_REF(_f_name)[];         \
           perf_info(_file, "%-12s\t", #_f_name);                 \
           perf_info(_file,"%s\t", ANNOTATE_FUN_REF(_f_name));    \
           perf_info(_file, "%-16s\t",_info);                     \
      } else                                                      \
      {                                                           \
           perf_info(_file, "%-12s\t%-16s\t", #_f_name,_info);    \
      }                                                           \
      if (IS_PRESENT(_f_name))                                    \
      {                                                           \
          PROFILE(_f_name _f_arg);                                \
          __norm=(float32_t)(_norm)/profiler_clk;                 \
          perf_info(_file, _fmt, profiler_clk,__norm);            \
          perf_info(_file, "\n");                                 \
      }                                                           \
      else                                                        \
      {                                                           \
           perf_info(_file, "NOT TESTED\n");                      \
      }                                                           \
  }                                                               \
}

#define PROFILE_NORMALIZED(_cond, _verbose, _f_name, _f_arg, _file, _info, _fmt, _norm) \
{                                                                 \
  if (_cond)                                                      \
  {                                                               \
      float32_t __norm;                                           \
      if (_verbose)                                               \
      {                                                           \
           extern const char ANNOTATE_FUN_REF(_f_name)[];         \
           perf_info(_file, "%-12s\t", #_f_name);                 \
           perf_info(_file,"%s\t", ANNOTATE_FUN_REF(_f_name));    \
           perf_info(_file, "%-16s\t",_info);                     \
      } else                                                      \
      {                                                           \
           perf_info(_file, "%-12s\t%-16s\t", #_f_name,_info);    \
      }                                                           \
      if (IS_PRESENT(_f_name))                                    \
      {                                                           \
           PROFILE(_f_name _f_arg);                               \
           __norm=profiler_clk/(float32_t)(_norm);                \
           perf_info(_file, _fmt, profiler_clk,__norm);           \
           perf_info(_file, "\n");                                \
      }                                                           \
      else                                                        \
      {                                                           \
           perf_info(_file, "NOT TESTED\n");                      \
      }                                                           \
  }                                                               \
}

#define OBJ_PROFILE_INVERTED(_cond, _verbose,_objname, _a_arg, _i_arg, _p_arg, _file,_info_,_fmt, _norm)\
{                                                                                  \
  int isPresent;                                                                   \
  _objname##_handle_t _objname=NULL;                                               \
  isPresent =IS_PRESENT(_objname##_alloc);                         \
  isPresent|=IS_PRESENT(_objname##_init);                          \
  isPresent|=IS_PRESENT(_objname##_process);                       \
  if (isPresent )     _objname = _objname##_init _i_arg;                           \
  if (_objname==NULL) isPresent=0;                                                 \
  PROFILE_INVERTED(_cond, _verbose,_objname##_process, _p_arg, _file, _info_, _fmt, _norm)         \
}


#define PROFILE_INVERTED_FFT(_cond, _verbose,_f_name,N, _f_arg, _file, _fmt, _norm) \
 if (IS_PRESENT(_f_name))                         \
 {                                                                \
    PROFILE_INVERTED(_cond, _verbose, _f_name, _f_arg, _file, "N="#N, _fmt, _norm) \
 }                                                                \
 else                                                             \
 {                                                                \
  perf_info(_file, #_f_name"\t\tNOT TESTED\n");                   \
 }

#define PROFILE_INVERTED_FFT_SC(_cond, _verbose,_f_name,N, i, _f_arg, _file, _fmt, _norm)           \
if (IS_PRESENT(_f_name))                                            \
 {                                                                                  \
    PROFILE_INVERTED(_cond, _verbose,_f_name, _f_arg, _file, "N=" #N ", scaling=" #i , _fmt, _norm) \
 }                                                                                  \
 else                                                                               \
 {                                                                                  \
 perf_info(_file, #_f_name"\t\tNOT TESTED\n");                                      \
 }

/* perf_info
  Works as printf() but duplicates output into a file (if f is not NULL)
 */
int perf_info(FILE * f, const char * fmt, ...);

/* global constants */
extern const char* prf_ptscycle;
extern const char* prf_ptscycle2;
extern const char* prf_ptscycle3;
extern const char* prf_maccycle;
extern const char* prf_cyclesmtx;
extern const char* prf_cyclesblk;
extern const char* prf_cycle;
extern const char* prf_cyclespts;
extern const char* prf_bytescycle;
extern const char* prf_bitscycle;
extern const char* prf_cyclesbqd;
extern const char* prf_cyclessampleM;
#if 0
extern const char* prf_samplecycle;
extern const char* prf_cyclessample;
#endif

#endif
