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
 * common variables/functions for cycle measurements
 */

#include "mips.h"

#if defined COMPILER_MSVC
#define MEM_ALIGNER  __declspec(align(128))
#else
#define MEM_ALIGNER __attribute__ ((aligned(128)))
#endif

#if defined (COMPILER_XTENSA) && defined(MEM_MODEL)
  #if MEM_MODEL==1
    #if (XCHAL_DATARAM0_SIZE <= (4*1024*1024))
    #error  The local memory size is too low. Please choose The Active Build Target with MEM_MODEL=2 (Release/Debug) to build with system memory.
    #endif	
  #endif
	/* for MEM_MODEL=2 place to sram, otherwise place to dram0/dram0  */
	#if !defined (DRAM0)
		#if MEM_MODEL==2
		#define DRAM0 ".sram.data"
		#elif MEM_MODEL==1
		  #define DRAM0 ".dram0.data"
    #else
		#error unsupported MEM_MODEL
		#endif
	#endif
	#if !defined (DRAM1)
		#if MEM_MODEL==2
		#define DRAM1 ".sram.data"
		#elif MEM_MODEL==1
		#define DRAM1 ".dram0.data"
		#else
		#error unsupported MEM_MODEL
		#endif
	#endif
	#define PLACE_IN_DRAM0 __attribute__ ((section(DRAM0)))
	#define PLACE_IN_DRAM1 __attribute__ ((section(DRAM1)))
#else
	#define PLACE_IN_DRAM0
	#define PLACE_IN_DRAM1
#endif

uint8_t           MEM_ALIGNER objinstance_memory[OBJINSTANCE_SIZE] PLACE_IN_DRAM1;
/*
tProfiler_dataInp MEM_ALIGNER inp0                                 PLACE_IN_DRAM0;
tProfiler_dataInp MEM_ALIGNER inp1                                 PLACE_IN_DRAM1;
tProfiler_dataInp MEM_ALIGNER inp2                                 PLACE_IN_DRAM1;
tProfiler_dataOut MEM_ALIGNER out0                                 PLACE_IN_DRAM0;
tProfiler_dataOut MEM_ALIGNER out1                                 PLACE_IN_DRAM1;
tProfiler_dataOut MEM_ALIGNER out2                                 PLACE_IN_DRAM1;
tProfiler_scratch MEM_ALIGNER scratch0                             PLACE_IN_DRAM0;
*/
tProfilerData MEM_ALIGNER mips PLACE_IN_DRAM0;
/* perf_info
  Works as printf() but duplicates output into a file (if f is not NULL)
 */
int perf_info(FILE * f, const char * fmt, ...)
{
  int res;
  va_list arg;
  va_start(arg, fmt);
  if (f) 
  {
      vfprintf(f, fmt, arg);
      fflush(f);
  }
  res=vprintf(fmt, arg);
  fflush(stdout);
  va_end(arg);
  return res;
}

int main_mips( int phaseNum, int isFull, int isVerbose, int breakOnError )
{
    (void)phaseNum; 
    (void)isFull;
    (void)isVerbose;
    (void)breakOnError;
#ifdef COMPILER_XTENSA
#if XCHAL_HAVE_CCOUNT==0 || !defined(XT_RSR_CCOUNT) || defined(NO_XT_RSR_CCOUNT)
    printf("CCOUNT register not present in the configuration or XT_RSR_CCOUNT not defined! Cycle measurements might be inaccurate!\n");
#endif
#if !defined(MEM_MODEL)
    printf("MEM_MODEL is not defined! Cycle measurements might be inaccurate!\n");
#endif
#endif
    Rand_reset(RAND_RESET_A, RAND_RESET_B);

    Rand_i32(mips.inp0.i32, sizeof(mips.inp0.i32)/sizeof(*mips.inp0.i32));
    Rand_i32(mips.inp1.i32, sizeof(mips.inp1.i32)/sizeof(*mips.inp1.i32));
#if 1
    Rand_i32(mips.inp2.i32, sizeof(mips.inp2.i32)/sizeof(*mips.inp2.i32));
#endif
    return (1);
} /* main_mips() */

/* global constants */
const char* prf_ptscycle      = "%d (%.1f pts/cycle)";
const char* prf_ptscycle2     = "%d (%.2f pts/cycle)";
const char* prf_ptscycle3     = "%d (%.3f pts/cycle)";
const char* prf_maccycle      = "%d (%.1f MACs/cycle)";
const char* prf_cyclesmtx     = "%d (%.1f cycles/matrix)";
const char* prf_cycle         = "%d (cycles)";
const char* prf_cyclespts     = "%d (%.1f cycles/pts)";
const char* prf_bytescycle     ="%d (%.1f bytes/cycle)";
const char* prf_bitscycle      ="%d (%.1f bits/cycle)";
const char* prf_cyclesblk     = "%d (%.1f cycles/block)";
const char* prf_cyclesbqd     =" %d (%.1f cycles/(biquad*pts)";
const char* prf_cyclessampleM =" %d (%.1f cycles/(sample*M)";

#if 0
const char* prf_samplecycle = "%d (%.1f samples/cycle)";
const char* prf_cyclessample = "%d (%.1f cycles/sample)";
#endif

/* ------------------------------------------------------------------------ */
