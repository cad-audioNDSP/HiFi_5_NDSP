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
 * floating point extra statistics
 */
#ifndef FPSTAT_H__
#define FPSTAT_H__
#include <stdio.h>
#include "fpstat.h"
#include "types.h"
#include "float16.h"

struct tagStatList;
typedef struct
{
    struct tagStatList* list;
}
tFPStat;

// open statistics
void FPStatCreate(tFPStat* pFPStat);

/*-----------------------------------------
add data to the statistics
Input:
patName    - name of pattern
refData[N] - reference data
data[N]    - data
N          - number of data points
-----------------------------------------*/
void FPStatAdd(tFPStat* pFPStat, const char* patName, const float64_t* refData, const float32_t* data, int N);
void FPStatAdd_fp16(tFPStat* pFPStat, const char* patName, const float64_t* refData, const float16_t* data, int N);

// print statistics to the file
void FPStatPrint(tFPStat* pFPStat, FILE* f);
// close statisctics
void FPStatClose (tFPStat* pFPStat);

#endif
