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
* Test-engine add-on for arithmetic and logic functions on data vectors
*/
#include <string.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test engine API. */
#include "testeng_vector.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

/****************************/
/* Apply the target function to the test case data set:
* vector X (in), scalar F0 (in), scalar F1 (in), scalar Y (in), vector Z (out) */
void te_vector_processFxn_vXsF0sF1sYvZ(tTestEngContext * context)
{
	typedef void tFxn_fr16(const fract16         * x, fract16         y, fract16         * z, int N);
	typedef void tFxn_fr32(const fract32         * x, fract32         y, fract32         * z, int N);
	typedef void tFxn_fl32(float32_t       * z, const float32_t       * x, float32_t       y, float32_t       f0, float32_t       f1, int N);
	typedef void tFxn_fl64(const float64_t       * x, float64_t       y, float64_t       * z, int N);
	typedef void tFxn_fr16c(const complex_fract16 * x, complex_fract16 y, complex_fract16 * z, int N);
	typedef void tFxn_fr32c(const complex_fract32 * x, complex_fract32 y, complex_fract32 * z, int N);
	typedef void tFxn_fl32c(const complex_float   * x, complex_float   y, complex_float   * z, int N);
	typedef void tFxn_fl64c(const complex_double  * x, complex_double  y, complex_double  * z, int N);

	tTestEngTarget   fxn;
	void *X, *Y, *Z, *F0, *F1;
	int N;

	ASSERT(context && context->target.fut);
	te_vReportStd(context);
	X = vecGetElem(&context->dataSet.X, 0);
	Y = vecGetElem(&context->dataSet.Y, 0);
	F0 = vecGetElem(&context->dataSet.Wlo, 0);
	F1 = vecGetElem(&context->dataSet.Whi, 0);
	Z = vecGetElem(&context->dataSet.Z, 0);


	N = context->args.dim[0];
	fxn = context->target.fut;

	switch (context->desc->fmt)
	{
	case FMT_REAL | FMT_FRACT16: ((tFxn_fr16 *)fxn)((const fract16        *)X, *(fract16        *)Y, (fract16         *)Z, N); break;
	case FMT_REAL | FMT_FRACT32: ((tFxn_fr32 *)fxn)((const fract32        *)X, *(fract32        *)Y, (fract32         *)Z, N); break;
	case FMT_REAL | FMT_FLOAT32: ((tFxn_fl32 *)fxn)((float32_t       *)Z, (const float32_t      *)X, *(float32_t      *)Y, *(float32_t      *)F0, *(float32_t      *)F1, N); break;
	case FMT_REAL | FMT_FLOAT64: ((tFxn_fl64 *)fxn)((const float64_t      *)X, *(float64_t      *)Y, (float64_t       *)Z, N); break;
	case FMT_CPLX | FMT_FRACT16: ((tFxn_fr16c*)fxn)((const complex_fract16*)X, *(complex_fract16*)Y, (complex_fract16 *)Z, N); break;
	case FMT_CPLX | FMT_FRACT32: ((tFxn_fr32c*)fxn)((const complex_fract32*)X, *(complex_fract32*)Y, (complex_fract32 *)Z, N); break;
	case FMT_CPLX | FMT_FLOAT32: ((tFxn_fl32c*)fxn)((const complex_float  *)X, *(complex_float  *)Y, (complex_float   *)Z, N); break;
	case FMT_CPLX | FMT_FLOAT64: ((tFxn_fl64c*)fxn)((const complex_double *)X, *(complex_double *)Y, (complex_double  *)Z, N); break;
	default: ASSERT(0);
	}

} /* te_vector_processFxn_vXsF0sF1sYvZ() */
/******************************/

/****************/
/* Allocate vectors and load the data set:
* vector X (in), scalar F0 (in), scalar F1 (in), scalar Y (in), vector Z (out) */
int te_vector_loadFxn_vXsF0sF1sYvZ(tTestEngContext * context)
{
	int M, N, L;
	int nElemX, nElemY, nElemZ, nElemW, res;

	ASSERT(context && context->seqFile);

	M = MAX(0, context->args.dim[1]);
	N = MAX(0, context->args.dim[0]);
	L = MAX(0, context->args.dim[2]);

	nElemX = M*N*L;
	nElemY = L;
	nElemW = L;
	nElemZ = M*N*L;

	memset(&context->dataSet, 0, sizeof(context->dataSet));

	/* Allocate data vectors memory. */
	res = (7 == vecsAlloc(context->desc->isAligned, context->desc->fmt,
		&context->dataSet.X, nElemX,
		&context->dataSet.Y, nElemY,
		&context->dataSet.Z, nElemZ,
		&context->dataSet.Wlo, nElemW,
		&context->dataSet.Whi, nElemW,
		&context->dataSet.Zlo, nElemZ,
		&context->dataSet.Zhi, nElemZ, 0));
	if (res)
	{
		/* Load vectors data from the SEQ-file. */
		if (!(res = seqFileReadVecs(context->seqFile,
			&context->dataSet.X,
			&context->dataSet.Y,
			&context->dataSet.Wlo,
			&context->dataSet.Whi,
			&context->dataSet.Zlo,
			&context->dataSet.Zhi, 0)))
		{
			printf("te_loadFxn_vXsF0sF1sYvZ(): failed to read vectors data; "
				"fmt = 0x%02x, nElemX = %d, nElemY = %d,  nElemW = %d, nElemZ = %d\n",
				(unsigned)context->desc->fmt, nElemX, nElemY, nElemW, nElemZ);
		}
	}
	else
	{
		printf("te_loadFxn_vXsF0sF1sYvZ(): failed to allocate vectors; "
			"fmt = 0x%02x, nElemX = %d, nElemY = %d, nElemW = %d, nElemZ = %d\n",
			(unsigned)context->desc->fmt, nElemX, nElemY, nElemW, nElemZ);
	}

	/* Free vectors data if failed. */
	if (!res) te_freeVectors(context);

	return (res);

} /* te_vector_loadFxn_vXsF0sF1sYvZ() */

void processFxn_scl_vXvZ32(tTestEngContext * context)
{
	typedef int32_t tFxn_fl32(float32_t);
	typedef int32_t tFxn_fl64(float64_t);
	typedef int32_t tFxn_fr16(fract16);
	typedef int32_t tFxn_fr32(fract32);

#define CALL_FXN( typeFxn, Fxn, typeXZ, X, Z ) \
		{ *(int32_t*)(Z) = ( (typeFxn*)(Fxn) )( *(typeXZ*)(X) ); \
        (X) = (typeXZ*)(X) + 1; (Z) = (typeXZ*)(Z) + 1; }

	tTestEngTarget  fxn;
	void *X, *Z;
	int n, N;

	X = vecGetElem(&context->dataSet.X, 0);
	Z = vecGetElem(&context->dataSet.Z, 0);

	fxn = context->target.fut;
	N = context->args.dim[0];
	te_vReportStd(context);
	switch (context->desc->fmt & FMT_DTYPE_MASK)
	{
	case FMT_FLOAT32: for (n = 0; n < N; n++) CALL_FXN(tFxn_fl32, fxn, float32_t, X, Z); break;
	case FMT_FLOAT64: for (n = 0; n < N; n++) CALL_FXN(tFxn_fl64, fxn, float64_t, X, Z); break;
	case FMT_FRACT16: for (n = 0; n < N; n++) CALL_FXN(tFxn_fr16, fxn, fract16, X, Z); break;
	case FMT_FRACT32: for (n = 0; n < N; n++) CALL_FXN(tFxn_fr32, fxn, fract32, X, Z); break;
	default: ASSERT(0);
	}

#undef CALL_FXN

} /* processFxn_scl_vXvZ32() */

