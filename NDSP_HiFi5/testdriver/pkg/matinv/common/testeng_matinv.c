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
 * Test-engine add-on for matrix inversion categories 
 */

#include <string.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test engine API. */
#include "testeng_matinv.h"
#include "malloc16.h"

#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )
#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )


int te_loadFxn_mtxinv(tTestEngContext * context)
{
    int ret,N;
    N=context->args.dim[0];
    if(N>0) context->args.dim[0]=N*N; /* load and process N^2 elements */
    ret=te_loadFxn_vXvZ(context);
    context->args.dim[0]=N; /* get back N */
    return ret;
}
int te_loadFxn_mtxgjelim(tTestEngContext * context)
{
  int  N,fmt;
  int nElemX,nElemY,nElemZ, res;

  ASSERT( context && context->seqFile );

  N = MAX( 0, context->args.dim[0] );

  nElemX = N*N;
  nElemY = N;
  nElemZ = N;
  fmt = context->desc->fmt;
  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  /* Allocate data vectors memory. */
  res = ( 1 == vecsAlloc( context->desc->isAligned, fmt,&context->dataSet.X  , nElemX, 0 ) );
  res&= ( 1 == vecsAlloc( context->desc->isAligned, fmt,&context->dataSet.Y  , nElemY, 0 ) );
  res&= ( 3 == vecsAlloc( context->desc->isAligned, fmt,&context->dataSet.Z  , nElemZ ,
                         &context->dataSet.Zlo, nElemZ,
                         &context->dataSet.Zhi, nElemZ, 0 ) );

  if ( res )
  {
    /* Load vectors data from the SEQ-file. */
    if ( !( res = seqFileReadVecs( context->seqFile,
                                  &context->dataSet.X,
                                  &context->dataSet.Y,
                                  &context->dataSet.Zlo,
                                  &context->dataSet.Zhi, 0 ) ) )
    {
      printf( "te_loadFxn_mtxgjelim(): failed to read vectors data; "
              "fmt = 0x%02x, nElemX = %d, nElemY = %d, nElemZ = %d\n",
              (unsigned)context->desc->fmt, nElemX, nElemY, nElemZ );
    }
  }
  else
  {
    printf( "te_loadFxn_mtxgjelim(): failed to allocate vectors; "
            "fmt = 0x%02x, nElemX = %d, nElemY = %d, nElemZ = %d\n",
            (unsigned)context->desc->fmt, nElemX, nElemY, nElemZ );
  }

  /* Free vectors data if failed. */
  if ( !res ) te_freeVectors(context);

  return (res);
}
#include <math.h>
// return representation
static int inplaceCopyFloatFract(float32_t *f, int N)
{
    int qA;
    int32_t * x=(int32_t * )f;
    int n;
    float32_t maxf;
    for (maxf=0.f, n=0; n<N; n++) maxf=MAX(maxf,fabsf(f[n]));
    qA=(int)floorf(log2f(maxf)); 
    qA=31-qA-1;
    for (n=0; n<N; n++)
    {
        int64_t t;
        t=ldexpf(f[n],qA);
        t=MAX(MIN_INT32,MIN(MAX_INT32,t));
        x[n]=(int32_t)t;
    }
    return qA;
}

static void inplaceCopyFractFloat(int32_t *x, int qA, int N)
{
    int n;
    float32_t *f=(float32_t *)x;
    for (n=0; n<N; n++)
    {
        f[n]=(float32_t)ldexp((double)x[n],-qA);
    }
}


/* Apply the target function to the test case data set:
 * vector X (in), vector Z (out) */
void te_processFxn_matinv( tTestEngContext * context )
{
    int N;
    typedef void tFxnFloat      (void * pScr, const void * X);
    typedef int  tFxnFixed      (void * pScr, const void * X,int qA);
    const tMtxInvApi* pApi=(const tMtxInvApi*)context->desc->extraPtr;
    void *X, *Z;
    void *pScr;
    int isFloatingPoint=pApi->isFloatingPoint;

    ASSERT( context && context->target.fut );
    X = vecGetElem( &context->dataSet.X, 0 );
    Z = vecGetElem( &context->dataSet.Z, 0 );
    N=context->args.dim[0];
    memcpy(Z,X,vecGetSize(&context->dataSet.X));
    te_vReportStd(context);
    pScr = mallocAlign(pApi->getScratchSize(),0);
    if (pScr==NULL )  
    {
        printf( "te_processFxn_matinv(): failed to allocate memory\n");
        return;
    }

    switch(context->desc->extraParam)
    {
    case MTXINV_PLAIN:    
        if (isFloatingPoint)
        {
            ((tFxnFloat*)context->target.fut)( pScr, Z ); 
        }
        else
        {
            // fixed point API: compute qA, make in-place conversion and inverse operation on the output
            float32_t *f=(float32_t *)Z;
            int32_t   *x=(int32_t   *)Z;
            int cplx=(context->desc->fmt & FMT_CPLX)?2:1;
            int qA,qB;
            qA=inplaceCopyFloatFract(f,N*N*cplx);
            qB=((tFxnFixed*)context->target.fut)( pScr, x,qA ); 
            inplaceCopyFractFloat(x,qB,N*N*cplx);
        }
        break;
    default:              ASSERT(0);
    }
    freeAlign(pScr);
} /* te_processFxn_matinv() */


/* Apply the target function to the test case data set for Gauss-Jordan elimination */
void te_processFxn_gjelim( tTestEngContext * context )
{
    int N;
    typedef void tFxnFloat      (void* pScr, float32_t *y, const float32_t* A,const float32_t * x);
    typedef int  tFxnFixed      (void* pScr, int32_t *y, const int32_t* A,const int32_t * x, int qA, int qX);
    const tMtxInvApi* pApi=(const tMtxInvApi*)context->desc->extraPtr;
    void *X, *Z, *Y;
    void *pScr;
    int isFloatingPoint=pApi->isFloatingPoint;

    ASSERT( context && context->target.fut );
    X = vecGetElem( &context->dataSet.X, 0 );
    Y = vecGetElem( &context->dataSet.Y, 0 );
    Z = vecGetElem( &context->dataSet.Z, 0 );
    N=context->args.dim[0];
    te_vReportStd(context);
    pScr = mallocAlign(pApi->getScratchSize(),0);
    if (pScr==NULL )  
    {
        printf( "te_processFxn_matinv(): failed to allocate memory\n");
        return;
    }

    switch(context->desc->extraParam)
    {
    case MTXINV_PLAIN:    
        if (isFloatingPoint)
        {
            float32_t *fA=(float32_t *)X;
            float32_t *fX=(float32_t *)Y;
            float32_t *fY=(float32_t *)Z;
           ((tFxnFloat*)context->target.fut)( pScr, fY, fA, fX);
        }
        else
        {
            // fixed point API: compute qA, make in-place conversion and inverse operation on the output
            float32_t *fA=(float32_t *)X;
            int32_t   *iA=(int32_t   *)X;
            float32_t *fX=(float32_t *)Y;
            int32_t   *iX=(int32_t   *)Y;
            int32_t   *iY=(int32_t   *)Z;
            int cplx=(context->desc->fmt & FMT_CPLX)?2:1;
            int qA,qX,qY;
            qA=inplaceCopyFloatFract(fA,N*N*cplx);
            qX=inplaceCopyFloatFract(fX,N*cplx);
            qY=((tFxnFixed*)context->target.fut)( pScr, iY, iA, iX, qA,qX ); 
            inplaceCopyFractFloat(iY,qY,N*cplx);
        }
        break;
    default:              ASSERT(0);
    }
    freeAlign(pScr);
} /* te_processFxn_matinv() */


/* Perform all tests for matrix inversion API functions. */
int te_ExecMtxInv(const tTestMtxinv* pTbl, int isFull, int isVerbose, int breakOnError )
{
    int res = 1;
    for (;pTbl->pFunDescr;pTbl++)
    {
        if (isFull | pTbl->isFull)
        {
            res &= (0!=TestEngRun(pTbl->fxns, pTbl->pFunDescr, pTbl->seqFile, isFull, isVerbose, breakOnError,0));
            if (res == 0 && breakOnError) break;
        }
    }
    return (res);
} /* te_ExecMtxInv() */

