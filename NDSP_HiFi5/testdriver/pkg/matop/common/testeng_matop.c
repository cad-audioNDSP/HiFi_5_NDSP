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
 * Test-engine add-on for matrix categories 
 */

#include <string.h>
#include <stdlib.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* DSP Library API. */
#include LIBRARY_HEADER(matop)
/* Test engine API. */
#include "testeng_matop.h"
/* Aligning memory allocator. */
#include "malloc16.h"
/* Test data vectors tools and SEQ-file reader */
#include "vectools.h"

#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )


 /* Apply the function under test to test case data set: matrix transpose */
void te_processFxn_mtx_transpose( tTestEngContext * context )
{
  typedef void tFxn( void * z, const void * x, int M, int N);
  void *X, *Z;
  int M, N;

  ASSERT( context && context->target.fut );

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  M   = context->args.dim[0];
  N   = context->args.dim[1];

  te_vReportStd(context);
  ((tFxn*)context->target.fut)(Z,  X, M, N    );

} /* te_processFxn_mtx_transpose() */

/* Allocate vectors and load the data set for [c]matmmlt:
 * X[M][N], Y[N][P], Z[M][P] */
int te_loadFxn_mmlt( tTestEngContext * context )
{
    tVec X, Y, Z, Zlo, Zhi;
    int M, N, P;
    int nElemX, nElemY, nElemZ, res = 0;

    ASSERT( context && context->seqFile );

    memset( &X  , 0, sizeof(X  ) );
    memset( &Y  , 0, sizeof(Y  ) );
    memset( &Z  , 0, sizeof(Z  ) );
    memset( &Zlo, 0, sizeof(Zlo) );
    memset( &Zhi, 0, sizeof(Zhi) );

    M = MAX( 0, context->args.dim[0] );
    N = MAX( 0, context->args.dim[1] );
    P = MAX(0, context->args.dim[2]);

    nElemX = M*N;
    nElemY = N*P;
    nElemZ = M*P;

    /* Allocate data vectors memory. */
    res=(3 == vecsAlloc( context->desc->isAligned,
                            context->desc->fmt,
                            &X, nElemX,
                            &Y, nElemY,
                            &Z, nElemZ, 0 ));
    res&=(2 == vecsAlloc( TE_ALIGN_NO,
                            context->desc->fmt,
                            &Zlo, nElemZ,
                            &Zhi, nElemZ, 0 ));
    if (!res)
    {
    printf( "loadFxn_mmlt(): failed to allocate vectors; "
            "fmt = 0x%02x, nElemX = %d, nElemY = %d, nElemZ = %d\n",
            (unsigned)context->desc->fmt, nElemX, nElemY, nElemZ );
    }
    /* Load vectors data from the SEQ-file. */
    else if ( !seqFileReadVecs( context->seqFile, &X, &Y, &Zlo, &Zhi, 0 ) )
    {
    printf( "loadFxn_mmlt(): failed to read vectors data; "
            "fmt = 0x%02x, nElemX = %d, nElemY = %d, nElemZ = %d\n",
            (unsigned)context->desc->fmt, nElemX, nElemY, nElemZ );
    }
    else
    {
        memset( &context->dataSet, 0, sizeof( context->dataSet ) );

        context->dataSet.X   = X;
        context->dataSet.Y   = Y;
        context->dataSet.Z   = Z;
        context->dataSet.Zlo = Zlo;
        context->dataSet.Zhi = Zhi;
        res = 1;
    }

  if ( !res ) te_freeVectors(context); /* Free vectors data if failed. */
  return (res);

} /* loadFxn_mmlt() */

static void inplace16x8Convert(void* Src, int N)
{
    int n;
    const int16_t *x=(const int16_t *)Src;
          int8_t  *y=(      int8_t  *)Src;
    for (n=0; n<N; n++) y[n]=(int8_t)(x[n]>>8);
}

/* Apply the function under test to test case data set: [c]matmmlt */
void te_processFxn_mmlt( tTestEngContext * context )
{
  typedef void tFxn      (void *pScr,  void * z,   const void * x, const void * y, int M, int N, int P        );
  typedef void tFxnlshscr(void *pScr, void * z,   const void * x, const void * y, int M, int N, int P, int lsh);
  void *X, *Y, *Z;
  int M, N, P, lsh;
  size_t scrSz;
  void *pScr;

  ASSERT( context && context->target.fut );

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  M   = context->args.dim[0];
  N   = context->args.dim[1];
  P   = context->args.dim[2];
  lsh = context->args.dim[3];

  scrSz=0;
  if(context->desc->extraPtr ) 
  {
      const tMatmulApi_MNP* api=(const tMatmulApi_MNP*)context->desc->extraPtr;
      scrSz=api->getScratchSize(M,N,P);
  }
  pScr=scrSz==0 ? NULL: mallocAlign(scrSz,0);
  te_vReportStd(context);
  switch(context->desc->extraParam & 3)
  {
  case MTX_PLAIN:     ((tFxn*   )context->target.fut)(pScr,Z, X, Y, M, N, P      );
     if (pScr) freeAlign(pScr);
     break;
  case MTX_PLAIN_LSH:
    {   
        ((tFxnlshscr*)context->target.fut)(pScr, Z, X, Y, M, N, P, lsh ); 
    }
    if(pScr) freeAlign(pScr);
    break;
  case MTX_8X16_LSH:
    {   
        inplace16x8Convert(X,(M>0 && N>0) ? M*N:0);
        ((tFxnlshscr*)context->target.fut)(pScr, Z, X, Y, M, N, P, lsh ); 
    }
    if(pScr) freeAlign(pScr);
    break;
  }

} /* processFxn_mmlt() */

/* Apply the function under test to test case data set: vecmpy */
void te_processFxn_vecmmlt( tTestEngContext * context )
{
  typedef void tFxn   (void * z,   const void * x, const void * y, int M, int N         );
  typedef void tFxnlsh(void * z,   const void * x, const void * y, int M, int N, int lsh);
  void *X, *Y, *Z;
  int M, N, P, lsh;

  ASSERT( context && context->target.fut );

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  M   = context->args.dim[0];
  N   = context->args.dim[1];
  P   = context->args.dim[2];
  lsh   = context->args.dim[3];
  ASSERT(P==1); (void)P;
  te_vReportStd(context);
  switch(context->desc->extraParam)
  {
  case MTX_PLAIN:     ((tFxn*   )context->target.fut)(Z, X, Y, M, N     ); break;
  case MTX_PLAIN_LSH: ((tFxnlsh*)context->target.fut)(Z, X, Y, M, N, lsh); break;
  case MTX_8X16_LSH:  inplace16x8Convert(X,(M>0 && N>0) ? M*N:0);
                      ((tFxnlsh*)context->target.fut)(Z, X, Y, M, N, lsh); break;
  }

} /* te_processFxn_vecmmlt() */

/* Allocate vectors and load the data set for [c]mtx_vecmpyt: X[N], Z[N][N] */
int  te_loadFxn_vecmpyt( tTestEngContext * context )
{
    tVec X, Z, Zlo, Zhi;
    int N;
    int nElemX, nElemZ, res = 0;

    ASSERT( context && context->seqFile );

    memset( &X  , 0, sizeof(X  ) );
    memset( &Z  , 0, sizeof(Z  ) );
    memset( &Zlo, 0, sizeof(Zlo) );
    memset( &Zhi, 0, sizeof(Zhi) );

    N = MAX( 0, context->args.dim[0] );
    nElemX = N;
    nElemZ = N*N;

    /* Allocate data vectors memory. */
    res = (4 == vecsAlloc( context->desc->isAligned,
                           context->desc->fmt,
                           &X  , nElemX,
                           &Z  , nElemZ,
                           &Zlo, nElemZ,
                           &Zhi, nElemZ, 0 ));
    if (!res)
    {
        printf( "te_loadFxn_vecmpyt(): failed to allocate vectors; "
                "fmt = 0x%02x, nElemX = %d, nElemZ = %d\n",
                (unsigned)context->desc->fmt, nElemX, nElemZ );
    }
    /* Load vectors data from the SEQ-file. */
    else if ( !seqFileReadVecs( context->seqFile, &X, &Zlo, &Zhi, 0 ) )
    {
        printf( "te_loadFxn_vecmpyt(): failed to read vectors data; "
                "fmt = 0x%02x, nElemX = %d, nElemZ = %d\n",
                (unsigned)context->desc->fmt, nElemX, nElemZ );
    }
    else
    {
        memset( &context->dataSet, 0, sizeof( context->dataSet ) );
        context->dataSet.X   = X;
        context->dataSet.Z   = Z;
        context->dataSet.Zlo = Zlo;
        context->dataSet.Zhi = Zhi;
        res = 1;
    }

    if ( !res ) te_freeVectors(context); /* Free vectors data if failed. */
    return (res);

} /* te_loadFxn_vecmpyt() */

/* Apply the function under test to test case data set: [c]mtx_vecmpyt */
void te_processFxn_vecmpyt( tTestEngContext * context )
{
    typedef void tFxn   (void * z, const void * x, int N         );
    typedef void tFxnlsh(void * z, const void * x, int N, int lsh);
    void *X, *Z;
    int N, lsh;
  
    ASSERT( context && context->target.fut );
  
    X = vecGetElem( &context->dataSet.X, 0 );
    Z = vecGetElem( &context->dataSet.Z, 0 );
  
    N   = context->args.dim[0];
    lsh = context->args.dim[1];
    te_vReportStd(context);
    switch(context->desc->extraParam)
    {
    case MTX_PLAIN:     ((tFxn*   )context->target.fut)(Z, X, N     ); break;
    case MTX_PLAIN_LSH: ((tFxnlsh*)context->target.fut)(Z, X, N, lsh); break;
    }
} /* te_processFxn_vecmpyt() */

/* Allocate vectors and load the data set for [c]mtx_lrmpy: X[N][N], Y[N][N], Z[N][N] */
int te_loadFxn_lrmpy( tTestEngContext * context )
{
    tVec X, Y, Z, Zlo, Zhi;
    int N;
    int nElemXYZ, res = 0;

    ASSERT( context && context->seqFile );

    memset( &X  , 0, sizeof(X  ) );
    memset( &Y  , 0, sizeof(Y  ) );
    memset( &Z  , 0, sizeof(Z  ) );
    memset( &Zlo, 0, sizeof(Zlo) );
    memset( &Zhi, 0, sizeof(Zhi) );

    N = MAX( 0, context->args.dim[0] );
    nElemXYZ = N*N;

    /* Allocate data vectors memory. */
    res = (5 == vecsAlloc( context->desc->isAligned,
                           context->desc->fmt,
                           &X  , nElemXYZ,
                           &Y  , nElemXYZ,
                           &Z  , nElemXYZ,
                           &Zlo, nElemXYZ,
                           &Zhi, nElemXYZ, 0 ));
    if (!res)
    {
        printf( "te_loadFxn_lrmpy(): failed to allocate vectors; "
                "fmt = 0x%02x, nElemXYZ = %d\n",
                (unsigned)context->desc->fmt, nElemXYZ );
    }
    /* Load vectors data from the SEQ-file. */
    else if ( !seqFileReadVecs( context->seqFile, &X, &Y, &Zlo, &Zhi, 0 ) )
    {
        printf( "te_loadFxn_lrmpy(): failed to read vectors data; "
                "fmt = 0x%02x, nElemXYZ = %d\n",
                (unsigned)context->desc->fmt, nElemXYZ );
    }
    else
    {
        memset( &context->dataSet, 0, sizeof( context->dataSet ) );
        context->dataSet.X   = X;
        context->dataSet.Y   = Y;
        context->dataSet.Z   = Z;
        context->dataSet.Zlo = Zlo;
        context->dataSet.Zhi = Zhi;
        res = 1;
    }

    if ( !res ) te_freeVectors(context); /* Free vectors data if failed. */
    return (res);

} /* te_loadFxn_lrmpy() */

/* Apply the function under test to test case data set: [c]mtx_lrmpy */
void te_processFxn_lrmpy( tTestEngContext * context )
{
    typedef void tFxn   (void * pScr, void * z, const void * x, const void * y, int N);
    typedef void tFxnlsh(void * pScr, void * z, const void * x, const void * y, int N, int rsh, int lsh);
    void *X, *Y, *Z, *pScr;
    size_t szScr;
    int N, rsh, lsh;
  
    ASSERT( context && context->target.fut );
  
    X = vecGetElem( &context->dataSet.X, 0 );
    Y = vecGetElem( &context->dataSet.Y, 0 );
    Z = vecGetElem( &context->dataSet.Z, 0 );
  
    N   = context->args.dim[0];
    rsh = context->args.dim[1];
    lsh = context->args.dim[2];

    szScr = 0;
    if (context->desc->extraPtr) 
    {
        const tMatmulApi_N* api=(const tMatmulApi_N*)context->desc->extraPtr;
        szScr = api->getScratchSize(N);
    }
    pScr = szScr==0 ? NULL : mallocAlign(szScr, 0);

    te_vReportStd(context);
    switch(context->desc->extraParam)
    {
    case MTX_PLAIN:         ((tFxn*   )context->target.fut)(pScr, Z, X, Y, N          ); break;
    case MTX_PLAIN_RSH_LSH: ((tFxnlsh*)context->target.fut)(pScr, Z, X, Y, N, rsh, lsh); break;
    }

    if (pScr) freeAlign(pScr);

} /* te_processFxn_lrmpy() */

