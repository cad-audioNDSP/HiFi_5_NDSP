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
 * Test procedures for image processing APIs.
 */

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API. */
#include LIBRARY_HEADER(img)
/* Test engine API. */
#include "testeng.h"
#include <stdlib.h>
#include <string.h>

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

#define MAX_FUNC_NUM   10
/* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }
/* Initializer for a test description structure. */
#define TEST_DESC( fmt, dimNum, align, loadFxn, procFxn, extraPtr ) { (fmt),0,extraPtr,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }

/* function copies image to output and fills real output frame with random data */
static void copyImgRandOutput(void* outImg,void* inImg, const imgsize_t* sz,size_t sz_element)
{
    int32_t r=233;
          uint8_t* dst=(      uint8_t*)outImg;
    const uint8_t* src=(const uint8_t*)inImg;
    int w,h,stride,n,m;
    stride=sz->stride*sz_element;
    w     =sz->width *sz_element;
    h     =sz->height;
    for (n=0; n<h; n++)
    {
        for (m=0; m<w; m++) 
        {
            r=(int32_t)(433494437*r+28657);
            *dst++=(uint8_t)r;
        }
        src+=w;
        for (;m<stride; m++) *dst++=*src++;
    }
}

/* Allocate vectors and load the data set for interleave*/
static int te_loadFxn_interleave( tTestEngContext * context )
{
    int res=0;
    int width,height,stride,N;
    ASSERT( context && context->seqFile );

    width  = MAX( 0, context->args.dim[0] );
    height = MAX( 0, context->args.dim[1] );
    stride = MAX( 0, context->args.dim[2] );
    (void)width;
    N=height*stride;

    /* Allocate output/reference data vectors memory. */
    if ( 6 != vecsAlloc( context->desc->isAligned, context->desc->fmt,
                        &context->dataSet.X  , N,
                        &context->dataSet.Y  , N,
                        &context->dataSet.U  , N,
                        &context->dataSet.Z  , 3*N,
                        &context->dataSet.Zlo, 3*N,
                        &context->dataSet.Zhi, 3*N,
                        0 ) )
    {
        printf( "te_loadFxn_interleave(): failed to allocate image data\n");
    }
    /* Load vectors data from the SEQ-file. */
    else if ( !seqFileReadVecs( context->seqFile,
                                &context->dataSet.X,
                                &context->dataSet.Y,
                                &context->dataSet.U,
                                &context->dataSet.Zlo,
                                &context->dataSet.Zhi, 
                                0 ) )
    {
        printf( "te_loadFxn_interleave(): failed to read image data\n");
        te_freeVectors(context);
    }
    else
    {
        res = 1;
        imgsize_t sz;
        sz.width =width *3;
        sz.height=height;
        sz.stride=stride*3;
        copyImgRandOutput(vecGetElem(&context->dataSet.Z,0),vecGetElem(&context->dataSet.Zlo,0),&sz,context->dataSet.Zlo.szElem);
    }
    return (res);
} 

/* test interleave functions */
static void te_processFxn_interleave( tTestEngContext * context )
{
    typedef void tFxn(      void * restrict outImg, 
                      const void * restrict inImgR, 
                      const void * restrict inImgG, 
                      const void * restrict inImgB, 
                      const imgsize_t* sz);
    imgsize_t sz;
    tTestEngTarget   fxn;

    sz.width  = MAX( 0, context->args.dim[0] );
    sz.height = MAX( 0, context->args.dim[1] );
    sz.stride = MAX( 0, context->args.dim[2] );

    ASSERT( context && context->target.fut );
    te_vReportStd(context);

    fxn = context->target.fut;

    ( (tFxn *)fxn )( vecGetElem( &context->dataSet.Z, 0 ), 
                     vecGetElem( &context->dataSet.X, 0 ), 
                     vecGetElem( &context->dataSet.Y, 0 ), 
                     vecGetElem( &context->dataSet.U, 0 ), 
                     &sz);
} /* te_processFxn_interleave() */

/* Allocate vectors and load the data set for convert*/
static int te_loadFxn_convert( tTestEngContext * context )
{
    int res=0,Nalloc;
    int width,height,stride,N;
    ASSERT( context && context->seqFile );

    width  = MAX( 0, context->args.dim[0] );
    height = MAX( 0, context->args.dim[1] );
    stride = MAX( 0, context->args.dim[2] );
    (void)width;
    N=height*stride;

    /* Allocate output/reference data vectors memory. */
    Nalloc=0;
    Nalloc += vecsAlloc( 0, FMT_REAL|FMT_INT32,
                        &context->dataSet.auxVec[0]  ,13,
                        0 );
    Nalloc += vecsAlloc( context->desc->isAligned, context->desc->fmt,
                        &context->dataSet.X  , N,
                        &context->dataSet.Y  , N,
                        &context->dataSet.U  , N,
                        &context->dataSet.Z  , 3*N,
                        &context->dataSet.Zlo, 3*N,
                        &context->dataSet.Zhi, 3*N,
                        0 );
    if (Nalloc!=7)
    {
        printf( "te_loadFxn_convert(): failed to allocate image data\n");
    }
    /* Load vectors data from the SEQ-file. */
    else if ( !seqFileReadVecs( context->seqFile,
                                &context->dataSet.auxVec[0],
                                &context->dataSet.X,
                                &context->dataSet.Y,
                                &context->dataSet.U,
                                &context->dataSet.Zlo,
                                &context->dataSet.Zhi, 
                                0 ) )
    {
        printf( "te_loadFxn_convert(): failed to read image data\n");
        te_freeVectors(context);
    }
    else
    {
        int off,szElem=context->dataSet.Zlo.szElem;
        imgsize_t sz;
        uint8_t *R,*G,*B;
        sz.width =width ;
        sz.height=height;
        sz.stride=stride;
        R=(uint8_t*)vecGetElem( &context->dataSet.Z, 0 );
        G=(uint8_t*)vecGetElem( &context->dataSet.Z, 1*sz.height*sz.stride );
        B=(uint8_t*)vecGetElem( &context->dataSet.Z, 2*sz.height*sz.stride );
        off=(uint8_t*)vecGetElem( &context->dataSet.Zlo, 0 )-R;
        copyImgRandOutput(R,R+off,&sz,szElem);
        copyImgRandOutput(G,G+off,&sz,szElem);
        copyImgRandOutput(B,B+off,&sz,szElem);
        res = 1;
    }
    return (res);
} 
/* test convert functions */
static void te_processFxn_convert( tTestEngContext * context )
{
    typedef void tFxn(      void * restrict outImgY, 
                            void * restrict outImgU, 
                            void * restrict outImgV, 
                      const void * restrict inImgR, 
                      const void * restrict inImgG, 
                      const void * restrict inImgB, 
                      const int32_t * coef,
                      const imgsize_t* sz);
    imgsize_t sz;
    tTestEngTarget   fxn;

    sz.width  = MAX( 0, context->args.dim[0] );
    sz.height = MAX( 0, context->args.dim[1] );
    sz.stride = MAX( 0, context->args.dim[2] );

    ASSERT( context && context->target.fut );
    te_vReportStd(context);

    fxn = context->target.fut;

    ( (tFxn *)fxn )( vecGetElem( &context->dataSet.Z,0),
                     vecGetElem( &context->dataSet.Z,1*sz.stride*sz.height),
                     vecGetElem( &context->dataSet.Z,2*sz.stride*sz.height),
                     vecGetElem( &context->dataSet.X, 0 ), 
                     vecGetElem( &context->dataSet.Y, 0 ), 
                     vecGetElem( &context->dataSet.U, 0 ), 
                     vecGetElem_i32( &context->dataSet.auxVec[0], 0 ),
                     &sz);
} /* te_processFxn_convert() */

/* Allocate vectors and load the data set for deinterleave*/
static int te_loadFxn_deinterleave( tTestEngContext * context )
{
    int res=0;
    int width,height,stride,N;
    ASSERT( context && context->seqFile );

    width  = MAX( 0, context->args.dim[0] );
    height = MAX( 0, context->args.dim[1] );
    stride = MAX( 0, context->args.dim[2] );
    (void)width;
    N=height*stride;

    /* Allocate output/reference data vectors memory. */
    if ( 4 != vecsAlloc( context->desc->isAligned, context->desc->fmt,
                        &context->dataSet.X  , 3*N,
                        &context->dataSet.Z  , 3*N,
                        &context->dataSet.Zlo, 3*N,
                        &context->dataSet.Zhi, 3*N,
                        0 ) )
    {
        printf( "te_loadFxn_deinterleave(): failed to allocate image data\n");
    }
    /* Load vectors data from the SEQ-file. */
    else if ( !seqFileReadVecs( context->seqFile,
                                &context->dataSet.X,
                                &context->dataSet.Zlo,
                                &context->dataSet.Zhi, 
                                0 ) )
    {
        printf( "te_loadFxn_deinterleave(): failed to read image data\n");
        te_freeVectors(context);
    }
    else
    {
        int off,szElem=context->dataSet.Zlo.szElem;
        imgsize_t sz;
        uint8_t *R,*G,*B;
        sz.width =width ;
        sz.height=height;
        sz.stride=stride;
        R=(uint8_t*)vecGetElem( &context->dataSet.Z, 0 );
        G=(uint8_t*)vecGetElem( &context->dataSet.Z, 1*sz.height*sz.stride);
        B=(uint8_t*)vecGetElem( &context->dataSet.Z, 2*sz.height*sz.stride);;
        off=(uint8_t*)vecGetElem( &context->dataSet.Zlo, 0 )-R;
        copyImgRandOutput(R,R+off,&sz,szElem);
        copyImgRandOutput(G,G+off,&sz,szElem);
        copyImgRandOutput(B,B+off,&sz,szElem);
        res = 1;
    }
    return (res);
} 

/* test interleave functions */
static void te_processFxn_deinterleave( tTestEngContext * context )
{
    typedef void tFxn(        void * restrict outImgR, 
                              void * restrict outImgG, 
                              void * restrict outImgB, 
                        const void * restrict inImg, 
                        const imgsize_t* sz);;
    imgsize_t sz;
    tTestEngTarget   fxn;
    sz.width  = MAX( 0, context->args.dim[0] );
    sz.height = MAX( 0, context->args.dim[1] );
    sz.stride = MAX( 0, context->args.dim[2] );

    ASSERT( context && context->target.fut );
    te_vReportStd(context);

    fxn = context->target.fut;

    ( (tFxn *)fxn )( 
            vecGetElem( &context->dataSet.Z, 0 ),
            vecGetElem( &context->dataSet.Z, 1*(sz.height*sz.stride)),
            vecGetElem( &context->dataSet.Z, 2*(sz.height*sz.stride)),
            vecGetElem( &context->dataSet.X, 0 ), &sz);
} /* te_processFxn_interleave() */


/* Allocate vectors and load the data set for histogram:
 * vector X (image), vector Z (histogram) */
static int te_loadFxn_hist( tTestEngContext * context )
{
    int32_t *ptr;
    int res,M,format,width,height,stride;
    int nAlloc;
    ASSERT( context && context->seqFile );

    format = MAX( 0, context->args.dim[0] );
    M      = MAX( 0, context->args.dim[1] );
    (void)format ;

    res=0;
    memset( &context->dataSet, 0, sizeof(context->dataSet) );
    nAlloc=vecsAlloc( context->desc->isAligned, FMT_INT32,
                        &context->dataSet.U  , 3, 
                        0) ;
    if (nAlloc!=1)
    {
        printf( "te_loadFxn_hist(): failed to allocate extra data\n");
        te_freeVectors(context);
        return 0;
    }
    if ( !seqFileReadVecs( context->seqFile,&context->dataSet.U,0))
    {
        printf( "te_loadFxn_hist(): failed to read extra data\n");
        te_freeVectors(context);
        return 0;
    }
    ptr=vecGetElem_i32(&context->dataSet.U,0);
    width =ptr[0];
    height=ptr[1];
    stride=ptr[2];
    (void)width;

    /* Allocate input data vector X. */
    if ( !vecAlloc( &context->dataSet.X, height*stride, context->desc->isAligned, context->desc->fmt, 0 ) )
    {
        printf( "te_loadFxn_hist(): failed to allocate image data\n");
    }
    /* Allocate output/reference data vectors memory. */
    else if ( 6 != vecsAlloc( 0, FMT_INT32,
                        &context->dataSet.Z  , M,
                        &context->dataSet.Zlo, M,
                        &context->dataSet.Zhi, M,
                        &context->dataSet.W  , 3,
                        &context->dataSet.Wlo, 3,
                        &context->dataSet.Whi, 3,
                        0 ) )
    {
        printf( "te_loadFxn_hist(): failed to allocate histogram\n");
    }
    /* Load vectors data from the SEQ-file. */
    else if ( !seqFileReadVecs( context->seqFile,
                                &context->dataSet.X,
                                &context->dataSet.Zlo,
                                &context->dataSet.Zhi, 
                                &context->dataSet.Wlo,
                                &context->dataSet.Whi, 
                                0 ) )
    {
        printf( "te_loadFxn_hist(): failed to read image data\n");
        te_freeVectors(context);
    }
    else
    {
        res = 1;
    }
    return (res);
} 

/* test histogram functions */
static void te_processFxn_hist( tTestEngContext * context )
{
    typedef void tFxn_gs  ( imghist_t * h, const void * inImg, const imgsize_t* sz, int M);
    int M;
    int32_t *ptr;
    imgsize_t sz;
    imghist_t histres;
    
    tTestEngTarget   fxn;
    void *X;
    int32_t *Z,*W;

    M         = MAX( 0, context->args.dim[1] );
    ptr=vecGetElem_i32(&context->dataSet.U,0);
    sz.width =ptr[0];
    sz.height=ptr[1];
    sz.stride=ptr[2];

    ASSERT( context && context->target.fut );
    te_vReportStd(context);

    X = vecGetElem( &context->dataSet.X, 0 );
    Z = (int32_t*)vecGetElem( &context->dataSet.Z, 0 );
    W = vecGetElem_i32( &context->dataSet.W, 0 );

    fxn = context->target.fut;
    histres.h=Z;

    ( (tFxn_gs *)fxn )( &histres, (const void*)X, &sz, M );
    W[0]=histres.mean;
    W[1]=histres.var ;
    W[2]=histres.M   ;
} /* te_processFxn_hist() */

/* Allocate vectors and load the data set for image normalization:
 * vector X (image), vector Z (histogram) */
static int te_loadFxn_norm( tTestEngContext * context )
{
    int nArgs,res,width,height,stride,sz;

    ASSERT( context && context->seqFile );

    width  = MAX( 0, context->args.dim[0] );
    height = MAX( 0, context->args.dim[1] );
    stride = MAX( 0, context->args.dim[2] );
    (void)width;

    res=0;
    memset( &context->dataSet, 0, sizeof(context->dataSet) );
    sz=height*stride;
    /* Allocate input data vector X. */
    nArgs=vecsAlloc( 1, FMT_INT16,
                        &context->dataSet.Y  , 2, 
                        &context->dataSet.U  , 64,
                        &context->dataSet.V,   3, 0 );
    nArgs+=vecsAlloc( context->desc->isAligned, context->desc->fmt,
                        &context->dataSet.X  , sz, 
                        &context->dataSet.Z  , sz,
                        &context->dataSet.Zlo, sz,
                        &context->dataSet.Zhi, sz, 0 );
    /* Allocate output/reference data vectors memory. */
    if ( 7 != nArgs )
    {
        printf( "te_loadFxn_norm(): failed to allocate images\n");
    }
    /* Load vectors data from the SEQ-file. */
    else if ( !seqFileReadVecs( context->seqFile,
                                &context->dataSet.Y,
                                &context->dataSet.U,
                                &context->dataSet.V,
                                &context->dataSet.X,
                                &context->dataSet.Zlo,
                                &context->dataSet.Zhi, 0 ) )
    {
        printf( "te_loadFxn_norm(): failed to read image data\n");
        te_freeVectors(context);
    }
    else
    {
        imgsize_t sz;
        sz.width  = width  ;
        sz.height = height ;
        sz.stride = stride ;
        copyImgRandOutput(vecGetElem(&context->dataSet.Z,0),vecGetElem(&context->dataSet.Zlo,0),&sz,context->dataSet.Zlo.szElem);
        res = 1;
    }
    return (res);
} 

/* test normalization functions */
static void te_processFxn_norm( tTestEngContext * context )
{
    typedef void tFxn_gs          ( void *  outImg, const void * inImg, const imgsize_t* sz, int minInt, int maxInt);
    typedef void tFxn_gs_nonlinear( void *  outImg, const void * inImg, const imgsize_t* sz, const int16_t* tbl);
    int16_t *ptr;
    int16_t *tbl;
    int minInt,maxInt,fast,method,format;
    imgsize_t sz;
    
    tTestEngTarget   fxn;
    void *X, *Z;
    (void)format,(void)fast;

    sz.width  = MAX( 0, context->args.dim[0] );
    sz.height = MAX( 0, context->args.dim[1] );
    sz.stride = MAX( 0, context->args.dim[2] );
    ptr=(int16_t*)vecGetElem( &context->dataSet.Y,0);
    minInt  = ptr[0];
    maxInt  = ptr[1];
    tbl=(int16_t*)vecGetElem( &context->dataSet.U,0);
    ptr=(int16_t*)vecGetElem( &context->dataSet.V,0);
    fast  =ptr[0];
    method=ptr[1];
    format=ptr[2];

    ASSERT( context && context->target.fut );
    te_vReportStd(context);

    X = vecGetElem( &context->dataSet.X, 0 );
    Z = vecGetElem( &context->dataSet.Z, 0 );

    fxn = context->target.fut;

    switch ( method )
    {
        case 0: ( (tFxn_gs *)fxn )( (void*)Z, (const void*)X, &sz, minInt,maxInt ); break;
        case 1: ( (tFxn_gs_nonlinear *)fxn )( (void*)Z, (const void*)X, &sz, tbl ); break;
        default: ASSERT( 0 );
    }
} /* te_processFxn_norm() */

typedef struct
{
    void   (*process)        ( void * restrict pScr, void * restrict outImg, const void * restrict inImg, const imgpad_params_t * params); /* resize images */ 
    size_t (*getScratchSize) ( const imgpad_params_t * params );/* Returns: size of scratch memory area, in bytes. */
}
tPadApi;

typedef struct
{
    imgpad_params_t params;
    int isFast;
    int format; // 0 - uint8, 1 - int16
}
tPadContext;

/* Allocate vectors and load the data set for image padding:*/
static int te_loadFxn_pad( tTestEngContext * context )
{
    int32_t *ptr;
    tPadContext padContext;
    tPadApi *pAPI;
    int res,szIn,szOut,scrSz,nAlloc;

    ASSERT( context && context->seqFile );
    pAPI=(tPadApi*)context->desc->extraPtr;
    padContext.params.in.width =context->args.dim[0];
    padContext.params.in.height=context->args.dim[1];
    padContext.params.in.stride=context->args.dim[2];

    res=0;
    memset( &context->dataSet, 0, sizeof(context->dataSet) );

    nAlloc=vecsAlloc( context->desc->isAligned, FMT_INT32,
                        &context->dataSet.X  , 3, 
                        &context->dataSet.U  , 3, 
                        &context->dataSet.Z  , 2,
                        0) ;
    if (nAlloc!=3)
    {
        printf( "te_loadFxn_pad(): failed to allocate extra data\n");
        te_freeVectors(context);
        return 0;
    }
    if ( !seqFileReadVecs( context->seqFile,&context->dataSet.X,&context->dataSet.U,&context->dataSet.Z,0))
    {
        printf( "te_loadFxn_pad(): failed to read extra data\n");
        te_freeVectors(context);
        return 0;
    }
    ptr=vecGetElem_i32(&context->dataSet.X,0);
    padContext.params.out.width =ptr[0];
    padContext.params.out.height=ptr[1];
    padContext.params.out.stride=ptr[2];
    ptr=vecGetElem_i32(&context->dataSet.U,0);
    padContext.isFast=ptr[0];
    padContext.format=ptr[1];
    padContext.params.fill=ptr[2];
    ptr=vecGetElem_i32(&context->dataSet.Z,0);
    padContext.params.x=ptr[0];
    padContext.params.y=ptr[1];
    te_freeVectors(context);

    szIn =padContext.params.in.height*padContext.params.in.stride;
    szOut=padContext.params.out.height*padContext.params.out.stride;
    scrSz=pAPI->getScratchSize(&padContext.params);

    /* Allocate output/reference data vectors memory. */
    nAlloc=vecsAlloc( context->desc->isAligned, context->desc->fmt,
                        &context->dataSet.X  , szIn, 
                        0 ) ;
    nAlloc+=vecsAlloc( context->desc->isAligned, context->desc->fmt,
                        &context->dataSet.Z  , szOut,
                        &context->dataSet.Zlo, szOut,
                        &context->dataSet.Zhi, szOut, 
                        0 ) ;
    nAlloc+=vecsAlloc( 1, FMT_UINT8,
                        &context->dataSet.U  , sizeof(tPadContext), 
                        &context->dataSet.V  , scrSz, 
                        NULL) ;
    if ( nAlloc!=6 )
    {
        printf( "te_loadFxn_pad(): failed to allocate images\n");
        te_freeVectors(context);
    }
    /* Load vectors data from the SEQ-file. */
    else if ( !seqFileReadVecs( context->seqFile,
                                &context->dataSet.X,
                                &context->dataSet.Zlo,
                                &context->dataSet.Zhi, 0 ) )
    {
        printf( "te_loadFxn_pad(): failed to read vectors data\n");
        te_freeVectors(context);
    }
    else
    {
        tPadContext* pPadContext;
        pPadContext=(tPadContext*)vecGetElem(&context->dataSet.U,0);
        pPadContext[0]=padContext;
        copyImgRandOutput(vecGetElem(&context->dataSet.Z,0),vecGetElem(&context->dataSet.Zlo,0),&padContext.params.out,context->dataSet.Zlo.szElem);
        res = 1;
    }
    return (res);
} 

/* test padding functions */
static void te_processFxn_pad( tTestEngContext * context )
{
    tPadApi *pAPI;
    void *outImg, *inImg, *pScr;
    tPadContext* pPadContext;
    pPadContext=(tPadContext*)vecGetElem(&context->dataSet.U,0);
    inImg  = vecGetElem( &context->dataSet.X, 0 );
    outImg = vecGetElem( &context->dataSet.Z, 0 );
    pScr   = vecGetElem( &context->dataSet.V, 0 );
    ASSERT( context && context->target.fut );
    te_vReportStd(context);
    pAPI=(tPadApi*)context->desc->extraPtr;
    pAPI->process(pScr,outImg,inImg,&pPadContext->params);
}

/* vec API test definitions. */
static const tPadApi imgpad_gu8_api      = {imgpad_gu8       ,imgpad_gu8_getScratchSize      };
static const tPadApi imgfastpad_gu8_api  = {imgfastpad_gu8   ,imgfastpad_gu8_getScratchSize  };
static const tPadApi imgpad_gs8_api      = {imgpad_gs8       ,imgpad_gs8_getScratchSize      };
static const tPadApi imgfastpad_gs8_api  = {imgfastpad_gs8   ,imgfastpad_gs8_getScratchSize  };
static const tPadApi imgpad_gs16_api    = {imgpad_gs16     ,imgpad_gs16_getScratchSize    };
static const tPadApi imgfastpad_gs16_api= {imgfastpad_gs16 ,imgfastpad_gs16_getScratchSize};

typedef struct  
{
  tTestEngTarget   funcList[MAX_FUNC_NUM];
  tTestEngDesc     testDesc;
}
deftbl_t;

/* Perform all tests for image processing APIs. */
int func_imgmisc1(int isFull, int isVerbose, int breakOnError)
{
    int res = 1;
    static const deftbl_t testDefTblmisc[] =
    {
    {FUNC_LIST( (tTestEngTarget)&imgconvert_rgbyuv,(tTestEngTarget)&imgconvert_yuvrgb ),         TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_convert, &te_processFxn_convert, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastconvert_rgbyuv,(tTestEngTarget)&imgfastconvert_yuvrgb ), TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_convert, &te_processFxn_convert, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgconvert_rgbyuv16,(tTestEngTarget)&imgconvert_yuvrgb16 ),         TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_convert, &te_processFxn_convert, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastconvert_rgbyuv16,(tTestEngTarget)&imgfastconvert_yuvrgb16 ), TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_convert, &te_processFxn_convert, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imginterleave ),       TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_interleave, &te_processFxn_interleave, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastinterleave ),   TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_interleave, &te_processFxn_interleave, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgdeinterleave ),     TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_deinterleave, &te_processFxn_deinterleave, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastdeinterleave ), TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_deinterleave, &te_processFxn_deinterleave, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imginterleave16 ),       TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_interleave, &te_processFxn_interleave, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastinterleave16 ),   TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_interleave, &te_processFxn_interleave, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgdeinterleave16 ),     TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_deinterleave, &te_processFxn_deinterleave, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastdeinterleave16 ), TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_deinterleave, &te_processFxn_deinterleave, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imghist_gu8 ),        TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_hist, &te_processFxn_hist, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgfasthist_gu8 ),    TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_hist, &te_processFxn_hist, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imghist_gs8 ),        TEST_DESC( FMT_REAL|FMT_INT8, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_hist, &te_processFxn_hist, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgfasthist_gs8 ),    TEST_DESC( FMT_REAL|FMT_INT8, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_hist, &te_processFxn_hist, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imghist_gs16),       TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_hist, &te_processFxn_hist, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgfasthist_gs16 ),  TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_hist, &te_processFxn_hist, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgnorm_gu8,(tTestEngTarget)&imgnorm_gu8_nonlinear ),            TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_norm, &te_processFxn_norm, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastnorm_gu8,(tTestEngTarget)&imgfastnorm_gu8_nonlinear ),    TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_norm, &te_processFxn_norm, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgnorm_gs8,(tTestEngTarget)&imgnorm_gs8_nonlinear ),            TEST_DESC( FMT_REAL|FMT_INT8, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_norm, &te_processFxn_norm, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastnorm_gs8,(tTestEngTarget)&imgfastnorm_gs8_nonlinear ),    TEST_DESC( FMT_REAL|FMT_INT8, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_norm, &te_processFxn_norm, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgnorm_gs16,(tTestEngTarget)&imgnorm_gs16_nonlinear ),        TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_norm, &te_processFxn_norm, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastnorm_gs16,(tTestEngTarget)&imgfastnorm_gs16_nonlinear ),TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_norm, &te_processFxn_norm, NULL ) },
    {FUNC_LIST( (tTestEngTarget)&imgpad_gu8 ),        TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_pad, &te_processFxn_pad, &imgpad_gu8_api ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastpad_gu8 ),    TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_pad, &te_processFxn_pad, &imgfastpad_gu8_api ) },  
    {FUNC_LIST( (tTestEngTarget)&imgpad_gs8 ),        TEST_DESC( FMT_REAL|FMT_INT8, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_pad, &te_processFxn_pad, &imgpad_gs8_api ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastpad_gs8 ),    TEST_DESC( FMT_REAL|FMT_INT8, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_pad, &te_processFxn_pad, &imgfastpad_gs8_api ) },  
    {FUNC_LIST( (tTestEngTarget)&imgpad_gs16 ),      TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_pad, &te_processFxn_pad, &imgpad_gs16_api ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastpad_gs16 ),  TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_pad, &te_processFxn_pad, &imgfastpad_gs16_api ) },  
    { FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, NULL, NULL, NULL ) } /* End of table */
    };
    #define DO_TEST(fxn, seqFile)                                                                       \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTblmisc,                                  \
                                                       sizeof(testDefTblmisc)/sizeof(testDefTblmisc[0]),\
                                                       MAX_FUNC_NUM,                                    \
                                                       (tTestEngTarget)(fxn), "imgmisc1/" seqFile,           \
                                                       isFull, isVerbose, breakOnError ) )

    DO_TEST( &imghist_gu8               , "imghist_gu8.seq"             );
    DO_TEST( &imgfasthist_gu8           , "imgfasthist_gu8.seq"         );
    DO_TEST( &imghist_gs8               , "imghist_gs8.seq"             );
    DO_TEST( &imgfasthist_gs8           , "imgfasthist_gs8.seq"         );
    DO_TEST( &imghist_gs16              , "imghist_gs16.seq"            );
    DO_TEST( &imgfasthist_gs16          , "imgfasthist_gs16.seq"        );

    DO_TEST( &imgnorm_gu8               , "imgnorm_gu8.seq"               );
    DO_TEST( &imgfastnorm_gu8           , "imgfastnorm_gu8.seq"           );
    DO_TEST( &imgnorm_gs8               , "imgnorm_gs8.seq"               );
    DO_TEST( &imgfastnorm_gs8           , "imgfastnorm_gs8.seq"           );
    DO_TEST( &imgnorm_gs16              , "imgnorm_gs16.seq"              );
    DO_TEST( &imgfastnorm_gs16          , "imgfastnorm_gs16.seq"          );
    DO_TEST( &imgnorm_gu8_nonlinear     , "imgnorm_gu8_nonlinear.seq"     );
    DO_TEST( &imgfastnorm_gu8_nonlinear , "imgfastnorm_gu8_nonlinear.seq" );
    DO_TEST( &imgnorm_gs8_nonlinear     , "imgnorm_gs8_nonlinear.seq"     );
    DO_TEST( &imgfastnorm_gs8_nonlinear , "imgfastnorm_gs8_nonlinear.seq" );
    DO_TEST( &imgnorm_gs16_nonlinear    , "imgnorm_gs16_nonlinear.seq"    );
    DO_TEST( &imgfastnorm_gs16_nonlinear, "imgfastnorm_gs16_nonlinear.seq");

    DO_TEST( &imginterleave             , "imginterleave.seq"           ) ;
    DO_TEST( &imgfastinterleave         , "imgfastinterleave.seq"       );
    DO_TEST( &imgdeinterleave           , "imgdeinterleave.seq"         );
    DO_TEST( &imgfastdeinterleave       , "imgfastdeinterleave.seq"     );
    DO_TEST( &imginterleave16           , "imginterleave16.seq"         );
    DO_TEST( &imgfastinterleave16       , "imgfastinterleave16.seq"     );
    DO_TEST( &imgdeinterleave16         , "imgdeinterleave16.seq"       );
    DO_TEST( &imgfastdeinterleave16     , "imgfastdeinterleave16.seq"   );

    DO_TEST( &imgconvert_rgbyuv         , "imgconvert_rgbyuv.seq"      );
    DO_TEST( &imgfastconvert_rgbyuv     , "imgfastconvert_rgbyuv.seq"  );
    DO_TEST( &imgconvert_yuvrgb         , "imgconvert_yuvrgb.seq"      );
    DO_TEST( &imgfastconvert_yuvrgb     , "imgfastconvert_yuvrgb.seq"  );
    DO_TEST( &imgconvert_rgbyuv16       , "imgconvert_rgbyuv16.seq"    );
    DO_TEST( &imgfastconvert_rgbyuv16   , "imgfastconvert_rgbyuv16.seq");
    DO_TEST( &imgconvert_yuvrgb16       , "imgconvert_yuvrgb16.seq"    );
    DO_TEST( &imgfastconvert_yuvrgb16   , "imgfastconvert_yuvrgb16.seq");

    DO_TEST( &imgpad_gu8     , "imgpad_gu8.seq"           );
    DO_TEST( &imgfastpad_gu8 , "imgfastpad_gu8.seq"       );
    DO_TEST( &imgpad_gs8     , "imgpad_gs8.seq"           );
    DO_TEST( &imgfastpad_gs8 , "imgfastpad_gs8.seq"       );
    DO_TEST( &imgpad_gs16    , "imgpad_gs16.seq"          );
    DO_TEST( &imgfastpad_gs16, "imgfastpad_gs16.seq"      );

    return (res);
    #undef DO_TEST
}
