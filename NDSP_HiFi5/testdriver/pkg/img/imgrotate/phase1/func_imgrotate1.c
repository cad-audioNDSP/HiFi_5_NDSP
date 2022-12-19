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

typedef struct
{
    size_t (*alloc)          ( const imgrotate_params_t * params );   
    imgrotate_handle_t 
           (*init)           ( void * objmem, const imgrotate_params_t * params ); /* Returns: handle to the object, or NULL if initialization failed. */
    void   (*getOutSize)     (imgsize_t * restrict outSz, const imgrotate_params_t * params); /* return output image size */
    void   (*process)        ( imgrotate_handle_t handle, void * restrict pScr, void * restrict outImg, const void * restrict inImg, const imgsize_t* outSz); /* rotate images */ 
    size_t (*getScratchSize) ( const imgrotate_params_t * params );/* Returns: size of scratch memory area, in bytes. */
}
tRotateApi;

typedef struct
{
    imgrotate_params_t params;
    imgsize_t outSz;
    int isFast;
    int format; // 0 - -8bit, 1 -16bit
    imgrotate_handle_t handle;
}
tRotateContext;

/* Allocate vectors and load the data set for image rotation:
 * vector X (image), vector Y(angle, fillColor, isFast), vector Z (rotated image) */
static int te_loadFxn_rotate( tTestEngContext * context )
{
    tRotateContext rotateContext;
    tRotateApi *pAPI;
    int res,sz,scrSz,inst_sz,nAlloc;
    int32_t * ptr;
    imgsize_t outSz;

    ASSERT( context && context->seqFile );
    pAPI=(tRotateApi*)context->desc->extraPtr;
    rotateContext.handle=NULL;
    rotateContext.params.in.width  = MAX( 0, context->args.dim[0] );
    rotateContext.params.in.height = MAX( 0, context->args.dim[1] );
    rotateContext.params.in.stride = MAX( 0, context->args.dim[2] );

    res=0;
    memset( &context->dataSet, 0, sizeof(context->dataSet) );
    sz=rotateContext.params.in.height*rotateContext.params.in.stride;

    nAlloc=vecsAlloc( context->desc->isAligned, FMT_INT32,
                        &context->dataSet.U  , 3, 
                        &context->dataSet.V  , 4,
                        0) ;
    if (nAlloc!=2)
    {
        printf( "te_loadFxn_rotate(): failed to allocate extra data\n");
        te_freeVectors(context);
        return 0;
    }
    if ( !seqFileReadVecs( context->seqFile,&context->dataSet.U,&context->dataSet.V,0))
    {
        printf( "te_loadFxn_rotate(): failed to read extra data\n");
        te_freeVectors(context);
        return 0;
    }
    ptr=vecGetElem_i32(&context->dataSet.U,0);
    rotateContext.outSz.width =ptr[0];
    rotateContext.outSz.height=ptr[1];
    rotateContext.outSz.stride=ptr[2];
    ptr=vecGetElem_i32(&context->dataSet.V,0);
    rotateContext.params.angleQ15 =ptr[0];
    rotateContext.params.fill  =ptr[1];
    rotateContext.isFast       =ptr[2];
    rotateContext.format       =ptr[3];
    te_freeVectors(context);

    pAPI->getOutSize(&outSz,&rotateContext.params);
    NASSERT(outSz.width ==rotateContext.outSz.width);
    NASSERT(outSz.height==rotateContext.outSz.height);
    NASSERT(outSz.stride<=rotateContext.outSz.stride);
    outSz=rotateContext.outSz;
    scrSz=pAPI->getScratchSize(&rotateContext.params);
    inst_sz=pAPI->alloc(&rotateContext.params);

    /* Allocate output/reference data vectors memory. */
    nAlloc=vecsAlloc( context->desc->isAligned | ALLOC_DRAM0, context->desc->fmt,
                        &context->dataSet.X  , sz, 
                        &context->dataSet.Z  , outSz.stride*outSz.height,
                        &context->dataSet.Zlo, outSz.stride*outSz.height,
                        &context->dataSet.Zhi, outSz.stride*outSz.height, 0 ) ;
    nAlloc+=vecsAlloc( 1, FMT_UINT8,
                        &context->dataSet.U  , inst_sz+sizeof(tRotateContext), 
                        &context->dataSet.V  , scrSz, 
                        NULL) ;
    if ( nAlloc!=6 )
    {
        printf( "te_loadFxn_rotate(): failed to allocate images\n");
        te_freeVectors(context);
    }
    /* Load vectors data from the SEQ-file. */
    else if ( !seqFileReadVecs( context->seqFile,
                                &context->dataSet.X,
                                &context->dataSet.Zlo,
                                &context->dataSet.Zhi, 0 ) )
    {
        printf( "te_loadFxn_rotate(): failed to read vectors data\n");
        te_freeVectors(context);
    }
    else
    {
        tRotateContext* pRotateContext;
        pRotateContext=(tRotateContext*)vecGetElem(&context->dataSet.U,0);
        pRotateContext[0]=rotateContext;
        pRotateContext->handle=pAPI->init((void*)(pRotateContext+1),&rotateContext.params);
        copyImgRandOutput(vecGetElem(&context->dataSet.Z,0),vecGetElem(&context->dataSet.Zlo,0),&rotateContext.outSz,context->dataSet.Zlo.szElem);
        res = 1;
    }
    return (res);
} 

/* test rotation functions */
static void te_processFxn_rotate( tTestEngContext * context )
{
    tRotateApi *pAPI;
    void *outImg, *inImg, *pScr;
    tRotateContext* pRotateContext;
    pRotateContext=(tRotateContext*)vecGetElem(&context->dataSet.U,0);
    inImg  = vecGetElem( &context->dataSet.X, 0 );
    outImg = vecGetElem( &context->dataSet.Z, 0 );
    pScr   = vecGetElem( &context->dataSet.V, 0 );
    ASSERT( context && context->target.fut );
    te_vReportStd(context);
    pAPI=(tRotateApi*)context->desc->extraPtr;
    pAPI->process(pRotateContext->handle,pScr,outImg,inImg,&pRotateContext->outSz);
}

/* rotate API test definitions. */
static const tRotateApi imgrotate_gu8_api      ={imgrotate_gu8_alloc      ,imgrotate_gu8_init      ,imgrotate_gu8_getOutSize      ,imgrotate_gu8_process    ,imgrotate_gu8_getScratchSize      };
static const tRotateApi imgfastrotate_gu8_api  ={imgfastrotate_gu8_alloc  ,imgfastrotate_gu8_init  ,imgfastrotate_gu8_getOutSize  ,imgfastrotate_gu8_process,imgfastrotate_gu8_getScratchSize  };
static const tRotateApi imgrotate_gs8_api      ={imgrotate_gs8_alloc      ,imgrotate_gs8_init      ,imgrotate_gs8_getOutSize      ,imgrotate_gs8_process    ,imgrotate_gs8_getScratchSize      };
static const tRotateApi imgfastrotate_gs8_api  ={imgfastrotate_gs8_alloc  ,imgfastrotate_gs8_init  ,imgfastrotate_gs8_getOutSize  ,imgfastrotate_gs8_process,imgfastrotate_gs8_getScratchSize  };
static const tRotateApi imgrotate_gs16_api    ={imgrotate_gs16_alloc    ,imgrotate_gs16_init    ,imgrotate_gs16_getOutSize    ,imgrotate_gu8_process    ,imgrotate_gs16_getScratchSize    };
static const tRotateApi imgfastrotate_gs16_api={imgfastrotate_gs16_alloc,imgfastrotate_gs16_init,imgfastrotate_gs16_getOutSize,imgfastrotate_gu8_process,imgfastrotate_gs16_getScratchSize};

typedef struct  
{
  tTestEngTarget   funcList[MAX_FUNC_NUM];
  tTestEngDesc     testDesc;
}
deftbl_t;

int func_imgrotate1(int isFull, int isVerbose, int breakOnError)
{
    int res = 1;

    // rotation subset
    static const deftbl_t testDefTblrotate[] =
    {
    { FUNC_LIST( (tTestEngTarget)&imgrotate_gu8_process ),      TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_rotate, &te_processFxn_rotate, &imgrotate_gu8_api ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastrotate_gu8_process ),   TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_rotate, &te_processFxn_rotate, &imgfastrotate_gu8_api ) },  
    {FUNC_LIST( (tTestEngTarget)&imgrotate_gs8_process ),      TEST_DESC( FMT_REAL|FMT_INT8,  TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_rotate, &te_processFxn_rotate, &imgrotate_gs8_api ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastrotate_gs8_process ),  TEST_DESC( FMT_REAL|FMT_INT8,  TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_rotate, &te_processFxn_rotate, &imgfastrotate_gs8_api ) },  
    {FUNC_LIST( (tTestEngTarget)&imgrotate_gs16_process ),     TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_rotate, &te_processFxn_rotate, &imgrotate_gs16_api ) },
    {FUNC_LIST( (tTestEngTarget)&imgfastrotate_gs16_process ), TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_rotate, &te_processFxn_rotate, &imgfastrotate_gs16_api ) },  
    {FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, NULL, NULL, NULL ) } /* End of table */    };
    #define DO_TEST(fxn, seqFile)                                                                           \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTblrotate,                                    \
                                                       sizeof(testDefTblrotate)/sizeof(testDefTblrotate[0]),\
                                                       MAX_FUNC_NUM,                                        \
                                                       (tTestEngTarget)(fxn), "imgrotate1/" seqFile,               \
                                                       isFull, isVerbose, breakOnError ) )
    
	vecInitRegion(ALLOC_DRAM0);
    DO_TEST( &imgrotate_gu8_process     , "imgrotate_gu8.seq"           );
    DO_TEST( &imgfastrotate_gu8_process , "imgfastrotate_gu8.seq"       );
    DO_TEST( &imgrotate_gs8_process     , "imgrotate_gs8.seq"           );
    DO_TEST( &imgfastrotate_gs8_process , "imgfastrotate_gs8.seq"       );
    DO_TEST( &imgrotate_gs16_process    , "imgrotate_gs16.seq"          );
    DO_TEST( &imgfastrotate_gs16_process, "imgfastrotate_gs16.seq"      );
    #undef DO_TEST
    return (res);
}
