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
static void copyImgRandOutput(void* outImg, void* inImg, const imgsize_t* sz, size_t sz_element)
{

	int32_t r = 233;
	uint8_t* dst = (uint8_t*)outImg;
	const uint8_t* src = (const uint8_t*)inImg;
	int w, h, stride, n, m;
	stride = sz->stride*sz_element;
	w = sz->width *sz_element;
	h = sz->height;
	for (n = 0; n<h; n++)
	{
		for (m = 0; m<w; m++)
		{
			r = (int32_t)(433494437 * r + 28657);
			*dst++ = (uint8_t)r;
		}
		src += w;
		for (; m<stride; m++) *dst++ = *src++;
	}
}

typedef struct
{
    size_t (*alloc)          ( const imgresize_params_t * params );   
    imgresize_handle_t 
           (*init)           ( void * objmem, const imgresize_params_t * params ); /* Returns: handle to the object, or NULL if initialization failed. */
    void   (*process)        ( imgresize_handle_t handle, void * restrict pScr, void * restrict outImg, const void * restrict inImg); /* resize images */ 
    size_t (*getScratchSize) ( const imgresize_params_t * params );/* Returns: size of scratch memory area, in bytes. */
}
tResizeApi;

typedef struct
{
    imgresize_params_t params;
    int isFast;
    int format; // 0 - uint8, 1 - int16, 2 -int8
    imgresize_handle_t handle;
}
tResizeContext;

/* Allocate vectors and load the data set for image rotation:
 * vector X (image), vector Y(angle, fillColor, isFast), vector Z (resized image) */
static int te_loadFxn_resize( tTestEngContext * context )
{
    int32_t *ptr;
    tResizeContext resizeContext;
    tResizeApi *pAPI;
    int res,szIn,szOut,scrSz,inst_sz,nAlloc;

    ASSERT( context && context->seqFile );
    pAPI=(tResizeApi*)context->desc->extraPtr;
    resizeContext.handle=NULL;
    resizeContext.params.in.width =context->args.dim[0];
    resizeContext.params.in.height=context->args.dim[1];
    resizeContext.params.in.stride=context->args.dim[2];

    res=0;
    memset( &context->dataSet, 0, sizeof(context->dataSet) );

    nAlloc=vecsAlloc( context->desc->isAligned, FMT_INT32,
                        &context->dataSet.U  , 3, 
                        &context->dataSet.V  , 3,
                        0) ;
    if (nAlloc!=2)
    {
        printf( "te_loadFxn_resize(): failed to allocate extra data\n");
        te_freeVectors(context);
        return 0;
    }
    if ( !seqFileReadVecs( context->seqFile,&context->dataSet.U,&context->dataSet.V,0))
    {
        printf( "te_loadFxn_resize(): failed to read extra data\n");
        te_freeVectors(context);
        return 0;
    }
    ptr=vecGetElem_i32(&context->dataSet.U,0);
    resizeContext.params.out.width =ptr[0];
    resizeContext.params.out.height=ptr[1];
    resizeContext.params.out.stride=ptr[2];
    ptr=vecGetElem_i32(&context->dataSet.V,0);
    resizeContext.isFast=ptr[0];
    switch(ptr[1])
    {
    case 0: resizeContext.params.method=imgresize_method_nearest; break;
    case 1: resizeContext.params.method=imgresize_method_bilinear; break;
    case 2: resizeContext.params.method=imgresize_method_bicubic; break;
    default: ASSERT("unsupported resize method");
    }
    resizeContext.format=ptr[2];
    te_freeVectors(context);

    szIn =resizeContext.params.in.height*resizeContext.params.in.stride;
    szOut=resizeContext.params.out.height*resizeContext.params.out.stride;
    scrSz=pAPI->getScratchSize(&resizeContext.params);
    inst_sz=pAPI->alloc(&resizeContext.params);

    /* Allocate output/reference data vectors memory. */
    nAlloc=vecsAlloc( context->desc->isAligned | ALLOC_DRAM0, context->desc->fmt,
                        &context->dataSet.X  , szIn, 
                        &context->dataSet.Z  , szOut,
                        &context->dataSet.Zlo, szOut,
                        &context->dataSet.Zhi, szOut, 0 ) ;
    nAlloc+=vecsAlloc( 1, FMT_UINT8,
                        &context->dataSet.U  , inst_sz+sizeof(tResizeContext), 
                        &context->dataSet.V  , scrSz, 
                        NULL) ;
    if ( nAlloc!=6 )
    {
        printf( "te_loadFxn_resize(): failed to allocate images\n");
        te_freeVectors(context);
    }
    /* Load vectors data from the SEQ-file. */
    else if ( !seqFileReadVecs( context->seqFile,
                                &context->dataSet.X,
                                &context->dataSet.Zlo,
                                &context->dataSet.Zhi, 0 ) )
    {
        printf( "te_loadFxn_resize(): failed to read vectors data\n");
        te_freeVectors(context);
    }
    else
    {
        tResizeContext* pResizeContext;
        pResizeContext=(tResizeContext*)vecGetElem(&context->dataSet.U,0);
        pResizeContext[0]=resizeContext;
        copyImgRandOutput(vecGetElem(&context->dataSet.Z,0),vecGetElem(&context->dataSet.Zlo,0),&resizeContext.params.out,context->dataSet.Zlo.szElem);
        pResizeContext->handle=pAPI->init((void*)(pResizeContext+1),&resizeContext.params);
        res = 1;
    }
    return (res);
} 

/* test resize functions */
static void te_processFxn_resize( tTestEngContext * context )
{
    tResizeApi *pAPI;
    void *outImg, *inImg, *pScr;
    tResizeContext* pResizeContext;
    pResizeContext=(tResizeContext*)vecGetElem(&context->dataSet.U,0);
    inImg  = vecGetElem( &context->dataSet.X, 0 );
    outImg = vecGetElem( &context->dataSet.Z, 0 );
    pScr   = vecGetElem( &context->dataSet.V, 0 );
    ASSERT( context && context->target.fut );
    te_vReportStd(context);
    pAPI=(tResizeApi*)context->desc->extraPtr;
    pAPI->process(pResizeContext->handle,pScr,outImg,inImg);
}

/* vec API test definitions. */
static const tResizeApi imgresize_gu8_api    =
{imgresize_gu8_alloc    ,
imgresize_gu8_init    ,
imgresize_gu8_process    ,
imgresize_gu8_getScratchSize    };
static const tResizeApi imgfastresize_gu8_api=
{imgfastresize_gu8_alloc,
imgfastresize_gu8_init,
imgfastresize_gu8_process,
imgfastresize_gu8_getScratchSize};

static const tResizeApi imgresize_gs8_api    =
{imgresize_gs8_alloc    ,
imgresize_gs8_init    ,
imgresize_gs8_process    ,
imgresize_gs8_getScratchSize    };
static const tResizeApi imgfastresize_gs8_api=
{imgfastresize_gs8_alloc,
imgfastresize_gs8_init,
imgfastresize_gs8_process,
imgfastresize_gs8_getScratchSize};

static const tResizeApi imgresize_gs16_api    ={
imgresize_gs16_alloc   ,
imgresize_gs16_init    ,
imgresize_gs16_process ,
imgresize_gs16_getScratchSize    };
static const tResizeApi imgfastresize_gs16_api={
imgfastresize_gs16_alloc  ,
imgfastresize_gs16_init   ,
imgfastresize_gs16_process,
imgfastresize_gs16_getScratchSize};

typedef struct  
{
  tTestEngTarget   funcList[MAX_FUNC_NUM];
  tTestEngDesc     testDesc;
}
deftbl_t;

/* Perform all tests for image processing APIs. */
int func_imgresize1(int isFull, int isVerbose, int breakOnError)
{
	int res = 1;
    static const deftbl_t testDefTblresize[] =
    {
     { FUNC_LIST( (tTestEngTarget)&imgresize_gu8_process ),      TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_resize, &te_processFxn_resize, &imgresize_gu8_api ) },
     { FUNC_LIST( (tTestEngTarget)&imgresize_gs8_process ),      TEST_DESC( FMT_REAL|FMT_INT8, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_resize, &te_processFxn_resize, &imgresize_gs8_api ) },
     { FUNC_LIST( (tTestEngTarget)&imgresize_gs16_process ),     TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_3, TE_ALIGN_NO, &te_loadFxn_resize, &te_processFxn_resize, &imgresize_gs16_api ) },
     { FUNC_LIST( (tTestEngTarget)&imgfastresize_gu8_process),   TEST_DESC( FMT_REAL|FMT_UINT8, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_resize, &te_processFxn_resize, &imgfastresize_gu8_api ) },  
     { FUNC_LIST( (tTestEngTarget)&imgfastresize_gs8_process),   TEST_DESC( FMT_REAL|FMT_INT8, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_resize, &te_processFxn_resize, &imgfastresize_gs8_api ) },  
     { FUNC_LIST( (tTestEngTarget)&imgfastresize_gs16_process),  TEST_DESC( FMT_REAL|FMT_INT16, TE_DIM_NUM_3, TE_ALIGN_YES, &te_loadFxn_resize, &te_processFxn_resize, &imgfastresize_gs16_api ) },  
     { FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, NULL, NULL, NULL ) } /* End of table */
    };
    
    // resize subset
    #define DO_TEST(fxn, seqFile)                                                                           \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTblresize,                                    \
                                                       sizeof(testDefTblresize)/sizeof(testDefTblresize[0]),\
                                                       MAX_FUNC_NUM,                                        \
                                                       (tTestEngTarget)(fxn), "imgresize1/" seqFile,               \
                                                       isFull, isVerbose, breakOnError ) )
	vecInitRegion(ALLOC_DRAM0);
    DO_TEST( &imgresize_gu8_process     , "imgresize_gu8_bilinear.seq"       );
    DO_TEST( &imgfastresize_gu8_process , "imgfastresize_gu8_bilinear.seq"   );
    DO_TEST( &imgresize_gs8_process     , "imgresize_gs8_bilinear.seq"       );
    DO_TEST( &imgfastresize_gs8_process , "imgfastresize_gs8_bilinear.seq"   );
    DO_TEST( &imgresize_gs16_process    , "imgresize_gs16_bilinear.seq"      );
    DO_TEST( &imgfastresize_gs16_process, "imgfastresize_gs16_bilinear.seq"  );
    DO_TEST( &imgresize_gu8_process     , "imgresize_gu8_bicubic.seq"        );
    DO_TEST( &imgfastresize_gu8_process , "imgfastresize_gu8_bicubic.seq"    );
    DO_TEST( &imgresize_gs8_process     , "imgresize_gs8_bicubic.seq"        );
    DO_TEST( &imgfastresize_gs8_process , "imgfastresize_gs8_bicubic.seq"    );
    DO_TEST( &imgresize_gs16_process    , "imgresize_gs16_bicubic.seq"       );
    DO_TEST( &imgfastresize_gs16_process, "imgfastresize_gs16_bicubic.seq"   );
    DO_TEST( &imgresize_gu8_process     , "imgresize_gu8_nearest.seq"        );
    DO_TEST( &imgfastresize_gu8_process , "imgfastresize_gu8_nearest.seq"    );
    DO_TEST( &imgresize_gs8_process     , "imgresize_gs8_nearest.seq"        );
    DO_TEST( &imgfastresize_gs8_process , "imgfastresize_gs8_nearest.seq"    );
    DO_TEST( &imgresize_gs16_process    , "imgresize_gs16_nearest.seq"       );
    DO_TEST( &imgfastresize_gs16_process, "imgfastresize_gs16_nearest.seq"   );
    #undef DO_TEST
    return (res);
}
