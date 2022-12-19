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
* Test module for testing cycle performance (image processing APIs.)
*/

#include <string.h>
#include <math.h>
/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API. */
#include LIBRARY_HEADER(img)
#include "NatureDSP_Signal_img.h"
/* MIPS measurement means. */
#include "mips.h"
/* Utility functions and macros. */
#include "utils.h"
#include <stdlib.h>
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

#define TEST_ODD_ANGLES 1

#define PROFILE_IMG_ROTATE(isFull,isVerbose,basename,fast,_file,resolution,w,h,angle,szpixel) {  \
        size_t szObj, szScr;                                                \
        imgsize_t szOut;                                                    \
        int bytesIn,bytesOut;                                               \
        void* memObj;                                                       \
        void* pScr;                                                         \
        void* pImg;                                                         \
        int stride = fast?((w+7)&~7):w+1;                                   \
        imgrotate_params_t params;                                          \
        imgrotate_handle_t handle;                                          \
        params.in.width=w;                                                  \
        params.in.height=h;                                                 \
        params.in.stride=stride;                                            \
        params.fill=0;                                                      \
        params.angleQ15=(int16_t)(((int64_t)11930464L*angle+32768)>>16);    \
        szObj=basename##_alloc  (&params);                                  \
        szScr=basename##_getScratchSize (&params );                         \
        basename##_getOutSize (&szOut,&params );                            \
        szScr=(szScr+15)&~15;                                               \
        bytesIn=szpixel*w*stride;                                           \
        bytesOut=szOut.width*szOut.stride;                                  \
        NASSERT(szScr+MAX(bytesIn,bytesOut)<sizeof(mips));                  \
        (void)szObj,(void)szScr,(void)bytesIn,(void)bytesOut;               \
        pScr=(void*)&mips;                                                  \
        pImg=(void*)(((uintptr_t)pScr)+szScr);                              \
        memObj=(void*)objinstance_memory;                                   \
        NASSERT(szObj<sizeof(objinstance_memory));                          \
        handle=basename##_init(memObj,&params);                             \
        PROFILE_SIMPLE(isFull,isVerbose,basename##_process,                 \
                       (handle,pScr,pImg,pImg,&szOut),                      \
                       _file,resolution " " #angle " degrees",prf_cycle);   \
}

#define PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull,isVerbose,basename,fast,_file,resolution,w,h,angle,szpixel) {  \
        size_t szObj, szScr;                                                \
        imgsize_t szOut;                                                    \
        int bytesIn,bytesOut;                                               \
        void* memObj;                                                       \
        void* pScr;                                                         \
        void* pImg;                                                         \
        void* pOutImg;                                                      \
        int stride = fast?((w+7)&~7):w+1;                                   \
        imgrotate_params_t params;                                          \
        imgrotate_handle_t handle;                                          \
        params.in.width=w;                                                  \
        params.in.height=h;                                                 \
        params.in.stride=stride;                                            \
        params.fill=0;                                                      \
        params.angleQ15=(int16_t)(((int64_t)11930464L*angle+32768)>>16);    \
        szObj=basename##_alloc  (&params);                                  \
        szScr=basename##_getScratchSize (&params );                         \
        basename##_getOutSize (&szOut,&params );                            \
        szScr=(szScr+15)&~15;                                               \
        bytesIn=szpixel*w*stride;                                           \
        bytesOut=szOut.width*szOut.stride;                                  \
        NASSERT(szScr + (bytesIn + bytesOut)<sizeof(mips));                 \
        (void)szObj,(void)szScr,(void)bytesIn,(void)bytesOut;               \
        pScr=(void*)&mips;                                                  \
        pImg=(void*)(((uintptr_t)pScr)+szScr);                              \
        pOutImg=(void*)(((uintptr_t)pImg)+bytesIn);                         \
        memObj=(void*)objinstance_memory;                                   \
        NASSERT(szObj<sizeof(objinstance_memory));                          \
        handle=basename##_init(memObj,&params);                             \
        PROFILE_SIMPLE(isFull,isVerbose,basename##_process,                 \
                       /*(handle,pScr,pImg,pImg,&szOut), */                 \
                        (handle,pScr,pOutImg,pImg,&szOut),                  \
                       _file,resolution " " #angle " degrees",prf_cycle);   \
}

static void mips_imgrotate_gu8(int isFull, int isVerbose, FILE * fout)
{
    (void)isFull;
    //8-bit unsigned
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gu8,0,fout,"SQCIF[128x96]",128,96,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gu8,0,fout,"SQCIF[128x96]",128,96, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gu8,0,fout,"SQCIF[128x96]",128,96,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gu8,0,fout,"SQCIF[128x96]",128,96,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gu8,0,fout,"SQCIF[128x96]",128,96, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"SQCIF[128x96]",128,96,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"SQCIF[128x96]",128,96,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"SQCIF[128x96]",128,96,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gu8,0,fout,"QCIF[176x144]",176,144,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gu8,0,fout,"QCIF[176x144]",176,144, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gu8,0,fout,"QCIF[176x144]",176,144,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gu8,0,fout,"QCIF[176x144]",176,144,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gu8,0,fout,"QCIF[176x144]",176,144, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"QCIF[176x144]",176,144,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"QCIF[176x144]",176,144,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"QCIF[176x144]",176,144,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gu8,0,fout,"CIF[352x288]",352,288,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gu8,0,fout,"CIF[352x288]",352,288, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gu8,0,fout,"CIF[352x288]",352,288,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gu8,0,fout,"CIF[352x288]",352,288,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gu8,0,fout,"CIF[352x288]",352,288, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"CIF[352x288]",352,288,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"CIF[352x288]",352,288,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"CIF[352x288]",352,288,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gu8,0,fout,"QVGA[320x240]",320,240,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gu8,0,fout,"QVGA[320x240]",320,240, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gu8,0,fout,"QVGA[320x240]",320,240,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gu8,0,fout,"QVGA[320x240]",320,240,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gu8,0,fout,"QVGA[320x240]",320,240, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"QVGA[320x240]",320,240,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"QVGA[320x240]",320,240,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"QVGA[320x240]",320,240,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gu8,0,fout,"VGA[640x480]",640,480,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gu8,0,fout,"VGA[640x480]",640,480, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gu8,0,fout,"VGA[640x480]",640,480,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gu8,0,fout,"VGA[640x480]",640,480,270,1);
#if TEST_ODD_ANGLES
//    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gu8,0,fout,"VGA[640x480]",640,480, 45,1);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"VGA[640x480]",640,480,135,1);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"VGA[640x480]",640,480,225,1);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gu8,0,fout,"VGA[640x480]",640,480,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gu8,1,fout,"SQCIF[128x96]",128,96,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gu8,1,fout,"SQCIF[128x96]",128,96, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gu8,1,fout,"SQCIF[128x96]",128,96,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gu8,1,fout,"SQCIF[128x96]",128,96,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gu8,1,fout,"SQCIF[128x96]",128,96, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"SQCIF[128x96]",128,96,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"SQCIF[128x96]",128,96,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"SQCIF[128x96]",128,96,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gu8,1,fout,"QCIF[176x144]",176,144,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gu8,1,fout,"QCIF[176x144]",176,144, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gu8,1,fout,"QCIF[176x144]",176,144,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gu8,1,fout,"QCIF[176x144]",176,144,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gu8,1,fout,"QCIF[176x144]",176,144, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"QCIF[176x144]",176,144,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"QCIF[176x144]",176,144,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"QCIF[176x144]",176,144,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gu8,1,fout,"CIF[352x288]",352,288,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gu8,1,fout,"CIF[352x288]",352,288, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gu8,1,fout,"CIF[352x288]",352,288,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gu8,1,fout,"CIF[352x288]",352,288,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gu8,1,fout,"CIF[352x288]",352,288, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"CIF[352x288]",352,288,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"CIF[352x288]",352,288,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"CIF[352x288]",352,288,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gu8,1,fout,"QVGA[320x240]",320,240,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gu8,1,fout,"QVGA[320x240]",320,240, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gu8,1,fout,"QVGA[320x240]",320,240,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gu8,1,fout,"QVGA[320x240]",320,240,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gu8,1,fout,"QVGA[320x240]",320,240, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"QVGA[320x240]",320,240,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"QVGA[320x240]",320,240,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"QVGA[320x240]",320,240,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gu8,1,fout,"VGA[640x480]",640,480,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gu8,1,fout,"VGA[640x480]",640,480, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gu8,1,fout,"VGA[640x480]",640,480,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gu8,1,fout,"VGA[640x480]",640,480,270,1);
#if TEST_ODD_ANGLES
//    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gu8,1,fout,"VGA[640x480]",640,480, 45,1);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"VGA[640x480]",640,480,135,1);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"VGA[640x480]",640,480,225,1);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gu8,1,fout,"VGA[640x480]",640,480,315,1);
#endif
} /* mips_imgrotate_gu8() */

static void mips_imgrotate_gs8(int isFull, int isVerbose, FILE * fout)
{
    (void)isFull;
 //8-bit signed
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs8,0,fout,"SQCIF[128x96]",128,96,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs8,0,fout,"SQCIF[128x96]",128,96, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs8,0,fout,"SQCIF[128x96]",128,96,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs8,0,fout,"SQCIF[128x96]",128,96,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gs8,0,fout,"SQCIF[128x96]",128,96, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"SQCIF[128x96]",128,96,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"SQCIF[128x96]",128,96,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"SQCIF[128x96]",128,96,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs8,0,fout,"QCIF[176x144]",176,144,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs8,0,fout,"QCIF[176x144]",176,144, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs8,0,fout,"QCIF[176x144]",176,144,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs8,0,fout,"QCIF[176x144]",176,144,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gs8,0,fout,"QCIF[176x144]",176,144, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"QCIF[176x144]",176,144,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"QCIF[176x144]",176,144,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"QCIF[176x144]",176,144,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs8,0,fout,"CIF[352x288]",352,288,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs8,0,fout,"CIF[352x288]",352,288, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs8,0,fout,"CIF[352x288]",352,288,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs8,0,fout,"CIF[352x288]",352,288,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gs8,0,fout,"CIF[352x288]",352,288, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"CIF[352x288]",352,288,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"CIF[352x288]",352,288,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"CIF[352x288]",352,288,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs8,0,fout,"QVGA[320x240]",320,240,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs8,0,fout,"QVGA[320x240]",320,240, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs8,0,fout,"QVGA[320x240]",320,240,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs8,0,fout,"QVGA[320x240]",320,240,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gs8,0,fout,"QVGA[320x240]",320,240, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"QVGA[320x240]",320,240,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"QVGA[320x240]",320,240,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"QVGA[320x240]",320,240,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs8,0,fout,"VGA[640x480]",640,480,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs8,0,fout,"VGA[640x480]",640,480, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs8,0,fout,"VGA[640x480]",640,480,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs8,0,fout,"VGA[640x480]",640,480,270,1);
#if TEST_ODD_ANGLES
//    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gs8,0,fout,"VGA[640x480]",640,480, 45,1);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"VGA[640x480]",640,480,135,1);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"VGA[640x480]",640,480,225,1);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs8,0,fout,"VGA[640x480]",640,480,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs8,1,fout,"SQCIF[128x96]",128,96,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs8,1,fout,"SQCIF[128x96]",128,96, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs8,1,fout,"SQCIF[128x96]",128,96,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs8,1,fout,"SQCIF[128x96]",128,96,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gs8,1,fout,"SQCIF[128x96]",128,96, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"SQCIF[128x96]",128,96,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"SQCIF[128x96]",128,96,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"SQCIF[128x96]",128,96,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs8,1,fout,"QCIF[176x144]",176,144,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs8,1,fout,"QCIF[176x144]",176,144, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs8,1,fout,"QCIF[176x144]",176,144,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs8,1,fout,"QCIF[176x144]",176,144,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gs8,1,fout,"QCIF[176x144]",176,144, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"QCIF[176x144]",176,144,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"QCIF[176x144]",176,144,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"QCIF[176x144]",176,144,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs8,1,fout,"CIF[352x288]",352,288,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs8,1,fout,"CIF[352x288]",352,288, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs8,1,fout,"CIF[352x288]",352,288,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs8,1,fout,"CIF[352x288]",352,288,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gs8,1,fout,"CIF[352x288]",352,288, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"CIF[352x288]",352,288,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"CIF[352x288]",352,288,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"CIF[352x288]",352,288,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs8,1,fout,"QVGA[320x240]",320,240,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs8,1,fout,"QVGA[320x240]",320,240, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs8,1,fout,"QVGA[320x240]",320,240,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs8,1,fout,"QVGA[320x240]",320,240,270,1);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gs8,1,fout,"QVGA[320x240]",320,240, 45,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"QVGA[320x240]",320,240,135,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"QVGA[320x240]",320,240,225,1);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"QVGA[320x240]",320,240,315,1);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs8,1,fout,"VGA[640x480]",640,480,  0,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs8,1,fout,"VGA[640x480]",640,480, 90,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs8,1,fout,"VGA[640x480]",640,480,180,1);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs8,1,fout,"VGA[640x480]",640,480,270,1);
#if TEST_ODD_ANGLES
//    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gs8,1,fout,"VGA[640x480]",640,480, 45,1);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"VGA[640x480]",640,480,135,1);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"VGA[640x480]",640,480,225,1);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs8,1,fout,"VGA[640x480]",640,480,315,1);
#endif
} /* mips_imgrotate_gs8() */

static void mips_imgrotate_gs16(int isFull, int isVerbose, FILE * fout)
{
    (void)isFull;
    // 16-bit
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs16,0,fout,"SQCIF[128x96]",128,96,  0,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs16,0,fout,"SQCIF[128x96]",128,96, 90,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs16,0,fout,"SQCIF[128x96]",128,96,180,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs16,0,fout,"SQCIF[128x96]",128,96,270,2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gs16,0,fout,"SQCIF[128x96]",128,96, 45,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"SQCIF[128x96]",128,96,135,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"SQCIF[128x96]",128,96,225,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"SQCIF[128x96]",128,96,315,2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs16,0,fout,"QCIF[176x144]",176,144,  0,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs16,0,fout,"QCIF[176x144]",176,144, 90,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs16,0,fout,"QCIF[176x144]",176,144,180,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs16,0,fout,"QCIF[176x144]",176,144,270,2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gs16,0,fout,"QCIF[176x144]",176,144, 45,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"QCIF[176x144]",176,144,135,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"QCIF[176x144]",176,144,225,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"QCIF[176x144]",176,144,315,2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs16,0,fout,"CIF[352x288]",352,288,  0,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs16,0,fout,"CIF[352x288]",352,288, 90,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs16,0,fout,"CIF[352x288]",352,288,180,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs16,0,fout,"CIF[352x288]",352,288,270,2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gs16,0,fout,"CIF[352x288]",352,288, 45,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"CIF[352x288]",352,288,135,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"CIF[352x288]",352,288,225,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"CIF[352x288]",352,288,315,2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs16,0,fout,"QVGA[320x240]",320,240,  0,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs16,0,fout,"QVGA[320x240]",320,240, 90,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs16,0,fout,"QVGA[320x240]",320,240,180,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs16,0,fout,"QVGA[320x240]",320,240,270,2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gs16,0,fout,"QVGA[320x240]",320,240, 45,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"QVGA[320x240]",320,240,135,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"QVGA[320x240]",320,240,225,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"QVGA[320x240]",320,240,315,2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs16,0,fout,"VGA[640x480]",640,480,  0,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgrotate_gs16,0,fout,"VGA[640x480]",640,480, 90,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs16,0,fout,"VGA[640x480]",640,480,180,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgrotate_gs16,0,fout,"VGA[640x480]",640,480,270,2);
#if TEST_ODD_ANGLES
//    PROFILE_IMG_ROTATE(1,      isVerbose,imgrotate_gs16,0,fout,"VGA[640x480]",640,480, 45,2);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"VGA[640x480]",640,480,135,2);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"VGA[640x480]",640,480,225,2);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgrotate_gs16,0,fout,"VGA[640x480]",640,480,315,2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs16,1,fout,"SQCIF[128x96]",128,96,  0,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs16,1,fout,"SQCIF[128x96]",128,96, 90,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs16,1,fout,"SQCIF[128x96]",128,96,180,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs16,1,fout,"SQCIF[128x96]",128,96,270,2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gs16,1,fout,"SQCIF[128x96]",128,96, 45,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"SQCIF[128x96]",128,96,135,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"SQCIF[128x96]",128,96,225,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"SQCIF[128x96]",128,96,315,2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs16,1,fout,"QCIF[176x144]",176,144,  0,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs16,1,fout,"QCIF[176x144]",176,144, 90,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs16,1,fout,"QCIF[176x144]",176,144,180,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs16,1,fout,"QCIF[176x144]",176,144,270,2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gs16,1,fout,"QCIF[176x144]",176,144, 45,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"QCIF[176x144]",176,144,135,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"QCIF[176x144]",176,144,225,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"QCIF[176x144]",176,144,315,2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs16,1,fout,"CIF[352x288]",352,288,  0,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs16,1,fout,"CIF[352x288]",352,288, 90,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs16,1,fout,"CIF[352x288]",352,288,180,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs16,1,fout,"CIF[352x288]",352,288,270,2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gs16,1,fout,"CIF[352x288]",352,288, 45,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"CIF[352x288]",352,288,135,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"CIF[352x288]",352,288,225,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"CIF[352x288]",352,288,315,2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs16,1,fout,"QVGA[320x240]",320,240,  0,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs16,1,fout,"QVGA[320x240]",320,240, 90,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs16,1,fout,"QVGA[320x240]",320,240,180,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs16,1,fout,"QVGA[320x240]",320,240,270,2);
#if TEST_ODD_ANGLES
    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gs16,1,fout,"QVGA[320x240]",320,240, 45,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"QVGA[320x240]",320,240,135,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"QVGA[320x240]",320,240,225,2);
    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"QVGA[320x240]",320,240,315,2);
#endif
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs16,1,fout,"VGA[640x480]",640,480,  0,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(1,      isVerbose,imgfastrotate_gs16,1,fout,"VGA[640x480]",640,480, 90,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs16,1,fout,"VGA[640x480]",640,480,180,2);
    PROFILE_IMG_ROTATE_EVEN_ANGLES(isFull, isVerbose,imgfastrotate_gs16,1,fout,"VGA[640x480]",640,480,270,2);
#if TEST_ODD_ANGLES
//    PROFILE_IMG_ROTATE(1,      isVerbose,imgfastrotate_gs16,1,fout,"VGA[640x480]",640,480, 45,2);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"VGA[640x480]",640,480,135,2);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"VGA[640x480]",640,480,225,2);
//    PROFILE_IMG_ROTATE(isFull, isVerbose,imgfastrotate_gs16,1,fout,"VGA[640x480]",640,480,315,2);
#endif
} /* mips_imgrotate_gs16() */

void mips_imgrotate1(int isFull, int isVerbose, FILE * fout)
{
  mips_imgrotate_gu8(isFull, isVerbose, fout);
  mips_imgrotate_gs8(isFull, isVerbose, fout);
  mips_imgrotate_gs16(isFull, isVerbose, fout);
}
