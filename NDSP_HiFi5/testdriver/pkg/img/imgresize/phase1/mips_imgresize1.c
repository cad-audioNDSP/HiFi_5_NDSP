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

#define PROFILE_IMG_RESIZE(isFull,isVerbose,basename,fast,method,_file,resolution,win,hin,wout,hout,method_str,szpixel) {  \
        size_t szObj, szScr;                                         \
        void* memObj;                                                \
        void* memScr;                                                \
        int istride = fast?(win+7)&~7:win+1;                         \
        int ostride = fast?(wout+7)&~7:wout+1;                       \
        imgresize_params_t params={{win,hin,istride},                \
                                   {wout,hout,ostride},method};      \
        imgresize_handle_t handle;                                   \
        szObj=basename##_alloc  (&params);                           \
        szScr=basename##_getScratchSize (&params );                  \
        NASSERT(hin*istride*szpixel<sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2)); \
        NASSERT(hout*ostride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2));\
        /*NASSERT(szScr<sizeof(tProfiler_scratch)+sizeof(mips.inp0)); */ \
        (void)szScr;                                                   \
        memObj=malloc(szObj);                                          \
        memScr=mallocAlign(szScr, 4 * NatureDSP_Signal_get_isa_opt(NATUREDSP_ISA_OPT_INT16_SIMD_WIDTH));  \
        handle=basename##_init(memObj,&params);                        \
        PROFILE_SIMPLE(isFull,isVerbose,basename##_process,            \
                       (handle,memScr,(void*)&mips.out0,               \
                       (const void*)&mips.inp0),                       \
                       _file,resolution ", " method_str,prf_cycle);    \
        freeAlign(memScr);                                             \
        free(memObj);                                                  \
}

static void mips_resize_nearest_gu8(int isFull, int isVerbose, FILE * fout)
{
    // 8-bit unsigned
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "nearest",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"nearest",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"nearest",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"nearest",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_nearest,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"nearest",1);

    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "nearest",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"nearest",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"nearest",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"nearest",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"nearest",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_nearest,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"nearest",1);
}

static void mips_resize_nearest_gs8(int isFull, int isVerbose, FILE * fout)
{
    // 8-bit signed
    PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "SQCIF[128x96]->QCIF[176x144]", 128, 96, 176, 144, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "SQCIF[128x96]->CIF [352x288]", 128, 96, 352, 288, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "SQCIF[128x96]->QVGA[320x240]", 128, 96, 320, 240, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "SQCIF[128x96]->VGA [640x480]", 128, 96, 640, 480, "nearest", 1);
    PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QCIF[176x144]->SQCIF[128x96]", 176, 144, 128, 96, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QCIF[176x144]->CIF [352x288]", 176, 144, 352, 288, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QCIF[176x144]->QVGA[320x240]", 176, 144, 320, 240, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QCIF[176x144]->VGA [640x480]", 176, 144, 640, 480, "nearest", 1);
    PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "CIF [352x288]->SQCIF[128x96]", 352, 288, 128, 96, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "CIF [352x288]->QCIF[176x144]", 352, 288, 176, 144, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "CIF [352x288]->QVGA[320x240]", 352, 288, 320, 240, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "CIF [352x288]->VGA [640x480]", 352, 288, 640, 480, "nearest", 1);
    PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QVGA[320x240]->SQCIF[128x96]", 320, 240, 128, 96, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QVGA[320x240]->QCIF[176x144]", 320, 240, 176, 144, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QVGA[320x240]->CIF [352x288]", 320, 240, 352, 288, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "QVGA[320x240]->VGA [640x480]", 320, 240, 640, 480, "nearest", 1);
    PROFILE_IMG_RESIZE(1, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "VGA [640x480]->SQCIF[128x96]", 640, 480, 128, 96, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "VGA [640x480]->QCIF[176x144]", 640, 480, 176, 144, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "VGA [640x480]->CIF [352x288]", 640, 480, 352, 288, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgresize_gs8, 0, imgresize_method_nearest, fout, "VGA [640x480]->QVGA[320x240]", 640, 480, 320, 240, "nearest", 1);

    PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "SQCIF[128x96]->QCIF[176x144]", 128, 96, 176, 144, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "SQCIF[128x96]->CIF [352x288]", 128, 96, 352, 288, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "SQCIF[128x96]->QVGA[320x240]", 128, 96, 320, 240, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "SQCIF[128x96]->VGA [640x480]", 128, 96, 640, 480, "nearest", 1);
    PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QCIF[176x144]->SQCIF[128x96]", 176, 144, 128, 96, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QCIF[176x144]->CIF [352x288]", 176, 144, 352, 288, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QCIF[176x144]->QVGA[320x240]", 176, 144, 320, 240, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QCIF[176x144]->VGA [640x480]", 176, 144, 640, 480, "nearest", 1);
    PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "CIF [352x288]->SQCIF[128x96]", 352, 288, 128, 96, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "CIF [352x288]->QCIF[176x144]", 352, 288, 176, 144, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "CIF [352x288]->QVGA[320x240]", 352, 288, 320, 240, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "CIF [352x288]->VGA [640x480]", 352, 288, 640, 480, "nearest", 1);
    PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QVGA[320x240]->SQCIF[128x96]", 320, 240, 128, 96, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QVGA[320x240]->QCIF[176x144]", 320, 240, 176, 144, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QVGA[320x240]->CIF [352x288]", 320, 240, 352, 288, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "QVGA[320x240]->VGA [640x480]", 320, 240, 640, 480, "nearest", 1);
    PROFILE_IMG_RESIZE(1, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "VGA [640x480]->SQCIF[128x96]", 640, 480, 128, 96, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "VGA [640x480]->QCIF[176x144]", 640, 480, 176, 144, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "VGA [640x480]->CIF [352x288]", 640, 480, 352, 288, "nearest", 1);
    PROFILE_IMG_RESIZE(isFull, isVerbose, imgfastresize_gs8, 1, imgresize_method_nearest, fout, "VGA [640x480]->QVGA[320x240]", 640, 480, 320, 240, "nearest", 1);
}

static void mips_resize_nearest_gs16(int isFull, int isVerbose, FILE * fout)
{
    // 16-bit 
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_nearest,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"nearest",2);

    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"nearest",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"nearest",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_nearest,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"nearest",2);
}

static void mips_resize_nearest(int isFull, int isVerbose, FILE * fout)
{
  mips_resize_nearest_gu8(isFull, isVerbose, fout);
  mips_resize_nearest_gs8(isFull, isVerbose, fout);
  mips_resize_nearest_gs16(isFull, isVerbose, fout);
}

static void mips_resize_bilinear_gu8(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bilinear,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"bilinear",1);

    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bilinear,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"bilinear",1);
}

static void mips_resize_bilinear_gs8(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bilinear,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"bilinear",1);

    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"bilinear",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"bilinear",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bilinear,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"bilinear",1);
}

static void mips_resize_bilinear_gs16(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "bilinear",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"bilinear",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"bilinear",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"bilinear",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bilinear,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"bilinear",2);

    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "bilinear",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"bilinear",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"bilinear",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"bilinear",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"bilinear",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bilinear,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"bilinear",2);
}

static void mips_resize_bilinear(int isFull, int isVerbose, FILE * fout)
{
  mips_resize_bilinear_gu8(isFull, isVerbose, fout);
  mips_resize_bilinear_gs8(isFull, isVerbose, fout);
  mips_resize_bilinear_gs16(isFull, isVerbose, fout);
}

static void mips_resize_bicubic_gu8(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gu8,0,imgresize_method_bicubic,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"bicubic",1);

    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gu8,1,imgresize_method_bicubic,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"bicubic",1);
}

static void mips_resize_bicubic_gs8(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_IMG_RESIZE(1,    isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs8,0,imgresize_method_bicubic,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"bicubic",1);

    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"bicubic",1);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"bicubic",1);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs8,1,imgresize_method_bicubic,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"bicubic",1);
}

static void mips_resize_bicubic_gs16(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgresize_gs16,0,imgresize_method_bicubic,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"bicubic",2);

    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"SQCIF[128x96]->QCIF[176x144]",128,96, 176,144, "bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"SQCIF[128x96]->CIF [352x288]",128,96, 352,288, "bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"SQCIF[128x96]->QVGA[320x240]",128,96, 320,240, "bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"SQCIF[128x96]->VGA [640x480]",128,96, 640,480, "bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QCIF[176x144]->SQCIF[128x96]",176,144, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QCIF[176x144]->CIF [352x288]",176,144, 352,288,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QCIF[176x144]->QVGA[320x240]",176,144, 320,240,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QCIF[176x144]->VGA [640x480]",176,144, 640,480,"bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"CIF [352x288]->SQCIF[128x96]",352,288, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"CIF [352x288]->QCIF[176x144]",352,288, 176,144,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"CIF [352x288]->QVGA[320x240]",352,288, 320,240,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"CIF [352x288]->VGA [640x480]",352,288, 640,480,"bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QVGA[320x240]->SQCIF[128x96]",320,240, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QVGA[320x240]->QCIF[176x144]",320,240, 176,144,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QVGA[320x240]->CIF [352x288]",320,240, 352,288,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"QVGA[320x240]->VGA [640x480]",320,240, 640,480,"bicubic",2);
    PROFILE_IMG_RESIZE(1,     isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"VGA [640x480]->SQCIF[128x96]",640,480, 128, 96,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"VGA [640x480]->QCIF[176x144]",640,480, 176,144,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"VGA [640x480]->CIF [352x288]",640,480, 352,288,"bicubic",2);
    PROFILE_IMG_RESIZE(isFull,isVerbose,imgfastresize_gs16,1,imgresize_method_bicubic,fout,"VGA [640x480]->QVGA[320x240]",640,480, 320,240,"bicubic",2);
}

static void mips_resize_bicubic(int isFull, int isVerbose, FILE * fout)
{
  mips_resize_bicubic_gu8(isFull, isVerbose, fout);
  mips_resize_bicubic_gs8(isFull, isVerbose, fout);
  mips_resize_bicubic_gs16(isFull, isVerbose, fout);
}

void mips_imgresize1(int isFull, int isVerbose, FILE * fout)
{
    mips_resize_nearest(isFull,isVerbose,fout);
    mips_resize_bilinear(isFull,isVerbose,fout);
    mips_resize_bicubic(isFull,isVerbose,fout);
}
