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

static const int32_t rgbyuv_transform[]={160524403,315143225,61203284,-78989817,-155080532,234075718,329638740,-276483151,-53692460,611941572,-211876105,-311755570,1090980749};

#define PROFILE_IMG_HIST(isFull,isVerbose,f_name,fast,_file,resolution,w,h) {    \
        int stride = fast?(w+7)&~7:w+1;                                          \
        imgsize_t sz={w,h,stride};                                               \
        imghist_t hist={0,0,0,mips.out1.i32};                                    \
        NASSERT(h*stride<sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2)); \
        PROFILE_SIMPLE(isFull,isVerbose,f_name,(&hist,                           \
                       (void*)mips.inp0.u8,&sz,256),_file,resolution,prf_cycle); \
        }

#define PROFILE_IMG_NORM(isFull,isVerbose,f_name,fast,_file,resolution,w,h,szpixel) {   \
        int stride = fast?(w+7)&~7:w+1;                                                 \
        imgsize_t sz={w,h,stride};                                                      \
        NASSERT(h*stride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2));\
        NASSERT(h*stride*szpixel<sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2));\
        PROFILE_SIMPLE(isFull,isVerbose,f_name,((void*)mips.out0.u8,                    \
                       (void*)mips.inp0.u8,&sz,0,255),_file,resolution,prf_cycle);      \
        }
#define PROFILE_IMG_NORM_NONLINEAR(isFull,isVerbose,f_name,fast,_file,resolution,w,h,szpixel) {\
        int stride = fast?(w+7)&~7:w+1;                                                        \
        imgsize_t sz={w,h,stride};                                                             \
        NASSERT(h*stride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2));       \
        NASSERT(h*stride*szpixel<sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2));       \
        PROFILE_SIMPLE(isFull,isVerbose,f_name,((void*)mips.out0.u8,                           \
                       (void*)mips.inp0.u8,&sz,mips.inp1.i16),_file,resolution,prf_cycle);     \
        }

#define PROFILE_IMG_INTERLEAVE(isFull,isVerbose,f_name,fast,_file,resolution,w,h,szpixel) { \
        int stride = fast?(w+7)&~7:w+1;                                                     \
        imgsize_t sz={w,h,stride};                                                          \
        NASSERT(3*h*stride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2)+   \
                                   sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2));  \
        NASSERT(  h*stride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2)+   \
                                   sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2));  \
        PROFILE_SIMPLE(isFull,isVerbose,f_name,((void*)mips.out0.u8,                        \
                                                (const void*)mips.inp0.u8,                  \
                                                (const void*)mips.inp0.u8,                  \
                                                (const void*)mips.inp0.u8,&sz),             \
                                               _file,resolution,prf_cycle);                 \
        }
#define PROFILE_IMG_DEINTERLEAVE(isFull,isVerbose,f_name,fast,_file,resolution,w,h,szpixel) {\
        int stride = fast?(w+7)&~7:w+1;                                                      \
        imgsize_t sz={w,h,stride};                                                           \
        NASSERT(3*h*stride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2)+    \
                                   sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2));   \
        NASSERT(  h*stride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2)+    \
                                   sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2));   \
        PROFILE_SIMPLE(isFull,isVerbose,f_name,((void*)mips.out0.u8,                         \
                                                (void*)mips.out0.u8,                         \
                                                (void*)mips.out0.u8,                         \
                                                (const void*)mips.inp0.u8, &sz),             \
                                               _file,resolution,prf_cycle);                  \
        }

#define PROFILE_IMG_CONVERT(isFull,isVerbose,f_name,fast,_file,resolution,w,h,szpixel) { \
        int stride = fast?(w+7)&~7:w+1;                                                  \
        imgsize_t sz={w,h,stride};                                                       \
        NASSERT(h*stride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2)); \
        NASSERT(h*stride*szpixel<sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2)); \
        PROFILE_SIMPLE(isFull,isVerbose,f_name,((void*)mips.out0.u8,                     \
                                                (void*)mips.out0.u8,                     \
                                                (void*)mips.out0.u8,                     \
                                                (const void*)mips.inp0.u8,               \
                                                (const void*)mips.inp1.u8,               \
                                                (const void*)mips.inp2.u8,               \
                                                rgbyuv_transform,&sz),                   \
                                               _file,resolution,prf_cycle);              \
        }

#define PROFILE_IMG_PAD(isFull,isVerbose,f_name,fast,_file,resolution,win,hin,wout,hout,szpixel) \
{                                                                           \
        int istride = fast?(win+7)&~7:win+1;                                \
        int ostride = fast?(wout+7)&~7:wout+1;                              \
        imgsize_t szin ={win,hin,istride};                                  \
        imgsize_t szout={wout,hout,ostride};                                \
        int x=(win-wout)>>1;                                                \
        int y=(hin-hout)>>1;                                                \
        imgpad_params_t params={szin,szout,x,y,IMGPAD_EDGE};                \
        size_t sz=f_name##_getScratchSize(&params);                         \
        NASSERT(sz<=sizeof(mips.scratch0));                                 \
        NASSERT(hout*ostride*szpixel<sizeof(mips.out0)+sizeof(mips.out1)+sizeof(mips.out2));\
        NASSERT(hin*istride *szpixel<sizeof(mips.inp0)+sizeof(mips.inp1)+sizeof(mips.inp2));\
        (void)sz;                                                           \
        PROFILE_SIMPLE(isFull,isVerbose,f_name,((void*)mips.scratch0.u8,    \
                                                (void*)mips.out0.u8,        \
                                                (const void*)mips.inp0.u8,  \
                                                &params),                   \
                                               _file,resolution,prf_cycle); \
}

static void mips_hist(int isFull, int isVerbose, FILE * fout)
{
    (void)isFull;
    PROFILE_IMG_HIST(1, isVerbose,imghist_gu8,0,fout,"SQCIF[128x96]",128,96);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gu8,0,fout,"QCIF[176x144]",176,144);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gu8,0,fout,"CIF [352x288]",352,288);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gu8,0,fout,"QVGA[320x240]",320,240);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gu8,0,fout,"VGA [640x480]",640,480);

    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gu8,1,fout,"SQCIF[128x96]",128,96);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gu8,1,fout,"QCIF[176x144]",176,144);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gu8,1,fout,"CIF [352x288]",352,288);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gu8,1,fout,"QVGA[320x240]",320,240);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gu8,1,fout,"VGA [640x480]",640,480);

    PROFILE_IMG_HIST(1, isVerbose, imghist_gs8, 0, fout, "SQCIF[128x96]", 128, 96);
    PROFILE_IMG_HIST(1, isVerbose, imghist_gs8, 0, fout, "QCIF[176x144]", 176, 144);
    PROFILE_IMG_HIST(1, isVerbose, imghist_gs8, 0, fout, "CIF [352x288]", 352, 288);
    PROFILE_IMG_HIST(1, isVerbose, imghist_gs8, 0, fout, "QVGA[320x240]", 320, 240);
    PROFILE_IMG_HIST(1, isVerbose, imghist_gs8, 0, fout, "VGA [640x480]", 640, 480);

    PROFILE_IMG_HIST(1, isVerbose, imgfasthist_gs8, 1, fout, "SQCIF[128x96]", 128, 96);
    PROFILE_IMG_HIST(1, isVerbose, imgfasthist_gs8, 1, fout, "QCIF[176x144]", 176, 144);
    PROFILE_IMG_HIST(1, isVerbose, imgfasthist_gs8, 1, fout, "CIF [352x288]", 352, 288);
    PROFILE_IMG_HIST(1, isVerbose, imgfasthist_gs8, 1, fout, "QVGA[320x240]", 320, 240);
    PROFILE_IMG_HIST(1, isVerbose, imgfasthist_gs8, 1, fout, "VGA [640x480]", 640, 480);

    PROFILE_IMG_HIST(1, isVerbose,imghist_gs16,0,fout,"SQCIF[128x96]",128,96);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gs16,0,fout,"QCIF[176x144]",176,144);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gs16,0,fout,"CIF [352x288]",352,288);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gs16,0,fout,"QVGA[320x240]",320,240);
    PROFILE_IMG_HIST(1, isVerbose,imghist_gs16,0,fout,"VGA [640x480]",640,480);

    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gs16,1,fout,"SQCIF[128x96]",128,96);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gs16,1,fout,"QCIF[176x144]",176,144);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gs16,1,fout,"CIF [352x288]",352,288);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gs16,1,fout,"QVGA[320x240]",320,240);
    PROFILE_IMG_HIST(1, isVerbose,imgfasthist_gs16,1,fout,"VGA [640x480]",640,480);
}

static void mips_interleave(int isFull, int isVerbose, FILE * fout)
{
    (void)isFull;
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave  ,0,fout,"SQCIF[128x96]",128,96,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave  ,0,fout,"QCIF[176x144]",176,144,1);
    PROFILE_IMG_INTERLEAVE(1     , isVerbose,imginterleave  ,0,fout,"CIF [352x288]",352,288,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave  ,0,fout,"QVGA[320x240]",320,240,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave  ,0,fout,"VGA [640x480]",640,480,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave16,0,fout,"SQCIF[128x96]",128,96,2);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave16,0,fout,"QCIF[176x144]",176,144,2);
    PROFILE_IMG_INTERLEAVE(1     , isVerbose,imginterleave16,0,fout,"CIF [352x288]",352,288,2);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imginterleave16,0,fout,"QVGA[320x240]",320,240,2);

    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave  ,1,fout,"SQCIF[128x96]",128,96,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave  ,1,fout,"QCIF[176x144]",176,144,1);
    PROFILE_IMG_INTERLEAVE(1     , isVerbose,imgfastinterleave  ,1,fout,"CIF [352x288]",352,288,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave  ,1,fout,"QVGA[320x240]",320,240,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave  ,1,fout,"VGA [640x480]",640,480,1);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave16,1,fout,"SQCIF[128x96]",128,96,2);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave16,1,fout,"QCIF[176x144]",176,144,2);
    PROFILE_IMG_INTERLEAVE(1     , isVerbose,imgfastinterleave16,1,fout,"CIF [352x288]",352,288,2);
    PROFILE_IMG_INTERLEAVE(isFull, isVerbose,imgfastinterleave16,1,fout,"QVGA[320x240]",320,240,2);

    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave  ,0,fout,"SQCIF[128x96]",128,96,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave  ,0,fout,"QCIF[176x144]",176,144,1);
    PROFILE_IMG_DEINTERLEAVE(1     , isVerbose,imgdeinterleave  ,0,fout,"CIF [352x288]",352,288,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave  ,0,fout,"QVGA[320x240]",320,240,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave  ,0,fout,"VGA [640x480]",640,480,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave16,0,fout,"SQCIF[128x96]",128,96,2);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave16,0,fout,"QCIF[176x144]",176,144,2);
    PROFILE_IMG_DEINTERLEAVE(1     , isVerbose,imgdeinterleave16,0,fout,"CIF [352x288]",352,288,2);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgdeinterleave16,0,fout,"QVGA[320x240]",320,240,2);

    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave  ,1,fout,"SQCIF[128x96]",128,96,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave  ,1,fout,"QCIF[176x144]",176,144,1);
    PROFILE_IMG_DEINTERLEAVE(1     , isVerbose,imgfastdeinterleave  ,1,fout,"CIF(352x288)",352,288,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave  ,1,fout,"QVGA[320x240]",320,240,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave  ,1,fout,"VGA [640x480]",640,480,1);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave16,1,fout,"SQCIF[128x96]",128,96,2);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave16,1,fout,"QCIF[176x144]",176,144,2);
    PROFILE_IMG_DEINTERLEAVE(1     , isVerbose,imgfastdeinterleave16,1,fout,"CIF [352x288]",352,288,2);
    PROFILE_IMG_DEINTERLEAVE(isFull, isVerbose,imgfastdeinterleave16,1,fout,"QVGA[320x240]",320,240,2);
}

static void mips_convert(int isFull, int isVerbose, FILE * fout)
{
    (void)isFull;
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv  ,0,fout,"SQCIF[128x96]",128,96,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv  ,0,fout,"QCIF[176x144]",176,144,1);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgconvert_rgbyuv  ,0,fout,"CIF [352x288]",352,288,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv  ,0,fout,"QVGA[320x240]",320,240,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv  ,0,fout,"VGA [640x480]",640,480,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv16,0,fout,"SQCIF[128x96]",128,96,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv16,0,fout,"QCIF[176x144]",176,144,2);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgconvert_rgbyuv16,0,fout,"CIF [352x288]",352,288,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv16,0,fout,"QVGA[320x240]",320,240,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_rgbyuv16,0,fout,"VGA [640x480]",640,480,2);

    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv  ,1,fout,"SQCIF[128x96]",128,96,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv  ,1,fout,"QCIF[176x144]",176,144,1);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgfastconvert_rgbyuv  ,1,fout,"CIF [352x288]",352,288,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv  ,1,fout,"QVGA[320x240]",320,240,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv  ,1,fout,"VGA [640x480]",640,480,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv16,1,fout,"SQCIF[128x96]",128,96,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv16,1,fout,"QCIF[176x144]",176,144,2);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgfastconvert_rgbyuv16,1,fout,"CIF [352x288]",352,288,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv16,1,fout,"QVGA[320x240]",320,240,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_rgbyuv16,1,fout,"VGA [640x480]",640,480,2);

    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb  ,0,fout,"SQCIF[128x96]",128,96,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb  ,0,fout,"QCIF[176x144]",176,144,1);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgconvert_yuvrgb  ,0,fout,"CIF [352x288]",352,288,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb  ,0,fout,"QVGA[320x240]",320,240,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb  ,0,fout,"VGA [640x480]",640,480,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb16,0,fout,"SQCIF[128x96]",128,96,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb16,0,fout,"QCIF[176x144]",176,144,2);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgconvert_yuvrgb16,0,fout,"CIF [352x288]",352,288,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb16,0,fout,"QVGA[320x240]",320,240,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgconvert_yuvrgb16,0,fout,"VGA [640x480]",640,480,2);

    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb  ,1,fout,"SQCIF[128x96]",128,96,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb  ,1,fout,"QCIF[176x144]",176,144,1);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgfastconvert_yuvrgb  ,1,fout,"CIF [352x288]",352,288,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb  ,1,fout,"QVGA[320x240]",320,240,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb  ,1,fout,"VGA [640x480]",640,480,1);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb16,1,fout,"SQCIF[128x96]",128,96,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb16,1,fout,"QCIF[176x144]",176,144,2);
    PROFILE_IMG_CONVERT(1     , isVerbose,imgfastconvert_yuvrgb16,1,fout,"CIF [352x288]",352,288,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb16,1,fout,"QVGA[320x240]",320,240,2);
    PROFILE_IMG_CONVERT(isFull, isVerbose,imgfastconvert_yuvrgb16,1,fout,"VGA [640x480]",640,480,2);
}

static void mips_norm(int isFull, int isVerbose, FILE * fout)
{
    (void)isFull;
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gu8,0,fout,"SQCIF[128x96]",128,96,1);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gu8,0,fout,"QCIF[176x144]",176,144,1);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gu8,0,fout,"CIF [352x288]",352,288,1);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gu8,0,fout,"QVGA[320x240]",320,240,1);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gu8,0,fout,"VGA [640x480]",640,480,1);

    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gu8,1,fout,"SQCIF[128x96]",128,96,1);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gu8,1,fout,"QCIF[176x144]",176,144,1);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gu8,1,fout,"CIF [352x288]",352,288,1);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gu8,1,fout,"QVGA[320x240]",320,240,1);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gu8,1,fout,"VGA [640x480]",640,480,1);

    PROFILE_IMG_NORM(1, isVerbose, imgnorm_gs8, 0, fout, "SQCIF[128x96]", 128, 96, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgnorm_gs8, 0, fout, "QCIF[176x144]", 176, 144, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgnorm_gs8, 0, fout, "CIF [352x288]", 352, 288, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgnorm_gs8, 0, fout, "QVGA[320x240]", 320, 240, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgnorm_gs8, 0, fout, "VGA [640x480]", 640, 480, 1);

    PROFILE_IMG_NORM(1, isVerbose, imgfastnorm_gs8, 1, fout, "SQCIF[128x96]", 128, 96, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgfastnorm_gs8, 1, fout, "QCIF[176x144]", 176, 144, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgfastnorm_gs8, 1, fout, "CIF [352x288]", 352, 288, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgfastnorm_gs8, 1, fout, "QVGA[320x240]", 320, 240, 1);
    PROFILE_IMG_NORM(1, isVerbose, imgfastnorm_gs8, 1, fout, "VGA [640x480]", 640, 480, 1);

    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gs16,0,fout,"SQCIF[128x96]",128,96,2);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gs16,0,fout,"QCIF[176x144]",176,144,2);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gs16,0,fout,"CIF [352x288]",352,288,2);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gs16,0,fout,"QVGA[320x240]",320,240,2);
    PROFILE_IMG_NORM(1, isVerbose,imgnorm_gs16,0,fout,"VGA [640x480]",640,480,2);

    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gs16,1,fout,"SQCIF[128x96]",128,96,2);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gs16,1,fout,"QCIF[176x144]",176,144,2);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gs16,1,fout,"CIF [352x288]",352,288,2);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gs16,1,fout,"QVGA[320x240]",320,240,2);
    PROFILE_IMG_NORM(1, isVerbose,imgfastnorm_gs16,1,fout,"VGA [640x480]",640,480,2);
    // nonlinear
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gu8_nonlinear,0,fout,"SQCIF[128x96]",128,96,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gu8_nonlinear,0,fout,"QCIF[176x144]",176,144,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gu8_nonlinear,0,fout,"CIF [352x288]",352,288,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gu8_nonlinear,0,fout,"QVGA[320x240]",320,240,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gu8_nonlinear,0,fout,"VGA [640x480]",640,480,1);

    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gu8_nonlinear,1,fout,"SQCIF[128x96]",128,96,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gu8_nonlinear,1,fout,"QCIF[176x144]",176,144,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gu8_nonlinear,1,fout,"CIF [352x288]",352,288,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gu8_nonlinear,1,fout,"QVGA[320x240]",320,240,1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gu8_nonlinear,1,fout,"VGA [640x480]",640,480,1);

    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgnorm_gs8_nonlinear, 0, fout, "SQCIF[128x96]", 128, 96, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgnorm_gs8_nonlinear, 0, fout, "QCIF[176x144]", 176, 144, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgnorm_gs8_nonlinear, 0, fout, "CIF [352x288]", 352, 288, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgnorm_gs8_nonlinear, 0, fout, "QVGA[320x240]", 320, 240, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgnorm_gs8_nonlinear, 0, fout, "VGA [640x480]", 640, 480, 1);

    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgfastnorm_gs8_nonlinear, 1, fout, "SQCIF[128x96]", 128, 96, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgfastnorm_gs8_nonlinear, 1, fout, "QCIF[176x144]", 176, 144, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgfastnorm_gs8_nonlinear, 1, fout, "CIF [352x288]", 352, 288, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgfastnorm_gs8_nonlinear, 1, fout, "QVGA[320x240]", 320, 240, 1);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose, imgfastnorm_gs8_nonlinear, 1, fout, "VGA [640x480]", 640, 480, 1);

    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gs16_nonlinear,0,fout,"SQCIF[128x96]",128,96,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gs16_nonlinear,0,fout,"QCIF[176x144]",176,144,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gs16_nonlinear,0,fout,"CIF [352x288]",352,288,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gs16_nonlinear,0,fout,"QVGA[320x240]",320,240,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgnorm_gs16_nonlinear,0,fout,"VGA [640x480]",640,480,2);

    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gs16_nonlinear,1,fout,"SQCIF[128x96]",128,96,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gs16_nonlinear,1,fout,"QCIF[176x144]",176,144,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gs16_nonlinear,1,fout,"CIF [352x288]",352,288,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gs16_nonlinear,1,fout,"QVGA[320x240]",320,240,2);
    PROFILE_IMG_NORM_NONLINEAR(1, isVerbose,imgfastnorm_gs16_nonlinear,1,fout,"VGA [640x480]",640,480,2);
}

static void mips_pad(int isFull, int isVerbose, FILE * fout)
{
    // pad
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gu8,0,fout,"padding, SQCIF[128x96]->QCIF[176x144]",128,96, 176,144  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, SQCIF[128x96]->CIF [352x288]",128,96, 352,288  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, SQCIF[128x96]->QVGA[320x240]",128,96, 320,240  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, SQCIF[128x96]->VGA [640x480]" ,128,96, 640,480 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, QCIF[176x144]->CIF [352x288]",176,144, 352,288 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, QCIF[176x144]->QVGA[320x240]",176,144, 320,240 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, QCIF[176x144]->VGA [640x480]",176,144, 640,480 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, QVGA[320x240]->CIF [352x288]",320,240, 352,288 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, QVGA[320x240]->VGA [640x480]",320,240, 640,480 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"padding, CIF [352x288]->VGA [640x480]",352,288, 640,480 ,1);

    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs8,0,fout,"padding, SQCIF[128x96]->QCIF[176x144]",128,96, 176,144  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, SQCIF[128x96]->CIF [352x288]",128,96, 352,288  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, SQCIF[128x96]->QVGA[320x240]",128,96, 320,240  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, SQCIF[128x96]->VGA [640x480]" ,128,96, 640,480 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, QCIF[176x144]->CIF [352x288]",176,144, 352,288 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, QCIF[176x144]->QVGA[320x240]",176,144, 320,240 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, QCIF[176x144]->VGA [640x480]",176,144, 640,480 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, QVGA[320x240]->CIF [352x288]",320,240, 352,288 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, QVGA[320x240]->VGA [640x480]",320,240, 640,480 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"padding, CIF [352x288]->VGA [640x480]",352,288, 640,480 ,1);

    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs16,0,fout,"padding, SQCIF[128x96]->QCIF[176x144]",128,96, 176,144 , 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, SQCIF[128x96]->CIF [352x288]",128,96, 352,288 , 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, SQCIF[128x96]->QVGA[320x240]",128,96, 320,240 , 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, SQCIF[128x96]->VGA [640x480]" ,128,96, 640,480, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, QCIF[176x144]->CIF [352x288]",176,144, 352,288, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, QCIF[176x144]->QVGA[320x240]",176,144, 320,240, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, QCIF[176x144]->VGA [640x480]",176,144, 640,480, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, QVGA[320x240]->CIF [352x288]",320,240, 352,288, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, QVGA[320x240]->VGA [640x480]",320,240, 640,480, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"padding, CIF [352x288]->VGA [640x480]",352,288, 640,480, 2);

    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gu8,1,fout,"padding, SQCIF[128x96]->QCIF[176x144]",128,96, 176,144 , 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, SQCIF[128x96]->CIF [352x288]",128,96, 352,288 , 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, SQCIF[128x96]->QVGA[320x240]",128,96, 320,240 , 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, SQCIF[128x96]->VGA [640x480]" ,128,96, 640,480, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, QCIF[176x144]->CIF [352x288]",176,144, 352,288, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, QCIF[176x144]->QVGA[320x240]",176,144, 320,240, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, QCIF[176x144]->VGA [640x480]",176,144, 640,480, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, QVGA[320x240]->CIF [352x288]",320,240, 352,288, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, QVGA[320x240]->VGA [640x480]",320,240, 640,480, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"padding, CIF [352x288]->VGA [640x480]",352,288, 640,480, 1);

    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs8,1,fout,"padding, SQCIF[128x96]->QCIF[176x144]",128,96, 176,144 , 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, SQCIF[128x96]->CIF [352x288]",128,96, 352,288 , 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, SQCIF[128x96]->QVGA[320x240]",128,96, 320,240 , 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, SQCIF[128x96]->VGA [640x480]" ,128,96, 640,480, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, QCIF[176x144]->CIF [352x288]",176,144, 352,288, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, QCIF[176x144]->QVGA[320x240]",176,144, 320,240, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, QCIF[176x144]->VGA [640x480]",176,144, 640,480, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, QVGA[320x240]->CIF [352x288]",320,240, 352,288, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, QVGA[320x240]->VGA [640x480]",320,240, 640,480, 1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"padding, CIF [352x288]->VGA [640x480]",352,288, 640,480, 1);

    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs16,1,fout,"padding, SQCIF[128x96]->QCIF[176x144]",128,96 , 176,144, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, SQCIF[128x96]->CIF [352x288]",128,96 , 352,288, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, SQCIF[128x96]->QVGA[320x240]",128,96 , 320,240, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, SQCIF[128x96]->VGA [640x480]",128,96 , 640,480, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, QCIF[176x144]->CIF [352x288]",176,144, 352,288, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, QCIF[176x144]->QVGA[320x240]",176,144, 320,240, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, QCIF[176x144]->VGA [640x480]",176,144, 640,480, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, QVGA[320x240]->CIF [352x288]",320,240, 352,288, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, QVGA[320x240]->VGA [640x480]",320,240, 640,480, 2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"padding, CIF [352x288]->VGA [640x480]",352,288, 640,480, 2);
    // crop
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gu8,0,fout,"cropping, QCIF[176x144]->SQCIF[128x96]",176,144,128,96 ,  1);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gu8,0,fout,"cropping, CIF [352x288]->SQCIF[128x96]",352,288,128,96 ,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"cropping, QVGA[320x240]->SQCIF[128x96]",320,240,128,96 ,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"cropping, VGA [640x480]->SQCIF[128x96]",640,480,128,96 ,  1);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gu8,0,fout,"cropping, CIF [352x288]->QCIF[176x144]",352,288,176,144,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"cropping, QVGA[320x240]->QCIF[176x144]",320,240,176,144,  1);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gu8,0,fout,"cropping, VGA [640x480]->QCIF[176x144]",640,480,176,144,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"cropping, CIF [352x288]->QVGA[320x240]",352,288,320,240,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"cropping, VGA [640x480]->QVGA[320x240]",640,480,320,240,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gu8,0,fout,"cropping, VGA [640x480]->CIF [352x288]",640,480,352,288,  1);

    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs8,0,fout,"cropping, QCIF[176x144]->SQCIF[128x96]",176,144,128,96 ,  1);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs8,0,fout,"cropping, CIF [352x288]->SQCIF[128x96]",352,288,128,96 ,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"cropping, QVGA[320x240]->SQCIF[128x96]",320,240,128,96 ,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"cropping, VGA [640x480]->SQCIF[128x96]",640,480,128,96 ,  1);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs8,0,fout,"cropping, CIF [352x288]->QCIF[176x144]",352,288,176,144,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"cropping, QVGA[320x240]->QCIF[176x144]",320,240,176,144,  1);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs8,0,fout,"cropping, VGA [640x480]->QCIF[176x144]",640,480,176,144,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"cropping, CIF [352x288]->QVGA[320x240]",352,288,320,240,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"cropping, VGA [640x480]->QVGA[320x240]",640,480,320,240,  1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs8,0,fout,"cropping, VGA [640x480]->CIF [352x288]",640,480,352,288,  1);

    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs16,0,fout,"cropping, QCIF[176x144]->SQCIF[128x96]",176,144,128,96 ,2);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs16,0,fout,"cropping, CIF [352x288]->SQCIF[128x96]",352,288,128,96 ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"cropping, QVGA[320x240]->SQCIF[128x96]",320,240,128,96 ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"cropping, VGA [640x480]->SQCIF[128x96]",640,480,128,96 ,2);
    PROFILE_IMG_PAD(1,     isVerbose,imgpad_gs16,0,fout,"cropping, CIF [352x288]->QCIF[176x144]",352,288,176,144,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"cropping, QVGA[320x240]->QCIF[176x144]",320,240,176,144,2);
    PROFILE_IMG_PAD(1     ,isVerbose,imgpad_gs16,0,fout,"cropping, VGA [640x480]->QCIF[176x144]",640,480,176,144,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"cropping, CIF [352x288]->QVGA[320x240]",352,288,320,240,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"cropping, VGA [640x480]->QVGA[320x240]",640,480,320,240,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgpad_gs16,0,fout,"cropping, VGA [640x480]->CIF [352x288]",640,480,352,288,2);

    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gu8,1,fout,"cropping, QCIF[176x144]->SQCIF[128x96]",176,144,128,96  ,1);
    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gu8,1,fout,"cropping, CIF [352x288]->SQCIF[128x96]",352,288,128,96  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"cropping, QVGA[320x240]->SQCIF[128x96]",320,240,128,96  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"cropping, VGA [640x480]->SQCIF[128x96]",640,480,128,96  ,1);
    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gu8,1,fout,"cropping, CIF [352x288]->QCIF[176x144]",352,288,176,144 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"cropping, QVGA[320x240]->QCIF[176x144]",320,240,176,144 ,1);
    PROFILE_IMG_PAD(1     ,isVerbose,imgfastpad_gu8,1,fout,"cropping, VGA [640x480]->QCIF[176x144]",640,480,176,144 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"cropping, CIF [352x288]->QVGA[320x240]",352,288,320,240 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"cropping, VGA [640x480]->QVGA[320x240]",640,480,320,240 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gu8,1,fout,"cropping, VGA [640x480]->CIF [352x288]",640,480,352,288 ,1);

    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs8,1,fout,"cropping, QCIF[176x144]->SQCIF[128x96]",176,144,128,96  ,1);
    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs8,1,fout,"cropping, CIF [352x288]->SQCIF[128x96]",352,288,128,96  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"cropping, QVGA[320x240]->SQCIF[128x96]",320,240,128,96  ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"cropping, VGA [640x480]->SQCIF[128x96]",640,480,128,96  ,1);
    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs8,1,fout,"cropping, CIF [352x288]->QCIF[176x144]",352,288,176,144 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"cropping, QVGA[320x240]->QCIF[176x144]",320,240,176,144 ,1);
    PROFILE_IMG_PAD(1     ,isVerbose,imgfastpad_gs8,1,fout,"cropping, VGA [640x480]->QCIF[176x144]",640,480,176,144 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"cropping, CIF [352x288]->QVGA[320x240]",352,288,320,240 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"cropping, VGA [640x480]->QVGA[320x240]",640,480,320,240 ,1);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs8,1,fout,"cropping, VGA [640x480]->CIF [352x288]",640,480,352,288 ,1);

    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs16,1,fout,"cropping, QCIF[176x144]->SQCIF[128x96]",176,144,128,96  ,2);
    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs16,1,fout,"cropping, CIF [352x288]->SQCIF[128x96]",352,288,128,96  ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"cropping, QVGA[320x240]->SQCIF[128x96]",320,240,128,96  ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"cropping, VGA [640x480]->SQCIF[128x96]",640,480,128,96  ,2);
    PROFILE_IMG_PAD(1,     isVerbose,imgfastpad_gs16,1,fout,"cropping, CIF [352x288]->QCIF[176x144]",352,288,176,144 ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"cropping, QVGA[320x240]->QCIF[176x144]",320,240,176,144 ,2);
    PROFILE_IMG_PAD(1     ,isVerbose,imgfastpad_gs16,1,fout,"cropping, VGA [640x480]->QCIF[176x144]",640,480,176,144 ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"cropping, CIF [352x288]->QVGA[320x240]",352,288,320,240 ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"cropping, VGA [640x480]->QVGA[320x240]",640,480,320,240 ,2);
    PROFILE_IMG_PAD(isFull,isVerbose,imgfastpad_gs16,1,fout,"cropping, VGA [640x480]->CIF [352x288]",640,480,352,288 ,2);
}

void mips_imgmisc1(int isFull, int isVerbose, FILE * fout)
{
    mips_hist(isFull,isVerbose,fout);
    mips_norm(isFull,isVerbose,fout);
    mips_interleave(isFull,isVerbose,fout);
    mips_convert(isFull,isVerbose,fout);
    mips_pad(isFull,isVerbose,fout);
} /* mips_imgmisc() */
