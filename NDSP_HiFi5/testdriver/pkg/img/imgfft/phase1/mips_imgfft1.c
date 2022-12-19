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

#define PROFILE_IMG_FFT(isFull,isVerbose,basename,_file,resolution,w,h) { \
        size_t szScr;                                                     \
        size_t bytesOut;                                                  \
        void* pScr;                                                       \
        void* pImg;                                                       \
        void* pOutImg;                                                    \
        int istride = (w+7)&~7;                                           \
        imgsize_t sz;                                                     \
        sz.width=w;                                                       \
        sz.height=h;                                                      \
        sz.stride=istride;                                                \
        szScr=basename##_getScratchSize (&sz );                           \
        szScr=(szScr+15)&~15;                                             \
        /* bytesOut=h*w*sizeof(complex_fract16); */                       \
        bytesOut=h*w*sizeof(fract16);                             		  \
        NASSERT(szScr+3 * bytesOut<=sizeof(mips));                        \
        (void)(bytesOut);                                                 \
        pScr=(void*)&mips;                                                \
        pImg=(void*)(((uintptr_t)pScr)+szScr);                            \
        pOutImg=(void*)(((uintptr_t)pImg)+bytesOut);                      \
        PROFILE_SIMPLE(isFull,isVerbose,basename,                         \
                      /* (pScr,(complex_fract16*)pImg,pImg,128<<16,&sz), */ \
                       (pScr,(complex_fract16*)pOutImg,pImg,128<<16,&sz),   \
                       _file,resolution,prf_cycle);                       	\
}

#define PROFILE_IMG_IFFT(isFull,isVerbose,basename,_file,resolution,w,h) {\
        size_t szScr;                                                     \
        size_t bytesIn;                                                   \
        void* pScr;                                                       \
        void* pImg;                                                       \
        void* pOutImg;                                                    \
        int istride = (w+7)&~7;                                           \
        imgsize_t sz;                                                     \
        sz.width=w;                                                       \
        sz.height=h;                                                      \
        sz.stride=istride;                                                \
        szScr=basename##_getScratchSize (&sz );                           \
        szScr=(szScr+15)&~15;                                             \
        /*bytesIn=h*w*sizeof(complex_fract16);  */                        \
        bytesIn=h*w*sizeof(fract16);                              		  \
        (void)(bytesIn);                                                  \
        NASSERT(szScr+3 * bytesIn<=sizeof(mips));                         \
        pScr=(void*)&mips;                                                \
        pImg=(void*)(((uintptr_t)pScr)+szScr);                            \
        pOutImg=(void*)(((uintptr_t)pImg)+2 * bytesIn);                   \
        PROFILE_SIMPLE(isFull,isVerbose,basename,                         \
                      /* (pScr, pImg, (complex_fract16*)pImg,128<<16,&sz,0), */ \
                       (pScr, pOutImg, (complex_fract16*)pImg,128<<16,&sz,0),	\
                       _file,resolution,prf_cycle);                       		\
}

void mips_imgfft1(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_IMG_FFT(1,     isVerbose,imgfft_gu8,fout,"64x64",   64,  64);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gu8,fout,"128x128",128, 128);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gu8,fout,"256x256",256, 256);
    PROFILE_IMG_FFT(1,     isVerbose,imgfft_gu8,fout,"512x512",512, 512);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gu8,fout,"SQCIF[128x96]", 128, 96);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gu8,fout,"QCIF[176x144]", 176,144);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gu8,fout,"QVGA[320x240]", 320,240);
    PROFILE_IMG_FFT(1     ,isVerbose,imgfft_gu8,fout,"CIF [352x288]", 352,288);
    PROFILE_IMG_FFT(1     ,isVerbose,imgfft_gu8,fout,"VGA [640x480]", 640,480);

    PROFILE_IMG_FFT(1,     isVerbose,imgfft_gs8,fout,"64x64",   64,  64);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs8,fout,"128x128",128, 128);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs8,fout,"256x256",256, 256);
    PROFILE_IMG_FFT(1,     isVerbose,imgfft_gs8,fout,"512x512",512, 512);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs8,fout,"SQCIF[128x96]", 128, 96);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs8,fout,"QCIF[176x144]", 176,144);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs8,fout,"QVGA[320x240]", 320,240);
    PROFILE_IMG_FFT(1     ,isVerbose,imgfft_gs8,fout,"CIF [352x288]", 352,288);
    PROFILE_IMG_FFT(1     ,isVerbose,imgfft_gs8,fout,"VGA [640x480]", 640,480);

    PROFILE_IMG_FFT(1,     isVerbose,imgfft_gs16,fout,"64x64",   64,  64);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs16,fout,"128x128",128, 128);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs16,fout,"256x256",256, 256);
    PROFILE_IMG_FFT(1,     isVerbose,imgfft_gs16,fout,"512x512",512, 512);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs16,fout,"SQCIF[128x96]", 128, 96);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs16,fout,"QCIF[176x144]", 176,144);
    PROFILE_IMG_FFT(isFull,isVerbose,imgfft_gs16,fout,"QVGA[320x240]", 320,240);
    PROFILE_IMG_FFT(1     ,isVerbose,imgfft_gs16,fout,"CIF [352x288]", 352,288);
    PROFILE_IMG_FFT(1     ,isVerbose,imgfft_gs16,fout,"VGA [640x480]", 640,480);

    PROFILE_IMG_IFFT(1,     isVerbose,imgifft_gu8,fout,"64x64",   64,  64);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gu8,fout,"128x128",128, 128);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gu8,fout,"256x256",256, 256);
    PROFILE_IMG_IFFT(1,     isVerbose,imgifft_gu8,fout,"512x512",512, 512);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gu8,fout,"SQCIF[128x96]", 128, 96);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gu8,fout,"QCIF[176x144]", 176,144);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gu8,fout,"QVGA[320x240]", 320,240);
    PROFILE_IMG_IFFT(1     ,isVerbose,imgifft_gu8,fout,"CIF [352x288]", 352,288);
    PROFILE_IMG_IFFT(1     ,isVerbose,imgifft_gu8,fout,"VGA [640x480]", 640,480);

    PROFILE_IMG_IFFT(1,     isVerbose,imgifft_gs8,fout,"64x64",   64,  64);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs8,fout,"128x128",128, 128);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs8,fout,"256x256",256, 256);
    PROFILE_IMG_IFFT(1,     isVerbose,imgifft_gs8,fout,"512x512",512, 512);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs8,fout,"SQCIF[128x96]", 128, 96);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs8,fout,"QCIF[176x144]", 176,144);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs8,fout,"QVGA[320x240]", 320,240);
    PROFILE_IMG_IFFT(1     ,isVerbose,imgifft_gs8,fout,"CIF [352x288]", 352,288);
    PROFILE_IMG_IFFT(1     ,isVerbose,imgifft_gs8,fout,"VGA [640x480]", 640,480);

    PROFILE_IMG_IFFT(1,     isVerbose,imgifft_gs16,fout,"64x64",   64,  64);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs16,fout,"128x128",128, 128);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs16,fout,"256x256",256, 256);
    PROFILE_IMG_IFFT(1,     isVerbose,imgifft_gs16,fout,"512x512",512, 512);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs16,fout,"SQCIF[128x96]", 128, 96);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs16,fout,"QCIF[176x144]", 176,144);
    PROFILE_IMG_IFFT(isFull,isVerbose,imgifft_gs16,fout,"QVGA[320x240]", 320,240);
    PROFILE_IMG_IFFT(1     ,isVerbose,imgifft_gs16,fout,"CIF [352x288]", 352,288);
    PROFILE_IMG_IFFT(1     ,isVerbose,imgifft_gs16,fout,"VGA [640x480]", 640,480);
}
