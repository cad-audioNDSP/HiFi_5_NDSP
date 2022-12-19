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
* Test module for testing cycle performance (Filtering)
*/

#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(fir)
#include "mips.h"

#define BKFIR16X16_PROFILE(cond,verb,N,M)    OBJ_PROFILE_INVERTED(cond,verb,bkfir16x16,  (M),(objinstance_memory, M, 0, mips.inp2.i16),(bkfir16x16  , mips.out2.i16, mips.inp2.i16,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M);
#define BKFIR24X24_PROFILE(cond,verb,N,M)    OBJ_PROFILE_INVERTED(cond,verb,bkfir24x24,  (M),(objinstance_memory, M, 0, mips.inp2.i32),(bkfir24x24  , mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M);
#define BKFIR24X24P_PROFILE(cond,verb,N,M)   OBJ_PROFILE_INVERTED(cond,verb,bkfir24x24p, (M),(objinstance_memory, M, 0, mips.inp2.i32),(bkfir24x24p , mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M);
#define BKFIR32X16_PROFILE(cond,verb,N,M)    OBJ_PROFILE_INVERTED(cond,verb,bkfir32x16,  (M),(objinstance_memory, M, 0, mips.inp2.i16),(bkfir32x16  , mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M);
#define BKFIR32X32_PROFILE(cond,verb,N,M)    OBJ_PROFILE_INVERTED(cond,verb,bkfir32x32,  (M),(objinstance_memory, M, 0, mips.inp2.i32),(bkfir32x32  , mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M);
#define BKFIR32X32EP_PROFILE(cond,verb,N,M)  OBJ_PROFILE_INVERTED(cond,verb,bkfir32x32ep,(M),(objinstance_memory, M, 0, mips.inp2.i32),(bkfir32x32ep, mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M);

#define BKFIRA16X16_PROFILE(cond,verb,N,M)   OBJ_PROFILE_INVERTED(cond,verb,bkfira16x16,  (M),(objinstance_memory, M, mips.inp2.i16),(bkfira16x16,  mips.out2.i16, mips.inp2.i16,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M);
#define BKFIRA32X32_PROFILE(cond,verb,N,M)   OBJ_PROFILE_INVERTED(cond,verb,bkfira32x32,  (M),(objinstance_memory, M, mips.inp2.i32),(bkfira32x32,  mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M);
#define BKFIRA32X32EP_PROFILE(cond,verb,N,M) OBJ_PROFILE_INVERTED(cond,verb,bkfira32x32ep,(M),(objinstance_memory, M, mips.inp2.i32),(bkfira32x32ep,mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M);
#define BKFIRA32X16_PROFILE(cond,verb,N,M)   OBJ_PROFILE_INVERTED(cond,verb,bkfira32x16,  (M),(objinstance_memory, M, mips.inp2.i16),(bkfira32x16,  mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M);
#define BKFIRA24X24_PROFILE(cond,verb,N,M)   OBJ_PROFILE_INVERTED(cond,verb,bkfira24x24,  (M),(objinstance_memory, M, mips.inp2.i32),(bkfira24x24,  mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M);

#define CXFIR16X16_PROFILE(cond,verb,N,M)    OBJ_PROFILE_INVERTED(cond,verb,cxfir16x16,  (M),(objinstance_memory, M, 0, mips.inp2.ci16),(cxfir16x16  , mips.out2.ci16, mips.inp2.ci16,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M*4);
#define CXFIR32X16_PROFILE(cond,verb,N,M)    OBJ_PROFILE_INVERTED(cond,verb,cxfir32x16,  (M),(objinstance_memory, M, 0, mips.inp2.ci16),(cxfir32x16  , mips.out2.ci32, mips.inp2.ci32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M*4);
#define CXFIR24X24_PROFILE(cond,verb,N,M)    OBJ_PROFILE_INVERTED(cond,verb,cxfir24x24,  (M),(objinstance_memory, M, 0, mips.inp2.ci32),(cxfir24x24  , mips.out2.ci32, mips.inp2.ci32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M*4);
#define CXFIR32X32_PROFILE(cond,verb,N,M)    OBJ_PROFILE_INVERTED(cond,verb,cxfir32x32,  (M),(objinstance_memory, M, 0, mips.inp2.ci32),(cxfir32x32  , mips.out2.ci32, mips.inp2.ci32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M*4);
#define CXFIR32X32EP_PROFILE(cond,verb,N,M)  OBJ_PROFILE_INVERTED(cond,verb,cxfir32x32ep,(M),(objinstance_memory, M, 0, mips.inp2.ci32),(cxfir32x32ep, mips.out2.ci32, mips.inp2.ci32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M*4);

#define STEREO_BKFIR16X16_PROFILE(cond,verb,N,M)    OBJ_PROFILE_INVERTED(cond,verb,stereo_bkfir16x16,  (M),(objinstance_memory, M, 0, mips.inp2.i16, mips.inp2.i16),(stereo_bkfir16x16  , mips.out2.i16, mips.inp2.i16,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M*2);
#define STEREO_BKFIR32X32_PROFILE(cond,verb,N,M)    OBJ_PROFILE_INVERTED(cond,verb,stereo_bkfir32x32,  (M),(objinstance_memory, M, 0, mips.inp2.i32, mips.inp2.i32),(stereo_bkfir32x32  , mips.out2.i32, mips.inp2.i32,N),fout,"N: " #N "; M: " #M ,prf_maccycle,N*M*2);

void mips_bkfir16x16(int isFull, int isVerbose, FILE * fout)
{
    BKFIR16X16_PROFILE(     1, isVerbose,  80, 256)
    BKFIR16X16_PROFILE(isFull, isVerbose,2048,   8)
    BKFIR16X16_PROFILE(isFull, isVerbose, 160,   8)
    BKFIR16X16_PROFILE(isFull, isVerbose, 160,  16)
    BKFIR16X16_PROFILE(isFull, isVerbose,1024,  32)

    BKFIRA16X16_PROFILE(     1, isVerbose,  80, 256)
    BKFIRA16X16_PROFILE(isFull, isVerbose,2048,   8)
    BKFIRA16X16_PROFILE(isFull, isVerbose, 160,   8)
    BKFIRA16X16_PROFILE(isFull, isVerbose, 160,  16)
    BKFIRA16X16_PROFILE(isFull, isVerbose,1024,  32)
}
void mips_bkfir24x24(int isFull, int isVerbose, FILE * fout)
{
#if 0// for HiFi3/3z
    BKFIR24X24_PROFILE(     1, isVerbose,  80, 256)
    BKFIR24X24_PROFILE(isFull, isVerbose,2048,   8)
    BKFIR24X24_PROFILE(isFull, isVerbose, 160,   8)
    BKFIR24X24_PROFILE(isFull, isVerbose, 160,  16)
    BKFIR24X24_PROFILE(isFull, isVerbose,1024,  32)
#endif
    BKFIR24X24P_PROFILE(     1, isVerbose,  80, 256)
    BKFIR24X24P_PROFILE(isFull, isVerbose,  80, 512)
    BKFIR24X24P_PROFILE(isFull, isVerbose,2048, 4)
    BKFIR24X24P_PROFILE(isFull, isVerbose,2048, 8)
    BKFIR24X24P_PROFILE(isFull, isVerbose, 160, 8)
    BKFIR24X24P_PROFILE(isFull, isVerbose, 160, 16)
    BKFIR24X24P_PROFILE(isFull, isVerbose,  80, 16)
    BKFIR24X24P_PROFILE(isFull, isVerbose, 512, 32)
    BKFIR24X24P_PROFILE(isFull, isVerbose,1024, 32)
#if 0// for HiFi3/3z
    BKFIRA24X24_PROFILE(     1, isVerbose,80, 256)
    BKFIRA24X24_PROFILE(isFull, isVerbose,2048, 8)
    BKFIRA24X24_PROFILE(isFull, isVerbose,160, 8)
    BKFIRA24X24_PROFILE(isFull, isVerbose,160, 16)
    BKFIRA24X24_PROFILE(isFull, isVerbose,1024, 32)
#endif
}

void mips_bkfir32x32(int isFull, int isVerbose, FILE * fout)
{
    BKFIR32X16_PROFILE(     1, isVerbose,80, 256)
    BKFIR32X16_PROFILE(isFull, isVerbose,80, 512)
    BKFIR32X16_PROFILE(isFull, isVerbose,2048, 4)
    BKFIR32X16_PROFILE(isFull, isVerbose,2048, 8)
    BKFIR32X16_PROFILE(isFull, isVerbose,160, 8)
    BKFIR32X16_PROFILE(isFull, isVerbose,160, 16)
    BKFIR32X16_PROFILE(isFull, isVerbose,80 , 16)
    BKFIR32X16_PROFILE(isFull, isVerbose,512, 32)
    BKFIR32X16_PROFILE(isFull, isVerbose,1024, 32)

    BKFIR32X32_PROFILE(     1, isVerbose,80, 256)
    BKFIR32X32_PROFILE(isFull, isVerbose,80  , 512)
    BKFIR32X32_PROFILE(isFull, isVerbose,2048, 4)
    BKFIR32X32_PROFILE(isFull, isVerbose,2048, 8)
    BKFIR32X32_PROFILE(isFull, isVerbose,160, 8)
    BKFIR32X32_PROFILE(isFull, isVerbose,160, 16)
    BKFIR32X32_PROFILE(isFull, isVerbose, 80, 16)
    BKFIR32X32_PROFILE(isFull, isVerbose,512, 32)
    BKFIR32X32_PROFILE(isFull, isVerbose,1024, 32)

    BKFIR32X32EP_PROFILE(     1, isVerbose,80, 256)
    BKFIR32X32EP_PROFILE(isFull, isVerbose,80  ,512)
    BKFIR32X32EP_PROFILE(isFull, isVerbose,2048,4  )
    BKFIR32X32EP_PROFILE(isFull, isVerbose,2048, 8)
    BKFIR32X32EP_PROFILE(isFull, isVerbose,160, 8)
    BKFIR32X32EP_PROFILE(isFull, isVerbose,160, 16)
    BKFIR32X32EP_PROFILE(isFull, isVerbose,80 , 16)
    BKFIR32X32EP_PROFILE(isFull, isVerbose,512, 32)
    BKFIR32X32EP_PROFILE(isFull, isVerbose,1024, 32)

    BKFIRA32X16_PROFILE(     1, isVerbose, 80, 256)
    BKFIRA32X16_PROFILE(isFull, isVerbose,80  ,512)
    BKFIRA32X16_PROFILE(isFull, isVerbose,2048,4  )
    BKFIRA32X16_PROFILE(isFull, isVerbose,2048,   8)
    BKFIRA32X16_PROFILE(isFull, isVerbose, 160,   8)
    BKFIRA32X16_PROFILE(isFull, isVerbose, 160,  16)
    BKFIRA32X16_PROFILE(isFull, isVerbose, 80 ,16)
    BKFIRA32X16_PROFILE(isFull, isVerbose, 512,32)
    BKFIRA32X16_PROFILE(isFull, isVerbose,1024,  32)

    BKFIRA32X32_PROFILE(     1, isVerbose, 80, 256)
    BKFIRA32X32_PROFILE(isFull, isVerbose,80  ,512)
    BKFIRA32X32_PROFILE(isFull, isVerbose,2048,4  )
    BKFIRA32X32_PROFILE(isFull, isVerbose,2048,   8)
    BKFIRA32X32_PROFILE(isFull, isVerbose, 160,   8)
    BKFIRA32X32_PROFILE(isFull, isVerbose, 160,  16)
    BKFIRA32X32_PROFILE(isFull, isVerbose, 80 , 16)
    BKFIRA32X32_PROFILE(isFull, isVerbose, 512, 32)
    BKFIRA32X32_PROFILE(isFull, isVerbose,1024,  32)
    
    BKFIRA32X32EP_PROFILE(     1, isVerbose, 80, 256)
    BKFIRA32X32EP_PROFILE(isFull, isVerbose,80  ,512)
    BKFIRA32X32EP_PROFILE(isFull, isVerbose,2048,4  )
    BKFIRA32X32EP_PROFILE(isFull, isVerbose,2048,   8)
    BKFIRA32X32EP_PROFILE(isFull, isVerbose, 160,   8)
    BKFIRA32X32EP_PROFILE(isFull, isVerbose, 160,  16)
    BKFIRA32X32EP_PROFILE(isFull, isVerbose, 80 ,  16)
    BKFIRA32X32EP_PROFILE(isFull, isVerbose, 512,  32)
    BKFIRA32X32EP_PROFILE(isFull, isVerbose,1024,  32)
}

void mips_cxfir(int isFull, int isVerbose, FILE * fout)
{
    CXFIR16X16_PROFILE(     1, isVerbose, 80, 128)
    CXFIR16X16_PROFILE(isFull, isVerbose,2048,   8)
    CXFIR16X16_PROFILE(isFull, isVerbose, 160,   8)
    CXFIR16X16_PROFILE(isFull, isVerbose, 160,  16)
    CXFIR16X16_PROFILE(isFull, isVerbose,1024,  32)

    CXFIR32X16_PROFILE(     1, isVerbose, 80, 128)
    CXFIR32X16_PROFILE(isFull, isVerbose,80  ,512)
    CXFIR32X16_PROFILE(isFull, isVerbose,2048,4  )
    CXFIR32X16_PROFILE(isFull, isVerbose,2048,   8)
    CXFIR32X16_PROFILE(isFull, isVerbose, 160,   8)
    CXFIR32X16_PROFILE(isFull, isVerbose, 160,  16)
    CXFIR32X16_PROFILE(isFull, isVerbose, 80 ,  16)
    CXFIR32X16_PROFILE(isFull, isVerbose, 512,  32)
    CXFIR32X16_PROFILE(isFull, isVerbose,1024,  32)
#if 0 //HiFi3/3z API
    CXFIR24X24_PROFILE(     1, isVerbose, 80, 128)
    CXFIR24X24_PROFILE(isFull, isVerbose,2048,   8)
    CXFIR24X24_PROFILE(isFull, isVerbose, 160,   8)
    CXFIR24X24_PROFILE(isFull, isVerbose, 160,  16)
    CXFIR24X24_PROFILE(isFull, isVerbose,1024,  32)
#endif
    CXFIR32X32_PROFILE(     1, isVerbose, 80, 128)
    CXFIR32X32_PROFILE(isFull, isVerbose,80  ,512)
    CXFIR32X32_PROFILE(isFull, isVerbose,2048,4  )
    CXFIR32X32_PROFILE(isFull, isVerbose,2048,   8)
    CXFIR32X32_PROFILE(isFull, isVerbose, 160,   8)
    CXFIR32X32_PROFILE(isFull, isVerbose, 160,  16)
    CXFIR32X32_PROFILE(isFull, isVerbose, 80 ,  16)
    CXFIR32X32_PROFILE(isFull, isVerbose, 512,  32)
    CXFIR32X32_PROFILE(isFull, isVerbose,1024,  32)

    CXFIR32X32EP_PROFILE(     1, isVerbose, 80, 128)
    CXFIR32X32EP_PROFILE(isFull, isVerbose,80  ,512 )
    CXFIR32X32EP_PROFILE(isFull, isVerbose,2048,4   )
    CXFIR32X32EP_PROFILE(isFull, isVerbose,2048,   8)
    CXFIR32X32EP_PROFILE(isFull, isVerbose, 160,   8)
    CXFIR32X32EP_PROFILE(isFull, isVerbose, 160,  16)
    CXFIR32X32EP_PROFILE(isFull, isVerbose, 80 ,  16)
    CXFIR32X32EP_PROFILE(isFull, isVerbose, 512,  32)
    CXFIR32X32EP_PROFILE(isFull, isVerbose,1024,  32)
}

void mips_stereo_bkfir16x16(int isFull, int isVerbose, FILE * fout)
{
    STEREO_BKFIR16X16_PROFILE(     1, isVerbose,  80, 256)
    STEREO_BKFIR16X16_PROFILE(isFull, isVerbose,2048,   8)
    STEREO_BKFIR16X16_PROFILE(isFull, isVerbose, 160,   8)
    STEREO_BKFIR16X16_PROFILE(isFull, isVerbose, 160,  16)
    STEREO_BKFIR16X16_PROFILE(isFull, isVerbose,1024,  32)
}
void mips_stereo_bkfir32x32(int isFull, int isVerbose, FILE * fout)
{
    STEREO_BKFIR32X32_PROFILE(     1, isVerbose,80, 256)
    STEREO_BKFIR32X32_PROFILE(isFull, isVerbose,2048, 8)
    STEREO_BKFIR32X32_PROFILE(isFull, isVerbose,160, 8)
    STEREO_BKFIR32X32_PROFILE(isFull, isVerbose,160, 16)
    STEREO_BKFIR32X32_PROFILE(isFull, isVerbose,1024, 32)
}

void mips_firblk1(int isFull, int isVerbose, FILE * fout)
{
    mips_bkfir16x16(isFull, isVerbose, fout);
    mips_bkfir24x24(isFull, isVerbose, fout);
    mips_bkfir32x32(isFull, isVerbose, fout);
    mips_cxfir(isFull, isVerbose, fout);
    mips_stereo_bkfir16x16(isFull, isVerbose, fout);
    mips_stereo_bkfir32x32(isFull, isVerbose, fout);
}
