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
* Test module for testing cycle performance (Correlation, Convolution, 
* Dispreading, LMS)
*/

#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(fir)
#include "mips.h"

#define RPPAF_PROFILE(cond,verb,M,N)        OBJ_PROFILE_INVERTED(cond,verb,rppafirf,      (M,N),(objinstance_memory, M, N, mips.inp2.f32),(rppafirf  , mips.out2.f32, mips.inp2.f32),fout,"M: " #M "; N: " #N ,prf_maccycle,N*M);
#define RPPSF_PROFILE(cond,verb,M,N)        OBJ_PROFILE_INVERTED(cond,verb,rppsfirf,      (M,N),(objinstance_memory, M, N, mips.inp2.f32),(rppsfirf  , mips.out2.f32, mips.inp2.f32),fout,"M: " #M "; N: " #N ,prf_maccycle,N*M);
#define CPPAF_PROFILE(cond,verb,M,N)        OBJ_PROFILE_INVERTED(cond,verb,cppafirf,      (M,N),(objinstance_memory, M, N, mips.inp2.f32),(cppafirf  , mips.out2.cf32, mips.inp2.cf32),fout,"M: " #M "; N: " #N ,prf_maccycle,2*N*M);
#define CPPSF_PROFILE(cond,verb,M,N)        OBJ_PROFILE_INVERTED(cond,verb,cppsfirf,      (M,N),(objinstance_memory, M, N, mips.inp2.f32),(cppsfirf  , mips.out2.cf32, mips.inp2.cf32),fout,"M: " #M "; N: " #N ,prf_maccycle,2*N*M);

void mips_convolf(int isFull, int isVerbose, FILE * fout)
{
    void* pScr = (void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_convolf,(mips.out0.f32, mips.inp2.f32, mips.inp1.f32,    80, 56),fout,"N: 80; M: 56",prf_maccycle,     80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convolf,(mips.out0.f32, mips.inp2.f32, mips.inp1.f32,   256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );

    PROFILE_INVERTED(isFull, isVerbose,fir_convolaf,(pScr, mips.out0.f32, mips.inp2.f32, mips.inp1.f32,    80, 56),fout,"N: 80; M: 56",prf_maccycle,     80*56 );
    PROFILE_INVERTED(     1, isVerbose,fir_convolaf,(pScr, mips.out0.f32, mips.inp2.f32, mips.inp1.f32,   256, 80),fout,"N: 256; M: 80",prf_maccycle,   80*256 );
}

void mips_xcorrf(int isFull, int isVerbose, FILE * fout)
{
    void* pScr = (void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_xcorrf,  (mips.out0.f32, mips.inp2.f32, mips.inp1.f32, 80, 56    ),fout,"N: 80; M: 56" ,prf_maccycle,    80*56);
    PROFILE_INVERTED(     1, isVerbose,fir_xcorrf,  (mips.out0.f32, mips.inp2.f32, mips.inp1.f32, 256, 80   ),fout,"N: 256; M: 80",prf_maccycle,   80*256);
    PROFILE_INVERTED(isFull, isVerbose,cxfir_xcorrf,(mips.out0.cf32, mips.inp2.cf32, mips.inp1.cf32, 80, 56 ),fout,"N: 80; M: 56" ,prf_maccycle,  4*80*56);
    PROFILE_INVERTED(     1, isVerbose,cxfir_xcorrf,(mips.out0.cf32, mips.inp2.cf32, mips.inp1.cf32, 256, 80),fout,"N: 256; M: 80",prf_maccycle, 4*80*256);

    PROFILE_INVERTED(isFull, isVerbose,fir_xcorraf,  (pScr, mips.out0.f32,  mips.inp2.f32,  mips.inp1.f32, 80, 56  ),fout,"N: 80; M: 56" ,prf_maccycle,    80*56);
    PROFILE_INVERTED(     1, isVerbose,fir_xcorraf,  (pScr, mips.out0.f32,  mips.inp2.f32,  mips.inp1.f32, 256, 80 ),fout,"N: 256; M: 80",prf_maccycle,   80*256);
    PROFILE_INVERTED(isFull, isVerbose,cxfir_xcorraf,(pScr, mips.out0.cf32, mips.inp2.cf32, mips.inp1.cf32, 80, 56 ),fout,"N: 80; M: 56" ,prf_maccycle,  4*80*56);
    PROFILE_INVERTED(     1, isVerbose,cxfir_xcorraf,(pScr, mips.out0.cf32, mips.inp2.cf32, mips.inp1.cf32, 256, 80),fout,"N: 256; M: 80",prf_maccycle, 4*80*256);
}

void mips_acorrf(int isFull, int isVerbose, FILE * fout)
{
    void* pScr = (void*)mips.scratch0.u8;
    PROFILE_INVERTED(isFull, isVerbose,fir_acorrf,(mips.out0.f32, mips.inp2.f32, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorrf,(mips.out0.f32, mips.inp2.f32, 256),fout,"N: 256",prf_maccycle,  256*256 );
    PROFILE_INVERTED(isFull, isVerbose,fir_acorraf,(pScr, mips.out0.f32, mips.inp2.f32, 80),fout,"N: 80",prf_maccycle,    80*80   );
    PROFILE_INVERTED(     1, isVerbose,fir_acorraf,(pScr, mips.out0.f32, mips.inp2.f32, 256),fout,"N: 256",prf_maccycle,  256*256 );
}

void mips_gccphatf(int isFull, int isVerbose, FILE * fout)
{
    void* pScr=(void*)mips.scratch0.u8;
    PROFILE_INVERTED(     1, isVerbose,gccphatf,(pScr, mips.out1.f32, mips.inp0.f32, mips.inp1.f32, 64),fout,"N=64" ,prf_ptscycle3, 64);
    PROFILE_INVERTED(isFull, isVerbose,gccphatf,(pScr, mips.out1.f32, mips.inp0.f32, mips.inp1.f32,128),fout,"N=160",prf_ptscycle3,160);
    PROFILE_INVERTED(isFull, isVerbose,gccphatf,(pScr, mips.out1.f32, mips.inp0.f32, mips.inp1.f32,256),fout,"N=256",prf_ptscycle3,256);
    PROFILE_INVERTED(     1, isVerbose,gccphatf,(pScr, mips.out1.f32, mips.inp0.f32, mips.inp1.f32,320),fout,"N=320",prf_ptscycle3,320);
}

void mips_blmsf(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_INVERTED(isFull, isVerbose,fir_blmsf,(mips.out1.f32, mips.out2.f32, mips.inp1.f32, mips.inp0.f32, 0.1f, 1.1f, 80, 16),fout,"N: 80; M: 16",prf_maccycle,    80*16*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blmsf,(mips.out1.f32, mips.out2.f32, mips.inp1.f32, mips.inp0.f32, 0.1f, 1.1f, 64, 16),fout,"N: 64; M: 16",prf_maccycle,    64*16*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blmsf,(mips.out1.f32, mips.out2.f32, mips.inp1.f32, mips.inp0.f32, 0.1f, 1.1f, 64, 64),fout,"N: 64; M: 64",prf_maccycle,    64*64*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blmsf,(mips.out1.f32, mips.out2.f32, mips.inp1.f32, mips.inp0.f32, 0.1f, 1.1f, 80, 64),fout,"N: 80; M: 64",prf_maccycle,    80*64*2 );
    PROFILE_INVERTED(     1, isVerbose,fir_blmsf,(mips.out1.f32, mips.out2.f32, mips.inp2.f32, mips.inp0.f32, 0.1f, 1.1f, 80, 128),fout,"N: 80; M: 128",prf_maccycle, 80*128*2 );
    PROFILE_INVERTED(isFull, isVerbose,fir_blmsf,(mips.out1.f32, mips.out2.f32, mips.inp2.f32, mips.inp0.f32, 0.1f, 1.1f, 64, 128),fout,"N: 64; M: 128",prf_maccycle, 64*128*2 );

    PROFILE_INVERTED(isFull, isVerbose,cxfir_blmsf,(mips.out1.cf32, mips.out2.cf32, mips.inp1.cf32, mips.inp0.cf32, 0.1f, 1.1f, 80, 16),fout,"N: 80; M: 16",prf_maccycle,    80*16*2*4 );
    PROFILE_INVERTED(isFull, isVerbose,cxfir_blmsf,(mips.out1.cf32, mips.out2.cf32, mips.inp1.cf32, mips.inp0.cf32, 0.1f, 1.1f, 64, 16),fout,"N: 64; M: 16",prf_maccycle,    64*16*2*4 );
    PROFILE_INVERTED(isFull, isVerbose,cxfir_blmsf,(mips.out1.cf32, mips.out2.cf32, mips.inp1.cf32, mips.inp0.cf32, 0.1f, 1.1f, 64, 64),fout,"N: 64; M: 64",prf_maccycle,    64*64*2*4 );
    PROFILE_INVERTED(isFull, isVerbose,cxfir_blmsf,(mips.out1.cf32, mips.out2.cf32, mips.inp1.cf32, mips.inp0.cf32, 0.1f, 1.1f, 80, 64),fout,"N: 80; M: 64",prf_maccycle,    80*64*2*4 );
    PROFILE_INVERTED(     1, isVerbose,cxfir_blmsf,(mips.out1.cf32, mips.out2.cf32, mips.inp2.cf32, mips.inp0.cf32, 0.1f, 1.1f, 80, 128),fout,"N: 80; M: 128",prf_maccycle, 80*128*2*4 );
    PROFILE_INVERTED(isFull, isVerbose,cxfir_blmsf,(mips.out1.cf32, mips.out2.cf32, mips.inp2.cf32, mips.inp0.cf32, 0.1f, 1.1f, 64, 128),fout,"N: 64; M: 128",prf_maccycle, 64*128*2*4 );
}

void mips_ppfirf(int isFull, int isVerbose, FILE * fout)
{
    RPPAF_PROFILE(isFull, isVerbose,  32, 24);
    RPPAF_PROFILE(     1, isVerbose,  64, 20);
    RPPAF_PROFILE(isFull, isVerbose,  96, 16);
    RPPAF_PROFILE(isFull, isVerbose, 128, 12);

    RPPSF_PROFILE(isFull, isVerbose,  32, 24);
    RPPSF_PROFILE(     1, isVerbose,  64, 20);
    RPPSF_PROFILE(isFull, isVerbose,  96, 16);
    RPPSF_PROFILE(isFull, isVerbose, 128, 12);

    CPPAF_PROFILE(isFull, isVerbose,  32, 24);
    CPPAF_PROFILE(     1, isVerbose,  64, 20);
    CPPAF_PROFILE(isFull, isVerbose,  96, 16);
    CPPAF_PROFILE(isFull, isVerbose, 128, 12);

    CPPSF_PROFILE(isFull, isVerbose,  32, 24);
    CPPSF_PROFILE(     1, isVerbose,  64, 20);
    CPPSF_PROFILE(isFull, isVerbose,  96, 16);
    CPPSF_PROFILE(isFull, isVerbose, 128, 12);
}

void mips_firother2(int isFull, int isVerbose, FILE * fout)
{
    mips_convolf(isFull, isVerbose, fout);
    mips_xcorrf(isFull, isVerbose, fout);
    mips_acorrf(isFull, isVerbose, fout);
    mips_gccphatf(isFull, isVerbose, fout);
    mips_blmsf(isFull, isVerbose, fout);
    mips_ppfirf(isFull, isVerbose, fout);
}
