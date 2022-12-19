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
 * Test procedures for FIR
 */
#include "test_firother.h"
/*-----------------------------------------------
wrapper for LMS to test convergence:
input:
x[M+N*P-1]  far end
r[N*P]      near end
norm[P]     norm factor for P blocks
mu          adapt rate
M           IR length
N           block size
P           number of blocks
output:
h[M]        estimated impulse response
temporary
e[N]        error
-----------------------------------------------*/
void fir_blmsf_convergence( 
                      float32_t * e, float32_t * h, const float32_t * r,
                const float32_t * x, const float32_t* norm, float32_t mu, int N, int M, int P )
{
    int m,p;
    for(m=0; m<M; m++) h[m]=0;
    for (p=0; p<P; p++,x+=N,r+=N,norm++)
    {
        fir_blmsf(e,h,r,x,*norm,mu,N,M);
    }
}

void cxfir_blmsf_convergence( complex_float * e, complex_float * h, 
                const complex_float * r,
                const complex_float * x, 
                const float32_t* norm, float32_t mu, 
                int          N, int       M , int P)
{
    int m,p;
    for(m=0; m<2*M; m++) ((float32_t*)h)[m]=0;
    for (p=0; p<P; p++,x+=N,r+=N,norm++)
    {
        cxfir_blmsf(e,h,r,x,*norm,mu,N,M);
    }
}

void fir_blms16x16_convergence (  int16_t * e, int16_t *  h,
                const int16_t * r,
                const int16_t * x,
                const int16_t * norm,int16_t   mu,
                int       N,   int       M, int P)
{
    int m,p;
    for(m=0; m<M; m++) h[m]=0;
    for (p=0; p<P; p++,x+=N,r+=N,norm++)
    {
        fir_blms16x16(e,h,r,x,*norm,mu,N,M);
    }
}

void fir_blms16x32_convergence (  int32_t * e, int32_t *  h,
                const int16_t * r,
                const int16_t * x,
                const int32_t * norm,int16_t   mu,
                int       N,   int       M, int P)
{
    int m,p;
    for(m=0; m<M; m++) h[m]=0;
    for (p=0; p<P; p++,x+=N,r+=N,norm++)
    {
        fir_blms16x32(e,h,r,x,*norm,mu,N,M);
    }
}

void fir_blms32x32_convergence(  int32_t * e, int32_t *  h,
                const int32_t * r,
                const int32_t * x,
                const int32_t * norm, int32_t mu,
                int       N,   int       M, int P)
{
    int m,p;
    for(m=0; m<M; m++) h[m]=0;
    for (p=0; p<P; p++,x+=N,r+=N,norm++)
    {
        fir_blms32x32(e,h,r,x,*norm,mu,N,M);
    }
}

void fir_blms32x32ep_convergence(  int32_t * e, int32_t *  h,
                const int32_t * r,
                const int32_t * x,
                const int32_t * norm, int32_t mu,
                int       N,   int       M, int P)
{
    int m,p;
    for(m=0; m<M; m++) h[m]=0;
    for (p=0; p<P; p++,x+=N,r+=N,norm++)
    {
        fir_blms32x32ep(e,h,r,x,*norm,mu,N,M);
    }
}

void cxfir_blms32x32_convergence (complex_fract32 *  e, complex_fract32 *  h,
                const complex_fract32 *  r,
                const complex_fract32 *  x,
                const int32_t *  norm, int32_t mu,
                int       N,   int       M, int P)
{
    int m,p;
    for(m=0; m<2*M; m++) ((int32_t*)h)[m]=0;
    for (p=0; p<P; p++,x+=N,r+=N,norm++)
    {
        cxfir_blms32x32(e,h,r,x,*norm,mu,N,M);
    }
}

/* check if target function is not available in the target configuration already */
int te_create_convergence(tTestEngContext * context)
{
    return IS_PRESENT(context->desc->extraPtr) ? 1:-1;
}

int test_firother( int phaseNum, int isFull, int isVerbose, int breakOnError, const tTbl * tbl, int szTbl )
{
    int res = 1;

    int n;
    for (n=0; n<szTbl; n++)
    {
        if ( ( phaseNum == 0 || phaseNum == tbl[n].phaseNum ) && ( (isFull&1) || tbl[n].runAlways ) )
        {
            res &= (0!=TestEngRun(tbl[n].fxns, tbl[n].pFirDescr, tbl[n].seqFile, isFull, isVerbose, breakOnError,0));
            if (res == 0 && breakOnError) break;
        }
    }
    
    return res;
}
