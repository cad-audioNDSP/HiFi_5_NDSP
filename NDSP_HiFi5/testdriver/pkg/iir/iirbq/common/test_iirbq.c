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
 * Test procedures for IIR
 */
#include "test_iirbq.h"

/* Perform all tests for IIR API functions. */
int main_iirbq( int phaseNum, int isFull, int isVerbose, int breakOnError, const tTbl * tbl, int szTbl )
{
    int n;
    int res = 1;
    for (n=0; n<szTbl; n++)
    {
        if ( ( phaseNum == 0 || phaseNum == tbl[n].phaseNum ) && ( (isFull&1) || tbl[n].runAlways ) )
        {
            res &= (0!=TestEngRun(tbl[n].fxns, tbl[n].pIirDescr, tbl[n].seqFile, isFull, isVerbose, breakOnError, 0));
            if (res == 0 && breakOnError) break;
        }
    }
    return (res);
} /* main_iirbq() */
