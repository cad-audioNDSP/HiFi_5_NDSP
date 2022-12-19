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
* Test procedures for complex FFT functions.
*/
#include "test_rfftie.h"

/* Perform all tests for fft_realMxN_ie, ifft_realMxN_ie API functions. */
int main_rfftie(int phaseNum, int isFull, int isVerbose, int breakOnError, const TestDef_t *tbl)
{
    int ix, res;
    for (ix = 0, res = 1; tbl[ix].seqFilename && (res || !breakOnError); ix++)
    {
        if (phaseNum == 0 || phaseNum == tbl[ix].phaseNum)
        {
            tTestEngTarget target = tbl[ix].target;
            /* Make sure that all functions is present */
            if (!IS_PRESENT(tbl[ix].desc.frwTransFxn) &&
                !IS_PRESENT(tbl[ix].desc.invTransFxn))
            {
                target = (tTestEngTarget)tbl[ix].desc.frwTransFxn;
            }

            res &= (0 != TestEngRun(target,
                &tbl[ix].desc.desc,
                tbl[ix].seqFilename,
                isFull, isVerbose, breakOnError,0));
        }
    }
    return (res);
} /* main_rfftie() */
