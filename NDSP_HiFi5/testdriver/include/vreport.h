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
 * verification report utilities
 */

#ifndef __VREPORT_H__
#define __VREPORT_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*tReportFUT)() ;
/* initialize verification reporting */
void vReportInit(int turnOn, const char* exeName);

typedef enum
{
	evr_NOTTESTED=0,
	evr_OK=1,
	evr_FAIL=2
}
evrResult;

/* add to statistics 
   Input:
   fun_ptr_t[Nfun]	pointer to the functions under the test (FUT). 
   Nfun             number of FUT (several FUT might be tested together by one test case)
   extraParams	    extra parameters of function (NULL if there is no specialized variants of FUT)
   dataFile         filename with data
   caseType		    test case type 
   dataSize         total size of input/output data for FUT
*/
void vReportAdd( const tReportFUT fun_ptr_t[], int Nfun, const char* extraParams, const char* dataFile, int caseType, size_t dataSize);

/* set results for the lasting run test 
   testResult   result of data testing of FUT
   exceptResult result of exception handling testing of FUT
   memResult    result of testing the memory after the call of FUT

   set evr_NOTTESTED if specific test result is not known at a moment
*/
void vReportSetResult (evrResult testResult, evrResult exceptResult, evrResult memResult);

/* print report */
void vReportPrint();

/* deallocate resources */
void vReportClose();

/*
	just stop reporting by vReportSetResult() until next calls of vReportAdd()
	call this function upon the end of data file
*/
void vReportFinish();

#ifdef __cplusplus
};
#endif
#endif/*__VREPORT_H__*/
