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
  Test module/example for testing NatureDSP_Signal library
  Integrit, 2006-2019
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* Cycles measurement API. */
#include "mips.h"
#include "vreport.h"
#include "package.h"
int main_mips     ( int phaseNum, int isFull, int isVerbose, int breakOnError );

int main(int argc, char** argv)
{
    eTestType testType;
    eTestFlags testFlags;
    int m,retval;
  /*----------------------------------------------------------------------------*
   *                            Show Library info                               *
   *----------------------------------------------------------------------------*/
  {
    char lib_version[32], api_version[32];

    GET_LIBRARY_VERSION    ( lib_version );
    GET_LIBRARY_API_VERSION( api_version );

    printf( "%s library version: %s\n"
            "%s API version:     %s\n",
            LIBRARY_NAME, lib_version, LIBRARY_NAME, api_version );
  }

    for (m=1; m<argc; m++)
    {
        if (PackageProcessOpt(argv[m])==0)
        {
          printf( "Invalid or not supported command line option: %s\n", argv[m] );
          return (-1);
        }
    }
    testType=PackageGetTestType();
    testFlags =PackageGetTestFlags();
    // cache warnimg support for profiler
    noWarmup=(testFlags & eTESTFLAG_NOWARMUP)?1:0;
    #if defined (COMPILER_XTENSA) && (MEM_MODEL!=2)
    if (noWarmup)
    {
        printf("-nowarmup option is ignored. It should be used if code is built with MEM_MODEL=2\n");
    }
    #endif

    switch(testType)
    {
    case eTEST_ACCURACY:break;
    case eTEST_MIPS:        /* Initialize mips testing. */
            main_mips( 0, 0, 0, 0 );
            break;
    case eTEST_FUNC:
            vReportInit( (testFlags & eTESTFLAG_VERIFY) ? 1:0,argv[0] );
            break;
    default: break;
    }
    retval=PackagePerformTest(perf_info);
    switch(testType)
    {
    case eTEST_ACCURACY:break;
    case eTEST_MIPS:    break;
    case eTEST_FUNC:    
            vReportPrint( );
            vReportClose( );
            break;
    default: break;
    }
    return retval;
}
