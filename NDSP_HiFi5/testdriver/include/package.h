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
plug and play package API.
Each package is registering by declaration of an class instance
once linked this instance will be available to use in the main()
for calling from the command line
*/

#ifndef __PACKAGE_H__
#define __PACKAGE_H__
#include <stddef.h>
#include <stdio.h>

/* -----------------------------------------------
global functions to use in main()
------------------------------------------------*/
#ifdef __cplusplus
extern "C"
{
#endif

	typedef enum
	{
		eTEST_FUNC,
		eTEST_MIPS,
		eTEST_ACCURACY,
		eTEST_UNKNOWN = -1
	} eTestType;

	// flags derived from the command line options
	typedef enum
	{
		eTESTFLAG_MIPS		= 0x0800,
		eTESTFLAG_NOABORT	= 0x1000,
		eTESTFLAG_VERBOSE	= 0x2000,
		eTESTFLAG_VERIFY	= 0x0080,
		eTESTFLAG_FUNC		= 0x0040,
		eTESTFLAG_FULL		= 0x4000,
		eTESTFLAG_BRIEF		= 0x0020,
		eTESTFLAG_SANITY	= 0x0002,
		eTESTFLAG_NOWARMUP	= 0x0010,
		eTESTFLAG_ACCURACY	= 0x0008,
		eTESTFLAG_EXCEPTION	= 0x0004,
		eTESTFLAG_HELP		= 0x8000
	} eTestFlags;

	int  PackageProcessOpt(const char * cmd);                           // process command line option, returns 0 on error
	eTestType PackageGetTestType();
	eTestFlags PackageGetTestFlags();
	// start all tests
	int  PackagePerformTest(int(*perf_info)(FILE * f, const char * fmt, ...));

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
/* -----------------------------------------------
just class declaration for c++ stubs
------------------------------------------------*/
typedef  int(*func_test)(int isFull, int isVerbose, int breakOnError);
typedef void(*mips_test)(int isFull, int isVerbose, FILE* f);
typedef  int( *acc_test)(int isFull, int isVerbose, int breakOnError, int optAccuracy, int optException);

class Package
{
public:
	Package(
		const char *cmd,     /* command line string to call                      */
		int group,           /* package group for which it belongs to            */
		int number,          /* number inside the group which controls the order
							 zero indicates the declaration of the main item
							 in the group                                     */
							 int phase,           /* phase number. Controls the datatype supported
												  in the routines in this specific package         */
												  const char* help,    /* help string                                      */
												  func_test   func,    /* functional testing (can be NULL)                 */
												  mips_test   mips,    /* cycle testing (can be NULL)                      */
												  acc_test    acc      /* accuracy testing (can be NULL)                   */
												  );
protected:
private:
	friend class PackageContainer;
	friend int  PackageProcessOpt(const char * cmd);
	friend int  PackagePerformTest(int(*perf_info)(FILE * f, const char * fmt, ...));
	static Package* first;
	int compare(const Package* ref)
	{
		if (group != ref->group) return group >ref->group ? 1 : -1;
		if (number != ref->number) return number > ref->number ? 1 : -1;
		if (phase != ref->phase) return phase > ref->phase ? 1 : -1;
		return 0;
	}
	void insert(Package* c); // insert package to the list in the proper order
	Package* prev;
	Package* next;

	/* copies of original data from initialization */
	const char *cmd;     /* command line string to call                      */
	int group;           /* package group for which it belongs to            */
	int number;          /* number inside the group which controls the order
						 zero indicates the declaration of the main item
						 in the group                                     */
	int phase;           /* phase number. Controls the datatype supported
						 in the routines in this specific package         */
	const char* help;    /* help string                                      */
	func_test   func;    /* functional testing (can be NULL)                 */
	mips_test   mips;    /* cycle testing (can be NULL)                      */
	acc_test    acc;     /* accuracy testing (can be NULL)                   */
	// operational variables
	int execFlags;  /* flags set by PackageProcessOpt() to signify execution of this module */
};
#endif

#endif /* __PACKAGE_H__ */
