/* ------------------------------------------------------------------------ */
/* Copyright (c) 2020 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
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
/*          Copyright (C) 2009-2020 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
plug and play package API.
Each package is registering by declaration of an class instance
once linked this instance will be available to use in the main()
for calling from the command line
*/
#include <stddef.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "package.h"

class PackageContainer;

Package* Package::first = NULL;

Package::Package(
	const char *cmd,     /* command line string to call                      */
	int group,           /* package group for which it belongs to            */
	int number,          /* number inside the group ehich controls the order
						 zero indicates the declaration of the main item
						 in the group                                     */
						 int phase,           /* phase number. Controls the datatype supported
											  in the routines in this specific package         */
											  const char* help,    /* help string                                      */
											  func_test   func,    /* functional testing (can be NULL)                 */
											  mips_test   mips,    /* cycle testing (can be NULL)                      */
											  acc_test    acc      /* accuracy testing (can be NULL)                   */
											  )
{
	this->cmd = cmd;
	this->group = group;
	this->number = number;
	this->phase = phase;
	this->help = help;
	this->func = func;
	this->mips = mips;
	this->acc = acc;
	insert(this);
}

/*
inserts package to the list in acsending order:
group
number
phase
*/
void Package::insert(Package* c)
{
	Package* prev, *oldprev = NULL;
	if (first == NULL)
	{   // first package
		first = c;
		c->prev = NULL;
		c->next = NULL;
		return;
	}
	// find a point for insertion
	prev = first;
	if (c->compare(first)<0)
	{
		// insert before the first
		c->prev = NULL;
		c->next = first;
		first->prev = c;
		first = c;
		return;
	}
	while (prev)
	{
		if (c->compare(prev)<0)
		{
			oldprev = prev->prev;
			// insert before prev
			assert(oldprev != NULL);
			oldprev->next = c;
			prev->prev = c;
			c->prev = oldprev;
			c->next = prev;
			return;
		}
		oldprev = prev;
		prev = prev->next;
	}
	// insert at the end
	oldprev->next = c;
	c->prev = oldprev;
	c->next = NULL;
}


/*
Package container:

Provides all necessary functions for searching inside the packages
and calls functional/cycle performance utilities according to the
command line
*/
class PackageContainer
{
public:
	PackageContainer() { phaseFlags = 0; flags = eTESTFLAG_MIPS | eTESTFLAG_SANITY; };
	static void printHelp();
	int phaseFlags;
	int flags;
private:
};
static PackageContainer myPackageContainer;

/*------------------------------------------------------
Print help content for all plugged packages
------------------------------------------------------*/
void PackageContainer::printHelp()
{
	Package* next = Package::first;
	while (next)
	{
		printf("%-8s  %s\n", next->cmd, next->help);
		next = next->next;
	}
}


typedef struct
{
	const char * cmd;
	const char * help;
	void(*fn)();
	int maskSet, maskClr;
}
tParamTbl;

static void main_help();

static const tParamTbl tbl[] =
{
	{ "-func", "functional testing", NULL, eTESTFLAG_FUNC, eTESTFLAG_MIPS | eTESTFLAG_ACCURACY | eTESTFLAG_EXCEPTION },
	{ "-mips", "test for performance", NULL, eTESTFLAG_MIPS, eTESTFLAG_FUNC | eTESTFLAG_ACCURACY | eTESTFLAG_EXCEPTION },
	{ "-accuracy", "testing for accuracy", NULL, eTESTFLAG_ACCURACY, eTESTFLAG_MIPS | eTESTFLAG_FUNC },
	{ "-exception", "testing for exceptions", NULL, eTESTFLAG_EXCEPTION, eTESTFLAG_MIPS | eTESTFLAG_FUNC },
	{ "-noabort", "do not stop after the first failure", NULL, eTESTFLAG_NOABORT, 0 },
	{ "-verbose", "print test progess info", NULL, eTESTFLAG_VERBOSE, 0 },
	{ "-vreport", "print verification test report", NULL, eTESTFLAG_VERIFY, 0 },
	{ "-full", "use full-size test vectors", NULL, eTESTFLAG_FULL, eTESTFLAG_BRIEF | eTESTFLAG_SANITY },
	{ "-sanity", "use sanity test vectors", NULL, eTESTFLAG_SANITY, eTESTFLAG_BRIEF | eTESTFLAG_FULL },
	{ "-brief", "use shorter test vectors", NULL, eTESTFLAG_BRIEF, eTESTFLAG_SANITY | eTESTFLAG_FULL },
	{ "-nowarmup", "invalidate cache before each cycle measurement", NULL, eTESTFLAG_NOWARMUP, 0 },
	{ "-help", "this help", main_help, eTESTFLAG_HELP, 0 },
	{ "-h", 0, main_help, eTESTFLAG_HELP, 0 }
};

static void main_help()
{
	int k;
	printf("General options\n");
	for (k = 0; k<(int)(sizeof(tbl) / sizeof(tbl[0])); k++)
	{
		if (tbl[k].help == NULL) continue;
		printf("%-8s  %s\n", tbl[k].cmd, tbl[k].help);
	}
	printf("Library specific package options:\n");
	printf("%-8s  %s\n", "-phase<n>", "limit test coverage to specific datatype (1...5)");
	myPackageContainer.printHelp();
}


/*-------------------------------------------------
process command line option, returns 0 on error
-------------------------------------------------*/
int PackageProcessOpt(const char* cmd)
{
	int bFound, group;
	Package* list;
	int k;
	// scan table first!
	for (k = 0; k<(int)(sizeof(tbl) / sizeof(tbl[0])); k++)
	{
		if (strcmp(cmd, tbl[k].cmd) == 0)
		{
			myPackageContainer.flags |= tbl[k].maskSet;
			myPackageContainer.flags &= ~tbl[k].maskClr;
			if (tbl[k].fn) tbl[k].fn();
			return 1;
		}
	}
	if (strlen(cmd) == 7 && memcmp("-phase", cmd, 6) == 0 && isdigit(cmd[6]))
	{
		myPackageContainer.phaseFlags |= 1 << (cmd[6] - '0');
		return 1;
	}
	list = Package::first;
	bFound = 0; group = 0;
	while (list != NULL)
	{
		if (strcmp(cmd, list->cmd) == 0)
		{
			list->execFlags |= 1;
			if (list->number == 0) group = list->group; // later, mark all members of this group for execution
			bFound = 1;
		}
		list = list->next;
	}
	if (bFound == 0) return 0;
	if (group != 0)
	{
		// mark all members of this group for execution
		list = Package::first;
		while (list != NULL)
		{
			if (list->group == group) list->execFlags |= 1;
			list = list->next;
		}
	}
	return 1;
}

eTestType PackageGetTestType()
{
	int isAccuracy = (myPackageContainer.flags & eTESTFLAG_ACCURACY) ? 1 : 0;
	int isException = (myPackageContainer.flags & eTESTFLAG_EXCEPTION) ? 1 : 0;
	int isMips = (myPackageContainer.flags & eTESTFLAG_MIPS) ? 1 : 0;
	int isFunc = (myPackageContainer.flags & eTESTFLAG_FUNC) ? 1 : 0;
	if (isAccuracy || isException) return eTEST_ACCURACY;
	if (isMips) return eTEST_MIPS;
	if (isFunc) return eTEST_FUNC;
	assert(0);
	return eTEST_UNKNOWN;
}

eTestFlags PackageGetTestFlags()
{
	return (eTestFlags)myPackageContainer.flags;
}


/*-------------------------------------------------
perform tests
function perf_info is used only for cycle
measurement
-------------------------------------------------*/
int  PackagePerformTest(int(*perf_info)(FILE * f, const char * fmt, ...))
{
	FILE *fout = NULL;
	Package* list, *groupMain;
	int bPrintedGroup;
	int  group = -1, num = -1, phase = -1, phaseFl = myPackageContainer.phaseFlags;
	if (phaseFl == 0) phaseFl = ~0;
	int isFull = (myPackageContainer.flags & eTESTFLAG_FULL) ? 1 : ( (myPackageContainer.flags & eTESTFLAG_BRIEF) ? 0 : 2 );
	int isVerbose = (0 != (myPackageContainer.flags & eTESTFLAG_VERBOSE));
	int isBreak = (0 == (myPackageContainer.flags & eTESTFLAG_NOABORT));
	int isAccuracy = (myPackageContainer.flags & eTESTFLAG_ACCURACY) ? 1 : 0;
	int isException = (myPackageContainer.flags & eTESTFLAG_EXCEPTION) ? 1 : 0;
	int isMips = (myPackageContainer.flags & eTESTFLAG_MIPS) ? 1 : 0;
	int isFunc = (myPackageContainer.flags & eTESTFLAG_FUNC) ? 1 : 0;

	if (isMips)
  {
    myPackageContainer.flags &= ~eTESTFLAG_SANITY; isFull &=~2;    
		fout = fopen("performance.txt", "wb");
		perf_info(fout, "Function Name\tCycles Measurements\t\n\tInvocation Parameters\tCycles\n");
	}
	else
	{
		fout = stdout;
		perf_info = fprintf;
	}

	// first, recheck execFlags. if no one is set, set all of them 
	{
		int bSet = 0;
		for (list = Package::first; list != NULL; list = list->next) bSet |= list->execFlags;
		if (bSet == 0)
		{
			for (list = Package::first; list != NULL; list = list->next) list->execFlags |= 1;
		}
	}

	for (bPrintedGroup = 0, groupMain = NULL, list = Package::first; list != NULL; list = list->next)
	{
		if (group != list->group)
		{
			if (list->number == 0)
			{
				groupMain = list;
				bPrintedGroup = 0;
			}
			phase = -1;
			num = -1;
			group = list->group;
		}
		if (list->execFlags == 0) continue;
		if (num != list->number) phase = -1;
		num = list->number;

		if (list->func == NULL && list->mips == NULL && list->acc == NULL) continue;
		if ((1 << list->phase) & phaseFl)
		{
			if (bPrintedGroup == 0 && groupMain != NULL && groupMain->group == list->group)
			{
				perf_info(fout, "\n%s\n", groupMain->help);
				bPrintedGroup = 1;
			}

			if (phase == -1) perf_info(fout, "\n%s\n", list->help);
			phase = list->phase;
			if (isFunc && list->func)
			{
				if (!list->func(isFull, isVerbose, isBreak) && isBreak) break;
			}
			if (isMips && list->mips)
			{
				list->mips(isFull, isVerbose, fout);
			}
			if ((isAccuracy || isException) && list->acc)
			{
				if (!list->acc(isFull, isVerbose, isBreak, isAccuracy, isException) && isBreak) break;
			}
		}
	}
	perf_info(fout, "================= test completed =================\n");
	if (isMips) fclose(fout);
	return 1;
}
