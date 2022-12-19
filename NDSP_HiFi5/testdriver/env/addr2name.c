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
 * convert address to name using map-file
 */
#include "types.h"
#include "addr2name.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAXLEN 200

static uintptr_t findLoadOffset(const char *mapfile)
#ifdef COMPILER_MSVC
{
    char str[MAXLEN];
    char fun[MAXLEN];
    int sect,haddr;
    FILE * f;
    f=fopen(mapfile,"rt");
    if (f==NULL) { return 0; }
    while(fgets(str,MAXLEN,f))
    {
        str[MAXLEN-1]=0;
        if (sscanf(str,"%d:%x %s",&sect,&haddr,fun)==3)
        {
            if (sect==1 && strcmp(fun,"_addr2name")==0)
            {
                fclose(f);
                return ((uintptr_t)&addr2name)-(uintptr_t)(haddr);
            }
        }
    }
    fclose(f);
    return 0;
}
#else
{
    char str[MAXLEN];
    char fun[MAXLEN];
    int haddr;
    FILE * f;
    f=fopen(mapfile,"rt");
    if (f==NULL) { fclose(f); return 0; }
    while(fgets(str,MAXLEN,f))
    {
        str[MAXLEN-1]=0;
        if (sscanf(str," 0x%x %s",&haddr,fun)==2)
        {
            if (strcmp(fun,"addr2name")==0)
            {
                fclose(f);
                return ((uintptr_t)&addr2name)-(uintptr_t)(haddr);
            }
        }
    }
    fclose(f);
    return 0;
}
#endif

/*
   returns the name of function or NULL if not found
   Input:
   addr    address
   mapfile name of map-file
   Returned value - pointer to ASCII string (have to be 
   freed after use). NULL if not found
*/
char* addr2name(const char *mapfile, uintptr_t addr)
#ifdef COMPILER_MSVC
{
    //   version for Visual Studio
    uintptr_t loadOff;
    char str[MAXLEN];
    char *fun;
    int sect,haddr;
    FILE * f;
    loadOff=findLoadOffset(mapfile);
    addr-=loadOff;
    f=fopen(mapfile,"rt");
    if (f==NULL) { return NULL; }
    fun=(char*)malloc(MAXLEN);
    while(fgets(str,MAXLEN,f))
    {
        str[MAXLEN-1]=0;
        if (sscanf(str,"%d:%x %s",&sect,&haddr,fun)==3)
        {
            if (sect==1 && (uintptr_t)haddr==addr) 
            {
                fclose(f);
                return fun;
            }
        }
    }
    fclose(f);
    free(fun);
    return NULL;
}
#else
{
    //   version for gcc linker
    uintptr_t loadOff;
    char str[MAXLEN];
    char *fun;
    int haddr;
    FILE * f;
    loadOff=findLoadOffset(mapfile);
    addr-=loadOff;
    f=fopen(mapfile,"rt");
    if (f==NULL) { fclose(f); return NULL; }
    fun=(char*)malloc(MAXLEN);
    while(fgets(str,MAXLEN,f))
    {
        str[MAXLEN-1]=0;
        if (sscanf(str," 0x%x %s",&haddr,fun)==2)
        {
            if ((uintptr_t)haddr==addr) 
            {
                fclose(f);
                return fun;
            }
        }
    }
    fclose(f);
    free(fun);
    return NULL;
}
#endif


