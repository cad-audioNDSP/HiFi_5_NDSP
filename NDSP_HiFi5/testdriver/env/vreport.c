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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "types.h"
#include "vreport.h"
#include "addr2name.h"
#include "testcase.h"

/* list of data files */
typedef struct tagDataFileList
{
    char*                   dataFileName;
    struct tagDataFileList* next;
}
tDataFileList;
/* array of FUT */
typedef struct 
{
    tReportFUT *fut;
    int N;
}
tReportFUTArray;

typedef struct tagReportList
{
    tReportFUTArray       futarr;
    char*                 extraParams;
    tDataFileList*        list;
    int32_t               numCalls[MAXTESTCASE];    /* number of calls for up to 16 possible test case types */
    size_t                dataSize[MAXTESTCASE]; /* size of data for each kind of test case */
    evrResult             testResult;   /* final test results */
    evrResult             exceptResult;
    evrResult             memResult;
    struct tagReportList* next;            /* pointer to the next element */
}
tReportList;

typedef struct
{
    int turnOn;
    const char* exeName;
          char* mapfile;
    tReportList *lastTest;    /* last running test or NULL */

    tReportList *list;
}
tReport;

static tReport vReport;


/* datafile list management */
static tDataFileList* DataFileListAdd(tDataFileList* list, const char* dataFileName)
{
    tDataFileList* last=NULL;
    tDataFileList* ptr;
    /* find and insert non-duplicating name */
    ptr=list;
    while(ptr)
    {
        last=ptr;
        if(strcmp(ptr->dataFileName,dataFileName)==0) return list;
        ptr=ptr->next;
    }
    /* not found so add it */
    ptr=(tDataFileList*)calloc(1,sizeof(tDataFileList));
    ptr->dataFileName=(char*)malloc(strlen(dataFileName)+1);
    strcpy(ptr->dataFileName,dataFileName);
    ptr->next=NULL;
    if (last) last->next=ptr;
    if (list==NULL) list = ptr;
    return list;
}

static tDataFileList* DataFileListRemove(tDataFileList* list)
{
    while(list)
    {
        tDataFileList* ptr=list;
        list=list->next;
        free(ptr->dataFileName);
        free(ptr);
    }
    return list;
}

/* FUT array management */
static void ReportFUTArrayCreate(tReportFUTArray* arr, const tReportFUT fut[], int Nfun)
{
    int n;
    arr->fut=(tReportFUT*)calloc(Nfun,sizeof(tReportFUT));
    arr->N=Nfun;
    for (n=0; n<Nfun; n++)
    {
        arr->fut[n]=fut[n];
    }
}

static void ReportFUTArrayRemove(tReportFUTArray* arr)
{
    arr->N=0;
    free(arr->fut);
}

/* compare arrays, returns 0 if they are the same */
int ReportFUTArrayCmp(const tReportFUTArray* arr, const tReportFUT fut[], int Nfun)
{
    int n,eq=1;
    if (arr->N!=Nfun) return -1;
    eq=1;
    for (n=0; n<Nfun; n++) 
    {
        eq &= (arr->fut[n]==fut[n]);
    }
    return eq ? 0:1;
}

/* initialize verification reporting */
void vReportInit(int turnOn, const char* exeName)
{
    vReport.turnOn=turnOn;
    vReport.exeName=exeName;
    vReport.mapfile=(char*)malloc(strlen(exeName)+5); sprintf(vReport.mapfile,"%s.map",vReport.exeName);
    vReport.lastTest=NULL;
    vReport.list=NULL;
    /* verify presence and validness of existing map-file */
    if (turnOn)
    {
        static uintptr_t knownFunctons[]={
            (uintptr_t)&vReportInit,
            (uintptr_t)&vReportAdd,
            (uintptr_t)&vReportSetResult,
            (uintptr_t)&vReportPrint,
            (uintptr_t)&vReportClose,
            (uintptr_t)&vReportFinish
        };
        int n,ok=1;
        for (n=0; n<(int)(sizeof(knownFunctons)/sizeof(knownFunctons[0])); n++)
        {
            char* funName;
            funName=addr2name(vReport.mapfile,knownFunctons[n]);
            ok &= funName!=NULL;
            free(funName);
        }
        if (!ok)
        {
            printf("Valid MAP-file is required for verification report! Function names might be wrong or replaced with their addresses!\n");
        }
    }

}

/* the same as strcmp but may accept NULL pointers */
static int __strcmp(const char* str1,const char* str2)
{
    if (str1==NULL) return str2==NULL ? 0:1;
    if (str2==NULL) return str1==NULL ? 0:1;
    return strcmp(str1,str2);
}

/* management of report list */
static tReportList* ReportFindInsert(tReportList* list,const tReportFUT *fut,int Nfun, const char* extraParams)
{
    tReportList *last=NULL;
    tReportList *ptr=list;
    while(ptr)
    {
        last=ptr;
        if (ReportFUTArrayCmp(&ptr->futarr,fut,Nfun)==0 && __strcmp(ptr->extraParams,extraParams)==0) return ptr;
        ptr=ptr->next;
    }
    /* not found - allocate empty entry */
    ptr=(tReportList*)calloc(1,sizeof(tReportList));
    memset(ptr,0,sizeof(tReportList));
    ReportFUTArrayCreate(&ptr->futarr,fut,Nfun);
    if (extraParams) 
    {
        strcpy(ptr->extraParams=(char*)malloc(strlen(extraParams)+1),extraParams);
    }
    if (last) last->next=ptr;
    return ptr;
}


static tReportList* ReportListRemove(tReportList* list)
{
    while(list)
    {
        tReportList* ptr=list;
        list=list->next;
        if (ptr->extraParams) free(ptr->extraParams);
        ReportFUTArrayRemove(&ptr->futarr);
        DataFileListRemove(ptr->list);
        free(ptr);
    }
    return list;
}

/* add to statistics 
   Input:
   fun_ptr_t[Nfun]    pointer to the functions under the test (FUT). 
   Nfun             number of FUT (several FUT might be tested together by one test case)
   extraParams        extra parameters of function (NULL if there is no specialized variants of FUT)
   dataFile         filename with data
   caseType            test case type 
   dataSize         total size of input/output data for FUT
*/
void vReportAdd( const tReportFUT fun_ptr_t[], int Nfun, const char* extraParams, const char* dataFile, int caseType, size_t dataSize)
{
    tReportList* ptr;
    if (vReport.turnOn==0) return;

    /* find and insert stat for specific function */
    ptr=ReportFindInsert(vReport.list,fun_ptr_t,Nfun,extraParams);
    if (vReport.list==NULL) vReport.list=ptr;
    ptr->list=DataFileListAdd(ptr->list,dataFile);
    if (caseType>=0 && caseType<(int)(sizeof(ptr->numCalls)/sizeof(ptr->numCalls[0])))
    {
        ptr->numCalls[caseType]++;
        ptr->dataSize[caseType]+=dataSize;
    }
    vReport.lastTest=ptr;
}

/* set results for the lasting run test 
   testResult   result of data testing of FUT
   exceptResult result of exception handling testing of FUT
   memResult    result of testing the memory after the call of FUT

   set evr_NOTTESTED if specific test result is not known at a moment
*/
void vReportSetResult (evrResult testResult, evrResult exceptResult, evrResult memResult)
{
    tReportList* ptr;
    ptr=vReport.lastTest;
    if (ptr==NULL) return;
    /* ignore OK if already failed and so on */
    if ((int)ptr->testResult  <=(int)testResult)   ptr->testResult  =testResult;
    if ((int)ptr->exceptResult<=(int)exceptResult) ptr->exceptResult=exceptResult;
    if ((int)ptr->memResult   <=(int)memResult)    ptr->memResult   =memResult;
}

/*
    just stop reporting by vReportSetResult() until next calls of vReportAdd()
    call this function upon the end of data file
*/
void vReportFinish()
{
    vReport.lastTest=NULL;
}


/* print report */
void vReportPrint()
{
    tReportList* ptr;
    if (vReport.turnOn==0) return;
    printf("==================================================\n");
    printf("Verification report\n");
    printf("==================================================\n");
    for(ptr=vReport.list;ptr;ptr=ptr->next)
    {
        int n;
        tDataFileList*        dataFileList;
        printf("Function(s) : ");
        for (n=0; n<ptr->futarr.N; n++)
        {
            char* funName;
            funName=addr2name(vReport.mapfile,(uintptr_t)ptr->futarr.fut[n]);
            if (funName)
            {
                printf("%s%s",funName,n==ptr->futarr.N-1 ? "":",");
                free(funName);
            }
            else
            {
                printf("0x%08x%s",(uintptr_t)ptr->futarr.fut[n],n==ptr->futarr.N-1 ? "":", ");
            }
        }
        printf("\n");
        if (ptr->extraParams && ptr->extraParams[0]!=0) printf("Parameters  : %s\n", ptr->extraParams );
        dataFileList=ptr->list;
        printf("Data file(s):");
        while(dataFileList)
        {
            printf("%s%s",dataFileList->dataFileName , dataFileList->next?", ":"");
            dataFileList=dataFileList->next;
        }
        printf("\n");
        for (n=0; n<(int)(sizeof(ptr->numCalls)/sizeof(ptr->numCalls[0])); n++)
        {
            if (ptr->numCalls[n]) printf("Test case type %d (%s): %d (%lu bytes)\n",n,testCaseTypeStr[n],ptr->numCalls[n],(unsigned long)ptr->dataSize[n]);
        }
        printf("Functional tests            : %s\n",ptr->testResult==evr_NOTTESTED   ? "NOT TESTED" :(ptr->testResult==evr_OK   ?"OK":"FAIL"));
        printf("FPU Exception handling tests: %s\n",ptr->exceptResult==evr_NOTTESTED ? "NOT TESTED" :(ptr->exceptResult==evr_OK ?"OK":"FAIL"));
        printf("Memory tests                : %s\n",ptr->memResult==evr_NOTTESTED    ? "NOT TESTED" :(ptr->memResult==evr_OK    ?"OK":"FAIL"));
    }
    
    {
        FILE *f;
        int n;
        f=fopen("verification_summary.csv","wt");
        if (f)
        {
            tDataFileList*        dataFileList;
            size_t totalDataSize;
            fprintf(f,"Function(s) name;");
            fprintf(f,"Data file(s);");
            fprintf(f,"Parameters;");
            fprintf(f,"Functional tests;");
            fprintf(f,"FPU Exception handling tests;");
            fprintf(f,"Memory tests;");
            for (n=0; n<(int)(sizeof(ptr->numCalls)/sizeof(ptr->numCalls[0])); n++)
            {
                fprintf(f,"%s;",testCaseTypeStr[n]);
            }
            fprintf(f,"Data size, bytes\n");
            for(ptr=vReport.list;ptr;ptr=ptr->next)
            {
                for (n=0; n<ptr->futarr.N; n++)
                {
                    char* funName;
                    funName=addr2name(vReport.mapfile,(uintptr_t)ptr->futarr.fut[n]);
                    if (funName)
                    {
                        fprintf(f,"%s%s",funName,n==ptr->futarr.N-1 ? "":",");
                        free(funName);
                    }
                    else
                    {
                        fprintf(f,"0x%08x%s",(uintptr_t)ptr->futarr.fut[n],n==ptr->futarr.N-1 ? "":", ");
                    }
                }
                fprintf(f,";");
                dataFileList=ptr->list;
                while(dataFileList)
                {
                    fprintf(f,"%s%s",dataFileList->dataFileName , dataFileList->next?", ":"");
                    dataFileList=dataFileList->next;
                }
                fprintf(f,";");
                fprintf(f,"%s;", (ptr->extraParams && ptr->extraParams[0]!=0) ? ptr->extraParams : "" );
                fprintf(f,"%c;",ptr->testResult==evr_NOTTESTED   ? ' ' :(ptr->testResult==evr_OK   ?'v':'x'));
                fprintf(f,"%c;",ptr->exceptResult==evr_NOTTESTED ? ' ' :(ptr->exceptResult==evr_OK ?'v':'x'));
                fprintf(f,"%c;",ptr->memResult==evr_NOTTESTED    ? ' ' :(ptr->memResult==evr_OK    ?'v':'x'));
                for (totalDataSize=0,n=0; n<(int)(sizeof(ptr->numCalls)/sizeof(ptr->numCalls[0])); n++)
                {
                    fprintf(f,"%c;",ptr->numCalls[n]?'v':' ');
                    totalDataSize+=ptr->dataSize[n];
                }
                fprintf(f,"%lu\n",(unsigned long)totalDataSize);
            }
            fclose(f);
        }
        else
        {
            printf("File verification_summary.csv can not be open for writing\n");
        }
    }
}

/* deallocate resources */
void vReportClose()
{
    ReportListRemove(vReport.list);
    free(vReport.mapfile);
}
