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
 * floating point extra statistics
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "fpstat.h"
#include "types.h"
#include "float16.h"

typedef struct tagStatList
{
    char*     patName;  // pattern name
    double    minErr;   // min peak err (negative diff)
    double    maxErr;   // max peak err (pos diff)
    double    avgErr;   // average absolute error
    int       Npts;     // number of points
    struct tagStatList* next;            /* pointer to the next element */
}
tStat;

/* stat list management */
static tStat* StatAdd(tStat* list, const char* patName)
{
    tStat* last=NULL;
    tStat* ptr;
    /* find and insert non-duplicating name */
    ptr=list;
    while(ptr)
    {
        last=ptr;
        if(strcmp(ptr->patName,patName)==0) return ptr;
        ptr=ptr->next;
    }
    /* not found so add it */
    ptr=(tStat*)calloc(1,sizeof(tStat));
    memset(ptr,0,sizeof(tStat));
    ptr->patName=(char*)malloc(strlen(patName)+1);
    strcpy(ptr->patName,patName);
    ptr->next=NULL;
    if (last) last->next=ptr;
    if (list==NULL) list = ptr;
    return list;
}

static tStat* StatRemove(tStat* list)
{
    while(list)
    {
        tStat* ptr=list;
        list=list->next;
        free(ptr->patName);
        free(ptr);
    }
    return list;
}


// open statistics
void FPStatCreate(tFPStat* pFPStat)
{
    pFPStat->list=NULL;
}

/*-----------------------------------------
add data to the statistics
Input:
patName    - name of pattern
refData[N] - reference data
data[N]    - data
N          - number of data points
-----------------------------------------*/
void FPStatAdd(tFPStat* pFPStat, const char* patName, const float64_t* refData, const float32_t* data, int N)
{
    static const union ufloat32uint32 realmaxf ={0x7f7fffff};
    static const union ufloat32uint32 plusInff ={0x7f800000};
    int n,M;
    tStat* list;
    list=StatAdd(pFPStat->list,patName);
    if (pFPStat->list==NULL) pFPStat->list=list;
    M=list->Npts;
    for (n=0; n<N; n++)
    {
        union ufloat32uint32 tmp;
        float64_t rdata;
        int binf, bnan;
        float64_t ulp;
        tmp.f = data[n]; 
        tmp.u ^= 1;
        ulp = fabsf(data[n]-tmp.f);
        rdata=refData[n];
        ulp = ((float64_t)data[n]-rdata)/ulp;
        binf = (fabs(rdata)>(float64_t)realmaxf.f && fabsf(data[n])==plusInff.f); // both Inf
        bnan = !(refData[n]==refData[n]) && !(data[n]==data[n]);     // both NaNs
        if (bnan) continue; // omit NaNs from statistics
        if (binf) ulp=0.f;  // assume 0 on Inf
        if (ulp<list->minErr) list->minErr=ulp;
        if (ulp>list->maxErr) list->maxErr=ulp;
        list->avgErr = (list->avgErr*M+fabs(ulp))/(M+1);  // update average
        M++;
    }
    list->Npts=M;
}
void FPStatAdd_fp16(tFPStat* pFPStat, const char* patName, const float64_t* refData, const float16_t* data, int N)
{
    static const double realmax_fp16 =65504.;
    int n,M;
    tStat* list;
    list=StatAdd(pFPStat->list,patName);
    if (pFPStat->list==NULL) pFPStat->list=list;
    M=list->Npts;
    for (n=0; n<N; n++)
    {
        union ufloat16uint16 tmp;
        float64_t rdata;
        int binf, bnan;
        float64_t ulp,datan=conv_f16_to_f64(data[n]);
        tmp.f = data[n]; 
        tmp.u ^= 1;
        ulp = fabs(datan-conv_f16_to_f64(tmp.f));
        rdata=refData[n];
        ulp = (datan-rdata)/ulp;
        binf = (fabs(rdata)>realmax_fp16 && isinf_f16(data[n])); // both Inf
        bnan = !(refData[n]==refData[n]) && !eq_f16(data[n],data[n]);     // both NaNs
        if (bnan) continue; // omit NaNs from statistics
        if (binf) ulp=0.f;  // assume 0 on Inf
        if (ulp<list->minErr) list->minErr=ulp;
        if (ulp>list->maxErr) list->maxErr=ulp;
        list->avgErr = (list->avgErr*M+fabs(ulp))/(M+1);  // update average
        M++;
    }
    list->Npts=M;
}

// print statistics to the file
void FPStatPrint(tFPStat* pFPStat, FILE* f)
{
    tStat* list=pFPStat->list;
    while (list)
    {
        fprintf(f,"\t%s\tmin ULP err\t%5.2lf\tmax ULP err\t%5.2lf\tavg abs ULP err\t%5.2lf\n",
            list->patName,
            list->minErr,
            list->maxErr,
            list->avgErr);

        list=list->next;
    }
}

// close statisctics
void FPStatClose (tFPStat* pFPStat)
{
    pFPStat->list=StatRemove(pFPStat->list);
}
