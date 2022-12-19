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
 * Aligning wrapper over standard malloc/free.
 */

#include <stdlib.h>
#include <stdio.h>
/* Portable data types. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
/* Aligning wrapper over standard malloc/free. */
#include "malloc16.h"

#define MAGIC   0xf2a3
#define MAGIC_REGION   0xf3b2
static void simpleAlloc_init(void* region, size_t sz);
static void* simpleAlloc_alloc(size_t sz);
static void simpleAlloc_free(void* ptr);

/* Allocation header that is stored immediately before an aligned memory area.
 * We use 8-bit fields to avoid any padding and/or alignment restrictions for the
 * header structure. */
typedef struct alloc_hdr_tag {
  uint8_t magic0, magic1;
  uint8_t padSz0, padSz1;
} alloc_hdr_t;

typedef void* (*fnmalloc)(size_t);
static void * _alloc(size_t sz, size_t align,fnmalloc pfn,uint16_t magic );

/*
 Allocate and align a memory area.
 Input:
   sz     Area size, in bytes
   align  A power of 2 number specifying the alignment boundary, in bytes.
          Set it to zero to select the default alignment, that is the SIMD
          vector width for the target platform.
 Returns:
   Pointer to the aligned memory area if succeeded, or null pointer
   if failed to allocate the memory. 
*/
void * mallocAlign(size_t sz, size_t align)
{
    return _alloc(sz,align,malloc,MAGIC);
}
void * mallocAlignRegion(size_t sz, size_t align)
{
    return _alloc(sz,align,simpleAlloc_alloc,MAGIC_REGION);
}


static void * _alloc(size_t sz, size_t align,fnmalloc pfnmalloc,uint16_t magic )
{
    size_t allocSz, padSz, hdrSz = sizeof(alloc_hdr_t);
    alloc_hdr_t * hdr;
    void * ptr;
    if (!align) {
        align = GET_ISA_OPT(INT16_SIMD_WIDTH)*sizeof(int16_t);
    }
    ASSERT(0==(align&(align-1)));
    ASSERT(0<align && align<=65536);
    allocSz = align-1 + hdrSz + sz;
    if (!(ptr = pfnmalloc(allocSz))) return (void*)0;
    padSz = -(int)((uintptr_t)ptr+hdrSz) & (align-1);
    hdr = (alloc_hdr_t*)((uintptr_t)ptr+padSz);
    hdr->magic0 = (uint8_t)(magic>>8);
    hdr->magic1 = (uint8_t)magic;
    hdr->padSz0 = (uint8_t)(padSz>>8);
    hdr->padSz1 = (uint8_t)padSz;
    return (void*)(hdr+1);
} /* _align() */

/*
 Free an aligned memory area.
 Input:
   ptr  Pointer to a memory area allocated with mallocAlign()
 Output:
   None
*/
void freeAlign(void * ptr)
{
    uint16_t magic, padSz;
    alloc_hdr_t * hdr;
    hdr = (alloc_hdr_t*)ptr-1;
    magic = ((uint16_t)hdr->magic0<<8) | (uint16_t)hdr->magic1;
    padSz = ((uint16_t)hdr->padSz0<<8) | (uint16_t)hdr->padSz1;
    switch(magic)
    {
    case MAGIC:         free((void*)((uintptr_t)hdr-padSz)); break;
    case MAGIC_REGION:  simpleAlloc_free((void*)((uintptr_t)hdr-padSz)); break;
    default:            ASSERT(!"freeAlign(): bad argument");
    }
} /* freeAlign() */

/*-------------------------------------------------------------------------------
additional functions allowing to allocate memory not from heap but from 
specified region
For allocation, user should map the region by mallocAlignInitRegion() and 
call mallocAlignRegion() instead of mallocAlign
-------------------------------------------------------------------------------*/
void mallocAlignInitRegion(void* region, size_t sz)
{
    simpleAlloc_init(region,sz);
}

/*-----------------------------------------------------------------
very simple memory allocator which uses specified region for 
allocation - use with care!
-----------------------------------------------------------------*/
#define MAGIC_FREE 0x376257ab
#define MAGIC_BUSY 0xfa35120d
#define MAGIC_END  0xfefefefe
/* 2way list of allocated blocks:
   1st element has prev==NULL
   last element has next==NULL and magic==MAGIC_END
*/
typedef struct taglist
{
    uint32_t magic;       /* MAGIC_FREE/MAGIC_BUSY */
    struct taglist *prev;
    struct taglist *next;
}
tlist;

typedef struct
{
    void* regionBegin;
    void* regionEnd;
}
tSimpleAlloc;

static tSimpleAlloc simpleAlloc={NULL,NULL};

typedef struct 
{
    int nAllocated;
    int nFreed;
}
tstat;
/* simple validation function */
static tstat simpleAlloc_validate()
{
    tstat stat={0,0};
    tlist* list=(tlist *)(simpleAlloc.regionBegin);
    int bFreed=0;
    while(list->magic!=MAGIC_END)
    {
        //if (list->magic==MAGIC_FREE || list->magic==MAGIC_BUSY);
        NASSERT(list->next!=NULL);
        NASSERT(list->next->prev==list);
        if (list->prev) NASSERT(list->prev->next==list);
        stat.nAllocated+=(list->magic==MAGIC_BUSY)?1:0;
        stat.nFreed     +=(list->magic==MAGIC_FREE)?1:0;
        if (bFreed)
        {   // allocator always binds successive freed blocks !
            NASSERT(list->magic==MAGIC_BUSY);
        }
        NASSERT(list->magic==MAGIC_FREE || list->magic==MAGIC_BUSY);
        bFreed=list->magic==MAGIC_FREE;
        list=list->next;
    }
    return stat;
}

/* initialization of region */
static void simpleAlloc_init(void* region, size_t sz)
{
    tlist *last,*first;
    simpleAlloc.regionBegin=region;
    simpleAlloc.regionEnd  =(void*)((uintptr_t)region+sz);
    last=((tlist *)(simpleAlloc.regionEnd))-1;
    first=(tlist *)(simpleAlloc.regionBegin);
    last->magic=MAGIC_END;
    last->next =NULL;
    last->prev=first;
    first->magic=MAGIC_FREE;
    first->next =last;
    first->prev=NULL;
    simpleAlloc_validate();
}

/* combine successive freed block */
static void simpleAlloc_combineFreed(tlist* list)
{
    if (list->magic==MAGIC_BUSY) return;
    NASSERT(list->magic==MAGIC_FREE);
    /* combine with next freed block */
    while (list->next && list->next->magic==MAGIC_FREE)
    {
        list=list->next;
        list->prev->next=list->next;
        list->next->prev=list->prev;
    }
    /* combine with previous freed block */
    while (list->prev && list->prev->magic==MAGIC_FREE)
    {
        list->prev->next=list->next;
        list->next->prev=list->prev;
        list=list->prev;
    }
}

/* allocate from region */
static void* simpleAlloc_alloc(size_t sz)
{
    NASSERT(simpleAlloc.regionBegin);
    sz=(sz+sizeof(tlist)+7)&~7;
    /* find first candidate */
    tlist* list=(tlist *)(simpleAlloc.regionBegin);
    do
    {
        tlist *newEntry,*next;
        next=list->next;
        if (next==NULL) break;
        NASSERT(list->magic==MAGIC_FREE || list->magic==MAGIC_BUSY);
        size_t size=(uintptr_t)next-(uintptr_t)list;
        size-=sizeof(tlist);
        if (size>=sz && list->magic==MAGIC_FREE)
        {
            /* split current block into 2 segments */
            newEntry=(tlist *)(((uintptr_t)list)+sz);
            newEntry->magic=MAGIC_FREE;
            newEntry->next=list->next;
            newEntry->prev=list;
            list->next=newEntry;
            list->magic=MAGIC_BUSY;
            next->prev=newEntry;
            simpleAlloc_combineFreed(newEntry);
            simpleAlloc_validate();
            return (void*)(list+1);
        }
        list=list->next;
    }
    while(list!=NULL);
    printf("memory block %d can not be allocated in memory region %08x \n",sz,(uintptr_t)simpleAlloc.regionBegin);
    return NULL;
}

/* free from region */
static void simpleAlloc_free(void* ptr)
{
    tlist * list;
    NASSERT(simpleAlloc.regionBegin);
    NASSERT(ptr);
    list=((tlist *)ptr)-1;
    if (list->magic!=MAGIC_BUSY)
    {
        NASSERT(list->magic==MAGIC_BUSY);
    }
    NASSERT(list->magic==MAGIC_BUSY);
    list->magic=MAGIC_FREE;
    simpleAlloc_combineFreed(list);
    simpleAlloc_validate();
}

#if 0
/* test utility */
#include <stdlib.h>
#include <string.h>
uint64_t tmp[20000];
#define NMEM 32
uint8_t *pMem[NMEM];
void test_alloc()
{
    int n,nAllocated=0;
    tstat stat;
    for (n=0; n<NMEM; n++) pMem[n]=NULL;
    simpleAlloc_init((void*)tmp,sizeof(tmp));
    for (n=0; n<10000; n++)
    {
        int bAlloc=rand()&1;
        int idx   =rand()&(NMEM-1);
        int nBytes=rand()&127;
        int rndFill=rand()&255;
        if (bAlloc) 
        {
            if (pMem[idx]!=NULL) 
            {
                simpleAlloc_free(pMem[idx]);
                pMem[idx]=NULL;
                nAllocated--;
                stat=simpleAlloc_validate();
                nAllocated=nAllocated;
            }
            pMem[idx]=simpleAlloc_alloc(nBytes);
            memset(pMem[idx],rndFill,nBytes);
            nAllocated++;
        }
        else
        {
            if (pMem[idx]!=NULL) 
            {
                simpleAlloc_free(pMem[idx]);
                pMem[idx]=NULL;
                nAllocated--;
            }
        }
        stat=simpleAlloc_validate();
        printf("%d\n",n);
    }
    for (n=0; n<NMEM; n++)
    {
        if (pMem[n]!=NULL) 
        {
            simpleAlloc_free(pMem[n]);
            pMem[n]=NULL;
        }
    }
}
#endif
