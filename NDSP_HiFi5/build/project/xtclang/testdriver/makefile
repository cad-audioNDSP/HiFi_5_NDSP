# makefile for building IntegrIT(R) DSP Library
# Version 0.34 tuned for Cadence HiFi cores
# ---------------------------------------------------
# Target platform: gcc/xcc
#
#  make invocation:
#  make clean   - clean up all temporary files
#  make depends - refresh dependencies (makedeps.mk)
#  make all     - make depends and release
#-----------------------------------------------------
#  build under xcc
#  make sure that all neccesary Xtensa system variables 
#  and pathes configured properly. makefile uses 
#  XTENSA_CORE system variable to detect the name of core
#  so, if makefile is running from the command-line shell
#  executed from proper processor config of Xtensa Exmplorer 
#  it does not require setting of any additional parameters
#
#  System variables XTENSA_CORE,XTENSA_SYSTEM should.
#  be properly configured for correct procedure.
#  if XTENSA_SYSTEM is undefined, g++ build procedure
#  will use cstub files and headers from local
#  directory, not from system path
#
#  Other variables:
#  LANG=C++    build under g++, otherwise xcc used
#  LANG=LLVM   build under llvm Xtensa compiler
#  SAVETEMPS=1 save temporary files
#  MEM_MODEL=1,2 enables memory modeling
#-----------------------------------------------------  
#  build under gcc require unsetting of variable XTENSA_SYSTEM
#  (either by command line unset XTENSA_SYSTEM or via make 
#  command line option -e XTENSA_SYSTEM=)
#  and selecting proper core by XTENSA_CORE variable. 
#  By default, it use C language, C++ will be selected by 
#  using LANG=C++ system variable.
# 
#  examples:
#  1. build with with cstubs :
#  make -e XTENSA_SYSTEM= XTENSA_CORE=XRC_B64EP_PM_20apr
#  2. build under g++ with cstubs:
#  make -e XTENSA_SYSTEM= XTENSA_CORE=XRC_B64EP_PM_20apr LANG=C++
#-----------------------------------------------------

SRCDIRS =	\
  ../../../../testdriver \
  ../../../../testdriver/env \
  ../../../../testdriver/include \
  ../../../../testdriver/pkg \
  ../../../../testdriver/pkg/audio \
  ../../../../testdriver/pkg/audio/common \
  ../../../../testdriver/pkg/audio/phase1 \
  ../../../../testdriver/pkg/audio/phase2 \
  ../../../../testdriver/pkg/complex \
  ../../../../testdriver/pkg/complex/common \
  ../../../../testdriver/pkg/complex/complexs \
  ../../../../testdriver/pkg/complex/complexs/common \
  ../../../../testdriver/pkg/complex/complexs/phase2 \
  ../../../../testdriver/pkg/complex/complexv \
  ../../../../testdriver/pkg/complex/complexv/common \
  ../../../../testdriver/pkg/complex/complexv/phase2 \
  ../../../../testdriver/pkg/ef \
  ../../../../testdriver/pkg/ef/common \
  ../../../../testdriver/pkg/ef/phase1 \
  ../../../../testdriver/pkg/fft \
  ../../../../testdriver/pkg/fft/cfft \
  ../../../../testdriver/pkg/fft/cfft/common \
  ../../../../testdriver/pkg/fft/cfft/phase1 \
  ../../../../testdriver/pkg/fft/cfftie \
  ../../../../testdriver/pkg/fft/cfftie/common \
  ../../../../testdriver/pkg/fft/cfftie/phase1 \
  ../../../../testdriver/pkg/fft/cfftie/phase2 \
  ../../../../testdriver/pkg/fft/cnfft \
  ../../../../testdriver/pkg/fft/cnfft/common \
  ../../../../testdriver/pkg/fft/cnfft/phase1 \
  ../../../../testdriver/pkg/fft/common \
  ../../../../testdriver/pkg/fft/dct \
  ../../../../testdriver/pkg/fft/dct/common \
  ../../../../testdriver/pkg/fft/dct/phase1 \
  ../../../../testdriver/pkg/fft/dct/phase2 \
  ../../../../testdriver/pkg/fft/rfft \
  ../../../../testdriver/pkg/fft/rfft/common \
  ../../../../testdriver/pkg/fft/rfft/phase1 \
  ../../../../testdriver/pkg/fft/rfftie \
  ../../../../testdriver/pkg/fft/rfftie/common \
  ../../../../testdriver/pkg/fft/rfftie/phase1 \
  ../../../../testdriver/pkg/fft/rfftie/phase2 \
  ../../../../testdriver/pkg/fft/rnfft \
  ../../../../testdriver/pkg/fft/rnfft/common \
  ../../../../testdriver/pkg/fft/rnfft/phase1 \
  ../../../../testdriver/pkg/fft/spectrum \
  ../../../../testdriver/pkg/fft/spectrum/common \
  ../../../../testdriver/pkg/fft/spectrum/phase1 \
  ../../../../testdriver/pkg/fft/spectrum/phase2 \
  ../../../../testdriver/pkg/fir \
  ../../../../testdriver/pkg/fir/common \
  ../../../../testdriver/pkg/fir/conv2d \
  ../../../../testdriver/pkg/fir/conv2d/common \
  ../../../../testdriver/pkg/fir/conv2d/phase1 \
  ../../../../testdriver/pkg/fir/conv2d/phase2 \
  ../../../../testdriver/pkg/fir/conv2d/phase3 \
  ../../../../testdriver/pkg/fir/firblk \
  ../../../../testdriver/pkg/fir/firblk/common \
  ../../../../testdriver/pkg/fir/firblk/phase1 \
  ../../../../testdriver/pkg/fir/firblk/phase2 \
  ../../../../testdriver/pkg/fir/firdec \
  ../../../../testdriver/pkg/fir/firdec/common \
  ../../../../testdriver/pkg/fir/firdec/phase1 \
  ../../../../testdriver/pkg/fir/firdec/phase2 \
  ../../../../testdriver/pkg/fir/firint \
  ../../../../testdriver/pkg/fir/firint/common \
  ../../../../testdriver/pkg/fir/firint/phase1 \
  ../../../../testdriver/pkg/fir/firint/phase2 \
  ../../../../testdriver/pkg/fir/firother \
  ../../../../testdriver/pkg/fir/firother/common \
  ../../../../testdriver/pkg/fir/firother/phase1 \
  ../../../../testdriver/pkg/fir/firother/phase2 \
  ../../../../testdriver/pkg/fit \
  ../../../../testdriver/pkg/fit/common \
  ../../../../testdriver/pkg/fit/pfit \
  ../../../../testdriver/pkg/fit/pfit/common \
  ../../../../testdriver/pkg/fit/pfit/phase1 \
  ../../../../testdriver/pkg/fit/pfit/phase2 \
  ../../../../testdriver/pkg/iir \
  ../../../../testdriver/pkg/iir/common \
  ../../../../testdriver/pkg/iir/iirbq \
  ../../../../testdriver/pkg/iir/iirbq/common \
  ../../../../testdriver/pkg/iir/iirbq/phase1 \
  ../../../../testdriver/pkg/iir/iirbq/phase2 \
  ../../../../testdriver/pkg/iir/iirkal \
  ../../../../testdriver/pkg/iir/iirkal/common \
  ../../../../testdriver/pkg/iir/iirkal/phase1 \
  ../../../../testdriver/pkg/iir/iirkal/phase2 \
  ../../../../testdriver/pkg/iir/iirlt \
  ../../../../testdriver/pkg/iir/iirlt/common \
  ../../../../testdriver/pkg/iir/iirlt/phase1 \
  ../../../../testdriver/pkg/iir/iirlt/phase2 \
  ../../../../testdriver/pkg/iir/iirbqnd \
  ../../../../testdriver/pkg/iir/iirbqnd/common \
  ../../../../testdriver/pkg/iir/iirbqnd/phase1 \
  ../../../../testdriver/pkg/iir/iirbqnd/phase2 \
  ../../../../testdriver/pkg/img \
  ../../../../testdriver/pkg/img/common \
  ../../../../testdriver/pkg/img/imgfft \
  ../../../../testdriver/pkg/img/imgfft/common \
  ../../../../testdriver/pkg/img/imgfft/phase1 \
  ../../../../testdriver/pkg/img/imgmisc \
  ../../../../testdriver/pkg/img/imgmisc/common \
  ../../../../testdriver/pkg/img/imgmisc/phase1 \
  ../../../../testdriver/pkg/img/imgresize \
  ../../../../testdriver/pkg/img/imgresize/common \
  ../../../../testdriver/pkg/img/imgresize/phase1 \
  ../../../../testdriver/pkg/img/imgrotate \
  ../../../../testdriver/pkg/img/imgrotate/common \
  ../../../../testdriver/pkg/img/imgrotate/phase1 \
  ../../../../testdriver/pkg/math \
  ../../../../testdriver/pkg/math/common \
  ../../../../testdriver/pkg/math/maths \
  ../../../../testdriver/pkg/math/maths/common \
  ../../../../testdriver/pkg/math/maths/phase1 \
  ../../../../testdriver/pkg/math/maths/phase2 \
  ../../../../testdriver/pkg/math/maths/phase3 \
  ../../../../testdriver/pkg/math/mathv \
  ../../../../testdriver/pkg/math/mathv/common \
  ../../../../testdriver/pkg/math/mathv/phase1 \
  ../../../../testdriver/pkg/math/mathv/phase2 \
  ../../../../testdriver/pkg/math/mathv/phase3 \
  ../../../../testdriver/pkg/math/mathvf \
  ../../../../testdriver/pkg/math/mathvf/common \
  ../../../../testdriver/pkg/math/mathvf/phase1 \
  ../../../../testdriver/pkg/matinv \
  ../../../../testdriver/pkg/matinv/chol \
  ../../../../testdriver/pkg/matinv/chol/common \
  ../../../../testdriver/pkg/matinv/chol/phase1 \
  ../../../../testdriver/pkg/matinv/chol/phase2 \
  ../../../../testdriver/pkg/matinv/common \
  ../../../../testdriver/pkg/matinv/gj \
  ../../../../testdriver/pkg/matinv/gj/common \
  ../../../../testdriver/pkg/matinv/gj/phase1 \
  ../../../../testdriver/pkg/matinv/gj/phase2 \
  ../../../../testdriver/pkg/matop \
  ../../../../testdriver/pkg/matop/common \
  ../../../../testdriver/pkg/matop/phase1 \
  ../../../../testdriver/pkg/matop/phase2 \
  ../../../../testdriver/pkg/vector \
  ../../../../testdriver/pkg/vector/common \
  ../../../../testdriver/pkg/vector/phase1 \
  ../../../../testdriver/pkg/vector/phase2 \
  ../../../../testdriver/pkg/vector/phase3

OUTDIR  =	../../../bin

INCDIRS =	../../../../library/include		\
		      ../../../../testdriver/include

PACKAGE_SUFFIX = ""

TARGET := testdriver$(PACKAGE_SUFFIX)
LIBS   := m
DLIBS  := ../../../bin/NatureDSP_Signal$(PACKAGE_SUFFIX)

BUILD_TYPE		:= EXECUTABLE
INTEGRIT_OLD_SUFFIXES	:= 1
COMPILE_CSTUB		:= 0
LSP			:=
EXTERNAL_PROJECTS = ../library

# check if c++ compilation required
ifeq (C++,$(LANG))
CPP = 1
HAVE_XTENSA = 0
else
CPP = 0
HAVE_XTENSA = 1
endif  

ifeq (,$(XTENSA_SYSTEM))
XTCORE := $(notdir $(wildcard ../../../../cstub/*-0.cpp))
XTCORE := $(subst -0.cpp,,$(XTCORE))
XTCORE := $(subst cstub-Xm_,,$(XTCORE))
else
XTCORE := $(XTENSA_CORE)
endif

EXTERNAL_PROJECTS.0 = ../cstub
EXTERNAL_PROJECTS += $(EXTERNAL_PROJECTS.$(HAVE_XTENSA))

ifeq (0,$(HAVE_XTENSA))
USE_PATCHES=1  # use patched cstubs for gcc, otherwise use REFERENCE variable
else
ifeq (1,$(REFERENCE))
USE_PATCHES=0
REF_SUFFIX=_reference
else
USE_PATCHES=1
REF_SUFFIX=
endif
endif

ifeq (1,$(MEM_MODEL))
MEM_SUFFIX=_mem_model
else
ifeq (2,$(MEM_MODEL))
MEM_SUFFIX=_mem_model2
else
MEM_SUFFIX=
endif
endif


MODULES = \
  main \
  addr2name \
  float16 \
  fpstat \
  malloc16 \
  mips \
  package \
  profiler \
  rms \
  testcase \
  testeng \
  testeng_errh \
  testeng_load_fxn \
  testeng_process_fxn \
  utils \
  vectools \
  vreport \
  mfcc \
  testeng_logmel \
  testeng_mfcc \
  func_mfcc1 \
  mfcc1 \
  mips_mfcc1 \
  func_mfcc2 \
  mfcc2 \
  mips_mfcc2 \
  complex \
  complexs \
  complexs2 \
  func_complexs2 \
  mips_complexs2 \
  complexv \
  complexv2 \
  func_complexv2 \
  mips_complexv2 \
  ef \
  ef1 \
  func_ef1 \
  mips_ef1 \
  cfft \
  test_cfft \
  cfft1 \
  func_cfft1 \
  mips_cfft1 \
  cfftie \
  test_cfftie \
  cfftie1 \
  func_cfftie1 \
  mips_cfftie1 \
  cfftie2 \
  func_cfftie2 \
  mips_cfftie2 \
  cnfft \
  cnfft1 \
  func_cnfft1 \
  mips_cnfft1 \
  fft \
  testeng_fft \
  dct \
  testeng_dct \
  dct1 \
  func_dct1 \
  mips_dct1 \
  dct2 \
  func_dct2 \
  mips_dct2 \
  rfft \
  test_rfft \
  func_rfft1 \
  mips_rfft1 \
  rfft1 \
  rfftie \
  test_rfftie \
  func_rfftie1 \
  mips_rfftie1 \
  rfftie1 \
  func_rfftie2 \
  mips_rfftie2 \
  rfftie2 \
  rnfft \
  func_rnfft1 \
  mips_rnfft1 \
  rnfft1 \
  spectrum \
  test_fft_spectrum \
  func_spectrum1 \
  mips_spectrum1 \
  spectrum1 \
  func_spectrum2 \
  mips_spectrum2 \
  spectrum2 \
  fir \
  testeng_fir \
  testeng_fir_old \
  conv2d \
  conv2d1 \
  func_conv2d1 \
  mips_conv2d1 \
  conv2d2 \
  func_conv2d2 \
  mips_conv2d2 \
  conv2d3 \
  func_conv2d3 \
  mips_conv2d3 \
  firblk \
  test_firblk \
  firblk1 \
  func_firblk1 \
  mips_firblk1 \
  firblk2 \
  func_firblk2 \
  mips_firblk2 \
  firdec \
  test_firdec \
  firdec1 \
  func_firdec1 \
  mips_firdec1 \
  firdec2 \
  func_firdec2 \
  mips_firdec2 \
  firint \
  test_firint \
  firint1 \
  func_firint1 \
  mips_firint1 \
  firint2 \
  func_firint2 \
  mips_firint2 \
  firother \
  test_firother \
  firother1 \
  func_firother1 \
  mips_firother1 \
  firother2 \
  func_firother2 \
  mips_firother2 \
  fit \
  pfit \
  func_pfit1 \
  mips_pfit1 \
  pfit1 \
  func_pfit2 \
  mips_pfit2 \
  pfit2 \
  iir \
  testeng_iir_lat \
  testeng_iir_old \
  iirbq \
  test_iirbq \
  func_iirbq1 \
  iirbq1 \
  mips_iirbq1 \
  func_iirbq2 \
  iirbq2 \
  mips_iirbq2 \
  iirkal \
  test_iirkal \
  func_iirkal1 \
  iirkal1 \
  mips_iirkal1 \
  func_iirkal2 \
  iirkal2 \
  mips_iirkal2 \
  iirlt \
  test_iirlt \
  func_iirlt1 \
  iirlt1 \
  mips_iirlt1 \
  func_iirlt2 \
  iirlt2 \
  mips_iirlt2 \
  iirbqnd \
  func_iirbqnd1 \
  iirbqnd1 \
  mips_iirbqnd1 \
  func_iirbqnd2 \
  iirbqnd2 \
  mips_iirbqnd2 \
  img \
  imgfft \
  func_imgfft1 \
  imgfft1 \
  mips_imgfft1 \
  imgmisc \
  func_imgmisc1 \
  imgmisc1 \
  mips_imgmisc1 \
  imgresize \
  func_imgresize1 \
  imgresize1 \
  mips_imgresize1 \
  imgrotate \
  func_imgrotate1 \
  imgrotate1 \
  mips_imgrotate1 \
  math \
  testeng_math \
  maths \
  func_maths1 \
  maths1 \
  mips_maths1 \
  acc_maths2 \
  func_maths2 \
  maths2 \
  mips_maths2 \
  acc_maths3 \
  func_maths3 \
  maths3 \
  mips_maths3 \
  mathv \
  func_mathv1 \
  mathv1 \
  mips_mathv1 \
  func_mathv2 \
  mathv2 \
  mips_mathv2 \
  func_mathv3 \
  mathv3 \
  mips_mathv3 \
  mathvf \
  func_mathvf1 \
  mathvf1 \
  mips_mathvf1 \
  chol \
  test_chol \
  chol1 \
  func_chol1 \
  mips_chol1 \
  chol2 \
  func_chol2 \
  mips_chol2 \
  matinv \
  testeng_matinv \
  gj \
  func_gj1 \
  gj1 \
  mips_gj1 \
  func_gj2 \
  gj2 \
  mips_gj2 \
  matop \
  testeng_matop \
  func_matop1 \
  matop1 \
  mips_matop1 \
  func_matop2 \
  matop2 \
  mips_matop2 \
  testeng_vector \
  vector \
  func_vector1 \
  mips_vector1 \
  vector1 \
  func_vector2 \
  mips_vector2 \
  vector2 \
  func_vector3 \
  mips_vector3 \
  vector3

ifneq (0, $(HAVE_XTENSA))
#----------------------------------------------------------
# CFLAGS_FILE_module_name   += some compiler switches
# LLVMFLAGS_FILE_module_name+= some compiler switches
# select specific compiler switches for xcc and llvm(clang)
# compilers respectively
#----------------------------------------------------------
#CFLAGS_FILE_module_name   += some compiler switches

CFLAGS_FILE_mips_img += -mtext-section-literals 
LLVMFLAGS_FILE_mips_img += -mtext-section-literals
LLVMFLAGS_FILE_fpstat += -fno-unsafe-math-optimizations

ifneq (,$(DRAM0))
DRAM0_MEM="$(DRAM0).data"
CFLAGS_FILE_mips += '-DDRAM0=$(DRAM0_MEM)'
endif
ifneq (,$(DRAM1))
DRAM1_MEM="$(DRAM1).data"
CFLAGS_FILE_mips += '-DDRAM1=$(DRAM1_MEM)'
endif
else
CFLAGS_FILE_testeng_errh += -DXTENSA_CORE=Xm_$(XTCORE)
endif

######### special definition for proper build under xws ##############
ifeq (0,0)
CFLAGS_FILE_mips += 
LLVMFLAGS_FILE_mips += 
endif
 
default:	all

.PHONY:	all clean cleanest cleanlib debug ferret default profile release
.PHONY: $(foreach proj,$(EXTERNAL_PROJECTS), $(proj))
.PHONY: $(foreach proj,$(EXTERNAL_PROJECTS), $(proj)/clean)
.PHONY: $(foreach proj,$(EXTERNAL_PROJECTS), $(proj)/cleanest)

debug:		all
ferret:		all
profile:	all
release:	all 
 

ifeq (0,$(CPP)) ######################################
ifeq (LLVM,$(LANG))
PLATFORM     := $(XTCORE)$(REF_SUFFIX)$(MEM_SUFFIX)_llvm
PLATFORM_LIB := $(XTCORE)$(REF_SUFFIX)$(MEM_SUFFIX)_llvm
else
PLATFORM     := $(XTCORE)$(REF_SUFFIX)$(MEM_SUFFIX)
PLATFORM_LIB := $(XTCORE)$(REF_SUFFIX)$(MEM_SUFFIX)
endif
else #################################################
PLATFORM     := $(XTCORE)$(MEM_SUFFIX)-C++
PLATFORM_LIB := $(XTCORE)$(MEM_SUFFIX)-C++
endif ################################################

ifneq (, $(findstring cstub, $(PLATFORM))) ###########
SOLID_LIBRARY		:= 1
CCFLAGS_DIR_testdriver	+= -Wno-uninitialized
CCFLAGS_DIR_testdriver	+= -Wno-unused-result
CCFLAGS_DIR_cstub 		+= -Wno-uninitialized
CCFLAGS					+= -Wno-strict-aliasing
else #################################################
SOLID_LIBRARY		:= 0
CCFLAGS_DIR_fft_files += -OPT:alias=restrict
CCFLAGS_DIR_P2_files += -LNO:simd
endif ################################################

###############################################################
##                      Host OS detection                    ##
###############################################################


WINDIR+=$(windir)

ifeq (,$(WINDIR))
# Unix-like system
ULIS=1
override HOSTOS=$(shell uname)
else
ULIS=0
override HOSTOS=win32
SHELL=cmd
endif

###############################################################
##              Platform-specific configuration              ##
###############################################################
CSTUB_MODULES += cstub-Xm_$(XTCORE)

ifeq (,$(XTENSA_SYSTEM))
CSTUB_INCDIRS	+= ../../../../cstub
CSTUB_SRCDIRS	+= ../../../../cstub
ifeq (LLVM,$(LANG))
LLVMFLAGS_FILE_$(subst -,_,$(CSTUB_MODULES)) = -Dasm=__asm__ -O0
else
CFLAGS_FILE_$(subst -,_,$(CSTUB_MODULES)) = -Dasm=__asm__ -O0
endif
else
CSTUB_INCDIRS += $(XTENSA_SYSTEM)/../src/cstub $(XTENSA_SYSTEM)/../src/cstub/xtensa/tie
CSTUB_SRCDIRS += $(XTENSA_SYSTEM)/../src/cstub
endif

ifneq (0, $(HAVE_XTENSA)) ############################
##############################################################
VARPREFIX.$(PLATFORM)          = $(PLATFORM)
ifeq (LLVM,$(LANG))
CC.$(PLATFORM)      		     = xt-clang
LD.$(PLATFORM)      		     = xt-clang
else
CC.$(PLATFORM)      		     = xt-xcc 
LD.$(PLATFORM)      		     = xt-xcc
endif
AR.$(PLATFORM)		       = xt-ar
MAKE.$(PLATFORM)	       = xt-make
OBJCOPY.$(PLATFORM)	       = xt-objcopy
RUN.$(PLATFORM)		       = xt-run
RUNMEM.$(PLATFORM)	       = xt-run --mem_model
RUNFAST.$(PLATFORM)	       = xt-run --turbo
PROF.$(PLATFORM)	       = xt-prof
ARCH.$(PLATFORM)	       = $(PLATFORM)
CFLAGS.$(PLATFORM)	       += -fmessage-length=0 -mcoproc -mlongcalls  -fno-builtin
CFLAGS.profile.$(PLATFORM)     += -O3 
CFLAGS.release.$(PLATFORM)     += -O3 
ifeq (LLVM,$(LANG))
else
CFLAGS.$(PLATFORM)	       += -INLINE:requested
CFLAGS.profile.$(PLATFORM)     += -LNO:simd
CFLAGS.release.$(PLATFORM)     += -LNO:simd
endif
CFLAGS.$(PLATFORM)     += $(CFLAGS_EXT)
LFLAGS.profile.$(PLATFORM)     += --profile
RUNFLAGS.debug.$(PLATFORM)     = --turbo
RUNFLAGS.release.$(PLATFORM)   = --turbo
RUNFLAGS.ferret.$(PLATFORM)   = --turbo
TARGETOS.$(PLATFORM)	       = Xtensa
else
VARPREFIX.$(PLATFORM)		 = $(PLATFORM)
ifeq (1, $(CPP))
CC.$(PLATFORM)      		 = g++
LD.$(PLATFORM)      		 = g++
LFLAGS.$(PLATFORM) 		 += -Wl,--section-start=unused_section=0
CFLAGS.$(PLATFORM)		+=-D__USE_XOPEN2K8
else
VARPREFIX.$(PLATFORM)		 = $(PLATFORM)
CC.$(PLATFORM)      		 = gcc
LD.$(PLATFORM)      		 = gcc
LFLAGS.$(PLATFORM) 		 +=  -Wl,--section-start=unused_section=0
endif
AR.$(PLATFORM)      		 = ar
#MAKE.$(PLATFORM)      		 = make
OBJCOPY.$(PLATFORM)              = objcopy
RUN.$(PLATFORM)     		 = 
PROF.$(PLATFORM)      		 = 
CFLAGS.$(PLATFORM)         	 += -m32
CFLAGS.$(PLATFORM) 		 += -DXTENSA_CORE=Xm_$(XTCORE)
ifeq (XRC_B64_LP,$(XTCORE))
CFLAGS.$(PLATFORM)         	 += -DB64_LP
endif
INCDIRS			 	 += $(CSTUB_INCDIRS)
#SRCDIRS			 	 += $(CSTUB_SRCDIRS)
ifeq (1, $(COMPILE_CSTUB))
MODULES				 += $(CSTUB_MODULES)
endif
CFLAGS.profile.$(PLATFORM)       += -O3 
CFLAGS.release.$(PLATFORM)	 += -O3 
RUNFLAGS.debug.$(PLATFORM)	 = 
RUNFLAGS.release.$(PLATFORM)	 = 
RUNFLAGS.ferret.$(PLATFORM)	 = 
TARGETOS.$(PLATFORM)		 = Linux
LFLAGS.$(PLATFORM) 		 += -m32
DLIBS  += ../../../bin/cstub
endif

CFLAGS.release.$(PLATFORM)	 += -DUSE_ISA_PATCHES=$(USE_PATCHES)
###############################################################
##                     Well-known suffixes                   ##
###############################################################
#POSTFIX.$(PLATFORM)-Xtensa-release.1=_$(PLATFORM)
#POSTFIX.$(PLATFORM)-Xtensa-debug.1=_$(PLATFORM)
#POSTFIX.$(PLATFORM)-Xtensa-ferret.1=_$(PLATFORM)
#POSTFIX.$(PLATFORM)-Xtensa-profile.1=_$(PLATFORM)

#POSTFIX_DLC.$(PLATFORM)-Xtensa-release.1=_$(PLATFORM_LIB)
#POSTFIX_DLC.$(PLATFORM)-Xtensa-debug.1=_$(PLATFORM_LIB)
#POSTFIX_DLC.$(PLATFORM)-Xtensa-profile.1=_$(PLATFORM_LIB)

###############################################################
##                 Host-specific configuration               ##
###############################################################

ifeq (1,$(ULIS)) ##############################################
ECHO=echo
override RM=rm -f
CP=cp
PSEP=/
NULLDEV=/dev/null
MKDIR=mkdir
JOIN=&&
else ##########################################################
ECHO=@echo
override RM=del
CP=copy
PSEP=\\\\#
NULLDEV=NUL
MKDIR=mkdir
JOIN=&
endif #########################################################

###############################################################
##                     Target OS definition                  ##
###############################################################

TARGETOS=$(firstword $(TARGETOS.$(PLATFORM)) $(HOSTOS))

VARPREFIX.Darwin=MAC
VARPREFIX.iOS=MAC
VARPREFIX.linux=LINUX
VARPREFIX.Linux=LINUX
VARPREFIX.unix-like.1=UNIX_LIKE
VARPREFIX.win32=WIN32

###############################################################
##                  Select build configuration               ##
###############################################################

CONFIG=$(firstword $(findstring release, $(MAKECMDGOALS)) $(findstring debug,\
 $(MAKECMDGOALS)) $(findstring ferret, $(MAKECMDGOALS)) $(findstring profile, $(MAKECMDGOALS)) release)

CONFIG_ID	= $(PLATFORM)-$(TARGETOS)-$(CONFIG)
DLIBS_CONFIG_ID	= $(PLATFORM_LIB)-$(TARGETOS)-$(firstword \
		  $(DLIBS_CONFIG) $(CONFIG))

# Intermediate directory

INTERMEDIATE_DIR= $(CONFIG_ID)
USE_SUFFIXES	= $(firstword $(INTEGRIT_SUFFIXES) $(INTEGRIT_OLD_SUFFIXES))
POSTFIX		= $(firstword \
	$(POSTFIX.$(CONFIG_ID).$(USE_SUFFIXES)) -$(CONFIG_ID))
POSTFIX_DLC	= $(firstword \
	$(POSTFIX_DLC.$(DLIBS_CONFIG_ID).$(USE_SUFFIXES)) \
	$(POSTFIX.$(DLIBS_CONFIG_ID).$(USE_SUFFIXES)) -$(DLIBS_CONFIG_ID))

$(INTERMEDIATE_DIR)/dirinfo.txt:
	-$(MKDIR) $(INTERMEDIATE_DIR)
	-$(ECHO) Intermediate directory.>$(INTERMEDIATE_DIR)/dirinfo.txt
	
###############################################################
##                Creating devtools options                  ##
###############################################################

IDIRS_OPT = $(foreach f,$(INCDIRS),-I $(f))
LDIRS_OPT = $(foreach f,$(LIBDIRS),-L $(f))

CC			= $(firstword $(CC.$(PLATFORM)) gcc)
LD			= $(firstword $(LD.$(PLATFORM)) gcc)
AR			= $(firstword $(AR.$(PLATFORM)) ar)
#MAKE			= $(firstword $(MAKE.$(PLATFORM)) make)
OBJCOPY			= $(firstword $(OBJCOPY.$(PLATFORM)) objcopy)
RUN			= $(RUN.$(PLATFORM))
PROF			= $(firstword $(PROF.$(PLATFORM)) gprof)

CFLAGS_DEPENDENT	+= $(CFLAGS.$(PLATFORM))

CFLAGS.default.debug.	= $(CFLAGS.default.debug) -g -O0 -D_DEBUG
CFLAGS.default.profile.	= $(CFLAGS.default.profile) -g -pg -O3
CFLAGS.default.release.	= $(CFLAGS.default.release) -O3
CFLAGS.default.ferret.	= $(CFLAGS.default.ferret) -g -O0 -D_DEBUG
LFLAGS.default.debug.	= $(LFLAGS.default.debug)
LFLAGS.default.ferret.	= $(LFLAGS.default.ferret)
LFLAGS.default.profile.	= $(LFLAGS.default.profile) -g -pg -O
LFLAGS.default.release.EXECUTABLE = -s 
LFLAGS.default.ferret.EXECUTABLE =  
LFLAGS.default.release.	= $(LFLAGS.default.release) -O \
			$(LFLAGS.default.release.$(BUILD_TYPE))
#ifeq (1,$(MAP))
LFLAGS += -Wl,-Map,$(TARGET_SPEC.EXECUTABLE).map
#endif

# support memory modeling
ifeq (1,$(HAVE_XTENSA))
ifeq (1,$(MEM_MODEL))
LFLAGS += -mlsp=sim-local
CFLAGS += -DMEM_MODEL=1
else
ifeq (2,$(MEM_MODEL))
LFLAGS += -mlsp=sim
CFLAGS += -DMEM_MODEL=2
else
endif
endif
endif

CFLAGS.debug		= $(CFLAGS.debug.$(PLATFORM)) \
			$(CFLAGS.default.debug.$(CFLAGS.debug.$(PLATFORM)))
CFLAGS.ferret		= $(CFLAGS.ferret.$(PLATFORM)) \
			$(CFLAGS.default.ferret.$(CFLAGS.ferret.$(PLATFORM)))			
CFLAGS.profile		= $(CFLAGS.profile.$(PLATFORM)) \
			$(CFLAGS.default.profile.$(CFLAGS.profile.$(PLATFORM)))
CFLAGS.release		= $(CFLAGS.release.$(PLATFORM)) \
			$(CFLAGS.default.release.$(CFLAGS.release.$(PLATFORM)))

CFLAGADDINCDIR		= $(firstword $(CFLAGADDINCDIR.$(PLATFORM)) -I)
CFLAGMAKEDEP		= $(firstword $(CFLAGMAKEDEP.$(PLATFORM)) -MM)
CFLAGMAKEOBJ		= $(firstword $(CFLAGMAKEDEP.$(PLATFORM)) -c)
CFLAGSAVETEMP.1		= $(firstword $(CFLAGSAVETEMP.$(PLATFORM)) -save-temps)

LFLAGS_DEPENDENT	+= $(LFLAGS.$(PLATFORM))
LFLAGS.debug		= $(LFLAGS.debug.$(PLATFORM)) \
			$(LFLAGS.default.debug.$(LFLAGS.debug.$(PLATFORM)))
LFLAGS.profile		= $(LFLAGS.profile.$(PLATFORM)) \
			$(LFLAGS.default.profile.$(LFLAGS.profile.$(PLATFORM)))
LFLAGS.release		= $(LFLAGS.release.$(PLATFORM)) \
			$(LFLAGS.default.release.$(LFLAGS.release.$(PLATFORM)))
LFLAGS.ferret		= $(LFLAGS.ferret.$(PLATFORM)) \
			$(LFLAGS.default.ferret.$(LFLAGS.ferret.$(PLATFORM))) -Wl,-u -Wl,malloc -lferret
LFLAGS.LIBRARY.1	+=--no-standard-libraries
LFLAGS.LIBRARY.1	+=-Wl,-r,--retain-symbols-file,$(TARGET_SYMFILE)

ARFLAGS			= $(firstword $(ARFLAGS.$(PLATFORM)) rcs)

RUNFLAGS		+= $(RUNFLAGS.$(PLATFORM))
RUNFLAGS		+= $(RUNFLAGS.$(CONFIG).$(PLATFORM))

PROFFLAGS		+= $(PROFFLAGS.$(PLATFORM))

LFLAGADDLIBDIR		= $(firstword $(LFLAGADDLIBDIR.$(PLATFORM)) -L)

INCDIRS			+= $($(VARPREFIX.$(PLATFORM))_INCDIRS)
SRCDIRS			+= $($(VARPREFIX.$(PLATFORM))_SRCDIRS)
LIBDIRS			+= $($(VARPREFIX.$(PLATFORM))_LIBDIRS)
SRCDIRS			+= $(SRCDIR)

# Header search options
CFLAGS_DEPENDENT += $(foreach idir, $(INCDIRS), $(CFLAGADDINCDIR) $(idir))

CFLAGS_DEPENDENT	+= $(CFLAGSAVETEMP.$(SAVETEMPS))

# Library search options
LFLAGS_DEPENDENT += $(foreach ldir, $(LIBDIRS), $(LFLAGADDLIBDIR) $(ldir))

# Configuration-dependent options
CFLAGS_DEPENDENT += $(CFLAGS.$(CONFIG))
LFLAGS_DEPENDENT += $(LFLAGS.$(CONFIG)) $(LFLAGS.$(BUILD_TYPE).$(SOLID_LIBRARY))
ifeq (0, $(CPP))
#CFLAGS_DEPENDENT += -std=c99
endif
CFLAGS_DEPENDENT += $(IDIRS_OPT) -c -Dasm=__asm__
CFLAGS_DEPENDENT += -Wall -Wno-error 
CFLAGS_DEPENDENT += '-DFULL_VECTOR_DIR="./../../vectors_full"'
CFLAGS_DEPENDENT += '-DBRIEF_VECTOR_DIR="./../../vectors_brief"'
CFLAGS_DEPENDENT += '-DSANITY_VECTOR_DIR="./../../testdriver/vectors_sanity"'
CFLAGS_DEPENDENT += '-DPACKAGE_SUFFIX=$(PACKAGE_SUFFIX)'
###############################################################
##                     Create target spec                    ##
###############################################################

#TARGET_SUFFIX.win32	= .exe

TARGET_SUFFIX		= $(firstword $(TARGET_SUFFIX.$(BUILD_TYPE))\
			  $(TARGET_SUFFIX.$(HOSTOS)))

TARGET_SPEC.LIBRARY	= $(OUTDIR)/$(TARGET)$(POSTFIX).a
TARGET_SPEC.EXECUTABLE	= $(OUTDIR)/$(TARGET)$(POSTFIX)$(TARGET_SUFFIX)

TARGET_INTERMEDIATE.LIBRARY.1.1	= $(INTERMEDIATE_DIR)/$(TARGET).solid.1.o
TARGET_INTERMEDIATE.LIBRARY.1.2	= $(INTERMEDIATE_DIR)/$(TARGET).solid.2.o
TARGET_INTERMEDIATE.LIBRARY.1.3	= $(INTERMEDIATE_DIR)/$(TARGET).solid.o


###############################################################
##                    Complete module list                   ##
###############################################################

MODULES += $($(VARPREFIX.$(PLATFORM))_MODULES)
MODULES += $($(VARPREFIX.$(TARGETOS))_MODULES)
MODULES += $($(VARPREFIX.unix-like.$(ULIS))_MODULES)
# Keep compatibility with IntegrIT makefiles below 0.30
MODULES += $(OBJS:.o=)
MODULES += $($(VARPREFIX.$(PLATFORM))_OBJS:.o=)
MODULES += $($(VARPREFIX.$(TARGETOS))_OBJS:.o=)
MODULES += $($(VARPREFIX.unix-like.$(ULIS))_OBJS:.o=)

###############################################################
##           Generate list of objects to be built            ##
###############################################################

ALL_OBJECTS = $(foreach module, $(MODULES), $(INTERMEDIATE_DIR)/$(module).o)

###############################################################
##        Generate list of extra libraries to be linked      ##
###############################################################

ALL_DLIBS=$(foreach lib,$(DLIBS),$(lib)$(POSTFIX_DLC).a)

###############################################################
##                     Dependency generator                  ##
###############################################################

ifneq (1,$(CPP))
vpath %.c $(SRCDIRS) $(CSTUB_SRCDIRS)
vpath %.cpp $(SRCDIRS) $(CSTUB_SRCDIRS)
else
vpath %.c $(SRCDIRS)
vpath %.cpp $(SRCDIRS) $(CSTUB_SRCDIRS)
endif

COMMAND_DEPGEN =  $(CC)
COMMAND_DEPGEN += $(CFLAGS_DEPENDENT)
COMMAND_DEPGEN += $(CFLAGS_DIR_$(subst -,_,$(notdir $(dir $(1)))))
ifeq (LLVM,$(LANG))
COMMAND_DEPGEN += $(LLVMFLAGS_FILE_$(subst -,_,$(notdir $(basename $(1)))))
else
COMMAND_DEPGEN += $(CFLAGS_FILE_$(subst -,_,$(notdir $(basename $(1)))))
endif
COMMAND_DEPGEN += $(CFLAGS)
COMMAND_DEPGEN += $(CFLAGMAKEDEP)
COMMAND_DEPGEN += -MQ $(2)
COMMAND_DEPGEN += -MQ $(2:.d=.o)
COMMAND_DEPGEN += $(1)
COMMAND_DEPGEN += -o $(2)

$(INTERMEDIATE_DIR)/%.d: %.c $(INTERMEDIATE_DIR)/dirinfo.txt
	$(call COMMAND_DEPGEN,$<,$@)

$(INTERMEDIATE_DIR)/%.d: %.cpp $(INTERMEDIATE_DIR)/dirinfo.txt
	$(call COMMAND_DEPGEN,$<,$@)

# Include dependencies
-include $(addprefix $(INTERMEDIATE_DIR)/, $(addsuffix .d, $(MODULES)))

###############################################################
##                        Compiling rules                    ##
###############################################################

COMMAND_COMPILE =  $(CC)
COMMAND_COMPILE += $(CFLAGS_DEPENDENT)
COMMAND_COMPILE += $(CFLAGS_DIR_$(subst -,_,$(notdir $(dir $(1)))))
ifeq (LLVM,$(LANG))
COMMAND_COMPILE += $(LLVMFLAGS_FILE_$(subst -,_,$(notdir $(basename $(1)))))
else
COMMAND_COMPILE += $(CFLAGS_FILE_$(subst -,_,$(notdir $(basename $(1)))))
endif
COMMAND_COMPILE += $(CFLAGS)
COMMAND_COMPILE += $(CFLAGMAKEOBJ)
COMMAND_COMPILE += $(1)
COMMAND_COMPILE += -o $(2) 

$(INTERMEDIATE_DIR)/%.o: %.c $(INTERMEDIATE_DIR)/dirinfo.txt
	$(call COMMAND_COMPILE,$<,$@)

$(INTERMEDIATE_DIR)/%.o: %.cpp $(INTERMEDIATE_DIR)/dirinfo.txt
	$(call COMMAND_COMPILE,$<,$@)

###############################################################
##                    Linking/archiving rules                ##
###############################################################

$(foreach lib,$(ALL_DLIBS),$(lib)): $(EXTERNAL_PROJECTS)

ALL_LIBS = $(ALL_DLIBS) $(foreach lib,$(LIBS),-l$(lib))

LD_LIST.  = $(ALL_OBJECTS)
LD_LIST.0 = $(ALL_OBJECTS)
LD_LIST.1 = $(TARGET_INTERMEDIATE.$(BUILD_TYPE).$(SOLID_LIBRARY).3)

$(TARGET_INTERMEDIATE.LIBRARY.1.1): $(ALL_OBJECTS) $(ALL_DLIBS)
	$(LD) $(LFLAGS_DEPENDENT) $(LFLAGS) -o \
	$(TARGET_INTERMEDIATE.LIBRARY.1.1) $(ALL_OBJECTS) $(ALL_LIBS)

$(TARGET_INTERMEDIATE.LIBRARY.1.2):	$(TARGET_INTERMEDIATE.LIBRARY.1.1)
	$(OBJCOPY) --keep-global-symbols=$(TARGET_SYMFILE) \
	$(TARGET_INTERMEDIATE.LIBRARY.1.1) $(TARGET_INTERMEDIATE.LIBRARY.1.2)

$(TARGET_INTERMEDIATE.LIBRARY.1.3):	$(TARGET_INTERMEDIATE.LIBRARY.1.2)
	$(OBJCOPY) --strip-unneeded \
	$(TARGET_INTERMEDIATE.LIBRARY.1.2) $(TARGET_INTERMEDIATE.LIBRARY.1.3)

$(INTERMEDIATE_DIR)/$(TARGET).ld: $(LD_LIST.$(SOLID_LIBRARY))
	$(shell $(RM) $@) $(foreach mod,$(LD_LIST.$(SOLID_LIBRARY)),$(shell $(ECHO) $(mod)>>$@))

$(TARGET_SPEC.LIBRARY): $(EXTERNAL_PROJECTS) $(LD_LIST.$(SOLID_LIBRARY)) $(INTERMEDIATE_DIR)/$(TARGET).ld
	$(RM) $@
	$(AR) rcs $@ @$(INTERMEDIATE_DIR)/$(TARGET).ld

$(TARGET_SPEC.EXECUTABLE): $(EXTERNAL_PROJECTS) $(ALL_OBJECTS) $(ALL_DLIBS)
	$(LD) $(LFLAGS_DEPENDENT) $(LFLAGS) -o \
	$(TARGET_SPEC.EXECUTABLE) $(ALL_OBJECTS) $(ALL_LIBS)

all:    $(TARGET_SPEC.$(BUILD_TYPE))

$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)):
	+$(MAKE) -C $@ CONFIG=$(CONFIG) all

###############################################################
##                         Clean rules                       ##
###############################################################
$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)/clean):
	$(MAKE) -C $(dir $@) CONFIG=$(CONFIG) clean
$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)/cleanest):
	$(MAKE) -C $(dir $@) CONFIG=$(CONFIG) cleanest

clean:	$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)/clean)
#clean:	
	-$(RM) $(subst /,$(PSEP),$(INTERMEDIATE_DIR)/*.o)
	-$(RM) $(subst /,$(PSEP),$(TARGET_SPEC.$(BUILD_TYPE)))

cleanest:	$(foreach proj,$(EXTERNAL_PROJECTS),$(proj)/cleanest)
	-$(RM) $(subst /,$(PSEP),$(INTERMEDIATE_DIR)/*.d)
	-$(RM) $(subst /,$(PSEP),$(INTERMEDIATE_DIR)/*.o)
	-$(RM) $(subst /,$(PSEP),$(TARGET_SPEC.$(BUILD_TYPE)))

cleanlib:
	-$(RM) $(subst /,$(PSEP),$(TARGET_SPEC.LIBRARY)) 2>$(NULLDEV)
