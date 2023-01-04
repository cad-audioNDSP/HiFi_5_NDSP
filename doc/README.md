# How to Build and Run the Source Code in Linux environment
	Get the latest or required version of NDSP HiFi5 Code from GitHub 
  https://github.com/cad-audioNDSP/HiFi_5_NDSP/tree/main/NDSP_HiFi5

# The source code is organized as follows.
  build - contains the make file 
  library - contains the optimized kernel functions for the HiFi core 
  testdriver - contains the demo driver code to tun the library   

# It is assumed that the required HiFi core configurations and the Xtensa toolchain are set up in the Linux environment. 
 	An example .cshrc file  that sets up the build environment accordingly is provided for reference.   

# Setting up the environment 
  A typical way is to place this .cshrc file in your home directory. 
  Source ~/.cshrc 
  setenv XTENSA_CORE AE_HiFi5e_LE5_AO_FP

# Compiling the Source Code: 
  •	Navigate to the testdriver directory: 
  …/ NDSP_HiFi5/build/project/xtclang/testdriver
  •	make clean -j -e LANG=LLVM  
  •	make all -j -e LANG=LLVM 


# Running the executable: 
•	Navigate to the bin directory: 
…/ NDSP_HiFi5/build/bin
•	Performance tests:
  o	xt-run testdriver-AE_HiFi5e_LE5_AO_FP_llvm-Xtensa-release -mips -sanity
  o	xt-run testdriver-AE_HiFi5e_LE5_AO_FP_llvm-Xtensa-release -mips -brief 
  o	xt-run testdriver-AE_HiFi5e_LE5_AO_FP_llvm-Xtensa-release -mips -full   
•	Functional tests:
  o	xt-run testdriver-AE_HiFi5e_LE5_AO_FP_llvm-Xtensa-release -func -sanity
  o	xt-run testdriver-AE_HiFi5e_LE5_AO_FP_llvm-Xtensa-release -func -brief
  o	xt-run testdriver-AE_HiFi5e_LE5_AO_FP_llvm-Xtensa-release -func -full
  o	xt-run testdriver-AE_HiFi5e_LE5_AO_FP_llvm-Xtensa-release -func -sanity -verbose 
  o	xt-run testdriver-AE_HiFi5e_LE5_AO_FP_llvm-Xtensa-release -func -sanity -fir -verbose 
  o	xt-run testdriver-AE_HiFi5e_LE5_AO_FP_llvm-Xtensa-release -func -brief -fir -iir -fft
