#!/bin/bash

scons gpu compute_level='sm_52' \
    sundials_inc_dir='/home/almousa/install/sundials-6.2/include' \
	sundials_lib_dir='/home/almousa/install/sundials-6.2/lib' \
	fftw3_inc_dir='/home/almousa/install/fftw/include' \
	fftw3_lib_dir='/home/almousa/install/fftw/lib' \
	mechanism_dir='/home/almousa/TFM/hpc_cvode/accelerInt/accelerInt-master/mechanism' \
	NVCCFLAGS='-m64 -Xptxas -v -Xptxas --disable-optimizer-constants -I/home/almousa/install/cuda-samples/Common --expt-relaxed-constexpr' \
	blas_lapack_libs='lapack,blas'