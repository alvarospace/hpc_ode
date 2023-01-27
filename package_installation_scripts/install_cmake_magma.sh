#!/bin/bash

cmake -DMAGMA_ENABLE_CUDA=ON \
	-DMAGMA_ENABLE_HIP=OFF \
	-DCMAKE_INSTALL_PREFIX=/home/almousa/install/magma \
	-DGPU_TARGET='Ampere' \
	..
