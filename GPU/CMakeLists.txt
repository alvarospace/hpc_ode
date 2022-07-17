cmake_minimum_required(VERSION 3.22.1)

# Project name
project(GPU_Integrator LANGUAGES CXX CUDA)

set(CMAKE_VERBOSE_MAKEFILE ON)


# Set SUNDIALS_DIR variable pointing to config file to be able to use "find_package"
#set(SUNDIALS_EXTERN_ENABLE OFF CACHE BOOL "Set ON to specify SUNDIALS install path")
set(SUNDIALS_EXTERN_ENABLE OFF CACHE BOOL "Set ON to specify SUNDIALS cmake config file")
set(SUNDIALS_EXTERN_CMAKE_FILE "/usr/lib/sundials/lib/cmake/sundials" CACHE PATH "SUNDIALS cmake config path")

if(SUNDIALS_EXTERN_ENABLE)
    set(SUNDIALS_DIR ${SUNDIALS_EXTERN_CMAKE_FILE} CACHE PATH "Location of SUNDIALSConfig.cmake")
else()
    set(SUNDIALS_DIR /home/almousa/install/sundials-6.2/lib/cmake/sundials CACHE PATH "Location of SUNDIALSConfig.cmake")
endif()

message("SUNDIALS EXTERN: " ${SUNDIALS_EXTERN_ENABLE})
message("SUNDIALS CMAKE PATH: " ${SUNDIALS_DIR})

#set(SUNDIALS_DIR /home/almousa/install/sundials-6.2/lib/cmake/sundials CACHE PATH "Location of SUNDIALSConfig.cmake")


# Pyjac CUDA Library
add_subdirectory(libpyjacGPU)

# CUDA Integrator
add_subdirectory(app_magma)

# CUDA Integrator without pyjac library
add_subdirectory(app_magma_diff_jac)
