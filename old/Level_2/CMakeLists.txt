cmake_minimum_required(VERSION 3.22.1)

# Project name
project(Standalone_Integrator LANGUAGES C CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)
# enable_language(CXX)
# set(CMAKE_CXX_FLAGS "-std=c++0x -O3")
# set(CMAKE_CXX_STANDARD 11)

# Configuration variable for sundials path
set(SUNDIALS_EXTERN_ENABLE OFF CACHE BOOL "Set ON to specify SUNDIALS install path")
set(SUNDIALS_EXTERN_LIBS "/usr/lib/sundials/lib/" CACHE PATH "SUNDIALS library path") 
set(SUNDIALS_EXTERN_INCLUDE_PATH "/usr/lib/sundials/include/" CACHE PATH "SUNDIALS include path")

message("SUNDIALS EXTERN: " ${SUNDIALS_EXTERN_ENABLE})
message("SUNDIALS EXTERN LIBS: " ${SUNDIALS_EXTERN_LIBS})
message("SUNDIALS EXTERN INCLUDE PATH: " ${SUNDIALS_EXTERN_INCLUDE_PATH})

# Pyjac Library
add_subdirectory(libpyjac)

# Serial app
add_subdirectory(app_serial)

# Serial app alvaro
add_subdirectory(app_serial_enhanced)

# OpenMP version
add_subdirectory(app_openmp)

# Installation
# install(TARGETS app_serial DESTINATION bin)
# install(TARGETS libpyjac DESTINATION lib)
