cmake_minimum_required(VERSION 3.1)
project (run)

set(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(CXX)
set(CMAKE_CXX_FLAGS "-std=c++0x -O3")
set(CMAKE_CXX_STANDARD 11)

find_package(Threads REQUIRED)

#Sundials
include_directories("/gpfs/projects/bsc21/bsc21575/02_Sources/sundials_6.1/build/include/")
link_directories("/gpfs/projects/bsc21/bsc21575/02_Sources/sundials_6.1/build/lib64/")
link_libraries(-L/gpfs/projects/bsc21/bsc21575/02_Sources/sundials_6.1/build/lib64/ -lsundials_cvodes -lsundials_nvecserial)

#Pyjac
include_directories("/gpfs/projects/bsc21/bsc21575/03_Random_tests/08_Random_Scripts/21_PyJac/libs_pyjac/out_est2/")
link_libraries(/gpfs/projects/bsc21/bsc21575/03_Random_tests/08_Random_Scripts/21_PyJac/libs_pyjac/out_est2/libc_pyjac.a)

add_executable(run_pyjac main.cpp)
