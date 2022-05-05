cmake_minimum_required(VERSION 3.1)
project (run)

set(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(CXX)
set(CMAKE_CXX_FLAGS "-std=c++0x -O3")
set(CMAKE_CXX_STANDARD 11)

find_package(Threads REQUIRED)

#Sundials Include 
include_directories("/path_to_sundials/sundials_6.1/build/include/")

#Sundials library 
link_libraries(-L/path_to_sundials/sundials_6.1/build/lib64/ -lsundials_cvodes -lsundials_nvecserial)
link_directories("/path_to_sundials/sundials_6.1/build/lib64/")

#Pyjac Includes 
include_directories("./out/")
#Pyjac library 
link_libraries(./libc_pyjac.a)

add_executable(run_pyjac main.cpp)
