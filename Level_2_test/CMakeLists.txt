cmake_minimum_required(VERSION 3.22.1)

# Project name
project(Standalone_Integrator LANGUAGES C CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)
# enable_language(CXX)
# set(CMAKE_CXX_FLAGS "-std=c++0x -O3")
# set(CMAKE_CXX_STANDARD 11)


# Pyjac Library
add_subdirectory(libpyjac)

# Serial app
add_subdirectory(app_serial)

# Installation
# install(TARGETS app_serial DESTINATION bin)
# install(TARGETS libpyjac DESTINATION lib)
