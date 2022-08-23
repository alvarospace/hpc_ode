cmake_minimum_required(VERSION 3.0.0)
project(HPC_ODE VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

########### Dependencies #########
include(ExternalProject)

ExternalProject_Add(cantera
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/ext/cantera
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND  ""
    BUILD_COMMAND scons build prefix=../../ext_install/cantera python_package=none
    INSTALL_COMMAND scons install
)

########### MAIN PROJECT #########

add_subdirectory(src)



############ TESTS ##########
include(CTest)
enable_testing()

add_subdirectory(tests)


########## PACK ##########
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
