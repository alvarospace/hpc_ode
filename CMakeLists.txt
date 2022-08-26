cmake_minimum_required(VERSION 3.0.0)
project(HPC_ODE VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

########### Dependencies #########
include(ExternalProject)


# Cantera
set(CANTERA_PREFIX ${CMAKE_SOURCE_DIR}/ext_install/cantera CACHE PATH "Path where cantera will be installed")
ExternalProject_Add(CanteraProject
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/ext/cantera
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND  ""
    BUILD_COMMAND scons build prefix=${CANTERA_PREFIX} python_package=none
    INSTALL_COMMAND scons install
)
set(CANTERA_LIB_DIR ${CANTERA_PREFIX}/lib)
set(CANTERA_INCLUDE_DIR ${CANTERA_PREFIX}/include)

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
