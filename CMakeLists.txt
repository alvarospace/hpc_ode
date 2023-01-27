cmake_minimum_required(VERSION 3.0.0)
project(HPC_ODE VERSION 0.1.0 LANGUAGES C CXX CUDA)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
enable_language(Fortran)

########### Dependencies #########
include(ExternalProject)
set(EXT_DIR ${CMAKE_SOURCE_DIR}/ext)
set(EXT_INSTALL_DIR ${CMAKE_SOURCE_DIR}/ext_install CACHE PATH "Path where dependencies will be installed")


# Cuda Samples
ExternalProject_Add(cuda-samplesProject
    SOURCE_DIR ${EXT_DIR}/cuda-samples
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND  ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
set(CUDA_SAMPLES_INCLUDE_DIR ${EXT_DIR}/cuda-samples/Common)

# Magma v2.6.2
set(GPU_TARGET Ampere CACHE STRING "GPU target device such as Fermi, Kepler, Maxwell, Pascal, Volta, Turing, Ampere")
ExternalProject_Add(MagmaProject
    SOURCE_DIR ${EXT_DIR}/magma-2.6.2
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${EXT_INSTALL_DIR}/magma
        -DGPU_TARGET=${GPU_TARGET}
        -DMAGMA_ENABLE_CUDA=ON
        -DMAGMA_ENABLE_HIP=OFF
        -DBUILD_SHARED_LIBS=ON
        -DUSE_FORTRAN=ON
)
set(MAGMA_LIB_DIR ${EXT_INSTALL_DIR}/magma/lib)
set(MAGMA_INCLUDE_DIR ${EXT_INSTALL_DIR}/magma/include)
# ExternalProject_Add_Step(MagmaProject prebuild
#     COMMAND ${CMAKE_COMMAND} -E echo "BACKEND=cuda\nGPU_TARGET=${GPU_TARGET}" > make.inc
#     COMMAND ${CMAKE_COMMAND} -E make generate
#     COMMENT "Preconfiguring magma with make.inc file -> GPU_TARGET = ${GPU_TARGET}"
#     DEPENDERS configure
#     WORKING_DIRECTORY ${EXT_DIR}/magma
# )
# ExternalProject_Add_StepTargets(MagmaProject prebuild)

# Sundials v6.3.0
set(CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda CACHE PATH "Path where cuda toolkit is installed")
ExternalProject_Add(SundialsProject
    SOURCE_DIR ${EXT_DIR}/sundials
    DEPENDS MagmaProject
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${EXT_INSTALL_DIR}/sundials
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DENABLE_CUDA=ON
        -DCMAKE_CUDA_ARCHITECTURES=80
        -DCUDA_TOOLKIT_ROOT_DIR=${CUDA_TOOLKIT_ROOT_DIR}
        -DEXAMPLES_ENABLE_C=ON
        -DEXAMPLES_ENABLE_CXX=ON
        -DEXAMPLES_ENABLE_CUDA=ON
        -DENABLE_LAPACK=ON
        -DENABLE_MAGMA=ON
        -DMAGMA_DIR=${EXT_INSTALL_DIR}/magma
        -DSUNDIALS_MAGMA_BACKENDS=CUDA
        -DSUNDIALS_INDEX_SIZE=32
        -DENABLE_OPENMP=ON
        -DENABLE_PTHREAD=ON
        -DSUNDIALS_BUILD_WITH_PROFILING=ON
        -DSUNDIALS_LOGGING_LEVEL=4
        -DBUILD_SHARED_LIBS=ON
)
set(SUNDIALS_LIB_DIR ${EXT_INSTALL_DIR}/sundials/lib)
set(SUNDIALS_INCLUDE_DIR ${EXT_INSTALL_DIR}/sundials/include)

# Cantera v2.6.0
set(CANTERA_PREFIX ${EXT_INSTALL_DIR}/cantera)
ExternalProject_Add(CanteraProject
    SOURCE_DIR ${EXT_DIR}/cantera
    DEPENDS SundialsProject
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND  ""
    BUILD_COMMAND scons build prefix=${CANTERA_PREFIX} python_package=none sundials_include=${SUNDIALS_INCLUDE_DIR} sundials_libdir=${SUNDIALS_LIB_DIR}
    INSTALL_COMMAND scons install
)
set(CANTERA_LIB_DIR ${CANTERA_PREFIX}/lib)
set(CANTERA_INCLUDE_DIR ${CANTERA_PREFIX}/include)

########### MAIN PROJECT #########

add_subdirectory(src)
add_subdirectory(app)

############ TESTS ##########
include(CTest)
enable_testing()

add_subdirectory(tests)

########## PACK ##########
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
