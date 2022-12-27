# CPU
add_library(IntegratorsCPU STATIC 
    CanteraIntegrator.cpp
    CanteraIntegratorOMP.cpp
    CVodeIntegrator.cpp
    CVodeIntegratorOMP.cpp
    Drivers/CVodeSerialDriver.cpp
)

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

add_dependencies(IntegratorsCPU CanteraProject SundialsProject)

target_link_directories(IntegratorsCPU PUBLIC ${CANTERA_LIB_DIR})
target_link_directories(IntegratorsCPU PUBLIC ${SUNDIALS_LIB_DIR})

target_link_libraries(IntegratorsCPU PUBLIC Logger Mesh Reader mechanismCPU cantera fmt yaml-cpp Threads::Threads OpenMP::OpenMP_CXX sundials_cvode.a)

target_include_directories(IntegratorsCPU PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(IntegratorsCPU PUBLIC ${CANTERA_INCLUDE_DIR})
target_include_directories(IntegratorsCPU PUBLIC ${SUNDIALS_INCLUDE_DIR})
target_include_directories(IntegratorsCPU PUBLIC ${CMAKE_SOURCE_DIR}/src)


# GPU
add_library(IntegratorsGPU STATIC
    CVodeIntegratorGPU.cpp
    Drivers/CVodeGPUDriver.cu
)

add_dependencies(IntegratorsGPU SundialsProject MagmaProject)

target_link_directories(IntegratorsGPU PUBLIC ${SUNDIALS_LIB_DIR})

target_link_libraries(IntegratorsGPU PUBLIC Mesh Reader mechanismGPU sundials_cvode.a sundials_nveccuda.a sundials_sunlinsolmagmadense.so sundials_sunmatrixmagmadense.so)

target_include_directories(IntegratorsGPU PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(IntegratorsGPU PUBLIC ${SUNDIALS_INCLUDE_DIR})
target_include_directories(IntegratorsGPU PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_include_directories(IntegratorsGPU PUBLIC ${MAGMA_INCLUDE_DIR})

set_target_properties(IntegratorsGPU PROPERTIES CUDA_SEPARABLE_COMPILATION ON)