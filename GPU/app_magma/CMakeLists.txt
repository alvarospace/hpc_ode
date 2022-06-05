# Find OpenMP
#find_package(OpenMP REQUIRED)

# Target executable
add_executable(app_magma main.cu utils.cpp utils.hpp cvode_user.cu cvode_user.cuh)

# Find Sundials
find_package(SUNDIALS COMPONENTS cvode nveccuda sunlinsolmagmadense REQUIRED)

set(SUNDIALS_TARGETS SUNDIALS::cvode SUNDIALS::nveccuda SUNDIALS::sunlinsolmagmadense)

# CUDA OPTIONS
set_property(TARGET app_magma PROPERTY CUDA_SEPARABLE_COMPILATION ON)

# Link libraries
target_link_libraries(app_magma PRIVATE ${SUNDIALS_TARGETS})
target_link_libraries(app_magma PRIVATE pyjacGPU)
#target_link_libraries(app_magma PRIVATE OpenMP::OpenMP_CXX) 

# Include headers
target_include_directories(app_magma PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(app_magma PRIVATE ${SUNDIALS_INCLUDE_DIRS})


