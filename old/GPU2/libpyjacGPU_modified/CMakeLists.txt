# Wrap source code
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.cu ${CMAKE_CURRENT_LIST_DIR}/jacobs/*.cu ${CMAKE_CURRENT_LIST_DIR}/rates/*.cu)
file(GLOB HEADER_FILES ${CMAKE_CURRENT_LIST_DIR}/*.cuh ${CMAKE_CURRENT_LIST_DIR}/jacobs/*.cuh ${CMAKE_CURRENT_LIST_DIR}/rates/*.cuh)


# Create static libpyjacGPU.a
add_library(pyjacGPU_modified ${SRC_FILES} ${HEADER_FILES})

target_include_directories(pyjacGPU_modified PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(pyjacGPU_modified PUBLIC ${CMAKE_CURRENT_LIST_DIR}/jacobs/)
target_include_directories(pyjacGPU_modified PUBLIC ${CMAKE_CURRENT_LIST_DIR}/rates/)

# We need to explicitly state that we need all CUDA files in the 
# particle library to be built with -dc as the member functions 
# could be called by other libraries and executables
set_target_properties(pyjacGPU_modified PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#set_target_properties(pyjacGPU_modified PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_options(pyjacGPU_modified PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xptxas --disable-optimizer-constants>)

