project(ODEIntegrator LANGUAGES CUDA CXX)

set(CMAKE_CUDA_ARCHITECTURES 80)

add_subdirectory(Logger)
add_subdirectory(Timer)
add_subdirectory(Mesh)
add_subdirectory(InputOutput/Reader)
add_subdirectory(InputOutput/Writer)
add_subdirectory(Mechanism)
add_subdirectory(Integrators)

add_library(${PROJECT_NAME} STATIC)
target_link_libraries(${PROJECT_NAME} PUBLIC
    Logger
    Timer
    Mesh
    Reader
    Writer
    IntegratorsCPU
    IntegratorsGPU
)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
set_target_properties(${PROJECT_NAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)