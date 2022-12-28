project(ODEIntegrator LANGUAGES CUDA CXX)

set(CMAKE_CUDA_ARCHITECTURES 80)

add_subdirectory(Logger)
add_subdirectory(Timer)
add_subdirectory(Mesh)
add_subdirectory(InputOutput/Reader)
add_subdirectory(InputOutput/Writer)
add_subdirectory(Mechanism)
add_subdirectory(Integrators)
