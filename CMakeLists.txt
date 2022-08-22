cmake_minimum_required(VERSION 3.0.0)
project(HPC_ODE VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(src)



############ TESTS ##########
include(CTest)
enable_testing()

add_subdirectory(tests)


########## PACK ##########
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
