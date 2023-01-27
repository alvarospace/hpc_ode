find_package(OpenMP REQUIRED)

add_executable(testLogger testLogger.cpp)
target_link_libraries(testLogger PRIVATE Logger OpenMP::OpenMP_CXX)

add_test(NAME testLogger COMMAND testLogger)