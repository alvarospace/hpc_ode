add_library(Logger STATIC
    ThreadSafeQueue.cpp
    BaseLogger.cpp
)

find_package(Threads REQUIRED)

target_include_directories(Logger PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(Logger PRIVATE Threads::Threads)