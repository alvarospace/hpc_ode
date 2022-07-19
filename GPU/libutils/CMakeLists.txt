# Create static libutils.a
add_library(utils utils.cpp utils.hpp)

target_include_directories(utils PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(utils PRIVATE ${CMAKE_SOURCE_DIR})
