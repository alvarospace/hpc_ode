add_library(Timer STATIC 
    Timer.cpp
)

target_include_directories(Timer PUBLIC ${CMAKE_SOURCE_DIR}/include)