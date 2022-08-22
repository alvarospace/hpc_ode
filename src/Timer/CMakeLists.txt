add_library(Timer STATIC 
    Timer.cpp
)

target_include_directories(Timer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)