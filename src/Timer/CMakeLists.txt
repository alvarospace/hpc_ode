add_library(Timer OBJECT 
    Timer.cpp
)
target_include_directories(Timer PUBLIC ${CMAKE_SOURCE_DIR}/include)