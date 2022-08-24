add_library(Mesh STATIC 
    Mesh.cpp
    Point.cpp
    PointBuilder.cpp
)

target_include_directories(Mesh PUBLIC ${CMAKE_SOURCE_DIR}/include)