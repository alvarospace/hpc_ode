add_library(Mesh OBJECT
    Mesh.cpp
    Point.cpp
    PointBuilder.cpp
)
target_include_directories(Mesh PUBLIC ${CMAKE_SOURCE_DIR}/include)