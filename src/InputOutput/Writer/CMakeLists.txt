add_library(Writer STATIC 
    csvWriter.cpp
)

target_link_libraries(Writer PUBLIC Logger Mesh)
target_include_directories(Writer PUBLIC ${CMAKE_SOURCE_DIR}/include)
