add_library(Reader STATIC 
    csvReader.cpp
)

target_link_libraries(Reader PUBLIC Logger Mesh)
target_include_directories(Reader PUBLIC ${CMAKE_SOURCE_DIR}/include)