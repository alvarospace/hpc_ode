add_library(Reader STATIC 
    csvReader.cpp
)

target_link_libraries(Reader PUBLIC Mesh)
target_include_directories(Reader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)