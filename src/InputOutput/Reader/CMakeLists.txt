add_library(Reader OBJECT
    csvReader.cpp
)
target_link_libraries(Reader PUBLIC Mesh Logger)
target_include_directories(Reader PUBLIC ${CMAKE_SOURCE_DIR}/include)