add_library(Writer OBJECT
    csvWriter.cpp
)
target_link_libraries(Writer PUBLIC Mesh Logger)
target_include_directories(Writer PUBLIC ${CMAKE_SOURCE_DIR}/include)