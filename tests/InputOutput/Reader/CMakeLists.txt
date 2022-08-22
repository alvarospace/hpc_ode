add_executable(testCsvReader testCsvReader.cpp)
target_link_libraries(testCsvReader PRIVATE Reader)

add_custom_command(TARGET testCsvReader 
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/data/ $<TARGET_FILE_DIR:testCsvReader>/data/
)   

add_test(NAME testCsvReader COMMAND testCsvReader)