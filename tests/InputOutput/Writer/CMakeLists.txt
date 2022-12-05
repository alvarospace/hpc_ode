add_executable(testCsvWriter testCsvWriter.cpp)
target_link_libraries(testCsvWriter PRIVATE Writer)

add_custom_command(TARGET testCsvWriter 
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/data/ $<TARGET_FILE_DIR:testCsvWriter>/data/
)   

add_test(NAME testCsvWriter COMMAND testCsvWriter)