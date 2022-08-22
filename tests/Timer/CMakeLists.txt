add_executable(testTimer testTimer.cpp)
target_link_libraries(testTimer PRIVATE Timer)

add_test(NAME testTimer COMMAND testTimer)