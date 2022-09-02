add_executable(testCanteraIntegrator testCanteraIntegrator.cpp)

target_link_libraries(testCanteraIntegrator PRIVATE Integrators Timer)
add_custom_command(TARGET testCanteraIntegrator
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/data/ $<TARGET_FILE_DIR:testCanteraIntegrator>/data/
)

add_test(NAME testCanteraIntegrator COMMAND testCanteraIntegrator)

##########

add_executable(testCVodeIntegrator testCVodeIntegrator.cpp)
target_link_libraries(testCVodeIntegrator PRIVATE Integrators Timer)
add_test(NAME testCVodeIntegrator COMMAND testCVodeIntegrator)
