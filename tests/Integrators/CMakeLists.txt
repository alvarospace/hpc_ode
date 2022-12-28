add_executable(testCanteraIntegrator testCanteraIntegrator.cpp)

target_link_libraries(testCanteraIntegrator PRIVATE Mesh Logger Timer Reader IntegratorsCPU)
add_custom_command(TARGET testCanteraIntegrator
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/data/ $<TARGET_FILE_DIR:testCanteraIntegrator>/data/
)

add_test(NAME testCanteraIntegrator COMMAND testCanteraIntegrator)

##########

add_executable(testCVodeIntegrator testCVodeIntegrator.cpp)
target_link_libraries(testCVodeIntegrator PRIVATE Mesh Logger Timer Reader IntegratorsCPU)
add_test(NAME testCVodeIntegrator COMMAND testCVodeIntegrator)

#########

add_executable(testCVodeIntegratorGPU testCVodeIntegratorGPU.cpp)
target_link_libraries(testCVodeIntegratorGPU PRIVATE Mesh Logger Timer Reader Writer IntegratorsGPU IntegratorsCPU)
add_test(NAME testCVodeIntegratorGPU COMMAND testCVodeIntegratorGPU)
