
# Target executable
add_executable(app_serial main.cpp)

# Find Sundials
find_package(SUNDIALS COMPONENTS sundials_cvode)

if(SUNDIALS_FOUND)
    target_link_libraries(app_serial PRIVATE ${SUNDIALS_LIBRARIES})
    target_include_directories(app_serial PRIVATE ${SUNDIALS_INCLUDE_DIRS})

elseif(SUNDIALS_EXTERN_ENABLE)
    target_link_libraries(app_serial PRIVATE ${SUNDIALS_EXTERN_LIBS}/libsundials_cvode.a)
    target_include_directories(app_serial PRIVATE ${SUNDIALS_EXTERN_INCLUDE_PATH})

endif()


# Link libraries
target_link_libraries(app_serial PRIVATE pyjac)

# Includes headers
target_include_directories(app_serial PRIVATE ${CMAKE_SOURCE_DIR})
