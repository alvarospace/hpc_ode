
# Find OpenMP
find_package(OpenMP REQUIRED)

# Target executable
add_executable(app_openmp main.cpp utils.cpp utils.hpp cvode_user.cpp cvode_user.hpp)

# Find Sundials
find_package(SUNDIALS COMPONENTS sundials_cvode)

if(SUNDIALS_FOUND)
    target_link_libraries(app_openmp PRIVATE ${SUNDIALS_LIBRARIES})
    target_include_directories(app_openmp PRIVATE ${SUNDIALS_INCLUDE_DIRS})

elseif(SUNDIALS_EXTERN_ENABLE)
    target_link_libraries(app_openmp PRIVATE ${SUNDIALS_EXTERN_LIBS}/libsundials_cvode.a)
    target_include_directories(app_openmp PRIVATE ${SUNDIALS_EXTERN_INCLUDE_PATH})

endif()


# Link libraries
target_link_libraries(app_openmp PRIVATE pyjac)
target_link_libraries(app_openmp PRIVATE OpenMP::OpenMP_CXX) 

# Includes headers
target_include_directories(app_openmp PRIVATE ${CMAKE_SOURCE_DIR})