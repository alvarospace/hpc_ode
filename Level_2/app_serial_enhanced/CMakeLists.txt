# Find Sundials Library
find_library(SUNDIALS_LIBRARY NAMES sundials_cvode 
            PATHS "/home/almousa/install/sundials-6.2/lib")


# Find OpenMP
find_package(OpenMP REQUIRED)

# Target executable
add_executable(app_serial_enhanced main.cpp utils.cpp utils.hpp cvode_user.cpp cvode_user.hpp)

# Link libraries
target_link_libraries(app_serial_enhanced PRIVATE pyjac)
target_link_libraries(app_serial_enhanced PRIVATE ${SUNDIALS_LIBRARY})
target_link_libraries(app_serial_enhanced PRIVATE OpenMP::OpenMP_CXX) 

# Includes headers
target_include_directories(app_serial_enhanced PRIVATE "/home/almousa/install/sundials-6.2/include")
target_include_directories(app_serial_enhanced PRIVATE ${CMAKE_SOURCE_DIR})
