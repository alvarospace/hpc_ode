# Find Sundials Library
find_library(SUNDIALS_LIBRARY NAMES sundials_cvode 
            PATHS "/home/almousa/install/sundials-6.2/lib")

# Target executable
add_executable(app_serial main.cpp)

# Link libreries
target_link_libraries(app_serial PRIVATE pyjac)
target_link_libraries(app_serial PRIVATE ${SUNDIALS_LIBRARY}) 

# Includes headers
target_include_directories(app_serial PRIVATE "/home/almousa/install/sundials-6.2/include")
target_include_directories(app_serial PRIVATE ${CMAKE_SOURCE_DIR})
