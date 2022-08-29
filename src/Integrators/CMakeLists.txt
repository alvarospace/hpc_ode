add_library(Integrators STATIC 
    CanteraIntegrator.cpp
    CanteraIntegratorOMP.cpp
)

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

add_dependencies(Integrators CanteraProject)

target_link_directories(Integrators PUBLIC ${CANTERA_LIB_DIR})

target_link_libraries(Integrators PUBLIC Mesh Reader cantera fmt yaml-cpp Threads::Threads OpenMP::OpenMP_CXX)

target_include_directories(Integrators PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(Integrators PUBLIC ${CANTERA_INCLUDE_DIR})