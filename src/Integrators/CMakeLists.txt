add_library(Integrators STATIC 
    CanteraIntegrator.cpp
)

target_include_directories(Integrators PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(Integrators PUBLIC ${CANTERA_PREFIX}/include)