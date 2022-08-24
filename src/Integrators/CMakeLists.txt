add_library(Integrators STATIC 
    CanteraIntegrator.cpp
)

target_include_directories(Integrators PUBLIC ${CMAKE_SOURCE_DIR}/include)