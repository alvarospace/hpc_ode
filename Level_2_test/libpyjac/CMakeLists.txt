# Wrap source code
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.c ${CMAKE_CURRENT_LIST_DIR}/jacobs/*.c)
file(GLOB HEADER_FILES ${CMAKE_CURRENT_LIST_DIR}/*.h ${CMAKE_CURRENT_LIST_DIR}/jacobs/*.h)

# Create static libpyjac.a
add_library(pyjac ${SRC_FILES} ${HEADER_FILES})

target_include_directories(pyjac PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(pyjac PUBLIC ${CMAKE_CURRENT_LIST_DIR}/jacobs/)
