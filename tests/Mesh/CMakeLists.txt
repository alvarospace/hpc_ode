add_executable(testPoint testPoint.cpp)
target_link_libraries(testPoint PRIVATE Mesh)

add_executable(testMesh testMesh.cpp)
target_link_libraries(testMesh PRIVATE Mesh)

add_test(NAME testPoint COMMAND testPoint)
add_test(NAME testMesh  COMMAND testMesh)