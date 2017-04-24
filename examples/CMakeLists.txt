FILE(
  GLOB_RECURSE EXAMPLE_SMOKE
  "${CMAKE_CURRENT_SOURCE_DIR}/smoke/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/smoke/*.cpp"
)

FILE(
  GLOB_RECURSE EXAMPLE_UNIT_TEST
  "${CMAKE_CURRENT_SOURCE_DIR}/general_test/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/general_test/*.cpp"
)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

#Homemade public headers and libs
include_directories(${CMAKE_SOURCE_DIR}/include)
link_directories(${CMAKE_SOURCE_DIR}/vfxepoch_libs)

message("WILL INSTALL INCLUDE PATH = ${CMAKE_SOURCE_DIR}/include")
message("WILL INSTALL LIB PATH = ${CMAKE_SOURCE_DIR}/vfxepoch_libs")

add_executable(smoke ${EXAMPLE_SMOKE})
add_executable(general_test ${EXAMPLE_UNIT_TEST})

target_link_libraries(smoke VFXEpoch ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
target_link_libraries(general_test VFXEpoch ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
