FILE(
  GLOB EXAMPLE_SMOKE
  "${CMAKE_CURRENT_SOURCE_DIR}/Common/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/Common/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/smoke/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/smoke/*.cpp"
)

FILE(
  GLOB_RECURSE EXAMPLE_UNIT_TEST
  "${CMAKE_CURRENT_SOURCE_DIR}/Common/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/Common/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/examiner/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/examiner/*.cpp"
)

FILE(
  GLOB_RECURSE 2DPARTICLE_VISUALIZER
  "${CMAKE_CURRENT_SOURCE_DIR}/Common/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/Common/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/viz2D/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/viz2D/*.cpp"
)

FILE(
  GLOB_RECURSE VORTEX_RINGS_2D
  "${CMAKE_CURRENT_SOURCE_DIR}/Common/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/Common/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/vortex_rings_2d/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/vortex_rings_2d/*.cpp"
)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

#Homemade public headers and libs
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/examples/Common)
link_directories(${CMAKE_SOURCE_DIR}/vfxepoch_libs)
link_directories(${CMAKE_SOURCE_DIR}/external_libs/OpenEXR/libs)

message("WILL INSTALL INCLUDE PATH = ${CMAKE_SOURCE_DIR}/include")
message("WILL INSTALL LIB PATH = ${CMAKE_SOURCE_DIR}/vfxepoch_libs")

add_executable(smoke ${EXAMPLE_SMOKE})
add_executable(examiner ${EXAMPLE_UNIT_TEST})
add_executable(2D_particle_visualizer ${2DPARTICLE_VISUALIZER})
add_executable(vortex_rings_2d ${VORTEX_RINGS_2D})

target_link_libraries(smoke VFXEpoch ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
target_link_libraries(examiner VFXEpoch ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
target_link_libraries(2D_particle_visualizer VFXEpoch ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} IlmImf-2_2 Half)
target_link_libraries(vortex_rings_2d VFXEpoch ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
