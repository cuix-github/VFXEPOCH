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

FILE(
  GLOB_RECURSE WRITE_ALEMBIC_POINTS
  "${CMAKE_CURRENT_SOURCE_DIR}/Common/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/Common/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/write_alembic/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/write_alembic/*.cpp"
)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

#Homemade public headers and libs
include_directories(/usr/local/include)
include_directories(/usr/local/include/OpenEXR)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/examples/Common)
link_directories(${CMAKE_SOURCE_DIR}/vfxepoch_libs)
link_directories(${CMAKE_SOURCE_DIR}/external_libs/Alembic/libs)
link_directories(/usr/local/lib)

message("-- INSTALL INCLUDE PATH = ${CMAKE_SOURCE_DIR}/include")
message("-- INSTALL LIB PATH = ${CMAKE_SOURCE_DIR}/vfxepoch_libs")

add_executable(smoke ${EXAMPLE_SMOKE})
add_executable(examiner ${EXAMPLE_UNIT_TEST})
add_executable(2D_particle_visualizer ${2DPARTICLE_VISUALIZER})
add_executable(vortex_rings_2d ${VORTEX_RINGS_2D})
add_executable(write_alembic_points ${WRITE_ALEMBIC_POINTS})

target_link_libraries(smoke VFXEpoch ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
target_link_libraries(examiner VFXEpoch ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
target_link_libraries(2D_particle_visualizer VFXEpoch ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} IlmImf-2_2 Half)
target_link_libraries(vortex_rings_2d VFXEpoch ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} IlmImf-2_2 Imath-2_2 Half Alembic)
target_link_libraries(write_alembic_points VFXEpoch ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} IlmImf-2_2 Imath-2_2 Half Alembic)
