file(GLOB_RECURSE VFXEpoch_SRC "fluids/*.h" "fluids/*.cpp" "utl/*.h" "utl/*.cpp")
add_library(VFXEpoch STATIC ${VFXEpoch_SRC})
file(GLOB_RECURSE VFXEpoch_HEADERS "fluids/*.h" "utl/*.h")
INSTALL(FILES ${VFXEpoch_HEADERS} DESTINATION ../include)
