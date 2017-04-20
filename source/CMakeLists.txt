file(GLOB_RECURSE VFXEpoch_SRC "fluids/*.h" "fluids/*.cpp" "utl/*.h" "utl/*.cpp")
add_library(VFXEpoch STATIC ${VFXEpoch_SRC})
