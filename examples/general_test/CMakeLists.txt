File(GLOB GENERAL_TEST_SRC
 "*.h"
 "*.cpp"
 )

add_executable(general_test ${GENERAL_TEST_SRC})

target_link_libraries (general_test LINK_PUBLIC VFXEpoch)
