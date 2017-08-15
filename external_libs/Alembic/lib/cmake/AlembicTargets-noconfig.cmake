#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Alembic::Alembic" for configuration ""
set_property(TARGET Alembic::Alembic APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(Alembic::Alembic PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_NOCONFIG "/usr/local/lib/libImath-2_2.so;/usr/local/lib/libIlmThread-2_2.so;/usr/local/lib/libIex-2_2.so;/usr/local/lib/libIexMath-2_2.so;/usr/local/lib/libHalf.so;-lpthread;-lm"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libAlembic.so.1.7.3"
  IMPORTED_SONAME_NOCONFIG "libAlembic.so.1.7"
  )

list(APPEND _IMPORT_CHECK_TARGETS Alembic::Alembic )
list(APPEND _IMPORT_CHECK_FILES_FOR_Alembic::Alembic "${_IMPORT_PREFIX}/lib/libAlembic.so.1.7.3" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
