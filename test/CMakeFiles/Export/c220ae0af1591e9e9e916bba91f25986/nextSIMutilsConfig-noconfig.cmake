#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "nextSIMutils::domain_decomp" for configuration ""
set_property(TARGET nextSIMutils::domain_decomp APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(nextSIMutils::domain_decomp PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libdomain_decomp.so.1.0"
  IMPORTED_SONAME_NOCONFIG "libdomain_decomp.so.1"
  )

list(APPEND _cmake_import_check_targets nextSIMutils::domain_decomp )
list(APPEND _cmake_import_check_files_for_nextSIMutils::domain_decomp "${_IMPORT_PREFIX}/lib/libdomain_decomp.so.1.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
