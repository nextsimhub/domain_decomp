cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Hide symbols by default
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)

set(PROJECT_NAME neXtSIMutils)
set(LIB_NAME domain_decomp)
set(EXEC_NAME decomp)
set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
project(${PROJECT_NAME} VERSION ${VERSION_MAJOR}.${VERSION_MINOR} LANGUAGES CXX)

# Set RPATH behaviour
set(CMAKE_SKIP_RPATH FALSE)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
# Embed absolute paths to external libraries that are not part of
# the project, (they are expected to be at the same location on all
# machines the project will be deployed to
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
include(GNUInstallDirs)
# Define RPATH for executables via a relative expression to enable a
# fully relocatable package
file(RELATIVE_PATH relDir1
 ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
 ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
file(RELATIVE_PATH relDir2
 ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/examples
 ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_RPATH $ORIGIN/${relDir1}:$ORIGIN/${relDir2})

# Dependencies
find_package(Catch2 REQUIRED)
find_package(MPI REQUIRED COMPONENTS CXX)
find_package(netCDF REQUIRED)
find_package(Zoltan REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)

# Generic options for compilation
add_compile_options(-Wall -Wextra -pedantic -fPIC)

# Library
add_library(${LIB_NAME} SHARED
  ${PROJECT_SOURCE_DIR}/Grid.cpp
  ${PROJECT_SOURCE_DIR}/Partitioner.cpp
  ${PROJECT_SOURCE_DIR}/ZoltanPartitioner.cpp
  )
set_target_properties(${LIB_NAME} PROPERTIES
  VERSION ${VERSION_MAJOR}.${VERSION_MINOR}
  SOVERSION ${VERSION_MAJOR}
  )
# Create an export header to control symbol visibility
include(GenerateExportHeader)
generate_export_header(${LIB_NAME}
  EXPORT_FILE_NAME ${LIB_NAME}_export.hpp
  EXPORT_MACRO_NAME LIB_EXPORT
  )
# Create an alias to be used within the project and by other projects
add_library(${PROJECT_NAME}::${LIB_NAME} ALIAS ${LIB_NAME})
# Set library dependencies
target_link_libraries(${LIB_NAME}
  PUBLIC MPI::MPI_CXX
  PRIVATE netcdf ${Zoltan_LIBRARIES}
  )
target_link_directories(${LIB_NAME}
  PRIVATE ${netCDF_LIB_DIR}
  )
target_include_directories(${LIB_NAME}
  PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE ${netCDF_INCLUDE_DIR} ${Zoltan_INCLUDE_DIRS}
  )
# Set the header files that need to be installed
set_target_properties(${LIB_NAME} PROPERTIES
  PUBLIC_HEADER "Grid.hpp;Partitioner.hpp"
  )
# Create import targets
install(TARGETS ${LIB_NAME}
  EXPORT ${PROJECT_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
install(EXPORT ${PROJECT_NAME}
  FILE ${PROJECT_NAME}Config.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake
  )

# Executable
add_executable(${EXEC_NAME}
  ${PROJECT_SOURCE_DIR}/Grid.cpp
  ${PROJECT_SOURCE_DIR}/Partitioner.cpp
  ${PROJECT_SOURCE_DIR}/ZoltanPartitioner.cpp
  ${PROJECT_SOURCE_DIR}/main.cpp)
# Create an alias to be used within the project and by other projects
add_executable(${PROJECT_NAME}::${EXEC_NAME} ALIAS ${EXEC_NAME})
# Set executable dependencies
target_include_directories(${EXEC_NAME}
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${netCDF_INCLUDE_DIR} ${Zoltan_INCLUDE_DIRS}
  )
target_link_libraries(${EXEC_NAME}
  PUBLIC MPI::MPI_CXX
  PRIVATE netcdf ${Zoltan_LIBRARIES} Boost::program_options
  )
target_link_directories(${EXEC_NAME}
  PRIVATE ${netCDF_LIB_DIR}
  )
install(TARGETS ${EXEC_NAME} DESTINATION ${CMAKE_INSTALL_BINDIR})

# Unit tests
enable_testing()
add_subdirectory(test)

# Examples
add_subdirectory(examples)
