cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(nextsim_domain_decomp CXX)

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(netCDF REQUIRED)
find_package(Zoltan REQUIRED)

# Generic options for compilation
add_compile_options(-Wall -Wextra -pedantic -fPIC)

# Executable
add_executable(partition grid.cpp partition.cpp)
target_link_libraries(partition PRIVATE MPI::MPI_CXX netcdf netcdf-cxx4 ${Zoltan_LIBRARIES})
target_link_directories(partition PRIVATE "${netCDF_LIB_DIR}")
target_include_directories(partition PRIVATE "${netCDF_INCLUDE_DIR}" "${Zoltan_INCLUDE_DIRS}")
