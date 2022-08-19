# Domain Decomposition for neXtSIM_DG

## Problem Statement

We address the problem of domain decomposition of numerical ocean and sea-ice models. Such models typically use a sea-land mask to omit unnecessary computations on land. A typical approach is to use a cartesian (or rectilinear or general block distribution) domain decomposition among processors, where the domain is divided in equally sized sub-domains, ignoring the sea-land mask. This, however, may lead to significant load imbalance which impedes scalability.

We have identified the following requirements for the domain decomposition algorithm:
 - produce balanced sub-domains in terms of work and communication
 - produce rectangular sub-domains as this will help in handling communication during halo exchanges between neighbouring processes
 - be static (computed either offline or during the initialiasation phase)
 - be scalable

The proposed approach is based on the Recursive Coordinate Bisection (RCB) geometric partitioning algorithm [^1]. Geometric coordinates are first partitioned into two balanced parts. Partitioning continues recursively in each part until the desired number of balanced parts has been created. The algorithm can be tuned to build rectilinear partitions. We are using the implementation of the RCB algorithm available in the [Zoltan](https://sandialabs.github.io/Zoltan/) library developed by Sandia National Laboratories.

[^1]: M. Berger and S. Bokhari. "A partitioning strategy for nonuniform problems on multiprocessors." IEEE Trans. Computers, C-36 (1987) 570-580.

## Getting Started

### Requirements
* A Unix-like operating system (e.g., Linux or Mac OS X)
* ANSI C++ compiler
* MPI library for message passing (e.g., MPICH or OpenMPI)
* CMake >= 3.10
* [netCDF-4 C](https://github.com/Unidata/netcdf-c/releases/tag/v4.8.1), built with parallel I/O support to netCDF-4 files through HDF5 and to classic files through PnetCDF
* Zoltan, built with CMake from the **[Trilinos](https://github.com/trilinos/Trilinos.git)** package
* [Catch2](https://github.com/catchorg/Catch2) for unit testing
* [Boost](https://www.boost.org/) program_options library

#### Building Zoltan from source

```
git clone https://github.com/trilinos/Trilinos.git
cd Trilinos
mkdir build && cd build
cmake \
 -DTPL_ENABLE_MPI:BOOL=ON \
 -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
 -DTrilinos_ENABLE_Zoltan:BOOL=ON \
 -DTrilinos_ENABLE_Fortran:BOOL=OFF \
 -DZoltan_ENABLE_EXAMPLES:BOOL=ON \
 -DZoltan_ENABLE_TESTS:BOOL=ON \
 -DBUILD_SHARED_LIBS=ON \
 -DCMAKE_INSTALL_PREFIX:FILEPATH=<path> \
 ..
make
make test
make install
```

If you want to use an MPI library installed in a non-standard location you will need to additionally set:
```
 -DCMAKE_CXX_COMPILER=<mpicxx> -DCMAKE_C_COMPILER=<mpicc>
```

#### Building PnetCDF from source
Download v1.12.3 of [PnetCDF](https://github.com/Parallel-NetCDF/PnetCDF/releases/tag/checkpoint.1.12.3).
```
cd PnetCDF
autoreconf -i
CXX=mpicxx CC=mpicc ./configure --disable-fortran --enable-shared --prefix=<installation_path>
make
make install
```

#### Building netCDF from source with parallel I/O
Download v4.8.1 of [netCDF-4 C](https://github.com/Unidata/netcdf-c/releases/tag/v4.8.1). Parallel I/O is implemented through HDF5 (built with parallel I/O).

For netCDF-C (for more details see [here](https://docs.unidata.ucar.edu/netcdf-c/current/netCDF-CMake.html) for requirements):
```
cd netcdf-c-4.8.1
mkdir build && cd build
cmake \
 -DCMAKE_C_COMPILER=mpicc \
 -DCMAKE_CXX_COMPILER=mpicxx \
 -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_PREFIX_PATH=<path_to_hdf5> \
 -DCMAKE_INSTALL_PREFIX=<installation_path> \
 ..
make
make test
make install
```

### Installation
It is recommended to build the code in a separate directory form the source directory. The basic steps for building with CMake are:
1. Create a build directory, outside of the source directory.
2. In your build directory run `cmake <path-to-src>`
3. It is recommended to set options by calling `ccmake` in your build directory. Alternatively you can use the `-DVARIABLE=value` syntax in the previous step.
4. Run `make` to build.
5. Run `make test` to run tests.
6. Run `make install` to install.

The project installs a shared library that can be imported by other CMake projects as neXtSIMutils::domain_decomp, as well as a binary `decomp` that be used to partition grids directly.

### How to run
The binary `decomp` can be used to partition a 2D grid with an optional land mask represented as a netCDF file. By default, the name of the dimensions in the netCDF file are `x` and `y` with the `y` dimension increasing the fastest, and the name of the variable representing the land mask is `mask`. For netCDF files using the enhanced model, we assume the group that contains the information of interest is named `data`. These can be overridden using command-line options. For example:
```
mpirun -n 2 ./decomp grid.nc --dim0 y --dim1 x --mask land_mask
```

The `decomp` tool produces two netCDF files named `partition_mask_<num_mpi_processes>.nc` and `partition_metadata_<num_mpi_processes>.nc` with the following layout:

```
netcdf partition_mask_2 {
dimensions:
	x = 30 ;
	y = 30 ;
variables:
	short pid(x, y) ;

// global attributes:
		:num_processes = 2s ;
data:

 pid = ...
}
```

The netCDF variable `pid` is defined as the process ID of each point in the grid.

```
netcdf partition_metadata_2 {
dimensions:
	P = 2 ;
variables:
	int global_x(P) ;
	int global_y(P) ;
	int local_extent_x(P) ;
	int local_extent_y(P) ;
data:

 global_x = 0, 16 ;
 global_y = 0, 0 ;
 local_extent_x = 16, 14 ;
 local_extent_y = 30, 30 ;
}

```

The netCDF variables `global_x/y` are defined as the coordinates of the upper left corner of the bounding box for each MPI process using zero-based indexing and `local_extent_x/y` are the extents in the corresponding dimensions of the bounding box for each MPI process.