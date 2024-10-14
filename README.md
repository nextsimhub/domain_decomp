# Domain Decomposition for neXtSIM_DG

## Coding conventions

We follow the same [criteria](https://github.com/nextsimhub/nextsimdg?tab=readme-ov-file#coding-conventions) as those used in the `neXtSIM_DG` project.

A [dedicated clang format file](https://github.com/nextsimhub/nextsimdg/blob/main/.clang-format) has been designed for the code. You may run it locally and manually with the command ```clang-format -i $yourfile``` or have a plugin with your favorite code editor or implement a git pre-commit hook locally by putting this pre-commit file in your .git/hooks/. An example pre-commit file can be found at [.pre-commit](https://github.com/nextsimhub/nextsimdg/blob/658_pre-commit-date/.pre-commit). This clang formatting will also be run each time a pull request is done as part of the continuous integration.

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
 -DENABLE_PNETCDF=ON \
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
It is recommended to build the code in a separate directory from the source directory. The basic steps for building with CMake are:
1. Create a build directory, outside of the source directory.
2. In your build directory run `cmake <path-to-src>`
3. It is recommended to set options by calling `ccmake` in your build directory. Alternatively you can use the `-DVARIABLE=value` syntax in the previous step.
4. Run `make` to build.
5. Run `make test` to run tests.
6. Run `make install` to install.

The project installs a shared library that can be imported by other CMake projects as `neXtSIMutils::domain_decomp,` as well as a binary `decomp` that be used to partition grids directly.

### How to run

Running the help of the `domain` tool gives the following:

```
$ ./build/decomp -h
Usage: ./build/decomp [options]
Options:
  -h [ --help ]             Display this help message
  -g [ --grid ] arg         NetCDF grid file
  -x [ --xdim ] arg (=x)    Name of x dimension in netCDF grid file
  -y [ --ydim ] arg (=y)    Name of y dimension in netCDF grid file
  -o [ --order ] arg (=yx)  Order of dimensions in netCDF grid file, e.g., 'yx'
                            or 'xy'
  -m [ --mask ] arg (=mask) Mask variable name in netCDF grid file
  -i [ --ignore-mask ]      Ignore mask in netCDF grid file

```

We can see that the most of the options have defaults (show in parentheses) e.g., if you do not specify the name of the land
mask `-m/--mask`, it will default to `"mask"`. To explain the other options it is instructive to first look at an example grid
file.

```
$ ncdump test_2.nc
netcdf test_2 {
dimensions:
        m = 6 ;
        n = 4 ;
variables:
        int land_mask(n, m) ;

// global attributes:
                :title = "Non-default dimension naming and ordering" ;
data:

 land_mask =
  0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1,
  0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1 ;
}
```

This is from a grid file `test_2.nc` where the `x` and `y` dimensions are called `m` and `n` respectively. Looking at the
`variables` sections we can see that the land mask `land_mask` is stored in the order `land_mask(n, m)`, where `n` is the first
dimension. No matter what the dimensions are called inside of the grid file, neXtSIM_DG will expect the metadata files (produced
by `decomp`) to have a specific output format. By default, the ordering of `decomp` assumes `'yx'` which (due to the names of our
dimensions) corresponds to `n` being the first dimension and `m` being the second. This can be changed to `'xy'` for example, if
the land mask was stored transposed i.e., `land_mask(m, n)`.

To run `decomp` on the `test_2.nc` file we can use the following command:
```
$ mpirun -n 2 ../decomp -g test_2.nc -x m -y n -o 'yx' -m "land_mask"
```

This should produce the following output (stdout) as well as two files, `partition_mask_2.nc` and `partition_metadata_2.nc`.

```
$ mpirun -n 2 ../decomp -g test_2.nc -x m -y n -o 'yx' -m "land_mask"
Partitioning total time: 9.2725e-05 (secs)
Partitioning Statistics:
 Total weight of dots = 12
 Weight on each proc: ave = 6, max = 6, min = 6
 Maximum weight of single dot = 1
 Send count: ave = 0, min = 0, max = 0
 Recv count: ave = 0, min = 0, max = 0
 Max dots: ave = 6, min = 6, max = 6
 Max memory: ave = 8, min = 8, max = 8
 # of OverAlloc: ave = 0, min = 0, max = 0
 Median find iteration counts:
   Serial iterations per process: avg 0.000000 min 0 max 0
   Parallel iterations:
     Per process: avg 1.000000 min 1 max 1
     Total for all parallel cuts: 1
     Detail of first cuts:
       Level 1 cut count: avg 1.000000 variance 0.000000 actual: 1
 Start-up time (secs): ave = 2.70065e-05, min = 2.6947e-05, max = 2.7066e-05
 Start-up time (%): ave = 29.1254, min = 29.0612, max = 29.1895
 Pre-median time (secs): ave = 4.311e-06, min = 4.034e-06, max = 4.588e-06
 Pre-median time (%): ave = 4.64923, min = 4.3505, max = 4.94796
 Median time (secs): ave = 8.2655e-06, min = 8.066e-06, max = 8.465e-06
 Median time (%): ave = 8.91399, min = 8.69884, max = 9.12915
 Comm time (secs): ave = 5.0859e-05, min = 5.0418e-05, max = 5.13e-05
 Comm time (%): ave = 54.8493, min = 54.3737, max = 55.3249
 STATS Runs 1  bal  CURRENT 1.000000  MAX 1.000000  MIN 1.000000  AVG 1.000000
 STATS Runs 1  moveCnt CURRENT 0.000000  MAX 0.000000  MIN 0.000000  AVG 0.000000
 STATS DETAIL count:  min 6  max 6  avg 6.000000  imbal 1.000000
 STATS DETAIL wdim = 0; no detail available
```

The `decomp` tool produces two netCDF-4 files (using the classic data model) named `partition_mask_<num_mpi_processes>.nc` and `partition_metadata_<num_mpi_processes>.nc` with the following layout:

- `partition_mask_<num_mpi_processes>.nc` - can be used to check the partitioning
- `partition_metadata_<num_mpi_processes>.nc` - is used by neXtSIM_DG to read in the domain decomposition information, such as
  domain sizes, neighbours, halos etc.

The following is an example of the partition mask file generated by running the command above.

```
netcdf partition_mask_2 {
dimensions:
y = 4 ;
x = 6 ;
variables:
int pid(y, x) ;
// global attributes:
:num_processes = 2 ;
data:
 pid =
  -1, -1, -1, -1, -1, -1,
  0, 0, 0, 1, 1, 1,
  -1, -1, -1, -1, -1, -1,
  0, 0, 0, 1, 1, 1 ;
}
```

The netCDF variable `pid` is defined as the process ID of each point in the grid. `pid=-1` notes regions where the land mask is
zero. Given an example partition mask (`pid`) layout below, the "domain" layout would be as follows:

```
        partition_mask_2                         domain layout
                                             ┌──────────────────┐
O┌───►X -1, -1, -1, -1, -1, -1               │         2        │
 │       0,  0,  0,  1,  1,  1               ├─────────┬────────┤
 │      -1, -1, -1, -1, -1, -1     Y▲        │         │        │
Y▼       0,  0,  0,  1,  1,  1      │        │    0    │    1   │
        -1, -1, -1, -1, -1, -1      │        │         │        │
         2,  2,  2,  2,  2,  2     0└───►X   └─────────┴────────┘
```

The following is an example of the partition metadata file generated by running the command above.

```
netcdf partition_metadata_2 {
dimensions:
NX = 6 ;
NY = 4 ;
P = 2 ;
L = 1 ;
R = 1 ;
B = UNLIMITED ; // (0 currently)
T = UNLIMITED ; // (0 currently)
group: bounding_boxes {
  variables:
  int domain_x(P) ;
  int domain_extent_x(P) ;
  int domain_y(P) ;
  int domain_extent_y(P) ;
  data:
   domain_x = 0, 3 ;
   domain_extent_x = 3, 3 ;
   domain_y = 0, 0 ;
   domain_extent_y = 4, 4 ;
  } // group bounding_boxes
group: connectivity {
  variables:
  int left_neighbours(P) ;
  int left_neighbour_ids(L) ;
  int left_neighbour_halos(L) ;
  int right_neighbours(P) ;
  int right_neighbour_ids(R) ;
  int right_neighbour_halos(R) ;
  int bottom_neighbours(P) ;
  int bottom_neighbour_ids(B) ;
  int bottom_neighbour_halos(B) ;
  int top_neighbours(P) ;
  int top_neighbour_ids(T) ;
  int top_neighbour_halos(T) ;
  data:
   left_neighbours = 0, 1 ;
   left_neighbour_ids = 0 ;
   left_neighbour_halos = 4 ;
   right_neighbours = 1, 0 ;
   right_neighbour_ids = 1 ;
   right_neighbour_halos = 4 ;
   bottom_neighbours = 0, 0 ;
   top_neighbours = 0, 0 ;
  } // group connectivity
}
```

The netCDF variables `domain_x/y` are defined as the coordinates of the bottom left corner of the bounding box for each MPI
process using zero-based indexing and `domain_extent_x/y` are the extents in the corresponding dimensions of the bounding box
for each MPI process. The file also defines the variables `X_neighbours(P)`, `X_neighbour_ids(X_dim)` and
`X_neighbour_halos(X_dim)`, where `X` is `top/bottom/left/right`, which correspond to the number of neighbours per process, the
neighbour IDs and halo sizes of each process sorted from lower to higher MPI rank.
