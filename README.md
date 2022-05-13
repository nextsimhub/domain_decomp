# Domain Decomposition for neXtSIM_DG

## Problem Statement

We address the problem of domain decomposition of numerical ocean and sea-ice models. Such models typically use a sea-land mask to omit unnecessary computations on land. A typical approach is to use a cartesian (or rectilinear or general block distribution) domain decomposition among processors, where the domain is divided in equally sized sub-domains, ignoring the sea-land mask. This, however, may lead to significant load imbalance which impedes scalability.

We have identyfied the following requirements for the domain decomposition algorithm:
 - produce balanced sub-domains in terms of work and communication
 - produce rectangular sub-domains as this will help in handling communication during halo exchanges between neighbouring processes
 - be static (computed either offline or during the initialiasation phase)
 - be scalable

The proposed approach is based on the Recursive Coordinate Bisection (RCB) geometric partitioning algorithm [^1]. Geometric coordinates are first partitioned into two balanced parts. Partitioning continues recursively in each part until the desired number of balanced parts has been created. The algorithm can be tuned to build rectilinear partitions.

[^1]: M. Berger and S. Bokhari. "A partitioning strategy for nonuniform problems on multiprocessors." IEEE Trans. Computers, C-36 (1987) 570-580.

### Example: partitioning the NEMO ORCA1 global domain

![NEMO ORCA1 decomposition to 16 processes](img/ORCA1_RCB_16.png)

## Getting Started

### Requirements
* CMake >= 3.10
* C++ compiler
* MPI
* netCDF-4 C, netCDF-4 C++, built with parallel I/O support
* Zoltan, built with CMake from the **[Trilinos](https://github.com/trilinos/Trilinos)** package

### Installation
It is recommended to build the code in a separate directory form the source directory. The basic steps for building with CMake are:
1. Create a build directory, outside of the source directory.
2. In your build directory run `cmake <path-to-src>`
3. It is recommended to set options by calling `ccmake` in your build directory. Alternatively you can use the `-DVARIABLE=value` syntax in the previous step.
4. Run `make` to build.
