# Domain Decomposition for neXtSIM_DG

## Problem Statement

We address the problem of domain decomposition of numerical ocean and sea-ice models. Such models typically use a sea-land mask to omit unnecessary computations on land. A typical approach is to use a Cartesian (or rectilinear or general block distribution) domain decomposition among processors, ignoring the sea-land mask. This, however, may lead to significant load imbalance which impedes scalability.

We have identyfied the following requirements for the domain decomposition algorithm:
 - produce balanced domains in terms of work and communication
 - produce rectangular domains as this will help in handling communication during halo exchanges between neighbouring processes
 - be static (computed either offline or during the initialiasation phase)
 - be scalable

The proposed approach is based on the Recursive Coordinate Bisection (RCB) geometric partitioining algorithm [^1]. Geometric coordinates are first partitioned into two balanced parts. Partitioning continues recursively in each part until the desired number of balanced parts has been created. The algorithm can be tuned to build rectilinear partitions.

[^1]: M. Berger and S. Bokhari. "A partitioning strategy for nonuniform problems on multiprocessors." IEEE Trans. Computers, C-36 (1987) 570-580.