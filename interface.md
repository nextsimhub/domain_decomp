# Description of the NetCDF-4 partition data format

In typical use `domain_decomp` produces two NetCDF-4 files that describe the
partition of the grid across the processors. These are:

- 'mask' file
- 'metadata' file

The purpose of this document is to give a brief specification of the format of
these files and to identify invariants that can be used for validation of the data.

## Abstract grid

The NextSimDG uses a 2D quadrilateral grid to represent the domain.
As a result, we can represent the grid as a 2D array of entries, where
each entry corresponds to a cell.

Every cell has exactly 4 neighbours (where the outer edge of a domain may be considered
a neighbour).

## Mask file

The mask file is currently not used to pass information to the NextSimDG and is
intended for debugging purposes.

The contents are as described below in the pseudo-CDL:
```
netcdf partition_mask {
dimensions:
    y = ... ; // Grid size in the y direction
    x = ... ; // Grid size in the x direction
variables:
    // Mask that assigns each grid point to a specific rank
    //
    // For: value = pid[i,j]
    //
    // - value in [0; ...) -> value is the rank of the processor in [0; num_processes-1)
    // - value in {-1} -> point is inactive ("Land")
    // - otherwise value is invalid
    //
    int pid(y, x) ;

    // Total number of ranks
    :num_processes = 3;
}
```

## Metadata file

### Assumptions about partition

The metadata NetCDF-4 file is the intended mechanism to communicate the partition
information from `domain_decomp` to the NextSimDG.

The partition will always have the following properties:

- Each rank is assigned a rectangular patch of the domain.
  (i.e., can be represented by integer coordinates of the lower left and upper right
  points)
- Domain edges may be periodic. If it is the case, neighbour information across
  the periodic edge is stored. For non-periodic edges, it is empty.
- Connectivity information is via:
  - 4 edges named (`left`, `right`, `top`, `bottom`)
  - 4 corners named (`top_left`, `top_right`, `bottom_left`, `bottom_right`)


### Dimensions

The NetCDF-4 file contains the following dimensions.
Note that in the description we reference necessary *variables* defined in the
file which will be described later:
```
dimensions:
    NX = ... ; // Number of grid points in the x direction
    NY = ... ; // Number of grid points in the y direction
    P  = ... ; // Number of partitions (ranks)

    // Sum of the number of 'left' neighbours across all partitions
    // Neighbour information is appended across all ranks into a single list
    // `L` is the size of this list.
    //
    // Satisfies: L == sum(left_neighbours)
    //
    // The index for each rank is computed by `cumsum`, e.g. first entry
    // for rank 3, in L-dimensioned variable `var(L)` by:
    //  var[cumsum(left_neighbours)[3]]
    //
    L  = [sum(left_neighbours)] ;

    // Same as `L` but for right edge
    R = [sum(right_neighbours)] ;

    // Same as `L` but for bottom edge
    B = [sum(bottom_neighbours)] ;

    // Same as `L` but for top edge
    T = [sum(top_neighbours)] ;

    // Same as `L` but for the top left corner
    // Since each rank can have at most 1 corner, the following invariant
    // holds: TL <= P
    TL = ...;

    // Other corners are similar to TL
    TR = ...;  // Top right corner
    BR = ...;  // Bottom right corner
    BL = ...;  // Bottom left corner
```

Note that NetCDF does not allow to define ordinary dimensions with size 0.
We do need to transfer empty lists in case, e.g., no partition has a left neighbour
(L == 0). In this case, the dimension will be defined as `UNLIMITED`, but will
never hold any entries.

We use the `UNLIMITED` dimension to represent dimensions with size 0 and empty lists.


### 'bounding_boxes' group

This group contains information about the patch assigned to each rank.

All entries for all variables in this group must be non-negative integers.

Each patch is stored by integer coordinates `(i_x, i_y)` in the grid and the
extent of the patch in the positive direction `(extent_x, extent_y)`.

```
group: bounding_boxes {
  variables:

    // 'i_x' for each rank
  	int domain_x(P) ;

    // 'extent_x' for each rank
  	int domain_extent_x(P) ;

    // 'i_y' for each rank
  	int domain_y(P) ;

    // 'extent_y' for each rank
  	int domain_extent_y(P) ;
}
```

The following invariants should hold:

    - `NX * NY == sum(domain_extent_x * domain_extent_y)` (area consistency)
    - `domain_x + domain_extent_x <= NX` (patch fits in the domain)
    - `domain_y + domain_extent_y <= NY` (patch fits in the domain)


### 'connectivity' group

This group contains information about the topology of the partition.

All entries for all variables in this group must be non-negative integers.

For conciseness, we will only describe the variables for the 'left' edge.
The other edges {'right', 'top', 'bottom' } and corners
{'top_left', 'top_right', 'bottom_left', 'bottom_right'} follow the same format:

The concept of a `perimeter buffer` is important. In our communication model,
each rank establishes two buffers, to be populated with the cell data 'send' and
'recv' (receive). Send is for *other ranks* to read from, recv is for *this rank*
to write to when fetching the data.

The buffers are addressed from the *lower left* in the anti-clockwise direction
(i.e., 'bottom' -> 'right' -> 'top' -> 'left').

When establishing the `recv` buffer, the corners are ignored.


```
group: connectivity {
  variables:

    // Number of left neighbours for each rank
    // i.e., number of ranks reachable through the left edge
    int left_neighbours(P) ;

    // List of left neighbour ids
    // Appended across all ranks
  	int left_neighbour_ids(L) ;

    // For each left neighbour, number of cells through which the neighbour
    // is connected
    // Must satisfy:
    //  sum(left_neighbour_halos) == domain_extent_y[rank]
    //
  	int left_neighbour_halos(L) ;

    // For each left neighbour, coordinate of the first point
    // in the target `send` buffer
    //
  	int left_neighbour_halo_send(L) ;


    // For each left neighbour, coordinate of the first point
    // in the 'own' `recv` buffer
    //
    // This variable is **not** defined for the corners!
    //
  	int left_neighbour_halo_recv(L) ;


    ...
}
```

Of course connectivity is bi-directional. If `n` is connected to `m` through
the left edge, then `m` is connected to `n` through the right edge.

