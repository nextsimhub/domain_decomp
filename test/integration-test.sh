#!/usr/bin/env bash

set -e

# run the domain decomp tool
${MPIEXEC} --oversubscribe ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ../decomp -g test_1.nc -x x -y y -m mask -o 'yx'

# check the two output files are identical to the reference files
for filename in partition_mask_3 partition_metadata_3
do
ncdump "${filename}.nc" > "${filename}.cdl"
diff "${filename}.cdl" "${CMAKE_CURRENT_SOURCE_DIR}/ref_${filename}.cdl"
done

