#!/usr/bin/env bash

set -e

# run the domain decomp tool on test_1.nc
${MPIEXEC} --oversubscribe ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ../decomp -g test_1.nc -x x -y y -m mask -o 'yx' > tmp.log

# check the two output files are identical to the reference files
for filename in partition_mask_3 partition_metadata_3
do
ncdump "${filename}.nc" > "${filename}.cdl"
diff "${filename}.cdl" "${CMAKE_CURRENT_SOURCE_DIR}/test_1/ref_${filename}.cdl"
done

# run the domain decomp tool on test_2.nc
${MPIEXEC} --oversubscribe ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} ../decomp -g test_2.nc -x m -y n -m land_mask -o 'yx' > tmp.log

# check the two output files are identical to the reference files
for filename in partition_mask_3 partition_metadata_3
do
ncdump "${filename}.nc" > "${filename}.cdl"
diff "${filename}.cdl" "${CMAKE_CURRENT_SOURCE_DIR}/test_2/ref_${filename}.cdl"
done
