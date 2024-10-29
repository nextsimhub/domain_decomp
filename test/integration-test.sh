#!/usr/bin/env bash

set -e

# set filename for each integration test
declare -A FNAMES=(
  ["test_1"]="test_1.nc"
  ["test_2"]="test_2.nc"
  ["test_1_p0"]="test_1.nc"
  ["test_1_p1"]="test_1.nc"
  ["test_1_p0_p1"]="test_1.nc"
)

# set flags for each integration test
declare -A FLAGS=(
  ["test_1"]="-x x -y y -m mask -o yx"
  ["test_2"]="-x m -y n -m land_mask -o yx"
  ["test_1_p0"]="-x x -y y -m mask -o yx --p0"
  ["test_1_p1"]="-x x -y y -m mask -o yx --p1"
  ["test_1_p0_p1"]="-x x -y y -m mask -o yx --p0 --p1"
)

# run the domain decomp tool for each test case
for TEST in test_1 test_2 test_1_p0 test_1_p1 test_1_p0_p1; do
  echo "Running integration test '${TEST}'"
  ${MPIEXEC} --oversubscribe ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS} \
    ../decomp -g ${FNAMES[${TEST}]} ${FLAGS[${TEST}]} >/dev/null

  # check the two output files are identical to the reference files
  for filename in partition_mask_3 partition_metadata_3; do
    ncdump "${filename}.nc" >"${filename}.cdl"
    diff "${filename}.cdl" "${CMAKE_CURRENT_SOURCE_DIR}/${TEST}/ref_${filename}.cdl"
  done
  echo -e "\033[0;32mTest passed\033[0m"
done
