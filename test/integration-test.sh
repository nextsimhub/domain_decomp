#!/usr/bin/env bash

set -eu

# set filename for each integration test
declare -A FNAMES=(
  ["test_1"]="test_1.nc"
  ["test_2"]="test_2.nc"
  ["test_3"]="test_3.nc"
  ["test_4"]="test_4.nc"
  ["test_1_px"]="test_1.nc"
  ["test_1_py"]="test_1.nc"
  ["test_1_px_py"]="test_1.nc"
  ["test_4_px"]="test_4.nc"
  ["test_4_py"]="test_4.nc"
  ["test_4_px_py"]="test_4.nc"
  ["test_tripolar_2"]="tripolar_grid.nc"
)

# set flags for each integration test
declare -A FLAGS=(
  ["test_1"]="-x x -y y -m mask -o yx"
  ["test_2"]="-x m -y n -m land_mask -o yx"
  ["test_3"]="-x x -y y -m mask -o yx"
  ["test_4"]="-x x -y y -m land_mask -o yx"
  ["test_1_px"]="-x x -y y -m mask -o yx --px"
  ["test_1_py"]="-x x -y y -m mask -o yx --py"
  ["test_1_px_py"]="-x x -y y -m mask -o yx --px --py"
  ["test_4_px"]="-x x -y y -m land_mask -o yx --px"
  ["test_4_py"]="-x x -y y -m land_mask -o yx --py"
  ["test_4_px_py"]="-x x -y y -m land_mask -o yx --px --py"
  ["test_tripolar_2"]="-x x -y y -m mask -o yx --px --py"
)

# run the domain decomp tool for each test case
for TEST in test_1 test_2 test_1_px test_1_py test_1_px_py; do
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

for TEST in test_3 test_4 test_4_px test_4_py test_4_px_py; do
  echo "Running integration corner test '${TEST}'"
  ${MPIEXEC} --oversubscribe ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} \
    ../decomp -g ${FNAMES[${TEST}]} ${FLAGS[${TEST}]} >/dev/null

  for filename in partition_mask_4 partition_metadata_4; do
    ncdump "${filename}.nc" >"${filename}.cdl"
    diff "${filename}.cdl" "${CMAKE_CURRENT_SOURCE_DIR}/${TEST}/ref_${filename}.cdl"
  done
  echo -e "\033[0;32mTest passed\033[0m"
done

for TEST in test_tripolar_2; do
  echo "Running integration tripolar test '${TEST}'"
  ${MPIEXEC} --oversubscribe ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} \
    ../decomp -g ${FNAMES[${TEST}]} ${FLAGS[${TEST}]} >/dev/null

  for filename in partition_mask_2 partition_metadata_2; do
    ncdump "${filename}.nc" >"${filename}.cdl"
    diff "${filename}.cdl" "${CMAKE_CURRENT_SOURCE_DIR}/${TEST}/ref_${filename}.cdl"
  done
  echo -e "\033[0;32mTest passed\033[0m"
done
