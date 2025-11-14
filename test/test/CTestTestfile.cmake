# CMake generated Testfile for 
# Source directory: /home/nvs31/nextsimdg_domain_decomp/domain_decomp/test
# Build directory: /home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_grid_0 "/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/view/bin/mpiexec" "--oversubscribe" "-n" "4" "./test_grid_0")
set_tests_properties(test_grid_0 PROPERTIES  _BACKTRACE_TRIPLES "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;24;add_test;/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;29;create_test;/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;0;")
add_test(test_grid_1 "/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/view/bin/mpiexec" "--oversubscribe" "-n" "4" "./test_grid_1")
set_tests_properties(test_grid_1 PROPERTIES  _BACKTRACE_TRIPLES "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;24;add_test;/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;30;create_test;/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;0;")
add_test(test_grid_2 "/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/view/bin/mpiexec" "--oversubscribe" "-n" "4" "./test_grid_2")
set_tests_properties(test_grid_2 PROPERTIES  _BACKTRACE_TRIPLES "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;24;add_test;/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;31;create_test;/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;0;")
add_test(test_zoltan_0 "/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/view/bin/mpiexec" "--oversubscribe" "-n" "4" "./test_zoltan_0")
set_tests_properties(test_zoltan_0 PROPERTIES  _BACKTRACE_TRIPLES "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;24;add_test;/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;32;create_test;/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;0;")
add_test(test_zoltan_1 "/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/view/bin/mpiexec" "--oversubscribe" "-n" "4" "./test_zoltan_1")
set_tests_properties(test_zoltan_1 PROPERTIES  _BACKTRACE_TRIPLES "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;24;add_test;/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;33;create_test;/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;0;")
add_test(test_zoltan_2 "/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/view/bin/mpiexec" "--oversubscribe" "-n" "4" "./test_zoltan_2")
set_tests_properties(test_zoltan_2 PROPERTIES  _BACKTRACE_TRIPLES "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;24;add_test;/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;34;create_test;/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;0;")
add_test(integration "bash" "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/integration-test.sh")
set_tests_properties(integration PROPERTIES  ENVIRONMENT "MPIEXEC=/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/view/bin/mpiexec;MPIEXEC_NUMPROC_FLAG=-n;MPIEXEC_PREFLAGS=;CMAKE_CURRENT_SOURCE_DIR=/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test" _BACKTRACE_TRIPLES "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;74;add_test;/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeLists.txt;0;")
