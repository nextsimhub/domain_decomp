# Install script for directory: /home/nvs31/nextsimdg_domain_decomp/domain_decomp

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/build")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdomain_decomp.so.1.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdomain_decomp.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "$ORIGIN/../lib:$ORIGIN/../../lib:/home/nvs31/spack/opt/spack/linux-skylake/netcdf-c-4.9.2-gnookzpsjrpbxdydol2k76nkl3szxnuc/lib:/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/view/lib")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/libdomain_decomp.so.1.0"
    "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/libdomain_decomp.so.1"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdomain_decomp.so.1.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdomain_decomp.so.1"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/nvs31/spack/opt/spack/linux-skylake/netcdf-c-4.9.2-gnookzpsjrpbxdydol2k76nkl3szxnuc/lib:/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/view/lib:::::::::::::::::::::::::::::::::"
           NEW_RPATH "$ORIGIN/../lib:$ORIGIN/../../lib:/home/nvs31/spack/opt/spack/linux-skylake/netcdf-c-4.9.2-gnookzpsjrpbxdydol2k76nkl3szxnuc/lib:/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/view/lib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/libdomain_decomp.so")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/Grid.hpp"
    "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/Partitioner.hpp"
    "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/DomainUtils.hpp"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/nextSIMutilsConfig.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/nextSIMutilsConfig.cmake"
         "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeFiles/Export/c220ae0af1591e9e9e916bba91f25986/nextSIMutilsConfig.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/nextSIMutilsConfig-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/nextSIMutilsConfig.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake" TYPE FILE FILES "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeFiles/Export/c220ae0af1591e9e9e916bba91f25986/nextSIMutilsConfig.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake" TYPE FILE FILES "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeFiles/Export/c220ae0af1591e9e9e916bba91f25986/nextSIMutilsConfig-noconfig.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/decomp" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/decomp")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/decomp"
         RPATH "$ORIGIN/../lib:$ORIGIN/../../lib:/home/nvs31/spack/opt/spack/linux-skylake/netcdf-c-4.9.2-gnookzpsjrpbxdydol2k76nkl3szxnuc/lib:/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/._view/kud5ek6g6cfs4qty3a3hqjumqtf26zd5/lib:/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/view/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/decomp")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/decomp" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/decomp")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/decomp"
         OLD_RPATH "/home/nvs31/spack/opt/spack/linux-skylake/netcdf-c-4.9.2-gnookzpsjrpbxdydol2k76nkl3szxnuc/lib:/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/._view/kud5ek6g6cfs4qty3a3hqjumqtf26zd5/lib:/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/view/lib:::::::::::::::::::::::::::::::::"
         NEW_RPATH "$ORIGIN/../lib:$ORIGIN/../../lib:/home/nvs31/spack/opt/spack/linux-skylake/netcdf-c-4.9.2-gnookzpsjrpbxdydol2k76nkl3szxnuc/lib:/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/._view/kud5ek6g6cfs4qty3a3hqjumqtf26zd5/lib:/home/nvs31/spack/var/spack/environments/sasip_env_test/.spack-env/view/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/decomp")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/CMakeFiles/decomp.dir/install-cxx-module-bmi-noconfig.cmake" OPTIONAL)
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/test/cmake_install.cmake")
  include("/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/examples/cmake_install.cmake")

endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/nvs31/nextsimdg_domain_decomp/domain_decomp/test/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
