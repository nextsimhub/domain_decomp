

add_executable(unit_tests main.cpp)
target_link_libraries(unit_tests PRIVATE doctest::doctest ${INTERNAL_CORE_NAME})
target_include_directories(
  unit_tests
  PRIVATE ${CMAKE_SOURCE_DIR} ${DOCTEST_INCLUDE_DIR})

add_test(NAME unit_tests COMMAND unit_tests)
