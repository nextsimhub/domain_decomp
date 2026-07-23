/*!
 * @file main.cpp
 *
 * Test driver for the lightweight unit tests 
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "DomainUtils.hpp"


TEST_CASE("Use symbol"){
    // Placeholder test to verify symbol visibility when compiling tests
    Domain d{{0,1}, {2,3}};

    CHECK(d.getWidth() == 2);
    CHECK(d.getHeight() == 2);
}
