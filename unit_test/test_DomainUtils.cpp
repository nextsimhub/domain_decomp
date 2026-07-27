#include "DomainUtils.hpp"

#include "doctest/doctest.h"

TEST_CASE("Point equality")
{
    const Point p1 { 1, 2 };
    const Point p2 { 1, 2 };
    const Point p3 { 3, 4 };

    CHECK(p1 == p2);
    CHECK(p1 != p3);
}

TEST_CASE("Point arithmetic")
{
    // This is not the most sophisticated test
    // It's purpose is to verify that the operators can be called and there is
    // nothing silly going on (e.g. addition or x,y components beeing flipped)
    // Single inputs should be enough to cover that
    //
    const int p_x = 1;
    const int p_y = 2;

    const Point p1 { p_x, p_y };

    SUBCASE("Point + Point")
    {
        const int x = 3;
        const int y = 4;

        const Point p2 { x, y };
        const Point p3 = p1 + p2;
        CHECK(p3 == Point { p_x + x, p_y + y });
    }

    SUBCASE("Point - Point")
    {
        const int x = 3;
        const int y = 4;

        const Point p2 { x, y };
        const Point p3 = p1 - p2;
        CHECK(p3 == Point { p_x - x, p_y - y });
    }

    SUBCASE("Point * scalar")
    {
        const int scalar = 3;
        const Point p2 = p1 * scalar;
        CHECK(p2 == Point { p_x * scalar, p_y * scalar });
    }

    SUBCASE("scalar * Point")
    {
        const int scalar = 3;
        const Point p2 = scalar * p1;
        CHECK(p2 == Point { p_x * scalar, p_y * scalar });
    }

    SUBCASE("Negation")
    {
        const Point p2 = -p1;
        CHECK(p2 == Point { -p_x, -p_y });
    }
}

TEST_SUITE("Domain construction")
{

    TEST_CASE("Points in Normalised order")
    {
        const Point p1 { 1, 2 };
        const Point p2 { 3, 4 };
        const Domain d = Domain::fromTwoPoints(p1, p2);
        CHECK(d.p1 == Point { 1, 2 });
        CHECK(d.p2 == Point { 3, 4 });
    }

    TEST_CASE("Point 1 right of Point 2")
    {
        const Point p1 { 3, 2 };
        const Point p2 { 1, 4 };
        const Domain d = Domain::fromTwoPoints(p1, p2);
        CHECK(d.p1 == Point { 1, 2 });
        CHECK(d.p2 == Point { 3, 4 });
    }

    TEST_CASE("Point 1 above Point 2")
    {
        const Point p1 { 1, 4 };
        const Point p2 { 3, 2 };
        const Domain d = Domain::fromTwoPoints(p1, p2);
        CHECK(d.p1 == Point { 1, 2 });
        CHECK(d.p2 == Point { 3, 4 });
    }

    TEST_CASE("Point 1 right and above Point 2")
    {
        const Point p1 { 3, 4 };
        const Point p2 { 1, 2 };
        const Domain d = Domain::fromTwoPoints(p1, p2);
        CHECK(d.p1 == Point { 1, 2 });
        CHECK(d.p2 == Point { 3, 4 });
    }
}

TEST_CASE("Domain Getters")
{
    const Point p1 { 1, 2 };
    const Point p2 { 4, 4 };
    const Domain d { p1, p2 };

    CHECK(d.getWidth() == 3);
    CHECK(d.getHeight() == 2);
}

TEST_SUITE("Domain status")
{

    TEST_CASE("Empty domain")
    {
        const Point p1 { 0, 0 };
        const Point p2 { -1, -1 };
        const Domain d { p1, p2 };

        CHECK(d.isEmpty());
        CHECK(!d.isPoint());
        CHECK(!d.isLine());
        CHECK(!d.isArea());
    }

    TEST_CASE("Point domain")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 0, 0 };
        const Domain d { p1, p2 };

        CHECK(!d.isEmpty());
        CHECK(d.isPoint());
        CHECK(!d.isLine());
        CHECK(!d.isArea());
    }

    TEST_CASE("Line domain - vertical")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 0, 4 };
        const Domain d { p1, p2 };

        CHECK(!d.isEmpty());
        CHECK(!d.isPoint());
        CHECK(d.isLine());
        CHECK(!d.isArea());
    }

    TEST_CASE("Line domain - horizontal")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 4, 0 };
        const Domain d { p1, p2 };

        CHECK(!d.isEmpty());
        CHECK(!d.isPoint());
        CHECK(d.isLine());
        CHECK(!d.isArea());
    }

    TEST_CASE("Normal domain")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 4, 4 };
        const Domain d { p1, p2 };

        CHECK(!d.isEmpty());
        CHECK(!d.isPoint());
        CHECK(!d.isLine());
        CHECK(d.isArea());
    }
}

TEST_SUITE("Domain - Point symmetry reflection")
{

    TEST_CASE("Full domain in positive quadrant")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 2, 1 };
        const Domain d { p1, p2 };

        const Point p { 2, 2 };
        const Domain reflected = pointReflection(p, d);

        CHECK(reflected.p1 == Point { 2, 3 });
        CHECK(reflected.p2 == Point { 4, 4 });
    }

    TEST_CASE("Full domain across quadrants")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 2, 1 };
        const Domain d { p1, p2 };

        const Point p { 0, 0 };
        const Domain reflected = pointReflection(p, d);

        CHECK(reflected.p1 == Point { -2, -1 });
        CHECK(reflected.p2 == Point { 0, 0 });
    }

    TEST_CASE("Identity map")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 4, 2 };
        const Domain d { p1, p2 };

        const Point p { 2, 1 };
        const Domain reflected = pointReflection(p, d);

        CHECK(reflected.p1 == Point { 0, 0 });
        CHECK(reflected.p2 == Point { 4, 2 });
    }

    TEST_CASE("Point domain")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 0, 0 };
        const Domain d { p1, p2 };

        const Point p { 2, 1 };
        const Domain reflected = pointReflection(p, d);

        CHECK(reflected.p1 == Point { 4, 2 });
        CHECK(reflected.p2 == Point { 4, 2 });
    }

    TEST_CASE("Line domain - vertical")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 0, 4 };
        const Domain d { p1, p2 };

        const Point p { 1, 1 };
        const Domain reflected = pointReflection(p, d);

        CHECK(reflected.p1 == Point { 2, -2 });
        CHECK(reflected.p2 == Point { 2, 2 });
    }

    TEST_CASE("Line domain - horizontal")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 4, 0 };
        const Domain d { p1, p2 };

        const Point p { 3, 0 };
        const Domain reflected = pointReflection(p, d);

        CHECK(reflected.p1 == Point { 2, 0 });
        CHECK(reflected.p2 == Point { 6, 0 });
    }

    TEST_CASE("Empty domain")
    {
        const Point p1 { 0, 0 };
        const Point p2 { -1, -1 };
        const Domain d { p1, p2 };

        const Point p { 3, 0 };
        const Domain reflected = pointReflection(p, d);

        CHECK(reflected.isEmpty());
    }
}

TEST_SUITE("Domain equality")
{
    TEST_CASE("Normal domains")
    {
        const Domain d1 { Point { 0, 0 }, Point { 4, 4 } };
        const Domain d2 { Point { 0, 0 }, Point { 4, 4 } };
        const Domain d3 { Point { 1, 1 }, Point { 4, 4 } };

        CHECK(d1 == d2);
        CHECK(d1 != d3);
    }

    TEST_CASE("Empty domains")
    {
        const Domain empty { Point { 0, 0 }, Point { -1, -1 } };
        const Domain empty2 { Point { 0, 0 }, Point { -2, -2 } };
        const Domain d2 { Point { 1, 1 }, Point { 2, 2 } };

        CHECK(empty == empty2);
        CHECK(empty != d2);
    }
}

TEST_SUITE("Domain intersection")
{

    TEST_CASE("2D intersection")
    {
        const Domain d1 { Point { 0, 0 }, Point { 6, 4 } };
        const Domain d2 { Point { 2, 2 }, Point { 6, 6 } };

        const Domain expectedIntersection { Point { 2, 2 }, Point { 6, 4 } };
        const Domain intersectionDomain = intersection(d1, d2);

        CHECK(intersectionDomain == expectedIntersection);
    }

    TEST_CASE("2D disjoint intersection")
    {
        const Domain d1 { Point { 0, 0 }, Point { 6, 4 } };
        const Domain d2 { Point { 7, 2 }, Point { 9, 6 } };

        const Domain intersectionDomain = intersection(d1, d2);

        CHECK(intersectionDomain.isEmpty());
    }

    TEST_CASE("Intersection with an empty set")
    {
        const Domain d1 { Point { 0, 0 }, Point { 6, 4 } };
        const Domain d2 { Point { 7, 2 }, Point { 6, 6 } };

        REQUIRE(d2.isEmpty());

        const Domain intersectionDomain = intersection(d1, d2);

        CHECK(intersectionDomain.isEmpty());
    }

    TEST_CASE("Intersection with line overlap - vertical")
    {
        const Domain d1 { Point { 0, 0 }, Point { 6, 4 } };
        const Domain d2 { Point { 6, 2 }, Point { 7, 6 } };

        const Domain expectedIntersection { Point { 6, 2 }, Point { 6, 4 } };
        const Domain intersectionDomain = intersection(d1, d2);

        CHECK(intersectionDomain.isLine());
        CHECK(intersectionDomain == expectedIntersection);
    }

    TEST_CASE("Intersection with line overlap - horizontal")
    {
        const Domain d1 { Point { 0, 0 }, Point { 5, 2 } };
        const Domain d2 { Point { 1, -2 }, Point { 3, 0 } };

        const Domain expectedIntersection { Point { 1, 0 }, Point { 3, 0 } };
        const Domain intersectionDomain = intersection(d1, d2);

        CHECK(intersectionDomain.isLine());
        CHECK(intersectionDomain == expectedIntersection);
    }

    TEST_CASE("Intersection with point overlap")
    {
        const Domain d1 { Point { 0, 0 }, Point { 6, 4 } };
        const Domain d2 { Point { -3, -2 }, Point { 0, 0 } };

        const Domain expectedIntersection { Point { 0, 0 }, Point { 0, 0 } };
        const Domain intersectionDomain = intersection(d1, d2);

        CHECK(intersectionDomain.isPoint());
        CHECK(intersectionDomain == expectedIntersection);
    }
}

TEST_SUITE("Domain Connections")
{
    TEST_CASE("Valid connections - Edges")
    {
        const Domain d1 { Point {
                              0,
                              0,
                          },
            Point { 4, 4 } };

        SUBCASE("Bottom edge")
        {
            const Domain shared { Point { 1, 0 }, Point { 3, 0 } };
            const Connection connection = identifyNeighbourConnectionType(d1, shared);
            CHECK(connection.getType() == ConnectionType::EDGE);
            CHECK(connection.getEdge() == BOTTOM);
        }

        SUBCASE("Top edge")
        {
            const Domain shared { Point { 1, 4 }, Point { 3, 4 } };
            const Connection connection = identifyNeighbourConnectionType(d1, shared);
            CHECK(connection.getType() == ConnectionType::EDGE);
            CHECK(connection.getEdge() == TOP);
        }

        SUBCASE("Left edge")
        {
            const Domain shared { Point { 0, 1 }, Point { 0, 3 } };
            const Connection connection = identifyNeighbourConnectionType(d1, shared);
            CHECK(connection.getType() == ConnectionType::EDGE);
            CHECK(connection.getEdge() == LEFT);
        }

        SUBCASE("Right edge")
        {
            const Domain shared { Point { 4, 1 }, Point { 4, 3 } };
            const Connection connection = identifyNeighbourConnectionType(d1, shared);
            CHECK(connection.getType() == ConnectionType::EDGE);
            CHECK(connection.getEdge() == RIGHT);
        }
    }

    TEST_CASE("Valid connections - Corners")
    {
        const Domain d1 { Point {
                              0,
                              0,
                          },
            Point { 4, 4 } };

        SUBCASE("Bottom left corner")
        {
            const Domain shared { Point { 0, 0 }, Point { 0, 0 } };
            const Connection connection = identifyNeighbourConnectionType(d1, shared);
            CHECK(connection.getType() == ConnectionType::CORNER);
            CHECK(connection.getCorner() == BOTTOM_LEFT);
        }

        SUBCASE("Bottom right corner")
        {
            const Domain shared { Point { 4, 0 }, Point { 4, 0 } };
            const Connection connection = identifyNeighbourConnectionType(d1, shared);
            CHECK(connection.getType() == ConnectionType::CORNER);
            CHECK(connection.getCorner() == BOTTOM_RIGHT);
        }

        SUBCASE("Top left corner")
        {
            const Domain shared { Point { 0, 4 }, Point { 0, 4 } };
            const Connection connection = identifyNeighbourConnectionType(d1, shared);
            CHECK(connection.getType() == ConnectionType::CORNER);
            CHECK(connection.getCorner() == TOP_LEFT);
        }

        SUBCASE("Top right corner")
        {
            const Domain shared { Point { 4, 4 }, Point { 4, 4 } };
            const Connection connection = identifyNeighbourConnectionType(d1, shared);
            CHECK(connection.getType() == ConnectionType::CORNER);
            CHECK(connection.getCorner() == TOP_RIGHT);
        }
    }

    TEST_CASE("Points on edges that are not corners")
    {
        // This should not appear. We identify it as a no connection
        const Domain d1 { Point {
                              0,
                              0,
                          },
            Point { 4, 4 } };

        SUBCASE("Bottom edge")
        {
            const Domain shared { Point { 1, 0 }, Point { 1, 0 } };
            const Connection connection = identifyNeighbourConnectionType(d1, shared);
            CHECK(connection.getType() == ConnectionType::INVALID);
        }
        SUBCASE("Top edge")
        {
            const Domain shared { Point { 1, 4 }, Point { 1, 4 } };
            const Connection connection = identifyNeighbourConnectionType(d1, shared);
            CHECK(connection.getType() == ConnectionType::INVALID);
        }
        SUBCASE("Left edge")
        {
            const Domain shared { Point { 0, 1 }, Point { 0, 1 } };
            const Connection connection = identifyNeighbourConnectionType(d1, shared);
            CHECK(connection.getType() == ConnectionType::INVALID);
        }
        SUBCASE("Right edge")
        {
            const Domain shared { Point { 4, 1 }, Point { 4, 1 } };
            const Connection connection = identifyNeighbourConnectionType(d1, shared);
            CHECK(connection.getType() == ConnectionType::INVALID);
        }
    }

    TEST_CASE("Shared edge extending outside domain")
    {
        const Domain d1 { Point {
                              0,
                              0,
                          },
            Point { 4, 4 } };

        // We test only bottom edge
        const Domain shared { Point { 1, 0 }, Point { 5, 0 } };
        const Connection connection = identifyNeighbourConnectionType(d1, shared);
        CHECK(connection.getType() == ConnectionType::INVALID);
    }

    TEST_CASE("Shared domain is an area")
    {
        const Domain d1 { Point {
                              0,
                              0,
                          },
            Point { 4, 4 } };
        const Domain shared { Point { 1, 1 }, Point { 3, 3 } };
        const Connection connection = identifyNeighbourConnectionType(d1, shared);
        CHECK(connection.getType() == ConnectionType::INVALID);
    }

    TEST_CASE("Domain is not an area")
    {
        const Domain shared { Point { 0, 0 }, Point { 1, 1 } };

        SUBCASE("Empty")
        {
            const Domain d1 { Point {
                                  5,
                                  5,
                              },
                Point { 4, 4 } };
            CHECK_THROWS(identifyNeighbourConnectionType(d1, shared));
        }
        SUBCASE("Point")
        {
            const Domain d1 { Point {
                                  0,
                                  0,
                              },
                Point { 0, 0 } };
            CHECK_THROWS(identifyNeighbourConnectionType(d1, shared));
        }
        SUBCASE("Line")
        {
            const Domain d1 { Point {
                                  0,
                                  0,
                              },
                Point { 0, 4 } };
            CHECK_THROWS(identifyNeighbourConnectionType(d1, shared));
        }
    }

    TEST_CASE("Shared domain is empty")
    {
        const Domain d1 { Point {
                              0,
                              0,
                          },
            Point { 4, 4 } };
        const Domain shared { Point { 5, 5 }, Point { 4, 4 } };

        const Connection connection = identifyNeighbourConnectionType(d1, shared);
        CHECK(connection.getType() == ConnectionType::INVALID);
    }
}
