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

TEST_CASE("Domain construction")
{

    SUBCASE("Points in Normalised order")
    {
        const Point p1 { 1, 2 };
        const Point p2 { 3, 4 };
        const Domain d = Domain::fromTwoPoints(p1, p2);
        CHECK(d.p1 == Point { 1, 2 });
        CHECK(d.p2 == Point { 3, 4 });
    }

    SUBCASE("Point 1 right of Point 2")
    {
        const Point p1 { 3, 2 };
        const Point p2 { 1, 4 };
        const Domain d = Domain::fromTwoPoints(p1, p2);
        CHECK(d.p1 == Point { 1, 2 });
        CHECK(d.p2 == Point { 3, 4 });
    }

    SUBCASE("Point 1 above Point 2")
    {
        const Point p1 { 1, 4 };
        const Point p2 { 3, 2 };
        const Domain d = Domain::fromTwoPoints(p1, p2);
        CHECK(d.p1 == Point { 1, 2 });
        CHECK(d.p2 == Point { 3, 4 });
    }

    SUBCASE("Point 1 right and above Point 2")
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

TEST_CASE("Domain status")
{
    SUBCASE("Empty domain")
    {
        const Point p1 { 0, 0 };
        const Point p2 { -1, -1 };
        const Domain d { p1, p2 };

        CHECK(d.isEmpty());
        CHECK(!d.isPoint());
        CHECK(!d.isLine());
        CHECK(!d.isArea());
    }

    SUBCASE("Point domain")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 0, 0 };
        const Domain d { p1, p2 };

        CHECK(!d.isEmpty());
        CHECK(d.isPoint());
        CHECK(!d.isLine());
        CHECK(!d.isArea());
    }

    SUBCASE("Line domain - vertical")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 0, 4 };
        const Domain d { p1, p2 };

        CHECK(!d.isEmpty());
        CHECK(!d.isPoint());
        CHECK(d.isLine());
        CHECK(!d.isArea());
    }

    SUBCASE("Line domain - horizontal")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 4, 0 };
        const Domain d { p1, p2 };

        CHECK(!d.isEmpty());
        CHECK(!d.isPoint());
        CHECK(d.isLine());
        CHECK(!d.isArea());
    }

    SUBCASE("Normal domain")
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

TEST_CASE("Domain - Point symmetry reflection")
{

    SUBCASE("Full domain in positive quadrant")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 2, 1 };
        const Domain d { p1, p2 };

        const Point p { 2, 2 };
        const Domain reflected = pointReflection(p, d);

        CHECK(reflected.p1 == Point { 2, 3 });
        CHECK(reflected.p2 == Point { 4, 4 });
    }

    SUBCASE("Full domain across quadrants")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 2, 1 };
        const Domain d { p1, p2 };

        const Point p { 0, 0 };
        const Domain reflected = pointReflection(p, d);

        CHECK(reflected.p1 == Point { -2, -1 });
        CHECK(reflected.p2 == Point { 0, 0 });
    }

    SUBCASE("Identity map")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 4, 2 };
        const Domain d { p1, p2 };

        const Point p { 2, 1 };
        const Domain reflected = pointReflection(p, d);

        CHECK(reflected.p1 == Point { 0, 0 });
        CHECK(reflected.p2 == Point { 4, 2 });
    }

    SUBCASE("Point domain")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 0, 0 };
        const Domain d { p1, p2 };

        const Point p { 2, 1 };
        const Domain reflected = pointReflection(p, d);

        CHECK(reflected.p1 == Point { 4, 2 });
        CHECK(reflected.p2 == Point { 4, 2 });
    }

    SUBCASE("Line domain - vertical")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 0, 4 };
        const Domain d { p1, p2 };

        const Point p { 1, 1 };
        const Domain reflected = pointReflection(p, d);

        CHECK(reflected.p1 == Point { 2, -2 });
        CHECK(reflected.p2 == Point { 2, 2 });
    }

    SUBCASE("Line domain - horizontal")
    {
        const Point p1 { 0, 0 };
        const Point p2 { 4, 0 };
        const Domain d { p1, p2 };

        const Point p { 3, 0 };
        const Domain reflected = pointReflection(p, d);

        CHECK(reflected.p1 == Point { 2, 0 });
        CHECK(reflected.p2 == Point { 6, 0 });
    }

    SUBCASE("Empty domain")
    {
        const Point p1 { 0, 0 };
        const Point p2 { -1, -1 };
        const Domain d { p1, p2 };

        const Point p { 3, 0 };
        const Domain reflected = pointReflection(p, d);

        CHECK(reflected.isEmpty());
    }
}

namespace {
bool equalDomains(const Domain& d1, const Domain& d2)
{

    // Empty domains are equal to each other
    // and never equal to a non empty domain
    if (d1.isEmpty() || d2.isEmpty()) {
        return d1.isEmpty() && d2.isEmpty();
    }

    // Representation is unique so we can just compare the points
    return d1.p1 == d2.p1 && d1.p2 == d2.p2;
}
}

TEST_CASE("Domain intersection")
{
    SUBCASE("2D intersection")
    {
        const Domain d1 { Point { 0, 0 }, Point { 6, 4 } };
        const Domain d2 { Point { 2, 2 }, Point { 6, 6 } };

        const Domain expectedIntersection { Point { 2, 2 }, Point { 6, 4 } };
        const Domain intersectionDomain = intersection(d1, d2);

        CHECK(equalDomains(intersectionDomain, expectedIntersection));
    }

    SUBCASE("2D disjoint intersection")
    {
        const Domain d1 { Point { 0, 0 }, Point { 6, 4 } };
        const Domain d2 { Point { 7, 2 }, Point { 9, 6 } };

        const Domain intersectionDomain = intersection(d1, d2);

        CHECK(intersectionDomain.isEmpty());
    }

    SUBCASE("Intersection with an empty set")
    {
        const Domain d1 { Point { 0, 0 }, Point { 6, 4 } };
        const Domain d2 { Point { 7, 2 }, Point { 6, 6 } };

        REQUIRE(d2.isEmpty());

        const Domain intersectionDomain = intersection(d1, d2);

        CHECK(intersectionDomain.isEmpty());
    }

    SUBCASE("Intersection with line overlap - vertical")
    {
        const Domain d1 { Point { 0, 0 }, Point { 6, 4 } };
        const Domain d2 { Point { 6, 2 }, Point { 7, 6 } };

        const Domain expectedIntersection { Point { 6, 2 }, Point { 6, 4 } };
        const Domain intersectionDomain = intersection(d1, d2);

        CHECK(intersectionDomain.isLine());
        CHECK(equalDomains(intersectionDomain, expectedIntersection));
    }

    SUBCASE("Intersection with line overlap - horizontal")
    {
        const Domain d1 { Point { 0, 0 }, Point { 5, 2 } };
        const Domain d2 { Point { 1, -2 }, Point { 3, 0 } };

        const Domain expectedIntersection { Point { 1, 0 }, Point { 3, 0 } };
        const Domain intersectionDomain = intersection(d1, d2);

        CHECK(intersectionDomain.isLine());
        CHECK(equalDomains(intersectionDomain, expectedIntersection));
    }

    SUBCASE("Intersection with point overlap")
    {
        const Domain d1 { Point { 0, 0 }, Point { 6, 4 } };
        const Domain d2 { Point { -3, -2 }, Point { 0, 0 } };

        const Domain expectedIntersection { Point { 0, 0 }, Point { 0, 0 } };
        const Domain intersectionDomain = intersection(d1, d2);

        CHECK(intersectionDomain.isPoint());
        CHECK(equalDomains(intersectionDomain, expectedIntersection));
    }
}
