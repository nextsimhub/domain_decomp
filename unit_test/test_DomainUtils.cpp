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
}
