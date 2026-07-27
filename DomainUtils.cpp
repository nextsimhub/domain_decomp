/*!
 * @file DomainUtils.cpp
 *
 * @date Oct 07, 2024
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#include "DomainUtils.hpp"
#include <algorithm>
#include <iostream>
#include <tuple>

int Domain::getWidth() const { return p2.x - p1.x; }
int Domain::getHeight() const { return p2.y - p1.y; }

int domainOverlap(const Domain d1, const Domain d2, const Edge edge)
{
    int overlap = 0;
    if (edge == TOP || edge == BOTTOM) {
        overlap = std::min(d1.p2.x, d2.p2.x) - std::max(d1.p1.x, d2.p1.x);
    } else if (edge == LEFT || edge == RIGHT) {
        overlap = std::min(d1.p2.y, d2.p2.y) - std::max(d1.p1.y, d2.p1.y);
    } else {
        std::cerr << "ERROR: edge must be LEFT, RIGHT, BOTTOM, TOP." << std::endl;
        exit(EXIT_FAILURE);
    }
    return overlap;
}

Domain Domain::fromTwoPoints(const Point& p1, const Point& p2)
{
    Domain d;
    // Lower left corner
    d.p1.x = std::min(p1.x, p2.x);
    d.p1.y = std::min(p1.y, p2.y);

    // Upper right corner
    d.p2.x = std::max(p1.x, p2.x);
    d.p2.y = std::max(p1.y, p2.y);
    return d;
}

Domain pointReflection(Point c, Domain d)
{
    if (d.isEmpty()) {
        return d;
    }

    const auto p1_reflected = 2 * c - d.p1;
    const auto p2_reflected = 2 * c - d.p2;
    return Domain::fromTwoPoints(p1_reflected, p2_reflected);
}

namespace {

/**
 * @brief Compute the intersection of two 1D intervals
 *
 * Note that the interval is represented as [start, end]:
 *   * start < end is a valid line segment
 *   * start == end is a valid point
 *   * start > end is an empty interval
 *
 * @param a First interval
 * @param b Second interval
 * @return The intersection of the two intervals.
 */
std::pair<int, int> intervalIntersection(const std::pair<int, int>& a, const std::pair<int, int>& b)
{
    const int start = std::max(a.first, b.first);
    const int end = std::min(a.second, b.second);
    return { start, end };
}
}

Domain intersection(Domain d1, Domain d2)
{
    // Intersection in x-axis
    // Make a structured binding in C++17
    int x_start, x_end;
    std::tie(x_start, x_end) = intervalIntersection({ d1.p1.x, d1.p2.x }, { d2.p1.x, d2.p2.x });

    // Intersection in y-axis
    // Make a structured binding in C++17
    int y_start, y_end;
    std::tie(y_start, y_end) = intervalIntersection({ d1.p1.y, d1.p2.y }, { d2.p1.y, d2.p2.y });

    return { Point { x_start, y_start }, Point { x_end, y_end } };
}
