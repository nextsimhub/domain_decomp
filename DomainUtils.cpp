/*!
 * @file DomainUtils.cpp
 *
 * @date Oct 07, 2024
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#include "DomainUtils.hpp"
#include <algorithm>
#include <iostream>

int Domain::get_width() const { return p2.x - p1.x; }
int Domain::get_height() const { return p2.y - p1.y; }

int domain_overlap(const Domain d1, const Domain d2, const Edge edge)
{
    int overlap = 0;
    if (edge == TOP || edge == BOTTOM) {
        // check that at least one of the points in domain 1 overlap with domain 2 in the
        // x-direction
        if (d1.p2.x >= d2.p1.x && d1.p1.x <= d2.p2.x) {
            overlap = std::min(d1.p2.x, d2.p2.x) - std::max(d1.p1.x, d2.p1.x);
        }
    } else if (edge == LEFT || edge == RIGHT) {
        // check that at least one of the points in domain 1 overlap with domain 2 in the
        // y-direction
        if (d1.p2.y >= d2.p1.y && d1.p1.y <= d2.p2.y) {
            overlap = std::min(d1.p2.y, d2.p2.y) - std::max(d1.p1.y, d2.p1.y);
        }
    } else {
        std::cerr << "ERROR: edge must be LEFT, RIGHT, BOTTOM, TOP." << std::endl;
        exit(EXIT_FAILURE);
    }
    return overlap;
}
