/*!
 * @file DomainUtils.hpp
 *
 * @date Oct 07, 2024
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#ifndef DOMAINUTILS_HPP
#define DOMAINUTILS_HPP

/*!
 * @brief 2D point structs. Points can be used to construct Domains.
 */
struct Point {

    int x, y;
};

/*!
 * @brief After grid decomposition, we are left with 2D domains.
 *
 * Each Domain is described by 2 points as follows
 *          p2 = (x2, y2)
 *   ┌────────────x
 *   │            │
 *   │            │
 *   │            │
 *   x────────────┘
 * P1 = (x1, y1)
 */
struct Domain {
    Point p1, p2;
};

/*!
 * @brief Compute the overlap between two domains. Will return zero if domains
 * do not overlap.
 *
 * Note that this doesn't check if two domains are neighbours, just that they
 * overlap in the x or y direction.
 *
 * @param d1 First Domain
 * @param d2 Second Domain
 * @param dir direction to find overlap ('x' or 'y')
 */
int domainOverlap(const Domain d1, const Domain d2, const char dir);

#endif /* DOMAINUTILS_HPP */