/*!
 * @file DomainUtils.hpp
 *
 * @date Oct 07, 2024
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#ifndef DOMAINUTILS_HPP
#define DOMAINUTILS_HPP

#include <array>
#include <stdexcept>

enum Edge { LEFT, RIGHT, BOTTOM, TOP, N_EDGE };
static constexpr std::array<Edge, N_EDGE> edges = { LEFT, RIGHT, BOTTOM, TOP };

enum Corner { TOP_LEFT, TOP_RIGHT, BOTTOM_RIGHT, BOTTOM_LEFT, N_CORNER };
static constexpr std::array<Corner, N_CORNER> corners
    = { TOP_LEFT, TOP_RIGHT, BOTTOM_RIGHT, BOTTOM_LEFT };

/*!
 * @brief 2D point structs. Points can be used to construct Domains.
 */
struct Point {

    int x, y;

    // Define arithmetic on Points
    Point operator+(const Point& other) const { return { x + other.x, y + other.y }; }
    Point operator-(const Point& other) const { return { x - other.x, y - other.y }; }
    Point operator*(const int scalar) const { return { x * scalar, y * scalar }; }
    Point operator-() const { return { -x, -y }; }

    // Implement equality
    bool operator==(const Point& other) const { return x == other.x && y == other.y; }
    bool operator!=(const Point& other) const { return !(*this == other); }
};
inline Point operator*(const int scalar, const Point& p) { return p * scalar; }

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
    /*!
     * @brief return width of domain
     */
    int getWidth() const;
    /*!
     * @brief return height of domain
     */
    int getHeight() const;

    /**
     * @brief Create a domain from two points
     *
     * Unlike the constructor, the points do not need to be in particular order.
     * i.e p1.y > p2.y is OK.
     *
     * Two points uniquely define a patch, the function will return a 'normalised'
     * representation of the domain where the invariant:
     *  `p1.x <= p2.x and p1.y <= p2.y` holds
     *
     * Degenerate domains (Lines and points) are allowed
     *
     * @param p1 First point
     * @param p2 Second point
     * @return Domain defined by the two points
     */
    static Domain fromTwoPoints(const Point& p1, const Point& p2);

    /**
     * @brief Non-normalised domain represents an empty domain (empty set)
     */
    bool isEmpty() const { return p1.x > p2.x || p1.y > p2.y; }

    /**
     *  @brief We support degenerate domains that collapse to a point
     */
    bool isPoint() const { return p1.x == p2.x && p1.y == p2.y; }

    /**
     * @brief We support degenerate domains that collapse to a line
     */
    bool isLine() const
    {
        const bool isVerticalLine = p1.x == p2.x && p1.y != p2.y;
        const bool isHorizontalLine = p1.x != p2.x && p1.y == p2.y;

        return isVerticalLine || isHorizontalLine;
    }
    /**
     * @brief For completeness we check is a domain represents an area (patch)
     */
    bool isArea() const { return !isEmpty() && !isPoint() && !isLine(); }

    /**
     * @brief Define domain equality
     *
     * Empty domains are always equal to each other
     */
    bool operator==(const Domain& other) const
    {
        // Empty domains are equal to each other
        // and never equal to a non empty domain
        if (this->isEmpty() || other.isEmpty()) {
            return this->isEmpty() && other.isEmpty();
        }

        // Representation is unique so we can just compare the points
        return this->p1 == other.p1 && this->p2 == other.p2;
    }
    bool operator!=(const Domain& other) const { return !(*this == other); }
};

/*!
 * @brief Compute the overlap between two domains. Will return zero if domains do not overlap.
 *
 * Note that this doesn't check if two domains are neighbours, just that they overlap in the x or y
 * direction.
 *
 * @param d1 First Domain
 * @param d2 Second Domain
 * @param dir direction to find overlap ('x' or 'y')
 */
int domainOverlap(const Domain d1, const Domain d2, const Edge edge);

/*!
 * @brief Applies point symmetry to a domain
 *
 * For a degenerate domain, returns a degenerate domain.
 * Makes no guarantees about the points inside.
 * Reflection of an empty set is an empty set is our logic here.
 *
 * @param p Point of the symmetry
 * @param d Domain to be reflected
 * @return Domain reflected about point p
 */
Domain pointReflection(Point p, Domain d);

/*!
 *
 * @brief Domain intersection
 *
 * Takes two domains and returns the intersection of the two.
 * The domains include their boundaries, so two domains that share
 * a section of an edge or a point will return a degenerate "line"
 * or "point" domain.
 *
 * The intersection of disjoint domains returns a degenerate domain that
 * represents an empty set.
 *
 * The intention is to use this function to compute neighbourhood relation
 * between the domains.
 *
 * The operation should be commutative and associative, but we don't test
 * this property at the moment (albeit we have a bug if it not holds).
 *
 * @param d1 First domain
 * @param d2 Second domain
 *
 */
Domain intersection(Domain d1, Domain d2);

enum class ConnectionType {
    EDGE,
    CORNER,
    INVALID,
};

/**
 *  @brief Utility struct to return a type of a connection between two domains
 *
 * Since CORNERS and EDGES are represented by different enums we tie
 * the together into a single return type.
 *
 */
class Connection {
    ConnectionType type;
    Edge edge;
    Corner corner;

public:
    static Connection viaEdge(Edge edge)
    {
        Connection c;
        c.type = ConnectionType::EDGE;
        c.edge = edge;
        return c;
    }
    static Connection viaCorner(Corner corner)
    {
        Connection c;
        c.type = ConnectionType::CORNER;
        c.corner = corner;
        return c;
    }
    static Connection invalid()
    {
        Connection c;
        c.type = ConnectionType::INVALID;
        return c;
    }

    ConnectionType getType() const { return type; }
    Edge getEdge() const
    {
        if (type != ConnectionType::EDGE) {
            throw std::runtime_error("Connection is not an EDGE type");
        }
        return edge;
    }
    Corner getCorner() const
    {
        if (type != ConnectionType::CORNER) {
            throw std::runtime_error("Connection is not a CORNER type");
        }
        return corner;
    }
};

/**
 * @brief Identify the CORNER/EDGE relationship
 *
 *
 * We determine the connection relative to `domain` (e.g. for a connectivity)
 *   ┌───┬───┐
 *   │ 1 │ 2 │
 *   └───┴───┘
 * If domain = 1, the result is RIGHT EDGE
 * If domain = 1, the result is LEFT EDGE
 *
 *
 * If a connection is through CANDIDATE, the following holds:
 *  CANDIDATE ∩ SHARED == SHARED and CANDIDATE and SHARED have same dimension
 *
 * Where:
 *   CANDIDATE ∈ {TOP_LEFT, TOP_RIGHT, BOTTOM_RIGHT, BOTTOM_LEFT, LEFT, RIGHT, BOTTOM, TOP}
 *   SHARED beeing non-empty intersection of the two domains
 *
 * @param domain The domain for which we want to identify the connection. It must be an 'Area'.
 * @param shared The shared domain between the two domains.
 *
 *  @result The connection descriptor. Allows to identify EDGE/CORNER.
 *     Invalid connection for:
 *      * Shared domain is empty
 *      * Shared domain is a point on an edge (not a corner)
 *      * Shared domain is an area (not a line or point)
 *      * Shared domain is a line that extends outside the domain
 *
 * @throws If the `domain` is not an area.
 *
 */
Connection identifyNeighbourConnectionType(const Domain domain, const Domain shared);


#endif /* DOMAINUTILS_HPP */
