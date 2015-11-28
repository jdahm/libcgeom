#include <cassert>
#include <cmath>

#include "cgl/geom2d.hpp"

namespace cgl
{

/* Point2d */

Point2d::Point2d() : Point2d::container_type({0, 0}) { }

Point2d::Point2d(real x, real y) : Point2d::container_type({x, y}) { }

bool Point2d::operator==(const Point2d& b) const {
        const Point2d& a = *this;
        for (size_type i=0; i<dim; i++)
                if (std::abs(a[i] - b[i]) > real_eps) return false;
        return true;
}

bool Point2d::operator!=(const Point2d& b) const {
        const Point2d& a = *this;
        return ! (a == b);
}

Point2d Point2d::operator+(const Point2d& b) const {
        const Point2d& a = *this;
        return Point2d(a[0] + b[0], a[1] + b[1]);
}

Point2d Point2d::operator-(const Point2d& b) const {
        const Point2d& a = *this;
        return Point2d(a[0] - b[0], a[1] - b[1]);
}

Point2d::value_type Point2d::distance() const {
        const Point2d& a = *this;
        return std::sqrt(a[0]*a[0] + a[1]*a[1]);
}

Point2d::value_type Point2d::distance(const Point2d& p) const {
        const Point2d& a = *this;
        real dist = 0;
        for (size_type i=0; i<dim; i++)
                dist += (a[i] - p[i]) * (a[i] - p[i]);
        return std::sqrt(dist);
}

std::istream& operator>> (std::istream& is, Point2d& p) {
        is >> p[0];
        is >> p[1];
        return is;
}

std::ostream& operator<< (std::ostream& os, const Point2d& p) {
        os << "(" << p[0] << "," << p[1] << ")";
        return os;
}

/* Line2d */

Line2d::Line2d(const Point2d& p, const Point2d& q)
{
        // Computes the normalized line equation through the points p
        // and q.
        Point2d t = q - p;
        real len = t.distance();
        assert(len != 0);
        a =  t[1] / len;
        b = -t[0] / len;
        c = -(a * p[0] + b * p[1]);
}


real Line2d::eval(const Point2d& p) const
{
        // Plugs point p into the line equation.
        return a * p[0] + b * p[1] + c;
}


Line2d::classify_type Line2d::classify(const Point2d& p) const
{
        // Returns -1, 0, or 1, if p is to the left of, on,
        // or right of the line, respectively.
        real d = eval(p);
        return (d < -real_eps) ? -1 : (d > real_eps ? 1 : 0);
}

Point2d intersect(const Line2d& k, const Line2d &l)
{
        // intersect: point of intersection p of two lines k, l
        real den = k.a * l.b - k.b * l.a;
        assert(den != 0);
        Point2d inter = {(k.b * l.c - k.c * l.b) / den, (k.c * l.a - k.a * l.c) / den};
        return inter;
}

std::ostream& operator<< (std::ostream &os, const Line2d& k)
{
        os << "Line2d(a=" << k.a << " b=" << k.b << " c=" << k.c << ")";
        return os;
}

// // Initialize the bounding box with zero area
// AABBB2d::AABB2d(const Point2d& p) : left_bottom(p), right_top(p) { }

// void AABBB2d::add_point(const Point2d& p)
// {
        
// }


// Returns twice the area of the oriented triangle (a, b, c), i.e., the
// area is positive if the triangle is oriented counterclockwise.
real tri_area(const Point2d& a, const Point2d& b, const Point2d& c)
{
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}


// Returns true if the point d is inside the circle defined by the
// points a, b, c. See Guibas and Stolfi (1985) p.107.
bool in_circle(const Point2d& a, const Point2d& b,
               const Point2d& c, const Point2d& d)
{
        return (a[0] * a[0] + a[1] * a[1]) * tri_area(b, c, d) -
               (b[0] * b[0] + b[1] * b[1]) * tri_area(a, c, d) +
               (c[0] * c[0] + c[1] * c[1]) * tri_area(a, b, d) -
               (d[0] * d[0] + d[1] * d[1]) * tri_area(a, b, c) > real_eps;
}


// Returns true if the points a, b, c are in a counterclockwiseorder
bool ccw(const Point2d& a, const Point2d& b, const Point2d& c)
{
        return tri_area(a, b, c) > 0;
}

} // namespace cgl
