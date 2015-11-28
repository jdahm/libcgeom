#ifndef CGL_GEOM2D_HPP
#define CGL_GEOM2D_HPP

/*
  Original author: Dani Lischinski
  Source: Graphics Gems IV
*/

#include <iostream>
#include <array>

#include "cgl/cgl.hpp"

namespace cgl
{

class Point2d : private std::array<real, 2> {
public:
        typedef real value_type;
        typedef std::array<real, 2> container_type;
        typedef unsigned int size_type;

        static constexpr size_type dim = 2;

        Point2d();

        Point2d(real x, real y);

        using container_type::begin;
        using container_type::end;
        using container_type::operator[];

        bool operator==(const Point2d& b) const;

        bool operator!=(const Point2d& b) const;

        Point2d operator+(const Point2d& b) const;

        Point2d operator-(const Point2d& b) const;

        value_type distance() const;

        value_type distance(const Point2d& p) const;

        friend std::istream& operator>>(std::istream&, Point2d&);

        friend std::ostream& operator<<(std::ostream&, const Point2d&);
};


class Line2d {
private:
        real a, b, c;
public:
        typedef int classify_type;

        Line2d(const Point2d&, const Point2d&);
        real eval(const Point2d&) const;

        classify_type classify(const Point2d&) const;

        friend Point2d intersect(const Line2d&, const Line2d&);

        friend std::ostream& operator<<(std::ostream&, const Line2d&);
};

// Might want this later
// class AABB2d {
// private:
//         Point2d left_bottom, top_right;
// public:
//         AABB2d(const Point2d&);

//         void add_point(const Point2d&);
// };


real tri_area(const Point2d& a, const Point2d& b,
                            const Point2d& c);

bool in_circle(const Point2d& a, const Point2d& b,
               const Point2d& c, const Point2d& d);

bool ccw(const Point2d& a, const Point2d& b, const Point2d& c);

} // namespace cgl

#endif
