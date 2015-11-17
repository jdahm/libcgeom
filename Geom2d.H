#ifndef GEOM2D_H
#define GEOM2D_H

////////////////////////////////////////////////////////////////////////
// This code is a modified version of the Delaunay triangulator
// written by Dani Lischinski
// in Graphics Gems IV
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include <cassert>
#include <limits>

typedef double Real;
static constexpr Real RealEps = std::numeric_limits<Real>::epsilon();

class Point2d {
public:
        Real x, y;
        Point2d() : x(0.0), y(0.0) { }
        Point2d(Real a, Real b) : x(a), y(b) { }
        Point2d(const Point2d &p) {
                *this = p;
        }

        Real norm() const;

        void normalize();

        Point2d operator + (const Point2d&) const;
        Point2d operator - (const Point2d&) const;
        bool operator==(const Point2d&) const;
};

class Line {
private:
        Real a, b, c;
public:
        Line() { }
        Line(const Point2d &, const Point2d &);
        Real eval(const Point2d&) const;

        int classify(const Point2d&) const;

        friend void intersect(const Line &k, const Line &l, Point2d &p);

        friend std::ostream& operator << (std::ostream &, const Line &);
};


#endif
