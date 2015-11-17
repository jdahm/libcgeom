#include "Geom2d.H"
#include <iostream>
#include <fstream>
/******** Point2d *********/

Real Point2d::norm() const {
        return sqrt(x * x + y * y);
}


void Point2d::normalize() {
        Real len;

        if ((len = sqrt(x * x + y * y)) == 0.0) {
                std::cerr << "Point2d::normalize: Division by 0\n";
        } else {
                x /= len;
                y /= len;
        }
}


Point2d operator* (Real c, const Point2d& v) {
        return Point2d(c * v.x, c * v.y);
}


// static Real dot(const Point2d& u, const Point2d& v) {
//         return u.x * v.x + u.y * v.y;
// }


Point2d Point2d::operator+ (const Point2d& v) const {
        return Point2d(x + v.x, y + v.y);
}


Point2d Point2d::operator- (const Point2d& p) const {
        return Point2d(x - p.x, y - p.y);
}


bool Point2d::operator== (const Point2d& p) const {
        return (*this - p).norm() < RealEps;
}


std::istream& operator>> (std::istream& is, Point2d& p) {
        is >> p.x >> p.y;
        return is;
}


std::ostream& operator<< (std::ostream& os, const Point2d& p) {
        os << '(' << p.x << "," << p.y << ')';
        Real fx = p.x - floor(p.x + .5);
        Real fy = p.y - floor(p.y + .5);
        if ((fx!=0) && (fy!=0) && (fabs(fx)<1e-3) && (fabs(fy)<1e-3))
                os << "=("
                   << floor(p.x + .5) << "+"[fx<0] << fx << ","
                   << floor(p.y + .5) << "+"[fy<0] << fy << ")";

        return os;
}


/******** Line *********/

Line::Line(const Point2d& p, const Point2d& q) {
// Computes the normalized line equation through the
// points p and q.
        Point2d t = q - p;
        Real len = t.norm();
        assert(len!=0);
        a =  t.y / len;
        b = -t.x / len;
        c = -(a * p.x + b * p.y);
}


// Plugs point p into the line equation.
Real Line::eval(const Point2d& p) const {
        return a * p.x + b * p.y + c;
}


// Returns -1, 0, or 1, if p is to the left of, on,
// or right of the line, respectively.
int Line::classify(const Point2d& p) const {
        Real d = eval(p);
        return (d < -RealEps) ? -1 : (d > RealEps ? 1 : 0);
}


std::ostream& operator<< (std::ostream &os, const Line &k) {
        os << "Line(a=" << k.a << " b=" << k.b << " c=" << k.c << ")";
        return os;
}
