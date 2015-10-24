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

using std::istream; using std::ostream; using std::cerr; using std::cout;

#ifndef ABS
#define ABS(a)	((a) >= 0 ? (a) : -(a))
#endif

#ifndef MAX
#define MAX(a, b)       ((a) >= (b) ? (a) : (b))
#define MIN(a, b)       ((a) <= (b) ? (a) : (b))
#endif

#ifndef TRUE
#define FALSE 0
#define TRUE  1
#endif

#define EPS 1e-6

typedef double  Real;

class Point2d {
public:
	Real x, y;
	Point2d()					{ x = 0; y = 0; }
	Point2d(Real a, Real b)		{ x = a; y = b; }
	Point2d(const Point2d& p)	{ *this = p; }
	Real norm() const;
	void normalize();
	Point2d operator+(const Point2d&) const;
	Point2d operator-(const Point2d&) const;
	int operator==(const Point2d&) const;
	friend Point2d operator*(Real, const Point2d&);
	friend Real dot(const Point2d&, const Point2d&);
	friend istream& operator>>(istream&, Point2d&);
	friend ostream& operator<<(ostream&, const Point2d&);
};

class Line {
public:
	Line()	{}
	Line(const Point2d&, const Point2d&);
	Real eval(const Point2d&) const;
	int classify(const Point2d&) const;
	friend void intersect(const Line &k, const Line &l, Point2d &p);
	friend ostream& operator<<(ostream&, const Line&);
private:
	Real a, b, c;
};

inline Real Point2d::norm() const
{
	return sqrt(x * x + y * y);
}

inline void Point2d::normalize()
{
	Real len;

	if ((len = sqrt(x * x + y * y)) == 0.0)
		cerr << "Point2d::normalize: Division by 0\n";
	else {
		x /= len;
		y /= len;
	}
}

inline Point2d operator*(Real c, const Point2d& v)
{
	return Point2d(c * v.x, c * v.y);
}

inline Real dot(const Point2d& u, const Point2d& v)
{
	return u.x * v.x + u.y * v.y;
}

inline Point2d Point2d::operator+(const Point2d& v) const
{
	return Point2d(x + v.x, y + v.y);
}

inline Point2d Point2d::operator-(const Point2d& p) const
{
	return Point2d(x - p.x, y - p.y);
}

inline int Point2d::operator==(const Point2d& p) const
{
	return ((*this - p).norm() < EPS);
}

inline istream& operator>>(istream& is, Point2d& p)
{
	is >> p.x >> p.y;
	return is;
}

inline ostream& operator<<(ostream& os, const Point2d& p)
{
	os << '(' << p.x << "," << p.y << ')';
	//??following lines to help with debugging
	Real fx = p.x-floor(p.x+.5);
	Real fy = p.y-floor(p.y+.5);
	if (fx!=0 && fy!=0 && fabs(fx)<1e-3 && fabs(fy)<1e-3)
		os << "=("
		    << floor(p.x+.5) << "+"[fx<0] << fx << ","
		    << floor(p.y+.5) << "+"[fy<0] << fy << ")";
	return os;
}

// Line:

inline Line::Line(const Point2d& p, const Point2d& q)
// Computes the normalized line equation through the
// points p and q.
{
	Point2d t = q - p;
	Real len = t.norm();
	assert(len!=0);
	a =   t.y / len;
	b = - t.x / len;
	c = -(a*p.x + b*p.y);
}

inline Real Line::eval(const Point2d& p) const
// Plugs point p into the line equation.
{
	return (a * p.x + b* p.y + c);
}

inline int Line::classify(const Point2d& p) const
// Returns -1, 0, or 1, if p is to the left of, on,
// or right of the line, respectively.
{
	Real d = eval(p);
	return (d < -EPS) ? -1 : (d > EPS ? 1 : 0);
}

inline ostream& operator<<(ostream &os, const Line &k)
{
	os << "Line(a=" << k.a << " b=" << k.b << " c=" << k.c << ")";
	return os;
}

#endif
