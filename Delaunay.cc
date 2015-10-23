#include "Delaunay.hh"
#include "QuadEdge.hh"

/*************** Geometric Predicates for Delaunay Diagrams *****************/

inline Real TriArea(const Point2d& a, const Point2d& b, const Point2d& c)
// Returns twice the area of the oriented triangle (a, b, c), i.e., the
// area is positive if the triangle is oriented counterclockwise.
{
  return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
}

int InCircle(const Point2d& a, const Point2d& b, const Point2d& c, const Point2d& d)
// Returns TRUE if the point d is inside the circle defined by the
// points a, b, c. See Guibas and Stolfi (1985) p.107.
{
  return (a.x*a.x + a.y*a.y) * TriArea(b, c, d) -
    (b.x*b.x + b.y*b.y) * TriArea(a, c, d) +
    (c.x*c.x + c.y*c.y) * TriArea(a, b, d) -
    (d.x*d.x + d.y*d.y) * TriArea(a, b, c) > EPS;
}

int ccw(const Point2d& a, const Point2d& b, const Point2d& c)
// Returns TRUE if the points a, b, c are in a counterclockwiseorder
{
  return (TriArea(a, b, c) > 0);
}

int RightOf(const Point2d& x, Edge* e)
{
  return ccw(x, e->Dest2d(), e->Org2d());
}

int LeftOf(const Point2d& x, Edge* e)
{
  return ccw(x, e->Org2d(), e->Dest2d());
}

int OnEdge(const Point2d& x, Edge* e)
// A predicate that determines if the point x is on the edge e.
// The point is considered on if it is in the EPS-neighborhood
// of the edge.
{
  Real t1, t2, t3;
  t1 = (x - e->Org2d()).norm();
  t2 = (x - e->Dest2d()).norm();
  if (t1 < EPS || t2 < EPS)
    return TRUE;
  t3 = (e->Org2d() - e->Dest2d()).norm();
  if (t1 > t3 || t2 > t3)
    return FALSE;
  Line line(e->Org2d(), e->Dest2d());
  return (fabs(line.eval(x)) < EPS);
}

