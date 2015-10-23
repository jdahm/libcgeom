#ifndef DELAUNAY_HH
#define DELAUNAY_HH

#include "Geom2D.hh"
#include "QuadEdge.hh"

/************* Topological Operations for Delaunay Diagrams *****************/

Edge* Connect(Edge* a, Edge* b);

void Swap(Edge* e);

/*************** Geometric Predicates for Delaunay Diagrams *****************/

Real TriArea(const Point2d& a, const Point2d& b, const Point2d& c);

int InCircle(const Point2d& a, const Point2d& b, const Point2d& c, const Point2d& d);

int ccw(const Point2d& a, const Point2d& b, const Point2d& c);

int RightOf(const Point2d& x, Edge* e);

int LeftOf(const Point2d& x, Edge* e);

int OnEdge(const Point2d& x, Edge* e);

#endif
