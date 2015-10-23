#ifndef SUBDIVISION_HH
#define SUBDIVISION_HH

#include "QuadEdge.hh"
#include "Geom2D.hh"
#include <list>

class Subdivision {
private:
  // List types
  typedef std::list<QuadEdge*> QuadEdgeList;
  typedef std::list<Point2d*> PointList;

  QuadEdgeList qeList;
  PointList pointList;
  Edge *startingEdge;

  // Operators to add points and edges updating pointList and qeList above
  Point2d *AddPoint(const Point2d&);
  Edge *AddEdge();

  // Operators to remove points and edges, updating pointList and qeList above
  void RemoveEdge(Edge*);

  PointList::size_type containsPoint(const Point2d&);

  Edge *Locate(const Point2d&);

public:
  Subdivision(const Point2d&, const Point2d&, const Point2d&);
  ~Subdivision();

  // Insert a new point and update Delaunay Diagram
  void InsertSite(const Point2d&);

  // Topological Operators for Delaunay Diagram
  Edge* Connect(Edge* a, Edge* b);
  void Swap(Edge* e);


};

#endif
