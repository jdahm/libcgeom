#ifndef DELAUNAY_HH
#define DELAUNAY_HH

#include "Geom2d.hh"
#include "QuadEdge.hh"
#include "Subdivision.hh"
#include <list>

class Delaunay : public Subdivision {
private:
  // List types
  typedef std::list<QuadEdge*> QuadEdgeList;
  typedef std::list<Point2d*> PointList;

  // Data members
  // Global lists of primal points and edges
  QuadEdgeList qeList;
  PointList pointList;
  // Keep a pointer to an edge from QuadEdgeList from whick to begin walking
  Edge *startingEdge;

  // Operators to add points and edges updating pointList and qeList above
  Point2d *AddPoint(const Point2d&);
  Edge *AddEdge();

  // Operators to locate a point or edge given a pointer
  // Output iterator
  PointList::iterator    LocatePointIter(Point2d* p );
  QuadEdgeList::iterator LocateEdgeIter (Edge *e    );
  // Output index
  PointList::size_type    LocatePointIndex(Point2d* p);
  QuadEdgeList::size_type LocateEdgeIndex (Edge *e   );

  // Operators to remove points and edges, updating pointList and qeList above
  void RemoveEdge(Edge*);

  // Locate an edge based on coordinates. Attempts to do a smart
  // search based on direction while walking the subdivision
  Edge *Locate(const Point2d&) const;

public:
  Delaunay(const Point2d&, const Point2d&, const Point2d&);
  ~Delaunay();

  // Insert a new point and update Delaunay Diagram
  void InsertSite(const Point2d&);

  // Topological Operators for Delaunay Diagram
  Edge* Connect(Edge*, Edge*);
  void Swap(Edge*);

  // Implementation of Subdivision interface
  void Write(const std::string&);
};

// /*************** Geometric Predicates for Delaunay Diagrams *****************/

// Real TriArea(const Point2d& a, const Point2d& b, const Point2d& c);

// int InCircle(const Point2d& a, const Point2d& b, const Point2d& c, const Point2d& d);

// bool ccw(const Point2d& a, const Point2d& b, const Point2d& c);

// bool RightOf(const Point2d& x, Edge* e);

// bool LeftOf(const Point2d& x, Edge* e);

// bool OnEdge(const Point2d& x, Edge* e);

#endif
