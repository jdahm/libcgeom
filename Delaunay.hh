#ifndef DELAUNAY_HH
#define DELAUNAY_HH

#include "Geom2d.hh"
#include "QuadEdge.hh"
#include <list>
#include <tuple>

class Delaunay {
private:
  // List types
  typedef std::list<QuadEdge*> QEdgeList;
  typedef std::list<Point2d*> PointList;

  // Data members
  // Global lists of primal points and edges
  QEdgeList qeList;
  PointList pointList;
  // Keep a pointer to an edge from QuadEdgeList from whick to begin walking
  Edge *startingEdge;

  // Operators to add points and edges updating pointList and qeList above
  Point2d *AddPoint(const Point2d&);
  Edge *AddEdge();

  // Operators to locate a point or edge given a pointer
  // Output iterator
  PointList::iterator LocatePointIter(Point2d* p );
  QEdgeList::iterator LocateEdgeIter (Edge *e    );
  // Output index
  PointList::size_type  LocatePointIndex(Point2d* p);
  QEdgeList::size_type LocateEdgeIndex (Edge *e   );

  // Operators to remove points and edges, updating pointList and qeList above
  void RemoveEdge(Edge*);

  // Locate an edge based on coordinates. Attempts to do a smart
  // search based on direction while walking the subdivision
  Edge *Locate(const Point2d&) const;

  void Init2(const Point2d&, const Point2d&, Edge*&, Edge*&);
  void Init3(const Point2d&, const Point2d&, const Point2d&, Edge*&, Edge*&);
  void InitDC(std::list<Point2d>&, Edge*&, Edge*&);

public:
  Delaunay(std::list<Point2d>&);
  Delaunay(const Point2d&, const Point2d&, const Point2d&);

  ~Delaunay();

  // Insert a new point and update Delaunay Diagram
  void InsertSite(const Point2d&);

  // Topological Operators for Delaunay Diagram
  Edge* Connect(Edge*, Edge*);
  void Swap(Edge*);

  // Simple text file output
  void WriteTxt(const std::string&);

  // Paraview friendly output
  void WriteVtu(const std::string&);
};


#endif
