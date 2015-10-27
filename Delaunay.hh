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
  Point2d *addPoint(const Point2d&);
  Edge *addEdge();

  // Operators to remove points and edges, updating pointList and qeList above
  void removeEdge(Edge*);

  // Operators to locate a point or edge given a pointer
  // Output iterator
  PointList::iterator locatePointIter(Point2d* p );
  QEdgeList::iterator locateEdgeIter (Edge *e    );
  // Output index
  PointList::size_type  locatePointIndex(Point2d* p);
  QEdgeList::size_type locateEdgeIndex (Edge *e   );

  // Locate an edge based on coordinates. Attempts to do a smart
  // search based on direction while walking the subdivision
  Edge *locate(const Point2d&) const;

  // Initialize a triangulation (used by constructors)
  void init2(const Point2d&, const Point2d&, Edge*&, Edge*&);
  void init3(const Point2d&, const Point2d&, const Point2d&, Edge*&, Edge*&);
  void initDC(std::list<Point2d>&, Edge*&, Edge*&);

  // Topological Operators for Delaunay Diagram
  Edge* connect(Edge*, Edge*);
  void swap(Edge*);

public:
  // Constructors
  Delaunay(std::list<Point2d>&);
  Delaunay(const Point2d&, const Point2d&, const Point2d&);

  ~Delaunay();

  // Insert a new point and update Delaunay Diagram
  void insertSite(const Point2d&);

  // Simple text file output
  void writeTxt(const std::string&);

  // Paraview friendly output
  void writeVtu(const std::string&);
};


#endif
