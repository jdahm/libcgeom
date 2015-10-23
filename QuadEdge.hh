#ifndef QUADEDGE_HH
#define QUADEDGE_HH

#include "Geom2D.hh"

class QuadEdge;

class Edge {
  friend QuadEdge;
  friend void Splice(Edge*, Edge*);
private:
  int num;
  Edge *next;
  Point2d *data;
public:
  Edge() { }

  /************** Edge Algebra ****************/

  inline Edge* Rot()
  // Return the dual of the current edge, directed from its right to its left.
  { return (num < 3) ? this + 1 : this - 3; }

  inline Edge* invRot()
  // Return the dual of the current edge, directed from its left to its right.
  { return (num > 0) ? this - 1 : this + 3; }

  inline Edge* Sym()
  // Return the edge from the destination to the origin of the current edge.
  { return (num < 2) ? this + 2 : this - 2; }

  inline Edge* Onext()
  // Return the next ccw edge around (from) the origin of the current edge.
  { return next; }

  inline Edge* Oprev()
  // Return the next cw edge around (from) the origin of the current edge.
  { return Rot()->Onext()->Rot(); }

  inline Edge* Dnext()
  // Return the next ccw edge around (into) the destination of the current edge.
  { return Sym()->Onext()->Sym(); }

  inline Edge* Dprev()
  // Return the next cw edge around (into) the destination of the current edge.
  { return invRot()->Onext()->invRot(); }

  inline Edge* Lnext()
  // Return the ccw edge around the left face following the current edge.
  { return invRot()->Onext()->Rot(); }

  inline Edge* Lprev()
  // Return the ccw edge around the left face before the current edge.
  { return Onext()->Sym(); }

  inline Edge* Rnext()
  // Return the edge around the right face ccw following the current edge.
  { return Rot()->Onext()->invRot(); }

  inline Edge* Rprev()
  // Return the edge around the right face ccw before the current edge.
  { return Sym()->Onext(); }

  /************** Access to data pointers ****************/
  inline Point2d* Org()
  // Origin
  { return data; }

  inline Point2d* Dest()
  // Destination
  { return Sym()->data; }

  inline const Point2d& Org2d() const
  { return *data; }

  inline const Point2d& Dest2d() const
  { return (num < 2) ? *((this + 2)->data) : *((this - 2)->data); }

  inline void  EndPoints(Point2d* o, Point2d* d)
  {
    data = o;
    Sym()->data = d;
  }

  inline QuadEdge* Qedge() { return (QuadEdge *)(this - num); }
};

class QuadEdge {
  friend Edge* MakeEdge();
private:
  Edge e[4];
public:
  QuadEdge();
};


/*********************** Basic Topological Operators ****************/
Edge* MakeEdge();

void Splice(Edge* a, Edge* b);

void DeleteEdge(Edge* e);

Edge* Connect(Edge* a, Edge* b);

#endif
