#ifndef SUBDIVISION_HH
#define SUBDIVISION_HH

#include "QuadEdge.hh"
#include "Geom2D.hh"

class Subdivision {
private:
  Edge *startingEdge;
  Edge *Locate(const Point2d&);
public:
  Subdivision(const Point2d&, const Point2d&, const Point2d&);
  void InsertSite(const Point2d&);
};

#endif
