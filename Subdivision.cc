#include "Subdivision.hh"
#include "QuadEdge.hh"
#include "Delaunay.hh"


Subdivision::Subdivision(const Point2d& a, const Point2d& b, const Point2d& c)
// Initialize a subdivision to the triangle defined by the points a, b, c.
{
  Point2d *da, *db, *dc;
  da = new Point2d(a), db = new Point2d(b), dc = new Point2d(c);
  Edge* ea = MakeEdge();
  ea->EndPoints(da, db);
  Edge* eb = MakeEdge();
  Splice(ea->Sym(), eb);
  eb->EndPoints(db, dc);
  Edge* ec = MakeEdge();
  Splice(eb->Sym(), ec);
  ec->EndPoints(dc, da);
  Splice(ec->Sym(), ea);
  startingEdge = ea;
}

Edge* Subdivision::Locate(const Point2d& x)
// Returns an edge e, s.t. either x is on e, or e is an edge of
// a triangle containing x. The search starts from startingEdge
// and proceeds in the general direction of x. Based on the
// pseudocode in Guibas and Stolfi (1985) p.121.
{
  Edge* e = startingEdge;
  while (TRUE) {
    if (x == e->Org2d() || x == e->Dest2d())
      return e;
    else if (RightOf(x, e))
      e = e->Sym();
    else if (!RightOf(x, e->Onext()))
      e = e->Onext();
    else if (!RightOf(x, e->Dprev()))
      e = e->Dprev();
    else
      return e;
  }
}

void Subdivision::InsertSite(const Point2d& x)
// Inserts a new point into a subdivision representing a Delaunay
// triangulation, and fixes the affected edges so that the result
// is still a Delaunay triangulation. This is based on the
// pseudocode from Guibas and Stolfi (1985) p.120, with slight
// modifications and a bug fix.
{
  Edge* e = Locate(x);
  if ((x == e->Org2d()) || (x == e->Dest2d()))  // point is already in
    return;
  else if (OnEdge(x, e)) {
    e = e->Oprev();
    DeleteEdge(e->Onext());
  }
  // Connect the new point to the vertices of the containing
  // triangle (or quadrilateral, if the new point fell on an
  // existing edge.)
  Edge* base = MakeEdge();
  base->EndPoints(e->Org(), new Point2d(x));
  Splice(base, e);
  startingEdge = base;
  do {
    base = Connect(e, base->Sym());
    e = base->Oprev();
  } while (e->Lnext() != startingEdge);
  // Examine suspect edges to ensure that the Delaunay condition
  // is satisfied.
  do {
    Edge* t = e->Oprev();
    if (RightOf(t->Dest2d(), e) &&
        InCircle(e->Org2d(), t->Dest2d(), e->Dest2d(), x)) {
      Swap(e);
      e = e->Oprev();
    }
    else if (e->Onext() == startingEdge)  // no more suspect edges
      return;
    else  // pop a suspect edge
      e = e->Onext()->Lprev();
  } while (TRUE);
}
