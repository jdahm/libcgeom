#include "Subdivision.hh"
#include "QuadEdge.hh"
#include "Delaunay.hh"
#include <algorithm>

Subdivision::Subdivision(const Point2d& a, const Point2d& b, const Point2d& c)
// Initialize a subdivision to the triangle defined by the points a, b, c.
{
  Point2d *da = AddPoint(a), *db = AddPoint(b), *dc = AddPoint(c);
  Edge* ea = AddEdge();
  ea->EndPoints(da, db);
  Edge* eb = AddEdge();
  Splice(ea->Sym(), eb);
  eb->EndPoints(db, dc);
  Edge* ec = AddEdge();
  Splice(eb->Sym(), ec);
  ec->EndPoints(dc, da);
  Splice(ec->Sym(), ea);
  startingEdge = ea;
}

Subdivision::~Subdivision()
{
  for(Point2d* p : pointList)
    delete p;

  for(QuadEdge* qe : qeList)
    delete qe;
}

Point2d* Subdivision::AddPoint(const Point2d& a)
// Returns a pointer to the point at index in pointList
{
  Point2d* da = new Point2d(a);
  pointList.push_back(da);
  // return pointList.back();
  return da;
}

Edge* Subdivision::AddEdge()
{
  Edge* eb = MakeEdge();
  qeList.push_back(eb->Qedge());
  // return eqList.back();
  return eb;
}

void Subdivision::RemoveEdge(Edge* e)
{
  // This step is expensive: O(n)
  QuadEdgeList::iterator qeIter = std::find(qeList.begin(), qeList.end(), e->Qedge());
  DeleteEdge(e); // Deallocates the pointer
  qeList.erase(qeIter); // Removes the element from the list
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


/************* Topological Operations for Delaunay Diagrams *****************/

Edge* Subdivision::Connect(Edge* a, Edge* b)
// Add a new edge e connecting the destination of a to the
// origin of b, in such a way that all three have the same
// left face after the connection is complete.
// Additionally, the data pointers of the new edge are set.
{
  Edge* e = AddEdge();
  Splice(e, a->Lnext());
  Splice(e->Sym(), b);
  e->EndPoints(a->Dest(), b->Org());
  return e;
}

void Subdivision::Swap(Edge* e)
// Essentially turns edge e counterclockwise inside its enclosing
// quadrilateral. The data pointers are modified accordingly.
{
  Edge* a = e->Oprev();
  Edge* b = e->Sym()->Oprev();
  Splice(e, a);
  Splice(e->Sym(), b);
  Splice(e, a->Lnext());
  Splice(e->Sym(), b->Lnext());
  e->EndPoints(a->Dest(), b->Dest());
}


/************* Iterative Delaunay Diagram *****************/

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
    RemoveEdge(e->Onext());
  }
  // Connect the new point to the vertices of the containing
  // triangle (or quadrilateral, if the new point fell on an
  // existing edge.)
  Edge* base = AddEdge();
  Point2d *dx = AddPoint(x);
  base->EndPoints(e->Org(), dx);
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
