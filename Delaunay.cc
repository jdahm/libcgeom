#include "Delaunay.hh"
#include "QuadEdge.hh"
#include <algorithm>
#include <fstream>
#include <iomanip>

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


/************* Topological Operations for Delaunay Diagrams *****************/

Edge* Delaunay::Connect(Edge* a, Edge* b)
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

void Delaunay::Swap(Edge* e)
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

Delaunay::Delaunay(const Point2d& a, const Point2d& b, const Point2d& c)
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

Delaunay::~Delaunay()
{
  for(Point2d* p : pointList)
    delete p;

  for(QuadEdge* qe : qeList)
    delete qe;
}

Point2d* Delaunay::AddPoint(const Point2d& a)
// Returns a pointer to the point at index in pointList
{
  Point2d* da = new Point2d(a);
  pointList.push_back(da);
  // return pointList.back();
  return da;
}

Edge* Delaunay::AddEdge()
{
  Edge* eb = MakeEdge();
  qeList.push_back(eb->Qedge());
  // return eqList.back();
  return eb;
}

Delaunay::PointList::iterator Delaunay::LocatePointIter(Point2d* p)
// This step is relatively expensive: O(nEdge-1)
{
  return std::find(pointList.begin(), pointList.end(), p);
}

Delaunay::QuadEdgeList::iterator Delaunay::LocateEdgeIter(Edge *e)
// This step is relatively expensive: O(nEdge-1)
{
  return std::find(qeList.begin(), qeList.end(), e->Qedge());
}

Delaunay::PointList::size_type Delaunay::LocatePointIndex(Point2d* p)
{
  return std::distance(pointList.begin(), LocatePointIter(p));
}

Delaunay::QuadEdgeList::size_type Delaunay::LocateEdgeIndex(Edge *e)
{
  return std::distance(qeList.begin(), LocateEdgeIter(e));
}

void Delaunay::RemoveEdge(Edge* e)
{
  QuadEdgeList::iterator qeIter = LocateEdgeIter(e);
  DeleteEdge(e); // Deallocates the pointer
  qeList.erase(qeIter); // Removes the element from the list
}


Edge* Delaunay::Locate(const Point2d& x) const
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


/************* Iterative Delaunay Diagram *****************/

void Delaunay::InsertSite(const Point2d& x)
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

void Delaunay::Write(const std::string& fileName)
// Outputs the delaunay subdivision to fileName
// Format:
// nVerts nEdges
// for i in [1,nVerts]:
//   Vert_i->x Vert_i->y
// for i in [1,nEdges]:
//   Edge_i->Org Edge_i->Dest
{
  std::ofstream fst(fileName, std::ios::out);

  fst << pointList.size() << " " << qeList.size() << std::endl;

  for (Point2d* p : pointList)
    fst << std::scientific << std::setprecision(16) << p->x << " "
        << std::scientific << std::setprecision(16) << p->y << std::endl;

  for (QuadEdge* qe : qeList){
    Edge* e = reinterpret_cast<Edge*>(qe);
    PointList::size_type oi = LocatePointIndex(e->Org()), di = LocatePointIndex(e->Dest());
    fst << oi << " " << di << std::endl;
  }

  fst.close();
}
