#include "Delaunay.hh"
#include "QuadEdge.hh"
#include <algorithm>
#include <fstream>
#include <iomanip>


/*************** Geometric Predicates for Delaunay Diagrams *****************/

static inline Real triArea(const Point2d& a, const Point2d& b, const Point2d& c)
// Returns twice the area of the oriented triangle (a, b, c), i.e., the
// area is positive if the triangle is oriented counterclockwise.
{
  return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
}

static inline bool inCircle(const Point2d& a, const Point2d& b, const Point2d& c, const Point2d& d)
// Returns TRUE if the point d is inside the circle defined by the
// points a, b, c. See Guibas and Stolfi (1985) p.107.
{
  return (a.x*a.x + a.y*a.y) * triArea(b, c, d) -
    (b.x*b.x + b.y*b.y) * triArea(a, c, d) +
    (c.x*c.x + c.y*c.y) * triArea(a, b, d) -
    (d.x*d.x + d.y*d.y) * triArea(a, b, c) > RealEps;
}

static inline bool ccw(const Point2d& a, const Point2d& b, const Point2d& c)
// Returns TRUE if the points a, b, c are in a counterclockwiseorder
{
  return (triArea(a, b, c) > 0);
}

static inline bool rightOf(const Point2d& x, Edge* e)
{
  return ccw(x, e->Dest2d(), e->Org2d());
}

static inline bool leftOf(const Point2d& x, Edge* e)
{
  return ccw(x, e->Org2d(), e->Dest2d());
}

static bool onEdge(const Point2d& x, Edge* e)
// A predicate that determines if the point x is on the edge e.
// The point is considered on if it is in the EPS-neighborhood
// of the edge.
{
  Real t1, t2, t3;
  t1 = (x - e->Org2d()).norm();
  t2 = (x - e->Dest2d()).norm();
  if (t1 < RealEps || t2 < RealEps)
    return true;
  t3 = (e->Org2d() - e->Dest2d()).norm();
  if (t1 > t3 || t2 > t3)
    return false;
  Line line(e->Org2d(), e->Dest2d());
  return (fabs(line.eval(x)) < RealEps);
}


/************* Topological Operations for Delaunay Diagrams *****************/

Edge* Delaunay::connect(Edge* a, Edge* b)
// Add a new edge e connecting the destination of a to the
// origin of b, in such a way that all three have the same
// left face after the connection is complete.
// Additionally, the data pointers of the new edge are set.
{
  Edge* e = addEdge();
  Splice(e, a->Lnext());
  Splice(e->Sym(), b);
  e->EndPoints(a->Dest(), b->Org());
  return e;
}

void Delaunay::swap(Edge* e)
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

void Delaunay::init2(const Point2d& a, const Point2d& b,
                     Edge* &el, Edge* &er)
// Initialize a subdivision to the line defined by the points a, b.
{
  Point2d *da = addPoint(a), *db = addPoint(b);
  Edge* ea = addEdge();
  ea->EndPoints(da, db);
  startingEdge = ea;
  el = ea;
  er = ea->Sym();
}

void Delaunay::init3(const Point2d& a, const Point2d& b, const Point2d& c,
                     Edge* &el, Edge* &er)
// Initialize a subdivision to the triangle defined by the points a, b, c.
{
  Point2d *da = addPoint(a), *db = addPoint(b), *dc = addPoint(c);
  // Assume a, b, and c are in sorted order
  // Create edges ea and eb connecting da to db and db to dc
  Edge *ea = addEdge(), *eb = addEdge();
  Splice(ea->Sym(), eb);
  ea->EndPoints(da, db);
  eb->EndPoints(ea->Dest(), dc);
  // Now close the triangle
  Edge *ec = connect(eb, ea);
  el = ea;
  er = eb->Sym();
  if (ccw(a, c, b)) {
    el = ec->Sym();
    er = ec;
  }
}

Delaunay::Delaunay(const Point2d& a, const Point2d& b, const Point2d& c)
{
  Edge *el, *er;
  init3(a, b, c, el, er);
  startingEdge = el;
}

Delaunay::Delaunay(std::list<Point2d>& point)
{
  Edge *el, *er;
  initDC(point, el, er);
  startingEdge = el;
}

Delaunay::~Delaunay()
{
  for(Point2d* p : pointList)
    delete p;

  for(QuadEdge* qe : qeList)
    delete qe;
}

Point2d* Delaunay::addPoint(const Point2d& a)
// Allocate a Point2d and add it to pointList
{
  Point2d* da = new Point2d(a);
  pointList.push_back(da);
  return da;
}

Edge* Delaunay::addEdge()
// Allocate an edge and add it to qeList
{
  Edge* eb = MakeEdge();
  qeList.push_back(eb->Qedge());
  return eb;
}

void Delaunay::removeEdge(Edge* e)
// Remove an edge and deallocate the memory
{
  QEdgeList::iterator qeIter = locateEdgeIter(e);
  DeleteEdge(e); // Deallocates the pointer
  qeList.erase(qeIter); // Removes the element from the list
}

Delaunay::PointList::iterator Delaunay::locatePointIter(Point2d* p)
// This step is relatively expensive: O(nEdge-1)
{
  return std::find(pointList.begin(), pointList.end(), p);
}

Delaunay::QEdgeList::iterator Delaunay::locateEdgeIter(Edge *e)
// This step is relatively expensive: O(nEdge-1)
{
  return std::find(qeList.begin(), qeList.end(), e->Qedge());
}

Delaunay::PointList::size_type Delaunay::locatePointIndex(Point2d* p)
// Convert pointer distance to size_type
{
  return std::distance(pointList.begin(), locatePointIter(p));
}

Delaunay::QEdgeList::size_type Delaunay::locateEdgeIndex(Edge *e)
// Convert pointer distance to size_type
{
  return std::distance(qeList.begin(), locateEdgeIter(e));
}

Edge* Delaunay::locate(const Point2d& x) const
// Returns an edge e, s.t. either x is on e, or e is an edge of
// a triangle containing x. The search starts from startingEdge
// and proceeds in the general direction of x. Based on the
// pseudocode in Guibas and Stolfi (1985) p.121.
{
  Edge* e = startingEdge;
  while (true) {
    if (x == e->Org2d() || x == e->Dest2d())
      return e;
    else if (rightOf(x, e))
      e = e->Sym();
    else if (!rightOf(x, e->Onext()))
      e = e->Onext();
    else if (!rightOf(x, e->Dprev()))
      e = e->Dprev();
    else
      return e;
  }
}


/************* Divide and Conquer (DC) Delaunay Diagram *****************/

static inline bool valid(Edge* e, Edge* basel)
// Tests whether the edge e is above basel
{ return rightOf(*e->Dest(), basel); }

void Delaunay::initDC(std::list<Point2d>& point, Edge* &el, Edge* &er)
// Creates a Delaunay triangulation from the points listed in 'point',
// using the Divide and Conquer algorithm from Guibas
// and Stolfi (1985) Fig. 23 p. 114.
{
  if (point.size() == 2) {
    // Create an edge from one to the other point
    init2(point.front(), point.back(), el, er);
  }
  else if (point.size() == 3) {
    // Create a triangle
    std::list<Point2d>::iterator pi = point.begin();
    init3(*pi, *std::next(pi), *std::next(pi, 2), el, er);
  }
  else{
    // Split into two lists
    std::list<Point2d> pointL, pointR;
    std::list<Point2d>::iterator mpit = std::next(point.begin(), (point.size()+1)/2);
    pointL.splice(pointL.begin(), point, point.begin(), mpit);
    pointR.splice(pointR.begin(), point, mpit, point.end());

    // Compute Delaunay recusively until one of the conditions above is met
    Edge *ldo, *ldi, *rdi, *rdo;
    initDC(pointL, ldo, ldi);
    initDC(pointR, rdi, rdo);

    // Compute the lower tangent of L and R
    while (true) {
      if      (leftOf (*rdi->Org(), ldi)) ldi = ldi->Lnext();
      else if (rightOf(*ldi->Org(), rdi)) rdi = rdi->Rprev();
      else break;
    }

    // Create a first cross edge basel from rdi.Org to ldi.Org
    Edge *basel = connect(rdi->Sym(), ldi);
    if (ldi->Org() == ldo->Org()) ldo = basel->Sym();
    if (rdi->Org() == rdo->Org()) rdo = basel;

    // This is the merge loop
    while (true) {
      // Locate the first L point (lcand.Dest) to be encountered by
      // the rising bubble, and delete L edges out of basel.Dest that
      // fail the circle test.
      Edge *lcand = basel->Sym()->Onext();
      if (valid(lcand, basel))
        while(inCircle(*basel->Dest(), *basel->Org(), *lcand->Dest(), *lcand->Onext()->Dest())) {
          Edge *t = lcand->Onext();
          DeleteEdge(lcand);
          lcand = t;
        }
      // Symmetrically, locate the first R point to be hit, and delete
      // R edges
      Edge *rcand = basel->Oprev();
      if (valid(rcand, basel))
        while(inCircle(*basel->Dest(), *basel->Org(), *rcand->Dest(), *rcand->Oprev()->Dest())) {
          Edge *t = rcand->Oprev();
          DeleteEdge(rcand);
          rcand = t;
        }

      // If both lcand and rcand are invalid, then basel is the upper
      // common tangent
      if (!valid(lcand, basel) && !valid(rcand, basel)) break;

      // The next cross edge is to be connected to either lcand.Dest
      // or rcand.Dest If both are valid, then choose the appropriate
      // one using the InCircle test.
      if (!valid(lcand, basel) ||
          (valid(rcand, basel) && inCircle(*lcand->Dest(), *lcand->Org(),
                                           *rcand->Org(), *rcand->Dest()))) {
        // Add cross edge basel from rcand.Dest to basel.Dest
        basel = connect(rcand, basel->Sym());
      }
      else {
        // Add cross edge basel from basel.Org to lcand.Dest
        basel = connect(basel->Sym(), lcand->Sym());
      }
    } // Merge loop

    // Return [ldo, rdo]
    el = ldo;
    er = rdo;
  }
}

/************* Iterative Delaunay Diagram *****************/

void Delaunay::insertSite(const Point2d& x)
// Inserts a new point into a subdivision representing a Delaunay
// triangulation, and fixes the affected edges so that the result
// is still a Delaunay triangulation. This is based on the
// pseudocode from Guibas and Stolfi (1985) p.120, with slight
// modifications and a bug fix.
{
  Edge* e = locate(x);
  if ((x == e->Org2d()) || (x == e->Dest2d()))  // point is already in
    return;
  else if (onEdge(x, e)) {
    e = e->Oprev();
    removeEdge(e->Onext());
  }
  // Connect the new point to the vertices of the containing
  // triangle (or quadrilateral, if the new point fell on an
  // existing edge.)
  Edge* base = addEdge();
  Point2d *dx = addPoint(x);
  base->EndPoints(e->Org(), dx);
  Splice(base, e);
  startingEdge = base;
  do {
    base = connect(e, base->Sym());
    e = base->Oprev();
  } while (e->Lnext() != startingEdge);
  // Examine suspect edges to ensure that the Delaunay condition
  // is satisfied.
  do {
    Edge* t = e->Oprev();
    if (rightOf(t->Dest2d(), e) &&
        inCircle(e->Org2d(), t->Dest2d(), e->Dest2d(), x)) {
      swap(e);
      e = e->Oprev();
    }
    else if (e->Onext() == startingEdge)  // no more suspect edges
      return;
    else  // pop a suspect edge
      e = e->Onext()->Lprev();
  } while (true);
}
