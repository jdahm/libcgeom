#include "Delaunay.H"
#include "QuadEdge.H"
#include <algorithm>
#include <fstream>
#include <iomanip>


/*************** Geometric Predicates for Delaunay Diagrams *****************/

// Returns twice the area of the oriented triangle (a, b, c), i.e., the
// area is positive if the triangle is oriented counterclockwise.
static inline Real tri_area(const Point2d& a, const Point2d& b, const
                            Point2d& c) {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}


// Returns true if the point d is inside the circle defined by the
// points a, b, c. See Guibas and Stolfi (1985) p.107.
static inline bool in_circle(const Point2d& a, const Point2d& b, const
                             Point2d& c, const Point2d& d) {
        return (a.x * a.x + a.y * a.y) * tri_area(b, c, d) -
               (b.x * b.x + b.y * b.y) * tri_area(a, c, d) +
               (c.x * c.x + c.y * c.y) * tri_area(a, b, d) -
               (d.x * d.x + d.y * d.y) * tri_area(a, b, c) > RealEps;
}


// Returns true if the points a, b, c are in a counterclockwiseorder
static inline bool ccw(const Point2d& a, const Point2d& b, const Point2d& c) {
        return tri_area(a, b, c) > 0;
}


static inline bool right_of(const Point2d& x, Edge* e) {
        return ccw(x, e->Dest2d(), e->Org2d());
}


static inline bool left_of(const Point2d& x, Edge* e) {
        return ccw(x, e->Org2d(), e->Dest2d());
}


// A predicate that determines if the point x is on the edge e.
// The point is considered on if it is in the EPS-neighborhood
// of the edge.
static bool on_edge(const Point2d& x, Edge* e) {
        Real t1, t2, t3;
        t1 = (x - e->Org2d()).norm();
        t2 = (x - e->Dest2d()).norm();
        if ((t1 < RealEps) || (t2 < RealEps)) return true;

        t3 = (e->Org2d() - e->Dest2d()).norm();
        if ((t1 > t3) || (t2 > t3)) return false;

        Line line(e->Org2d(), e->Dest2d());
        return fabs(line.eval(x)) < RealEps;
}


/************* Topological Operations for Delaunay Diagrams *****************/

// Add a new edge e connecting the destination of a to the
// origin of b, in such a way that all three have the same
// left face after the connection is complete.
// Additionally, the data pointers of the new edge are set.
Edge* Delaunay::connect(Edge* a, Edge* b) {
        Edge* e = add_edge();
        splice(e, a->Lnext());
        splice(e->Sym(), b);
        e->end_points(a->Dest(), b->Org());
        return e;
}


// Essentially turns edge e counterclockwise inside its enclosing
// quadrilateral. The data pointers are modified accordingly.
void Delaunay::swap(Edge* e) {
        Edge* a = e->Oprev();
        Edge* b = e->Sym()->Oprev();
        splice(e, a);
        splice(e->Sym(), b);
        splice(e, a->Lnext());
        splice(e->Sym(), b->Lnext());
        e->end_points(a->Dest(), b->Dest());
}


// Initialize a subdivision to the line defined by the points a, b.
void Delaunay::init_2(const Point2d& a, const Point2d& b, Edge* &el,
                      Edge* &er) {
        Point2d *da = add_point(a), *db = add_point(b);
        Edge *ea = add_edge();
        ea->end_points(da, db);
        startingEdge = ea;
        el = ea;
        er = ea->Sym();
}


// Initialize a subdivision to the triangle defined by the points a, b, c.
void Delaunay::init_3(const Point2d &a, const Point2d &b, const Point2d &c,
                      Edge *&el, Edge *&er) {
        Point2d *da = add_point(a), *db = add_point(b), *dc = add_point(c);
        // Assume a, b, and c are in sorted order
        // Create edges ea and eb connecting da to db and db to dc
        Edge *ea = add_edge(), *eb = add_edge();
        splice(ea->Sym(), eb);
        ea->end_points(da, db);
        eb->end_points(ea->Dest(), dc);
        // Now close the triangle
        Edge *ec = connect(eb, ea);
        el = ea;
        er = eb->Sym();
        if (ccw(a, c, b)) {
                el = ec->Sym();
                er = ec;
        }
}


Delaunay::Delaunay(const Point2d &a, const Point2d &b, const Point2d &c) {
        Edge *el, *er;
        init_3(a, b, c, el, er);
        startingEdge = el;
}


Delaunay::Delaunay(std::list<Point2d>& point) {
        Edge *el, *er;
        init_dc(point, el, er);
        startingEdge = el;
}


Delaunay::~Delaunay() {
        for (Point2d *p : pointList) delete p;
        for (QuadEdge *qe : qeList) delete qe;
}


// Allocate a Point2d and add it to pointList
Point2d* Delaunay::add_point(const Point2d& a) {
        Point2d *da = new Point2d(a);
        pointList.push_back(da);
        return da;
}


// Allocate an edge and add it to qeList
Edge* Delaunay::add_edge() {
        Edge *eb = make_edge();
        qeList.push_back(eb->Qedge());
        return eb;
}


// Remove an edge and deallocate the memory
void Delaunay::remove_edge(Edge *e) {
        QEdgeList::iterator qeIter = locate_edge_iter(e);
        delete_edge(e);        // Deallocates the pointer
        qeList.erase(qeIter); // Removes the element from the list
}


// This step is relatively expensive: O(nEdge-1)
Delaunay::PointList::iterator Delaunay::locate_point_iter(Point2d* p) {
        return std::find(pointList.begin(), pointList.end(), p);
}


// This step is relatively expensive: O(nEdge-1)
Delaunay::QEdgeList::iterator Delaunay::locate_edge_iter(Edge *e) {
        return std::find(qeList.begin(), qeList.end(), e->Qedge());
}


// Convert pointer distance to size_type
Delaunay::PointList::size_type Delaunay::locate_point_index(Point2d* p) {
        return std::distance(pointList.begin(), locate_point_iter(p));
}


// Convert pointer distance to size_type
Delaunay::QEdgeList::size_type Delaunay::locate_edge_index(Edge *e) {
        return std::distance(qeList.begin(), locate_edge_iter(e));
}


// Returns an edge e, s.t. either x is on e, or e is an edge of
// a triangle containing x. The search starts from startingEdge
// and proceeds in the general direction of x. Based on the
// pseudocode in Guibas and Stolfi (1985) p.121.
Edge* Delaunay::locate(const Point2d& x) const {
        Edge *e = startingEdge;
        while (true) {
                if ((x == e->Org2d()) || (x == e->Dest2d())) return e;
                else if (right_of(x, e))
                        e = e->Sym();
                else if (!right_of(x, e->Onext()))
                        e = e->Onext();
                else if (!right_of(x, e->Dprev()))
                        e = e->Dprev();
                else return e;
        }
}


/************* Divide and Conquer (DC) Delaunay Diagram *****************/

// Tests whether the edge e is above basel
static inline bool valid(Edge *e, Edge *basel) {
        return right_of(*e->Dest(), basel);
}


// Creates a Delaunay triangulation from the points listed in 'point',
// using the Divide and Conquer algorithm from Guibas
// and Stolfi (1985) Fig. 23 p. 114.
void Delaunay::init_dc(std::list<Point2d> &point, Edge*& el, Edge*& er) {
        if (point.size() == 2) {
                // Create an edge from one to the other point
                init_2(point.front(), point.back(), el, er);
        } else if (point.size() == 3) {
                // Create a triangle
                std::list<Point2d>::iterator pi = point.begin();
                init_3(*pi, *std::next(pi), *std::next(pi, 2), el, er);
        } else {
                // Split into two lists
                std::list<Point2d>           pointL, pointR;
                std::list<Point2d>::iterator mpit;
                mpit = std::next(point.begin(), (point.size() + 1) / 2);

                pointL.splice(pointL.begin(), point, point.begin(), mpit);
                pointR.splice(pointR.begin(), point, mpit, point.end());

                // Compute Delaunay recusively until one of the conditions above is met
                Edge *ldo, *ldi, *rdi, *rdo;
                init_dc(pointL, ldo, ldi);
                init_dc(pointR, rdi, rdo);

                // Compute the lower tangent of L and R
                while (true) {
                        if      (left_of(*rdi->Org(), ldi))
                                ldi = ldi->Lnext();
                        else if (right_of(*ldi->Org(), rdi))
                                rdi = rdi->Rprev();
                        else
                                break;
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
                                while (in_circle(*basel->Dest(), *basel->Org(),
                                                 *lcand->Dest(),
                                                 *lcand->Onext()->Dest())) {
                                        Edge *t = lcand->Onext();
                                        delete_edge(lcand);
                                        lcand = t;
                                }

                        // Symmetrically, locate the first R point to be hit, and delete
                        // R edges
                        Edge *rcand = basel->Oprev();
                        if (valid(rcand, basel))
                                while (in_circle(*basel->Dest(), *basel->Org(),
                                                 *rcand->Dest(),
                                                 *rcand->Oprev()->Dest())) {
                                        Edge *t = rcand->Oprev();
                                        delete_edge(rcand);
                                        rcand = t;
                                }

                        // If both lcand and rcand are invalid, then basel is the upper
                        // common tangent
                        if (!valid(lcand, basel) && !valid(rcand, basel))
                                break;

                        // The next cross edge is to be connected to either lcand.Dest
                        // or rcand.Dest If both are valid, then choose the appropriate
                        // one using the InCircle test.
                        if (!valid(lcand, basel) ||
                            (valid(rcand, basel) &&
                             in_circle(*lcand->Dest(), *lcand->Org(),
                                       *rcand->Org(),
                                       *rcand->Dest())))
                                // Add cross edge basel from rcand.Dest to basel.Dest
                                basel = connect(rcand, basel->Sym());
                        else
                                // Add cross edge basel from basel.Org to lcand.Dest
                                basel = connect(basel->Sym(), lcand->Sym());
                } // Merge loop

                // Return [ldo, rdo]
                el = ldo;
                er = rdo;
        }
}


/************* Iterative Delaunay Diagram *****************/

// Inserts a new point into a subdivision representing a Delaunay
// triangulation, and fixes the affected edges so that the result
// is still a Delaunay triangulation. This is based on the
// pseudocode from Guibas and Stolfi (1985) p.120, with slight
// modifications and a bug fix.
void Delaunay::insert_point(const Point2d& x) {
        Edge *e = locate(x);
        if ((x == e->Org2d()) || (x == e->Dest2d())) {
                // point is already in
                return;
        } else if (on_edge(x, e)) {
                e = e->Oprev();
                remove_edge(e->Onext());
        }

        // Connect the new point to the vertices of the containing
        // triangle (or quadrilateral, if the new point fell on an
        // existing edge.)
        Edge    *base = add_edge();
        Point2d *dx = add_point(x);
        base->end_points(e->Org(), dx);
        splice(base, e);
        startingEdge = base;
        do {
                base = connect(e, base->Sym());
                e = base->Oprev();
        } while (e->Lnext() != startingEdge);

        // Examine suspect edges to ensure that the Delaunay condition
        // is satisfied.
        do {
                Edge* t = e->Oprev();
                if (right_of(t->Dest2d(), e) &&
                    in_circle(e->Org2d(), t->Dest2d(), e->Dest2d(), x)) {
                        swap(e);
                        e = e->Oprev();
                } else if (e->Onext() == startingEdge) {
                        // no more suspect edges
                        return;
                } else {
                        // pop a suspect edge
                        e = e->Onext()->Lprev();
                }
        } while (true);
}
