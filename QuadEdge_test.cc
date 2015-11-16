#include "QuadEdge.hh"
#include <iostream>

int main() {
  Point2d *a = new Point2d(0.0, 0.0);
  Point2d *b = new Point2d(1.0, 0.0);
  Point2d *c = new Point2d(1.0, 1.0);
  Point2d *d = new Point2d(0.0, 1.0);

  Edge *ea = makeEdge();
  Edge *eb = makeEdge();
  Edge *ec = makeEdge();
  Edge *ed = makeEdge();
  Edge *ee = makeEdge();

  ea->endPoints(a, b);
  // Connect eb to ea->Dest()
  splice(ea->Sym(), eb); eb->endPoints(ea->Dest(), c);
  // Connect ec to eb->Dest()
  splice(eb->Sym(), ec); ec->endPoints(eb->Dest(), d);
  // Connect ed to ec->Sym() and ea
  splice(ec->Sym(), ed); splice(ed->Sym(), ea); ed->endPoints(ec->Dest(), a);
  // Connect ee to ec and ea -- this is a Delaunay::connect procedure
  splice(ee, eb->Lnext()); splice(ee->Sym(), ea); ee->endPoints(eb->Dest(), ea->Org());

  // Traverse the convex hull
  if (ea->Dnext() != eb->Sym()) {
    std::cerr << "Error: ea->Dnext() != eb->Sym()" << std::endl;
    return 1;
  }

  if (eb->Dnext() != ec->Sym()) {
    std::cerr << "Error: eb->Dnext() != ec->Sym()" << std::endl;
    return 1;
  }

  if (ec->Dnext() != ed->Sym()) {
    std::cerr << "Error: ec->Dnext() != ed->Sym()" << std::endl;
    return 1;
  }

  if (ed->Dprev() != ee) {
    std::cerr << "Error: ed->Dprev() != ee" << std::endl;
    return 1;
  }

  // Cross linkage -- test algebra
  // Rot / invRot
  if (ee->Rot()->Org() != ed->invRot()->Org()) {
    std::cerr << "Error: ee->Rot()->Org() != ed->invRot()->Org()" << std::endl;
    return 1;
  }

  // Sym
  if (ee->Sym()->Org() != a) {
    std::cerr << "Error: ee->Sym()->Org() != a" << std::endl;
    return 1;
  }

  // Onext
  if (ee->Onext() != eb->Sym()) {
    std::cerr << "Error: ee->Onext() != eb->Sym()" << std::endl;
    return 1;
  }

  // Oprev
  if (ee->Oprev() != ec) {
    std::cerr << "Error: ee->Oprev() != ec" << std::endl;
    return 1;
  }

  // Dnext
  if (ee->Dnext() != ed) {
    std::cerr << "Error: ee->Dnext() != ed" << std::endl;
    return 1;
  }

  // Dprev
  if (ee->Dprev() != ea->Sym()) {
    std::cerr << "Error: ee->Dprev() != ea->Sym()" << std::endl;
    return 1;
  }

  // Lnext
  if (ee->Lnext() != ea) {
    std::cerr << "Error: ee->Lnext() != ea" << std::endl;
    return 1;
  }

  // Lprev
  if (ee->Lprev() != eb) {
    std::cerr << "Error: ee->Lprev() != eb" << std::endl;
    return 1;
  }

  // Rnext
  if (ee->Rnext() != ec->Sym()) {
    std::cerr << "Error: ee->Rnext() != ec->Sym()" << std::endl;
    return 1;
  }

  // Rprev
  if (ee->Rprev() != ed->Sym()) {
    std::cerr << "Error: ee->Rprev() != ed->Sym()" << std::endl;
    return 1;
  }

  // Org / Dest
  if (ed->Org()  != d) {
    std::cerr << "Error: ed->Org() != d" << std::endl;
    return 1;
  }

  if (ed->Dest() != a) {
    std::cerr << "Error: ed->Dest() != a" << std::endl;
    return 1;
  }

  if (ee->Org()  != c) {
    std::cerr << "Error: ee->Org() != c" << std::endl;
    return 1;
  }

  if (ee->Dest()  != a) {
    std::cerr << "Error: ee->Dest() != a" << std::endl;
    return 1;
  }

  // Let the system clean up the memory

  return 0;
}
