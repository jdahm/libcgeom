#include "QuadEdge.hh"
#include <iostream>

int main()
{
  Point2d *a = new Point2d(0.0, 0.0);
  Point2d *b = new Point2d(1.0, 0.0);
  Point2d *c = new Point2d(1.0, 1.0);
  Point2d *d = new Point2d(0.0, 1.0);

  Edge *ea = MakeEdge();
  Edge *eb = MakeEdge();
  Edge *ec = MakeEdge();
  Edge *ed = MakeEdge();
  Edge *ee = MakeEdge();

  ea->EndPoints(a, b);
  eb->EndPoints(b, c);
  ec->EndPoints(c, d);
  ed->EndPoints(ec->Dest(), ea->Org());
  ee->EndPoints(eb->Dest(), ea->Org());

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

  delete a;
  delete b;
  delete c;
  delete d;

  return 0;
}
