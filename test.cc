#include "Geom2D.hh"
#include "Subdivision.hh"

int main()
{

  Point2d a(0.0, 0.0), b(0.0, 10.0), c(10.0, 0.0);

  Subdivision S(a, b, c);

  return 0;
}
