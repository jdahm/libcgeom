#include "Geom2d.hh"
#include "Delaunay.hh"

int main()
{

  Point2d a(0.0, 0.0), b(0.0, 10.0), c(10.0, 0.0);

  Delaunay DT(a, b, c);

  DT.InsertSite(Point2d(2.0, 2.0));

  DT.Write("out.txt");

  DT.WriteVtuFiles("out");
  
  return 0;
}
