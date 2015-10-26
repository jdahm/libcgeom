#include "Geom2d.hh"
#include "Delaunay.hh"
#include <fstream>
#include <tuple>
#include <iostream>
#include <random>
#include <cmath>

typedef std::list<Point2d> PointList;
typedef std::tuple<PointList::size_type, PointList::size_type> PointConnect;
typedef std::list<PointConnect> PointConnectList;

void Read(const std::string& fileName, PointList& point, PointConnectList& connect)
{
  // Open the file
  std::ifstream fst(fileName, std::ios::in);

  // Number of points and lines
  PointList::size_type nPoint;
  PointConnectList::size_type nLine;
  fst >> nPoint;
  fst >> nLine;

  // Read in points
  {
    Real x, y;
    for (PointList::size_type i=0; i<nPoint; i++) {
      fst >> x;
      fst >> y;
      point.push_back(Point2d(x, y));
    }
  }

  // Read in edges
  {
    PointList::size_type j, k;
    for (PointConnectList::size_type i=0; i<nLine; i++) {
      fst >> j;
      fst >> k;
      connect.push_back(PointConnect(j, k));
    }
  }

  fst.close();
}

std::list<Point2d> generate_random_points(unsigned long int n, unsigned int prec)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(-prec, prec);
  unsigned long int np = 0;
  std::list<Point2d> point;

  while (np < n) {
    Real x = static_cast<Real>(dis(gen))/static_cast<Real>(prec),
      y = static_cast<Real>(dis(gen))/static_cast<Real>(prec);
    if (std::sqrt(x*x + y*y) < 1) {
      point.push_back(Point2d(x, y));
      np++;
    }
  }

  return point;
}


int main()
{
  PointList point;
  PointConnectList connect;

  // // Read the points
  // Read("test.txt", point, connect);
  point = generate_random_points(1e3, 1e9);

  // Sort
  point.sort([](const Point2d& a, const Point2d& b)
             {
               // WARNING: Floating point equivalence check here
               // Not very safe, but seems to work for now
               if (a.x == b.x) return (a.y < b.y);
               else return (a.x < b.x);
             });

  // Delaunay (DD)
  Delaunay DT(point);

  // Point2d a(0.0, 0.0), b(0.0, 10.0), c(10.0, 0.0);

  // Delaunay DT(a, b, c);

  // DT.InsertSite(Point2d(2.0, 2.0));

  DT.Write("out2.txt");

  return 0;
}
