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

void ReadTxt(const std::string& fileName, PointList& point, PointConnectList& connect)
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

std::list<Point2d> genRandomPoints(unsigned long int n, unsigned int prec)
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
  point = genRandomPoints(1e3, 1e9);

  // Sort
  point.sort([](const Point2d& a, const Point2d& b)
             {
               // Maybe this is slightly safer...
               if (std::abs(a.x - b.x) < RealEps) return (a.y < b.y);
               else return (a.x < b.x);
             });

  // Remove duplicates (has to be already sorted)
  point.unique([](const Point2d& a, const Point2d& b)
               { return ((std::abs(a.x - b.x) < RealEps) &&
                         (std::abs(a.y - b.y) < RealEps)); });

  // Delaunay (DC)
  Delaunay DT(point);

  // DT.Write("out2.txt");
  DT.WriteVtu("out");
  
  return 0;
}
