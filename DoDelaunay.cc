#include "Parallel.hh"
#include "Geom2d.hh"
#include "Delaunay.hh"
#include <fstream>
#include <tuple>
#include <iostream>
#include <random>
#include <cmath>
#include <sstream>

typedef std::list<Point2d> PointList;
typedef std::tuple<PointList::size_type, PointList::size_type> PointConnect;
typedef std::list<PointConnect> PointConnectList;

void readTxt(const std::string& fileName, PointList& point, PointConnectList& connect)
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
  std::list<Point2d> point;

  for (unsigned long int np = 0; np < n; np++) {
    Real x = static_cast<Real>(dis(gen))/static_cast<Real>(prec);
    Real y = static_cast<Real>(dis(gen))/static_cast<Real>(prec);
    point.push_back(Point2d(x, y));
  }

  return point;
}

template<class Cont>
void queueTransfer(Cont& v, unsigned int rankDest)
{
  // Nothing here, for now.
}


int main(int argc, char *argv[])
{
  unsigned int myRank, nProc;
  PointList point;
  PointConnectList connect;

  initParallel(&argc, &argv);

  myRank = rank();
  nProc = numProcs();

  if (myRank == 0) {
    if (argc < 3) {
      throw std::runtime_error("USAGE: DoDelaunay [-f | -n] arg");
    }
    std::string sw(argv[1]);
    if (sw == "-f") {
      readTxt(std::string(argv[2]), point, connect);
    }
    else if (sw == "-n") {
      std::stringstream ss;
      ss << argv[2];
      unsigned long int n;
      ss >> n;
      point = genRandomPoints(n, 1e8);
    }

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

    for (unsigned int i=1; i<nProc; i++){
      std::list<Point2d> outpoint;
      std::list<Point2d>::iterator startit, endit;
      startit = std::next(point.begin(), i*(point.size()+1)/nProc);
      endit   = std::next(point.begin(), std::min((i+1)*(point.size()+1)/nProc, point.size()));
      outpoint.splice(outpoint.begin(), point, startit, endit);
      queueTransfer(outpoint, i);
    }
  }

  // Delaunay (DC)
  Delaunay DT(point);

  // DT.writeTxt("out2.txt");
  DT.writeVtu("out");

  finalizeParallel();

  return 0;
}
