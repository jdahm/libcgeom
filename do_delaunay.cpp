#include <string>
#include <sstream>
#include <iostream>
#include <list>
#include <vector>
#include <random>

#include "par/environment.hpp"
#include "cgl/point_set.hpp"
#include "cgl/triangulation.hpp"

using cgl::ProcTopology;
using cgl::real;
using cgl::Point2d;
using cgl::PointSet;

typedef std::list<Point2d> PointList;
typedef std::tuple<PointList::size_type, PointList::size_type> PointConnect;
typedef std::list<PointConnect> PointConnectList;

struct ProgramArgs {
        enum class InitType { File, Random };
        enum class SaveType { Txt, Vtu };
        InitType init;
        // Could be a union, but this deletes the default constructor
        unsigned int nPoints;
        std::string inFileName;
        SaveType save;
        std::string outFilePrefix;
};

ProgramArgs parse_args(int argc, char *argv[]) {
        static const std::string usage("USAGE: DoDelaunay [-f | -n] arg [-t | -v] out");

        ProgramArgs args;

        const par::communicator& comm_world = par::comm_world();

        // Root proc parses the arguments
        if (argc != 5) comm_world.abort(usage, 1);

        const std::string isw(argv[1]);
        std::stringstream iss(argv[2]);
        if (isw == "-f") {
                args.init = ProgramArgs::InitType::File;
                iss >> args.inFileName;
        }
        else if (isw == "-n") {
                args.init = ProgramArgs::InitType::Random;
                iss >> args.nPoints;
        }
        else {
                comm_world.abort("Command line parsing error. " + usage, 1);
        }

        const std::string osw(argv[3]);
        if (osw == "-t") args.save = ProgramArgs::SaveType::Txt;
        else if (osw == "-v") args.save = ProgramArgs::SaveType::Vtu;

        // Output file name
        args.outFilePrefix = argv[4];

        return args;
}

void read_txt(const std::string& fileName, PointSet& ps,
              PointConnectList& connect) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(-1e8, 1e8);

        // Open the file
        std::ifstream fst(fileName, std::ios::in);
        const real delta_angle = 0.01;

        // Number of points and lines
        PointList::size_type nPoint;
        PointConnectList::size_type nLine;
        fst >> nPoint;
        fst >> nLine;

        // Read in points
        for (PointList::size_type i = 0; i<nPoint; i++) {
                real x, y;
                fst >> x;
                fst >> y;
                // real xp, yp;
                // xp = x * std::cos(delta_angle) - y * std::sin(delta_angle);
                // yp = x * std::sin(delta_angle) + y * std::cos(delta_angle);
                // x = std::sqrt(x*x + y*y) * std::cos(std::atan(y/x) + angle);
                // y = std::sqrt(x*x + y*y) * std::sin(std::atan(y/x) + angle);
                ps.add(Point2d(x, y));
        }

        // Read in edges
        for (PointConnectList::size_type i = 0; i<nLine; i++) {
                PointList::size_type j, k;
                fst >> j;
                fst >> k;
                connect.push_back(PointConnect(j, k));
        }

        fst.close();
}

bool file_exists(const std::string& fileName)
{
        std::ifstream fst(fileName);
        if (fst) return true;
        else return false;
}


PointSet generate_random_pointset(unsigned long int n, unsigned int prec) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(-prec, prec);
        std::list<Point2d> point;
        PointSet ps;

        for (unsigned long int np = 0; np < n; np++) {
                real x = static_cast<real>(dis(gen)) / static_cast<real>(prec);
                real y = static_cast<real>(dis(gen)) / static_cast<real>(prec);
                ps.add(Point2d(x, y));
        }

        return ps;
}

int main(int argc, char *argv[]) {
        const par::communicator& comm_world = par::comm_world();

        PointSet ps;
        PointList point;
        PointConnectList connect;
        ProgramArgs args;

        // Root proc parses the arguments
        args = parse_args(argc, argv);

        if (args.init == ProgramArgs::InitType::File) {
                std::string file(args.inFileName);
                std::stringstream ss;
                ss << comm_world.rank();
                std::string rank;
                ss >> rank;
                if (!file_exists(file)) {
                        file = file + "_" + rank + ".txt";
                        read_txt(file, ps, connect);
                }
                else {
                        if (comm_world.rank() == 0)
                                read_txt(file, ps, connect);
                }
        }
        else if (args.init == ProgramArgs::InitType::Random)
                ps = generate_random_pointset(args.nPoints, 1e8);
        else
                throw std::runtime_error("Unknown input type");

        // Distribute the points so they are evenly distributed among the
        // processors in a line
        ps.distribute(ProcTopology::Line);

        // Delaunay (DC)
        cgl::Delaunay DT(ps);
        
        if (args.save == ProgramArgs::SaveType::Txt)
                DT.write_txt(args.outFilePrefix);
        // else if (args.save == ProgramArgs::SaveType::Vtu)
        //         DT.write_vtu(args.outFilePrefix);
        else {
                comm_world.abort("Unknown output format", 1);
        }

        return 0;
}
