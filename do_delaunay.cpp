#include <string>
#include <sstream>
#include <iostream>
#include <list>
#include <vector>
#include <random>

#include "par/environment.hpp"
#include "cgl/point_set.hpp"
#include "cgl/triangulation.hpp"

using cgl::real;
using cgl::Point2d;
using cgl::PointSet;

typedef std::list<Point2d> PointList;
typedef std::tuple<PointList::size_type, PointList::size_type> PointConnect;
typedef std::list<PointConnect> PointConnectList;

struct ProgramArgs {
        enum class SaveType { None, Txt, Vtu };
        std::string infile;
        SaveType save;
        std::string outfile;
};

ProgramArgs parse_args(int argc, char *argv[]) {
        static const std::string usage("USAGE: DoDelaunay in_file ([-t | -v] out)");

        ProgramArgs args;

        const par::communicator& comm_world = par::comm_world();

        // Root proc parses the arguments
        if (argc < 2) comm_world.abort(usage, 1);

        {
                std::stringstream ss(argv[1]);
                ss >> args.infile;
        }

        args.save = ProgramArgs::SaveType::None;
        if (argc > 2) {
                const std::string osw(argv[2]);
                if (osw == "-t") args.save = ProgramArgs::SaveType::Txt;
                else if (osw == "-v") args.save = ProgramArgs::SaveType::Vtu;
        }

        // Output file name
        {
                std::stringstream ss(argv[3]);
                ss >> args.outfile;
        }

        std::cout << args.infile << " " << args.outfile << std::endl;

        return args;
}

void read_txt(const std::string& fileName, PointSet& ps,
              PointConnectList& connect) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(-1e8, 1e8);

        // Open the file
        std::ifstream fst(fileName, std::ios::in);

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

int main(int argc, char *argv[]) {
        const par::communicator& comm_world = par::comm_world();

        PointSet ps;
        PointList point;
        PointConnectList connect;
        ProgramArgs args;

        // Root proc parses the arguments
        args = parse_args(argc, argv);

        if (file_exists(args.infile)) {
                if (comm_world.rank() == 0)
                        read_txt(args.infile, ps, connect);
        }
        else {
                std::stringstream ss;
                ss << comm_world.rank();
                std::string rank;
                ss >> rank;
                std::string file(args.infile);
                file = file + "_" + rank + ".txt";
                read_txt(file, ps, connect);
        }

        std::cout << ps.size() << std::endl;

        // Distribute the points so they are evenly distributed among the
        // processors in a line
        ps.distribute(cgl::PSTopology::Unary, cgl::LBMethod::S2A2, 0);

        // Delaunay (DC)
        cgl::Delaunay DT(ps);

        if (args.save != ProgramArgs::SaveType::None) {
                if (args.save == ProgramArgs::SaveType::Txt)
                        DT.write_txt(args.outfile);
                // else if (args.save == ProgramArgs::SaveType::Vtu)
                //         DT.write_vtu(args.outFilePrefix);
                else
                        comm_world.abort("Unknown output format", 1);
        }

        return 0;
}
