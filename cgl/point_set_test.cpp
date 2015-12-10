#include "par/environment.hpp"
#include "cgl/point_set.hpp"
#include <random>

using cgl::real;
using cgl::PointSet;
using cgl::Point2d;

typename PointSet::container_type
generate_random_points(unsigned long int n, unsigned int prec)
{
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(-prec, prec);
        typename PointSet::container_type point;

        for (unsigned long int np = 0; np < n; np++) {
                real x = static_cast<real>(dis(gen)) / static_cast<real>(prec);
                real y = static_cast<real>(dis(gen)) / static_cast<real>(prec);
                point.push_back(Point2d(x, y));
        }

        return point;
}

int main()
{
        // TODO: Add a real test

        PointSet::container_type point = generate_random_points(15, 1e8);
        PointSet ps(std::move(point));

        ps.distribute(cgl::PSTopology::Unary, cgl::LBMethod::S2A2, 0);

        write_txt(ps, "t1s");

        return 0;
}
