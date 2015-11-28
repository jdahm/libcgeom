#include <iostream>
#include <cmath>

#include "cgl/geom2d.hpp"

int main()
{
        using cgl::Point2d;
        using cgl::Line2d;

        Point2d p;

        if (p.distance() != 0.) {
                std::cerr << "Origin distance failed." << std::endl;
                return 1;
        }

        if (p.distance(Point2d(1, 2)) != std::sqrt(5)) {
                std::cerr << "Point distance failed." << std::endl;
                return 1;
        }

        if (p == Point2d({0.0 + 2*cgl::real_eps, 0.0})) {
                std::cerr << "Inequality test failed." << std::endl;
                return 1;
        }

        if (p != Point2d({0.0, 0.0})) {
                std::cerr << "Inequality test failed." << std::endl;
                return 1;
        }

        Line2d l1({0, 1}, {1, 1});
        Line2d l2({1, 1}, {0, 0});

        if (intersect(l1, l2) != Point2d(1, 1)) {
                std::cerr << "Line intersect test failed." << std::endl;
                return 1;
        }

        return 0;
}
