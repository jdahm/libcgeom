#include <iostream>

#include "cgl/subdivision.hpp"
#include "cgl/geom2d.hpp"

using cgl::Subdivision;
using cgl::Point2d;

int main()
{
        Subdivision s;

        Point2d p1(0, 0), p2(1, 0), p3(0, 10);

        Subdivision::edge_type a = s.add_edge(p1, p2);

        Subdivision::edge_type b = s.extend_edge(a, p3);

        Subdivision::edge_type c = s.connect_edges(b, a);

        s.remove_edge(c);

        if (s.num_points() != 3 || s.num_edges() != 2) {
                std::cerr << "Failed num_points() and num_edges() test." << std::endl;
                return 1;
        }

        return 0;
}
