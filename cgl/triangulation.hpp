#ifndef CGL_TRIANGULATION_HPP
#define CGL_TRIANGULATION_HPP

#include "cgl/subdivision.hpp"
#include "cgl/point_set.hpp"

namespace cgl
{

class Triangulation : public Subdivision {
public:
        void swap(const edge_type&);

        void write_txt(const std::string& filePrefix);
};

class Delaunay : public Triangulation {
public:
        Delaunay(PointSet&);

private:
        void merge(edge_type&, const edge_type&, const edge_type&, edge_type&);

        void init_dc(PointSet&, edge_type&, edge_type&);
};

} // namespace cgl


#endif
