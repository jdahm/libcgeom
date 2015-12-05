#ifndef CGL_TRIANGULATION_HPP
#define CGL_TRIANGULATION_HPP

#include <stack>

#include "cgl/geom2d.hpp"
#include "cgl/subdivision.hpp"
#include "cgl/point_set.hpp"

namespace cgl
{

struct MergeInfo {
        int neighbor;
        Direction neighbor_dir;
};

class Triangulation : public Subdivision {
public:
        void swap(const edge_type&);
};

class Delaunay : public Triangulation {
public:
        Delaunay(PointSet&);

private:
        void create_merge_stack();

        void merge(edge_type&, const edge_type&, const edge_type&, edge_type&);

        void init_dc(PointSet&, edge_type&, edge_type&, int);

        void proc_merge(unsigned int, Direction, edge_type&);

        std::stack<MergeInfo> merge_stack;
        AABB2d bounding_box;
};

} // namespace cgl


#endif
