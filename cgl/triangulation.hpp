#ifndef CGL_TRIANGULATION_HPP
#define CGL_TRIANGULATION_HPP

#include <stack>

#include "cgl/geom2d.hpp"
#include "cgl/subdivision.hpp"
#include "cgl/point_set.hpp"

namespace cgl
{

struct MergeInfo {
        par::communicator comm;
        unsigned int neighbor;
        Direction neighbor_dir;
        bool active;
};

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

        void init_dc(PointSet&, edge_type&, edge_type&, int);

        void create_merge_stack(const PointSet&, const par::communicator&);

        void proc_merge(const par::communicator&, unsigned int, Direction);

        std::stack<MergeInfo> merge_stack;
        AABB2d bounding_box;
};

} // namespace cgl


#endif
