#ifndef CGL_TRIANGULATION_HPP
#define CGL_TRIANGULATION_HPP

#include <list>
#include <utility> // for std::pair
#include <vector>
#include <stack>

#include "cgl/geom2d.hpp"
#include "cgl/subdivision.hpp"
#include "cgl/point_set.hpp"

namespace cgl
{

class Hull {
public:
        typedef Point2d point_type;
        typedef std::pair<point_type, point_type> data_type;
        typedef std::list<data_type> container_type;
        typedef typename container_type::reverse_iterator reverse_iterator;
        typedef typename container_type::iterator iterator;

        Hull(std::vector<real> raw);

        iterator swap(iterator it, std::vector<real> raw);

        iterator begin();
        iterator end();
        reverse_iterator rbegin();
        reverse_iterator rend();

        const point_type& org() const;

private:
        point_type base;
        container_type hlist;
};


struct Neighbor {
        int neighbor;
        Direction dir;
};

class Triangulation : public Subdivision {
public:
        void swap(const edge_type&);
};

class Delaunay : public Triangulation {
public:
        Delaunay(PointSet&);

private:
        // Divide and conquer
        void init_dc(PointSet&, edge_type&, edge_type&, int);

        // Create a stack of merges that have to occur
        void create_merge_stack();

        // Merge within the processor
        void merge(edge_type&, const edge_type&, const edge_type&, edge_type&);

        // Merge with another processor
        void merge_proc(unsigned int, Direction, edge_type&, edge_type&);
        void merge_hull(const par::communicator&, unsigned int, edge_type&, const edge_type&, Hull&);
        void merge_hull(const par::communicator&, unsigned int, Hull&, const edge_type&, edge_type&);

        std::stack< std::vector<Neighbor> > merge_stack;
};

} // namespace cgl


#endif
