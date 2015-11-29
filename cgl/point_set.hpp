#ifndef CGL_POINT_SET_HPP
#define CGL_POINT_SET_HPP

#include <list>
#include <vector>
#include <sstream>
#include <fstream>
#include <zoltan.h>

#include "cgl/geom2d.hpp"
#include "par/environment.hpp"

namespace cgl
{

class PointSet {
        /*
          PointSet wraps a std::list<Point2d>. It uses this data structure
          because a list can be split in constant time and without additional
          copies. This speeds up divide and conquer algorithms.
        */

public:
        typedef Point2d point_type;
        typedef std::list<point_type> container_type;
        typedef typename container_type::size_type size_type;
        typedef typename container_type::iterator iterator;
        typedef typename container_type::const_iterator const_iterator;

        static constexpr size_type dim() { return 2; }

        PointSet();

        PointSet(const container_type &other);

        PointSet(container_type &&other);

        ~PointSet();

        PointSet(const PointSet &other);

        PointSet& operator=(const PointSet &other);

        // Only have a ForwardIterator, but this requires a RandomAccess iterator
        // const point_type& operator[](size_type i) const { return point[i]; }
        // point_type& operator[](size_type i) { return point[i]; }

        size_type size() const;

        PointSet split(iterator);

        PointSet halve();

        void sort(unsigned int);

        void distribute(ProcTopology);

        iterator begin();
        iterator end();
        const_iterator begin() const;
        const_iterator end() const;

        const point_type& front() const;
        const point_type& back() const;

        void add(const point_type&);

        iterator erase(iterator);

        size_type global_offset();

private:
        // Interactions with zoltan
        void z_init();
        void z_balance(ProcTopology);
        void z_destroy();

        void recompute_global_offset();

        void partition_1d(unsigned int);

        container_type point;
        bool global_valid;
        size_type global_offset_;
        struct Zoltan_Struct *zz;
        bool z_initialized;
        int sorted_dimension;
};

void write_csv(const PointSet &ps, const std::string &prefix);

template <typename T>
std::vector< std::list<T> > split_multi(std::list<T>& original_list, std::vector<T>& bound);

} // namespace cgl

#endif
