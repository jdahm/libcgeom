#ifndef CGL_POINT_SET_HPP
#define CGL_POINT_SET_HPP

#include <list>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>

#include "cgl/geom2d.hpp"
#include "cgl/load_balance.hpp"
#include "par/communicator.hpp"

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

        friend void read_csv(PointSet& ps, const std::string& prefix);

        PointSet();

        PointSet(const container_type&);

        PointSet(container_type&&);

        PointSet(const PointSet &other);

        PointSet& operator=(const PointSet &other);

        // Only have a ForwardIterator, but this requires a RandomAccess iterator
        // const point_type& operator[](size_type i) const { return point[i]; }
        // point_type& operator[](size_type i) { return point[i]; }

        size_type size() const;

        PointSet split(iterator);

        PointSet halve();

        void sort(unsigned int);

        void distribute(PSTopology, LBMethod, unsigned int);

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
        void recompute_global_offset();

        container_type point;
        int sorted_dimension;
};

void write_txt(const PointSet&, const std::string&);
void read_txt(PointSet&, const std::string&);

} // namespace cgl

#endif
