// #include <cmath>
// #include <vector>
// #include <array>
// #include <numeric>
#include <algorithm>
#include <iomanip>

#include "cgl/point_set.hpp"
#include "cgl/geom2d.hpp"
#include "par/environment.hpp"

namespace cgl
{

PointSet::PointSet() : point(), sorted_dimension(-1) { }

PointSet::PointSet(const container_type& p) :
        point(p), sorted_dimension(-1) { }

PointSet::PointSet(container_type&& p) :
        point(p), sorted_dimension(-1) { }

PointSet::PointSet(const PointSet& other) :
        point(other.point), sorted_dimension(-1) { }

PointSet& PointSet::operator=(const PointSet& other) {
        point = other.point;
        sorted_dimension = other.sorted_dimension;
        return *this;
}

PointSet PointSet::split(PointSet::iterator it) {
        // Create a new container
        container_type pointR;

        // Take the right side of the list
        pointR.splice(pointR.begin(), point, it, point.end());

        // Return this->begin() in it
        it = begin();

        // Return a PointSet
        return PointSet(pointR);
}

PointSet PointSet::halve() {
        return split(std::next(begin(), (size() + 1) / 2));
}

void PointSet::sort(unsigned int dimen)
{
        if (sorted_dimension != static_cast<int>(dimen))
                sort_point2d_list(dimen, point);
        sorted_dimension = dimen;
}

void PointSet::distribute(PSTopology top, LBMethod method, unsigned int dimen)
{
        // Balance the set
        balance_set(top, method, dimen, point);
}

PointSet::size_type PointSet::size() const { return point.size(); }

PointSet::iterator PointSet::begin() { return point.begin(); }
PointSet::iterator PointSet::end() { return point.end(); }
PointSet::const_iterator PointSet::begin() const { return point.begin(); }
PointSet::const_iterator PointSet::end() const { return point.end(); }

const PointSet::point_type& PointSet::front() const { return point.front(); }
const PointSet::point_type& PointSet::back() const { return point.back(); }

void PointSet::add(const point_type& p) { point.push_back(p); sorted_dimension = -1; }

PointSet::iterator PointSet::erase(iterator it) { return point.erase(it); }

void write_txt(const PointSet &ps, const std::string &prefix)
{
        const par::communicator &comm_world = par::comm_world();

        std::stringstream ss;
        ss << comm_world.rank();

        std::string rank;
        ss >> rank;
        std::ofstream fst(prefix + "_" + rank + ".txt", std::ios::out);

        fst << ps.size() << " 0" << std::endl;
        for (PointSet::const_iterator it=ps.begin(); it!=ps.end(); ++it) {
                const PointSet::point_type &a = *it;
                fst << a[0] << " " << a[1] << std::endl;
        }

        fst.close();
}

void read_txt(PointSet& ps, const std::string& prefix)
{
        const par::communicator &comm_world = par::comm_world();

        std::stringstream ss;
        ss << comm_world.rank();

        std::string rank;
        ss >> rank;
        std::ifstream fst(prefix + "_" + rank + ".txt", std::ios::in);

        PointSet::size_type num_point;
        fst >> num_point;

        PointSet::size_type i = 0;
        while (fst.good()) {
                real x, y;
                fst >> x;
                fst.ignore();
                fst >> y;
                if (fst.eof() || i >= num_point) break;
                ps.add(PointSet::point_type(x, y));
        }

        fst.close();
}

} // namespace cgl
