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

template <typename T1, typename T2, typename T3>
std::vector< std::list<T1> > split_multi(std::list<T1>& original_list, std::vector<T2>& bound, T3 pred)
{
        // Concept from @Yakk answer here:
        // http://stackoverflow.com/questions/23794679/splitting-an-stl-list-based-on-a-condition
        std::vector< std::list<T1> > result;

        std::list<T1> current;
        for (T2& b : bound) {
                do
                        current.splice(current.end(), original_list, original_list.begin());
                while (pred(b, current.back()));
                result.emplace_back(std::move(current));                
        }
        return result;
}

int main()
{
        PointSet ps(generate_random_points(10, 1e8));
        // // std::list<int> l = {1,2,3,4,5,6,7,8,9,10};
        // std::list<Point2d> l;
        // l.emplace_back(Point2d(0.1, 1.0));
        // l.emplace_back(Point2d(0.2, 1.0));
        // l.emplace_back(Point2d(0.3, 1.0));
        // l.emplace_back(Point2d(0.4, 1.0));
        // l.emplace_back(Point2d(0.5, 1.0));
        // l.emplace_back(Point2d(0.6, 1.0));
        // l.emplace_back(Point2d(0.7, 1.0));
        // l.emplace_back(Point2d(0.8, 1.0));
        // l.emplace_back(Point2d(0.9, 1.0));
        // l.emplace_back(Point2d(1.0, 1.0));
        // l.emplace_back(Point2d(1.1, 1.0));
        // l.emplace_back(Point2d(1.2, 1.0));
        // l.emplace_back(Point2d(1.3, 1.0));
        // std::vector<real> b = {0.3, 0.7, 0.9, 1.2};

        // auto a = split_multi(l, b, [](real b, const Point2d& p) { return b > p[0]; });
        // for (unsigned int i=0; i<a.size(); i++) {
        //         for (auto& b : a[i]) std::cout << b << " ";
        //         std::cout << std::endl;
        // }

        if (ps.size() != 10) {
                std::cerr << "Size test failed." << std::endl;
                return 1;
        }

        ps.partition_1d(0);
        // write_csv(ps, "point_out");

        // TODO: Add a real test

        // write_csv(ps, "out1");
        // ps.balance();
        // write_csv(ps, "out2");

        return 0;
}
