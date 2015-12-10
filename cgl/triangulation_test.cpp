#include "cgl/point_set.hpp"
#include "cgl/triangulation.hpp"
#include <random>

using cgl::real;
using cgl::PointSet;
using cgl::Point2d;
using cgl::Delaunay;

#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <cstdio>

// Source: http://stackoverflow.com/questions/15118661/in-c-whats-the-fastest-way-to-tell-whether-two-string-or-binary-files-are-di
template <typename InputIterator1, typename InputIterator2>
bool range_equal(InputIterator1 first1, InputIterator1 last1,
                 InputIterator2 first2, InputIterator2 last2) {
        while (first1 != last1 && first2 != last2) {
                if (*first1 != *first2) return false;
                ++first1;
                ++first2;
        }
        return (first1 == last1) && (first2 == last2);
}


bool compare_files(const std::string& filename1, const std::string& filename2) {
        std::ifstream file1(filename1);
        std::ifstream file2(filename2);

        std::istreambuf_iterator<char> begin1(file1);
        std::istreambuf_iterator<char> begin2(file2);

        std::istreambuf_iterator<char> end;

        return range_equal(begin1, end, begin2, end);
}


bool test_txt(const std::string &fileName) {
        typedef std::list<Point2d> PointList;
        typedef std::tuple<PointList::size_type,
                           PointList::size_type> PointConnect;
        typedef std::list<PointConnect> PointConnectList;

        PointSet ps;

        // Open the file
        std::ifstream fst(fileName, std::ios::in);

        // Number of points and lines
        PointList::size_type nPoint;
        PointConnectList::size_type nLine;
        fst >> nPoint;
        fst >> nLine;

        // Read in points
        for (PointList::size_type i=0; i<nPoint; i++) {
                real x, y;
                fst >> x;
                fst >> y;
                ps.add(Point2d(x, y));
        }

        Delaunay dt(ps);

        // Write a temporary file
        dt.write_txt("._out");

        // Compare it
        bool equal = compare_files("._out_0.txt", fileName);

        // Remove the temporary file
        std::remove("._out_0.txt");

        return equal;
}


int main() {
        if (!test_txt("test.txt")) {
                std::cerr << "test.txt Delaunay test failed" << std::endl;
                return 1;
        }

        return 0;
}
