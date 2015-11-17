#ifndef DELAUNAY_HH
#define DELAUNAY_HH

#include "Geom2d.H"
#include "QuadEdge.H"
#include <list>
#include <tuple>

class Delaunay {
        // List types
        typedef std::list<QuadEdge*> QEdgeList;
        typedef std::list<Point2d*> PointList;

        // Data members
        // Global lists of primal points and edges
        QEdgeList qeList;
        PointList pointList;
        // Keep a pointer to an edge from QuadEdgeList from whick to begin walking
        Edge      *startingEdge;

        // Operators to add points and edges updating pointList and qeList above
        Point2d *add_point(const Point2d&);

        Edge *add_edge();

        // Operators to remove points and edges, updating pointList and qeList above
        void remove_edge(Edge*);

        // Operators to locate a point or edge given a pointer
        // Output iterator
        PointList::iterator locate_point_iter(Point2d*);

        QEdgeList::iterator locate_edge_iter (Edge*);

        // Output index
        PointList::size_type locate_point_index(Point2d*);

        QEdgeList::size_type locate_edge_index (Edge*);

        // Locate an edge based on coordinates. Attempts to do a smart
        // search based on direction while walking the subdivision
        Edge *locate(const Point2d&) const;

        // Initialize a triangulation (used by constructors)
        void init_2(const Point2d&, const Point2d&, Edge*&, Edge*&);

        void init_3(const Point2d&,
                    const Point2d&,
                    const Point2d&,
                    Edge*&,
                    Edge*&);

        void init_dc(std::list<Point2d>&, Edge*&, Edge*&);

        // Topological Operators for Delaunay Diagram
        Edge* connect(Edge*, Edge*);

        void swap(Edge*);

public:
        // Constructors
        Delaunay(std::list<Point2d>&);
        Delaunay(const Point2d &, const Point2d &, const Point2d &);

        ~Delaunay();

        // Insert a new point and update Delaunay Diagram
        void insert_point(const Point2d&);

        // Simple text file output
        void write_txt(const std::string&);

        // Paraview friendly output
        void write_vtu(const std::string&);

};


#endif
